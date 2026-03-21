"""全局阻尼 Newton–Raphson 求解器（Hamel 1999 §2.3 Eq. 2.9）。

将所有 cell 状态打包为全局向量，用有限差分构造稀疏 Jacobian，
阻尼 NR 一次性更新整个反应器。边界条件（含循环返料）在残差前施加，
故 FD 会体现顶格→底格的耦合（Nebenelemente）。

残差行缩放采用**按方程类型的物理参考量**（非当前 ‖F‖ 动态缩放），
使气相/固相/能量方程在 ‖F̂‖ 中可比，避免 NR 只“看见”能量方程。

Source: Hamel (1999) §2.3, p.20; Bild 2.2
"""

from __future__ import annotations

from collections.abc import Callable
from typing import List

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import N_GAS

_FD_EPS = 1e-6


def n_var(cell: Cell) -> int:
    """每格未知量：N_d + N_b + m_solid(nk*4) + T。"""
    return 2 * N_GAS + cell.solid.n_size_classes * N_SOLID_COMP + 1


def pack_cell(cell: Cell) -> np.ndarray:
    return np.concatenate([
        cell.N_d,
        cell.N_b,
        cell.m_solid.flatten(),
        np.array([cell.T], dtype=np.float64)
    ])


def unpack_cell(x: np.ndarray, cell: Cell) -> None:
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    np.maximum(x[:N_GAS], 0.0, out=cell.N_d)
    np.maximum(x[N_GAS : 2 * N_GAS], 0.0, out=cell.N_b)
    m_flat = np.maximum(x[2 * N_GAS : 2 * N_GAS + nv_sol], 0.0)
    cell.m_solid[:] = m_flat.reshape((nk, N_SOLID_COMP))
    cell.T = float(np.clip(x[-1], 300.0, 2500.0))


def pack_reactor(cells: List[Cell]) -> np.ndarray:
    return np.concatenate([pack_cell(c) for c in cells])


def unpack_reactor(x: np.ndarray, cells: List[Cell]) -> None:
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        unpack_cell(x[offset : offset + nv], cell)
        offset += nv


def cell_offsets(cells: List[Cell]) -> List[int]:
    offs: List[int] = [0]
    for c in cells:
        offs.append(offs[-1] + n_var(c))
    return offs


def global_residual(
    x: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
) -> np.ndarray:
    """F(x)：解包 → 边界条件（含循环）→ 各 cell 残差拼接。"""
    unpack_reactor(x, cells)
    apply_bc_fn()
    return np.concatenate([c.residuals() for c in cells])


def build_equation_scales(
    cells: List[Cell],
    ref_gas_mol_s: float,
    ref_solid_kg_s: float,
    ref_energy_W: float,
) -> np.ndarray:
    """按方程类型构造静态行缩放（物理参考量级，非当前残差）。

    - 气相守恒 [mol/s]：ref_gas_mol_s
    - 固相守恒 [kg/s]：ref_solid_kg_s
    - 能量 [W]：ref_energy_W
    """
    rg = max(float(ref_gas_mol_s), 1e-9)
    rs = max(float(ref_solid_kg_s), 1e-12)
    re = max(float(ref_energy_W), 1.0)

    n_total = sum(n_var(c) for c in cells)
    scale = np.empty(n_total, dtype=np.float64)
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        nk = cell.solid.n_size_classes
        nv_sol = nk * N_SOLID_COMP
        scale[offset : offset + N_GAS] = rg
        scale[offset + N_GAS : offset + 2 * N_GAS] = rg
        scale[offset + 2 * N_GAS : offset + 2 * N_GAS + nv_sol] = rs
        scale[offset + 2 * N_GAS + nv_sol] = re
        offset += nv
    return scale


def build_jacobian_fd(
    x: np.ndarray,
    F0: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    eq_scale: np.ndarray,
) -> sp.csr_matrix:
    """列有限差分构造稀疏 Jacobian（行与 F̂=F/eq_scale 一致）。"""
    n_total = len(x)
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []

    F0s = F0 / eq_scale

    for j in range(n_total):
        # 修正步长逻辑：对于微量组分，使用其绝对值相关的步长，而非硬编码的 max(1, abs(x))
        xj = float(x[j])
        h = _FD_EPS * (abs(xj) + 1e-12)
        
        x_p = x.copy()
        x_p[j] += h
        
        # 只在扰动后更新状态并计算残差
        unpack_reactor(x_p, cells)
        apply_bc_fn()
        F_p = np.concatenate([c.residuals() for c in cells])
        
        dF = (F_p / eq_scale - F0s) / h
        
        # 快速寻找非零元并记录（稀疏矩阵构造）
        indices = np.where(np.abs(dF) > 1e-14)[0]
        for i in indices:
            rows.append(int(i))
            cols.append(j)
            vals.append(float(dF[i]))

    # 计算结束后恢复原始状态
    unpack_reactor(x, cells)
    apply_bc_fn()
    return sp.csr_matrix((vals, (rows, cols)), shape=(n_total, n_total))


def _rms_norm(F_hat: np.ndarray) -> float:
    """缩放残差向量的 RMS（用于收敛判据）。"""
    if len(F_hat) == 0:
        return 0.0
    return float(np.sqrt(np.mean(F_hat**2)))


def _clip_dx(
    dx: np.ndarray,
    cells: List[Cell],
    ref_gas_mol_s: float,
    ref_solid_kg_s: float,
) -> np.ndarray:
    """严格限制 Newton 步长，防止物理量跳变导致发散。"""
    rg_limit = 0.2 * max(float(ref_gas_mol_s), 1.0)
    rs_limit = 0.2 * max(float(ref_solid_kg_s), 0.1)
    dT_limit = 50.0  # 单次迭代温度变化限制在 50K 以内

    out = dx.copy()
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        nk = cell.solid.n_size_classes
        nv_sol = nk * N_SOLID_COMP
        # 气相
        out[offset : offset + 2 * N_GAS] = np.clip(
            out[offset : offset + 2 * N_GAS], -rg_limit, rg_limit
        )
        # 固相
        out[offset + 2 * N_GAS : offset + 2 * N_GAS + nv_sol] = np.clip(
            out[offset + 2 * N_GAS : offset + 2 * N_GAS + nv_sol], -rs_limit, rs_limit
        )
        # 温度
        out[offset + nv - 1] = float(np.clip(out[offset + nv - 1], -dT_limit, dT_limit))
        offset += nv
    return out


def solve_global_nr(
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    ref_gas_mol_s: float = 80.0,
    ref_solid_kg_s: float = 1.0,
    ref_energy_W: float = 2e7,
    max_iter: int = 25,
    tol_rms: float = 0.01,
    lambda_init: float = 0.1,  # 初始步长更为保守
    n_damp_halvings: int = 8,
    lambda_min: float = 1.0 / 1024.0,
    verbose: bool = False,
) -> dict:
    """阻尼 Newton–Raphson：ĴΔx = −F̂（F̂ = F / 静态 eq_scale），再线搜索阻尼。

    Hamel Eq. 2.9（gedämpftes Newton-Verfahren）的工程实现。
    """
    x = pack_reactor(cells)
    apply_bc_fn()
    F = global_residual(x, cells, apply_bc_fn)
    scale = build_equation_scales(cells, ref_gas_mol_s, ref_solid_kg_s, ref_energy_W)

    n_total = len(F)
    n_gas = 2 * N_GAS * len(cells)
    n_sol = sum(c.solid.n_size_classes * N_SOLID_COMP for c in cells)
    n_ene = len(cells)

    F_hat = F / scale
    norm_F = _rms_norm(F_hat)
    history: List[float] = [norm_F]
    converged = False

    for it in range(max_iter):
        if verbose:
            gas_norm = float(np.sqrt(np.mean((F_hat[:n_gas]) ** 2)))
            sol_norm = float(np.sqrt(np.mean((F_hat[n_gas : n_gas + n_sol]) ** 2)))
            ene_norm = float(np.sqrt(np.mean((F_hat[n_gas + n_sol :]) ** 2)))
            print(
                f"  NR iter {it}: RMS={norm_F:.3e}  gas={gas_norm:.3e}  "
                f"solid={sol_norm:.3e}  energy={ene_norm:.3e}"
            )

        if norm_F < tol_rms:
            converged = True
            break

        # 构造 Jacobian。注意：这里 F 必须是与当前 x 对应的最新残差
        J = build_jacobian_fd(x, F, cells, apply_bc_fn, scale)

        try:
            # 使用稀疏解法求解 Newton 方向 Δx
            lu = splu(J.tocsc())
            dx = lu.solve(-F_hat)
        except Exception:
            # 退回到最小二乘法处理奇异矩阵
            Jd = J.toarray()
            dx, _, _, _ = np.linalg.lstsq(Jd, -F_hat, rcond=None)

        # 对 dx 进行物理约束下的步长剪切
        dx_clipped = _clip_dx(dx, cells, ref_gas_mol_s, ref_solid_kg_s)

        # --- 阻尼线搜索 (Damped Newton) ---
        lam = float(lambda_init)
        nt = norm_F
        found_step = False
        
        for damp_step in range(n_damp_halvings):
            x_trial = x + lam * dx_clipped
            
            # 对尝试步进行物理边界保护（必须正值，且 T 在范围内）
            offset = 0
            for cell in cells:
                nv = n_var(cell)
                # 气/固质量流率 >= 0
                x_trial[offset : offset + nv - 1] = np.maximum(x_trial[offset : offset + nv - 1], 0.0)
                # 温度保护
                x_trial[offset + nv - 1] = float(np.clip(x_trial[offset + nv - 1], 300.0, 2500.0))
                offset += nv

            # 计算尝试步的残差范数
            F_trial = global_residual(x_trial, cells, apply_bc_fn)
            F_hat_trial = F_trial / scale
            nt = _rms_norm(F_hat_trial)
            
            if nt < norm_F:
                if verbose and damp_step > 0:
                    print(f"    Line search OK at step {damp_step}: lam={lam:.4f}, norm={nt:.3e}")
                found_step = True
                x = x_trial
                F = F_trial
                F_hat = F_hat_trial
                norm_F = nt
                break
            
            lam *= 0.5
            if lam <= lambda_min:
                break

        if not found_step:
            if verbose:
                print(f"  NR iter {it}: Line search failed to reduce norm. Stopping.")
            break

        history.append(norm_F)

    # 结果解包回反应器对象
    unpack_reactor(x, cells)
    apply_bc_fn()
    F_final = global_residual(x, cells, apply_bc_fn)
    final_norm = float(np.linalg.norm(F_final))

    return {
        "converged": converged,
        "n_iter": len(history),
        "residual": final_norm,
        "rms_scaled_final": norm_F,
        "norm_history": history,
    }
