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

from src.core.cell import Cell
from src.core.species import N_GAS

_FD_EPS = np.sqrt(np.finfo(float).eps)


def n_var(cell: Cell) -> int:
    """每格未知量：N_d + N_b + m_solid + T。"""
    return 2 * N_GAS + cell.solid.n_size_classes + 1


def pack_cell(cell: Cell) -> np.ndarray:
    return np.concatenate([cell.N_d, cell.N_b, cell.m_solid, np.array([cell.T], dtype=np.float64)])


def unpack_cell(x: np.ndarray, cell: Cell) -> None:
    nk = cell.solid.n_size_classes
    np.maximum(x[:N_GAS], 0.0, out=cell.N_d)
    np.maximum(x[N_GAS : 2 * N_GAS], 0.0, out=cell.N_b)
    np.maximum(x[2 * N_GAS : 2 * N_GAS + nk], 0.0, out=cell.m_solid)
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
        scale[offset : offset + N_GAS] = rg
        scale[offset + N_GAS : offset + 2 * N_GAS] = rg
        scale[offset + 2 * N_GAS : offset + 2 * N_GAS + nk] = rs
        scale[offset + 2 * N_GAS + nk] = re
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

    x_base = x.copy()
    F0s = F0 / eq_scale

    for j in range(n_total):
        h = _FD_EPS * max(1.0, abs(x_base[j]))
        x_p = x_base.copy()
        x_p[j] += h
        F_p = global_residual(x_p, cells, apply_bc_fn)
        dF = (F_p / eq_scale - F0s) / h
        unpack_reactor(x_base, cells)
        apply_bc_fn()
        for i in range(n_total):
            v = float(dF[i])
            if abs(v) > 1e-14:
                rows.append(i)
                cols.append(j)
                vals.append(v)

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
    """用绝对参考界限制 Newton 步，避免“小流量相对界”导致无法爬升。"""
    rg = max(float(ref_gas_mol_s), 1e-9)
    rs = max(float(ref_solid_kg_s), 1e-12)
    out = dx.copy()
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        nk = cell.solid.n_size_classes
        out[offset : offset + 2 * N_GAS] = np.clip(
            out[offset : offset + 2 * N_GAS],
            -rg,
            rg,
        )
        out[offset + 2 * N_GAS : offset + 2 * N_GAS + nk] = np.clip(
            out[offset + 2 * N_GAS : offset + 2 * N_GAS + nk],
            -rs,
            rs,
        )
        out[offset + nv - 1] = float(np.clip(out[offset + nv - 1], -200.0, 200.0))
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
    lambda_init: float = 0.7,
    n_damp_halvings: int = 8,
    lambda_min: float = 1.0 / 256.0,
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
    n_sol = sum(c.solid.n_size_classes for c in cells)
    n_ene = len(cells)
    assert n_gas + n_sol + n_ene == n_total, "残差维数与方程类型分块不一致"

    F_hat = F / scale
    norm0 = _rms_norm(F_hat)
    history: List[float] = [norm0]
    converged = False

    for it in range(max_iter):
        F_hat = F / scale
        norm_F = _rms_norm(F_hat)
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

        J = build_jacobian_fd(x, F, cells, apply_bc_fn, scale)

        try:
            lu = splu(J.tocsc())
            dx = lu.solve(-F_hat)
        except Exception:
            Jd = J.toarray()
            try:
                dx = np.linalg.solve(Jd, -F_hat)
            except np.linalg.LinAlgError:
                dx, _, _, _ = np.linalg.lstsq(Jd, -F_hat, rcond=None)

        dx = _clip_dx(dx, cells, ref_gas_mol_s, ref_solid_kg_s)

        lam = float(lambda_init)
        x_trial = x.copy()
        F_trial = F.copy()
        nt = norm_F

        for _ in range(n_damp_halvings):
            x_trial = x + lam * dx
            offset = 0
            for cell in cells:
                nv = n_var(cell)
                sl = slice(offset, offset + nv - 1)
                x_trial[sl] = np.maximum(x_trial[sl], 0.0)
                x_trial[offset + nv - 1] = float(
                    np.clip(x_trial[offset + nv - 1], 300.0, 2500.0)
                )
                offset += nv

            F_trial = global_residual(x_trial, cells, apply_bc_fn)
            F_hat_trial = F_trial / scale
            nt = _rms_norm(F_hat_trial)
            if nt < norm_F or lam <= lambda_min:
                break
            lam *= 0.5

        x = x_trial.copy()
        F = F_trial.copy()
        history.append(nt)

    unpack_reactor(x, cells)
    apply_bc_fn()
    F_final = global_residual(x, cells, apply_bc_fn)
    final_norm = float(np.linalg.norm(F_final))
    rms_scaled_final = _rms_norm(F_final / scale)

    return {
        "converged": converged,
        "n_iter": len(history),
        "residual": final_norm,
        "rms_scaled_final": rms_scaled_final,
        "norm_history": history,
    }
