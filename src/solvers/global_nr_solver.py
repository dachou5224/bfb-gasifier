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
from time import perf_counter
from typing import List

import numpy as np
import scipy.sparse as sp

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import N_GAS
from src.solvers.convergence import InnerConvergence
from src.solvers.nr_indexing import (
    MAIN_NR_GAS_COUNT,
    MAIN_NR_GAS_IDX,
    affected_residual_cells_for_var_cell,
    build_var_index_to_cell_map,
    cell_offsets,
    n_solid_var,
    n_var,
    pack_cell,
    pack_reactor,
    solid_flat_indices_for_nr,
    unpack_cell,
    unpack_reactor,
)
from src.solvers.structured_jacobian import (
    assemble_jacobian_fd_structured,
    build_jacobian_structure,
    build_solver_graph,
    solve_sparse_direct_linear_step,
    solve_structured_linear_step,
    structured_strategy_name_from_cells,
)

_FD_EPS = 1e-6

# 本模块内历史命名（与 nr_indexing 常量一致）
_MAIN_NR_GAS_COUNT = MAIN_NR_GAS_COUNT
_MAIN_NR_GAS_IDX = MAIN_NR_GAS_IDX
_build_var_index_to_cell_map = build_var_index_to_cell_map
_affected_residual_cells_for_var_cell = affected_residual_cells_for_var_cell


def global_residual(
    x: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
) -> np.ndarray:
    """F(x)：解包 → 边界条件（含循环）→ 各 cell 残差拼接。"""
    unpack_reactor(x, cells)
    apply_bc_fn()
    return np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])


def _project_cell_residual_to_main_nr(cell: Cell) -> np.ndarray:
    """Project full cell residual vector to main-NR unknown subset (exclude minor gas species)."""
    full = cell.residuals()
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    gas_d = full[:N_GAS]
    gas_b = full[N_GAS : 2 * N_GAS]
    solid = full[2 * N_GAS : 2 * N_GAS + nv_sol][solid_flat_indices_for_nr(cell)]
    energy = full[-1]
    return np.concatenate(
        (
            gas_d[_MAIN_NR_GAS_IDX],
            gas_b[_MAIN_NR_GAS_IDX],
            solid,
            np.array([energy], dtype=np.float64),
        )
    )


def _solid_reference_for_cell(cell: Cell, ref_solid_kg_s: float) -> float:
    """Cell-aware solid reference scale for residual normalization and step clipping."""
    if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
        local_ref = float(
            np.sum(
                np.maximum(
                    cell.m_solid
                    + cell.m_solid_zu
                    + cell.m_solid_rez
                    + cell.m_solid_in
                    + cell.m_solid_auf_in
                    + cell.m_solid_ab_in,
                    0.0,
                )
            )
        )
        return max(local_ref, 1e-8)
    return max(float(ref_solid_kg_s), 1e-12)


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
    re = max(float(ref_energy_W), 1.0)

    n_total = sum(n_var(c) for c in cells)
    scale = np.empty(n_total, dtype=np.float64)
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        nv_sol = n_solid_var(cell)
        rs = _solid_reference_for_cell(cell, ref_solid_kg_s)
        scale[offset : offset + _MAIN_NR_GAS_COUNT] = rg
        scale[offset + _MAIN_NR_GAS_COUNT : offset + 2 * _MAIN_NR_GAS_COUNT] = rg
        scale[offset + 2 * _MAIN_NR_GAS_COUNT : offset + 2 * _MAIN_NR_GAS_COUNT + nv_sol] = rs
        scale[offset + 2 * _MAIN_NR_GAS_COUNT + nv_sol] = re
        offset += nv
    return scale


def build_jacobian_fd_dense(
    x: np.ndarray,
    F0: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    eq_scale: np.ndarray,
    verbose: bool = False,
) -> tuple[sp.csr_matrix, dict]:
    """列有限差分构造稀疏 Jacobian（行与 F̂=F/eq_scale 一致）。"""
    n_total = len(x)
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    residual_cell_calls = 0

    F0s = F0 / eq_scale

    if verbose and len(x) > 100:
        print(f"DEBUG: Building Jacobian for {n_total} variables...")

    x_p = x.copy()
    for j in range(n_total):
        if verbose and j > 0 and j % 100 == 0:
            print(f"DEBUG:   Jacobian progress: {j}/{n_total}")
        # FD 步长：相对步长 + 绝对下限，避免 x≈0 时 h 过小导致差分噪声
        xj = float(x[j])
        h = max(_FD_EPS * max(abs(xj), 1.0), 1e-8)
        x_p[j] = xj
        x_p[j] += h
        
        # 只在扰动后更新状态并计算残差
        unpack_reactor(x_p, cells)
        apply_bc_fn()
        F_p = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
        residual_cell_calls += len(cells)
        
        dF = (F_p / eq_scale - F0s) / h
        
        # 快速寻找非零元并记录（稀疏矩阵构造）
        indices = np.where(np.abs(dF) > 1e-14)[0]
        for i in indices:
            rows.append(int(i))
            cols.append(j)
            vals.append(float(dF[i]))
        x_p[j] = xj

    # 计算结束后恢复原始状态
    unpack_reactor(x, cells)
    apply_bc_fn()
    J = sp.csr_matrix((vals, (rows, cols)), shape=(n_total, n_total))
    zero_cols = int(np.count_nonzero(np.diff(J.tocsc().indptr) == 0))
    zero_rows = int(np.count_nonzero(np.diff(J.indptr) == 0))
    return J, {
        "strategy": "dense_fd",
        "residual_cell_calls": int(residual_cell_calls),
        "nnz": int(J.nnz),
        "zero_cols": zero_cols,
        "zero_rows": zero_rows,
    }


def build_jacobian_fd_block_tridiag(
    x: np.ndarray,
    F0: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    eq_scale: np.ndarray,
    apply_local_bc_fn: Callable[[int], None] | None = None,
    affected_residual_cells_fn: Callable[[int], tuple[int, ...]] | None = None,
    verbose: bool = False,
) -> tuple[sp.csr_matrix, dict]:
    """按 block-tridiagonal 依赖结构仅重算受影响 residual 块。"""
    n_total = len(x)
    n_cells = len(cells)
    offsets = cell_offsets(cells)
    var_to_cell = _build_var_index_to_cell_map(offsets)
    affected_cache = {
        cell_idx: (
            affected_residual_cells_fn(cell_idx)
            if affected_residual_cells_fn is not None
            else _affected_residual_cells_for_var_cell(cell_idx, n_cells)
        )
        for cell_idx in range(n_cells)
    }

    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    residual_cell_calls = 0

    F0s = F0 / eq_scale

    if verbose and len(x) > 100:
        print(f"DEBUG: Building block-tridiagonal Jacobian for {n_total} variables...")

    x_p = x.copy()
    for j in range(n_total):
        if verbose and j > 0 and j % 100 == 0:
            print(f"DEBUG:   Jacobian progress: {j}/{n_total}")

        xj = float(x[j])
        h = max(_FD_EPS * max(abs(xj), 1.0), 1e-8)
        x_p[j] = xj
        x_p[j] += h

        unpack_reactor(x_p, cells)
        if apply_local_bc_fn is not None:
            apply_local_bc_fn(int(var_to_cell[j]))
        else:
            apply_bc_fn()

        for affected_cell_idx in affected_cache[int(var_to_cell[j])]:
            row_start = offsets[affected_cell_idx]
            row_end = offsets[affected_cell_idx + 1]
            F_p_local = _project_cell_residual_to_main_nr(cells[affected_cell_idx])
            residual_cell_calls += 1
            dF_local = (F_p_local / eq_scale[row_start:row_end] - F0s[row_start:row_end]) / h
            indices = np.flatnonzero(np.abs(dF_local) > 1e-14)
            if indices.size == 0:
                continue
            rows.extend((row_start + indices).tolist())
            cols.extend([j] * int(indices.size))
            vals.extend(dF_local[indices].astype(float).tolist())
        x_p[j] = xj

    unpack_reactor(x, cells)
    apply_bc_fn()
    J = sp.csr_matrix((vals, (rows, cols)), shape=(n_total, n_total))
    zero_cols = int(np.count_nonzero(np.diff(J.tocsc().indptr) == 0))
    zero_rows = int(np.count_nonzero(np.diff(J.indptr) == 0))
    return J, {
        "strategy": "block_tridiag_fd",
        "residual_cell_calls": int(residual_cell_calls),
        "nnz": int(J.nnz),
        "zero_cols": zero_cols,
        "zero_rows": zero_rows,
    }


def build_jacobian_fd(
    x: np.ndarray,
    F0: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    eq_scale: np.ndarray,
    *,
    strategy: str = "block_tridiag_structured",
    apply_local_bc_fn: Callable[[int], None] | None = None,
    affected_residual_cells_fn: Callable[[int], tuple[int, ...]] | None = None,
    verbose: bool = False,
) -> tuple[sp.csr_matrix, dict]:
    if strategy == "dense_fd":
        return build_jacobian_fd_dense(
            x=x,
            F0=F0,
            cells=cells,
            apply_bc_fn=apply_bc_fn,
            eq_scale=eq_scale,
            verbose=verbose,
        )
    if strategy in {"block_tridiag_structured", "band_plus_side_elements_structured", "side_elements_sparse_fd"}:
        resolved_strategy = (
            "band_plus_side_elements_structured"
            if strategy == "side_elements_sparse_fd"
            else str(strategy)
        )
        graph = build_solver_graph(cells)
        structure = build_jacobian_structure(graph, cells)
        structured = assemble_jacobian_fd_structured(
            x=x,
            F0=F0,
            cells=cells,
            apply_bc_fn=apply_bc_fn,
            eq_scale=eq_scale,
            structure=structure,
            verbose=verbose,
        )
        return structured.matrix, {
            "strategy": resolved_strategy,
            "residual_cell_calls": int(structured.residual_cell_calls),
            "nnz": int(structured.matrix.nnz),
            "zero_cols": int(structured.zero_cols),
            "zero_rows": int(structured.zero_rows),
            "structured_jacobian": structured,
            "jacobian_structure": {
                "main_chain_blocks": int(len(structure.main_chain)),
                "band_block_count": int(len(structure.band_block_pairs)),
                "side_element_count": int(len(structured.side_blocks)),
                "side_tail_block_count": int(len(structure.side_tail_blocks)),
                "side_head_block_count": int(len(structure.side_head_blocks)),
                "unexpected_block_count": int(len(structured.unexpected_blocks)),
                "unexpected_blocks": [list(p) for p in structured.unexpected_blocks],
                "structure_validation_ok": bool(structured.structure_validation_ok),
            },
        }
    if strategy == "block_tridiag_fd":
        return build_jacobian_fd_block_tridiag(
            x=x,
            F0=F0,
            cells=cells,
            apply_bc_fn=apply_bc_fn,
            eq_scale=eq_scale,
            apply_local_bc_fn=apply_local_bc_fn,
            affected_residual_cells_fn=affected_residual_cells_fn,
            verbose=verbose,
        )
    raise ValueError(f"Unsupported Jacobian strategy: {strategy!r}")


def _rms_norm(F_hat: np.ndarray) -> float:
    """缩放残差向量的 RMS（用于收敛判据）。"""
    if len(F_hat) == 0:
        return 0.0
    return float(np.sqrt(np.mean(F_hat**2)))


def _state_rms_norm(x_vec: np.ndarray) -> float:
    """全局状态向量 RMS，用于 InnerConvergence 中 rtol·‖x‖ 项（Check1）。"""
    if len(x_vec) == 0:
        return 0.0
    return float(np.sqrt(np.mean(x_vec**2)))


def _clip_dx(
    dx: np.ndarray,
    cells: List[Cell],
    ref_gas_mol_s: float,
    ref_solid_kg_s: float,
) -> np.ndarray:
    """严格限制 Newton 步长，防止物理量跳变导致发散。"""
    rg_limit = 0.2 * max(float(ref_gas_mol_s), 1.0)
    dT_limit = 35.0  # Historical damping guardrail: 35 K keeps LU convergence behavior stable
                      # while reducing outer-loop temperature stair-stepping vs tighter caps.

    out = dx.copy()
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        nv_sol = n_solid_var(cell)
        rs_limit = 0.2 * max(_solid_reference_for_cell(cell, ref_solid_kg_s), 0.1)
        # 气相
        out[offset : offset + 2 * _MAIN_NR_GAS_COUNT] = np.clip(
            out[offset : offset + 2 * _MAIN_NR_GAS_COUNT], -rg_limit, rg_limit
        )
        # 固相
        out[offset + 2 * _MAIN_NR_GAS_COUNT : offset + 2 * _MAIN_NR_GAS_COUNT + nv_sol] = np.clip(
            out[offset + 2 * _MAIN_NR_GAS_COUNT : offset + 2 * _MAIN_NR_GAS_COUNT + nv_sol], -rs_limit, rs_limit
        )
        # 温度
        out[offset + nv - 1] = float(np.clip(out[offset + nv - 1], -dT_limit, dT_limit))
        offset += nv
    return out


def _clip_dx_diagnostics(dx: np.ndarray, dx_clipped: np.ndarray, cells: List[Cell]) -> dict:
    """Summarize how strongly the Newton step was clipped this iteration."""
    clipped_mask = np.abs(dx_clipped - dx) > 1e-12
    temp_clipped_cells = 0
    max_abs_dT_step = 0.0
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        temp_idx = offset + nv - 1
        if bool(clipped_mask[temp_idx]):
            temp_clipped_cells += 1
        max_abs_dT_step = max(max_abs_dT_step, abs(float(dx_clipped[temp_idx])))
        offset += nv
    return {
        "clipped_fraction": float(np.mean(clipped_mask)) if clipped_mask.size else 0.0,
        "clipped_components": int(np.count_nonzero(clipped_mask)),
        "temp_clipped_cells": int(temp_clipped_cells),
        "max_abs_dT_step": float(max_abs_dT_step),
    }


def _next_lambda_seed(lambda_cap: float, accepted_lambda: float, trial_evals: int) -> float:
    """Adapt the next line-search start from the last accepted damping factor."""
    lam_cap = float(max(lambda_cap, 1e-12))
    lam_acc = float(np.clip(accepted_lambda, 1e-12, lam_cap))
    if int(trial_evals) <= 1:
        return float(min(lam_cap, 2.0 * lam_acc))
    return lam_acc


def _solve_linear_step(
    J: sp.csr_matrix,
    rhs: np.ndarray,
    *,
    zero_rows: int = 0,
    zero_cols: int = 0,
    structured_jacobian=None,
    linear_solver_backend: str = "sparse_direct_fallback",
) -> tuple[np.ndarray, str, bool, dict[str, float | int | bool]]:
    """Solve linear Newton step with optional diagonal regularization for rank-deficient Jacobian."""
    use_regularization = int(zero_rows) > 0 or int(zero_cols) > 0
    if linear_solver_backend == "structured_direct" and structured_jacobian is not None:
        dx, method, regularized, diag = solve_structured_linear_step(
            structured_jacobian,
            rhs,
            regularize=use_regularization,
        )
        return dx, method, regularized, diag
    dx, method, regularized = solve_sparse_direct_linear_step(
        J,
        rhs,
        regularize=use_regularization,
    )
    return dx, method, regularized, {
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "schur_size": 0,
        "fallback_used": method != "splu",
    }


def solve_global_nr(
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    apply_local_bc_fn: Callable[[int], None] | None = None,
    affected_residual_cells_fn: Callable[[int], tuple[int, ...]] | None = None,
    ref_gas_mol_s: float = 80.0,
    ref_solid_kg_s: float = 1.0,
    ref_energy_W: float = 2e7,
    max_iter: int = 25,
    tol_rms: float = 0.01,
    inner_convergence: InnerConvergence | None = None,
    lambda_init: float = 0.5,  # 0.1 过于保守；配合 _clip_dx 会把有效温度步长压到约 5 K/iter
    n_damp_halvings: int = 8,
    lambda_min: float = 1.0 / 1024.0,
    jacobian_strategy: str = "block_tridiag_structured",
    linear_solver_backend: str | None = None,
    jacobian_lag_steps: int = 1,
    verbose: bool = False,
) -> dict:
    """阻尼 Newton–Raphson：ĴΔx = −F̂（F̂ = F / 静态 eq_scale），再线搜索阻尼。

    Hamel Eq. 2.9（gedämpftes Newton-Verfahren）的工程实现。

    inner_convergence
        若给定，Check1 使用 ``atol + rtol * ‖x‖_RMS``（‖x‖ 为打包状态向量 RMS）；
        若为 ``None``，则沿用 ``norm_F < tol_rms``（与历史行为一致）。
    """
    solve_started = perf_counter()
    x = pack_reactor(cells)
    timing = {
        "initial_residual_s": 0.0,
        "jacobian_build_s": 0.0,
        "linear_solve_s": 0.0,
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "line_search_residual_s": 0.0,
        "final_residual_s": 0.0,
        "total_s": 0.0,
    }
    counts = {
        "global_residual_calls": 0,
        "jacobian_cell_residual_calls": 0,
        "jacobian_nnz_last": 0,
        "jacobian_zero_cols_last": 0,
        "jacobian_zero_rows_last": 0,
        "jacobian_rebuilds": 0,
        "jacobian_reuses": 0,
        "line_search_evaluations": 0,
        "line_search_backtracks": 0,
        "line_search_failures": 0,
        "line_search_retries": 0,
        "linear_regularized_solves": 0,
        "linear_lstsq_fallbacks": 0,
        "structured_fallbacks": 0,
        "nr_schur_size_last": 0,
    }
    accepted_lambda_history: list[float | None] = []
    line_search_trial_counts: list[int] = []
    clip_history: list[dict] = []

    t0 = perf_counter()
    F = global_residual(x, cells, apply_bc_fn)
    timing["initial_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] += 1
    x = pack_reactor(cells)
    scale = build_equation_scales(cells, ref_gas_mol_s, ref_solid_kg_s, ref_energy_W)

    n_total = len(F)
    n_gas = 2 * _MAIN_NR_GAS_COUNT * len(cells)
    n_sol = sum(n_solid_var(c) for c in cells)
    n_ene = len(cells)

    F_hat = F / scale
    norm_F = _rms_norm(F_hat)
    history: List[float] = [norm_F]
    converged = False
    lambda_seed = float(lambda_init)
    consecutive_line_search_failures = 0
    iter_attempts = 0
    lag_steps = max(int(jacobian_lag_steps), 1)
    cached_J: sp.csr_matrix | None = None
    cached_jac_meta: dict | None = None
    cached_jac_iter = -1
    linear_solver_backend_last = ""
    if linear_solver_backend is None:
        linear_solver_backend = (
            "structured_direct"
            if str(jacobian_strategy) in {"block_tridiag_structured", "band_plus_side_elements_structured"}
            else "sparse_direct_fallback"
        )

    for it in range(max_iter):
        if verbose:
            gas_norm = float(np.sqrt(np.mean((F_hat[:n_gas]) ** 2)))
            sol_norm = float(np.sqrt(np.mean((F_hat[n_gas : n_gas + n_sol]) ** 2)))
            ene_norm = float(np.sqrt(np.mean((F_hat[n_gas + n_sol :]) ** 2)))
            print(
                f"  NR iter {it}: RMS={norm_F:.3e}  gas={gas_norm:.3e}  "
                f"solid={sol_norm:.3e}  energy={ene_norm:.3e}"
            )

        x_rms = _state_rms_norm(x)
        if inner_convergence is not None:
            if inner_convergence.is_converged(norm_F, x_rms):
                converged = True
                break
        elif norm_F < tol_rms:
            converged = True
            break
        iter_attempts += 1

        # 构造 Jacobian。注意：这里 F 必须是与当前 x 对应的最新残差
        # 平稳阶段自适应 Jacobian lag：
        # 当上一轮 line-search 一次命中且阻尼不小，允许多复用一次 Jacobian
        # 以减少 FD 构造开销（不改变物理方程，仅改变线性化频率）。
        smooth_step = (
            len(line_search_trial_counts) > 0
            and line_search_trial_counts[-1] == 1
            and len(accepted_lambda_history) > 0
            and accepted_lambda_history[-1] is not None
            and float(accepted_lambda_history[-1]) >= 0.5
            and norm_F < max(10.0 * tol_rms, 0.2)
        )
        lag_steps_eff = max(lag_steps, 2) if smooth_step else lag_steps

        rebuild_jacobian = (
            cached_J is None
            or (it - cached_jac_iter) >= lag_steps_eff
        )
        if rebuild_jacobian:
            t0 = perf_counter()
            J, jac_meta = build_jacobian_fd(
                x,
                F,
                cells,
                apply_bc_fn,
                scale,
                strategy=jacobian_strategy,
                apply_local_bc_fn=apply_local_bc_fn,
                affected_residual_cells_fn=affected_residual_cells_fn,
                verbose=verbose,
            )
            timing["jacobian_build_s"] += perf_counter() - t0
            cached_J = J
            cached_jac_meta = jac_meta
            cached_jac_iter = it
            counts["jacobian_rebuilds"] += 1
            counts["jacobian_cell_residual_calls"] += int(jac_meta.get("residual_cell_calls", 0))
            counts["jacobian_nnz_last"] = int(jac_meta.get("nnz", J.nnz))
            counts["jacobian_zero_cols_last"] = int(jac_meta.get("zero_cols", 0))
            counts["jacobian_zero_rows_last"] = int(jac_meta.get("zero_rows", 0))
        else:
            J = cached_J
            jac_meta = cached_jac_meta or {}
            counts["jacobian_reuses"] += 1

        # 使用稀疏解法求解 Newton 方向 Δx；对秩亏 Jacobian 做轻度对角正则化。
        t0 = perf_counter()
        dx, linear_solver_used, linear_regularized, linear_diag = _solve_linear_step(
            J,
            -F_hat,
            zero_rows=int(jac_meta.get("zero_rows", 0)),
            zero_cols=int(jac_meta.get("zero_cols", 0)),
            structured_jacobian=jac_meta.get("structured_jacobian"),
            linear_solver_backend=str(linear_solver_backend),
        )
        timing["linear_solve_s"] += perf_counter() - t0
        timing["band_lu_s"] = timing.get("band_lu_s", 0.0) + float(linear_diag.get("band_lu_s", 0.0))
        timing["side_update_s"] = timing.get("side_update_s", 0.0) + float(linear_diag.get("side_update_s", 0.0))
        if linear_regularized:
            counts["linear_regularized_solves"] += 1
        if linear_solver_used == "lstsq":
            counts["linear_lstsq_fallbacks"] += 1
        linear_solver_backend_last = str(linear_solver_used)
        counts["structured_fallbacks"] = int(counts.get("structured_fallbacks", 0)) + int(
            bool(linear_diag.get("fallback_used", False) and str(linear_solver_backend) == "structured_direct")
        )
        counts["nr_schur_size_last"] = int(linear_diag.get("schur_size", 0))

        # 对 dx 进行物理约束下的步长剪切
        dx_clipped = _clip_dx(dx, cells, ref_gas_mol_s, ref_solid_kg_s)
        clip_diag = _clip_dx_diagnostics(dx, dx_clipped, cells)

        # --- 阻尼线搜索 (Damped Newton) ---
        lam = float(np.clip(lambda_seed, lambda_min, max(float(lambda_init), lambda_min)))
        nt = norm_F
        found_step = False
        trial_evals = 0
        best_trial_nt = float(np.inf)
        best_trial_lambda: float | None = None
        
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
            t0 = perf_counter()
            F_trial = global_residual(x_trial, cells, apply_bc_fn)
            timing["line_search_residual_s"] += perf_counter() - t0
            counts["global_residual_calls"] += 1
            trial_evals += 1
            F_hat_trial = F_trial / scale
            nt = _rms_norm(F_hat_trial)
            if np.isfinite(nt) and nt < best_trial_nt:
                best_trial_nt = float(nt)
                best_trial_lambda = float(lam)
            
            if nt < norm_F:
                if verbose and damp_step > 0:
                    print(f"    Line search OK at step {damp_step}: lam={lam:.4f}, norm={nt:.3e}")
                found_step = True
                consecutive_line_search_failures = 0
                x = pack_reactor(cells)
                F = F_trial
                F_hat = F_hat_trial
                norm_F = nt
                counts["line_search_evaluations"] += int(trial_evals)
                counts["line_search_backtracks"] += max(int(trial_evals) - 1, 0)
                accepted_lambda_history.append(float(lam))
                line_search_trial_counts.append(int(trial_evals))
                clip_history.append(
                    {
                        "iter": int(it + 1),
                        "accepted_lambda": float(lam),
                        "line_search_trials": int(trial_evals),
                        "line_search_failed": False,
                        "best_trial_rms_scaled": float(best_trial_nt),
                        "best_trial_lambda": best_trial_lambda,
                        **clip_diag,
                    }
                )
                lambda_seed = _next_lambda_seed(float(lambda_init), float(lam), int(trial_evals))
                break
            
            lam *= 0.5
            if lam <= lambda_min:
                break

        if not found_step:
            counts["line_search_evaluations"] += int(trial_evals)
            counts["line_search_backtracks"] += max(int(trial_evals) - 1, 0)
            counts["line_search_failures"] += 1
            accepted_lambda_history.append(None)
            line_search_trial_counts.append(int(trial_evals))
            clip_history.append(
                {
                    "iter": int(it + 1),
                    "accepted_lambda": None,
                    "line_search_trials": int(trial_evals),
                    "line_search_failed": True,
                    "best_trial_rms_scaled": float(best_trial_nt),
                    "best_trial_lambda": best_trial_lambda,
                    **clip_diag,
                }
            )
            if verbose:
                print(f"  NR iter {it}: Line search failed to reduce norm. Stopping.")
            consecutive_line_search_failures += 1
            retry_allowed = consecutive_line_search_failures <= 1 and (it + 1) < max_iter
            if retry_allowed:
                counts["line_search_retries"] += 1
                lam_retry_base = float(best_trial_lambda) if best_trial_lambda is not None else float(lambda_seed)
                lambda_seed = float(np.clip(0.5 * lam_retry_base, lambda_min, max(float(lambda_init), lambda_min)))
                cached_J = None
                cached_jac_meta = None
                cached_jac_iter = -1
                continue
            break

        history.append(norm_F)

    # 结果解包回反应器对象
    t0 = perf_counter()
    F_final = global_residual(x, cells, apply_bc_fn)
    timing["final_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] += 1
    x = pack_reactor(cells)
    final_norm = float(np.linalg.norm(F_final))
    timing["total_s"] = perf_counter() - solve_started

    return {
        "converged": converged,
        "n_iter": len(history),
        "residual": final_norm,
        "rms_scaled_final": norm_F,
        "norm_history": history,
        "n_newton_iters_attempted": int(iter_attempts),
        "jacobian_strategy": jacobian_strategy,
        "jacobian_structure": (
            cached_jac_meta.get("jacobian_structure")
            if cached_jac_meta is not None
            else {
                "main_chain_blocks": int(len(cells)),
                "band_block_count": 0,
                "side_element_count": 0,
                "side_tail_block_count": 0,
                "side_head_block_count": 0,
                "unexpected_block_count": 0,
                "unexpected_blocks": [],
                "structure_validation_ok": True,
            }
        ),
        "linear_solver_backend": str(linear_solver_backend),
        "linear_solver_backend_last": linear_solver_backend_last,
        "jacobian_lag_steps": int(lag_steps),
        "timing": timing,
        "counts": counts,
        "accepted_lambda_history": accepted_lambda_history,
        "line_search_trial_counts": line_search_trial_counts,
        "clip_history": clip_history,
    }
