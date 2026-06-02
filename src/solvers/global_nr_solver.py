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
    gas_mode_for_nr,
    has_temperature_var,
    n_gas_var,
    n_solid_var,
    n_var,
    pack_cell,
    pack_reactor,
    solid_offset,
    solid_flat_indices_for_nr,
    temperature_offset,
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


def _is_bed_two_phase_gas_cell(cell: Cell) -> bool:
    """Cells whose NR layout keeps distinct dense and bubble gas balances."""
    return gas_mode_for_nr(cell) == "two_phase"


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
    gas_mode = gas_mode_for_nr(cell)
    parts: list[np.ndarray] = []
    if gas_mode == "two_phase":
        parts.extend([gas_d[_MAIN_NR_GAS_IDX], gas_b[_MAIN_NR_GAS_IDX]])
    elif gas_mode == "single":
        parts.append((gas_d + gas_b)[_MAIN_NR_GAS_IDX])
    if solid.size:
        parts.append(solid)
    if has_temperature_var(cell):
        parts.append(np.array([energy], dtype=np.float64))
    out = np.concatenate(parts)
    expected = n_var(cell)
    if out.size != expected:
        raise AssertionError(
            f"NR residual layout mismatch for {getattr(cell, 'cell_type', 'bed')}: "
            f"projected={out.size}, n_var={expected}"
        )
    return out


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
                + np.maximum(cell.R_solid, 0.0)
            )
        )
        # Use at least 10% of feed rate as minimum scale to prevent 1e-8 collapse
        # in cells with near-zero holdup (e.g., upper bed cells at initialization).
        # A 1e-8 minimum causes F_hat explosion when any solid residual perturbs
        # away from zero during line search (even tiny changes → F_hat ~ 1e8).
        return max(local_ref, max(float(ref_solid_kg_s) * 0.1, 1e-6))
    return max(float(ref_solid_kg_s), 1e-12)


def _solid_equation_scales_for_cell(cell: Cell, ref_solid_kg_s: float) -> np.ndarray:
    """Per-solid-DOF residual scales for the packed NR solid balances."""
    nv_sol = n_solid_var(cell)
    if nv_sol <= 0:
        return np.zeros(0, dtype=np.float64)
    floor = max(float(ref_solid_kg_s) * 1e-3, 1e-8)
    if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
        support = (
            np.maximum(
                cell.m_solid
                + cell.m_solid_zu
                + cell.m_solid_rez
                + cell.m_solid_in
                + cell.m_solid_auf_in
                + cell.m_solid_ab_in,
                0.0,
            )
            + np.abs(cell.R_solid)
        )
        return np.maximum(support.reshape(-1)[solid_flat_indices_for_nr(cell)], floor)
    return np.full(nv_sol, max(float(ref_solid_kg_s), 1e-12), dtype=np.float64)


def _is_solid_empty_cell(cell: Cell, ref_solid_kg_s: float = 1.0) -> bool:
    """Return True if a holdup_transport cell has effectively zero solid inventory AND
    zero solid inflow from ALL sources.  Perturbing such cells creates artificial
    outflow with no inflow → residual spike.  Used to pin solid dx to 0 in those cells.

    Threshold is set relative to the total feed rate so that cells that received a
    tiny amount of solid from the first accepted NR step (e.g. 1e-4 kg/s when feed
    is 1 kg/s) are still treated as "empty" and not perturbed further.
    """
    if str(getattr(cell, "solid_state_model", "legacy_stream")) not in {"holdup_transport", "freeboard_closure"}:
        return False
    raw = float(
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
                + np.maximum(cell.R_solid, 0.0)
            )
        )
    # Use 0.1% of the solid feed rate (or at least 1e-5) as the "empty" threshold.
    threshold = max(float(ref_solid_kg_s) * 1e-3, 1e-5)
    return raw < threshold


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
        ng = n_gas_var(cell)
        nv_sol = n_solid_var(cell)
        rs = _solid_reference_for_cell(cell, ref_solid_kg_s)
        if ng > 0:
            scale[offset : offset + ng] = rg
        sol_start = offset + solid_offset(cell)
        scale[sol_start : sol_start + nv_sol] = (
            _solid_equation_scales_for_cell(cell, ref_solid_kg_s)
            if nv_sol > 0
            else rs
        )
        if has_temperature_var(cell):
            scale[offset + temperature_offset(cell)] = re
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

    unpack_reactor(x, cells)
    apply_bc_fn()
    F0 = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
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

    unpack_reactor(x, cells)
    apply_bc_fn()
    F0 = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
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

        affected_cells = tuple(affected_cache[int(var_to_cell[j])])
        if len(affected_cells) == n_cells:
            F_p_all = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
            residual_cell_calls += len(cells)
        else:
            F_p_all = None
        for affected_cell_idx in affected_cells:
            row_start = offsets[affected_cell_idx]
            row_end = offsets[affected_cell_idx + 1]
            if F_p_all is None:
                F_p_local = _project_cell_residual_to_main_nr(cells[affected_cell_idx])
                residual_cell_calls += 1
            else:
                F_p_local = F_p_all[row_start:row_end]
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
        use_local_callbacks = resolved_strategy != "band_plus_side_elements_structured"
        # Band+side structured Jacobian is recycle-coupled and sensitive to nonlocal
        # boundary propagation. Local-BC/affected-row pruning can miss these terms
        # (both closure-owned and explicit-freeboard graphs), so always assemble with
        # full-BC callbacks disabled for this strategy.
        structured = assemble_jacobian_fd_structured(
            x=x,
            F0=F0,
            cells=cells,
            apply_bc_fn=apply_bc_fn,
            eq_scale=eq_scale,
            structure=structure,
            apply_local_bc_fn=(apply_local_bc_fn if use_local_callbacks else None),
            affected_residual_cells_fn=(affected_residual_cells_fn if use_local_callbacks else None),
            verbose=verbose,
        )
        return structured.matrix, {
            "strategy": resolved_strategy,
            "residual_cell_calls": int(structured.residual_cell_calls),
            "nnz": int(structured.matrix.nnz),
            "zero_cols": int(structured.zero_cols),
            "zero_rows": int(structured.zero_rows),
            "structured_jacobian": structured,
            "fd_row_support_by_col": structured.row_support_by_col,
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


def _residual_group_metrics(F_hat: np.ndarray, cells: List[Cell]) -> dict[str, float]:
    """Compute gas/solid/energy RMS and max metrics from the real DOF layout."""
    gas_terms: list[np.ndarray] = []
    sol_terms: list[np.ndarray] = []
    ene_terms: list[float] = []
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        ng = n_gas_var(cell)
        ns = n_solid_var(cell)
        gas_terms.append(F_hat[offset : offset + ng])
        sol_start = offset + solid_offset(cell)
        sol_terms.append(F_hat[sol_start : sol_start + ns])
        ene_terms.append(float(F_hat[offset + nv - 1]))
        offset += nv
    gas = np.concatenate(gas_terms) if gas_terms else np.zeros(0, dtype=np.float64)
    solid = np.concatenate(sol_terms) if sol_terms else np.zeros(0, dtype=np.float64)
    energy = np.asarray(ene_terms, dtype=np.float64)
    max_abs = lambda arr: float(np.max(np.abs(arr))) if arr.size else 0.0
    return {
        "gas_rms": float(np.sqrt(np.mean(gas**2))) if gas.size else 0.0,
        "solid_rms": float(np.sqrt(np.mean(solid**2))) if solid.size else 0.0,
        "energy_rms": float(np.sqrt(np.mean(energy**2))) if energy.size else 0.0,
        "gas_max_abs": max_abs(gas),
        "solid_max_abs": max_abs(solid),
        "energy_max_abs": max_abs(energy),
    }


def _residual_gas_phase_metrics(F_hat: np.ndarray, cells: List[Cell]) -> dict[str, float]:
    """Separate total gas conservation from dense/bubble phase-split residuals."""
    combined_terms: list[np.ndarray] = []
    split_terms: list[np.ndarray] = []
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        gas_mode = gas_mode_for_nr(cell)
        if gas_mode == "two_phase":
            dense = F_hat[offset : offset + _MAIN_NR_GAS_COUNT]
            bubble = F_hat[offset + _MAIN_NR_GAS_COUNT : offset + 2 * _MAIN_NR_GAS_COUNT]
            if dense.size and bubble.size:
                combined_terms.append(dense + bubble)
                split_terms.append(dense - bubble)
        elif gas_mode == "single":
            gas = F_hat[offset : offset + _MAIN_NR_GAS_COUNT]
            if gas.size:
                combined_terms.append(gas)
        offset += nv
    combined = np.concatenate(combined_terms) if combined_terms else np.zeros(0, dtype=np.float64)
    split = np.concatenate(split_terms) if split_terms else np.zeros(0, dtype=np.float64)
    max_abs = lambda arr: float(np.max(np.abs(arr))) if arr.size else 0.0
    return {
        "gas_combined_rms": float(np.sqrt(np.mean(combined**2))) if combined.size else 0.0,
        "gas_phase_split_rms": float(np.sqrt(np.mean(split**2))) if split.size else 0.0,
        "gas_combined_max_abs": max_abs(combined),
        "gas_phase_split_max_abs": max_abs(split),
    }


def _residual_group_norms(F_hat: np.ndarray, cells: List[Cell]) -> dict[str, float]:
    """Compute gas/solid/energy RMS using the real cell-by-cell DOF layout."""
    metrics = _residual_group_metrics(F_hat, cells)
    return {
        "gas": metrics["gas_rms"],
        "solid": metrics["solid_rms"],
        "energy": metrics["energy_rms"],
    }


def _nr_layout_audit(cells: List[Cell], projected_residual: np.ndarray | None = None) -> list[dict[str, int | str]]:
    """Describe the per-cell NR row/column layout for result diagnostics."""
    out: list[dict[str, int | str]] = []
    offsets = cell_offsets(cells)
    for idx, cell in enumerate(cells):
        nv = n_var(cell)
        if projected_residual is None:
            residual_len = int(_project_cell_residual_to_main_nr(cell).size)
        else:
            residual_len = int(offsets[idx + 1] - offsets[idx])
        out.append(
            {
                "index": int(idx),
                "cell_type": str(getattr(cell, "cell_type", "bed")),
                "gas_mode": gas_mode_for_nr(cell),
                "gas_unknown_count": int(n_gas_var(cell)),
                "solid_unknown_count": int(n_solid_var(cell)),
                "unknown_count": int(nv),
                "residual_count": int(residual_len),
            }
        )
    return out


def _clip_dx(
    dx: np.ndarray,
    cells: List[Cell],
    ref_gas_mol_s: float,
    ref_solid_kg_s: float,
    zero_empty_solid: bool = False,
    t_step_limit_K: float | None = None,
) -> np.ndarray:
    """严格限制 Newton 步长，防止物理量跳变导致发散。

    对物理上"空"的 holdup_transport 单元（m_solid ≈ 0 且无任何固相进流），
    仅当 zero_empty_solid=True 时才将固相 dx 清零（GD fallback 专用）。
    NR 方向不应清零：Newton 步由 Jacobian 计算，天然保证 descent，
    强制清零会使最大残差 DOF 无法更新，导致线搜索反下降。
    """
    rg_limit = 0.2 * max(float(ref_gas_mol_s), 1.0)

    out = dx.copy()
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        ng = n_gas_var(cell)
        nv_sol = n_solid_var(cell)
        rs_limit = 0.2 * max(_solid_reference_for_cell(cell, ref_solid_kg_s), 0.1)
        # 气相
        if ng > 0:
            out[offset : offset + ng] = np.clip(out[offset : offset + ng], -rg_limit, rg_limit)
        # 固相：对物理上空的单元，清零固相步长，防止无进流时产生残差尖峰
        sol_start = offset + solid_offset(cell)
        sol_slice = slice(sol_start, sol_start + nv_sol)
        if _is_solid_empty_cell(cell, ref_solid_kg_s):
            if zero_empty_solid:
                # GD fallback：允许负步（排干虚假固相库存），阻止正步（防止向空单元充填）。
                # 完全清零会让 GD 无法减小空单元的固相残差，导致求解器停滞。
                out[sol_slice] = np.minimum(out[sol_slice], 0.0)
            else:
                # NR 方向：用紧致限制（0.1% ref），允许 NR 缓慢建立固相库存，
                # 同时防止过大步长导致线搜索失败
                tight_limit = 1e-3 * max(_solid_reference_for_cell(cell, ref_solid_kg_s), 1.0)
                out[sol_slice] = np.clip(out[sol_slice], -tight_limit, tight_limit)
        else:
            out[sol_slice] = np.clip(out[sol_slice], -rs_limit, rs_limit)
        if nv_sol > 0:
            # Keep the line-search positivity projection from silently deleting a
            # source-supported solid inventory.  The projected trial state clips
            # negative masses to zero, so an over-large negative Newton component
            # can otherwise erase char/ash holdup and still reduce gas RMS.
            solid_state = cell.m_solid.reshape(-1)[solid_flat_indices_for_nr(cell)]
            out[sol_slice] = np.maximum(out[sol_slice], -0.8 * np.maximum(solid_state, 0.0))
        if t_step_limit_K is not None and has_temperature_var(cell):
            # Structured solves bypass the sparse-direct Tikhonov fallback, so keep
            # the same physical trust region here as a final common guardrail.
            temp_idx = offset + temperature_offset(cell)
            limit = float(max(t_step_limit_K, 1e-12))
            out[temp_idx] = float(np.clip(out[temp_idx], -limit, limit))
        offset += nv
    return out


def _temperature_bounds_for_cell(cell: Cell) -> tuple[float, float]:
    """Return physical/Vorabrechnung temperature bounds for one NR cell."""
    lo = float(getattr(cell, "_nr_temperature_min_K", 300.0))
    hi = float(getattr(cell, "_nr_temperature_max_K", 2500.0))
    if not np.isfinite(lo):
        lo = 300.0
    if not np.isfinite(hi):
        hi = 2500.0
    lo = max(lo, 300.0)
    hi = min(max(hi, lo + 1.0), 2500.0)
    return lo, hi


def _clip_dx_diagnostics(dx: np.ndarray, dx_clipped: np.ndarray, cells: List[Cell]) -> dict:
    """Summarize how strongly the Newton step was clipped this iteration."""
    clipped_mask = np.abs(dx_clipped - dx) > 1e-12
    temp_clipped_cells = 0
    max_abs_dT_step = 0.0
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        if has_temperature_var(cell):
            temp_idx = offset + temperature_offset(cell)
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


def _next_lambda_seed(lambda_cap: float, accepted_lambda: float, trial_evals: int, reduction_ratio: float) -> float:
    """Adapt the next line-search start from the last accepted damping factor."""
    lam_cap = float(max(lambda_cap, 1e-12))
    lam_acc = float(np.clip(accepted_lambda, 1e-12, lam_cap))
    rr = float(np.clip(reduction_ratio, 0.0, 10.0))
    if int(trial_evals) <= 1:
        if rr < 0.4:
            return float(min(lam_cap, 2.0 * lam_acc))
        if rr < 0.7:
            return float(min(lam_cap, 1.5 * lam_acc))
        return float(min(lam_cap, 1.2 * lam_acc))
    if int(trial_evals) <= 2 and rr < 0.8:
        return float(min(lam_cap, 1.2 * lam_acc))
    return float(max(0.8 * lam_acc, 1e-12))


def _solve_linear_step(
    J: sp.csr_matrix,
    rhs: np.ndarray,
    *,
    zero_rows: int = 0,
    zero_cols: int = 0,
    structured_jacobian=None,
    linear_solver_backend: str = "sparse_direct_fallback",
    t_dof_indices: np.ndarray | None = None,
    t_max_step_K: float = 600.0,
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
        t_dof_indices=t_dof_indices,
        t_max_step_K=t_max_step_K,
    )
    return dx, method, regularized, {
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "schur_size": 0,
        "fallback_used": method != "splu",
    }


def _solve_lm_step(
    J: sp.csr_matrix,
    F_hat: np.ndarray,
    *,
    lm_mu: float,
    t_dof_indices: np.ndarray | None = None,
    t_max_step_K: float = 600.0,
) -> tuple[np.ndarray, str, bool, dict[str, float | int | bool]]:
    """Levenberg–Marquardt step: solve ``(J^T J + μI) dx = -J^T F_hat``."""
    mu = float(max(lm_mu, 1.0e-16))
    H = (J.T @ J).tocsr()
    g = J.T @ F_hat
    H_mu = H + sp.eye(H.shape[0], format="csr") * mu
    dx, method, regularized = solve_sparse_direct_linear_step(
        H_mu,
        -g,
        regularize=False,
        t_dof_indices=t_dof_indices,
        t_max_step_K=t_max_step_K,
    )
    return dx, f"lm_{method}", regularized, {
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "schur_size": 0,
        "fallback_used": method != "splu",
        "lm_mu": mu,
    }


def _solve_ptc_step(
    J: sp.csr_matrix,
    F_hat: np.ndarray,
    *,
    ptc_alpha: float,
    t_dof_indices: np.ndarray | None = None,
    t_max_step_K: float = 600.0,
) -> tuple[np.ndarray, str, bool, dict[str, float | int | bool]]:
    """Pseudo-transient Newton step: solve ``(J + αI) dx = -F_hat``."""
    alpha = float(max(ptc_alpha, 1.0e-16))
    J_eff = J + sp.eye(J.shape[0], format="csr") * alpha
    dx, method, regularized = solve_sparse_direct_linear_step(
        J_eff,
        -F_hat,
        regularize=False,
        t_dof_indices=t_dof_indices,
        t_max_step_K=t_max_step_K,
    )
    return dx, f"ptc_{method}", regularized, {
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "schur_size": 0,
        "fallback_used": method != "splu",
        "ptc_alpha": alpha,
    }


def _solve_equilibrated_newton_step(
    J: sp.csr_matrix,
    rhs: np.ndarray,
    *,
    zero_rows: int = 0,
    zero_cols: int = 0,
    equil_iters: int = 3,
    scale_clip: float = 1.0e3,
    t_dof_indices: np.ndarray | None = None,
    t_max_step_K: float = 600.0,
) -> tuple[np.ndarray, str, bool, dict[str, float | int | bool]]:
    """Solve Newton step after Ruiz-style row/column equilibration."""
    n = int(J.shape[0])
    if n == 0:
        return np.zeros(0, dtype=np.float64), "equil_empty", False, {
            "band_lu_s": 0.0,
            "side_update_s": 0.0,
            "schur_size": 0,
            "fallback_used": False,
        }
    if int(equil_iters) <= 0:
        dx, method, regularized = solve_sparse_direct_linear_step(
            J,
            rhs,
            regularize=(int(zero_rows) > 0 or int(zero_cols) > 0),
            t_dof_indices=t_dof_indices,
            t_max_step_K=t_max_step_K,
        )
        return dx, f"equil0_{method}", regularized, {
            "band_lu_s": 0.0,
            "side_update_s": 0.0,
            "schur_size": 0,
            "fallback_used": method != "splu",
            "equil_iters": 0,
        }

    J_eq = J.tocsr().astype(np.float64)
    rhs_eq = np.asarray(rhs, dtype=np.float64).copy()
    d_col = np.ones(n, dtype=np.float64)
    eps = 1.0e-12
    sclip = float(max(scale_clip, 1.0))

    for _ in range(int(equil_iters)):
        row_norm = np.sqrt(np.maximum(np.asarray(J_eq.power(2).sum(axis=1)).ravel(), eps))
        d_row = np.clip(1.0 / row_norm, 1.0 / sclip, sclip)
        J_eq = sp.diags(d_row, format="csr") @ J_eq
        rhs_eq *= d_row

        col_norm = np.sqrt(np.maximum(np.asarray(J_eq.power(2).sum(axis=0)).ravel(), eps))
        d_c = np.clip(1.0 / col_norm, 1.0 / sclip, sclip)
        J_eq = J_eq @ sp.diags(d_c, format="csr")
        d_col *= d_c

    dx_eq, method, regularized = solve_sparse_direct_linear_step(
        J_eq,
        rhs_eq,
        regularize=(int(zero_rows) > 0 or int(zero_cols) > 0),
        t_dof_indices=t_dof_indices,
        t_max_step_K=t_max_step_K,
    )
    dx = d_col * dx_eq
    return dx, f"equil_{method}", regularized, {
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "schur_size": 0,
        "fallback_used": method != "splu",
        "equil_iters": int(equil_iters),
        "equil_scale_clip": sclip,
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
    n_damp_halvings: int = 12,
    lambda_min: float = 1.0 / 1024.0,
    jacobian_strategy: str = "block_tridiag_structured",
    linear_solver_backend: str | None = None,
    jacobian_lag_steps: int = 1,
    allow_gd_fallback: bool = False,
    prefer_full_step: bool = False,
    step_model: str = "newton",
    lm_mu0: float = 1.0e-4,
    lm_mu_growth: float = 10.0,
    ptc_alpha0: float = 1.0,
    ptc_alpha_growth: float = 2.0,
    equil_iters: int = 0,
    equil_scale_clip: float = 1.0e3,
    line_search_max_trials: int = 0,
    nonmonotone_enabled: bool = False,
    nonmonotone_window: int = 5,
    nonmonotone_relax: float = 1.0,
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
    group_norms = {"gas": float("inf"), "solid": float("inf"), "energy": float("inf")}
    group_max_abs = {"gas": float("inf"), "solid": float("inf"), "energy": float("inf")}

    t0 = perf_counter()
    F = global_residual(x, cells, apply_bc_fn)
    timing["initial_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] += 1
    x = pack_reactor(cells)
    scale = build_equation_scales(cells, ref_gas_mol_s, ref_solid_kg_s, ref_energy_W)

    n_total = len(F)

    F_hat = F / scale
    norm_F = _rms_norm(F_hat)
    history: List[float] = [norm_F]
    converged = False
    lambda_seed = 1.0 if bool(prefer_full_step) else float(lambda_init)
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
    step_model_resolved = str(step_model).strip().lower()
    if step_model_resolved not in {"newton", "lm", "ptc", "equil_newton"}:
        raise ValueError(f"Unsupported step_model={step_model!r}; expected 'newton', 'lm', 'ptc', or 'equil_newton'")
    lm_mu0_resolved = float(max(lm_mu0, 1.0e-12))
    lm_mu_growth_resolved = float(max(lm_mu_growth, 1.0))
    ptc_alpha0_resolved = float(max(ptc_alpha0, 1.0e-12))
    ptc_alpha_growth_resolved = float(max(ptc_alpha_growth, 1.0))
    equil_iters_resolved = int(max(equil_iters, 0))
    equil_scale_clip_resolved = float(max(equil_scale_clip, 1.0))
    line_search_max_trials_resolved = int(max(line_search_max_trials, 0))
    nonmonotone_enabled_resolved = bool(nonmonotone_enabled)
    nonmonotone_window_resolved = int(max(nonmonotone_window, 1))
    nonmonotone_relax_resolved = float(max(nonmonotone_relax, 1.0))

    # Pre-compute T DOF indices (last DOF of each cell block).
    # Used by _solve_linear_step to pass Tikhonov regularisation info.
    _t_dof_idx: np.ndarray | None = None
    _t_dof_list: list[int] = []
    _off = 0
    for _c in cells:
        _nv = n_var(_c)
        if has_temperature_var(_c):
            _t_dof_list.append(_off + temperature_offset(_c))
        _off += _nv
    if _t_dof_list:
        _t_dof_idx = np.array(_t_dof_list, dtype=np.intp)

    for it in range(max_iter):
        # 总是按真实 cell-by-cell 布局计算三类残差 RMS，用于自适应 T 步长判断。
        group_metrics = _residual_group_metrics(F_hat, cells)
        group_norms = {
            "gas": group_metrics["gas_rms"],
            "solid": group_metrics["solid_rms"],
            "energy": group_metrics["energy_rms"],
        }
        group_max_abs = {
            "gas": group_metrics["gas_max_abs"],
            "solid": group_metrics["solid_max_abs"],
            "energy": group_metrics["energy_max_abs"],
        }
        gas_norm = group_norms["gas"]
        sol_norm = group_norms["solid"]
        ene_norm = group_norms["energy"]
        if verbose:
            print(
                f"  NR iter {it}: RMS={norm_F:.3e}  gas={gas_norm:.3e}  "
                f"solid={sol_norm:.3e}  energy={ene_norm:.3e}"
            )

        # 自适应 T 步长上限：当能量方程已充分满足而气/固残差仍大时，
        # 收紧 Tikhonov 的 T 限制，避免 T 扰动破坏已收敛的能量平衡。
        # 仅在能量残差比收敛目标低 3 个数量级时（< 0.1% of tol_rms）才生效，
        # 避免在气相迭代中过早固定 T（导致气相 → 能量方程违反的振荡）。
        _dominant_non_energy = max(gas_norm, sol_norm, 1e-30)
        _energy_frac = ene_norm / _dominant_non_energy
        _t_max_K_adaptive: float
        if ene_norm < tol_rms * 1e-3 and _energy_frac < 1e-3:
            # 能量深度收敛（<0.1% tol_rms）：将 T 步长限制到 1K
            _t_max_K_adaptive = 1.0
        else:
            _t_max_K_adaptive = 600.0

        x_rms = _state_rms_norm(x)
        component_norm = max(gas_norm, sol_norm, ene_norm)
        component_max_abs = max(group_max_abs.values())
        if inner_convergence is not None:
            if (
                inner_convergence.is_converged(norm_F, x_rms)
                and component_norm <= tol_rms
                and component_max_abs <= tol_rms
            ):
                converged = True
                break
        elif norm_F < tol_rms and component_norm <= tol_rms and component_max_abs <= tol_rms:
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
            # 新 Jacobian 意味着全新牛顿方向；重置 lambda_seed 到初始值，
            # 以便线搜索从大步长开始探索，不沿用上次成功的小步长。
            lambda_seed = 1.0 if bool(prefer_full_step) else float(lambda_init)
        else:
            J = cached_J
            jac_meta = cached_jac_meta or {}
            counts["jacobian_reuses"] += 1

        # 使用稀疏解法求解 Newton 方向 Δx；对秩亏 Jacobian 做轻度对角正则化。
        t0 = perf_counter()
        if step_model_resolved == "lm":
            lm_mu_iter = lm_mu0_resolved * (lm_mu_growth_resolved ** max(consecutive_line_search_failures, 0))
            dx, linear_solver_used, linear_regularized, linear_diag = _solve_lm_step(
                J,
                F_hat,
                lm_mu=lm_mu_iter,
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
            counts["lm_steps"] = int(counts.get("lm_steps", 0)) + 1
        elif step_model_resolved == "ptc":
            ptc_alpha_iter = ptc_alpha0_resolved * (ptc_alpha_growth_resolved ** max(consecutive_line_search_failures, 0))
            dx, linear_solver_used, linear_regularized, linear_diag = _solve_ptc_step(
                J,
                F_hat,
                ptc_alpha=ptc_alpha_iter,
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
            counts["ptc_steps"] = int(counts.get("ptc_steps", 0)) + 1
        elif step_model_resolved == "equil_newton":
            dx, linear_solver_used, linear_regularized, linear_diag = _solve_equilibrated_newton_step(
                J,
                -F_hat,
                zero_rows=int(jac_meta.get("zero_rows", 0)),
                zero_cols=int(jac_meta.get("zero_cols", 0)),
                equil_iters=equil_iters_resolved,
                scale_clip=equil_scale_clip_resolved,
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
            counts["equil_steps"] = int(counts.get("equil_steps", 0)) + 1
        else:
            dx, linear_solver_used, linear_regularized, linear_diag = _solve_linear_step(
                J,
                -F_hat,
                zero_rows=int(jac_meta.get("zero_rows", 0)),
                zero_cols=int(jac_meta.get("zero_cols", 0)),
                structured_jacobian=jac_meta.get("structured_jacobian"),
                linear_solver_backend=str(linear_solver_backend),
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
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
        dx_clipped = _clip_dx(dx, cells, ref_gas_mol_s, ref_solid_kg_s, t_step_limit_K=_t_max_K_adaptive)
        clip_diag = _clip_dx_diagnostics(dx, dx_clipped, cells)

        # --- 阻尼线搜索 (Damped Newton) ---
        lambda_cap = 1.0 if bool(prefer_full_step) else max(float(lambda_init), lambda_min)
        # Hamel-style damped Newton starts fresh Newton directions from full step.
        # A line-search retry, however, is the same accepted state and the same
        # linearisation; restart it from the best smaller lambda instead of
        # rebuilding an identical FD Jacobian and trying lambda=1 again.
        if bool(prefer_full_step) and consecutive_line_search_failures == 0:
            lam = 1.0
        else:
            lam = float(np.clip(lambda_seed, lambda_min, lambda_cap))
        nt = norm_F
        found_step = False
        trial_evals = 0
        best_trial_nt = float(np.inf)
        best_trial_lambda: float | None = None
        nm_recent_max = float(max(history[-nonmonotone_window_resolved:])) if history else float(norm_F)
        nm_accept_upper = float(nm_recent_max * nonmonotone_relax_resolved)
        ls_trial_cap = (
            min(int(n_damp_halvings), line_search_max_trials_resolved)
            if line_search_max_trials_resolved > 0
            else int(n_damp_halvings)
        )

        for damp_step in range(ls_trial_cap):
            x_trial = x + lam * dx_clipped
            temp_bound_hits = 0
            # 对尝试步进行物理边界保护（必须正值，且 T 在范围内）
            offset = 0
            for cell in cells:
                nv = n_var(cell)
                if has_temperature_var(cell):
                    temp_local = temperature_offset(cell)
                    # 气/固质量流率 >= 0
                    x_trial[offset : offset + temp_local] = np.maximum(
                        x_trial[offset : offset + temp_local],
                        0.0,
                    )
                    # 温度保护
                    raw_T_trial = float(x_trial[offset + temp_local])
                    t_min, t_max = _temperature_bounds_for_cell(cell)
                    clipped_T_trial = float(np.clip(raw_T_trial, t_min, t_max))
                    temp_bound_hits += int(abs(clipped_T_trial - raw_T_trial) > 1e-12)
                    x_trial[offset + temp_local] = clipped_T_trial
                else:
                    x_trial[offset : offset + nv] = np.maximum(x_trial[offset : offset + nv], 0.0)
                offset += nv

            # Evaluate projected trial states even when temperature bounds are hit.
            # Rejecting all bound-hit trials can starve line-search (0 evaluations)
            # when one cell sits near a hard bound (e.g. T≈300 K).

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
            if verbose and trial_evals <= 3:
                print(f"    LS trial {damp_step}: lam={lam:.4e}  rms={nt:.4e}  (norm_F={norm_F:.4e})")
            
            accepted_nonmonotone = bool(
                nonmonotone_enabled_resolved
                and consecutive_line_search_failures > 0
                and nt >= norm_F
                and np.isfinite(nt)
                and nt <= nm_accept_upper
            )
            if nt < norm_F or accepted_nonmonotone:
                if verbose and damp_step > 0:
                    print(f"    Line search OK at step {damp_step}: lam={lam:.4f}, norm={nt:.3e}")
                found_step = True
                consecutive_line_search_failures = 0
                prev_norm_F = float(norm_F)
                x = pack_reactor(cells)
                F = F_trial
                F_hat = F_hat_trial
                norm_F = nt
                counts["line_search_evaluations"] += int(trial_evals)
                counts["line_search_backtracks"] += max(int(trial_evals) - 1, 0)
                if accepted_nonmonotone:
                    counts["nonmonotone_accepts"] = int(counts.get("nonmonotone_accepts", 0)) + 1
                accepted_lambda_history.append(float(lam))
                line_search_trial_counts.append(int(trial_evals))
                clip_history.append(
                    {
                        "iter": int(it + 1),
                        "accepted_lambda": float(lam),
                        "line_search_trials": int(trial_evals),
                        "line_search_failed": False,
                        "nonmonotone_accepted": bool(accepted_nonmonotone),
                        "best_trial_rms_scaled": float(best_trial_nt),
                        "best_trial_lambda": best_trial_lambda,
                        **clip_diag,
                    }
                )
                lambda_seed = _next_lambda_seed(
                    float(lambda_cap),
                    float(lam),
                    int(trial_evals),
                    float(nt / max(prev_norm_F, 1.0e-30)),
                )
                break
            
            lam *= 0.5
            if lam < lambda_min:
                break

        if not found_step:
            if line_search_max_trials_resolved > 0 and int(trial_evals) >= int(ls_trial_cap):
                counts["line_search_cap_hits"] = int(counts.get("line_search_cap_hits", 0)) + 1
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
                lambda_seed = float(np.clip(0.5 * lam_retry_base, lambda_min, lambda_cap))
                # 恢复 cells 到接受状态（line search trial 会把 cells 留在最后一次试探）
                _ = global_residual(x, cells, apply_bc_fn)
                continue

            if not bool(allow_gd_fallback):
                # Hamel-aligned damped Newton path: stop after line-search failure
                # (after one smaller-lambda retry) instead of switching to a
                # gradient-descent surrogate direction.
                _ = global_residual(x, cells, apply_bc_fn)
                counts["global_residual_calls"] += 1
                break

            # --- 梯度下降兜底 (Gradient-descent fallback) ---
            # 当 NR 线搜索彻底失败（lstsq 方向在当前状态为反下降方向）时，
            # 使用 dx_gd = -J^T F_hat 作为替代方向。
            # 该方向保证是 ||F_hat||^2 的下降方向：
            #   d/dλ ||F_hat(x + λ dx_gd)||^2|_{λ=0} = -2||J^T F_hat||^2 ≤ 0
            # 即使 Jacobian 近奇异（holdup_transport + frozen Vorabrechnung），
            # 只要 J^T F_hat ≠ 0，梯度方向就是下降方向。
            # 恢复 cells 到接受状态（line search trial 会把 cells 留在最后一次试探状态）
            _ = global_residual(x, cells, apply_bc_fn)
            J_arr = J.toarray() if sp.issparse(J) else np.asarray(J)
            dx_gd = -(J_arr.T @ F_hat)
            # 当能量已满足时，清零 T 分量：任何 T 扰动都会破坏已收敛的能量方程。
            _energy_satisfied_for_gd = ene_norm < tol_rms * 1e-3 and _energy_frac < 1e-3
            if _energy_satisfied_for_gd and _t_dof_idx is not None:
                dx_gd[_t_dof_idx] = 0.0
            gd_norm = float(np.linalg.norm(dx_gd))
            if verbose:
                print(f"  GD fallback: gd_norm={gd_norm:.3e}, norm_F={norm_F:.3e}, freeze_T={_energy_satisfied_for_gd}")
            if gd_norm > 1e-10:
                # 将方向归一化后按气相参考量缩放，确保步长物理合理
                gd_scale = 0.1 * float(ref_gas_mol_s) / gd_norm
                dx_gd_clipped = _clip_dx(
                    gd_scale * dx_gd,
                    cells,
                    ref_gas_mol_s,
                    ref_solid_kg_s,
                    zero_empty_solid=True,
                    t_step_limit_K=_t_max_K_adaptive,
                )
                if verbose:
                    print(f"    gd_scale={gd_scale:.3e}, max|dx_gd_clip|={np.max(np.abs(dx_gd_clipped)):.3e}")
                clip_diag_gd = _clip_dx_diagnostics(dx_gd_clipped, dx_gd_clipped, cells)
                lam_gd = float(lambda_init)
                found_gd = False
                gd_trial_evals = 0
                # 追踪线搜索中最优试探（用于非单调兜底）
                best_gd_nt = float(np.inf)
                best_gd_x: np.ndarray | None = None
                best_gd_F: np.ndarray | None = None
                best_gd_F_hat: np.ndarray | None = None
                best_gd_lam: float = float(lambda_init)

                for _ in range(n_damp_halvings * 2):
                    x_trial_gd = x + lam_gd * dx_gd_clipped
                    offset = 0
                    gd_temp_hits = 0
                    for cell in cells:
                        nv = n_var(cell)
                        if has_temperature_var(cell):
                            temp_local = temperature_offset(cell)
                            x_trial_gd[offset : offset + temp_local] = np.maximum(
                                x_trial_gd[offset : offset + temp_local], 0.0
                            )
                            raw_T = float(x_trial_gd[offset + temp_local])
                            clipped_T = float(np.clip(raw_T, 300.0, 2500.0))
                            gd_temp_hits += int(abs(clipped_T - raw_T) > 1e-12)
                            x_trial_gd[offset + temp_local] = clipped_T
                        else:
                            x_trial_gd[offset : offset + nv] = np.maximum(
                                x_trial_gd[offset : offset + nv], 0.0
                            )
                        offset += nv
                    if gd_temp_hits > 0:
                        if lam_gd <= lambda_min:
                            break
                        lam_gd *= 0.5
                        continue
                    F_gd = global_residual(x_trial_gd, cells, apply_bc_fn)
                    counts["global_residual_calls"] += 1
                    gd_trial_evals += 1
                    F_hat_gd = F_gd / scale
                    nt_gd = _rms_norm(F_hat_gd)
                    if np.isfinite(nt_gd) and nt_gd < best_gd_nt:
                        best_gd_nt = float(nt_gd)
                        best_gd_x = x_trial_gd.copy()
                        best_gd_F = F_gd.copy()
                        best_gd_F_hat = F_hat_gd.copy()
                        best_gd_lam = float(lam_gd)
                    if verbose:
                        print(f"    GD lam={lam_gd:.4f}  rms={nt_gd:.4e}  {'<<ACCEPT' if np.isfinite(nt_gd) and nt_gd < norm_F else ''}")
                        if gd_trial_evals == 1:
                            # First GD trial: diagnose which residual component is largest
                            big_idx = int(np.argmax(np.abs(F_hat_gd)))
                            print(
                                f"    [GD diag] biggest F_hat: idx={big_idx} val={F_hat_gd[big_idx]:.3e} "
                                f"x={x_trial_gd[big_idx]:.3e} prev_F_hat={F_hat[big_idx]:.3e}"
                            )
                            offsets_diag = cell_offsets(cells)
                            cell_idx = max(0, int(np.searchsorted(offsets_diag, big_idx, side="right") - 1)) if cells else 0
                            print(f"    [GD diag] cell≈{cell_idx}, change_at_idx={x_trial_gd[big_idx]-x[big_idx]:.3e}")
                    if np.isfinite(nt_gd) and nt_gd < norm_F:
                        # 严格下降
                        x = x_trial_gd
                        F = F_gd
                        F_hat = F_hat_gd
                        norm_F = nt_gd
                        found_gd = True
                        consecutive_line_search_failures = 0
                        counts.setdefault("gd_fallback_accepts", 0)
                        counts["gd_fallback_accepts"] += 1
                        accepted_lambda_history.append(float(lam_gd))
                        line_search_trial_counts.append(int(gd_trial_evals))
                        clip_history.append(
                            {
                                "iter": int(it + 1),
                                "accepted_lambda": float(lam_gd),
                                "line_search_trials": int(gd_trial_evals),
                                "line_search_failed": False,
                                "gd_fallback": True,
                                "best_trial_rms_scaled": float(nt_gd),
                                "best_trial_lambda": float(lam_gd),
                                **clip_diag_gd,
                            }
                        )
                        lambda_seed = float(lambda_min)
                        break
                    lam_gd *= 0.5
                    if lam_gd < lambda_min:
                        break

                if not found_gd and best_gd_x is not None:
                    # 非单调有界兜底：若 GD 最优试探比当前 rms 仅高 ≤5%，
                    # 则接受之（有界 non-monotone，防止发散），以便外层继续推进。
                    nm_rel_tol = 1.05
                    if best_gd_nt < nm_rel_tol * norm_F:
                        _old_norm_F = norm_F
                        x = best_gd_x
                        F = best_gd_F
                        F_hat = best_gd_F_hat  # type: ignore[assignment]
                        norm_F = best_gd_nt
                        found_gd = True
                        consecutive_line_search_failures = 0
                        counts.setdefault("gd_nm_accepts", 0)
                        counts["gd_nm_accepts"] += 1
                        accepted_lambda_history.append(float(best_gd_lam))
                        line_search_trial_counts.append(int(gd_trial_evals))
                        clip_history.append(
                            {
                                "iter": int(it + 1),
                                "accepted_lambda": float(best_gd_lam),
                                "line_search_trials": int(gd_trial_evals),
                                "line_search_failed": False,
                                "gd_fallback": True,
                                "gd_nonmonotone": True,
                                "best_trial_rms_scaled": float(best_gd_nt),
                                "best_trial_lambda": float(best_gd_lam),
                                **clip_diag_gd,
                            }
                        )
                        if verbose:
                            print(
                                f"  GD non-monotone accept: rms {norm_F:.4e} "
                                f"(ratio={norm_F / _old_norm_F:.3f})"
                            )
                        lambda_seed = float(lambda_min)
                if found_gd:
                    history.append(norm_F)
                    continue  # 继续下一个 NR 迭代
            break

        history.append(norm_F)

    # 结果解包回反应器对象
    t0 = perf_counter()
    F_final = global_residual(x, cells, apply_bc_fn)
    timing["final_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] += 1
    x = pack_reactor(cells)
    final_norm = float(np.linalg.norm(F_final))
    final_group_metrics = _residual_group_metrics(F_final / scale, cells)
    final_gas_phase_metrics = _residual_gas_phase_metrics(F_final / scale, cells)
    timing["total_s"] = perf_counter() - solve_started

    return {
        "converged": converged,
        "n_iter": len(history),
        "residual": final_norm,
        "rms_scaled_final": norm_F,
        "rms_scaled_gas_final": final_group_metrics["gas_rms"],
        "rms_scaled_solid_final": final_group_metrics["solid_rms"],
        "rms_scaled_energy_final": final_group_metrics["energy_rms"],
        "rms_scaled_gas_combined_final": final_gas_phase_metrics["gas_combined_rms"],
        "rms_scaled_gas_phase_split_final": final_gas_phase_metrics["gas_phase_split_rms"],
        "rms_scaled_component_max_final": max(
            final_group_metrics["gas_rms"],
            final_group_metrics["solid_rms"],
            final_group_metrics["energy_rms"],
        ),
        "max_abs_scaled_gas_final": final_group_metrics["gas_max_abs"],
        "max_abs_scaled_solid_final": final_group_metrics["solid_max_abs"],
        "max_abs_scaled_energy_final": final_group_metrics["energy_max_abs"],
        "max_abs_scaled_gas_combined_final": final_gas_phase_metrics["gas_combined_max_abs"],
        "max_abs_scaled_gas_phase_split_final": final_gas_phase_metrics["gas_phase_split_max_abs"],
        "max_abs_scaled_final": max(
            final_group_metrics["gas_max_abs"],
            final_group_metrics["solid_max_abs"],
            final_group_metrics["energy_max_abs"],
        ),
        "norm_history": history,
        "n_newton_iters_attempted": int(iter_attempts),
        "jacobian_strategy": jacobian_strategy,
        "nr_layout_audit": _nr_layout_audit(cells, F_final),
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
        "step_model": step_model_resolved,
        "lm_mu0": lm_mu0_resolved,
        "lm_mu_growth": lm_mu_growth_resolved,
        "ptc_alpha0": ptc_alpha0_resolved,
        "ptc_alpha_growth": ptc_alpha_growth_resolved,
        "equil_iters": equil_iters_resolved,
        "equil_scale_clip": equil_scale_clip_resolved,
        "line_search_max_trials": line_search_max_trials_resolved,
        "nonmonotone_enabled": nonmonotone_enabled_resolved,
        "nonmonotone_window": nonmonotone_window_resolved,
        "nonmonotone_relax": nonmonotone_relax_resolved,
        "jacobian_lag_steps": int(lag_steps),
        "timing": timing,
        "counts": counts,
        "accepted_lambda_history": accepted_lambda_history,
        "line_search_trial_counts": line_search_trial_counts,
        "clip_history": clip_history,
    }
