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
from typing import Any, List

import numpy as np
import scipy.sparse as sp

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.convergence import InnerConvergence
from src.solvers.nr_indexing import (
    affected_residual_cells_for_var_cell,
    build_var_index_to_cell_map,
    cell_offsets,
    gas_mode_for_nr,
    has_temperature_var,
    main_nr_gas_idx_for_cell,
    main_nr_gas_phase_split_positions_for_cell,
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
    assemble_jacobian_hybrid_structured,
    build_jacobian_structure,
    build_solver_graph,
    solve_sparse_direct_linear_step,
    solve_structured_linear_step,
    structured_strategy_name_from_cells,
)

_FD_EPS = 1e-6

_build_var_index_to_cell_map = build_var_index_to_cell_map
_affected_residual_cells_for_var_cell = affected_residual_cells_for_var_cell


def _is_bed_two_phase_gas_cell(cell: Cell) -> bool:
    """Cells whose NR layout keeps distinct dense and bubble gas balances."""
    return gas_mode_for_nr(cell) == "two_phase"


def snap_trace_o2_holdups(cells: List[Cell]) -> int:
    """将无入口支撑的微量 O₂ holdup 粘滞到 0（handoff §5.13）。

    同构裸动力学下，N_O2≈1e-16 即可把 bed 内 R1/R12 等速率推到 O(10²⁺) mol/s，
    而真实贫氧格点应为精确零分支。仅当 N_in/N_zu 亦低于阈值时才清零，
    避免挡住有入口 O₂ 的穿透。

    Returns
    -------
    int
        被清零的相持有量个数（dense/bubble 各计一次）。
    """
    if not cells:
        return 0
    j = int(GAS_SPECIES_INDEX["O2"])
    n_hit = 0
    for cell in cells:
        if not bool(getattr(cell, "_nr_snap_trace_o2", False)):
            continue
        atol = float(getattr(cell, "_nr_snap_trace_o2_atol_mol_s", 1.0e-12))
        atol = max(atol, 0.0)
        nin = float(
            cell.N_d_in[j]
            + cell.N_b_in[j]
            + cell.N_zu_d[j]
            + cell.N_zu_b[j]
        )
        if nin > atol:
            continue
        for arr in (cell.N_d, cell.N_b):
            val = float(arr[j])
            if val != 0.0 and val <= atol:
                arr[j] = 0.0
                n_hit += 1
    return int(n_hit)


def _install_bed_scoped_flag(
    cells: List[Cell],
    *,
    enabled: bool,
    atol: float,
    min_bed: int,
    flag_attr: str,
    atol_attr: str,
) -> None:
    """按 bed 序号挂布尔旗标 + atol（snap / rate-gate 共用）。"""
    atol = float(max(atol, 0.0))
    min_bed = int(max(min_bed, 0))
    bed_i = -1
    for cell in cells:
        setattr(cell, atol_attr, atol)
        apply = bool(enabled)
        if apply and str(getattr(cell, "cell_type", "bed")) == "bed":
            bed_i += 1
            apply = bed_i >= min_bed
        elif apply and min_bed > 0:
            apply = False
        setattr(cell, flag_attr, bool(apply))


def install_trace_o2_snap_on_cells(cells: List[Cell], config: Any | None) -> None:
    """把 trace-O₂ 粘滞开关挂到 solver cells（供 residual / clip 读取）。

    ``nr_snap_trace_o2_min_bed_index_thesis``：仅 bed 序号 ≥ 该值的床层启用
    （用于避开 bed0/1 氧化主战场，专治 bed2+ 贫氧速率面跳跃）。
    """
    enabled = bool(getattr(config, "nr_snap_trace_o2_holdup_thesis", False)) if config else False
    atol = (
        float(getattr(config, "nr_snap_trace_o2_atol_mol_s_thesis", 1.0e-12))
        if config
        else 1.0e-12
    )
    min_bed = (
        int(getattr(config, "nr_snap_trace_o2_min_bed_index_thesis", 0) or 0)
        if config
        else 0
    )
    _install_bed_scoped_flag(
        cells,
        enabled=enabled,
        atol=atol,
        min_bed=min_bed,
        flag_attr="_nr_snap_trace_o2",
        atol_attr="_nr_snap_trace_o2_atol_mol_s",
    )


def install_trace_o2_rate_gate_on_cells(cells: List[Cell], config: Any | None) -> None:
    """速率评价时把微量 O₂ 当零（不动 holdup；handoff §5.15）。

    与 snap 不同：Newton 状态可保留 1e-16 灰尘，但 R1/R5/R12 等走精确零分支。
    若 ``nr_trace_o2_rate_gate_after_y_co_thesis``：先只挂 atol/作用域，
    待 ``arm_trace_o2_rate_gate_if_y_co`` 在 CO 带宽后武装。
    """
    enabled = bool(getattr(config, "nr_trace_o2_rate_gate_thesis", False)) if config else False
    defer = bool(getattr(config, "nr_trace_o2_rate_gate_after_y_co_thesis", False)) if config else False
    atol = (
        float(getattr(config, "nr_trace_o2_rate_gate_atol_mol_s_thesis", 1.0e-12))
        if config
        else 1.0e-12
    )
    min_bed = (
        int(getattr(config, "nr_trace_o2_rate_gate_min_bed_index_thesis", 2) or 0)
        if config
        else 2
    )
    # 延迟模式：安装作用域元数据，但暂不启用门控
    active_now = bool(enabled) and (not defer)
    _install_bed_scoped_flag(
        cells,
        enabled=active_now,
        atol=atol,
        min_bed=min_bed,
        flag_attr="_nr_trace_o2_rate_gate",
        atol_attr="_nr_trace_o2_rate_gate_atol_mol_s",
    )
    # 记录延迟武装所需的作用域（即使当前未启用）
    bed_i = -1
    for cell in cells:
        in_scope = False
        if enabled and str(getattr(cell, "cell_type", "bed")) == "bed":
            bed_i += 1
            in_scope = bed_i >= int(max(min_bed, 0))
        cell._nr_trace_o2_rate_gate_in_scope = bool(in_scope)
        cell._nr_trace_o2_rate_gate_deferred = bool(enabled and defer)


def _bed_top_y_co(cells: List[Cell]) -> float:
    bed = [c for c in cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    if not bed:
        return float("nan")
    top = bed[-1]
    j_co = int(GAS_SPECIES_INDEX["CO"])
    n_co = float(top.N_d[j_co] + top.N_b[j_co])
    n_gas = float(sum(top.N_d) + sum(top.N_b))
    return float(n_co / max(n_gas, 1.0e-30))


def arm_trace_o2_rate_gate_if_y_co(cells: List[Cell], config: Any | None) -> dict[str, float | bool]:
    """CO 带宽后武装延迟 rate-gate；已武装则保持。

    Returns
    -------
    dict
        ``armed`` / ``y_co_top`` / ``newly_armed``
    """
    thr = (
        float(getattr(config, "nr_trace_o2_rate_gate_y_co_threshold_thesis", 0.14))
        if config
        else 0.14
    )
    y_co = _bed_top_y_co(cells)
    already = any(bool(getattr(c, "_nr_trace_o2_rate_gate", False)) for c in cells)
    deferred_pending = any(bool(getattr(c, "_nr_trace_o2_rate_gate_deferred", False)) for c in cells)
    if already:
        return {"armed": True, "newly_armed": False, "y_co_top": y_co}
    if (not deferred_pending) or (not np.isfinite(y_co)) or y_co < thr:
        return {"armed": False, "newly_armed": False, "y_co_top": y_co}
    for cell in cells:
        if bool(getattr(cell, "_nr_trace_o2_rate_gate_in_scope", False)):
            cell._nr_trace_o2_rate_gate = True
        cell._nr_trace_o2_rate_gate_deferred = False
    return {"armed": True, "newly_armed": True, "y_co_top": y_co}


def install_bed_top_solid_energy_passthrough_on_cells(
    cells: List[Cell], config: Any | None
) -> None:
    """床顶 Eq.2.7 固体焓透传旗标（§5.19 / §5.20）。

    ``after_y_co``：安装 deferred，待 ``arm_bed_top_solid_energy_passthrough_if_y_co``。
    """
    enabled = (
        bool(getattr(config, "nr_bed_top_solid_energy_passthrough_thesis", False))
        if config
        else False
    )
    defer = (
        bool(getattr(config, "nr_bed_top_solid_energy_passthrough_after_y_co_thesis", False))
        if config
        else False
    )
    bed = [c for c in cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    for c in cells:
        c._solid_energy_passthrough = False
        c._solid_energy_passthrough_deferred = False
        c._solid_energy_passthrough_is_top = False
    if not enabled or not bed:
        return
    top = bed[-1]
    top._solid_energy_passthrough_is_top = True
    if defer:
        top._solid_energy_passthrough_deferred = True
    else:
        top._solid_energy_passthrough = True


def arm_bed_top_solid_energy_passthrough_if_y_co(
    cells: List[Cell], config: Any | None
) -> dict[str, float | bool]:
    """CO 带宽后武装延迟床顶焓透传；已武装则保持。"""
    thr = (
        float(getattr(config, "nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis", 0.14))
        if config
        else 0.14
    )
    y_co = _bed_top_y_co(cells)
    already = any(bool(getattr(c, "_solid_energy_passthrough", False)) for c in cells)
    deferred_pending = any(
        bool(getattr(c, "_solid_energy_passthrough_deferred", False)) for c in cells
    )
    if already:
        return {"armed": True, "newly_armed": False, "y_co_top": y_co}
    if (not deferred_pending) or (not np.isfinite(y_co)) or y_co < thr:
        return {"armed": False, "newly_armed": False, "y_co_top": y_co}
    for cell in cells:
        if bool(getattr(cell, "_solid_energy_passthrough_is_top", False)):
            cell._solid_energy_passthrough = True
        cell._solid_energy_passthrough_deferred = False
    return {"armed": True, "newly_armed": True, "y_co_top": y_co}


def gate_trace_o2_for_rate_eval(
    cell: Cell,
    C_b: np.ndarray,
    C_d: np.ndarray,
    y_b: np.ndarray,
    y_d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """若 cell 启用 rate-gate 且局地贫氧，返回 O₂ 已置零的 C/y 副本。

    只看局地 ``N`` / ``C``，**不**因上游 ``N_in`` 放行——速率面取决于局地浓度；
    bed2 常有微量 ``N_in(O2)`` 但 holdup≈1e-16，按入口跳过会让门控失效（§5.15）。
    """
    if not bool(getattr(cell, "_nr_trace_o2_rate_gate", False)):
        return C_b, C_d, y_b, y_d
    j = int(GAS_SPECIES_INDEX["O2"])
    atol = float(max(getattr(cell, "_nr_trace_o2_rate_gate_atol_mol_s", 1.0e-12), 0.0))
    c_atol = float(max(getattr(cell, "_nr_trace_o2_rate_gate_C_atol_mol_m3", 1.0e-9), 0.0))
    n_hold = float(cell.N_d[j] + cell.N_b[j])
    c_hold = float(max(float(C_d[j]), 0.0) + max(float(C_b[j]), 0.0))
    if n_hold > atol and c_hold > c_atol:
        return C_b, C_d, y_b, y_d
    C_b_g = np.asarray(C_b, dtype=np.float64).copy()
    C_d_g = np.asarray(C_d, dtype=np.float64).copy()
    y_b_g = np.asarray(y_b, dtype=np.float64).copy()
    y_d_g = np.asarray(y_d, dtype=np.float64).copy()
    C_b_g[j] = 0.0
    C_d_g[j] = 0.0
    y_b_g[j] = 0.0
    y_d_g[j] = 0.0
    return C_b_g, C_d_g, y_b_g, y_d_g


def global_residual(
    x: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
) -> np.ndarray:
    """F(x)：解包 → 边界条件（含循环）→ 各 cell 残差拼接。"""
    unpack_reactor(x, cells)
    apply_bc_fn()
    if snap_trace_o2_holdups(cells) > 0:
        # 与速率评价所用 holdup 对齐，避免 x 残留 1e-16 而 F 走零分支
        x[:] = pack_reactor(cells)
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
    gas_idx = main_nr_gas_idx_for_cell(cell)
    parts: list[np.ndarray] = []
    if gas_mode == "two_phase":
        parts.extend([gas_d[gas_idx], gas_b[gas_idx]])
    elif gas_mode == "single":
        parts.append((gas_d + gas_b)[gas_idx])
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
        # Hamel coarse-PSD：单档进料时多数粒径 DOF 局部 support≈0，但 cell 总 holdup 仍 O(1)。
        # 仅用 per-DOF floor 会把 ~1e-2 kg/s 物理残差放大为 F_hat~10+，误触 merit gate。
        cell_ref = _solid_reference_for_cell(cell, ref_solid_kg_s)
        per_dof_floor = max(floor, cell_ref)
        return np.maximum(
            support.reshape(-1)[solid_flat_indices_for_nr(cell)],
            per_dof_floor,
        )
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

    from src.solvers.jacobian_transport_freeze import (
        enter_jacobian_transport_k_freeze,
        exit_jacobian_transport_k_freeze,
    )

    enter_jacobian_transport_k_freeze(cells)
    try:
        x_p = x.copy()
        for j in range(n_total):
            if verbose and j > 0 and j % 100 == 0:
                print(f"DEBUG:   Jacobian progress: {j}/{n_total}")

            unpack_reactor(x, cells)
            apply_bc_fn()

            xj = float(x[j])
            h = max(_FD_EPS * max(abs(xj), 1.0), 1e-8)
            x_p[:] = x
            x_p[j] = xj + h

            unpack_reactor(x_p, cells)
            apply_bc_fn()
            F_p = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
            residual_cell_calls += len(cells)

            dF = (F_p / eq_scale - F0s) / h

            indices = np.where(np.abs(dF) > 1e-14)[0]
            for i in indices:
                rows.append(int(i))
                cols.append(j)
                vals.append(float(dF[i]))
    finally:
        exit_jacobian_transport_k_freeze(cells)

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

    from src.solvers.jacobian_transport_freeze import (
        enter_jacobian_transport_k_freeze,
        exit_jacobian_transport_k_freeze,
    )

    enter_jacobian_transport_k_freeze(cells)
    try:
        x_p = x.copy()
        for j in range(n_total):
            if verbose and j > 0 and j % 100 == 0:
                print(f"DEBUG:   Jacobian progress: {j}/{n_total}")

            unpack_reactor(x, cells)
            apply_bc_fn()

            xj = float(x[j])
            h = max(_FD_EPS * max(abs(xj), 1.0), 1e-8)
            x_p[:] = x
            x_p[j] = xj + h

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
    finally:
        exit_jacobian_transport_k_freeze(cells)

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


def resolve_hybrid_jacobian_strategy(structured_strategy: str) -> str:
    """Map structured FD strategy name to its hybrid analytic overlay twin."""
    name = str(structured_strategy).strip()
    if name.startswith("hybrid_"):
        return name
    if name == "block_tridiag_structured":
        return "hybrid_block_tridiag_structured"
    if name == "band_plus_side_elements_structured":
        return "hybrid_band_plus_side_elements_structured"
    raise ValueError(f"Unsupported structured Jacobian strategy for hybrid mapping: {name!r}")


_STRUCTURED_DIRECT_JACOBIAN_STRATEGIES = frozenset(
    {
        "block_tridiag_structured",
        "band_plus_side_elements_structured",
        "hybrid_block_tridiag_structured",
        "hybrid_band_plus_side_elements_structured",
    }
)


def is_band_plus_side_elements_jacobian_strategy(strategy: str) -> bool:
    name = str(strategy).strip()
    return name in {
        "band_plus_side_elements_structured",
        "hybrid_band_plus_side_elements_structured",
    }


def structured_jacobian_uses_direct_linear_solver(strategy: str) -> bool:
    """Return True when the Jacobian strategy can use the structured direct LU path."""
    return str(strategy).strip() in _STRUCTURED_DIRECT_JACOBIAN_STRATEGIES


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
    connectivity_graph: Any | None = None,
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
    if strategy in {
        "block_tridiag_structured",
        "band_plus_side_elements_structured",
        "side_elements_sparse_fd",
        "hybrid_block_tridiag_structured",
        "hybrid_band_plus_side_elements_structured",
    }:
        hybrid = str(strategy).startswith("hybrid_")
        if strategy == "side_elements_sparse_fd":
            resolved_strategy = "band_plus_side_elements_structured"
        elif strategy in {"hybrid_band_plus_side_elements_structured"}:
            resolved_strategy = "band_plus_side_elements_structured"
        elif strategy in {"hybrid_block_tridiag_structured"}:
            resolved_strategy = "block_tridiag_structured"
        else:
            resolved_strategy = str(strategy)
        graph = build_solver_graph(cells)
        structure = build_jacobian_structure(graph, cells, connectivity_graph=connectivity_graph)
        use_local_callbacks = resolved_strategy != "band_plus_side_elements_structured"
        assemble_fn = assemble_jacobian_hybrid_structured if hybrid else assemble_jacobian_fd_structured
        structured = assemble_fn(
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
        reported_strategy = f"hybrid_{resolved_strategy}" if hybrid else resolved_strategy
        return structured.matrix, {
            "strategy": reported_strategy,
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
        ng = n_gas_var(cell)
        if gas_mode == "two_phase":
            half = ng // 2
            dense = F_hat[offset : offset + half]
            bubble = F_hat[offset + half : offset + ng]
            if dense.size and bubble.size:
                combined_terms.append(dense + bubble)
                split_pos = main_nr_gas_phase_split_positions_for_cell(cell)
                if split_pos.size:
                    split_terms.append(dense[split_pos] - bubble[split_pos])
        elif gas_mode == "single":
            gas = F_hat[offset : offset + ng]
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


def _line_search_merit(
    F_hat: np.ndarray,
    cells: List[Cell],
    *,
    use_gas_phase_split: bool,
    use_energy: bool = False,
) -> float:
    """线搜索标量 merit。

    Wirsum (1997/98) Eq. 2.23 / Hamel §2.3：默认仅为缩放残差 RMS（‖F̂‖₂ 族），
    试探接受要求 merit 严格下降。``use_gas_phase_split`` / ``use_energy`` 为工程桥，
    打开时与分量 RMS 取 max（非 Wirsum 原文）。

    Ref: docs/wirsum_1997_nr_solver_notes.md; handoff §6.6
    """
    merit = _rms_norm(F_hat)
    if use_gas_phase_split:
        split_rms = float(_residual_gas_phase_metrics(F_hat, cells)["gas_phase_split_rms"])
        merit = max(merit, split_rms)
    if use_energy:
        energy_rms = float(_residual_group_metrics(F_hat, cells)["energy_rms"])
        merit = max(merit, energy_rms)
    return merit


def _energy_dominates_merit(
    energy_F: float,
    merit_F: float,
    *,
    energy_dominate_ratio: float = 0.85,
) -> bool:
    """True when energy residual co-dominates the max-merit floor."""
    return float(energy_F) >= float(energy_dominate_ratio) * max(float(merit_F), 1.0e-30)


def _refresh_aware_joint_acceptable(
    *,
    energy_F: float,
    energy_trial: float,
    merit_F: float,
    merit_post: float,
    split_F: float,
    split_post: float,
    energy_dominate_ratio: float = 0.85,
) -> bool:
    """Refresh-aware joint accept (engineering bridge on Hamel Eq. 2.9).

    When energy dominates the current merit and the trial lowers energy but
    raises pre-refresh max-merit (split climb), accept only if after a hydro
    refresh the post merit and split are not worse than the current point.

    Bare energy-descent accept without this gate was falsified
    (``phase2_extent_raw_energy_desc_accept_ab_o8.json``).
    """
    if not _energy_dominates_merit(
        energy_F, merit_F, energy_dominate_ratio=energy_dominate_ratio
    ):
        return False
    if not (np.isfinite(energy_trial) and float(energy_trial) < float(energy_F)):
        return False
    if not (np.isfinite(merit_post) and float(merit_post) <= float(merit_F) + 1.0e-15):
        return False
    if float(split_post) > float(split_F) + 1.0e-15:
        return False
    return True


def _zero_gas_holdup_newton_step(dx: np.ndarray, cells: List[Cell]) -> np.ndarray:
    """Zero packed gas holdup DOFs (dense/bubble) for split-constrained LS.

    Energy-dominated plateau: freeze phase-split inventory updates and keep
    temperature/solid Newton components. Engineering bridge — not Hamel.
    """
    from src.solvers.nr_indexing import n_gas_var, n_var

    out = np.asarray(dx, dtype=np.float64).copy()
    offset = 0
    for cell in cells:
        nv = int(n_var(cell))
        ng = int(n_gas_var(cell))
        if ng > 0:
            out[offset : offset + ng] = 0.0
        offset += nv
    return out


_SAME_PHASE_SYNGAS_SPECIES: tuple[str, ...] = ("H2", "CO", "CH4")
_SAME_PHASE_OVERLAP_ATOL_MOL_S: float = 1.0e-12  # [mol/s]


def _zero_same_phase_oxidizer_fuel_newton_step(
    dx: np.ndarray,
    cells: List[Cell],
    *,
    atol_mol_s: float = _SAME_PHASE_OVERLAP_ATOL_MOL_S,
    allow_oxidizer_into_fuel: bool = False,
) -> np.ndarray:
    """Zero Newton components that would create same-phase O₂ ∩ syngas overlap.

    Two-phase Startwert keeps bubble O₂ and dense H₂/CO/CH₄ apart (handoff §5.49).
    Putting syngas into an O₂-bearing phase (e.g. bubble CO) then per-phase
    fastox burns it and knocks N_ex off the linearization (handoff §5.50–5.51).

    ``allow_oxidizer_into_fuel=False`` (default): also block O₂ entering a
    fuel-bearing phase. That keeps R1 dark because heterogeneous oxidation
    lives only in the emulsion.

    ``allow_oxidizer_into_fuel=True``: allow ``dN_O2>0`` into fuel. Per-phase
    fastox burns the overlap; leftover O₂ lights R1. Naked R12 without fastox
    still jumps F to 10⁸ — this path requires per-phase fastox.

    Step-manifold restriction on Hamel Eq. 2.9, not a change to R1/R2/K_bd.

    Source: Hamel (1999) Eq. 2.9; handoff §5.51–5.52.
    """
    atol = float(max(atol_mol_s, 0.0))
    out = np.asarray(dx, dtype=np.float64).copy()
    j_o2 = int(GAS_SPECIES_INDEX["O2"])
    syngas_idx = tuple(int(GAS_SPECIES_INDEX[name]) for name in _SAME_PHASE_SYNGAS_SPECIES)
    block_o2_into_fuel = not bool(allow_oxidizer_into_fuel)
    offset = 0
    for cell in cells:
        nv = int(n_var(cell))
        gas_idx = main_nr_gas_idx_for_cell(cell)
        mode = gas_mode_for_nr(cell)
        n_sp = int(gas_idx.size)
        if n_sp > 0 and mode in {"two_phase", "single"}:
            loc = {int(gi): k for k, gi in enumerate(gas_idx)}
            phases: list[tuple[int, np.ndarray]] = [
                (offset, np.asarray(cell.N_d, dtype=np.float64)),
            ]
            if mode == "two_phase":
                phases.append((offset + n_sp, np.asarray(cell.N_b, dtype=np.float64)))
            for phase_offset, N_phase in phases:
                fuel = float(sum(max(float(N_phase[j]), 0.0) for j in syngas_idx))
                o2 = float(max(float(N_phase[j_o2]), 0.0))
                if block_o2_into_fuel and j_o2 in loc and fuel > atol:
                    i = phase_offset + loc[j_o2]
                    if float(out[i]) > 0.0:
                        out[i] = 0.0
                if o2 > atol:
                    for j in syngas_idx:
                        if j not in loc:
                            continue
                        i = phase_offset + loc[j]
                        if float(out[i]) > 0.0:
                            out[i] = 0.0
        offset += nv
    return out


def _line_search_trial_acceptable(
    *,
    merit_trial: float,
    merit_F: float,
    split_trial: float,
    split_F: float,
    energy_F: float,
    use_gas_phase_split: bool,
    use_energy: bool,
    accepted_nonmonotone: bool = False,
    energy_dominate_ratio: float = 0.85,
    min_merit_rel_drop_for_split_rise: float = 0.05,
) -> bool:
    """Accept LS trial; block split inflation into the energy merit floor.

    Pathology (post-F3 O3): energy ~0.146 sets merit; tiny energy chips are
    accepted while phase-split climbs into co-dominance (0.098→0.146).  Allow
    mild split rises that stay below ``ratio * merit_trial``; require a
    meaningful relative merit drop only when the trial would make split
    co-dominant.

    Engineering bridge on Hamel Eq. 2.9 line search — not a thesis equation.
    """
    if bool(accepted_nonmonotone):
        if (
            use_gas_phase_split
            and use_energy
            and float(energy_F) >= float(energy_dominate_ratio) * max(float(merit_F), 1.0e-30)
            and float(split_trial) > float(split_F) + 1.0e-15
            and float(split_trial)
            >= float(energy_dominate_ratio) * max(float(merit_trial), 1.0e-30)
        ):
            return False
        return True
    if not (np.isfinite(merit_trial) and float(merit_trial) < float(merit_F)):
        return False
    if not (use_gas_phase_split and use_energy):
        return True
    if float(split_trial) <= float(split_F) + 1.0e-15:
        return True
    merit_floor = max(float(merit_F), 1.0e-30)
    if float(energy_F) < float(energy_dominate_ratio) * merit_floor:
        return True
    # Mild split rise still below co-dominance on the trial merit is OK.
    if float(split_trial) < float(energy_dominate_ratio) * max(float(merit_trial), 1.0e-30):
        return True
    rel_drop = (float(merit_F) - float(merit_trial)) / merit_floor
    return rel_drop >= float(min_merit_rel_drop_for_split_rise)


def update_t_cap_split_hysteresis(
    hold: bool,
    *,
    split_rms: float,
    merit_F: float,
    merit_without_energy: float | None = None,
    engage_ratio: float = 0.85,
    release_ratio: float = 0.70,
) -> bool:
    """Latch split-dominant T-cap across NR iters / outers (F3).

    Engage when ``split >= engage_ratio * merit_F``.  Release only when split
    is clearly secondary versus non-energy merit (``max(rms, split)``), so a
    transient energy spike that inflates ``merit_F`` cannot unlock 600 K steps.

    Engineering bridge around Hamel Eq. 2.9 damping — not a thesis equation.
    """
    split = float(split_rms)
    merit = max(float(merit_F), 1.0e-30)
    if split > 1.0e-12 and split >= float(engage_ratio) * merit:
        return True
    if not bool(hold):
        return False
    gate = float(merit_without_energy) if merit_without_energy is not None else merit
    gate = max(gate, 1.0e-30)
    # Stay latched while split remains co-dominant vs non-energy residuals.
    if split >= float(release_ratio) * gate:
        return True
    return False


def resolve_inner_t_max_step_K(
    *,
    ene_norm: float,
    tol_rms: float,
    split_rms: float,
    merit_F: float,
    dominant_non_energy: float,
    default_cap_K: float = 600.0,
    energy_deep_cap_K: float = 1.0,
    split_dominant_cap_K: float | None = 150.0,
    split_dominant_merit_ratio: float = 0.85,
    line_search_split_enabled: bool = False,
    hold_split_cap: bool = False,
) -> float:
    """Adaptive per-iteration temperature step cap for inner NR.

    Symmetric to the energy-deep branch: when phase-split merit co-dominates
    (``split_rms >= ratio * merit_F``), limit ``|ΔT|`` so Newton does not use
    temperature to chase split/solid residuals.  ``merit_F`` is already
    ``max(rms, split, energy)``; requiring ``split >= energy`` was too strict
    (outer-8 v8: split/merit≈0.86 with energy slightly above split).

    ``hold_split_cap`` (F3 hysteresis): keep the split-dominant cap even when
    a transient energy spike drops ``split/merit`` below the engage ratio.

    Ref: ``global_nr_iteration.execute_global_nr_solve``; Hamel Eq. 2.9 damped NR.
    """
    _energy_frac = float(ene_norm) / max(float(dominant_non_energy), 1.0e-30)
    if float(ene_norm) < float(tol_rms) * 1.0e-3 and _energy_frac < 1.0e-3:
        return float(energy_deep_cap_K)
    cap = float(split_dominant_cap_K) if split_dominant_cap_K is not None else 0.0
    if line_search_split_enabled and cap > 0.0:
        split_now = float(split_rms) > 1.0e-12 and float(split_rms) >= float(
            split_dominant_merit_ratio
        ) * float(merit_F)
        if bool(hold_split_cap) or split_now:
            return float(cap)
    return float(default_cap_K)


def resolve_t_step_with_abs_dT_budget(
    t_max_adaptive_K: float,
    *,
    dT_so_far_K: float,
    abs_dT_budget_K: float | None,
) -> float:
    """Clamp per-step |ΔT| by remaining absolute temperature budget for this inner solve.

    Prevents multi-iter accumulation (e.g. 4×150 K) from producing outer-level
    |ΔT|~600 K when a warmer bed0 seed unlocks large Newton steps.

    Engineering trust-region on Hamel Eq. 2.9 — not a thesis equation.
    """
    step = float(max(t_max_adaptive_K, 0.0))
    budget = float(abs_dT_budget_K) if abs_dT_budget_K is not None else 0.0
    if budget <= 0.0:
        return step
    remaining = budget - float(max(dT_so_far_K, 0.0))
    if remaining <= 1.0e-6:
        return 1.0e-6
    return float(min(step, remaining))


def _inner_nr_check1_satisfied(
    *,
    F_hat: np.ndarray,
    cells: List[Cell],
    norm_F: float,
    x_rms: float,
    tol_rms: float,
    inner_convergence: InnerConvergence | None,
    enforce_gas_phase_split_convergence: bool,
    line_search_merit_enabled: bool,
    line_search_energy_merit_enabled: bool = False,
    tol_gas_phase_split_max: float | None,
) -> bool:
    """Check1 gate aligned with thesis line-search merit when split merit is enabled."""
    group_metrics = _residual_group_metrics(F_hat, cells)
    component_norm = max(
        group_metrics["gas_rms"],
        group_metrics["solid_rms"],
        group_metrics["energy_rms"],
    )
    split_tol = float(tol_rms if tol_gas_phase_split_max is None else tol_gas_phase_split_max)
    use_merit_check1 = bool(line_search_merit_enabled and enforce_gas_phase_split_convergence)
    if use_merit_check1:
        merit = _line_search_merit(
            F_hat,
            cells,
            use_gas_phase_split=True,
            use_energy=bool(line_search_energy_merit_enabled),
        )
        residual_ok = (
            inner_convergence.is_converged(merit, x_rms)
            if inner_convergence is not None
            else merit < tol_rms
        )
        return bool(residual_ok and component_norm <= tol_rms)

    gas_phase_metrics = _residual_gas_phase_metrics(F_hat, cells)
    gas_split_max_abs = float(gas_phase_metrics["gas_phase_split_max_abs"])
    split_ok = (not enforce_gas_phase_split_convergence) or (gas_split_max_abs <= split_tol)
    component_max_abs = max(
        group_metrics["gas_max_abs"],
        group_metrics["solid_max_abs"],
        group_metrics["energy_max_abs"],
    )
    if inner_convergence is not None:
        return bool(
            inner_convergence.is_converged(norm_F, x_rms)
            and component_norm <= tol_rms
            and component_max_abs <= tol_rms
            and split_ok
        )
    return bool(
        norm_F < tol_rms
        and component_norm <= tol_rms
        and component_max_abs <= tol_rms
        and split_ok
    )


def compute_live_nr_scaled_metrics(
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    *,
    ref_gas_mol_s: float,
    ref_solid_kg_s: float,
    ref_energy_W: float,
    use_gas_phase_split_merit: bool,
    use_energy_merit: bool = False,
) -> dict[str, float]:
    """Recompute scaled NR metrics from the current packed cell state."""
    apply_bc_fn()
    F = global_residual(pack_reactor(cells), cells, apply_bc_fn)
    scale = build_equation_scales(cells, ref_gas_mol_s, ref_solid_kg_s, ref_energy_W)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    group_metrics = _residual_group_metrics(F_hat, cells)
    gas_phase_metrics = _residual_gas_phase_metrics(F_hat, cells)
    merit = _line_search_merit(
        F_hat,
        cells,
        use_gas_phase_split=use_gas_phase_split_merit,
        use_energy=use_energy_merit,
    )
    return {
        "rms_scaled_final": float(merit),
        "rms_scaled_gas_final": float(group_metrics["gas_rms"]),
        "rms_scaled_solid_final": float(group_metrics["solid_rms"]),
        "rms_scaled_energy_final": float(group_metrics["energy_rms"]),
        "rms_scaled_gas_combined_final": float(gas_phase_metrics["gas_combined_rms"]),
        "rms_scaled_gas_phase_split_final": float(gas_phase_metrics["gas_phase_split_rms"]),
        "rms_scaled_component_max_final": float(
            max(group_metrics["gas_rms"], group_metrics["solid_rms"], group_metrics["energy_rms"])
        ),
        "max_abs_scaled_gas_final": float(group_metrics["gas_max_abs"]),
        "max_abs_scaled_solid_final": float(group_metrics["solid_max_abs"]),
        "max_abs_scaled_energy_final": float(group_metrics["energy_max_abs"]),
        "max_abs_scaled_gas_combined_final": float(gas_phase_metrics["gas_combined_max_abs"]),
        "max_abs_scaled_gas_phase_split_final": float(gas_phase_metrics["gas_phase_split_max_abs"]),
        "max_abs_scaled_final": float(
            max(
                group_metrics["gas_max_abs"],
                group_metrics["solid_max_abs"],
                group_metrics["energy_max_abs"],
            )
        ),
        "residual": float(np.linalg.norm(F)),
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
    solid_ftb_relative_floor: bool = False,
) -> np.ndarray:
    """严格限制 Newton 步长，防止物理量跳变导致发散。

    对物理上"空"的 holdup_transport 单元（m_solid ≈ 0 且无任何固相进流），
    仅当 zero_empty_solid=True 时才将固相 dx 清零（GD fallback 专用）。
    NR 方向不应清零：Newton 步由 Jacobian 计算，天然保证 descent，
    强制清零会使最大残差 DOF 无法更新，导致线搜索反下降。

    ``solid_ftb_relative_floor``
        True 时固相 FTB 地板与气相一致：``max(1e-6, 1e-4·max m)``。
        微量库存（~1e-6 kg）负步改为局部清零，不进全局 α——避免掐死上段
        dense-CO 等气相 Newton 步（handoff §5.6）。
    """
    rg_limit = 0.2 * max(float(ref_gas_mol_s), 1.0)

    out = np.asarray(dx, dtype=np.float64).copy()

    # Holdups share one global α (trust + significant-inventory FTB). Temperature
    # is clipped component-wise and excluded from α — raw |dT|~1e5–1e6 K would
    # otherwise collapse the whole Newton step. Near-empty solid inventories may
    # be excluded from α via ``solid_ftb_relative_floor``. CO wipe false-accepts
    # are blocked by the LS collapse guard.
    alpha = 1.0
    offset = 0
    gas_state_max = 0.0
    gas_entries: list[tuple[int, float]] = []
    for cell in cells:
        nv = n_var(cell)
        ng = n_gas_var(cell)
        nv_sol = n_solid_var(cell)
        rs_limit = 0.2 * max(_solid_reference_for_cell(cell, ref_solid_kg_s), 0.1)

        if ng > 0:
            gas_idx = main_nr_gas_idx_for_cell(cell)
            gas_mode = gas_mode_for_nr(cell)
            if gas_mode == "two_phase":
                state = np.concatenate(
                    [
                        np.maximum(np.asarray(cell.N_d[gas_idx], dtype=np.float64), 0.0),
                        np.maximum(np.asarray(cell.N_b[gas_idx], dtype=np.float64), 0.0),
                    ]
                )
            elif gas_mode == "single":
                state = np.maximum(np.asarray(cell.N_d[gas_idx], dtype=np.float64), 0.0)
            else:
                state = np.zeros(0, dtype=np.float64)
            if state.size == ng:
                gas_state_max = max(gas_state_max, float(np.max(state)))
                snap_o2 = bool(getattr(cell, "_nr_snap_trace_o2", False))
                o2_atol = float(getattr(cell, "_nr_snap_trace_o2_atol_mol_s", 1.0e-12))
                j_o2 = int(GAS_SPECIES_INDEX["O2"])
                n_sp = int(gas_idx.size)
                nin_o2 = float(
                    cell.N_d_in[j_o2]
                    + cell.N_b_in[j_o2]
                    + cell.N_zu_d[j_o2]
                    + cell.N_zu_b[j_o2]
                )
                for i in range(ng):
                    gas_entries.append((offset + i, float(state[i])))
                    dxi = float(out[offset + i])
                    # 贫氧粘滞：holdup 与入口均低于阈值时禁止 Newton 引入微量 O₂
                    if snap_o2 and nin_o2 <= o2_atol and n_sp > 0:
                        gi = int(gas_idx[i if i < n_sp else i - n_sp])
                        if gi == j_o2 and float(state[i]) <= o2_atol:
                            out[offset + i] = 0.0
                            dxi = 0.0
                    if abs(dxi) > rg_limit > 0.0:
                        alpha = min(alpha, rg_limit / abs(dxi))

        sol_start = offset + solid_offset(cell)
        if nv_sol > 0:
            solid_state = cell.m_solid.reshape(-1)[solid_flat_indices_for_nr(cell)]
            s_max = float(np.max(solid_state)) if solid_state.size else 0.0
            if solid_ftb_relative_floor:
                # Align with gas x_floor: max(1e-3-ish floor scale, 1e-4·max).
                s_floor = max(1.0e-6, 1.0e-4 * max(s_max, 0.0))
            else:
                s_floor = max(1.0e-6, 1.0e-8 * max(s_max, 0.0)) if solid_state.size else 1.0e-6
            for j in range(nv_sol):
                idx = sol_start + j
                dxi = float(out[idx])
                if _is_solid_empty_cell(cell, ref_solid_kg_s):
                    if zero_empty_solid:
                        if dxi > 0.0:
                            out[idx] = 0.0
                    else:
                        tight_limit = 1e-3 * max(
                            _solid_reference_for_cell(cell, ref_solid_kg_s), 1.0
                        )
                        if abs(dxi) > tight_limit > 0.0:
                            alpha = min(alpha, tight_limit / abs(dxi))
                else:
                    if abs(dxi) > rs_limit > 0.0:
                        alpha = min(alpha, rs_limit / abs(dxi))
                xi = float(max(solid_state[j], 0.0)) if j < solid_state.size else 0.0
                if dxi < 0.0 and xi > s_floor:
                    alpha = min(alpha, 0.8 * xi / (-dxi))
                elif dxi < 0.0:
                    out[idx] = 0.0

        if t_step_limit_K is not None and has_temperature_var(cell):
            temp_idx = offset + temperature_offset(cell)
            limit = float(max(t_step_limit_K, 1e-12))
            out[temp_idx] = float(np.clip(out[temp_idx], -limit, limit))
        offset += nv

    x_floor = max(1.0e-3, 1.0e-4 * gas_state_max)
    for idx, xi in gas_entries:
        dxi = float(out[idx])
        if dxi >= 0.0:
            continue
        if xi <= x_floor:
            out[idx] = 0.0
            continue
        if (-dxi) <= 1.0e-12 * max(xi, 1.0):
            continue
        alpha = min(alpha, 0.8 * xi / (-dxi))

    if alpha < 1.0:
        offset = 0
        for cell in cells:
            nv = n_var(cell)
            if has_temperature_var(cell):
                temp_idx = offset + temperature_offset(cell)
                t_step = float(out[temp_idx])
                out[offset : offset + nv] *= alpha
                out[temp_idx] = t_step
            else:
                out[offset : offset + nv] *= alpha
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


def diagnose_newton_step_consistency(
    J: sp.spmatrix | np.ndarray,
    F_hat: np.ndarray,
    dx: np.ndarray,
    dx_clipped: np.ndarray | None = None,
) -> dict[str, float]:
    """P0 probe: relative linear defect of Newton / clipped Newton directions.

    For a consistent linearisation, ``||J dx + F_hat|| / ||F_hat||`` should be
    near machine/solver tolerance. Large growth after clipping means the accepted
    search direction has lost the local descent property of ``Jd=-F``.
    """
    F = np.asarray(F_hat, dtype=np.float64).reshape(-1)
    d = np.asarray(dx, dtype=np.float64).reshape(-1)
    if F.size == 0 or d.size == 0 or F.size != d.size:
        return {
            "f_hat_norm": 0.0,
            "linear_defect_rel": float("nan"),
            "dx_norm": 0.0,
        }
    f_norm = float(np.linalg.norm(F))
    denom = max(f_norm, 1.0e-30)
    lin = np.asarray(J @ d, dtype=np.float64).reshape(-1) + F
    out: dict[str, float] = {
        "f_hat_norm": f_norm,
        "linear_defect_rel": float(np.linalg.norm(lin) / denom),
        "dx_norm": float(np.linalg.norm(d)),
    }
    if dx_clipped is not None:
        dc = np.asarray(dx_clipped, dtype=np.float64).reshape(-1)
        lin_c = np.asarray(J @ dc, dtype=np.float64).reshape(-1) + F
        d_norm = max(float(np.linalg.norm(d)), 1.0e-30)
        out["clip_linear_defect_rel"] = float(np.linalg.norm(lin_c) / denom)
        out["clip_displacement_rel"] = float(np.linalg.norm(dc - d) / d_norm)
        out["clipped_components"] = float(np.count_nonzero(np.abs(dc - d) > 1e-12))
    return out


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
            t_dof_indices=t_dof_indices,
            t_max_step_K=float(t_max_step_K),
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
    enforce_gas_phase_split_convergence: bool = False,
    tol_gas_phase_split_max: float | None = None,
    line_search_use_gas_phase_split_merit: bool = False,
    line_search_use_energy_merit: bool = False,
    refresh_hydrodynamics_on_accepted_step: bool = False,
    line_search_fast_oxidation_projection: bool = False,
    line_search_per_phase_fastox: bool = False,
    line_search_per_phase_fastox_bed_limit: int = 0,
    clip_same_phase_oxidizer_fuel_step: bool = False,
    allow_oxidizer_into_fuel_step: bool = False,
    line_search_syngas_collapse_guard: bool = True,
    line_search_syngas_collapse_ratio: float = 0.05,
    line_search_syngas_floor_ratio: float = 0.0,
    line_search_refresh_aware_accept: bool = False,
    line_search_split_constrained_step: bool = False,
    clip_solid_ftb_relative_floor: bool = False,
    inner_t_step_default_cap_K: float = 600.0,
    inner_split_dominant_t_step_cap_K: float | None = 150.0,
    inner_split_dominant_merit_ratio: float = 0.85,
    inner_t_cap_hold_split_dominant: bool = False,
    inner_t_cap_split_release_ratio: float = 0.70,
    inner_abs_dT_budget_K: float | None = None,
    connectivity_graph: Any | None = None,
    record_inner_group_history: bool = False,
    verbose: bool = False,
) -> dict:
    """阻尼 Newton–Raphson：ĴΔx = −F̂（F̂ = F / 静态 eq_scale），再线搜索阻尼。

    Hamel Eq. 2.9（gedämpftes Newton-Verfahren）的工程实现。

    inner_convergence
        若给定，Check1 使用 ``atol + rtol * ‖x‖_RMS``（‖x‖ 为打包状态向量 RMS）；
        若为 ``None``，则沿用 ``norm_F < tol_rms``（与历史行为一致）。
    """
    from src.solvers.global_nr_iteration import execute_global_nr_solve

    return execute_global_nr_solve(
        cells=cells,
        apply_bc_fn=apply_bc_fn,
        apply_local_bc_fn=apply_local_bc_fn,
        affected_residual_cells_fn=affected_residual_cells_fn,
        ref_gas_mol_s=ref_gas_mol_s,
        ref_solid_kg_s=ref_solid_kg_s,
        ref_energy_W=ref_energy_W,
        max_iter=max_iter,
        tol_rms=tol_rms,
        inner_convergence=inner_convergence,
        lambda_init=lambda_init,
        n_damp_halvings=n_damp_halvings,
        lambda_min=lambda_min,
        jacobian_strategy=jacobian_strategy,
        linear_solver_backend=linear_solver_backend,
        jacobian_lag_steps=jacobian_lag_steps,
        allow_gd_fallback=allow_gd_fallback,
        prefer_full_step=prefer_full_step,
        step_model=step_model,
        lm_mu0=lm_mu0,
        lm_mu_growth=lm_mu_growth,
        ptc_alpha0=ptc_alpha0,
        ptc_alpha_growth=ptc_alpha_growth,
        equil_iters=equil_iters,
        equil_scale_clip=equil_scale_clip,
        line_search_max_trials=line_search_max_trials,
        nonmonotone_enabled=nonmonotone_enabled,
        nonmonotone_window=nonmonotone_window,
        nonmonotone_relax=nonmonotone_relax,
        enforce_gas_phase_split_convergence=enforce_gas_phase_split_convergence,
        tol_gas_phase_split_max=tol_gas_phase_split_max,
        line_search_use_gas_phase_split_merit=line_search_use_gas_phase_split_merit,
        line_search_use_energy_merit=line_search_use_energy_merit,
        refresh_hydrodynamics_on_accepted_step=refresh_hydrodynamics_on_accepted_step,
        line_search_fast_oxidation_projection=line_search_fast_oxidation_projection,
        line_search_per_phase_fastox=line_search_per_phase_fastox,
        line_search_per_phase_fastox_bed_limit=line_search_per_phase_fastox_bed_limit,
        clip_same_phase_oxidizer_fuel_step=clip_same_phase_oxidizer_fuel_step,
        allow_oxidizer_into_fuel_step=allow_oxidizer_into_fuel_step,
        line_search_syngas_collapse_guard=line_search_syngas_collapse_guard,
        line_search_syngas_collapse_ratio=line_search_syngas_collapse_ratio,
        line_search_syngas_floor_ratio=line_search_syngas_floor_ratio,
        line_search_refresh_aware_accept=line_search_refresh_aware_accept,
        line_search_split_constrained_step=line_search_split_constrained_step,
        clip_solid_ftb_relative_floor=clip_solid_ftb_relative_floor,
        inner_t_step_default_cap_K=inner_t_step_default_cap_K,
        inner_split_dominant_t_step_cap_K=inner_split_dominant_t_step_cap_K,
        inner_split_dominant_merit_ratio=inner_split_dominant_merit_ratio,
        inner_t_cap_hold_split_dominant=inner_t_cap_hold_split_dominant,
        inner_t_cap_split_release_ratio=inner_t_cap_split_release_ratio,
        inner_abs_dT_budget_K=inner_abs_dT_budget_K,
        connectivity_graph=connectivity_graph,
        record_inner_group_history=record_inner_group_history,
        verbose=verbose,
    )
