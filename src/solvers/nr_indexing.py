"""全局 NR 未知量布局与 pack/unpack 契约（与 Hamel 全局联立思路对齐）。

主 NR 向量默认排除 H2S（Gibbs 路径）；thesis 模式下 NH3/TAR1/TAR2 纳入 NR，
与干燥/热解 dense 源项闭合。Legacy 模式仍可将 NH3/TAR 排除在 NR 外。

所有索引与 ``global_nr_solver.solve_global_nr`` 使用的布局必须与本模块保持一致。
"""

from __future__ import annotations

from typing import Any, List

import numpy as np

from src.core.cell import Cell, N_SOLID_COMP, S_ASH, S_CHAR
from src.core.species import GAS_SPECIES_INDEX, N_GAS

# Legacy NR layout (pre-C-section): pyrolysis species fixed at vorab seed.
LEGACY_MAIN_NR_EXCLUDED_SPECIES: tuple[str, ...] = ("H2S", "NH3", "TAR1", "TAR2")
# Hamel C-section closure: only H2S stays Gibbs-only outside the NR Jacobian.
THESIS_MAIN_NR_EXCLUDED_SPECIES: tuple[str, ...] = ("H2S",)

PYROLYSIS_NR_GAS_SPECIES: tuple[str, ...] = ("NH3", "TAR1", "TAR2")
MAJOR_NR_GAS_SPECIES: tuple[str, ...] = ("CO", "CO2", "H2", "H2O", "CH4", "O2", "N2")

# Backward-compatible module constants (legacy 7-species layout).
MAIN_NR_EXCLUDED_SPECIES: tuple[str, ...] = LEGACY_MAIN_NR_EXCLUDED_SPECIES
_MAIN_EXCLUDED_IDX = {GAS_SPECIES_INDEX[s] for s in MAIN_NR_EXCLUDED_SPECIES}

MAIN_NR_GAS_IDX = np.array(
    [i for i in range(N_GAS) if i not in _MAIN_EXCLUDED_IDX],
    dtype=np.int64,
)
MAIN_NR_GAS_COUNT = int(MAIN_NR_GAS_IDX.size)


def _excluded_from_config(config: Any | None) -> tuple[str, ...]:
    if config is None:
        return LEGACY_MAIN_NR_EXCLUDED_SPECIES
    include_pyro = bool(getattr(config, "nr_include_pyrolysis_gas_species_in_unknowns_thesis", False))
    if include_pyro:
        return THESIS_MAIN_NR_EXCLUDED_SPECIES
    return LEGACY_MAIN_NR_EXCLUDED_SPECIES


def main_nr_excluded_species(config: Any | None = None) -> tuple[str, ...]:
    """Species excluded from the global NR unknown/residual vector."""
    return _excluded_from_config(config)


def build_main_nr_gas_idx(excluded: tuple[str, ...]) -> np.ndarray:
    excluded_idx = {GAS_SPECIES_INDEX[s] for s in excluded}
    return np.array([i for i in range(N_GAS) if i not in excluded_idx], dtype=np.int64)


def main_nr_gas_idx_for_config(config: Any | None = None) -> np.ndarray:
    return build_main_nr_gas_idx(main_nr_excluded_species(config))


def main_nr_gas_count_for_config(config: Any | None = None) -> int:
    return int(main_nr_gas_idx_for_config(config).size)


def main_nr_gas_idx_for_cell(cell: Cell) -> np.ndarray:
    custom = getattr(cell, "_main_nr_gas_idx", None)
    if custom is not None:
        return np.asarray(custom, dtype=np.int64)
    return MAIN_NR_GAS_IDX


def main_nr_gas_count_for_cell(cell: Cell) -> int:
    return int(main_nr_gas_idx_for_cell(cell).size)


def main_nr_gas_phase_split_positions_for_cell(cell: Cell) -> np.ndarray:
    """Packed two-phase block positions used for dense/bubble split merit/convergence.

    Pyrolysis species (NH3/TAR*) stay in combined gas NR unknowns but are excluded
    from phase-split gates because their direct source is dense-only.
    """
    gas_idx = main_nr_gas_idx_for_cell(cell)
    major_set = {int(GAS_SPECIES_INDEX[s]) for s in MAJOR_NR_GAS_SPECIES}
    positions = [i for i, gi in enumerate(gas_idx) if int(gi) in major_set]
    return np.asarray(positions, dtype=np.int64)


def nr_projection_gas_species(config: Any | None = None) -> tuple[str, ...]:
    """Bottom-joint / start-value LSQ species aligned with NR gas unknowns (excl. H2S)."""
    excluded = set(main_nr_excluded_species(config))
    out: list[str] = []
    for sp in (*MAJOR_NR_GAS_SPECIES, *PYROLYSIS_NR_GAS_SPECIES):
        if sp not in excluded and sp not in out:
            out.append(sp)
    return tuple(out)


def nr_dense_only_pyrolysis_support_species(config: Any | None = None) -> tuple[str, ...]:
    """Pyrolysis species kept as fixed dense support when not NR unknowns."""
    if bool(getattr(config, "nr_include_pyrolysis_gas_species_in_unknowns_thesis", False)):
        return ()
    return PYROLYSIS_NR_GAS_SPECIES


def install_main_nr_gas_indexing_on_cells(cells: List[Cell], config: Any | None) -> None:
    """Attach per-cell NR gas index map and NH3 Gibbs policy from reactor config."""
    gas_idx = main_nr_gas_idx_for_config(config)
    nh3_idx = int(GAS_SPECIES_INDEX["NH3"])
    nh3_in_nr = nh3_idx in set(int(i) for i in gas_idx)
    for cell in cells:
        cell._main_nr_gas_idx = gas_idx.copy()
        cell._nr_solves_nh3_transport = nh3_in_nr


def install_nr_jacobian_metadata_on_cells(cells: List[Cell], config: Any | None) -> None:
    """Attach per-cell metadata consumed by ``jacobian_analytic`` band overlays."""
    lower = float(getattr(config, "solid_lower_inlet_frac", 0.0)) if config is not None else 0.0
    lower = float(np.clip(lower, 0.0, 1.0))
    for cell in cells:
        if cell_kind(cell) == "bed":
            cell._nr_solid_lower_inlet_frac = lower


def cell_kind(cell: Cell) -> str:
    """Normalized NR layout kind for one solver cell."""
    return str(getattr(cell, "cell_type", "bed")).strip().lower() or "bed"


def gas_mode_for_nr(cell: Cell) -> str:
    """Gas unknown layout for the Hamel-aligned NR vector."""
    kind = cell_kind(cell)
    if kind in {"cyclone", "return_leg"}:
        return "none"
    if kind == "freeboard":
        return "single"
    return "two_phase"


def n_gas_var(cell: Cell) -> int:
    """Number of gas unknowns included in the global-NR vector for one cell."""
    count = main_nr_gas_count_for_cell(cell)
    mode = gas_mode_for_nr(cell)
    if mode == "none":
        return 0
    if mode == "single":
        return count
    return 2 * count


def solid_offset(cell: Cell) -> int:
    """Start offset of solid unknowns inside one packed cell block."""
    return n_gas_var(cell)


def temperature_offset(cell: Cell) -> int:
    """Offset of the cell temperature unknown inside one packed cell block."""
    if not has_temperature_var(cell):
        raise ValueError(f"{cell_kind(cell)} has no NR temperature unknown")
    return n_gas_var(cell) + n_solid_var(cell)


def has_temperature_var(cell: Cell) -> bool:
    """Whether this cell carries an independent NR temperature unknown.

    Cyclone / return-leg side elements carry a pass-through energy closure row
    (``Cell._calc_side_temperature_closure_residual``) in the global NR system.
    Ref: Hamel (1999) §2.3 Nebenelemente; ``docs/hamel_engineering_bridge_notes.md`` §5.
    """
    return cell_kind(cell) != "connection_pipe"


def n_var(cell: Cell) -> int:
    """Number of NR unknowns for this cell's physical layout."""
    return n_gas_var(cell) + n_solid_var(cell) + (1 if has_temperature_var(cell) else 0)


def solid_component_indices_for_nr(cell: Cell) -> np.ndarray:
    """Solid component indices kept in global-NR state vector for this cell."""
    if cell_kind(cell) == "freeboard":
        # Hamel freeboard alignment:
        # freeboard solids are computed in Vorabrechnung via analytical
        # trajectory closure and should not enter the inner global-NR Jacobian as
        # explicit Eulerian holdup unknowns. Only gas composition and T remain in
        # the freeboard cell state vector for Newton updates.
        return np.zeros(0, dtype=np.int64)
    # Global NR should only carry transported/resident solid holdups that remain
    # physical state variables throughout the bed. VM/moisture are handled as
    # fresh-feed release budgets / source supports in the lower cells; once they
    # lose inlet support, BC projection hard-resets them to zero. Keeping them in
    # the packed Newton vector therefore creates dead DOFs and Jacobian zero
    # rows/cols without adding physical information.
    return np.array([S_CHAR, S_ASH], dtype=np.int64)


def solid_flat_indices_for_nr(cell: Cell) -> np.ndarray:
    """Flattened m_solid indices included in global-NR vector for this cell."""
    comp_idx = solid_component_indices_for_nr(cell)
    nk = cell.solid.n_size_classes
    if nk <= 0:
        return np.zeros(0, dtype=np.int64)
    base = np.arange(nk, dtype=np.int64)[:, None] * int(N_SOLID_COMP)
    return (base + comp_idx[None, :]).reshape(-1)


def n_solid_var(cell: Cell) -> int:
    """Number of solid unknowns included in global-NR vector for one cell."""
    return int(cell.solid.n_size_classes) * int(solid_component_indices_for_nr(cell).size)


UNPACK_MODE_PRODUCTION_CLIP: str = "production_clip"
UNPACK_MODE_STRICT_DOMAIN_NO_PROJECTION: str = "strict_domain_no_projection"
# T>0 是闭定义域的开下界；不是生产 unpack 的 300 K fence。
CLOSED_DOMAIN_T_OPEN_LOWER_K: float = 0.0  # [K]
_NTOT_DENOM_FLOOR_MOL_S: float = 1.0e-30  # [mol/s]


class DomainViolation(ValueError):
    """打包状态不在闭定义域：拒绝写入，不替换为裁剪后的状态。"""


def _clip_production_temperature_K(T: float, cell: Cell) -> float:
    """生产 unpack 的温度 fence。原文无此式。"""
    t_min = float(getattr(cell, "_nr_temperature_min_K", 300.0))
    t_max = float(getattr(cell, "_nr_temperature_max_K", 2500.0))
    if not np.isfinite(t_min):
        t_min = 300.0
    if not np.isfinite(t_max):
        t_max = 2500.0
    t_min = max(t_min, 300.0)
    t_max = min(max(t_max, t_min + 1.0), 2500.0)
    if T < t_min:
        return t_min
    if T > t_max:
        return t_max
    return T


def _write_packed_into_cell(
    x: np.ndarray,
    cell: Cell,
    *,
    clip_nonnegative: bool,
    clip_temperature: bool,
) -> None:
    """把打包块写入 cell。``clip_*`` 仅生产路径为 True。"""
    gas_idx = main_nr_gas_idx_for_cell(cell)
    gas_count = int(gas_idx.size)
    gas_mode = gas_mode_for_nr(cell)
    nv_sol = n_solid_var(cell)
    offset = 0
    if gas_mode == "two_phase":
        nd = np.asarray(x[offset : offset + gas_count], dtype=np.float64).copy()
        offset += gas_count
        nb = np.asarray(x[offset : offset + gas_count], dtype=np.float64).copy()
        offset += gas_count
        if clip_nonnegative:
            nd = np.maximum(nd, 0.0)
            nb = np.maximum(nb, 0.0)
        cell.N_d[gas_idx] = nd
        cell.N_b[gas_idx] = nb
    elif gas_mode == "single":
        nd = np.asarray(x[offset : offset + gas_count], dtype=np.float64).copy()
        offset += gas_count
        if clip_nonnegative:
            nd = np.maximum(nd, 0.0)
        cell.N_d[gas_idx] = nd
        cell.N_b[gas_idx] = 0.0
    else:
        cell.N_d[gas_idx] = 0.0
        cell.N_b[gas_idx] = 0.0
    if nv_sol > 0:
        m_flat = np.asarray(x[offset : offset + nv_sol], dtype=np.float64).copy()
        if clip_nonnegative:
            m_flat = np.maximum(m_flat, 0.0)
        solid_flat = cell.m_solid.reshape(-1)
        solid_flat[solid_flat_indices_for_nr(cell)] = m_flat
        offset += nv_sol
    if has_temperature_var(cell):
        T = float(x[-1])  # [K]
        cell.T = _clip_production_temperature_K(T, cell) if clip_temperature else T


def packed_nonnegative_slice(cell: Cell) -> slice:
    """打包向量中必须 ``≥0`` 的前缀（气体+固体）；温度不在此切片。"""
    n_nn = n_gas_var(cell) + n_solid_var(cell)
    return slice(0, n_nn)


def packed_cell_in_closed_domain(x: np.ndarray, cell: Cell) -> tuple[bool, str]:
    """闭定义域：有限、物种/库存 ``≥0``（允许零）、``T>0``。不写入 cell。

    Source: Wirsum (1997) Eq. 2.21–2.23 无 clip；N≥0 是定义域不是互补方程。
    """
    z = np.asarray(x, dtype=np.float64)
    nv = n_var(cell)
    if z.shape != (nv,):
        return False, f"packed length {z.size} != n_var={nv}"
    if not np.all(np.isfinite(z)):
        return False, "packed state is not finite"
    nn = packed_nonnegative_slice(cell)
    if np.any(z[nn] < 0.0):
        return False, "gas or solid packed entry is negative"
    if has_temperature_var(cell):
        T = float(z[-1])  # [K]
        if T <= CLOSED_DOMAIN_T_OPEN_LOWER_K:
            return False, "T must be > 0 K"
    return True, ""


def packed_cell_safe_to_evaluate_F(x: np.ndarray, cell: Cell) -> tuple[bool, str]:
    """评价 F 前还要求各存在相的总摩尔流可作分母。零组分本身合法。"""
    ok, why = packed_cell_in_closed_domain(x, cell)
    if not ok:
        return False, why
    gas_mode = gas_mode_for_nr(cell)
    gas_count = main_nr_gas_count_for_cell(cell)
    if gas_mode == "two_phase":
        nd = np.asarray(x[:gas_count], dtype=np.float64)
        nb = np.asarray(x[gas_count : 2 * gas_count], dtype=np.float64)
        if float(np.sum(nd)) < _NTOT_DENOM_FLOOR_MOL_S:
            return False, "dense ntot too small for mole fractions"
        if float(np.sum(nb)) < _NTOT_DENOM_FLOOR_MOL_S:
            return False, "bubble ntot too small for mole fractions"
    elif gas_mode == "single":
        nd = np.asarray(x[:gas_count], dtype=np.float64)
        if float(np.sum(nd)) < _NTOT_DENOM_FLOOR_MOL_S:
            return False, "single-phase ntot too small for mole fractions"
    return True, ""


def packed_reactor_in_closed_domain(x: np.ndarray, cells: List[Cell]) -> tuple[bool, str]:
    z = np.asarray(x, dtype=np.float64)
    expected = int(sum(n_var(c) for c in cells))
    if z.size != expected:
        return False, f"reactor packed length {z.size} != {expected}"
    offset = 0
    for i, cell in enumerate(cells):
        nv = n_var(cell)
        ok, why = packed_cell_in_closed_domain(z[offset : offset + nv], cell)
        if not ok:
            return False, f"cell[{i}]: {why}"
        offset += nv
    return True, ""


def unpack_cell_no_projection(x: np.ndarray, cell: Cell) -> None:
    """原样写入打包值：无 ``max``、无温度 fence。"""
    _write_packed_into_cell(
        np.asarray(x, dtype=np.float64),
        cell,
        clip_nonnegative=False,
        clip_temperature=False,
    )


def unpack_cell_strict_domain_no_projection(x: np.ndarray, cell: Cell) -> None:
    """合法则原样写入；非法则拒绝且不改 cell。

    Source: docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q9
    """
    ok, why = packed_cell_in_closed_domain(x, cell)
    if not ok:
        raise DomainViolation(why)
    unpack_cell_no_projection(x, cell)


def unpack_reactor_no_projection(x: np.ndarray, cells: List[Cell]) -> None:
    offset = 0
    z = np.asarray(x, dtype=np.float64)
    for cell in cells:
        nv = n_var(cell)
        unpack_cell_no_projection(z[offset : offset + nv], cell)
        offset += nv


def unpack_reactor_strict_domain_no_projection(x: np.ndarray, cells: List[Cell]) -> None:
    """先检查全部 cell，再写入；任一非法则全部不改。"""
    ok, why = packed_reactor_in_closed_domain(x, cells)
    if not ok:
        raise DomainViolation(why)
    unpack_reactor_no_projection(x, cells)


def snapshot_cell_packed_fields(cell: Cell) -> dict[str, np.ndarray | float]:
    """用于断言非法步不改写。"""
    return {
        "N_d": np.array(cell.N_d, dtype=np.float64, copy=True),
        "N_b": np.array(cell.N_b, dtype=np.float64, copy=True),
        "m_solid": np.array(cell.m_solid, dtype=np.float64, copy=True),
        "T": float(cell.T),
    }


def cell_packed_fields_equal(
    a: dict[str, np.ndarray | float],
    b: dict[str, np.ndarray | float],
) -> bool:
    if float(a["T"]) != float(b["T"]):
        return False
    return (
        np.array_equal(a["N_d"], b["N_d"])
        and np.array_equal(a["N_b"], b["N_b"])
        and np.array_equal(a["m_solid"], b["m_solid"])
    )


def max_feasible_alpha_closed_domain(
    x: np.ndarray,
    dx: np.ndarray,
    cells: List[Cell],
    *,
    alpha0: float = 1.0,
) -> tuple[float, str]:
    """具名分数边界步长：最大 ``α∈(0,α0]`` 使 ``x+α dx`` 仍在闭定义域。

    若某零组分的方向指向负域，返回 ``α=0`` 并注明需要活跃集，不设正数 floor。
    """
    z = np.asarray(x, dtype=np.float64)
    d = np.asarray(dx, dtype=np.float64)
    assert z.shape == d.shape, "x 与 dx 维度不一致"
    alpha = float(alpha0)
    policy = "fractional_boundary_to_closed_domain"
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        xi = z[offset : offset + nv]
        di = d[offset : offset + nv]
        nn = packed_nonnegative_slice(cell)
        for val, step in zip(xi[nn], di[nn]):
            if step < 0.0:
                if float(val) <= 0.0:
                    return 0.0, "active_set_required_zero_component"
                alpha = min(alpha, float(val) / (-float(step)))
        if has_temperature_var(cell):
            T = float(xi[-1])  # [K]
            dT = float(di[-1])  # [K]
            if dT < 0.0:
                if T <= CLOSED_DOMAIN_T_OPEN_LOWER_K:
                    return 0.0, "active_set_required_T"
                # 保持 T+α dT > 0
                alpha = min(alpha, (T / (-dT)) * (1.0 - 1.0e-15))
        offset += nv
    if alpha <= 0.0:
        return 0.0, "active_set_required"
    return float(min(alpha, float(alpha0))), policy


def pack_cell(cell: Cell) -> np.ndarray:
    gas_idx = main_nr_gas_idx_for_cell(cell)
    gas_count = int(gas_idx.size)
    solid_flat_idx = solid_flat_indices_for_nr(cell)
    gas_mode = gas_mode_for_nr(cell)
    parts: list[np.ndarray] = []
    if gas_mode == "two_phase":
        parts.extend([cell.N_d[gas_idx], cell.N_b[gas_idx]])
    elif gas_mode == "single":
        parts.append(cell.N_d[gas_idx])
    if solid_flat_idx.size:
        parts.append(cell.m_solid.reshape(-1)[solid_flat_idx])
    if has_temperature_var(cell):
        parts.append(np.array([cell.T], dtype=np.float64))
    if not parts:
        return np.zeros(0, dtype=np.float64)
    return np.concatenate(parts)


def unpack_cell(x: np.ndarray, cell: Cell) -> None:
    """生产默认：负物流裁零、温度 fence。不是 Q9 无投影路径。"""
    _write_packed_into_cell(
        np.asarray(x, dtype=np.float64),
        cell,
        clip_nonnegative=True,
        clip_temperature=True,
    )


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


def build_var_index_to_cell_map(offsets: List[int]) -> np.ndarray:
    n_total = int(offsets[-1])
    out = np.empty(n_total, dtype=np.int32)
    for cell_idx in range(len(offsets) - 1):
        out[offsets[cell_idx] : offsets[cell_idx + 1]] = cell_idx
    return out


def affected_residual_cells_for_var_cell(cell_idx: int, n_cells: int) -> tuple[int, ...]:
    """默认 bed 链式邻接 + 顶格与底格循环耦合（无显式 freeboard/side 图时）。"""
    affected: list[int] = []
    if cell_idx == 0:
        affected.extend(range(n_cells))
    if cell_idx > 0:
        affected.append(cell_idx - 1)
    affected.append(cell_idx)
    if cell_idx + 1 < n_cells:
        affected.append(cell_idx + 1)
    if cell_idx == n_cells - 1 and 0 not in affected:
        affected.append(0)
    return tuple(sorted(set(affected)))
