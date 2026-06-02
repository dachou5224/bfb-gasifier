"""全局 NR 未知量布局与 pack/unpack 契约（与 Hamel 全局联立思路对齐）。

主 NR 向量排除微量气相组分（如 H2S/NH3），由 Gibbs/后处理路径处理。
所有索引与 `global_nr_solver.solve_global_nr` 使用的布局必须与本模块保持一致。
"""

from __future__ import annotations

from typing import List

import numpy as np

from src.core.cell import Cell, N_SOLID_COMP, S_ASH, S_CHAR
from src.core.species import GAS_SPECIES_INDEX, N_GAS

# Hamel thesis alignment: minor species not in global Jacobian unknown vector.
# TAR1/TAR2 are also excluded: at gasification temperatures (>800 K) they crack
# completely and their NR balance rows are near-zero, making the Jacobian
# ill-conditioned (condition number ~10^16).  They are kept fixed at their
# vorabrechnung-seeded (near-zero) values throughout the inner NR solve.
MAIN_NR_EXCLUDED_SPECIES: tuple[str, ...] = ("H2S", "NH3", "TAR1", "TAR2")
_MAIN_EXCLUDED_IDX = {GAS_SPECIES_INDEX[s] for s in MAIN_NR_EXCLUDED_SPECIES}

MAIN_NR_GAS_IDX = np.array(
    [i for i in range(N_GAS) if i not in _MAIN_EXCLUDED_IDX],
    dtype=np.int64,
)
MAIN_NR_GAS_COUNT = int(MAIN_NR_GAS_IDX.size)


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
    mode = gas_mode_for_nr(cell)
    if mode == "none":
        return 0
    if mode == "single":
        return MAIN_NR_GAS_COUNT
    return 2 * MAIN_NR_GAS_COUNT


def solid_offset(cell: Cell) -> int:
    """Start offset of solid unknowns inside one packed cell block."""
    return n_gas_var(cell)


def temperature_offset(cell: Cell) -> int:
    """Offset of the cell temperature unknown inside one packed cell block."""
    if not has_temperature_var(cell):
        raise ValueError(f"{cell_kind(cell)} has no NR temperature unknown")
    return n_gas_var(cell) + n_solid_var(cell)


def has_temperature_var(cell: Cell) -> bool:
    """Whether this cell carries an independent NR temperature unknown."""
    return cell_kind(cell) not in {"cyclone", "return_leg"}


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


def pack_cell(cell: Cell) -> np.ndarray:
    solid_flat_idx = solid_flat_indices_for_nr(cell)
    gas_mode = gas_mode_for_nr(cell)
    parts: list[np.ndarray] = []
    if gas_mode == "two_phase":
        parts.extend([cell.N_d[MAIN_NR_GAS_IDX], cell.N_b[MAIN_NR_GAS_IDX]])
    elif gas_mode == "single":
        parts.append(cell.N_d[MAIN_NR_GAS_IDX])
    if solid_flat_idx.size:
        parts.append(cell.m_solid.reshape(-1)[solid_flat_idx])
    if has_temperature_var(cell):
        parts.append(np.array([cell.T], dtype=np.float64))
    if not parts:
        return np.zeros(0, dtype=np.float64)
    return np.concatenate(parts)


def unpack_cell(x: np.ndarray, cell: Cell) -> None:
    gas_mode = gas_mode_for_nr(cell)
    nv_sol = n_solid_var(cell)
    offset = 0
    if gas_mode == "two_phase":
        nd = np.maximum(x[offset : offset + MAIN_NR_GAS_COUNT], 0.0)
        offset += MAIN_NR_GAS_COUNT
        nb = np.maximum(x[offset : offset + MAIN_NR_GAS_COUNT], 0.0)
        offset += MAIN_NR_GAS_COUNT
        cell.N_d[MAIN_NR_GAS_IDX] = nd
        cell.N_b[MAIN_NR_GAS_IDX] = nb
    elif gas_mode == "single":
        nd = np.maximum(x[offset : offset + MAIN_NR_GAS_COUNT], 0.0)
        offset += MAIN_NR_GAS_COUNT
        cell.N_d[MAIN_NR_GAS_IDX] = nd
        cell.N_b[MAIN_NR_GAS_IDX] = 0.0
    else:
        cell.N_d[MAIN_NR_GAS_IDX] = 0.0
        cell.N_b[MAIN_NR_GAS_IDX] = 0.0
    if nv_sol > 0:
        m_flat = np.maximum(x[offset : offset + nv_sol], 0.0)
        solid_flat = cell.m_solid.reshape(-1)
        solid_flat[solid_flat_indices_for_nr(cell)] = m_flat
        offset += nv_sol
    if has_temperature_var(cell):
        T = float(x[-1])
        t_min = float(getattr(cell, "_nr_temperature_min_K", 300.0))
        t_max = float(getattr(cell, "_nr_temperature_max_K", 2500.0))
        if not np.isfinite(t_min):
            t_min = 300.0
        if not np.isfinite(t_max):
            t_max = 2500.0
        t_min = max(t_min, 300.0)
        t_max = min(max(t_max, t_min + 1.0), 2500.0)
        if T < t_min:
            cell.T = t_min
        elif T > t_max:
            cell.T = t_max
        else:
            cell.T = T


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
