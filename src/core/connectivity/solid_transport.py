"""Connectivity submodule — see package ``__init__`` for public API."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from src.core.cell import Cell, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_gas_enthalpy_flow, calc_solid_enthalpy_flow
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX, cp_ash, cp_char, cp_sand
from src.core.connectivity_graph import (
    ConnectivityBlock,
    ConnectivityEdge,
    ConnectivityGraph,
    bed_top_freeboard_ab_source_block_id,
    block_id_bed_index,
    build_connectivity_graph_from_reactor_params,
    cyclone_inlet_source_block_id,
    axial_upstream_block_id,
    freeboard_chain_inlet_source_block_id,
    infer_stream_class,
    ordered_freeboard_block_ids,
    resolve_bed_neighbor_block_ids,
    resolve_registry_cell,
    solid_ab_freeboard_downstream_block_id,
    solid_ab_source_block_id,
    solid_auf_source_block_id,
)
from src.solvers.nr_indexing import (
    MAJOR_NR_GAS_SPECIES,
    nr_dense_only_pyrolysis_support_species,
    nr_projection_gas_species,
)

_F_W = 0.25
_ZETA_W = 0.40

SOLID_TRANSPORT_INVENTORY_ENGINEERING: str = "engineering_bottom_threshold"
SOLID_TRANSPORT_INVENTORY_SPECIES_OWN: str = "species_own_inventory"
BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING: float = 0.15
_INVENTORY_MODES = frozenset(
    {SOLID_TRANSPORT_INVENTORY_ENGINEERING, SOLID_TRANSPORT_INVENTORY_SPECIES_OWN}
)


def resolve_solid_transport_inventory_mode(configured: str) -> str:
    """校验 ``solid_transport_inventory_mode``。"""
    name = str(configured or SOLID_TRANSPORT_INVENTORY_ENGINEERING).strip().lower()
    if name not in _INVENTORY_MODES:
        raise ValueError(
            f"Unknown solid_transport_inventory_mode={configured!r}; "
            f"expected {SOLID_TRANSPORT_INVENTORY_ENGINEERING!r} or "
            f"{SOLID_TRANSPORT_INVENTORY_SPECIES_OWN!r}"
        )
    return name


def cell_active_char_ash_holdup(cell: Cell) -> float:
    """char+ash 库存 [kg]，不含 VM/水分/惰性几何床料。"""
    return float(
        np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))
        + np.sum(np.maximum(cell.m_solid[:, S_ASH], 0.0))
    )

def recycled_solid_stream(m_solid_top: np.ndarray, frac: float) -> np.ndarray:
    """Recycle stream for the current thesis-aligned bed model."""
    m_recycle = np.zeros_like(m_solid_top)
    m_recycle[:, S_CHAR] = np.maximum(m_solid_top[:, S_CHAR], 0.0) * float(frac)
    m_recycle[:, S_ASH] = np.maximum(m_solid_top[:, S_ASH], 0.0) * float(frac)
    return m_recycle

def propagated_solid_stream(
    m_solid: np.ndarray,
    frac: float = 1.0,
    *,
    include_reactive: bool = False,
) -> np.ndarray:
    """Axial solid propagation for the current thesis-aligned bed model."""
    out = np.zeros_like(m_solid)
    out[:, S_CHAR] = np.maximum(m_solid[:, S_CHAR], 0.0) * float(frac)
    out[:, S_ASH] = np.maximum(m_solid[:, S_ASH], 0.0) * float(frac)
    if include_reactive:
        out[:, S_VM] = np.maximum(m_solid[:, S_VM], 0.0) * float(frac)
        out[:, S_MOISTURE] = np.maximum(m_solid[:, S_MOISTURE], 0.0) * float(frac)
    return out

def propagated_solid_stream_from_cell(
    cell: Cell,
    frac: float = 1.0,
    *,
    include_reactive: bool = False,
) -> np.ndarray:
    """Project the cell's outgoing solid transport stream onto routed components."""
    upflow = cell._solid_upflow_rates()
    routed = propagated_solid_stream(upflow, frac=frac, include_reactive=False)
    if include_reactive:
        if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
            # Propagate the remaining resident reactive stock, not the original
            # fresh/source support. Reusing m_solid_zu + m_solid_in would count
            # the same VM/moisture budget again in every lower-bed cell.
            reactive_source = np.maximum(cell.m_solid, 0.0)
        else:
            reactive_source = np.maximum(upflow, 0.0)
        routed[:, S_VM] = reactive_source[:, S_VM] * float(frac)
        routed[:, S_MOISTURE] = reactive_source[:, S_MOISTURE] * float(frac)
    return routed

def recycled_solid_stream_from_cell(cell: Cell, frac: float) -> np.ndarray:
    """Recycle only the routed fraction of the cell's outgoing solid stream."""
    return recycled_solid_stream(cell._solid_upflow_rates(), frac)

def project_zero_source_solid_components(cell: Cell) -> None:
    """Project VM/moisture state onto the physically allowed inlet support."""
    vm_in = float(
        np.sum(
            np.maximum(
                cell.m_solid_zu[:, S_VM]
                + cell.m_solid_in[:, S_VM]
                + cell.m_solid_rez[:, S_VM],
                0.0,
            )
        )
    )
    if vm_in <= 1e-12:
        cell.m_solid[:, S_VM] = 0.0

    moist_in = float(
        np.sum(
            np.maximum(
                cell.m_solid_zu[:, S_MOISTURE]
                + cell.m_solid_in[:, S_MOISTURE]
                + cell.m_solid_rez[:, S_MOISTURE],
                0.0,
            )
        )
    )
    if moist_in <= 1e-12:
        cell.m_solid[:, S_MOISTURE] = 0.0

def cell_total_solid_holdup(cell: Cell) -> float:
    """Estimate total solid holdup [kg] from local hydrodynamics and geometry.

    Ref: Hamel (1999) Eq. 2.6 coupling interpretation; user-provided bed discretization
    for ``M_solid = A * dh * (1-eps_b) * rho_p * (1-eps_d)``.
    """
    area = float(np.pi * cell.geo.D_bed**2 / 4.0)
    return float(
        max(area, 0.0)
        * max(float(cell.geo.dh), 0.0)
        * max(1.0 - float(cell.eps_b), 0.0)
        * max(float(cell.solid.rho_s), 0.0)
        * max(1.0 - float(cell.eps_d_voidage), 0.0)
    )


def cell_resident_holdup_for_transport_closure(
    cell: Cell,
    *,
    is_bottom_bed_cell: bool = False,
    inventory_mode: str = SOLID_TRANSPORT_INVENTORY_ENGINEERING,
) -> float:
    """``K_ab`` / ``m_auf`` 连续式闭合质量 [kg]；见 ``update_bed_solid_transport_coefficients``。

    ``engineering_bottom_threshold``：上层用几何库存；床底在活性库存 < 15% 几何库存时
    改用 char+ash。这是工程开关，不是 Hamel Eq. (2-6)。
    ``species_own_inventory``：所有格都返回 char+ash，供诊断；实际 ``K`` 按物种分母重建。

    Ref: Hamel (1999) Eq.2.6; docs/audit_reports/secondary_diagnosis_2026-09-14.md §5
    """
    geo = cell_total_solid_holdup(cell)
    mode = resolve_solid_transport_inventory_mode(inventory_mode)
    resident = cell_active_char_ash_holdup(cell)
    if mode == SOLID_TRANSPORT_INVENTORY_SPECIES_OWN:
        return resident
    if not is_bottom_bed_cell:
        return geo
    if str(getattr(cell, "solid_state_model", "legacy_stream")) not in {
        "holdup_transport",
        "freeboard_closure",
    }:
        return geo
    if resident <= 1e-9 or resident >= BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING * geo:
        return geo
    return resident

def bed_upflow_coefficient(cell: Cell, *, is_top_bed_cell: bool) -> float:
    """Hamel-style wake carrying coefficient ``K_auf`` for a bed cell [1/s]."""
    denom = max((1.0 - float(cell.eps_b)) * float(cell.geo.dh), 1e-12)
    k_auf = _F_W * max(float(cell.eps_b), 0.0) * max(float(cell.u_b), 0.0) / denom
    if is_top_bed_cell:
        k_auf *= _ZETA_W
    return float(max(k_auf, 0.0))

def _update_bed_k_species_own(
    bed_cells: list[Cell],
    *,
    top_downflow_components: np.ndarray | None,
    bottom_no_penetration: bool,
) -> None:
    """按物种库存重建 ``K``：``m_up,j = K_up M_j``，``K_ab,j = m_ab,j / M_j``。"""
    n = len(bed_cells)
    k_auf_vals = np.array(
        [bed_upflow_coefficient(cell, is_top_bed_cell=(i == n - 1)) for i, cell in enumerate(bed_cells)],
        dtype=np.float64,
    )
    m_char = np.array(
        [float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) for cell in bed_cells],
        dtype=np.float64,
    )
    m_ash = np.array(
        [float(np.sum(np.maximum(cell.m_solid[:, S_ASH], 0.0))) for cell in bed_cells],
        dtype=np.float64,
    )
    m_up_char = k_auf_vals * m_char
    m_up_ash = k_auf_vals * m_ash
    from_above = np.zeros(2, dtype=np.float64)
    if isinstance(top_downflow_components, np.ndarray) and top_downflow_components.shape[0] >= 2:
        from_above[:] = np.maximum(top_downflow_components[:2], 0.0)

    for i in range(n - 1, -1, -1):
        cell = bed_cells[i]
        ext_char = float(
            np.sum(
                np.maximum(
                    cell.m_solid_zu[:, S_CHAR] + cell.m_solid_rez[:, S_CHAR] + cell.m_solid_in[:, S_CHAR],
                    0.0,
                )
            )
        )
        ext_ash = float(
            np.sum(
                np.maximum(
                    cell.m_solid_zu[:, S_ASH] + cell.m_solid_rez[:, S_ASH] + cell.m_solid_in[:, S_ASH],
                    0.0,
                )
            )
        )
        rxn_char = float(np.sum(np.asarray(cell.R_solid[:, S_CHAR], dtype=np.float64)))
        rxn_ash = float(np.sum(np.asarray(cell.R_solid[:, S_ASH], dtype=np.float64)))
        up_in_char = float(m_up_char[i - 1]) if i > 0 else 0.0
        up_in_ash = float(m_up_ash[i - 1]) if i > 0 else 0.0
        if i == 0 and bottom_no_penetration:
            m_ab_char = 0.0
            m_ab_ash = 0.0
        else:
            m_ab_char = max(up_in_char - float(m_up_char[i]) + from_above[0] + ext_char + rxn_char, 0.0)
            m_ab_ash = max(up_in_ash - float(m_up_ash[i]) + from_above[1] + ext_ash + rxn_ash, 0.0)
        k_ab_char = m_ab_char / m_char[i] if m_char[i] > 1.0e-12 else 0.0
        k_ab_ash = m_ab_ash / m_ash[i] if m_ash[i] > 1.0e-12 else 0.0
        cell.K_solid_auf.fill(0.0)
        cell.K_solid_ab.fill(0.0)
        cell.K_solid_auf[:, S_CHAR] = float(k_auf_vals[i])
        cell.K_solid_auf[:, S_ASH] = float(k_auf_vals[i])
        cell.K_solid_ab[:, S_CHAR] = float(k_ab_char)
        cell.K_solid_ab[:, S_ASH] = float(k_ab_ash)
        from_above[0] = m_ab_char
        from_above[1] = m_ab_ash


def update_bed_solid_transport_coefficients(
    bed_cells: list[Cell],
    cfg: Any,
    *,
    top_downflow_total: float = 0.0,
    top_downflow_components: np.ndarray | None = None,
) -> None:
    """更新床层 ``K_auf`` / ``K_ab``。

    ``engineering_bottom_threshold``（默认）：总闭合库存作分母，床底 15% 切换。
    ``species_own_inventory``：``K_ab,j = m_ab,j / M_j``，无 15% 跳变。
    ``solid_transport_bottom_no_penetration``：床底 ``K_ab=0``。

    Ref: Hamel (1999) Eq. (2-6); docs/audit_reports/secondary_diagnosis_2026-09-14.md §5
    """
    if not bed_cells:
        return

    inventory_mode = resolve_solid_transport_inventory_mode(
        getattr(cfg, "solid_transport_inventory_mode", SOLID_TRANSPORT_INVENTORY_ENGINEERING)
    )
    bottom_no_penetration = bool(getattr(cfg, "solid_transport_bottom_no_penetration", False))
    if inventory_mode == SOLID_TRANSPORT_INVENTORY_SPECIES_OWN:
        _update_bed_k_species_own(
            bed_cells,
            top_downflow_components=top_downflow_components,
            bottom_no_penetration=bottom_no_penetration,
        )
        return

    holdups = np.array(
        [
            cell_resident_holdup_for_transport_closure(
                cell,
                is_bottom_bed_cell=(i == 0),
                inventory_mode=inventory_mode,
            )
            for i, cell in enumerate(bed_cells)
        ],
        dtype=np.float64,
    )
    k_auf_vals = np.array(
        [bed_upflow_coefficient(cell, is_top_bed_cell=(i == len(bed_cells) - 1)) for i, cell in enumerate(bed_cells)],
        dtype=np.float64,
    )
    m_auf_totals = k_auf_vals * holdups

    m_ab_from_above_total = float(max(top_downflow_total, 0.0))
    m_ab_from_above_components = np.zeros(2, dtype=np.float64)
    if isinstance(top_downflow_components, np.ndarray) and top_downflow_components.shape[0] >= 2:
        m_ab_from_above_components[:] = np.maximum(top_downflow_components[:2], 0.0)
    elif m_ab_from_above_total > 0.0:
        m_ab_from_above_components[:] = 0.5 * m_ab_from_above_total

    for i in range(len(bed_cells) - 1, -1, -1):
        cell = bed_cells[i]
        up_in = float(m_auf_totals[i - 1]) if i > 0 else 0.0
        ext_in_char = float(np.sum(np.maximum(cell.m_solid_zu[:, S_CHAR] + cell.m_solid_rez[:, S_CHAR] + cell.m_solid_in[:, S_CHAR], 0.0)))
        ext_in_ash = float(np.sum(np.maximum(cell.m_solid_zu[:, S_ASH] + cell.m_solid_rez[:, S_ASH] + cell.m_solid_in[:, S_ASH], 0.0)))
        reaction_char = float(np.sum(np.asarray(cell.R_solid[:, S_CHAR], dtype=np.float64)))
        reaction_ash = float(np.sum(np.asarray(cell.R_solid[:, S_ASH], dtype=np.float64)))
        # Bed axial upflow remains hydrodynamics-driven (shared ``K_auf``), while the
        # top-down closure for ``K_ab`` keeps char/ash source terms separated.
        m_ab_total = max(
            up_in - float(m_auf_totals[i]) + m_ab_from_above_total + ext_in_char + ext_in_ash + reaction_char + reaction_ash,
            0.0,
        )
        m_ab_char = max(up_in - float(m_auf_totals[i]) + m_ab_from_above_components[0] + ext_in_char + reaction_char, 0.0)
        m_ab_ash = max(up_in - float(m_auf_totals[i]) + m_ab_from_above_components[1] + ext_in_ash + reaction_ash, 0.0)
        k_ab_char = m_ab_char / max(float(holdups[i]), 1e-12) if holdups[i] > 1e-12 else 0.0
        k_ab_ash = m_ab_ash / max(float(holdups[i]), 1e-12) if holdups[i] > 1e-12 else 0.0

        cell.K_solid_auf.fill(0.0)
        cell.K_solid_ab.fill(0.0)
        # Thesis holdup transport: only char/ash are transported across cells.
        cell.K_solid_auf[:, S_CHAR] = float(k_auf_vals[i])
        cell.K_solid_auf[:, S_ASH] = float(k_auf_vals[i])
        cell.K_solid_ab[:, S_CHAR] = float(k_ab_char)
        cell.K_solid_ab[:, S_ASH] = float(k_ab_ash)
        m_ab_from_above_total = m_ab_total
        m_ab_from_above_components[0] = m_ab_char
        m_ab_from_above_components[1] = m_ab_ash

def _zero_reactive_solid_transport_components(cell: Cell) -> None:
    """Thesis holdup transport: char/ash only in axial K channels."""
    cell.m_solid_auf_in[:, S_VM] = 0.0
    cell.m_solid_auf_in[:, S_MOISTURE] = 0.0
    cell.m_solid_ab_in[:, S_VM] = 0.0
    cell.m_solid_ab_in[:, S_MOISTURE] = 0.0

def _set_bed_solid_auf_in_from_source(cell: Cell, source: Cell) -> None:
    cell.m_solid_auf_in[:, :] = source.K_solid_auf * np.maximum(source.m_solid, 0.0)
    cell.T_solid_auf_in = float(source.T)

def _set_bed_solid_ab_in_from_bed_source(cell: Cell, source: Cell) -> None:
    cell.m_solid_ab_in[:, :] = source.K_solid_ab * np.maximum(source.m_solid, 0.0)
    cell.T_solid_ab_in = float(source.T)

def _set_bed_solid_ab_in_from_freeboard_source(cell: Cell, source: Cell) -> None:
    cell.m_solid_ab_in[:, :] = np.maximum(source._solid_downflow_rates(), 0.0)
    cell.T_solid_ab_in = float(source.T)

def update_bed_solid_transport_inflows(bed_cells: list[Cell]) -> None:
    """Map neighboring bed-cell holdup onto thesis-style axial solid inflow terms."""
    update_bed_solid_transport_inflows_from_above(bed_cells)

def update_bed_solid_transport_inflows_from_above(
    bed_cells: list[Cell],
    top_above_cell: Cell | None = None,
    *,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Map neighboring bed/freeboard holdup onto thesis-style bed inflow terms."""
    if not bed_cells:
        return
    for cell in bed_cells:
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)
        cell.T_solid_auf_in = float(cell.T)
        cell.T_solid_ab_in = float(cell.T)

    if connectivity_graph is not None and cell_registry is not None:
        for i, cell in enumerate(bed_cells):
            block_id = f"bed_cell_{i}"
            auf_src = resolve_registry_cell(
                connectivity_graph,
                cell_registry,
                solid_auf_source_block_id(connectivity_graph, block_id),
            )
            if auf_src is not None:
                _set_bed_solid_auf_in_from_source(cell, auf_src)

            ab_block = solid_ab_source_block_id(connectivity_graph, block_id)
            ab_src = resolve_registry_cell(connectivity_graph, cell_registry, ab_block)
            if ab_src is not None:
                _set_bed_solid_ab_in_from_bed_source(cell, ab_src)
            else:
                fb_block = bed_top_freeboard_ab_source_block_id(connectivity_graph)
                if fb_block is not None and block_id == f"bed_cell_{len(bed_cells) - 1}":
                    fb_src = resolve_registry_cell(connectivity_graph, cell_registry, fb_block)
                    if fb_src is not None:
                        _set_bed_solid_ab_in_from_freeboard_source(cell, fb_src)
                elif top_above_cell is not None and i == len(bed_cells) - 1:
                    _set_bed_solid_ab_in_from_freeboard_source(cell, top_above_cell)
            _zero_reactive_solid_transport_components(cell)
        return

    for i, cell in enumerate(bed_cells):
        if i > 0:
            below = bed_cells[i - 1]
            _set_bed_solid_auf_in_from_source(cell, below)
        if i + 1 < len(bed_cells):
            above = bed_cells[i + 1]
            _set_bed_solid_ab_in_from_bed_source(cell, above)
        elif top_above_cell is not None:
            _set_bed_solid_ab_in_from_freeboard_source(cell, top_above_cell)
        _zero_reactive_solid_transport_components(cell)

def propagate_explicit_freeboard_chain(
    bed_cells: list[Cell],
    freeboard_cells: list[Cell],
    *,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Propagate axial gas/solid inlets from bed top through explicit freeboard cells."""
    if not bed_cells or not freeboard_cells:
        return

    if connectivity_graph is not None and cell_registry is not None:
        fb_ids = ordered_freeboard_block_ids(connectivity_graph)
        if not fb_ids:
            return
        inlet_id = freeboard_chain_inlet_source_block_id(connectivity_graph)
        prev = cell_registry.get(inlet_id) if inlet_id else bed_cells[-1]
        if prev is None:
            prev = bed_cells[-1]
        for fb_id in fb_ids:
            cell = cell_registry.get(fb_id)
            if cell is None:
                continue
            _set_explicit_freeboard_inlet_from_prev(prev, cell)
            prev = cell
        return

    prev = bed_cells[-1]
    for cell in freeboard_cells:
        _set_explicit_freeboard_inlet_from_prev(prev, cell)
        prev = cell

def _set_degenerate_single_phase_gas_inlet(prev: Cell, cell: Cell) -> None:
    """Route upstream total gas into the suspension channel of degenerate segments."""
    prev_kind = str(getattr(prev, "cell_type", "bed")).strip().lower()
    prev_total = prev.N_d if prev_kind in {"freeboard", "cyclone", "return_leg"} else prev.N_d + prev.N_b
    cell.N_d_in[:] = np.maximum(prev_total, 0.0)
    cell.N_b_in.fill(0.0)

def _set_explicit_freeboard_inlet_from_prev(prev: Cell, cell: Cell) -> None:
    """Refresh one explicit freeboard cell inlet from its upstream source cell."""
    _set_degenerate_single_phase_gas_inlet(prev, cell)
    cell.T_in_gas = float(prev.T)
    cell.T_in_solid = float(prev.T)
    cell.m_solid_in.fill(0.0)
    # Freeboard-closure solids remain outside the NR solid residual, but their
    # entrained char/ash stream still carries sensible/formation enthalpy into
    # the freeboard energy residual.
    cell.m_solid_in[:] = propagated_solid_stream_from_cell(prev, include_reactive=False)
    active_mask = getattr(cell, "_freeboard_active_char_ash_mask", None)
    if (
        str(getattr(cell, "solid_state_model", "legacy_stream")) != "freeboard_closure"
        and isinstance(active_mask, np.ndarray)
        and active_mask.shape[0] == cell.solid.n_size_classes
    ):
        mask = np.maximum(active_mask, 0.0).reshape(-1)
        cell.m_solid_in[:, S_CHAR] *= mask
        cell.m_solid_in[:, S_ASH] *= mask
    if str(getattr(cell, "solid_state_model", "legacy_stream")) == "holdup_transport":
        vm_col = np.maximum(cell.m_solid_in[:, S_VM], 0.0)
        moist_col = np.maximum(cell.m_solid_in[:, S_MOISTURE], 0.0)
        cell.m_solid_in.fill(0.0)
        cell.m_solid_in[:, S_VM] = vm_col
        cell.m_solid_in[:, S_MOISTURE] = moist_col
    cell.T_in_solid = float(prev.T)

def update_freeboard_solid_transport_inflows(
    freeboard_cells: list[Cell],
    *,
    bottom_below_cell: Cell,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Project explicit freeboard K_auf/K_ab into neighboring inflow terms.

    Hamel-aligned ``freeboard_closure`` cells keep closure hold-up and closure
    outflow coefficients, but do not participate in explicit Eulerian solid
    inflow chaining inside the NR boundary path.
    """
    if not freeboard_cells:
        return
    for cell in freeboard_cells:
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)
        cell.T_solid_auf_in = float(cell.T)
        cell.T_solid_ab_in = float(cell.T)

    if connectivity_graph is not None and cell_registry is not None:
        for fb_id in ordered_freeboard_block_ids(connectivity_graph):
            cell = cell_registry.get(fb_id)
            if cell is None or str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
                continue
            auf_block = solid_auf_source_block_id(connectivity_graph, fb_id)
            if auf_block is None:
                auf_block = axial_upstream_block_id(connectivity_graph, fb_id)
            auf_src = resolve_registry_cell(connectivity_graph, cell_registry, auf_block)
            if auf_src is None and fb_id == "freeboard_cell_0":
                auf_src = bottom_below_cell
            if auf_src is not None:
                cell.m_solid_auf_in[:, :] = np.maximum(auf_src._solid_upflow_rates(), 0.0)
                cell.T_solid_auf_in = float(auf_src.T)
            ab_block = solid_ab_freeboard_downstream_block_id(connectivity_graph, fb_id)
            ab_src = resolve_registry_cell(
                connectivity_graph,
                cell_registry,
                ab_block,
                direction="downstream",
            )
            if ab_src is not None:
                cell.m_solid_ab_in[:, :] = np.maximum(ab_src._solid_downflow_rates(), 0.0)
                cell.T_solid_ab_in = float(ab_src.T)
            _zero_reactive_solid_transport_components(cell)
            _apply_freeboard_active_class_mask_to_transport_inflows(cell)
        return

    for i, cell in enumerate(freeboard_cells):
        if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
            continue
        if i == 0:
            cell.m_solid_auf_in[:, :] = np.maximum(bottom_below_cell._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(bottom_below_cell.T)
        else:
            below = freeboard_cells[i - 1]
            cell.m_solid_auf_in[:, :] = np.maximum(below._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(below.T)
        if i + 1 < len(freeboard_cells):
            above = freeboard_cells[i + 1]
            cell.m_solid_ab_in[:, :] = np.maximum(above._solid_downflow_rates(), 0.0)
            cell.T_solid_ab_in = float(above.T)
        _zero_reactive_solid_transport_components(cell)
        _apply_freeboard_active_class_mask_to_transport_inflows(cell)

def _refresh_freeboard_solid_transport_inflows_for_indices(
    freeboard_cells: list[Cell],
    *,
    bottom_below_cell: Cell,
    indices: set[int],
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Refresh only selected explicit freeboard transport inflow rows."""
    if not freeboard_cells:
        return
    if connectivity_graph is not None and cell_registry is not None:
        fb_ids = ordered_freeboard_block_ids(connectivity_graph)
        selected = {
            fb_ids[i]
            for i in indices
            if 0 <= int(i) < len(fb_ids)
        }
        for fb_id in selected:
            cell = cell_registry.get(fb_id)
            if cell is None:
                continue
            cell.m_solid_auf_in.fill(0.0)
            cell.m_solid_ab_in.fill(0.0)
            cell.T_solid_auf_in = float(cell.T)
            cell.T_solid_ab_in = float(cell.T)
            if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
                continue
            auf_block = solid_auf_source_block_id(connectivity_graph, fb_id) or axial_upstream_block_id(
                connectivity_graph, fb_id
            )
            auf_src = resolve_registry_cell(connectivity_graph, cell_registry, auf_block)
            if auf_src is None and fb_id == "freeboard_cell_0":
                auf_src = bottom_below_cell
            if auf_src is not None:
                cell.m_solid_auf_in[:, :] = np.maximum(auf_src._solid_upflow_rates(), 0.0)
                cell.T_solid_auf_in = float(auf_src.T)
            ab_block = solid_ab_freeboard_downstream_block_id(connectivity_graph, fb_id)
            ab_src = resolve_registry_cell(
                connectivity_graph,
                cell_registry,
                ab_block,
                direction="downstream",
            )
            if ab_src is not None:
                cell.m_solid_ab_in[:, :] = np.maximum(ab_src._solid_downflow_rates(), 0.0)
                cell.T_solid_ab_in = float(ab_src.T)
            _zero_reactive_solid_transport_components(cell)
            _apply_freeboard_active_class_mask_to_transport_inflows(cell)
        return

    n = len(freeboard_cells)
    for i in sorted(int(idx) for idx in indices if 0 <= int(idx) < n):
        cell = freeboard_cells[i]
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)
        cell.T_solid_auf_in = float(cell.T)
        cell.T_solid_ab_in = float(cell.T)
        if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
            continue
        if i == 0:
            cell.m_solid_auf_in[:, :] = np.maximum(bottom_below_cell._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(bottom_below_cell.T)
        else:
            below = freeboard_cells[i - 1]
            cell.m_solid_auf_in[:, :] = np.maximum(below._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(below.T)
        if i + 1 < n:
            above = freeboard_cells[i + 1]
            cell.m_solid_ab_in[:, :] = np.maximum(above._solid_downflow_rates(), 0.0)
            cell.T_solid_ab_in = float(above.T)
        cell.m_solid_auf_in[:, S_VM] = 0.0
        cell.m_solid_auf_in[:, S_MOISTURE] = 0.0
        cell.m_solid_ab_in[:, S_VM] = 0.0
        cell.m_solid_ab_in[:, S_MOISTURE] = 0.0
        _apply_freeboard_active_class_mask_to_transport_inflows(cell)

def _apply_freeboard_active_class_mask_to_transport_inflows(cell: Cell) -> None:
    """Mask char/ash transport inflows by closure-active freeboard size classes."""
    active_mask = getattr(cell, "_freeboard_active_char_ash_mask", None)
    if not (isinstance(active_mask, np.ndarray) and active_mask.shape[0] == cell.solid.n_size_classes):
        return
    mask = np.maximum(active_mask, 0.0).reshape(-1)
    cell.m_solid_auf_in[:, S_CHAR] *= mask
    cell.m_solid_auf_in[:, S_ASH] *= mask
    cell.m_solid_ab_in[:, S_CHAR] *= mask
    cell.m_solid_ab_in[:, S_ASH] *= mask

def _sideblock_capture_efficiency(cfg: Any, comp_idx: int) -> float:
    """Cyclone capture efficiency for routed solid components."""
    if comp_idx == S_CHAR:
        return float(np.clip(getattr(cfg, "freeboard_cyclone_capture_char_frac", 0.0), 0.0, 1.0))
    if comp_idx == S_ASH:
        return float(np.clip(getattr(cfg, "freeboard_cyclone_capture_ash_frac", 0.0), 0.0, 1.0))
    return 0.0

SIDEBLOCK_K_MODEL: str = "engineering_gas_residence_capture_split"


def _estimate_cyclone_residence_time_s(cyclone_cell: Cell) -> float:
    """工程气体停留 ``V/V̇``，裁到 0.1–1 s。不能由此推出颗粒停留或分离。"""
    gas_total = float(np.sum(np.maximum(cyclone_cell.N_b + cyclone_cell.N_d, 0.0)))
    if gas_total <= 1e-12:
        gas_total = float(np.sum(np.maximum(cyclone_cell.N_b_in + cyclone_cell.N_d_in, 0.0)))
    T_ref = float(cyclone_cell.T if gas_total > 1e-12 and cyclone_cell.T > 0.0 else cyclone_cell.T_in_gas)
    if gas_total <= 1e-12 or T_ref <= 0.0:
        return 1.0
    gas_vol_flow = gas_total * Rg * T_ref / max(float(cyclone_cell.P), 1e-12)
    V_cyc = max(float(cyclone_cell.V_d), float(np.pi * cyclone_cell.geo.D_bed**2 * cyclone_cell.geo.dh / 4.0), 1e-12)
    dt = V_cyc / max(gas_vol_flow, 1e-12)
    return float(np.clip(dt, 0.1, 1.0))

def update_side_block_solid_transport_coefficients(
    cyclone_cell: Cell,
    return_leg_cell: Cell,
    cfg: Any,
) -> None:
    """Update engineering_gas_residence_capture_split side-block ``K = 1/dt``.

    原文无旋风 K / η 公式。这里用气体 ``V/V̇`` 裁到 0.1–1 s，再按捕集分数
    拆 ``K_auf/K_ab``。Eq. 2.6 只提供 ``m = K M`` 拓扑，不是分离规律。
    """
    K_cyc = 1.0 / max(_estimate_cyclone_residence_time_s(cyclone_cell), 1e-12)
    cyclone_cell.K_solid_auf.fill(0.0)
    cyclone_cell.K_solid_ab.fill(0.0)
    for comp_idx in (S_CHAR, S_ASH):
        eta = _sideblock_capture_efficiency(cfg, comp_idx)
        cyclone_cell.K_solid_auf[:, comp_idx] = (1.0 - eta) * K_cyc
        cyclone_cell.K_solid_ab[:, comp_idx] = eta * K_cyc

    leg_capacity = float(
        max(return_leg_cell.V_d, np.pi * return_leg_cell.geo.D_bed**2 * return_leg_cell.geo.dh / 4.0)
        * max(1.0 - float(return_leg_cell.solid.eps_mf), 0.0)
        * max(float(return_leg_cell.solid.rho_s), 0.0)
    )
    leg_throughput = float(np.sum(np.maximum(cyclone_cell._solid_downflow_rates(), 0.0)))
    K_leg = leg_throughput / max(leg_capacity, 1e-12) if leg_capacity > 1e-12 else 0.0
    return_leg_cell.K_solid_auf.fill(0.0)
    return_leg_cell.K_solid_ab.fill(0.0)
    for comp_idx in (S_CHAR, S_ASH):
        return_leg_cell.K_solid_ab[:, comp_idx] = K_leg

def _transport_inflow_support_holdup(
    cell: Cell,
    *,
    include_negative_R: bool = False,
) -> np.ndarray:
    """Steady Eq.2.6 resident solid inventory implied by current inflows and ``K`` [kg].

    Default：只把正 ``R_solid``（如热解产炭）计入 support，与历史 init seed 一致。
    ``include_negative_R=True``：用 ``in + R``（反应汇也降低目标 holdup），再 clip 到 ≥0
    （handoff §5.22；由 ``nr_holdup_support_include_negative_R_thesis`` 打开）。

    Ref: Hamel (1999) Eq.2.6 holdup transport closure
    """
    inflow = np.maximum(
        cell.m_solid_zu
        + cell.m_solid_rez
        + cell.m_solid_in
        + cell.m_solid_auf_in
        + cell.m_solid_ab_in,
        0.0,
    )
    R = np.asarray(cell.R_solid, dtype=np.float64)
    if include_negative_R:
        # net source；净消耗时目标 holdup→0，避免贴「忽略汇」的偏高 support
        inflow = np.maximum(inflow + R, 0.0)
    else:
        # Positive local solid formation (e.g. char from devolatilization) also supports holdup.
        inflow = inflow + np.maximum(R, 0.0)
    K_total = np.maximum(cell.K_solid_auf + cell.K_solid_ab, 0.0)
    support = np.zeros_like(cell.m_solid, dtype=np.float64)
    mask = K_total > 1e-12
    support[mask] = inflow[mask] / K_total[mask]
    return support


def seed_holdup_from_inflows(cell: Cell) -> bool:
    """Seed one holdup state from current transport inflows and ``K`` coefficients.

    For thesis ``holdup_transport`` states, upper bed cells depend on lower-cell
    ``K_auf * m_solid`` transport. During initialization this creates a one-way
    dependency chain: bed0 is known from fresh feed, bed1 depends on bed0, bed2
    depends on bed1, and so on. This helper raises under-filled holdup states
    toward their current steady ``inflow / K`` support without lowering existing
    inventory, so it is safe as an initialization sweep before NR updates.
    """
    if str(cell.solid_state_model) not in {"holdup_transport", "freeboard_closure"}:
        return False
    support = _transport_inflow_support_holdup(cell)
    K_total = np.maximum(cell.K_solid_auf + cell.K_solid_ab, 0.0)
    mask = K_total > 1e-12
    if not np.any(mask):
        return False
    underfilled_mask = mask & (np.maximum(cell.m_solid, 0.0) < support) & (support > 0.0)
    if not np.any(underfilled_mask):
        return False
    cell.m_solid[underfilled_mask] = support[underfilled_mask]
    return True


def max_upper_bed_solid_auf_flux_ratio_error(
    bed_cells: list[Cell],
    *,
    min_bed_index: int = 1,
    transport_components: tuple[int, ...] = (S_CHAR, S_ASH),
) -> float:
    """Max |1 - auf_in/out_auf| on bed1+ char/ash Eq.2.6–2.7 transport closure.

    ``out_auf`` uses ``K_auf * m_solid``; ``auf_in`` uses ``m_solid_auf_in``.

    Ref: Hamel (1999) Eq.2.6–2.7
    """
    max_err = 0.0
    for bed_idx, cell in enumerate(bed_cells):
        if bed_idx < int(min_bed_index):
            continue
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        auf_in = float(
            np.sum(np.maximum(cell.m_solid_auf_in[:, transport_components], 0.0))
        )
        out_auf = float(
            np.sum(
                np.maximum(cell.K_solid_auf[:, transport_components], 0.0)
                * np.maximum(cell.m_solid[:, transport_components], 0.0)
            )
        )
        if auf_in <= 1.0e-12 and out_auf <= 1.0e-12:
            continue
        ref = max(auf_in, out_auf, 1.0e-12)
        ratio = min(auf_in, out_auf) / ref
        max_err = max(max_err, 1.0 - ratio)
    return float(max_err)


def max_upper_bed_solid_stream_imbalance(
    bed_cells: list[Cell],
    *,
    min_bed_index: int = 1,
    transport_components: tuple[int, ...] = (S_CHAR, S_ASH),
) -> float:
    """Max |total_in - total_out| / ref on bed1+ char/ash transport streams [−].

    ``total_in`` uses ``m_solid_auf_in + m_solid_ab_in``; ``total_out`` uses
    ``_solid_outflow_rates``. When backmix (``K_ab``) carries most of the solid
    flux, auf-only closure error can be large while this stream balance stays small.

    Ref: Hamel (1999) Eq.2.6–2.7
    """
    max_err = 0.0
    for bed_idx, cell in enumerate(bed_cells):
        if bed_idx < int(min_bed_index):
            continue
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        total_in = float(
            np.sum(
                np.maximum(cell.m_solid_auf_in[:, transport_components], 0.0)
                + np.maximum(cell.m_solid_ab_in[:, transport_components], 0.0)
            )
        )
        total_out = float(
            np.sum(np.maximum(cell._solid_outflow_rates()[:, transport_components], 0.0))
        )
        ref = max(total_in, total_out, 1.0e-12)
        max_err = max(max_err, abs(total_in - total_out) / ref)
    return float(max_err)


def align_bed_solid_holdup_to_transport_inflow_support(
    cell: Cell,
    *,
    transport_components: tuple[int, ...] = (S_CHAR, S_ASH),
    allow_reduce: bool = False,
    include_negative_R: bool = False,
) -> bool:
    """Align char/ash holdup toward ``inflow / K`` transport support for upper beds.

    By default only raises under-filled inventory (init seed path). When
    ``allow_reduce`` is True, also lowers over-filled char/ash toward the same
    Eq.2.6 support — needed after outer refresh lowers upstream ``m_solid_auf_in``
    while inner NR left downstream holdup too high.

    Ref: Hamel (1999) Eq.2.6–2.7
    """
    if str(cell.solid_state_model) not in {"holdup_transport", "freeboard_closure"}:
        return False
    support = _transport_inflow_support_holdup(
        cell, include_negative_R=bool(include_negative_R)
    )
    changed = False
    for comp_idx in transport_components:
        comp_support = support[:, comp_idx]
        comp_k = np.maximum(
            cell.K_solid_auf[:, comp_idx] + cell.K_solid_ab[:, comp_idx],
            0.0,
        )
        mask = (comp_k > 1e-12) & (comp_support > 0.0)
        underfilled = mask & (np.maximum(cell.m_solid[:, comp_idx], 0.0) < comp_support - 1e-9)
        if np.any(underfilled):
            cell.m_solid[underfilled, comp_idx] = comp_support[underfilled]
            changed = True
        if allow_reduce:
            overfilled = mask & (cell.m_solid[:, comp_idx] > comp_support + 1e-9)
            if np.any(overfilled):
                cell.m_solid[overfilled, comp_idx] = comp_support[overfilled]
                changed = True
    return changed


def align_bed_solid_holdup_to_auf_upflow_closure(
    cell: Cell,
    *,
    transport_components: tuple[int, ...] = (S_CHAR, S_ASH),
    allow_reduce: bool = True,
) -> bool:
    """Set char/ash holdup so ``K_auf * m_solid`` matches ``m_solid_auf_in``.

    The full Eq.2.6 support ``inflow / (K_auf + K_ab)`` can leave a finite gap to
    the Eq.2.7 upflow term when ``K_ab > 0``. Pre-inner / init reconcile uses this
    pass after inflow refresh to close that gap explicitly.

    Ref: Hamel (1999) Eq.2.6–2.7
    """
    if str(cell.solid_state_model) not in {"holdup_transport", "freeboard_closure"}:
        return False
    changed = False
    for comp_idx in transport_components:
        k_auf = np.maximum(cell.K_solid_auf[:, comp_idx], 0.0)
        auf_in = np.maximum(cell.m_solid_auf_in[:, comp_idx], 0.0)
        mask = k_auf > 1e-12
        if not np.any(mask):
            continue
        target = np.zeros_like(k_auf)
        target[mask] = auf_in[mask] / k_auf[mask]
        current = np.maximum(cell.m_solid[:, comp_idx], 0.0)
        if allow_reduce:
            mismatch = mask & (np.abs(current - target) > 1e-9)
        else:
            mismatch = mask & (current < target - 1e-9)
        if np.any(mismatch):
            cell.m_solid[mismatch, comp_idx] = target[mismatch]
            changed = True
    return changed


def seed_bed_holdup_chain_from_transport(
    bed_cells: list[Cell],
    *,
    top_above_cell: Cell | None = None,
) -> bool:
    """Propagate thesis bed holdup seeding bottom→top from current transport states.

    Hamel Eq. 2.6 transport uses neighboring solid holdup to define axial inflow
    terms. After switching x0 to active-solid-only semantics, upper bed cells can
    remain empty unless we explicitly walk this dependency chain once during
    initialization. This sweep seeds only empty cells and then refreshes inflow
    terms from the newly created resident holdups.
    """
    if not bed_cells:
        return False
    changed_any = False
    for _ in range(max(3 * len(bed_cells), 1)):
        update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
        changed = False
        for cell in bed_cells:
            changed |= seed_holdup_from_inflows(cell)
        changed_any |= changed
        if not changed:
            break
    update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
    return changed_any


def reconcile_upper_bed_solid_holdup_for_init(
    bed_cells: list[Cell],
    cfg: Any,
    *,
    top_above_cell: Cell | None = None,
    min_align_bed_index: int = 1,
    align_passes: int = 2,
    allow_reduce: bool | None = None,
    apply_auf_upflow_closure: bool = False,
) -> bool:
    """Re-align bed1+ char/ash holdup after later init passes disturb transport support.

    Unlike ``seed_and_align_bed_holdup_chain_for_init``, this does not rerun the
    bottom→top seed sweep and is safe to call after product seed / abgleich.

    ``apply_auf_upflow_closure`` closes ``K_auf * m = m_solid_auf_in`` explicitly
    and is intended for outer pre-inner only (not init seeding).

    Ref: Hamel (1999) Eq.2.6–2.7
    """
    if not bed_cells:
        return False
    if allow_reduce is None:
        allow_reduce = bool(getattr(cfg, "nr_init_align_upper_bed_holdup_allow_reduce_thesis", False))
    include_neg_R = bool(getattr(cfg, "nr_holdup_support_include_negative_R_thesis", False))
    changed_any = False
    for pass_idx in range(max(int(align_passes), 1)):
        if not apply_frozen_bed_solid_transport_coefficients(bed_cells):
            update_bed_solid_transport_coefficients(bed_cells, cfg)
        update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
        changed = False
        for bed_idx, cell in enumerate(bed_cells):
            if bed_idx >= int(min_align_bed_index):
                if include_neg_R:
                    # 刷新 R_solid，使负反应汇进入 support
                    cell.residuals()
                changed |= align_bed_solid_holdup_to_transport_inflow_support(
                    cell,
                    allow_reduce=allow_reduce,
                    include_negative_R=include_neg_R,
                )
                if apply_auf_upflow_closure:
                    changed |= align_bed_solid_holdup_to_auf_upflow_closure(
                        cell,
                        allow_reduce=allow_reduce,
                    )
        changed_any |= changed
        if not changed and pass_idx > 0:
            break
    update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
    return changed_any


def seed_and_align_bed_holdup_chain_for_init(
    bed_cells: list[Cell],
    cfg: Any,
    *,
    top_above_cell: Cell | None = None,
    min_align_bed_index: int = 1,
    align_passes: int = 2,
) -> bool:
    """Seed bed holdup chain and reconcile bed1+ char/ash with ``K_auf`` transport support.

    The bottom→top seed sweep can leave upper beds under-filled relative to the
    hydrodynamic ``K_auf`` already frozen from Vorabrechnung. Recompute ``K_ab`` from
    updated holdups, refresh inflows, and raise bed1+ char/ash inventory toward
    ``inflow / (K_auf + K_ab)`` so ``K_auf * m_solid`` matches ``m_solid_auf_in``.

    Ref: Hamel (1999) Eq.2.6–2.7
    """
    if not bed_cells:
        return False
    changed_any = seed_bed_holdup_chain_from_transport(bed_cells, top_above_cell=top_above_cell)
    changed_any |= reconcile_upper_bed_solid_holdup_for_init(
        bed_cells,
        cfg,
        top_above_cell=top_above_cell,
        min_align_bed_index=min_align_bed_index,
        align_passes=align_passes,
    )
    return changed_any


def apply_frozen_bed_solid_transport_coefficients(bed_cells: list[Cell]) -> bool:
    """Replay frozen Vorabrechnung hydrodynamics/transport coefficients if available."""
    if not bed_cells:
        return False
    if not all(bool(cell._freeze_vorabrechnung_inner_nr) and bool(cell._vorab_solid_transport_cache_valid) for cell in bed_cells):
        return False
    for cell in bed_cells:
        cell._apply_frozen_vorabrechnung_hydrodynamics()
    return True

def apply_frozen_side_block_solid_transport_coefficients(
    cyclone_cell: Cell,
    return_leg_cell: Cell,
) -> bool:
    """Replay frozen side-block transport coefficients inside inner NR."""
    side_cells = (cyclone_cell, return_leg_cell)
    if not all(bool(cell._freeze_vorabrechnung_inner_nr) and bool(cell._vorab_solid_transport_cache_valid) for cell in side_cells):
        return False
    for cell in side_cells:
        cell._apply_frozen_vorabrechnung_hydrodynamics()
    return True

