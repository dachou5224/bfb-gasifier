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

from src.core.connectivity.gas_inlet import apply_bottom_recycle, apply_bottom_recycle_from_return_leg, set_bottom_cell_feeds
from src.core.connectivity.solid_transport import (
    apply_frozen_bed_solid_transport_coefficients,
    apply_frozen_side_block_solid_transport_coefficients,
    cell_total_solid_holdup,
    project_zero_source_solid_components,
    propagate_explicit_freeboard_chain,
    propagated_solid_stream_from_cell,
    update_bed_solid_transport_coefficients,
    update_bed_solid_transport_inflows,
    update_bed_solid_transport_inflows_from_above,
    update_freeboard_solid_transport_inflows,
    _set_degenerate_single_phase_gas_inlet,
)

def _apply_bed_upstream_inlet(
    prev: Cell,
    curr: Cell,
    above: Cell | None,
    cfg: Any,
    *,
    is_top_bed: bool,
) -> None:
    """Apply gas/solid upstream inlet terms for one bed cell (Hamel ``auf/ab`` semantics)."""
    curr.N_b_in[:] = prev.N_b
    curr.N_d_in[:] = prev.N_d
    curr.T_in_gas = prev.T
    xi_curr = float(curr.geo.h_center / max(float(cfg.H_bed), 1e-12))
    reactive_inlet_allowed = bool(cfg.allow_reactive_solid_propagation) and (
        xi_curr <= float(np.clip(cfg.reactive_solid_cutoff_xi, 0.0, 1.0))
    )
    if not is_top_bed and above is not None:
        lower_frac = float(np.clip(cfg.solid_lower_inlet_frac, 0.0, 1.0))
        upper_frac = 1.0 - lower_frac
        above_solid = propagated_solid_stream_from_cell(above, include_reactive=reactive_inlet_allowed)
        lower_solid = propagated_solid_stream_from_cell(prev, include_reactive=reactive_inlet_allowed)
        if reactive_inlet_allowed:
            above_reactive_empty = float(
                np.sum(np.maximum(above_solid[:, S_VM] + above_solid[:, S_MOISTURE], 0.0))
            ) <= 1e-12
            if above_reactive_empty:
                above_solid[:, S_VM] = lower_solid[:, S_VM]
                above_solid[:, S_MOISTURE] = lower_solid[:, S_MOISTURE]
        above_char_ash = float(np.sum(np.maximum(above_solid[:, S_CHAR] + above_solid[:, S_ASH], 0.0)))
        above_T = above.T if above_char_ash > 1e-12 else prev.T
        curr.m_solid_in[:] = upper_frac * above_solid + lower_frac * lower_solid
        curr.T_in_solid = upper_frac * above_T + lower_frac * prev.T
    else:
        curr.m_solid_in[:] = propagated_solid_stream_from_cell(
            prev,
            frac=float(np.clip(cfg.top_solid_inlet_frac, 0.0, 1.0)),
            include_reactive=reactive_inlet_allowed,
        )
        curr.T_in_solid = prev.T

    if str(getattr(curr, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
        vm_col = np.maximum(curr.m_solid_in[:, S_VM], 0.0)
        moist_col = np.maximum(curr.m_solid_in[:, S_MOISTURE], 0.0)
        curr.m_solid_in.fill(0.0)
        curr.m_solid_in[:, S_VM] = vm_col
        curr.m_solid_in[:, S_MOISTURE] = moist_col

def propagate_upstream(
    cells: list[Cell],
    cfg: Any,
    i: int,
    *,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Propagate axial gas and suspension-phase solid coupling into cell ``i``."""
    if i == 0:
        return

    if connectivity_graph is not None and cell_registry is not None:
        upstream_id, downstream_id = resolve_bed_neighbor_block_ids(connectivity_graph, i)
        if upstream_id is None:
            return
        prev = cell_registry.get(upstream_id)
        curr = cells[i]
        if prev is None:
            return
        above = cell_registry.get(downstream_id) if downstream_id is not None else None
        _apply_bed_upstream_inlet(prev, curr, above, cfg, is_top_bed=above is None)
        return

    prev, curr = cells[i - 1], cells[i]
    above = cells[i + 1] if i < len(cells) - 1 else None
    _apply_bed_upstream_inlet(prev, curr, above, cfg, is_top_bed=i == len(cells) - 1)

def route_side_blocks_from_connectivity_graph(
    graph: ConnectivityGraph,
    cell_registry: dict[str, Cell],
    cfg: Any,
) -> None:
    """Route cyclone/return-leg inflows using the Verbindungsmatrix edge table."""
    from src.core.connectivity_graph import ordered_solver_block_ids

    cyclone_ids = ordered_solver_block_ids(graph, "cyclone")
    return_ids = ordered_solver_block_ids(graph, "return_leg")
    if not cyclone_ids:
        return
    source_id = cyclone_inlet_source_block_id(graph)
    source_cell = resolve_registry_cell(graph, cell_registry, source_id)
    if source_cell is None:
        return

    upstream = source_cell
    for cyclone_id in cyclone_ids:
        cyclone = cell_registry.get(cyclone_id)
        if cyclone is None:
            continue
        _route_cyclone_inlet_from_upstream(upstream, cyclone, cfg)
        upstream = cyclone

    for return_id in return_ids:
        return_leg = cell_registry.get(return_id)
        if return_leg is None:
            continue
        from src.core.connectivity_graph import return_leg_upstream_block_id

        cyc_id = return_leg_upstream_block_id(graph, return_id)
        cyclone = cell_registry.get(cyc_id) if cyc_id is not None else upstream
        if cyclone is not None:
            _route_return_leg_inlet_from_cyclone(cyclone, return_leg, cfg)

def apply_recirculation_from_connectivity_graph(
    graph: ConnectivityGraph,
    bed_cells: list[Cell],
    cell_registry: dict[str, Cell],
    cfg: Any,
    *,
    relax: float | None = None,
) -> bool:
    """Apply ``N_rez`` / ``m_rez`` injection for graph-marked recirculation edges."""
    if not bed_cells:
        return False
    applied = False
    for edge in graph.edges:
        stream_class = edge.stream_class or infer_stream_class(edge.kind, edge.direction)
        if stream_class != "rez":
            continue
        bed_idx = block_id_bed_index(edge.to_block)
        if bed_idx is None:
            continue
        return_leg = resolve_registry_cell(graph, cell_registry, edge.from_block)
        if return_leg is None or str(getattr(return_leg, "cell_type", "")) != "return_leg":
            continue
        apply_bottom_recycle_from_return_leg(
            bed_cells,
            return_leg,
            cfg,
            relax=relax,
            target_bed_index=bed_idx,
        )
        applied = True
    return applied

def _route_cyclone_inlet_from_upstream(
    upstream: Cell,
    cyclone_cell: Cell,
    cfg: Any,
) -> None:
    """Route one upstream outlet into a cyclone side block."""
    _set_degenerate_single_phase_gas_inlet(upstream, cyclone_cell)
    cyclone_cell.T_in_gas = float(upstream.T)
    cyclone_cell.m_solid_in.fill(0.0)
    cyclone_cell.m_solid_auf_in[:, :] = propagated_solid_stream_from_cell(upstream, include_reactive=False)
    cyclone_cell.m_solid_ab_in.fill(0.0)
    cyclone_cell.T_in_solid = float(upstream.T)
    cyclone_cell.T_solid_auf_in = float(upstream.T)
    cyclone_cell.T_solid_ab_in = float(cyclone_cell.T)
    cyclone_cell.N_zu_b.fill(0.0)
    cyclone_cell.N_zu_d.fill(0.0)
    cyclone_cell.N_rez_b.fill(0.0)
    cyclone_cell.N_rez_d.fill(0.0)
    cyclone_cell.m_solid_zu.fill(0.0)
    cyclone_cell.m_solid_rez.fill(0.0)
    cyclone_cell.T_zu_gas = float(upstream.T)
    cyclone_cell.T_rez_gas = float(upstream.T)
    cyclone_cell.T_zu_solid = float(upstream.T)
    cyclone_cell.T_rez_solid = float(upstream.T)

def _route_return_leg_inlet_from_cyclone(
    cyclone_cell: Cell,
    return_leg_cell: Cell,
    cfg: Any,
) -> None:
    """Route cyclone underflow into one return-leg side block."""
    return_leg_cell.N_b_in.fill(0.0)
    return_leg_cell.N_d_in.fill(0.0)
    return_leg_cell.T_in_gas = float(cyclone_cell.T)
    return_leg_cell.m_solid_in.fill(0.0)
    return_leg_cell.m_solid_auf_in.fill(0.0)
    return_leg_cell.m_solid_ab_in.fill(0.0)
    return_leg_cell.T_in_solid = float(cyclone_cell.T)
    return_leg_cell.T_solid_auf_in = float(return_leg_cell.T)
    return_leg_cell.T_solid_ab_in = float(cyclone_cell.T)
    return_leg_cell.N_zu_b.fill(0.0)
    return_leg_cell.N_zu_d.fill(0.0)
    return_leg_cell.N_rez_b.fill(0.0)
    return_leg_cell.N_rez_d.fill(0.0)
    return_leg_cell.m_solid_zu.fill(0.0)
    return_leg_cell.m_solid_rez.fill(0.0)
    return_leg_cell.T_zu_gas = float(cyclone_cell.T)
    return_leg_cell.T_rez_gas = float(cyclone_cell.T)
    return_leg_cell.T_zu_solid = float(cyclone_cell.T)
    return_leg_cell.T_rez_solid = float(cyclone_cell.T)

def route_auxiliary_side_blocks(
    source_cells: list[Cell],
    cyclone_cell: Cell,
    return_leg_cell: Cell,
    cfg: Any,
) -> None:
    """Route the current top source outlet into cyclone, then cyclone underflow into return leg."""
    if not source_cells:
        return
    top = source_cells[-1]
    _route_cyclone_inlet_from_upstream(top, cyclone_cell, cfg)
    _route_return_leg_inlet_from_cyclone(cyclone_cell, return_leg_cell, cfg)

