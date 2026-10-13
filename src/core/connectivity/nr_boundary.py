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
from src.core.connectivity.routing import (
    apply_recirculation_from_connectivity_graph,
    propagate_upstream,
    route_auxiliary_side_blocks,
    route_side_blocks_from_connectivity_graph,
)
from src.core.connectivity.solid_transport import (
    apply_frozen_bed_solid_transport_coefficients,
    apply_frozen_side_block_solid_transport_coefficients,
    project_zero_source_solid_components,
    propagate_explicit_freeboard_chain,
    seed_holdup_from_inflows,
    update_bed_solid_transport_coefficients,
    update_bed_solid_transport_inflows,
    update_bed_solid_transport_inflows_from_above,
    update_freeboard_solid_transport_inflows,
    update_side_block_solid_transport_coefficients,
    _refresh_freeboard_solid_transport_inflows_for_indices,
    _set_explicit_freeboard_inlet_from_prev,
)

def _bed_transport_coefficients_are_frozen(cells: list[Cell]) -> bool:
    """Return True when bed ``K_solid`` should not be recomputed from holdup."""
    if apply_frozen_bed_solid_transport_coefficients(cells):
        return True
    from src.solvers.jacobian_transport_freeze import restore_jacobian_frozen_transport_k

    return restore_jacobian_frozen_transport_k(cells)


def _side_transport_coefficients_are_frozen(side_cells: list[Cell]) -> bool:
    if len(side_cells) >= 2 and apply_frozen_side_block_solid_transport_coefficients(
        side_cells[0],
        side_cells[1],
    ):
        return True
    from src.solvers.jacobian_transport_freeze import restore_jacobian_frozen_transport_k

    return restore_jacobian_frozen_transport_k(side_cells)

def _refresh_bed_transport_boundary_data(
    cells: list[Cell],
    cfg: Any,
    *,
    freeboard_cells: list[Cell] | None = None,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Refresh bed-side transport coefficients and axial inflow terms."""
    top_downflow_total = 0.0
    top_downflow_components: np.ndarray | None = None
    if freeboard_cells:
        top_downflow_matrix = np.maximum(freeboard_cells[0]._solid_downflow_rates(), 0.0)
        top_downflow_total = float(np.sum(top_downflow_matrix))
        top_downflow_components = np.array(
            [
                float(np.sum(top_downflow_matrix[:, S_CHAR])),
                float(np.sum(top_downflow_matrix[:, S_ASH])),
            ],
            dtype=np.float64,
        )
    if not _bed_transport_coefficients_are_frozen(cells):
        update_bed_solid_transport_coefficients(
            cells,
            cfg,
            top_downflow_total=top_downflow_total,
            top_downflow_components=top_downflow_components,
        )
    update_bed_solid_transport_inflows_from_above(
        cells,
        freeboard_cells[0] if freeboard_cells else None,
        connectivity_graph=connectivity_graph,
        cell_registry=cell_registry,
    )

def _refresh_side_block_boundary_data(
    cells: list[Cell],
    cfg: Any,
    *,
    freeboard_cells: list[Cell],
    cyclone_cell: Cell,
    return_leg_cell: Cell,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Refresh explicit cyclone/return-leg routing and bottom recycle."""
    from src.core.connectivity_graph import ordered_solver_block_ids, return_leg_upstream_block_id

    if connectivity_graph is not None and cell_registry is not None:
        route_side_blocks_from_connectivity_graph(connectivity_graph, cell_registry, cfg)
        cyclone_ids = ordered_solver_block_ids(connectivity_graph, "cyclone")
        return_ids = ordered_solver_block_ids(connectivity_graph, "return_leg")
        cyclone_cells = [cell_registry[cid] for cid in cyclone_ids if cid in cell_registry]
        return_cells = [cell_registry[rid] for rid in return_ids if rid in cell_registry]
    else:
        cyclone_cells = [cyclone_cell]
        return_cells = [return_leg_cell]

    side_transport_frozen = _side_transport_coefficients_are_frozen(
        [cyclone_cells[0] if cyclone_cells else cyclone_cell, return_cells[0] if return_cells else return_leg_cell]
    )
    for cyc in cyclone_cells:
        cyc.N_d[:] = np.maximum(cyc.N_d_in, 0.0)
        cyc.N_b.fill(0.0)
        cyc.T = float(cyc.T_in_gas)

    last_cyclone = cyclone_cells[-1] if cyclone_cells else cyclone_cell
    if not side_transport_frozen:
        for leg in return_cells:
            update_side_block_solid_transport_coefficients(last_cyclone, leg, cfg)
        for cyc in cyclone_cells:
            seed_holdup_from_inflows(cyc)

    for leg in return_cells:
        cyc_src = last_cyclone
        if connectivity_graph is not None and cell_registry is not None:
            leg_id = getattr(leg, "_side_block_graph_id", None)
            if leg_id is None:
                for block_id, cell in cell_registry.items():
                    if cell is leg:
                        leg_id = block_id
                        break
            if leg_id is not None:
                upstream_id = return_leg_upstream_block_id(connectivity_graph, leg_id)
                resolved = cell_registry.get(upstream_id) if upstream_id else None
                if resolved is not None:
                    cyc_src = resolved
        if not side_transport_frozen:
            update_side_block_solid_transport_coefficients(cyc_src, leg, cfg)
        leg.m_solid_ab_in[:, :] = np.maximum(cyc_src._solid_downflow_rates(), 0.0)
        leg.T_in_solid = float(cyc_src.T)
        leg.T_solid_ab_in = float(cyc_src.T)
        leg.T = float(cyc_src.T)
        seed_holdup_from_inflows(leg)

    if connectivity_graph is not None and cell_registry is not None:
        if not apply_recirculation_from_connectivity_graph(
            connectivity_graph,
            cells,
            cell_registry,
            cfg,
            relax=None,
        ):
            apply_bottom_recycle_from_return_leg(cells, return_leg_cell, cfg, relax=None)
    else:
        apply_bottom_recycle_from_return_leg(cells, return_leg_cell, cfg, relax=None)
    for cyc in cyclone_cells:
        project_zero_source_solid_components(cyc)
    for leg in return_cells:
        project_zero_source_solid_components(leg)

def apply_all_nr_boundary_data(
    cells: list[Cell],
    cfg: Any,
    *,
    freeboard_cells: list[Cell] | None = None,
    cyclone_cell: Cell | None = None,
    return_leg_cell: Cell | None = None,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Apply the full thesis-mode axial routing / recycle boundary update."""
    if not cells:
        return
    freeboard_cells = freeboard_cells or []
    set_bottom_cell_feeds(cells, cfg)
    cells[-1].m_solid_in[:] = 0.0
    for i in range(len(cells)):
        propagate_upstream(
            cells,
            cfg,
            i,
            connectivity_graph=connectivity_graph,
            cell_registry=cell_registry,
        )
    if freeboard_cells:
        propagate_explicit_freeboard_chain(
            cells,
            freeboard_cells,
            connectivity_graph=connectivity_graph,
            cell_registry=cell_registry,
        )
    _refresh_bed_transport_boundary_data(
        cells,
        cfg,
        freeboard_cells=freeboard_cells,
        connectivity_graph=connectivity_graph,
        cell_registry=cell_registry,
    )
    if freeboard_cells:
        update_freeboard_solid_transport_inflows(
            freeboard_cells,
            bottom_below_cell=cells[-1],
            connectivity_graph=connectivity_graph,
            cell_registry=cell_registry,
        )
    if cyclone_cell is not None and return_leg_cell is not None:
        _refresh_side_block_boundary_data(
            cells,
            cfg,
            freeboard_cells=freeboard_cells,
            cyclone_cell=cyclone_cell,
            return_leg_cell=return_leg_cell,
            connectivity_graph=connectivity_graph,
            cell_registry=cell_registry,
        )
    else:
        apply_bottom_recycle(cells, cfg, relax=None)
    for cell in cells:
        project_zero_source_solid_components(cell)
    for cell in freeboard_cells:
        project_zero_source_solid_components(cell)

def apply_local_nr_boundary_data(
    cells: list[Cell],
    cfg: Any,
    changed_cell_idx: int,
    *,
    freeboard_cells: list[Cell] | None = None,
    cyclone_cell: Cell | None = None,
    return_leg_cell: Cell | None = None,
    connectivity_graph: ConnectivityGraph | None = None,
    cell_registry: dict[str, Cell] | None = None,
) -> None:
    """Apply only the BC updates touched by a single perturbed cell."""
    n_cells = len(cells)
    if n_cells == 0:
        return
    if freeboard_cells or cyclone_cell is not None or return_leg_cell is not None:
        freeboard_cells = freeboard_cells or []
        idx = int(changed_cell_idx)
        n_bed = len(cells)
        n_fb = len(freeboard_cells)
        touched_bed: set[int] = set()
        touched_fb: set[int] = set()
        touched_side = False

        set_bottom_cell_feeds(cells, cfg)
        cells[-1].m_solid_in[:] = 0.0

        if idx < n_bed:
            touched_bed.add(idx)
            if idx > 0:
                touched_bed.add(idx - 1)
            if idx + 1 < n_bed:
                propagate_upstream(
                    cells,
                    cfg,
                    idx + 1,
                    connectivity_graph=connectivity_graph,
                    cell_registry=cell_registry,
                )
                touched_bed.add(idx + 1)
            if idx == n_bed - 1:
                if n_fb > 0:
                    _set_explicit_freeboard_inlet_from_prev(cells[-1], freeboard_cells[0])
                    _refresh_freeboard_solid_transport_inflows_for_indices(
                        freeboard_cells,
                        bottom_below_cell=cells[-1],
                        indices={0},
                        connectivity_graph=connectivity_graph,
                        cell_registry=cell_registry,
                    )
                    touched_fb.add(0)
                elif cyclone_cell is not None and return_leg_cell is not None:
                    touched_side = True
            _refresh_bed_transport_boundary_data(
                cells,
                cfg,
                freeboard_cells=freeboard_cells,
                connectivity_graph=connectivity_graph,
                cell_registry=cell_registry,
            )
        elif idx < n_bed + n_fb:
            fb_idx = idx - n_bed
            touched_fb.add(fb_idx)
            if fb_idx > 0:
                _refresh_freeboard_solid_transport_inflows_for_indices(
                    freeboard_cells,
                    bottom_below_cell=cells[-1],
                    indices={fb_idx - 1},
                    connectivity_graph=connectivity_graph,
                    cell_registry=cell_registry,
                )
                touched_fb.add(fb_idx - 1)
            if fb_idx + 1 < n_fb:
                _set_explicit_freeboard_inlet_from_prev(freeboard_cells[fb_idx], freeboard_cells[fb_idx + 1])
                _refresh_freeboard_solid_transport_inflows_for_indices(
                    freeboard_cells,
                    bottom_below_cell=cells[-1],
                    indices={fb_idx + 1},
                    connectivity_graph=connectivity_graph,
                    cell_registry=cell_registry,
                )
                touched_fb.add(fb_idx + 1)
            if fb_idx == 0:
                _refresh_bed_transport_boundary_data(
                cells,
                cfg,
                freeboard_cells=freeboard_cells,
                connectivity_graph=connectivity_graph,
                cell_registry=cell_registry,
            )
                touched_bed.update(range(n_bed))
            if fb_idx == n_fb - 1 and cyclone_cell is not None and return_leg_cell is not None:
                touched_side = True
        elif cyclone_cell is not None and idx == n_bed + n_fb:
            touched_side = True
        elif return_leg_cell is not None and idx == n_bed + n_fb + 1:
            if connectivity_graph is not None and cell_registry is not None:
                apply_recirculation_from_connectivity_graph(
                    connectivity_graph,
                    cells,
                    cell_registry,
                    cfg,
                    relax=None,
                )
            else:
                apply_bottom_recycle_from_return_leg(cells, return_leg_cell, cfg, relax=None)
            touched_bed.add(0)
            if n_bed > 1:
                touched_bed.add(1)
            project_zero_source_solid_components(return_leg_cell)

        if n_fb > 0 and not touched_fb and idx < n_bed and idx != n_bed - 1:
            # Keep freeboard transport coefficients/inflows coherent only when explicit
            # freeboard cells are already active in the solver graph.
            _refresh_freeboard_solid_transport_inflows_for_indices(
                freeboard_cells,
                bottom_below_cell=cells[-1],
                indices=set(),
                connectivity_graph=connectivity_graph,
                cell_registry=cell_registry,
            )
        if touched_side and cyclone_cell is not None and return_leg_cell is not None:
            _refresh_side_block_boundary_data(
                cells,
                cfg,
                freeboard_cells=freeboard_cells,
                cyclone_cell=cyclone_cell,
                return_leg_cell=return_leg_cell,
                connectivity_graph=connectivity_graph,
                cell_registry=cell_registry,
            )
            touched_bed.add(0)
            if n_bed > 1:
                touched_bed.add(1)
            if n_fb > 0:
                touched_fb.add(n_fb - 1)

        for bed_idx in sorted(touched_bed):
            if 0 <= bed_idx < n_bed:
                project_zero_source_solid_components(cells[bed_idx])
        for fb_idx in sorted(touched_fb):
            if 0 <= fb_idx < n_fb:
                project_zero_source_solid_components(freeboard_cells[fb_idx])
        return

    set_bottom_cell_feeds(cells, cfg)
    cells[-1].m_solid_in[:] = 0.0
    apply_bottom_recycle(cells, cfg, relax=None)
    if int(changed_cell_idx) == 0:
        for i in range(len(cells)):
            propagate_upstream(
                cells,
                cfg,
                i,
                connectivity_graph=connectivity_graph,
                cell_registry=cell_registry,
            )
        if not _bed_transport_coefficients_are_frozen(cells):
            update_bed_solid_transport_coefficients(cells, cfg)
        update_bed_solid_transport_inflows(cells)
        for cell in cells:
            project_zero_source_solid_components(cell)
        return
    project_zero_source_solid_components(cells[0])

    touched = {int(changed_cell_idx)}
    if int(changed_cell_idx) - 1 > 0:
        touched.add(int(changed_cell_idx) - 1)
    if int(changed_cell_idx) + 1 < n_cells:
        touched.add(int(changed_cell_idx) + 1)
    if int(changed_cell_idx) == n_cells - 1:
        touched.add(0)

    for idx in sorted(touched):
        if idx > 0:
            propagate_upstream(
                cells,
                cfg,
                idx,
                connectivity_graph=connectivity_graph,
                cell_registry=cell_registry,
            )
        project_zero_source_solid_components(cells[idx])
    if not _bed_transport_coefficients_are_frozen(cells):
        update_bed_solid_transport_coefficients(cells, cfg)
    update_bed_solid_transport_inflows(cells)

