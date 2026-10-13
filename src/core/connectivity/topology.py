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

def build_thesis_connectivity_topology(
    *,
    n_bed_cells: int,
    n_freeboard_cells: int,
    H_bed: float,
    H_freeboard: float,
    recycle_gas: bool,
    recirculation_frac: float,
    cyclone_capture_char_frac: float,
    cyclone_capture_ash_frac: float,
    secondary_injection_applied: bool,
    secondary_injection_segment: int | None,
    explicit_side_blocks: bool = False,
) -> dict:
    """Build the minimal Hamel-style connectivity graph for thesis mode."""
    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=n_bed_cells,
        n_freeboard_cells=n_freeboard_cells,
        H_bed=H_bed,
        H_freeboard=H_freeboard,
        recycle_gas=recycle_gas,
        recirculation_frac=recirculation_frac,
        cyclone_capture_char_frac=cyclone_capture_char_frac,
        cyclone_capture_ash_frac=cyclone_capture_ash_frac,
        secondary_injection_applied=secondary_injection_applied,
        secondary_injection_segment=secondary_injection_segment,
        explicit_side_blocks=explicit_side_blocks,
        preset="hamel_minimal_connectivity",
    )
    return graph.to_dict()

