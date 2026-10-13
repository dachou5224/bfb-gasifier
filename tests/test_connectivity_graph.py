from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from src.core.connectivity import build_thesis_connectivity_topology
from src.core.connectivity_graph import (
    affected_nr_residual_cells,
    build_connectivity_from_reactor_config,
    build_connectivity_graph_from_reactor_params,
    connectivity_graph_from_dict,
    infer_stream_class,
    load_connectivity_graph_json,
)
from src.core.reactor import Reactor, ReactorConfig


def test_affected_cells_simple_bed_chain():
    assert affected_nr_residual_cells(
        changed_cell_idx=2,
        n_bed=4,
        n_freeboard=0,
        explicit_freeboard_graph=False,
        side_blocks_in_boundary_path=False,
        has_side_blocks=False,
    ) == (1, 2, 3)


def test_affected_cells_bed_with_side_blocks():
    assert affected_nr_residual_cells(
        changed_cell_idx=3,
        n_bed=4,
        n_freeboard=0,
        explicit_freeboard_graph=False,
        side_blocks_in_boundary_path=True,
        has_side_blocks=True,
    ) == (2, 3, 4)


def test_affected_cells_explicit_freeboard_graph_with_recycle():
    assert affected_nr_residual_cells(
        changed_cell_idx=5,
        n_bed=3,
        n_freeboard=1,
        explicit_freeboard_graph=True,
        side_blocks_in_boundary_path=True,
        has_side_blocks=True,
    ) == (0, 1, 2, 5)


def test_infer_stream_class_maps_recirculation():
    assert infer_stream_class("external_zirkulation", "downward_to_bottom") == "rez"
    assert infer_stream_class("side_feed", "into_reactor") == "zu"
    assert infer_stream_class("axial_upflow", "upward") == "auf"


def test_build_vtt_graph_has_cross_links_and_solver_indices():
    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=3,
        n_freeboard_cells=1,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
        preset="vtt_minimal_loop",
    )
    assert graph.n_cross_links == 2
    assert graph.freeboard_active is True
    assert graph.recycle_active is True
    solver_ids = graph.solver_cell_block_ids()
    assert solver_ids[0] == "bed_cell_0"
    assert "cyclone_block" in solver_ids
    assert "return_leg_block" in solver_ids


def test_load_vtt_json_roundtrip():
    path = Path(__file__).resolve().parents[1] / "data" / "connectivity" / "vtt_minimal_loop.json"
    graph = load_connectivity_graph_json(path)
    assert graph.mode == "vtt_minimal_loop"
    assert len(graph.blocks) >= 8
    assert graph.n_cross_links == 2
    reloaded = connectivity_graph_from_dict(json.loads(path.read_text(encoding="utf-8")))
    assert reloaded.to_dict()["counts"] == graph.to_dict()["counts"]


def test_build_thesis_topology_matches_graph_builder():
    topo = build_thesis_connectivity_topology(
        n_bed_cells=4,
        n_freeboard_cells=2,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
    )
    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=4,
        n_freeboard_cells=2,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
        preset="hamel_minimal_connectivity",
    )
    assert topo["mode"] == graph.mode
    assert topo["counts"]["n_cross_links"] == graph.n_cross_links
    assert len(topo["edges"]) == len(graph.edges)


def test_reactor_builds_connectivity_graph_on_init():
    cfg = ReactorConfig(n_cells=2, thesis_mode=True, n_freeboard_cells=0, H_freeboard=0.0)
    reactor = Reactor(cfg)
    assert reactor.connectivity_graph is not None
    assert reactor.connectivity_graph.mode == "hamel_minimal_connectivity"
    assert reactor.connectivity_graph.blocks[2].block_id == "bed_cell_0"


def test_reactor_loads_connectivity_from_json_path():
    path = Path(__file__).resolve().parents[1] / "data" / "connectivity" / "vtt_minimal_loop.json"
    cfg = ReactorConfig(
        n_cells=3,
        thesis_mode=True,
        n_freeboard_cells=1,
        H_freeboard=9.5,
        connectivity_graph_path=str(path),
    )
    reactor = Reactor(cfg)
    assert reactor.connectivity_graph.mode == "vtt_minimal_loop"
    assert reactor.connectivity_graph.n_cross_links == 2


def test_affected_cells_from_connectivity_graph():
    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=3,
        n_freeboard_cells=1,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
    )
    cyclone_idx = next(
        block.solver_cell_index
        for block in graph.blocks
        if block.block_id == "cyclone_block"
    )
    neighbors = affected_nr_residual_cells(
        changed_cell_idx=int(cyclone_idx),
        n_bed=3,
        n_freeboard=1,
        explicit_freeboard_graph=True,
        side_blocks_in_boundary_path=True,
        has_side_blocks=True,
        connectivity_graph=graph,
    )
    assert int(cyclone_idx) in neighbors
    assert 0 in neighbors
    assert 3 in neighbors  # freeboard_cell_0
    assert 5 in neighbors  # return_leg_block


def test_jacobian_side_pairs_from_graph_match_expected_vtt_wraparound():
    from src.core.connectivity_graph import jacobian_side_block_pairs_from_connectivity_graph

    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=3,
        n_freeboard_cells=2,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
        preset="vtt_minimal_loop",
    )
    pairs = jacobian_side_block_pairs_from_connectivity_graph(graph)
    assert (0, 3) in pairs
    assert (0, 6) in pairs
    assert (6, 4) in pairs
    assert graph.n_cross_links == 2


def test_build_jacobian_structure_side_pairs_match_graph_derivation():
    """P0-2: declared side blocks must equal connectivity-graph Nebenelemente pairs."""
    from src.core.connectivity_graph import jacobian_side_block_pairs_from_connectivity_graph
    from src.solvers.structured_jacobian import build_jacobian_structure, build_solver_graph

    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=3,
        n_freeboard_cells=1,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
        preset="vtt_minimal_loop",
    )
    cfg = ReactorConfig(
        n_cells=3,
        thesis_mode=True,
        n_freeboard_cells=1,
        H_freeboard=9.5,
        explicit_side_blocks_enabled=True,
        connectivity_graph_path=str(
            Path(__file__).resolve().parents[1] / "data" / "connectivity" / "vtt_minimal_loop.json"
        ),
    )
    reactor = Reactor(cfg)
    solver_cells = reactor._solver_cells_for_nr()
    solver_graph = build_solver_graph(solver_cells)
    structure = build_jacobian_structure(
        solver_graph,
        solver_cells,
        connectivity_graph=reactor.connectivity_graph,
    )
    expected = jacobian_side_block_pairs_from_connectivity_graph(reactor.connectivity_graph)
    assert structure.side_block_pairs == expected
    assert len(structure.side_block_pairs) > 0


def test_invalid_graph_rejects_unknown_edge_endpoint():
    with pytest.raises(ValueError, match="unknown block"):
        connectivity_graph_from_dict(
            {
                "blocks": [{"block_id": "a", "kind": "bed_cell"}],
                "edges": [
                    {
                        "from_block": "a",
                        "to_block": "missing",
                        "kind": "axial_upflow",
                        "direction": "upward",
                        "solver_coupling": "explicit_state",
                    }
                ],
            }
        )


def test_build_connectivity_from_reactor_config_default():
    cfg = ReactorConfig(n_cells=3, thesis_mode=True, explicit_side_blocks_enabled=True)
    graph = build_connectivity_from_reactor_config(cfg)
    assert graph.mode == "hamel_minimal_connectivity"


def test_bed_axial_neighbors_from_graph_match_linear_chain():
    from src.core.connectivity_graph import resolve_bed_neighbor_block_ids

    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=4,
        n_freeboard_cells=0,
        H_bed=5.0,
        H_freeboard=0.0,
        recycle_gas=False,
        recirculation_frac=0.0,
        cyclone_capture_char_frac=0.0,
        cyclone_capture_ash_frac=0.0,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=False,
    )
    assert resolve_bed_neighbor_block_ids(graph, 0) == (None, "bed_cell_1")
    assert resolve_bed_neighbor_block_ids(graph, 2) == ("bed_cell_1", "bed_cell_3")
    assert resolve_bed_neighbor_block_ids(graph, 3) == ("bed_cell_2", None)


def test_propagate_upstream_graph_path_matches_index_path():
    from src.core.connectivity import propagate_upstream
    from src.core.connectivity_graph import build_solver_cell_registry

    cfg = ReactorConfig(n_cells=3, solid_lower_inlet_frac=0.25, top_solid_inlet_frac=0.5)
    reactor = Reactor(cfg)
    registry = build_solver_cell_registry(reactor)
    graph = reactor.connectivity_graph
    assert graph is not None

    index_cells = [cell for cell in reactor.cells]
    graph_cells = [cell for cell in reactor.cells]
    for cell in index_cells + graph_cells:
        cell.m_solid.fill(0.0)
    for prev, above in ((index_cells[0], index_cells[2]), (graph_cells[0], graph_cells[2])):
        prev.m_solid[0, :] = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
        above.m_solid[0, :] = np.array([5.0, 6.0, 7.0, 8.0], dtype=float)

    propagate_upstream(index_cells, cfg, 1)
    propagate_upstream(
        graph_cells,
        cfg,
        1,
        connectivity_graph=graph,
        cell_registry=registry,
    )
    np.testing.assert_allclose(graph_cells[1].N_d_in, index_cells[1].N_d_in)
    np.testing.assert_allclose(graph_cells[1].N_b_in, index_cells[1].N_b_in)
    np.testing.assert_allclose(graph_cells[1].m_solid_in, index_cells[1].m_solid_in)


def test_bed_solid_transport_graph_path_matches_index_path():
    from src.core.cell import S_CHAR
    from src.core.connectivity import update_bed_solid_transport_inflows_from_above
    from src.core.connectivity_graph import build_solver_cell_registry

    cfg = ReactorConfig(
        n_cells=3,
        thesis_mode=True,
        H_freeboard=2.0,
        n_freeboard_cells=1,
    )
    reactor = Reactor(cfg)
    registry = build_solver_cell_registry(reactor)
    graph = reactor.connectivity_graph
    assert graph is not None

    beds = list(reactor.cells)
    fb = reactor.freeboard_cells[0]
    fb.solid_state_model = "freeboard_closure"
    fb.m_solid[0, S_CHAR] = 9.0
    fb.K_solid_ab[0, S_CHAR] = 0.35
    fb.T = 1180.0

    for i, cell in enumerate(beds):
        cell.solid_state_model = "holdup_transport"
        cell.m_solid[0, S_CHAR] = 1.0 + float(i)
        cell.K_solid_auf.fill(0.12)
        cell.K_solid_ab.fill(0.21)
        cell.T = 1100.0 + float(i)

    expected_auf = [np.zeros_like(cell.m_solid_auf_in) for cell in beds]
    expected_ab = [np.zeros_like(cell.m_solid_ab_in) for cell in beds]
    update_bed_solid_transport_inflows_from_above(beds, top_above_cell=fb)
    for i, cell in enumerate(beds):
        expected_auf[i][:, :] = cell.m_solid_auf_in
        expected_ab[i][:, :] = cell.m_solid_ab_in

    update_bed_solid_transport_inflows_from_above(
        beds,
        top_above_cell=fb,
        connectivity_graph=graph,
        cell_registry=registry,
    )
    for i, cell in enumerate(beds):
        np.testing.assert_allclose(cell.m_solid_auf_in, expected_auf[i])
        np.testing.assert_allclose(cell.m_solid_ab_in, expected_ab[i])


def test_freeboard_solid_transport_graph_path_matches_index_path():
    from src.core.cell import S_CHAR
    from src.core.connectivity import update_freeboard_solid_transport_inflows
    from src.core.connectivity_graph import build_solver_cell_registry

    cfg = ReactorConfig(
        n_cells=2,
        thesis_mode=True,
        H_freeboard=2.0,
        n_freeboard_cells=2,
    )
    reactor = Reactor(cfg)
    registry = build_solver_cell_registry(reactor)
    graph = reactor.connectivity_graph
    assert graph is not None

    bed_top = reactor.cells[-1]
    fb0, fb1 = reactor.freeboard_cells
    bed_top.solid_state_model = "holdup_transport"
    bed_top.m_solid[0, S_CHAR] = 1.0
    bed_top.K_solid_auf[0, S_CHAR] = 0.4
    fb0.solid_state_model = "holdup_transport"
    fb1.solid_state_model = "holdup_transport"
    fb0.m_solid[0, S_CHAR] = 2.0
    fb0.K_solid_ab[0, S_CHAR] = 0.3
    fb1.m_solid[0, S_CHAR] = 3.0
    fb1.K_solid_auf[0, S_CHAR] = 0.2

    expected_auf = [np.zeros_like(fb0.m_solid_auf_in), np.zeros_like(fb1.m_solid_auf_in)]
    expected_ab = [np.zeros_like(fb0.m_solid_ab_in), np.zeros_like(fb1.m_solid_ab_in)]
    update_freeboard_solid_transport_inflows([fb0, fb1], bottom_below_cell=bed_top)
    expected_auf[0][:, :] = fb0.m_solid_auf_in
    expected_ab[0][:, :] = fb0.m_solid_ab_in
    expected_auf[1][:, :] = fb1.m_solid_auf_in
    expected_ab[1][:, :] = fb1.m_solid_ab_in

    update_freeboard_solid_transport_inflows(
        [fb0, fb1],
        bottom_below_cell=bed_top,
        connectivity_graph=graph,
        cell_registry=registry,
    )
    np.testing.assert_allclose(fb0.m_solid_auf_in, expected_auf[0])
    np.testing.assert_allclose(fb0.m_solid_ab_in, expected_ab[0])
    np.testing.assert_allclose(fb1.m_solid_auf_in, expected_auf[1])
    np.testing.assert_allclose(fb1.m_solid_ab_in, expected_ab[1])


def test_connection_pipe_passthrough_resolves_cyclone_inlet_to_physical_freeboard():
    from src.core.connectivity_graph import (
        add_connection_pipe_passthrough,
        cyclone_inlet_source_block_id,
    )

    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=3,
        n_freeboard_cells=1,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
        preset="vtt_minimal_loop",
    )
    assert cyclone_inlet_source_block_id(graph) == "freeboard_cell_0"

    piped = add_connection_pipe_passthrough(graph)
    assert cyclone_inlet_source_block_id(piped) == "freeboard_cell_0"
    assert any(block.block_id == "connection_pipe_0" for block in piped.blocks)

    json_path = Path(__file__).resolve().parents[1] / "data" / "connectivity" / "vtt_connection_pipe_passthrough.json"
    loaded = load_connectivity_graph_json(json_path)
    assert cyclone_inlet_source_block_id(loaded) == "freeboard_cell_0"


def test_dual_cyclone_dual_return_graph_and_reactor_side_blocks():
    from src.core.connectivity_graph import ordered_solver_block_ids, solver_topology_from_graph

    graph = build_connectivity_graph_from_reactor_params(
        n_bed_cells=3,
        n_freeboard_cells=1,
        H_bed=5.0,
        H_freeboard=9.5,
        recycle_gas=True,
        recirculation_frac=0.1,
        cyclone_capture_char_frac=0.9,
        cyclone_capture_ash_frac=0.95,
        secondary_injection_applied=False,
        secondary_injection_segment=None,
        explicit_side_blocks=True,
        preset="vtt_dual_return",
        n_cyclone_stages=2,
        n_return_legs=2,
        return_leg_bed_targets=(0, 2),
    )
    assert ordered_solver_block_ids(graph, "cyclone") == ("cyclone_block_0", "cyclone_block_1")
    assert ordered_solver_block_ids(graph, "return_leg") == ("return_leg_block_0", "return_leg_block_1")
    topo = solver_topology_from_graph(graph)
    assert len(topo.cyclone_indices) == 2
    assert len(topo.return_leg_indices) == 2

    cfg = ReactorConfig(
        n_cells=3,
        thesis_mode=True,
        n_freeboard_cells=1,
        H_freeboard=9.5,
        n_cyclone_stages_thesis=2,
        n_return_legs_thesis=2,
        return_leg_bed_targets_thesis=(0, 2),
    )
    reactor = Reactor(cfg)
    assert len(reactor.cyclone_cells) == 2
    assert len(reactor.return_leg_cells) == 2
    assert len(reactor.side_cells) == 4
    assert reactor.side_block_cells_by_id["return_leg_block_1"] is reactor.return_leg_cells[1]


def test_dual_return_recycle_targets_top_bed_from_graph():
    from src.core.cell import S_ASH, S_CHAR
    from src.core.connectivity import apply_recirculation_from_connectivity_graph
    from src.core.connectivity_graph import build_solver_cell_registry

    cfg = ReactorConfig(
        n_cells=3,
        thesis_mode=True,
        n_freeboard_cells=0,
        n_return_legs_thesis=2,
        return_leg_bed_targets_thesis=(0, 2),
        recirculation_frac=0.2,
    )
    reactor = Reactor(cfg)
    registry = build_solver_cell_registry(reactor)
    graph = reactor.connectivity_graph
    assert graph is not None

    leg0, leg1 = reactor.return_leg_cells
    for leg in (leg0, leg1):
        leg.K_solid_ab[:, S_CHAR] = 0.5
        leg.m_solid[0, S_CHAR] = 2.0
        leg.m_solid[0, S_ASH] = 1.0

    apply_recirculation_from_connectivity_graph(graph, reactor.cells, registry, cfg, relax=None)
    assert float(reactor.cells[0].m_solid_rez[0, S_CHAR]) > 0.0
    assert float(reactor.cells[2].m_solid_rez[0, S_CHAR]) > 0.0


def test_wirsum_damping_params_from_config():
    from src.solvers.wirsum_damping import wirsum_damped_newton_params_from_config

    cfg = ReactorConfig(
        thesis_mode=True,
        nr_damping_halvings_thesis=18,
        nr_lambda_init_thesis=0.5,
        nr_lambda_min_thesis=1.0 / 2048.0,
        nr_prefer_full_step_thesis=False,
    )
    params = wirsum_damped_newton_params_from_config(cfg)
    assert params.n_damp_halvings == 18
    assert params.lambda_init == pytest.approx(0.5)
    assert params.lambda_min == pytest.approx(1.0 / 2048.0)
    assert params.prefer_full_step is False
