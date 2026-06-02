from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from src.core.connectivity import (
    _set_explicit_freeboard_inlet_from_prev,
    cell_total_solid_holdup,
    seed_holdup_from_inflows,
    update_bed_solid_transport_coefficients,
    update_bed_solid_transport_inflows,
    update_freeboard_solid_transport_inflows,
)
from src.core.cell import S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.freeboard_segment import simulate_freeboard
from src.core.reactor import GAS_SPECIES, GAS_SPECIES_INDEX, Reactor, ReactorConfig, _propagated_solid_stream, _recycled_solid_stream
from src.solvers.vorabrechnung import vorabrechnung_tau_for_cell


def _load_char_audit_summary_fn():
    path = Path(__file__).resolve().parent.parent / "scripts" / "audit_char_mass_conservation_lu.py"
    spec = importlib.util.spec_from_file_location("audit_char_mass_conservation_lu", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.summarize_reactor_char_ledger


def test_recycled_solid_stream_keeps_only_char_and_ash():
    m = np.array([[1.0, 2.0, 3.0, 4.0]], dtype=float)
    out = _recycled_solid_stream(m, frac=0.25)

    assert out[0, S_CHAR] == 0.25
    assert out[0, S_ASH] == 1.0
    assert out[0, S_VM] == 0.0
    assert out[0, S_MOISTURE] == 0.0


def test_propagated_solid_stream_keeps_only_char_and_ash():
    m = np.array([[1.0, 2.0, 3.0, 4.0]], dtype=float)
    out = _propagated_solid_stream(m, frac=0.5)

    assert out[0, S_CHAR] == 0.5
    assert out[0, S_ASH] == 2.0
    assert out[0, S_VM] == 0.0
    assert out[0, S_MOISTURE] == 0.0


def test_propagated_solid_stream_can_keep_reactive_components_when_enabled():
    m = np.array([[1.0, 2.0, 3.0, 4.0]], dtype=float)
    out = _propagated_solid_stream(m, frac=0.5, include_reactive=True)

    assert out[0, S_CHAR] == 0.5
    assert out[0, S_VM] == 1.0
    assert out[0, S_MOISTURE] == 1.5
    assert out[0, S_ASH] == 2.0


def test_propagate_upstream_filters_vm_and_moisture_from_solid_inlet():
    reactor = Reactor(ReactorConfig(n_cells=3, solid_lower_inlet_frac=0.25, top_solid_inlet_frac=0.5))
    prev = reactor.cells[0]
    curr = reactor.cells[1]
    above = reactor.cells[2]

    prev.m_solid[0, :] = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    above.m_solid[0, :] = np.array([5.0, 6.0, 7.0, 8.0], dtype=float)

    reactor._propagate_upstream(1)

    np.testing.assert_allclose(
        curr.m_solid_in[0, :],
        np.array(
            [
                0.75 * 5.0 + 0.25 * 1.0,
                0.0,
                0.0,
                0.75 * 8.0 + 0.25 * 4.0,
            ],
            dtype=float,
        ),
    )


def test_propagate_upstream_allows_reactive_components_in_lower_zone_only():
    reactor = Reactor(
        ReactorConfig(
            n_cells=4,
            H_bed=4.0,
            solid_lower_inlet_frac=0.0,
            allow_reactive_solid_propagation=True,
            reactive_solid_cutoff_xi=0.35,
        )
    )
    # cell 1 center xi=0.375 (>0.35), cell 0->1 should already be filtered.
    reactor.cells[0].m_solid[0, :] = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    reactor._propagate_upstream(1)
    assert reactor.cells[1].m_solid_in[0, S_VM] == 0.0
    assert reactor.cells[1].m_solid_in[0, S_MOISTURE] == 0.0

    # Raise cutoff so lower cells can carry reactive stock.
    reactor.config.reactive_solid_cutoff_xi = 0.5
    reactor._propagate_upstream(1)
    assert reactor.cells[1].m_solid_in[0, S_VM] > 0.0
    assert reactor.cells[1].m_solid_in[0, S_MOISTURE] > 0.0


def test_propagate_upstream_allows_reactive_components_in_holdup_mode_lower_zone():
    reactor = Reactor(
        ReactorConfig(
            n_cells=4,
            H_bed=4.0,
            thesis_mode=True,
            solid_lower_inlet_frac=0.0,
            reactive_solid_cutoff_xi=0.5,
        )
    )
    reactor.cells[0].m_solid[0, S_VM] = 2.0
    reactor.cells[0].m_solid[0, S_MOISTURE] = 3.0
    reactor.cells[2].m_solid[0, :] = np.array([5.0, 0.0, 0.0, 6.0], dtype=float)
    reactor.cells[2].K_solid_auf[0, S_CHAR] = 0.2
    reactor.cells[2].K_solid_auf[0, S_ASH] = 0.2

    reactor._propagate_upstream(1)

    np.testing.assert_allclose(
        reactor.cells[1].m_solid_in[0, :],
        np.array([0.0, 2.0, 3.0, 0.0], dtype=float),
    )


def test_holdup_mode_reactive_propagation_does_not_duplicate_fresh_support():
    reactor = Reactor(
        ReactorConfig(
            n_cells=3,
            H_bed=3.0,
            thesis_mode=True,
            solid_lower_inlet_frac=1.0,
            reactive_solid_cutoff_xi=1.0,
        )
    )
    reactor.cells[0].m_solid_zu[0, S_VM] = 2.0
    reactor.cells[0].m_solid_zu[0, S_MOISTURE] = 3.0
    reactor.cells[0].m_solid[0, S_VM] = 0.5
    reactor.cells[0].m_solid[0, S_MOISTURE] = 0.25

    reactor._propagate_upstream(1)

    assert reactor.cells[1].m_solid_in[0, S_VM] == 0.5
    assert reactor.cells[1].m_solid_in[0, S_MOISTURE] == 0.25


def test_reactor_init_applies_thesis_reactive_defaults_when_mode_enabled_late():
    cfg = ReactorConfig(n_cells=2)
    cfg.thesis_mode = True
    cfg.allow_reactive_solid_propagation = False
    cfg.reactive_solid_cutoff_xi = 1.2
    cfg.major_gibbs_solver_mode = "augmented"

    reactor = Reactor(cfg)

    assert reactor.config.allow_reactive_solid_propagation is True
    assert reactor.config.reactive_solid_cutoff_xi == 1.0
    assert reactor.config.major_gibbs_solver_mode == "hamel_reduced"


def test_propagate_upstream_uses_cell_upflow_rates_under_holdup_transport():
    reactor = Reactor(ReactorConfig(n_cells=2, top_solid_inlet_frac=1.0))
    prev = reactor.cells[0]
    top = reactor.cells[1]
    prev.solid_state_model = "holdup_transport"
    prev.m_solid[0, :] = np.array([2.0, 4.0, 6.0, 8.0], dtype=float)
    prev.K_solid_auf.fill(0.1)
    prev.K_solid_ab.fill(0.2)

    reactor._propagate_upstream(1)

    np.testing.assert_allclose(top.m_solid_in[0, :], np.array([0.2, 0.0, 0.0, 0.8], dtype=float))


def test_top_cell_propagation_filters_vm_and_moisture_from_solid_inlet():
    reactor = Reactor(ReactorConfig(n_cells=2, top_solid_inlet_frac=0.4))
    prev = reactor.cells[0]
    top = reactor.cells[1]
    prev.m_solid[0, :] = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)

    reactor._propagate_upstream(1)

    np.testing.assert_allclose(top.m_solid_in[0, :], np.array([0.4, 0.0, 0.0, 1.6], dtype=float))


def test_thesis_mode_bed_only_builds_explicit_cyclone_and_return_leg_cells():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, n_freeboard_cells=0))

    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None
    assert all(cell.solid_state_model == "holdup_transport" for cell in reactor.cells + reactor.side_cells)
    assert [cell.cell_type for cell in reactor.side_cells] == ["cyclone", "return_leg"]
    assert len(reactor._solver_cells_for_nr()) == len(reactor.cells) + 2


def test_thesis_mode_with_freeboard_builds_full_solver_graph_when_no_local_refine():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, H_freeboard=2.0, n_freeboard_cells=2, freeboard_secondary_local_refine=1))

    assert reactor._use_explicit_freeboard_solver_graph() is True
    assert reactor._use_side_blocks_in_nr_boundary_path() is True
    assert all(cell.solid_state_model == "holdup_transport" for cell in reactor.cells)
    assert all(cell.solid_state_model == "freeboard_closure" for cell in reactor.freeboard_cells)
    assert len(reactor._solver_cells_for_nr()) == len(reactor.cells) + len(reactor.freeboard_cells) + len(reactor.side_cells)


def test_closure_owned_freeboard_keeps_side_blocks_in_boundary_path():
    reactor = Reactor(
        ReactorConfig(
            n_cells=2,
            thesis_mode=True,
            H_freeboard=2.0,
            n_freeboard_cells=2,
            explicit_freeboard_solver_graph_enabled=False,
        )
    )

    assert reactor._use_explicit_freeboard_cells() is True
    assert reactor._use_explicit_freeboard_solver_graph() is False
    assert reactor._use_side_blocks_in_nr_boundary_path() is True
    assert len(reactor._solver_cells_for_nr()) == len(reactor.cells) + len(reactor.side_cells)


def test_explicit_return_leg_recycle_keeps_gas_zero_and_solid_only():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, n_freeboard_cells=0, recirculation_frac=0.25))
    for cell in reactor.cells:
        cell.eps_b = 0.2
        cell.u_b = 0.5
        cell.eps_d_voidage = 0.5
        cell.R_solid.fill(0.0)
    top = reactor.cells[-1]
    top.N_d[:] = np.array([1.0] * 11, dtype=float)
    top.N_b[:] = np.array([0.5] * 11, dtype=float)
    top.m_solid[0, :] = np.array([2.0, 1.0, 0.5, 4.0], dtype=float)
    top.T = 1180.0

    reactor._initialize_explicit_side_block_states()
    reactor._apply_all_bc_for_nr()
    bot = reactor.cells[0]

    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None
    np.testing.assert_allclose(reactor.cyclone_cell.N_d_in, top.N_d + top.N_b)
    np.testing.assert_allclose(reactor.cyclone_cell.N_b_in, 0.0)
    np.testing.assert_allclose(reactor.return_leg_cell.N_d_in, 0.0)
    np.testing.assert_allclose(reactor.return_leg_cell.N_b_in, 0.0)
    np.testing.assert_allclose(reactor.cyclone_cell.m_solid_auf_in[0, S_CHAR], top._solid_upflow_rates()[0, S_CHAR])
    np.testing.assert_allclose(
        reactor.return_leg_cell.m_solid_ab_in[0, S_CHAR],
        reactor.cyclone_cell._solid_downflow_rates()[0, S_CHAR],
    )
    np.testing.assert_allclose(bot.N_rez_d, 0.0)
    np.testing.assert_allclose(bot.N_rez_b, 0.0)
    np.testing.assert_allclose(bot.m_solid_rez, reactor.return_leg_cell._solid_downflow_rates())
    assert bot.m_solid_rez[0, S_CHAR] > 0.0
    assert bot.m_solid_rez[0, S_ASH] > 0.0
    assert bot.m_solid_rez[0, S_VM] == 0.0
    assert bot.m_solid_rez[0, S_MOISTURE] == 0.0


def test_apply_all_bc_with_explicit_freeboard_routes_freeboard_to_side_blocks_and_backmixes_into_bed_top():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, H_freeboard=2.0, n_freeboard_cells=2, recirculation_frac=0.25))
    bed_top = reactor.cells[-1]
    fb0, fb1 = reactor.freeboard_cells
    bed_top.solid_state_model = "holdup_transport"
    bed_top.m_solid[0, S_CHAR] = 1.0
    bed_top.K_solid_auf[0, S_CHAR] = 0.4
    fb0.solid_state_model = "freeboard_closure"
    fb1.solid_state_model = "freeboard_closure"
    fb0.m_solid[0, S_CHAR] = 2.0
    fb0.K_solid_auf[0, S_CHAR] = 0.1
    fb0.K_solid_ab[0, S_CHAR] = 0.3
    fb1.m_solid[0, S_CHAR] = 3.0
    fb1.K_solid_auf[0, S_CHAR] = 0.2
    fb1.N_d[:] = 1.0
    fb1.T = 1150.0

    reactor._apply_all_bc_for_nr()

    np.testing.assert_allclose(bed_top.m_solid_ab_in[0, S_CHAR], fb0._solid_downflow_rates()[0, S_CHAR])
    assert bed_top.T_solid_ab_in == pytest.approx(fb0.T)
    assert fb0.m_solid_auf_in[0, S_CHAR] == 0.0
    assert fb1.m_solid_auf_in[0, S_CHAR] == 0.0
    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None
    np.testing.assert_allclose(reactor.cyclone_cell.m_solid_auf_in[0, S_CHAR], fb1._solid_upflow_rates()[0, S_CHAR])
    assert reactor.cyclone_cell.T_solid_auf_in == pytest.approx(fb1.T)
    np.testing.assert_allclose(
        reactor.return_leg_cell.m_solid_ab_in[0, S_CHAR],
        reactor.cyclone_cell._solid_downflow_rates()[0, S_CHAR],
    )
    assert reactor.return_leg_cell.T_solid_ab_in == pytest.approx(reactor.cyclone_cell.T)


def test_explicit_freeboard_inlet_filters_classes_not_active_in_closure_mask():
    reactor = Reactor(
        ReactorConfig(
            n_cells=2,
            thesis_mode=True,
            H_freeboard=2.0,
            n_freeboard_cells=1,
            n_age_classes=2,
        )
    )
    prev = reactor.cells[-1]
    fb = reactor.freeboard_cells[0]

    prev.solid_state_model = "holdup_transport"
    prev.m_solid.fill(0.0)
    prev.K_solid_auf.fill(0.0)
    prev.m_solid[0, S_CHAR] = 1.0
    prev.m_solid[1, S_CHAR] = 2.0
    prev.m_solid[0, S_ASH] = 0.5
    prev.m_solid[1, S_ASH] = 1.0
    prev.K_solid_auf[:, S_CHAR] = 0.5
    prev.K_solid_auf[:, S_ASH] = 0.5

    fb._freeboard_active_char_ash_mask = np.array([0.0, 1.0], dtype=np.float64)
    prev.N_d[:] = np.arange(1, len(prev.N_d) + 1, dtype=np.float64)
    prev.N_b[:] = 0.25 * np.arange(1, len(prev.N_b) + 1, dtype=np.float64)
    _set_explicit_freeboard_inlet_from_prev(prev, fb)
    update_freeboard_solid_transport_inflows([fb], bottom_below_cell=prev)

    np.testing.assert_allclose(fb.N_d_in, prev.N_d + prev.N_b)
    np.testing.assert_allclose(fb.N_b_in, 0.0)
    assert fb.m_solid_auf_in[0, S_CHAR] == 0.0
    assert fb.m_solid_auf_in[0, S_ASH] == 0.0
    assert fb.m_solid_auf_in[1, S_CHAR] == 0.0
    assert fb.m_solid_auf_in[1, S_ASH] == 0.0


def test_side_block_transport_coefficients_freeze_inside_inner_nr_boundary_updates():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, n_freeboard_cells=0))
    for cell in reactor.cells:
        cell.eps_b = 0.2
        cell.u_b = 0.5
        cell.eps_d_voidage = 0.5
        cell.R_solid.fill(0.0)
    top = reactor.cells[-1]
    top.N_d[:] = np.array([1.0] * len(top.N_d), dtype=float)
    top.N_b[:] = np.array([0.25] * len(top.N_b), dtype=float)
    top.m_solid[0, S_CHAR] = 3.0
    top.m_solid[0, S_ASH] = 1.0
    top.T = 1175.0

    reactor._apply_all_bc_for_nr()

    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None
    for cell in (reactor.cyclone_cell, reactor.return_leg_cell):
        cell._snapshot_vorabrechnung_hydrodynamics()
        cell.enable_inner_nr_vorabrechnung_freeze()

    reactor.cyclone_cell._vorab_K_solid_auf.fill(0.0)
    reactor.cyclone_cell._vorab_K_solid_ab.fill(0.0)
    reactor.cyclone_cell._vorab_K_solid_auf[0, S_CHAR] = 1.23
    reactor.cyclone_cell._vorab_K_solid_ab[0, S_CHAR] = 4.56
    reactor.return_leg_cell._vorab_K_solid_auf.fill(0.0)
    reactor.return_leg_cell._vorab_K_solid_ab.fill(0.0)
    reactor.return_leg_cell._vorab_K_solid_ab[0, S_CHAR] = 0.78

    top.N_d[:] = np.array([2.0] * len(top.N_d), dtype=float)
    top.T = 1230.0
    reactor._apply_all_bc_for_nr()

    assert reactor.cyclone_cell.K_solid_auf[0, S_CHAR] == 1.23
    assert reactor.cyclone_cell.K_solid_ab[0, S_CHAR] == 4.56
    assert reactor.return_leg_cell.K_solid_ab[0, S_CHAR] == 0.78
    np.testing.assert_allclose(
        reactor.return_leg_cell.m_solid_ab_in[0, S_CHAR],
        reactor.cyclone_cell._solid_downflow_rates()[0, S_CHAR],
    )


def test_bed_solid_transport_coefficients_follow_hamel_wake_and_top_ejection_formula():
    reactor = Reactor(ReactorConfig(n_cells=2, H_bed=2.0, D_bed=1.0))
    lower, upper = reactor.cells

    for cell in reactor.cells:
        cell.eps_b = 0.2
        cell.u_b = 0.5
        cell.eps_d_voidage = 0.5
        cell.m_solid.fill(0.0)
        cell.m_solid_zu.fill(0.0)
        cell.m_solid_rez.fill(0.0)
        cell.R_solid.fill(0.0)

    update_bed_solid_transport_coefficients(reactor.cells, reactor.config)

    expected_lower = 0.25 * 0.2 * 0.5 / ((1.0 - 0.2) * 1.0)
    expected_upper = expected_lower * 0.40
    np.testing.assert_allclose(lower.K_solid_auf[:, S_CHAR], expected_lower)
    np.testing.assert_allclose(lower.K_solid_auf[:, S_ASH], expected_lower)
    np.testing.assert_allclose(lower.K_solid_auf[:, S_VM], 0.0)
    np.testing.assert_allclose(lower.K_solid_auf[:, S_MOISTURE], 0.0)
    np.testing.assert_allclose(upper.K_solid_auf[:, S_CHAR], expected_upper)
    np.testing.assert_allclose(upper.K_solid_auf[:, S_ASH], expected_upper)
    np.testing.assert_allclose(upper.K_solid_auf[:, S_VM], 0.0)
    np.testing.assert_allclose(upper.K_solid_auf[:, S_MOISTURE], 0.0)


def test_seed_holdup_from_inflows_uses_positive_local_char_source():
    reactor = Reactor(ReactorConfig(n_cells=1, thesis_mode=True))
    cell = reactor.cells[0]
    cell.m_solid.fill(0.0)
    cell.m_solid_zu.fill(0.0)
    cell.m_solid_rez.fill(0.0)
    cell.m_solid_in.fill(0.0)
    cell.m_solid_auf_in.fill(0.0)
    cell.m_solid_ab_in.fill(0.0)
    cell.R_solid.fill(0.0)
    cell.K_solid_auf.fill(0.0)
    cell.K_solid_ab.fill(0.0)
    cell.K_solid_auf[0, S_CHAR] = 2.0
    cell.R_solid[0, S_CHAR] = 0.4
    cell.R_solid[1, S_CHAR] = -0.8

    changed = seed_holdup_from_inflows(cell)

    assert changed is True
    assert cell.m_solid[0, S_CHAR] == 0.2
    assert cell.m_solid[1, S_CHAR] == 0.0


def test_seed_holdup_from_inflows_raises_underfilled_inventory_without_lowering():
    reactor = Reactor(ReactorConfig(n_cells=1, thesis_mode=True))
    cell = reactor.cells[0]
    cell.m_solid.fill(0.0)
    cell.m_solid[0, S_CHAR] = 0.1
    cell.m_solid[0, S_ASH] = 0.5
    cell.m_solid_zu.fill(0.0)
    cell.m_solid_zu[0, S_CHAR] = 0.6
    cell.m_solid_zu[0, S_ASH] = 0.2
    cell.K_solid_auf.fill(0.0)
    cell.K_solid_ab.fill(0.0)
    cell.K_solid_auf[0, S_CHAR] = 2.0
    cell.K_solid_auf[0, S_ASH] = 2.0

    changed = seed_holdup_from_inflows(cell)

    assert changed is True
    assert cell.m_solid[0, S_CHAR] == pytest.approx(0.3)
    assert cell.m_solid[0, S_ASH] == pytest.approx(0.5)


def test_bed_solid_transport_backmix_coefficients_follow_top_down_continuity():
    reactor = Reactor(ReactorConfig(n_cells=2, H_bed=2.0, D_bed=1.0))
    lower, upper = reactor.cells
    for cell in reactor.cells:
        cell.eps_b = 0.2
        cell.u_b = 0.5
        cell.eps_d_voidage = 0.5
        cell.m_solid.fill(0.0)
        cell.m_solid_zu.fill(0.0)
        cell.m_solid_rez.fill(0.0)
        cell.R_solid.fill(0.0)
    lower.m_solid_zu[0, S_CHAR] = 10.0

    update_bed_solid_transport_coefficients(reactor.cells, reactor.config)

    holdup_lower = cell_total_solid_holdup(lower)
    holdup_upper = cell_total_solid_holdup(upper)
    m_auf_lower = float(lower.K_solid_auf[0, S_CHAR]) * holdup_lower
    m_auf_upper = float(upper.K_solid_auf[0, S_CHAR]) * holdup_upper

    expected_upper_kab = max(m_auf_lower - m_auf_upper, 0.0) / holdup_upper
    expected_lower_kab_char = max(-m_auf_lower + max(m_auf_lower - m_auf_upper, 0.0) + 10.0, 0.0) / holdup_lower
    expected_lower_kab_ash = max(-m_auf_lower + max(m_auf_lower - m_auf_upper, 0.0), 0.0) / holdup_lower

    np.testing.assert_allclose(upper.K_solid_ab[:, S_CHAR], expected_upper_kab)
    np.testing.assert_allclose(upper.K_solid_ab[:, S_ASH], expected_upper_kab)
    np.testing.assert_allclose(upper.K_solid_ab[:, S_VM], 0.0)
    np.testing.assert_allclose(upper.K_solid_ab[:, S_MOISTURE], 0.0)
    np.testing.assert_allclose(lower.K_solid_ab[:, S_CHAR], expected_lower_kab_char)
    np.testing.assert_allclose(lower.K_solid_ab[:, S_ASH], expected_lower_kab_ash)
    np.testing.assert_allclose(lower.K_solid_ab[:, S_VM], 0.0)
    np.testing.assert_allclose(lower.K_solid_ab[:, S_MOISTURE], 0.0)


def test_bed_solid_transport_inflows_use_neighbor_holdup_times_frozen_k():
    reactor = Reactor(ReactorConfig(n_cells=3))
    lower, mid, upper = reactor.cells
    lower.m_solid[0, :] = np.array([2.0, 0.5, 0.1, 1.0], dtype=float)
    upper.m_solid[0, :] = np.array([3.0, 0.7, 0.2, 1.5], dtype=float)
    lower.K_solid_auf.fill(0.2)
    upper.K_solid_ab.fill(0.3)

    update_bed_solid_transport_inflows(reactor.cells)

    expected_auf = 0.2 * lower.m_solid
    expected_ab = 0.3 * upper.m_solid
    expected_auf[:, S_VM] = 0.0
    expected_auf[:, S_MOISTURE] = 0.0
    expected_ab[:, S_VM] = 0.0
    expected_ab[:, S_MOISTURE] = 0.0
    np.testing.assert_allclose(mid.m_solid_auf_in, expected_auf)
    np.testing.assert_allclose(mid.m_solid_ab_in, expected_ab)


def test_apply_all_bc_replays_frozen_bed_transport_coefficients_inside_inner_nr():
    reactor = Reactor(ReactorConfig(n_cells=2))
    lower, upper = reactor.cells
    lower.m_solid[0, :] = np.array([2.0, 0.5, 0.1, 1.0], dtype=float)
    upper.m_solid[0, :] = np.array([3.0, 0.7, 0.2, 1.5], dtype=float)
    for cell, k_auf, k_ab in ((lower, 0.2, 0.05), (upper, 0.3, 0.15)):
        cell.u_mf = 1.0
        cell.u_b = 2.0
        cell.d_b = 0.1
        cell.eps_b = 0.2
        cell.eps_d = 0.8
        cell.eps_d_voidage = 0.5
        cell.K_bd = 3.0
        cell.V_b = 0.1
        cell.V_d = 0.2
        cell.u0 = 1.1
        cell.u_d = 0.4
        cell.n_rz = 1.0
        cell.K_solid_auf.fill(k_auf)
        cell.K_solid_ab.fill(k_ab)
        cell._snapshot_vorabrechnung_hydrodynamics()
        cell.enable_inner_nr_vorabrechnung_freeze()
        cell.K_solid_auf.fill(-1.0)
        cell.K_solid_ab.fill(-1.0)

    reactor._apply_all_bc_for_nr()

    np.testing.assert_allclose(lower.K_solid_auf, 0.2)
    np.testing.assert_allclose(lower.K_solid_ab, 0.05)
    np.testing.assert_allclose(upper.K_solid_auf, 0.3)
    np.testing.assert_allclose(upper.K_solid_ab, 0.15)
    expected_auf = 0.2 * lower.m_solid
    expected_auf[:, S_VM] = 0.0
    expected_auf[:, S_MOISTURE] = 0.0
    np.testing.assert_allclose(upper.m_solid_auf_in, expected_auf)


def test_reactor_carbon_conversion_uses_cell_upflow_rates_under_holdup_transport():
    reactor = Reactor(ReactorConfig(n_cells=2))
    bot, top = reactor.cells
    bot.m_solid_zu[0, S_CHAR] = 10.0
    top.solid_state_model = "holdup_transport"
    top.m_solid[0, S_CHAR] = 10.0
    top.K_solid_auf.fill(0.1)
    top.K_solid_ab.fill(0.2)

    assert reactor._compute_carbon_conversion() == 0.9


def test_reactor_carbon_conversion_excludes_bottom_recycle_from_system_exit():
    reactor = Reactor(ReactorConfig(n_cells=2, recirculation_frac=0.25))
    bot, top = reactor.cells
    bot.m_solid_zu[0, S_CHAR] = 10.0
    top.solid_state_model = "holdup_transport"
    top.m_solid[0, S_CHAR] = 10.0
    top.K_solid_auf.fill(0.1)
    top.K_solid_ab.fill(0.0)
    bot.m_solid_rez[0, S_CHAR] = 0.25

    assert reactor._compute_carbon_conversion() == 0.925


def test_char_mass_audit_reuses_solved_reactor_state_without_solving():
    summarize_reactor_char_ledger = _load_char_audit_summary_fn()
    reactor = Reactor(ReactorConfig(n_cells=2))
    bot, top = reactor.cells
    bot.m_solid_zu[0, S_CHAR] = 10.0
    bot.R_solid[0, S_CHAR] = -4.0
    top.solid_state_model = "holdup_transport"
    top.m_solid[0, S_CHAR] = 8.0
    top.K_solid_auf[0, S_CHAR] = 0.5
    bot.m_solid_rez[0, S_CHAR] = 1.0

    payload = summarize_reactor_char_ledger(
        reactor,
        {"converged": True, "rms_scaled_final": 0.0, "carbon_conv": 0.7},
        mode="unit",
        solve_kwargs={},
    )

    summary = payload["system_char_summary"]
    assert summary["fresh_char_kg_s"] == 10.0
    assert summary["reaction_char_kg_s"] == -4.0
    assert summary["top_up_char_kg_s"] == 4.0
    assert summary["bottom_recycle_char_kg_s"] == 1.0
    assert summary["system_exit_proxy_top_up_minus_recycle_kg_s"] == 3.0
    assert summary["fresh_minus_reaction_minus_exit_proxy_kg_s"] == 3.0


def test_char_mass_audit_reports_bed_transport_profile_by_size_class():
    summarize_reactor_char_ledger = _load_char_audit_summary_fn()
    reactor = Reactor(ReactorConfig(n_cells=2, n_age_classes=2))
    bot, top = reactor.cells
    bot.solid_state_model = "holdup_transport"
    top.solid_state_model = "holdup_transport"
    bot.m_solid[:, S_CHAR] = [1.0, 2.0]
    top.m_solid[:, S_CHAR] = [3.0, 4.0]
    bot.K_solid_auf[:, S_CHAR] = [0.1, 0.2]
    top.K_solid_auf[:, S_CHAR] = [0.3, 0.4]
    top.m_solid_auf_in[:, S_CHAR] = bot._solid_upflow_rates()[:, S_CHAR]

    payload = summarize_reactor_char_ledger(
        reactor,
        {"converged": True, "rms_scaled_final": 0.0, "carbon_conv": 0.0},
        mode="unit",
        solve_kwargs={},
    )

    profile = payload["bed_transport_profile"]
    assert profile[0]["char_holdup_by_class_kg"] == [1.0, 2.0]
    assert profile[0]["char_up_out_by_class_kg_s"] == pytest.approx([0.1, 0.4])
    assert profile[1]["char_up_out_total_kg_s"] == pytest.approx(2.5)
    assert profile[1]["char_auf_in_gap_vs_below_kg_s"] == pytest.approx(0.0)
    assert payload["bed_transport_summary"]["max_abs_auf_in_gap_vs_below_kg_s"] == pytest.approx(0.0)


def test_char_mass_audit_reports_phase2_freeboard_closure_transport_summary():
    summarize_reactor_char_ledger = _load_char_audit_summary_fn()
    reactor = Reactor(ReactorConfig(n_cells=1, thesis_mode=True, H_freeboard=2.0, n_freeboard_cells=2, n_age_classes=1))
    bed = reactor.cells[0]
    fb0, fb1 = reactor.freeboard_cells
    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None
    bed.solid_state_model = "holdup_transport"
    bed.m_solid[0, S_CHAR] = 10.0
    bed.K_solid_auf[0, S_CHAR] = 0.3
    fb0.solid_state_model = "freeboard_closure"
    fb1.solid_state_model = "freeboard_closure"
    fb0.m_solid[0, S_CHAR] = 1.0
    fb1.m_solid[0, S_CHAR] = 0.5
    fb0.K_solid_ab[0, S_CHAR] = 0.2
    fb1.K_solid_auf[0, S_CHAR] = 0.4
    reactor.cyclone_cell.m_solid_auf_in[0, S_CHAR] = 0.2
    reactor.return_leg_cell.m_solid_ab_in[0, S_CHAR] = 0.18

    payload = summarize_reactor_char_ledger(
        reactor,
        {
            "converged": True,
            "rms_scaled_final": 0.0,
            "carbon_conv": 0.0,
            "freeboard_active": True,
            "freeboard_bed_top_up_char_kg_s": 3.0,
            "freeboard_entrained_char_kg_s": [0.8, 0.2],
            "freeboard_entrained_return_char_profile_kg_s": [0.2, 0.05],
            "freeboard_solid_holdup_char_profile_kg": [1.0, 0.5],
            "freeboard_solid_holdup_char_before_profile_kg": [1.2, 0.6],
            "freeboard_explicit_char_sink_applied_profile_kg": [0.1, 0.0],
            "freeboard_entrained_exit_char_kg_s": 0.2,
            "freeboard_entrained_return_char_kg_s": 0.2,
            "freeboard_cyclone_capture_char_kg_s": 0.18,
            "freeboard_axial_z_m": [7.0, 8.0],
            "freeboard_axial_xi": [0.7, 0.8],
            "freeboard_carry_ratio_profile": [1.1, 0.9],
        },
        mode="phase2-unit",
        solve_kwargs={},
    )

    summary = payload["freeboard_closure_transport_summary"]
    assert summary["freeboard_active"] is True
    assert summary["bed_top_up_char_kg_s"] == pytest.approx(3.0)
    assert summary["current_bed_top_up_char_kg_s"] == pytest.approx(3.0)
    assert summary["gap_current_bed_top_vs_closure_bed_top_kg_s"] == pytest.approx(0.0)
    assert summary["freeboard_return_char_kg_s"] == pytest.approx(0.2)
    assert summary["freeboard_exit_char_kg_s"] == pytest.approx(0.2)
    assert summary["system_escape_after_cyclone_char_kg_s"] == pytest.approx(0.02)
    assert summary["gap_result_exit_vs_explicit_top_kg_s"] == pytest.approx(0.0)
    assert summary["gap_cyclone_capture_vs_return_leg_in_kg_s"] == pytest.approx(0.0)
    assert payload["freeboard_closure_transport_profile"][0]["char_return_kg_s"] == pytest.approx(0.2)
    interface = payload["bed_freeboard_interface_component_balance"]
    char_row = next(row for row in interface if row["component"] == "char")
    total_row = next(row for row in interface if row["component"] == "char+ash")
    assert char_row["ab_in_kg_s"] == pytest.approx(0.0)
    assert char_row["freeboard_bottom_return_kg_s"] == pytest.approx(0.2)
    assert char_row["bed_top_ab_in_gap_vs_freeboard_return_kg_s"] == pytest.approx(-0.2)
    assert total_row["component"] == "char+ash"
    sensitivity = payload["top_bed_solid_local_sensitivity"]
    char_sens = next(row for row in sensitivity if row["component"] == "char")
    assert char_sens["d_residual_d_holdup_1_s"] == pytest.approx(-0.3)
    assert char_sens["minus_K_sum_1_s"] == pytest.approx(-0.3)


def test_sync_freeboard_cells_from_closure_loads_holdup_transport_coefficients():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, H_freeboard=2.0, n_freeboard_cells=2, n_age_classes=1))
    bed_top = reactor.cells[-1]
    bed_top.u_b = 1.2
    bed_top.d_b = 0.12
    bed_top.solid_state_model = "holdup_transport"
    bed_top.m_solid[0, S_CHAR] = 0.5
    bed_top.K_solid_auf[0, S_CHAR] = 0.2

    n_in = np.zeros(len(GAS_SPECIES), dtype=np.float64)
    n_in[GAS_SPECIES_INDEX["CO"]] = 1.0
    n_in[GAS_SPECIES_INDEX["H2"]] = 1.0
    n_in[GAS_SPECIES_INDEX["H2O"]] = 1.0
    n_in[GAS_SPECIES_INDEX["N2"]] = 3.0
    fb = simulate_freeboard(
        N_in=n_in,
        T_in=1100.0,
        P=2.5e6,
        D_bed=0.6,
        H_freeboard=2.0,
        n_cells=2,
        u_b_bed_top=1.2,
        d_b_bed_top=0.12,
        eps_b_bed_top=0.3,
        eps_d_void_bed_top=0.5,
        rho_solid_bed_top=1700.0,
        d_p_classes_bed_top=np.array([7e-4]),
        m_char_classes_bed_top=np.array([0.02]),
        m_ash_classes_bed_top=np.array([0.01]),
        fuel_type="coal",
        trajectory_model="analytical_wirsum",
    )

    reactor._sync_freeboard_cells_from_closure(fb)

    assert reactor.freeboard_cells[0].solid_state_model == "freeboard_closure"
    assert reactor.freeboard_cells[0].m_solid[0, S_CHAR] >= 0.0
    assert reactor.freeboard_cells[0].K_solid_auf[0, S_CHAR] >= 0.0
    assert reactor.freeboard_cells[0].K_solid_ab[0, S_ASH] >= 0.0
    np.testing.assert_allclose(reactor.freeboard_cells[1].m_solid_auf_in, 0.0)


def test_sync_freeboard_cells_from_closure_can_preserve_nr_gas_state():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, H_freeboard=2.0, n_freeboard_cells=2, n_age_classes=1))
    bed_top = reactor.cells[-1]
    bed_top.u_b = 1.2
    bed_top.d_b = 0.12

    n_in = np.zeros(len(GAS_SPECIES), dtype=np.float64)
    n_in[GAS_SPECIES_INDEX["CO"]] = 1.0
    n_in[GAS_SPECIES_INDEX["H2"]] = 1.0
    n_in[GAS_SPECIES_INDEX["H2O"]] = 1.0
    n_in[GAS_SPECIES_INDEX["N2"]] = 3.0
    fb = simulate_freeboard(
        N_in=n_in,
        T_in=1100.0,
        P=2.5e6,
        D_bed=0.6,
        H_freeboard=2.0,
        n_cells=2,
        u_b_bed_top=1.2,
        d_b_bed_top=0.12,
        eps_b_bed_top=0.3,
        eps_d_void_bed_top=0.5,
        rho_solid_bed_top=1700.0,
        d_p_classes_bed_top=np.array([7e-4]),
        m_char_classes_bed_top=np.array([0.02]),
        m_ash_classes_bed_top=np.array([0.01]),
        fuel_type="coal",
        trajectory_model="analytical_wirsum",
    )

    fb0 = reactor.freeboard_cells[0]
    fb0.T = 987.0
    fb0.N_b[:] = np.linspace(0.01, 0.02, len(GAS_SPECIES))
    fb0.N_d[:] = np.linspace(0.1, 0.2, len(GAS_SPECIES))
    t_before = float(fb0.T)
    n_b_before = np.array(fb0.N_b, copy=True)
    n_d_before = np.array(fb0.N_d, copy=True)

    reactor._sync_freeboard_cells_from_closure(fb, preserve_gas_state=True)

    assert reactor.freeboard_cells[0].solid_state_model == "freeboard_closure"
    assert reactor.freeboard_cells[0].m_solid[0, S_CHAR] >= 0.0
    assert reactor.freeboard_cells[0].K_solid_auf[0, S_CHAR] >= 0.0
    assert reactor.freeboard_cells[0].T == pytest.approx(t_before)
    np.testing.assert_allclose(reactor.freeboard_cells[0].N_b, n_b_before)
    np.testing.assert_allclose(reactor.freeboard_cells[0].N_d, n_d_before)


def test_sync_freeboard_cells_from_closure_reconstructs_componentwise_k():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, H_freeboard=1.0, n_freeboard_cells=1, n_age_classes=1))
    state = SimpleNamespace(
        T=1050.0,
        N=np.array([0.0] * len(GAS_SPECIES), dtype=np.float64),
        m_hold_char_classes=np.array([2.0], dtype=np.float64),
        m_hold_ash_classes=np.array([1.0], dtype=np.float64),
        m_dot_auf_char_classes=np.array([1.0], dtype=np.float64),
        m_dot_auf_ash_classes=np.array([0.25], dtype=np.float64),
        m_dot_ab_char_classes=np.array([0.2], dtype=np.float64),
        m_dot_ab_ash_classes=np.array([0.05], dtype=np.float64),
        m_char_sink_applied_classes=np.array([0.0], dtype=np.float64),
        segment_index=0,
    )

    reactor._sync_freeboard_cells_from_closure({"states": [state]})

    cell = reactor.freeboard_cells[0]
    assert cell.K_solid_auf[0, S_CHAR] == pytest.approx(0.5)
    assert cell.K_solid_auf[0, S_ASH] == pytest.approx(0.25)
    assert cell.K_solid_ab[0, S_CHAR] == pytest.approx(0.1)
    assert cell.K_solid_ab[0, S_ASH] == pytest.approx(0.05)


def test_sync_freeboard_cells_from_refined_closure_aggregates_back_to_nominal_cells():
    reactor = Reactor(
        ReactorConfig(
            n_cells=2,
            thesis_mode=True,
            H_freeboard=2.0,
            n_freeboard_cells=2,
            n_age_classes=1,
            freeboard_secondary_injection_xi=0.6,
            freeboard_secondary_O2_mol_s=0.05,
            freeboard_secondary_N2_mol_s=0.05 * 3.76,
            freeboard_secondary_local_refine=4,
        )
    )
    fb = simulate_freeboard(
        N_in=np.array([1.0 if sp == "N2" else 0.0 for sp in GAS_SPECIES], dtype=np.float64),
        T_in=1100.0,
        P=2.5e6,
        D_bed=0.6,
        H_freeboard=2.0,
        n_cells=2,
        u_b_bed_top=1.2,
        d_b_bed_top=0.12,
        eps_b_bed_top=0.3,
        eps_d_void_bed_top=0.5,
        rho_solid_bed_top=1700.0,
        d_p_classes_bed_top=np.array([7e-4]),
        m_char_classes_bed_top=np.array([0.02]),
        m_ash_classes_bed_top=np.array([0.01]),
        fuel_type="coal",
        trajectory_model="analytical_wirsum",
        secondary_injection_xi=0.6,
        secondary_O2_mol_s=0.05,
        secondary_N2_mol_s=0.05 * 3.76,
        secondary_local_refine=4,
    )

    reactor._sync_freeboard_cells_from_closure(fb)

    assert reactor._use_explicit_freeboard_solver_graph() is True
    assert len(fb["states"]) > len(reactor.freeboard_cells)
    assert reactor.freeboard_cells[0].solid_state_model == "freeboard_closure"
    assert reactor.freeboard_cells[0].m_solid[0, S_CHAR] >= 0.0
    assert reactor.freeboard_cells[1].K_solid_auf[0, S_CHAR] >= 0.0


def test_bed_dh_profile_builds_nonuniform_bed_mesh():
    cfg = ReactorConfig(
        n_cells=3,
        H_bed=1.0,
        bed_dh_profile=(0.2, 0.3, 0.5),
    )
    reactor = Reactor(cfg)
    dh = [c.geo.dh for c in reactor.cells]
    hc = [c.geo.h_center for c in reactor.cells]
    assert dh == pytest.approx([0.2, 0.3, 0.5])
    assert hc == pytest.approx([0.1, 0.35, 0.75])


def test_nonuniform_bed_vorabrechnung_tau_uses_total_bed_height():
    cfg = ReactorConfig(
        n_cells=3,
        H_bed=1.0,
        bed_dh_profile=(0.2, 0.3, 0.5),
    )
    reactor = Reactor(cfg)
    for cell in reactor.cells:
        cell.u_mf = 0.25
        assert vorabrechnung_tau_for_cell(cell) == pytest.approx(4.0)


def test_bed_dh_profile_rejects_inconsistent_height():
    with pytest.raises(ValueError, match="sum\\(bed_dh_profile\\)"):
        Reactor(ReactorConfig(n_cells=2, H_bed=1.0, bed_dh_profile=(0.2, 0.7)))
