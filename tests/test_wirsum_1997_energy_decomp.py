"""Eq. 2.7：bed0 −10 MW 来自 H_zu 进、H_sout=0（holdup 为零）。"""

from __future__ import annotations

import inspect

from src.solvers.wirsum_1997_energy_decomp import (
    build_eq27_energy_decomp_payload,
    decompose_cell_eq27,
)


def test_eq27_decomp_matches_live_residual() -> None:
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    reactor.cells[0].residuals()
    row = decompose_cell_eq27(reactor.cells[0], 0)
    assert row["decomp_matches_live"] is True
    assert row["H_in_parts_match_api"] is True


def test_bed0_minus_10mw_is_feed_enthalpy_with_zero_holdup() -> None:
    payload = build_eq27_energy_decomp_payload()
    before = payload["before"]["per_cell"][0]
    after = payload["after"]["per_cell"][0]
    assert before["dominant_H_in_term"] == "H_solid_zu"
    assert before["holdup_zero"] is True
    assert before["H_solid_out_zero"] is True
    assert before["vm_in_not_in_energy_out"] is True
    assert abs(float(before["residual_W"] or 0.0)) > 1.0e6
    assert float(before["H_zu_by_comp_W"]["vm"] or 0.0) < -1.0e6
    assert float(before["H_zu_by_comp_W"]["char"] or 0.0) < -1.0e6
    assert after["holdup_zero"] is True
    assert after["H_solid_out_zero"] is True
    assert abs(float(after["residual_W"] or 0.0)) > 1.0e6
    assert after["K_auf_max_1_per_s"]["char"] > 0.0
    assert after["K_auf_max_1_per_s"]["vm"] == 0.0


def test_energy_decomp_does_not_add_vm_outlet_or_t_fence() -> None:
    import src.solvers.wirsum_1997_energy_decomp as mod

    source = inspect.getsource(mod)
    assert "temperature_fence" in source
    assert "vm_in_eq27_outlet" in source
    dumped = str(build_eq27_energy_decomp_payload())
    assert "experiment" not in dumped.lower()
    assert "y_CO_target" not in dumped
