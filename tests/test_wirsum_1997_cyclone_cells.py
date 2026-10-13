"""Hamel Bild 7.7：旋风 cell 衡算，回料进残差。不用 Lapple。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.solvers.wirsum_1997_cyclone_cells import (
    CycloneCellK,
    CycloneLoopSpec,
    bed_residual_with_rez,
    build_stage_f_cyclone_payload,
    cyclone_k_implied,
    d_bed_d_Mret,
    loop_residual,
    n_unknowns,
    recycle_streams,
    solve_loop,
)


def test_module_does_not_use_lapple() -> None:
    import src.solvers.wirsum_1997_cyclone_cells as mod

    src = inspect.getsource(mod)
    assert "from src.physics.wirsum_1997_cyclone" not in src
    assert "grade_efficiency_lapple" not in src
    assert "from src.core.freeboard" not in src


def test_gas_passthrough_eq24_degenerate() -> None:
    spec = CycloneLoopSpec(n_cyc=2, recycle=False, N_fb=1.8, m_leave=0.5)
    out = solve_loop(spec)
    assert out["physical_convergence"] is True
    np.testing.assert_allclose(out["N_mol_s"], [1.8, 1.8])


def test_solids_leave_equals_fly_plus_under() -> None:
    spec = CycloneLoopSpec(n_cyc=2, recycle=False, m_leave=0.5)
    out = solve_loop(spec)
    assert out["leave_minus_fly_under"] == pytest.approx(0.0, abs=1e-12)
    assert out["m_rez_kg_s"] == pytest.approx(0.0)
    assert out["m_aus_kg_s"] == pytest.approx(out["m_underflow_kg_s"])


def test_recycle_enters_bed_residual() -> None:
    off = CycloneLoopSpec(n_cyc=2, recycle=False, m_leave=0.5, m_feed=0.32)
    on = CycloneLoopSpec(n_cyc=2, recycle=True, m_leave=0.5, m_feed=0.32)
    a = solve_loop(off)
    b = solve_loop(on)
    assert a["dF_bed_dMret"] == pytest.approx(0.0)
    assert b["dF_bed_dMret"] == pytest.approx(on.k.K_return)
    assert a["F_bed"] == pytest.approx(bed_residual_with_rez(off, 0.0))
    assert b["F_bed"] == pytest.approx(bed_residual_with_rez(on, b["m_rez_kg_s"]))
    assert b["m_rez_kg_s"] == pytest.approx(b["m_underflow_kg_s"])
    assert d_bed_d_Mret(on) != pytest.approx(d_bed_d_Mret(off))


def test_n_cyc_one_or_two() -> None:
    assert n_unknowns(1) == 3
    assert n_unknowns(2) == 5
    one = solve_loop(CycloneLoopSpec(n_cyc=1, recycle=False, m_leave=0.5))
    two = solve_loop(CycloneLoopSpec(n_cyc=2, recycle=False, m_leave=0.5))
    assert one["physical_convergence"] is True
    assert two["physical_convergence"] is True
    assert one["N_mol_s"].shape == (1,)
    assert two["N_mol_s"].shape == (2,)


def test_residual_zero_at_root() -> None:
    spec = CycloneLoopSpec(n_cyc=2, recycle=True, m_leave=0.5)
    out = solve_loop(spec)
    F = loop_residual(
        np.concatenate([out["N_mol_s"], out["M_cyc_kg"], np.array([out["M_ret_kg"]])]),
        spec,
    )
    assert float(np.max(np.abs(F))) < 1.0e-10
    rez, aus = recycle_streams(out["M_ret_kg"], spec)
    assert aus == pytest.approx(0.0)
    assert rez == pytest.approx(out["m_rez_kg_s"])


def test_scaled_k_keeps_split_inverts_inventory() -> None:
    k1 = CycloneCellK(K_auf=0.25, K_ab=0.75, K_return=1.0)
    k2 = CycloneCellK(K_auf=0.50, K_ab=1.50, K_return=2.0)
    a = solve_loop(CycloneLoopSpec(n_cyc=2, recycle=False, m_leave=0.5, k=k1))
    b = solve_loop(CycloneLoopSpec(n_cyc=2, recycle=False, m_leave=0.5, k=k2))
    assert a["implied_split_down"] == pytest.approx(b["implied_split_down"])
    assert a["implied_split_down"] == pytest.approx(0.75)
    assert a["tau_cyc_s"] == pytest.approx(2.0 * b["tau_cyc_s"])
    np.testing.assert_allclose(a["M_cyc_kg"], 2.0 * b["M_cyc_kg"], atol=1e-12)
    assert a["M_ret_kg"] == pytest.approx(2.0 * b["M_ret_kg"])
    d = cyclone_k_implied(k1)
    assert d["not_a_separation_law"] is True
    assert d["K_source"] == "frozen_context"


def test_payload_is_not_experiment_and_not_lapple_mainline() -> None:
    payload = build_stage_f_cyclone_payload()
    assert payload["solved"] is False
    assert payload["K_model"] == "frozen_context"
    assert "engineering_lapple_1951" in payload["not"]
    dumped = str(payload)
    assert "y_CO_target" not in dumped
    assert payload["recycle_closed"]["dF_bed_dMret"] == pytest.approx(0.0)
    assert payload["recycle_open"]["dF_bed_dMret"] > 0.0
    assert payload["recycle_closed"]["physical_convergence"] is True
