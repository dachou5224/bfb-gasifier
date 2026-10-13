"""审阅 WP-A / WP-B 的固定状态门。不重解长链。"""

import numpy as np
import pytest
from dataclasses import replace

from src.solvers.wesseling_case_contract import bed_phase_volumes_m3, wesseling1_conditional_boundary


def test_bed_phase_volumes_sum_to_the_geometric_bed() -> None:
    bound = wesseling1_conditional_boundary()
    volumes = bed_phase_volumes_m3(
        bound["area_m2"], bound["H_bed_m"], bound["eps_b"], bound["eps_d"],
    )
    assert volumes["V_b_m3"] + volumes["V_d_m3"] + volumes["V_solids_m3"] == pytest.approx(
        volumes["V_bed_m3"]
    )
    assert bound["V_d_m3"] + bound["V_b_m3"] < bound["V_bed_m3"]
    assert bound["voidage_status"] == "kernel_default_not_tabulated"
    assert (bound["V_d_m3"] + bound["V_b_m3"]) / bound["V_bed_m3"] < 1.0


def test_illegal_state_does_not_pass_with_a_zero_residual() -> None:
    from src.solvers.wirsum_1997_named_f_control import certify_controlled_state

    x = np.ones(50, dtype=np.float64)
    x[-2:] = 1100.0
    f = np.zeros(50, dtype=np.float64)
    x[-1] = 100.0
    cert = certify_controlled_state(x, f, hydro=None)
    assert cert["passed"] is False
    assert cert["temperature_domain"] == "temporary_solver_domain_250_2500_K"
    x[-1] = 1100.0
    x[9] = -1.0e-6
    assert certify_controlled_state(x, f, hydro=None)["passed"] is False


def test_joint_step_opens_bubble_h2_and_ch4_when_the_full_merit_drops(monkeypatch) -> None:
    from src.solvers.wirsum_1997_named_f_control import joint_open_zeros

    x = np.zeros(50, dtype=np.float64)
    x[12] = 0.4
    x[13] = 1.5
    x[-2:] = 1100.0
    f = np.zeros(50, dtype=np.float64)
    f[9] = 3.0e5
    f[11] = 0.179
    calls = {"n": 0}

    def fake_residual(z, q, stats=None):
        calls["n"] += 1
        row = np.zeros(50, dtype=np.float64)
        if float(z[9]) > 0.0 and float(z[11]) > 0.0:
            row[9] = 1.0
            row[11] = 1.0
        else:
            row[9] = 3.0e5
            row[11] = 0.179
        return row

    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control._count_residual", fake_residual,
    )
    opened_x, _f, opened = joint_open_zeros(x, f, [9, 11], q=None)
    assert calls["n"] > 0
    assert opened == [9, 11]
    assert float(opened_x[9]) > 0.0 and float(opened_x[11]) > 0.0


def test_rejected_reactivation_keeps_q_with_the_accepted_state(monkeypatch) -> None:
    from src.solvers.wirsum_1997_named_f import outer_P_named
    from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
    from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x0 = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"].copy()
    x0[-2] += 50.0
    seen = []

    def fake_open(q, x, stats=None):
        y = np.array(x, copy=True)
        y[-2] = float(x[-2]) + 80.0
        return y, [7]

    def fake_correct(q, x, **kwargs):
        seen.append((float(q.energy.species.geo_ex.T), float(x[-2])))
        return np.array(x, copy=True)

    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control.reactivate_zero_unknowns", fake_open,
    )
    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control._e8_log_scaled_correct", fake_correct,
    )
    out = solve_named_controlled(nq, x0, deadline_s=2.0)
    assert seen
    assert seen[0][0] == pytest.approx(seen[0][1])
    assert out["x"][-2] == pytest.approx(float(x0[-2]))
    _ = outer_P_named


def test_progress_allows_a_block_tradeoff_and_an_energy_only_step() -> None:
    from src.solvers.wirsum_1997_named_f_control import step_improves_progress

    old = np.zeros(50, dtype=np.float64)
    old[0] = 7.74468e-6
    old[14] = 1.40873e-6
    old[40] = 1.90527e-5
    old[-2] = -0.00109757
    old[-1] = 0.00013564
    traded = old.copy()
    traded[0] = 1.2598e-5
    traded[14] = 5.78e-7
    traded[40] = 8.073e-6
    traded[-2] = -0.0352
    traded[-1] = 0.0903
    assert step_improves_progress(old, traded)
    guarded = old.copy()
    guarded[14] = 1.0e-7
    lost = guarded.copy()
    lost[14] = 2.0e-6
    lost[40] = 1.0e-6
    assert not step_improves_progress(guarded, lost)
    energy_old = np.zeros(50, dtype=np.float64)
    energy_old[0] = 1.0e-7
    energy_old[14] = 1.0e-7
    energy_old[40] = 1.0e-7
    energy_old[-2] = -2.0
    energy_new = energy_old.copy()
    energy_new[-2] = -1.5
    assert step_improves_progress(energy_old, energy_new)


def test_backtrack_keeps_the_fractional_step_with_lower_composite_merit(monkeypatch) -> None:
    from src.solvers.wirsum_1997_named_f_control import composite_backtrack

    accepted = np.ones(50, dtype=np.float64)
    accepted[-2:] = 1100.0
    proposed = accepted.copy()
    proposed[0] = 2.0

    def fake_residual(z, q, stats=None):
        if stats is not None:
            stats["F_calls"] = float(stats.get("F_calls", 0.0)) + 1.0
        row = np.zeros(50, dtype=np.float64)
        row[0] = abs(float(z[0]) - 1.25)
        return row

    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control._count_residual", fake_residual,
    )
    class _Q:
        hydro = None

    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control.outer_P_named", lambda x, nq: _Q(),
    )
    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control.certify_controlled_state",
        lambda x, f, hydro=None: {"passed": False, "transport_closed": True, "blocks": {}},
    )
    f_accepted = np.zeros(50, dtype=np.float64)
    f_accepted[0] = 0.25
    chosen, rows = composite_backtrack(None, accepted, proposed, f_accepted, stats={})
    assert chosen is not None
    assert chosen["alpha"] == pytest.approx(0.25)
    assert float(chosen["x"][0]) == pytest.approx(2.0 ** 0.25)
    assert any(row["alpha"] == 1.0 or row["alpha"] == 0.5 for row in rows)


def test_rejected_candidate_is_followed_by_a_tighter_trust(monkeypatch) -> None:
    from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
    from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x0 = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"].copy()
    x0[-2] += 50.0
    spans = []

    def fake_open(q, x, stats=None):
        return np.array(x, copy=True), []

    def fake_correct(q, x, **kwargs):
        spans.append(float(kwargs["log_span"]))
        y = np.array(x, copy=True)
        y[0] = -abs(float(kwargs["log_span"]))
        return y

    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control.reactivate_zero_unknowns", fake_open,
    )
    monkeypatch.setattr(
        "src.solvers.wirsum_1997_named_f_control._e8_log_scaled_correct", fake_correct,
    )
    out = solve_named_controlled(nq, x0, deadline_s=2.0)
    assert len(spans) >= 2
    assert spans[1] < spans[0]
    assert out["passed"] is False
    assert out["x"][0] == pytest.approx(float(x0[0]))


def test_timeout_keeps_call_counts(monkeypatch) -> None:
    from src.solvers.wirsum_1997_named_f import outer_P_named
    from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x0 = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    q = outer_P_named(x0, nq)
    stats: dict = {}
    got = _e8_log_scaled_correct(q, x0, deadline_s=0.02, stats=stats, max_j=2)
    assert stats["F_calls"] > 0.0
    assert stats["inner_status"] in ("deadline", "linear_step", "jacobian_budget")
    assert np.all(np.isfinite(got))
    assert float(stats["inner_merit_best"]) <= float(stats["inner_merit_entrance"]) + 1.0e-9
    stopped: dict = {}
    back = _e8_log_scaled_correct(q, x0, deadline_s=0.0, max_j=1, stats=stopped)
    assert stopped["inner_status"] == "deadline"
    assert stopped["checkpoint"] == "entrance"
    from src.solvers.wirsum_1997_named_f import inner_hydro_frozen_residual
    from src.solvers.wirsum_1997_named_f_e8 import residual_block_merit
    assert np.allclose(back, x0)
    assert float(stopped["inner_merit_best"]) <= float(stopped["inner_merit_entrance"]) + 1.0e-12


def test_undecided_trajectory_is_not_a_microsecond_discharge() -> None:
    from src.solvers.wirsum_1997_named_f import fb_char_balance, outer_P_named
    from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
    from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    q = outer_P_named(x, nq)
    open_hydro = replace(
        q.hydro,
        w_undecided=np.ones(2, dtype=np.float64),
        tau_s=np.zeros(2, dtype=np.float64),
        w_leave=np.zeros(2, dtype=np.float64),
        w_return=np.zeros(2, dtype=np.float64),
    )
    q_open = replace(q, hydro=open_hydro)
    resid_open, _, _ = fb_char_balance(x, q_open, rate_scale=1.0)
    floor_hydro = replace(open_hydro, w_undecided=np.zeros(2), tau_s=np.full(2, 1.0e-6))
    resid_floor, _, _ = fb_char_balance(x, replace(q, hydro=floor_hydro), rate_scale=1.0)
    assert not np.allclose(resid_open, resid_floor)
    f = np.zeros_like(x)
    assert certify_controlled_state(x, f, open_hydro)["passed"] is False
    assert certify_controlled_state(x, f, open_hydro)["transport_closed"] is False
