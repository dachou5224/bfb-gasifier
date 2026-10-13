"""未判定轨迹的步长续积分。窗口长度不能当成停留时间。"""

import math

import pytest

from src.physics.wirsum_1997_freeboard import (
    accept_refined_event,
    continue_undecided_event,
    integrate_particle_fate_with_time,
)


def test_both_undecided_with_close_states_is_not_a_step_mismatch() -> None:
    from src.physics.wirsum_1997_freeboard import classify_refined_pair

    kind = classify_refined_pair(
        ("undecided", 5.46, 3.184560, 0.2915169),
        ("undecided", 5.46, 3.183977, 0.2915169),
    )
    assert kind == "both_window_exhausted"
    assert classify_refined_pair(
        ("reaches_next_cell", 0.030, 1.0, 0.3),
        ("reaches_next_cell", 0.033, 1.0, 0.3),
        tol_t_abs_s=0.004,
    ) == "accepted"
    assert classify_refined_pair(
        ("reaches_next_cell", 0.030, 1.0, 0.3),
        ("reaches_next_cell", 0.034, 1.0, 0.3),
    ) == "event_time_disagreement"
    assert classify_refined_pair(
        ("returns", 1.0, 0.0, -0.1),
        ("reaches_next_cell", 1.0, 1.0, 0.2),
    ) == "fate_disagreement"


def test_event_is_accepted_only_when_both_steps_agree() -> None:
    assert accept_refined_event("returns", 1.0, "returns", 1.01) == ("returns", 1.01)
    assert accept_refined_event("returns", 1.0, "reaches_next_cell", 1.0) is None
    assert accept_refined_event("returns", 1.0, "returns", 2.0) is None
    assert accept_refined_event("undecided", 5.0, "undecided", 5.0) is None


def test_internal_face_event_continues_into_the_next_cell(monkeypatch) -> None:
    import src.physics.wirsum_1997_freeboard as fb

    calls = []

    def initial(*args, **kwargs):
        calls.append(float(kwargs["cell_height"]))
        return ("undecided", 5.0, 0.9, 0.2)

    def continued(*args, **kwargs):
        return ("reaches_next_cell", 1.0, 1.05, 0.4, True, "accepted")

    monkeypatch.setattr(fb, "integrate_particle_fate_with_time", initial)
    monkeypatch.setattr(fb, "continue_undecided_event", continued)
    trace = []
    fate, _t, cell = fb.integrate_particle_through_cells(
        1.0, [1.0, 1.0],
        u_g=0.3, rho_g=1.0, rho_p=1000.0, d_p=0.001, mu_g=1.0e-5, z=1.0,
        t_max=5.0, continue_budget_s=15.0, trace=trace,
    )
    assert len(calls) >= 2
    assert trace[-1]["cell"] == 1
    assert trace[-1]["reason"] == "both_window_exhausted"
    assert "u_m_s" in trace[-1]
    assert fate == "undecided"


def test_runtime_cap_stops_without_using_height_over_gas_speed() -> None:
    from src.physics.wirsum_1997_freeboard import (
        drag_equilibrium_speed_m_s,
        drag_work_steps,
        integrate_particle_through_cells,
    )

    u_eq = drag_equilibrium_speed_m_s(
        u_g=0.5, u_p=0.2, rho_g=1.0, rho_p=1400.0, d_p=3.0e-4, mu_g=4.0e-5, z=0.5,
    )
    assert math.isfinite(u_eq)
    steps = drag_work_steps(2.0, 0.2, u_eq, 1.0e-3)
    assert isinstance(steps, int) and steps > 0
    trace: list = []
    fate, t_s, _cell = integrate_particle_through_cells(
        0.2, [5.0],
        u_g=0.5, rho_g=1.0, rho_p=1400.0, d_p=3.0e-4, mu_g=4.0e-5, z=0.5,
        dt=1.0e-3, t_max=0.01, continue_budget_s=0.0, extra_physical_s=0.02,
        trace=trace,
    )
    assert fate == "undecided"
    assert trace[-1]["reason"] == "physical_budget_incomplete"
    assert trace[-1]["used_physical_s"] > 0.0
    assert trace[-1]["used_work_steps"] > 0
    assert "wall_elapsed_s" in trace[-1]
    assert trace[-1]["t_s"] == pytest.approx(t_s)
    assert t_s < 0.05
    assert t_s != pytest.approx(5.0 / 0.5)
    assert "dt_mismatch" not in trace[-1]["reason"]
    assert trace[-1]["drag_work_steps"] > 0
    assert "coarse" in trace[-1] and "fine" in trace[-1]


def test_work_steps_and_wall_clock_stop_separately() -> None:
    from src.physics.wirsum_1997_freeboard import integrate_particle_through_cells

    kwargs = dict(
        u_g=0.5, rho_g=1.0, rho_p=1400.0, d_p=3.0e-4, mu_g=4.0e-5, z=0.5,
        dt=1.0e-3, t_max=0.01, continue_budget_s=0.0, extra_physical_s=5.0,
    )
    trace: list = []
    fate, _t, _cell = integrate_particle_through_cells(0.2, [5.0], max_work_steps=3, trace=trace, **kwargs)
    assert fate == "undecided"
    assert trace[-1]["reason"] == "work_budget_incomplete"
    assert trace[-1]["used_work_steps"] <= 3
    trace_wall: list = []
    fate_w, _tw, _c = integrate_particle_through_cells(
        0.2, [5.0], wall_s=0.0, max_work_steps=0, extra_physical_s=0.0, trace=trace_wall, **{k: kwargs[k] for k in kwargs if k != "extra_physical_s"},
    )
    assert trace_wall[-1]["reason"] == "both_window_exhausted"
    _ = fate_w


def test_agreed_rising_particle_reaches_the_top_inside_the_runtime_cap() -> None:
    from src.physics.wirsum_1997_freeboard import integrate_particle_through_cells

    trace: list = []
    fate, t_s, cell = integrate_particle_through_cells(
        2.0, [0.2],
        u_g=3.0, rho_g=1.0, rho_p=500.0, d_p=1.0e-4, mu_g=2.0e-5, z=0.44,
        dt=1.0e-3, t_max=0.01, continue_budget_s=0.0, extra_physical_s=2.0,
        trace=trace,
    )
    assert fate == "reaches_next_cell"
    assert trace[-1]["reason"] == "furnace_exit"
    assert cell == 0
    assert t_s > 0.01
    assert trace[-1]["height_m"] >= 0.2


def test_short_window_stays_undecided_until_a_refined_event() -> None:
    kwargs = dict(
        u_g=0.2, rho_g=1.0, rho_p=1400.0, d_p=3.0e-4, mu_g=4.0e-5, z=0.5, cell_height=2.0,
    )
    short = integrate_particle_fate_with_time(0.3, dt=1.0e-3, t_max=0.02, **kwargs)
    assert short[0] == "undecided"
    fate, t_add, z_end, u_end, ok, kind = continue_undecided_event(
        short[3], short[2], dt=1.0e-3, extra_s=2.0, **kwargs,
    )
    assert ok is True
    assert kind == "accepted"
    assert fate == "returns"
    assert t_add > 0.02
    assert t_add != pytest.approx(2.0)
    assert z_end <= 0.0
    assert u_end <= 0.0
