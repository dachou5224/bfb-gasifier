"""常物性自由空间轨迹的解析边界事件。

对照 docs/audit_reports/trajectory_hamel_review_2026-10-09.md 的独立积分。
窗口末时刻必须和状态属于同一步，不把 120 s 当成物理边界。
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from src.physics.wirsum_1997_freeboard import (
    constant_property_boundary_event,
    hamel_eq386_asymptotic_z,
    integrate_particle_fate_with_time,
    integrate_particle_through_cells,
    particle_acceleration_cw,
)

# 审计报告 §3 的冻结输入
_UG = 0.4550760  # [m/s]
_RHO = 5.51522  # [kg/m³]
_RP = 1400.0  # [kg/m³]
_D = 2.0e-4  # [m]
_MU = 4.0e-5  # [Pa·s]
_H = 4.26880  # [m]
_Z = float(hamel_eq386_asymptotic_z(0.8))


def _props(u_g: float = _UG) -> dict[str, float]:
    return dict(u_g=u_g, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z)


def _kw(u_g: float = _UG) -> dict[str, float]:
    return dict(**_props(u_g), cell_height=_H)


def test_time_window_returns_the_state_at_that_instant() -> None:
    """t_max 不是 dt 的整数倍时，两个窗口不能共用同一步的末态。"""
    short = integrate_particle_fate_with_time(1.34, dt=0.01, t_max=0.015, **_kw())
    long = integrate_particle_fate_with_time(1.34, dt=0.01, t_max=0.020, **_kw())
    assert short[0] == long[0] == "undecided"
    assert short[1] == pytest.approx(0.015)
    assert long[1] == pytest.approx(0.020)
    assert long[2] == pytest.approx(0.0212412032, abs=1.0e-9)
    assert long[3] == pytest.approx(0.8429500451, abs=1.0e-8)
    assert short[2] != pytest.approx(long[2])
    assert short[2] < long[2]


@pytest.mark.parametrize(
    ("u0", "t_ref"),
    (
        (0.89, 27.8459),
        (1.34, 41.3357),
        (2.24, 60.5344),
        (4.47, 87.4294),
    ),
)
def test_uniform_field_return_matches_audit_table(u0: float, t_ref: float) -> None:
    fate, t_ev, z_end, u_end, reason = constant_property_boundary_event(u0, **_kw())
    assert fate == "returns"
    assert reason == "bed_return"
    assert z_end == pytest.approx(0.0, abs=1.0e-12)
    assert u_end < 0.0
    assert t_ev == pytest.approx(t_ref, abs=2.0e-3)


def test_gas_speed_shift_keeps_finite_events_past_120_s() -> None:
    back_slow = constant_property_boundary_event(1.34, **_kw(_UG + 0.0009))
    back_long = constant_property_boundary_event(1.34, **_kw(_UG + 0.0013))
    leave = constant_property_boundary_event(1.34, **_kw(_UG + 0.0020))
    assert back_slow[0] == "returns"
    assert back_slow[1] == pytest.approx(114.120, abs=0.01)
    assert back_long[0] == "returns"
    assert back_long[1] == pytest.approx(526.470, abs=0.05)
    assert leave[0] == "reaches_next_cell"
    assert leave[4] == "upper_boundary"
    assert leave[1] == pytest.approx(7144.611, abs=0.1)
    assert leave[2] == pytest.approx(_H)
    windowed = integrate_particle_fate_with_time(
        1.34, dt=2.0e-3, t_max=5.0, resolve_events=False, **_kw(_UG + 0.0013),
    )
    assert windowed[0] == "undecided"
    assert windowed[1] == pytest.approx(5.0)
    resolved = integrate_particle_fate_with_time(
        1.34, dt=2.0e-3, t_max=5.0, resolve_events=True, **_kw(_UG + 0.0013),
    )
    assert resolved[0] == "returns"
    assert resolved[1] == pytest.approx(526.470, abs=0.05)


def test_return_time_matches_independent_force_integration() -> None:
    fate, t_ev, _z, _u, _reason = constant_property_boundary_event(1.34, **_kw())

    def rhs(_t: float, y: np.ndarray) -> list[float]:
        return [float(y[1]), particle_acceleration_cw(float(y[1]), **_props())]

    def bed(_t: float, y: np.ndarray) -> float:
        return float(y[0])

    bed.terminal = True
    bed.direction = -1
    sol = solve_ivp(
        rhs, (0.0, 80.0), np.array([0.0, 1.34]), method="Radau",
        events=bed, rtol=1.0e-8, atol=1.0e-10,
    )
    assert fate == "returns"
    assert sol.t_events[0].size == 1
    assert t_ev == pytest.approx(float(sol.t_events[0][0]), abs=1.0e-3)


def test_launch_above_gas_speed_crosses_branch_before_return() -> None:
    """初速高于气速时先过 u=u_g，回床仍是床面事件。"""
    started_below = constant_property_boundary_event(0.2, **_kw())
    started_above = constant_property_boundary_event(1.34, **_kw())
    assert started_below[0] == started_above[0] == "returns"
    assert started_below[1] < started_above[1]
    assert started_above[1] > 1.0


def test_zero_terminal_speed_stays_inside_the_cell() -> None:
    sample = constant_property_boundary_event(1.34, **_kw())
    u_eq = float(sample[3])
    stalled = constant_property_boundary_event(1.34, **_kw(_UG - u_eq))
    assert stalled[0] == "undecided"
    assert stalled[4] == "asymptotic_residence"
    assert 0.0 < stalled[2] < _H
    assert abs(stalled[3]) < 1.0e-6


def test_upper_cell_return_reaches_the_bed() -> None:
    one = constant_property_boundary_event(1.34, **_kw())
    chained, t_chain, cell = integrate_particle_through_cells(
        1.34, [_H / 2.0, _H / 2.0], resolve_events=True, **_props(),
    )
    assert chained == "returns"
    assert cell == 0
    assert t_chain == pytest.approx(one[1], abs=1.0e-9)
    crossed, t_cross, cell_cross = integrate_particle_through_cells(
        1.34, [0.02, 2.0], resolve_events=True, **_props(),
    )
    assert crossed == "returns"
    assert cell_cross == 0
    assert t_cross == pytest.approx(one[1], abs=1.0e-6)
