"""自由空间抛出、气速、逐格解析轨迹与事件。期望值按冻结稿和 ADR 手算。"""

from __future__ import annotations

import math

import numpy as np
import pytest
from dataclasses import replace

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError,
    advance_cell,
    approach_factor,
    cell_height_at_velocity,
    cell_mean_velocity,
    cell_state,
    coalescence_weight,
    discriminant,
    drag_coefficient,
    drag_limit_z,
    ejection_flux,
    ejection_mass_flow,
    ejection_velocity_samples,
    ghost_bubble_velocity,
    motion_coefficients,
    pair_coalescence_probability,
    particle_acceleration,
    spacing_ratio,
    terminal_settling_velocity,
    time_to_velocity,
    trace_until_fate,
)

G = 9.81
PHI = 0.8
RHO_G = 1.2
MU = 2.0e-5
RHO_P = 1200.0


def _hand_z(phi_s: float) -> float:
    return 73.69 * math.exp(-5.0784 * phi_s)


def _config(**changes) -> SolverConfig:
    return replace(SolverConfig(), **changes)


def _cellwise(*speeds: float):
    def gas(cell, _height):
        return float(speeds[cell])

    return gas


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def _rk4(a: float, b: float, c: float, u0: float, h0: float, duration: float, steps: int = 200000):
    """常系数参考积分。只作闭式的对照，不进入求解。求解手段：ADR 0047。"""
    step = duration / steps
    velocity, height = u0, h0
    f = lambda value: a * value * value + b * value + c
    for _ in range(steps):
        k1 = f(velocity)
        k2 = f(velocity + 0.5 * step * k1)
        k3 = f(velocity + 0.5 * step * k2)
        k4 = f(velocity + step * k3)
        height += step * (velocity + 2.0 * (velocity + 0.5 * step * k1) + 2.0 * (velocity + 0.5 * step * k2) + velocity + step * k3) / 6.0
        velocity += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return velocity, height


def test_spacing_and_approach_by_hand():
    """式 (3.57)。解释性选择：ADR 0028。不用式 (3.58)。"""
    fraction = 0.25
    ratio = (6.0 * fraction / math.pi) ** (-1.0 / 3.0)
    assert spacing_ratio(fraction) == pytest.approx(ratio)
    assert approach_factor() == pytest.approx(1.0 - math.cos(math.radians(48.0)))
    prime = 0.5 * (1.0 + math.erf((1.0 - ratio) / (0.4 * math.sqrt(2.0))))
    assert pair_coalescence_probability(fraction) == pytest.approx(approach_factor() * prime)


def test_coalescence_weight_and_ejection_flux_by_hand():
    fraction = 0.2
    probability = pair_coalescence_probability(fraction)
    weight = 0.0
    for count in range(2, 6):
        weight += (count - 1) * (1.0 - probability) ** (count - 1) * probability
    assert coalescence_weight(fraction) == pytest.approx(weight)
    bubble_velocity = 1.4
    bubble_diameter = 0.08
    mean_diameter = 0.7e-3
    voidage = 0.55
    density = 1100.0
    roof = 3.0 * mean_diameter / bubble_diameter
    bracket = roof + (0.4 * 0.25 - roof) * weight
    flux = fraction * bubble_velocity * density * (1.0 - voidage) * bracket
    assert ejection_flux(
        fraction, bubble_velocity, bubble_diameter, density, voidage, mean_diameter
    ) == pytest.approx(flux)
    area = 0.017671
    assert ejection_mass_flow(
        fraction, bubble_velocity, bubble_diameter, density, voidage, mean_diameter, area
    ) == pytest.approx(flux * area)


def test_negative_ejection_speed_is_direct_return():
    speeds, shares, direct = ejection_velocity_samples(1.2, 8)
    mean = 1.53 * 1.2
    sigma = 0.6 * 1.2
    expected = 0.5 * (1.0 + math.erf((0.0 - mean) / (sigma * math.sqrt(2.0))))
    assert direct == pytest.approx(expected)
    assert speeds.min() > 0.0
    assert shares.sum() + direct == pytest.approx(1.0)


def test_drag_and_acceleration_by_hand():
    z = _hand_z(PHI)
    assert drag_limit_z(PHI) == pytest.approx(z)
    assert drag_coefficient(12.5, PHI) == pytest.approx(24.0 / 12.5 + z)
    diameter = 4.0e-4
    gas = 2.0
    particle = 0.4
    relative = particle - gas
    stokes = 18.0 * MU / (RHO_P * diameter**2)
    quadratic = 3.0 * RHO_G * z / (4.0 * RHO_P * diameter)
    weight = (RHO_P - RHO_G) * G / RHO_P
    expected = -stokes * relative - quadratic * relative * abs(relative) - weight
    assert particle_acceleration(particle, gas, RHO_G, MU, diameter, RHO_P, PHI) == pytest.approx(expected)
    settling = (-stokes + math.sqrt(stokes**2 + 4.0 * quadratic * weight)) / (2.0 * quadratic)
    assert terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI) == pytest.approx(settling)
    assert particle_acceleration(-settling, 0.0, RHO_G, MU, diameter, RHO_P, PHI) == pytest.approx(0.0, abs=1e-9)


def test_polynomial_matches_acceleration_on_both_sides():
    diameter = 8.0e-4
    gas = 1.7
    for particle, faster in ((2.4, True), (0.2, False)):
        a, b, c = motion_coefficients(gas, faster, RHO_G, MU, diameter, RHO_P, PHI)
        polynomial = a * particle**2 + b * particle + c
        direct = particle_acceleration(particle, gas, RHO_G, MU, diameter, RHO_P, PHI)
        assert polynomial == pytest.approx(direct)
    # 两支在 u_P = u_g 处给出同一个加速度 −(ρ_P−ρ_g)g/ρ_P，逐格分界连续。
    slower = motion_coefficients(gas, False, RHO_G, MU, diameter, RHO_P, PHI)
    faster = motion_coefficients(gas, True, RHO_G, MU, diameter, RHO_P, PHI)
    at_gas = lambda coefficient: coefficient[0] * gas**2 + coefficient[1] * gas + coefficient[2]
    assert at_gas(slower) == pytest.approx(at_gas(faster))
    assert at_gas(slower) == pytest.approx(-(RHO_P - RHO_G) * G / RHO_P)


def test_cell_mean_velocity_is_the_height_average():
    def gas(_cell, height):
        return 1.5 + 0.8 * float(height)

    left, right = 0.2, 0.6
    assert cell_mean_velocity(gas, 0, left, right) == pytest.approx(1.5 + 0.8 * 0.5 * (left + right))


def test_three_closed_forms_match_rk4_within_1e_8():
    """Δ > 0、Δ < 0、Δ = 0 三种闭式对同一逐格常系数的高精度 RK4。"""
    cases = [
        (-2.0, -30.0, -150.0, 6.0, 0.02),
        (-2.0, -40.0, 50.0, 0.5, 0.2),
        (2.0, -40.0, -10.0, -0.5, 0.2),
    ]
    for a, b, c, u0, duration in cases:
        closed_u, closed_h = cell_state(duration, a, b, c, u0, 0.0)
        numeric_u, numeric_h = _rk4(a, b, c, u0, 0.0, duration)
        assert _relative(closed_u, numeric_u) < 1e-8
        assert _relative(closed_h, numeric_h) < 1e-8
    limit = cell_state(0.1, 2.0, 0.0, 0.0, 3.0, 0.0)
    numeric = _rk4(2.0, 0.0, 0.0, 3.0, 0.0, 0.1)
    assert _relative(limit[0], numeric[0]) < 1e-8
    assert _relative(limit[1], numeric[1]) < 1e-8


def test_near_zero_discriminant_uses_the_limit():
    a, b, c, u0, duration = 2.0, 40.0, 200.0 + 1.0e-12, -20.0, 0.05
    assert abs(discriminant(a, b, c)) <= 1.0e-10 * b * b
    closed_u, closed_h = cell_state(duration, a, b, c, u0, 0.0)
    numeric_u, numeric_h = _rk4(a, b, c, u0, 0.0, duration)
    assert _relative(closed_u, numeric_u) < 1e-6
    assert _relative(closed_h, numeric_h) < 1e-6


def test_velocity_crossing_times_are_analytic():
    a, b, c = motion_coefficients(1.2, True, RHO_G, MU, 2.0e-4, RHO_P, PHI)
    u0 = 2.0
    target = 0.3
    tau = time_to_velocity(a, b, c, u0, target)
    assert tau is not None and tau > 0.0
    assert cell_state(tau, a, b, c, u0, 0.0)[0] == pytest.approx(target, rel=1e-10)
    assert time_to_velocity(a, b, c, u0, 1.0e4) is None

    a, b, c = motion_coefficients(1.0, False, RHO_G, MU, 2.0e-4, RHO_P, PHI)
    tau = time_to_velocity(a, b, c, -0.5, 0.0)
    assert tau is not None and tau > 0.0
    assert cell_state(tau, a, b, c, -0.5, 0.0)[0] == pytest.approx(0.0, abs=1e-12)


def test_advance_cell_reports_the_earliest_event():
    a, b, c = motion_coefficients(1.2, True, RHO_G, MU, 2.0e-4, RHO_P, PHI)
    kind, tau, height, velocity = advance_cell(0.2, 2.0, 0.0, 1.0, a, b, c, 1.2, 80.0)
    assert kind == "gas"
    assert velocity == pytest.approx(1.2)
    assert 0.0 < height < 1.0 and tau > 0.0

    # Δ > 0、A > 0：负速度向上加速，穿过 0 的时刻早于到达格界。
    kind, tau, height, velocity = advance_cell(0.5, -0.5, 0.0, 1.0, 2.0, -40.0, 50.0, 10.0, 80.0)
    assert kind == "zero"
    assert velocity == 0.0 and tau > 0.0
    assert height < 0.5


def test_cell_boundary_root_is_within_1e_12_m():
    kind, _tau, height, _velocity = advance_cell(0.1, 0.5, 0.0, 1.0, -2.0, -40.0, 50.0, 1.2, 80.0)
    assert kind == "upper"
    assert abs(height - 1.0) <= 1.0e-12


def test_fast_gas_carries_the_particle_out_and_slow_gas_returns_it():
    edges = np.array([0.0, 0.5, 1.0])
    rho = np.full(2, RHO_G)
    mu = np.full(2, MU)
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)

    fast = _cellwise(terminal + 3.0, terminal + 3.0)
    fate = trace_until_fate(0.0, 1.0, edges, fast, rho, mu, diameter, RHO_P, PHI, SolverConfig())
    assert fate.fate == "elutriate"
    assert fate.height_m == pytest.approx(1.0)
    assert fate.velocity_m_s > 0.0

    slow = _cellwise(0.2 * terminal, 0.2 * terminal)
    fate = trace_until_fate(0.2, 0.3, edges, slow, rho, mu, diameter, RHO_P, PHI, SolverConfig())
    assert fate.fate == "return"
    assert fate.height_m == pytest.approx(0.0)
    assert fate.velocity_m_s < 0.0


def test_boundary_event_inside_the_time_limit_is_not_missed():
    """没有速度事件时，格界求根也必须受剩余期限约束，但期限内的真实事件不能漏。"""
    terminal = terminal_settling_velocity(RHO_G, MU, 2.0e-4, RHO_P, PHI)
    edges = np.array([0.0, 0.75])
    rho, mu = np.array([RHO_G]), np.array([MU])

    def gas(_cell, _height, terminal=terminal):
        return terminal + 1.0

    fate = trace_until_fate(
        0.0, 1.0, edges, gas, rho, mu, 2.0e-4, RHO_P, PHI, _config(freeboard_time_limit_s=1.0)
    )
    assert fate.fate == "elutriate"
    assert fate.time_s == pytest.approx(0.75, rel=1e-6)
    assert fate.height_m == pytest.approx(0.75, abs=1e-12)


def test_event_beyond_the_time_limit_is_not_accepted():
    """期限之后的格界事件不构成归宿；到时只报错，不改判为带出。"""
    edges = np.array([0.0, 0.01])
    rho, mu = np.array([RHO_G]), np.array([MU])

    def gas(_cell, _height):
        return 1.0

    with pytest.raises(FreeboardError, match="积分时间达到上限"):
        trace_until_fate(
            0.0, 2.0, edges, gas, rho, mu, 2.0e-4, RHO_P, PHI, _config(freeboard_time_limit_s=1.0e-6)
        )


def test_time_limit_is_an_error_not_a_fate():
    edges = np.array([0.0, 0.5, 1.0])
    rho = np.full(2, RHO_G)
    mu = np.full(2, MU)
    still = _cellwise(0.0, 0.0)
    with pytest.raises(FreeboardError, match="积分时间达到上限"):
        trace_until_fate(
            0.4, 0.0, edges, still, rho, mu, 2.0e-4, RHO_P, PHI, _config(freeboard_time_limit_s=1.0e-3)
        )


def test_doubling_ejection_samples_changes_the_mean_flight_time_by_less_than_1e_4():
    edges = np.array([0.0, 0.5, 1.0])
    rho = np.full(2, RHO_G)
    mu = np.full(2, MU)
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)
    gas = _cellwise(terminal + 3.0, terminal + 3.0)
    count = SolverConfig().freeboard_ejection_samples

    def mean_time(sample_count: int) -> float:
        speeds, shares, _direct = ejection_velocity_samples(1.1, sample_count)
        total = 0.0
        weight = 0.0
        for speed, share in zip(speeds, shares):
            fate = trace_until_fate(0.0, float(speed), edges, gas, rho, mu, diameter, RHO_P, PHI, SolverConfig())
            total += float(share) * fate.time_s
            weight += float(share)
        return total / weight

    assert _relative(mean_time(count), mean_time(2 * count)) < 1.0e-4


def test_stationary_boundary_uses_acceleration_and_not_every_zero_speed_leaves():
    edges = np.array([0.0, 1.0])
    rho = np.array([RHO_G])
    mu = np.array([MU])
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)
    still = _cellwise(0.0)
    lifting = _cellwise(terminal + 3.0)
    balanced = _cellwise(terminal)

    fate = trace_until_fate(0.0, 0.0, edges, still, rho, mu, diameter, RHO_P, PHI, _config(freeboard_time_limit_s=1.0e-6))
    assert (fate.fate, fate.time_s, fate.height_m, fate.velocity_m_s) == ("return", 0.0, 0.0, 0.0)
    fate = trace_until_fate(1.0, 0.0, edges, lifting, rho, mu, diameter, RHO_P, PHI, SolverConfig())
    assert (fate.fate, fate.time_s, fate.height_m, fate.velocity_m_s) == ("elutriate", 0.0, 1.0, 0.0)
    fate = trace_until_fate(0.0, 0.0, edges, lifting, rho, mu, diameter, RHO_P, PHI, SolverConfig())
    assert fate.fate == "elutriate"
    assert fate.time_s > 0.0
    assert fate.height_m == pytest.approx(1.0)
    fate = trace_until_fate(1.0, 0.0, edges, still, rho, mu, diameter, RHO_P, PHI, SolverConfig())
    assert fate.fate == "return"
    assert fate.time_s > 0.0
    assert fate.height_m == pytest.approx(0.0)
    assert fate.velocity_m_s < 0.0
    with pytest.raises(FreeboardError, match="中性"):
        trace_until_fate(
            0.0, 0.0, edges, balanced, rho, mu, diameter, RHO_P, PHI, _config(freeboard_time_limit_s=1.0e-3)
        )


def test_ghost_bubble_uses_local_flow_and_area():
    edges = np.array([0.0, 0.2875, 2.3])
    diameters = np.array([0.15, 0.25, 0.25])
    primary = 0.04
    total = 0.07
    flows = np.array([primary, total])
    low = ghost_bubble_velocity(0, 0.1, edges, diameters, flows, 1.6, 0.05)
    high = ghost_bubble_velocity(1, 1.0, edges, diameters, flows, 1.6, 0.05)
    area_low = math.pi * (0.15 + (0.25 - 0.15) * (0.1 / 0.2875)) ** 2 / 4.0
    approach_low = primary / area_low
    expected_low = approach_low + (1.6 - approach_low) * math.exp(-0.1 / (3.0 * 0.05))
    approach_high = total / (math.pi * 0.25**2 / 4.0)
    expected_high = approach_high + (1.6 - approach_high) * math.exp(-1.0 / (3.0 * 0.05))
    assert low == pytest.approx(expected_low)
    assert high == pytest.approx(expected_high)
