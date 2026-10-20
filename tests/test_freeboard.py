"""自由空间抛出、气速与轨迹。期望值按冻结稿和 ADR 手算。"""

from __future__ import annotations

import math

import numpy as np
import pytest
from dataclasses import replace

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError,
    approach_factor,
    closed_form_state,
    coalescence_weight,
    drag_coefficient,
    drag_limit_z,
    ejection_flux,
    ejection_mass_flow,
    ejection_velocity_samples,
    ghost_bubble_velocity,
    integrate_freeboard,
    integrate_motion,
    motion_coefficients,
    pair_coalescence_probability,
    particle_acceleration,
    spacing_ratio,
    terminal_settling_velocity,
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


def test_spacing_and_approach_by_hand():
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


def test_numerical_integral_matches_closed_form_when_delta_positive():
    # 大颗粒、静气、向上的初速：这一支的 Δ > 0，短时间内速度保持为正。
    diameter = 8.0e-3
    rho_g = 1.1
    mu = 1.8e-5
    rho_p = 900.0
    gas = 0.0
    initial = 2.0
    a, b, c = motion_coefficients(gas, True, rho_g, mu, diameter, rho_p, PHI)
    assert a * initial**2 + b * initial + c == pytest.approx(
        particle_acceleration(initial, gas, rho_g, mu, diameter, rho_p, PHI)
    )
    delta = 4.0 * a * c - b**2
    assert delta > 0.0
    duration = 0.02
    analytic_u, analytic_h = closed_form_state(duration, a, b, c, initial, 0.3)
    numeric_h, numeric_u = integrate_motion(
        0.3, initial, duration, gas, rho_g, mu, diameter, rho_p, PHI, 1.0e-5
    )
    assert numeric_u == pytest.approx(analytic_u, rel=1e-4, abs=1e-6)
    assert numeric_h == pytest.approx(analytic_h, rel=1e-4, abs=1e-6)
    assert numeric_u > gas



def _column():
    edges = np.array([0.0, 0.5, 1.0])
    rho = np.full(2, RHO_G)
    mu = np.full(2, MU)
    return edges, rho, mu


def test_fast_gas_carries_the_particle_out_and_slow_gas_returns_it():
    edges, rho, mu = _column()
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)

    def fast(_cell, _height, terminal=terminal):
        return terminal + 3.0

    fate, _time, height, velocity = integrate_freeboard(
        0.0, 1.0, edges, fast, rho, mu, diameter, RHO_P, PHI, SolverConfig()
    )
    assert fate == "elutriate"
    assert height == pytest.approx(1.0)
    assert velocity > 0.0

    def slow(_cell, _height, terminal=terminal):
        return 0.2 * terminal

    fate, _time, height, velocity = integrate_freeboard(
        0.2, 0.3, edges, slow, rho, mu, diameter, RHO_P, PHI, SolverConfig()
    )
    assert fate == "return"
    assert height == pytest.approx(0.0)
    assert velocity < 0.0


def test_time_limit_is_an_error_not_a_fate():
    edges, rho, mu = _column()

    def gas(_cell, _height):
        return 0.0

    with pytest.raises(FreeboardError, match="积分时间达到上限"):
        integrate_freeboard(
            0.4,
            0.0,
            edges,
            gas,
            rho,
            mu,
            2.0e-4,
            RHO_P,
            PHI,
            _config(freeboard_time_limit_s=1.0e-3),
        )


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def test_step_halving_changes_the_state_by_less_than_1e_4():
    edges, rho, mu = _column()
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)
    base_step = SolverConfig().freeboard_max_step_s

    def gas(_cell, _height, terminal=terminal):
        return terminal + 3.0

    base = integrate_freeboard(0.0, 2.0, edges, gas, rho, mu, diameter, RHO_P, PHI, SolverConfig())
    fine = integrate_freeboard(
        0.0,
        2.0,
        edges,
        gas,
        rho,
        mu,
        diameter,
        RHO_P,
        PHI,
        _config(freeboard_max_step_s=0.5 * base_step),
    )
    assert base[0] == fine[0] == "elutriate"
    assert _relative(base[1], fine[1]) < 1.0e-4
    assert _relative(base[3], fine[3]) < 1.0e-4


def test_doubling_ejection_samples_changes_the_mean_flight_time_by_less_than_1e_4():
    edges, rho, mu = _column()
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)
    count = SolverConfig().freeboard_ejection_samples

    def gas(_cell, _height, terminal=terminal):
        return terminal + 3.0

    def mean_time(sample_count: int) -> float:
        speeds, shares, _direct = ejection_velocity_samples(1.1, sample_count)
        total = 0.0
        weight = 0.0
        for speed, share in zip(speeds, shares):
            _fate, time, _height, _velocity = integrate_freeboard(
                0.0, float(speed), edges, gas, rho, mu, diameter, RHO_P, PHI, SolverConfig()
            )
            total += float(share) * time
            weight += float(share)
        return total / weight

    assert _relative(mean_time(count), mean_time(2 * count)) < 1.0e-4


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
