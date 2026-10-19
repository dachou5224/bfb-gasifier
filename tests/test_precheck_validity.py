# -*- coding: utf-8 -*-
"""T27 预检：悬停判据与释放时间积分。期望值是手算的解析式，不是工况 1 的输出。"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    path = ROOT / "tools" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


hover = _load("precheck_hover")
pyro = _load("precheck_pyrolysis")


def test_monotonic_decline_finds_the_stable_height():
    def velocity(height):
        return 0.4 + 2.0 * math.exp(-1.5 * height)

    grid = np.linspace(0.0, 3.0, 2001)
    found = hover.stable_heights(grid, [velocity(item) for item in grid], 1.0)
    expected = -math.log(0.3) / 1.5
    assert found == pytest.approx([expected], abs=1.0e-3)


def test_rising_crossing_is_not_stable():
    grid = np.linspace(0.0, 2.0, 401)
    velocity = 0.2 + 2.0 * grid
    assert hover.stable_heights(grid, velocity, 1.0) == []


def test_terminal_velocity_zeros_the_motion_equation():
    particle = hover.Particle(2.0e-4, 800.0, 1.2, 3.0e-5, 0.8)
    settling = hover.terminal_settling_velocity(particle)
    assert settling > 0.0
    assert hover.motion_residual(particle, -settling) == pytest.approx(0.0, abs=1.0e-8)


def test_bed_samples_miss_while_local_trajectory_approaches():
    particle = hover.Particle(2.0e-4, 800.0, 1.2, 3.0e-5, 0.8)
    settling = hover.terminal_settling_velocity(particle)

    def velocity(height):
        if height < 1.5:
            return 0.05 * settling
        slope = (0.2 * settling - 2.5 * settling) / 0.5
        return 2.5 * settling + slope * (height - 1.5)

    grid = np.linspace(0.0, 2.0, 2001)
    assessed = hover.assess_hover(
        velocity,
        grid,
        [velocity(item) for item in grid],
        [0.0, 0.5, 1.0, 1.5, 2.0],
        particle,
        [(0.02, 0.5), (0.05, 0.5)],
    )
    assert assessed.stable_heights_m
    assert assessed.bed.status == "不会"
    assert assessed.local.status == "会"
    assert assessed.overall == "会"
    assert assessed.unreachable_cells


def test_no_stable_height_is_checked_for_both_families():
    particle = hover.Particle(2.0e-4, 800.0, 1.2, 3.0e-5, 0.8)
    grid = np.linspace(0.0, 1.0, 21)
    assessed = hover.assess_hover(
        lambda _height: 3.0,
        grid,
        np.full(grid.shape, 3.0),
        [0.0, 1.0],
        particle,
        [(1.0, 1.0)],
    )
    assert assessed.bed.status == "不会"
    assert assessed.local.status == "不会"
    assert assessed.overall == "不会"


def test_release_time_matches_exponential_decay():
    moisture = 3.0
    volatile = 1.5
    grid = np.linspace(0.0, 40.0, 8001)
    water = np.exp(-grid / moisture)
    vapours = np.exp(-grid / volatile)
    water_time = pyro.release_time(grid, water, pyro.exponential_tail(water[-1], moisture))
    volatile_time = pyro.release_time(grid, vapours, pyro.exponential_tail(vapours[-1], volatile))
    assert water_time == pytest.approx(moisture, rel=1.0e-6)
    assert volatile_time == pytest.approx(volatile, rel=1.0e-6)
    assert volatile_time < water_time


def test_front_rate_at_the_surface_matches_the_flux_limit():
    particle = pyro.case1_properties()
    radius = 5.0e-4
    rate = pyro.front_rate(1.0, particle, radius)
    alpha = pyro.heat_transfer_coefficient(particle, 2.0 * radius)
    expected = -alpha * (particle.temperature_k - particle.evaporation_temperature_k) / (
        radius * particle.density_kg_m3 * particle.corrected_latent_heat * particle.moisture_dry_basis
    )
    assert rate == pytest.approx(expected, rel=1.0e-12)
    assert rate < 0.0


def test_uniform_remaining_fraction_averages_to_one():
    radius = np.linspace(0.0, 1.2e-3, 25)
    assert pyro.average_volatile(radius, np.ones_like(radius), radius[-1]) == pytest.approx(1.0, rel=1.0e-6)


def test_fuel_conversion_conserves_elements_and_gasifies_char():
    converted = hover.fuel_to_gas()
    assert converted["char_left_mol_s"] == pytest.approx(0.0, abs=1.0e-8)
    for name, expected in converted["expected_atoms"].items():
        assert converted["atoms"][name] == pytest.approx(expected, rel=1.0e-9, abs=1.0e-12)
    lower = hover.case1_gas()
    upper = hover.case1_gas_upper()
    assert upper["primary_mol_s"] > lower["primary_mol_s"] + 0.4
    assert hover.superficial_velocity(2.0, upper) > hover.superficial_velocity(2.0, lower)


def test_adr_0023_diameters_match_the_sauter_mean():
    assert hover.fine_lower_mm() == pytest.approx(0.0291, abs=5.0e-5)
    published = (
        (0.0291, 0.1, 0.0506),
        (0.1, 0.25, 0.1527),
        (0.25, 0.5, 0.3466),
        (0.5, 0.71, 0.5928),
        (0.71, 1.0, 0.8385),
        (1.0, 2.0, 1.3863),
        (2.0, 3.15, 2.4885),
    )
    for low, high, diameter in published:
        assert hover.sauter_diameter_mm(low, high) == pytest.approx(diameter, abs=5.0e-4)


def test_declining_profile_defines_a_hover_band():
    def velocity(height):
        return 0.4 + 2.0 * math.exp(-1.5 * height)

    spans = hover.hover_velocity_spans(velocity, [(0.0, 3.0)])
    assert len(spans) == 1
    low, high = spans[0]
    assert low == pytest.approx(velocity(3.0), rel=1.0e-3)
    assert high == pytest.approx(velocity(0.0), rel=1.0e-3)
    assert hover.hover_velocity_spans(lambda height: 0.2 + 2.0 * height, [(0.0, 2.0)]) == []
    particle = hover.Particle(2.0e-4, 800.0, 1.2, 3.0e-5, 0.8)
    target = hover.terminal_settling_velocity(particle)
    found = hover.diameter_for_settling(target, 800.0, 1.2, 3.0e-5, 0.8)
    assert found == pytest.approx(2.0e-4, rel=1.0e-4)


def test_each_density_combination_keeps_its_own_verdict():
    rows = [
        {
            "label": "未检验：没有下限",
            "diameter_m": None,
            "epsilon_mf": 0.4,
            "phi_s": 0.6,
            "rho_p": math.nan,
            "v_t": math.nan,
            "stable_m": [],
            "bed": "未检验",
            "local": "未检验",
            "overall": "未检验",
        },
        {
            "label": "筛孔 1 mm",
            "diameter_m": 1.0e-3,
            "epsilon_mf": 0.40,
            "phi_s": 0.60,
            "rho_p": 900.0,
            "v_t": 1.2,
            "stable_m": [],
            "bed": "不会",
            "local": "不会",
            "overall": "不会",
        },
        {
            "label": "筛孔 1 mm",
            "diameter_m": 1.0e-3,
            "epsilon_mf": 0.55,
            "phi_s": 0.90,
            "rho_p": 1400.0,
            "v_t": 1.5,
            "stable_m": [2.0],
            "bed": "会",
            "local": "不会",
            "overall": "会",
        },
    ]
    text = hover.format_combinations(hover.combination_lines(rows))
    assert text.count("筛孔 1 mm") == 2
    assert "900" in text and "1400" in text
    assert "| 会 |" in text and "| 不会 |" in text


def test_drying_shell_closes_the_energy_balance():
    """湿核推进期间，显热增量与蒸发潜热之和等于表面传入的热。"""
    particle = pyro.case1_properties()
    radius = 0.5e-3
    solved = pyro.integrate_drying(particle, radius)
    assert solved is not None
    balance = pyro.drying_energy_balance(particle, radius, solved)
    assert balance["entered_j"][-1] > 0.0
    assert balance["latent_j"][-1] > 0.0
    assert balance["relative"] < 0.02
    end = solved.sol(float(solved.t[-1]))
    assert float(end[0]) < 0.05 * radius
    core = pyro.temperature_on_radius(end, [0.0], particle, radius)
    assert float(core[0]) == pytest.approx(particle.evaporation_temperature_k)


def test_small_particle_times_are_finite_and_use_the_same_integral():
    particle = pyro.case1_properties()
    result = pyro.single_particle(0.1e-3, particle, radial_count=13)
    assert result["status"] in {"成立", "不成立"}
    assert result["T_w_s"] > 0.0
    assert result["T_v_s"] > 0.0
    assert math.isfinite(result["T_w_s"])
    assert math.isfinite(result["T_v_s"])
