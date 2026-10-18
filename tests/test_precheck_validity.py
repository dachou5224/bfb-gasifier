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


def test_small_particle_times_are_finite_and_use_the_same_integral():
    particle = pyro.case1_properties()
    result = pyro.single_particle(0.1e-3, particle, radial_count=13)
    assert result["status"] in {"成立", "不成立"}
    assert result["T_w_s"] > 0.0
    assert result["T_v_s"] > 0.0
    assert math.isfinite(result["T_w_s"])
    assert math.isfinite(result["T_v_s"])
