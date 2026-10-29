"""自由空间参考轨迹的归宿与输运系数。

格内气速为常数，稳定点只在界面（ADR 0047）。例子都用两三个格与给定的常气速，
不是 Hamel 的鬼泡剖面。工况 1 的冒烟另标 slow。
"""

from __future__ import annotations

import numpy as np
import pytest
from dataclasses import replace

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError,
    account_freeboard,
    ejection_velocity_samples,
    ghost_bubble_velocity,
    stable_equilibria,
    terminal_settling_velocity,
    trace_until_fate,
)

_GAS = 1.2
_MU = 2.0e-5
_DENSITY = 1200.0
_PHI = 0.8
_DIAMETER = 2.0e-4
# 界面上下气速：下格 1,0、上格 0,2；v_t ≈ 0,849 落在两者之间。
_LOWER = 1.0
_UPPER = 0.2


def _config(**changes) -> SolverConfig:
    return replace(SolverConfig(), **changes)


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def _column(count: int, length: float, rho_g: float, mu: float):
    edges = np.linspace(0.0, length, count + 1)
    return edges, np.full(count, rho_g), np.full(count, mu)


def _cellwise(*speeds: float):
    def gas(cell, _height):
        return float(speeds[cell])

    return gas


def _account(edges, gas, rho, mu, diameter, rho_p, phi, bubble_velocity, samples, scale=1.0, **config):
    speeds, shares, direct = ejection_velocity_samples(bubble_velocity, samples)
    result = account_freeboard(
        edges, gas, rho, mu, diameter, rho_p, phi, speeds, shares, direct, _config(**config), scale
    )
    return result, direct


def test_interface_capture_keeps_the_last_entry_into_the_capture_region():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    equilibria = stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI)
    assert len(equilibria) == 1
    assert equilibria[0].height_m == pytest.approx(1.0)
    assert equilibria[0].cells == (1, 2)
    assert equilibria[0].half_width_m == pytest.approx(0.25)
    traced = trace_until_fate(
        0.2, 2.0, edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, _config(), equilibria
    )
    assert traced.fate == "capture"
    assert not traced.captured_at_start
    assert traced.entered_cell == 1
    assert traced.capture_time_s is not None and traced.capture_time_s > 0.0
    assert traced.velocity_reversals > 2
    assert traced.residence_s[0] == pytest.approx(traced.capture_time_s)
    assert traced.residence_s[1] == 0.0 and traced.residence_s[2] == 0.0
    assert traced.upward_hover[0] == 1.0


def test_local_generation_is_captured_at_start_and_has_no_free_state():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    result, direct = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8)
    assert result.capture_cells == (1, 2)
    assert result.return_share + result.elutriation_share + result.capture_share == pytest.approx(1.0)
    assert result.ejection_free + result.ejection_hover == pytest.approx(1.0 - direct)
    lower, middle, top = result.cells
    assert lower.free_state == "存在"
    assert lower.theta == 0 and middle.theta == 1 and top.theta == 1
    assert lower.a_t > 0.0 and lower.a_f == 0.0 and lower.b_f == 0.0
    assert middle.free_state == "自由状态不存在" and middle.a_f is None
    assert top.free_state == "自由状态不存在" and top.a_f is None
    assert sum(result.local_hover_inflow) + result.local_return + result.local_elutriation == 3.0
    assert result.local_hover_inflow[0] == 0.0
    assert sum(result.hover_inflow) == pytest.approx(result.capture_share)
    assert result.hover_inflow[1] == pytest.approx(result.capture_share)


def test_halving_delta_leaves_the_books_unchanged():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    full, _direct = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8)
    half, _unused = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8, scale=0.5)
    assert half.cells == full.cells
    assert half.hover_inflow == full.hover_inflow
    assert half.return_share == full.return_share
    assert half.elutriation_share == full.elutriation_share
    assert half.capture_share == full.capture_share
    assert half.ejection_free == full.ejection_free
    assert half.ejection_hover == full.ejection_hover


def test_doubling_samples_keeps_the_split_coefficients_within_1e_4():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    base, _direct = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 48)
    fine, _unused = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 96)
    assert fine.capture_share == pytest.approx(base.capture_share, rel=1e-4, abs=0.0)
    assert fine.hover_inflow == pytest.approx(base.hover_inflow, rel=1e-4, abs=0.0)
    for left, right in zip(base.cells, fine.cells):
        assert left.free_state == right.free_state
        for name in ("a_f", "a_t", "b_f", "b_t"):
            first, second = getattr(left, name), getattr(right, name)
            if first is None:
                assert second is None
                continue
            assert _relative(first, second) < 1.0e-4


def test_monotonic_profiles_are_not_captured():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    for speed, expect in ((0.2, "return"), (1.0 + 4.0, "elutriate")):
        gas = _cellwise(speed, speed, speed)
        assert stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI) == ()
        result, direct = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, speed, 8)
        assert result.capture_cells == ()
        assert result.capture_share == 0.0
        assert result.ejection_hover == 0.0
        assert all(cell.theta == 0 and cell.a_t == 0.0 and cell.b_t == 0.0 for cell in result.cells)
        if expect == "return":
            assert result.return_share == pytest.approx(1.0)
            assert result.elutriation_share == pytest.approx(0.0)
        else:
            assert result.return_share == pytest.approx(direct)
            assert result.elutriation_share == pytest.approx(1.0 - direct)


def test_neutral_cell_reports_and_stops():
    edges, rho, mu = _column(2, 1.0, _GAS, _MU)
    settling = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI)
    gas = _cellwise(settling, settling)
    with pytest.raises(FreeboardError, match="中性"):
        stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI)


def test_reaction_free_solution_matches_the_reference_books():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    speed = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI) + 4.0
    gas = _cellwise(speed, speed, speed)
    result, direct = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 48)
    assert result.capture_share == 0.0
    assert result.return_share == pytest.approx(direct)
    assert result.elutriation_share == pytest.approx(1.0 - direct)
    count = len(result.cells)
    matrix = np.zeros((count, count))
    right = np.zeros(count)
    right[0] = result.ejection_free
    for index, cell in enumerate(result.cells):
        assert cell.free_state == "存在"
        upward = cell.a_f + cell.a_t
        downward = cell.b_f + cell.b_t
        matrix[index, index] = upward + downward
        if index:
            matrix[index, index - 1] = -(result.cells[index - 1].a_f + result.cells[index - 1].a_t)
        if index + 1 < count:
            matrix[index, index + 1] = -(result.cells[index + 1].b_f + result.cells[index + 1].b_t)
    inventory = np.linalg.solve(matrix, right)
    for solved, cell in zip(inventory, result.cells):
        assert solved == pytest.approx(cell.m_hat)
    returned = (result.cells[0].b_f + result.cells[0].b_t) * inventory[0]
    elutriated = (result.cells[-1].a_f + result.cells[-1].a_t) * inventory[-1]
    assert returned == pytest.approx(0.0)
    assert elutriated == pytest.approx(result.elutriation_share)
    assert returned + direct + elutriated == pytest.approx(1.0)


def test_bed_sample_inside_the_capture_region_hovers_from_the_start():
    edges, rho, mu = _column(2, 1.0, _GAS, _MU)
    gas = _cellwise(_LOWER, _UPPER)
    result, direct = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8)
    assert result.capture_cells == (0, 1)
    assert result.ejection_hover == pytest.approx(result.capture_share)
    assert result.hover_inflow[0] == pytest.approx(result.capture_share)
    assert result.cells[0].theta == 1 and result.cells[1].theta == 1
    assert result.cells[0].free_state == "自由状态不存在"


def test_oversized_gas_region_has_no_interface_even_if_speeds_differ():
    """上下两格都高于 v_t（或都低于）时没有稳定界面。"""
    edges, rho, mu = _column(2, 1.0, _GAS, _MU)
    slow = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI)
    for lower, upper in ((slow + 0.4, slow + 0.2), (0.3 * slow, 0.1 * slow)):
        gas = _cellwise(lower, upper)
        assert stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI) == ()


def test_ghost_bubble_velocity_feeds_the_cell_mean():
    edges = np.array([0.0, 0.2875, 2.3])
    diameters = np.array([0.15, 0.25, 0.25])
    flows = np.array([0.04, 0.07])

    def gas(cell, height):
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, 1.6, 0.05)

    from hamel.freeboard import cell_mean_velocity

    mean = cell_mean_velocity(gas, 0, 0.0, 0.2875)
    samples = [gas(0, 0.2875 * (index + 0.5) / 200) for index in range(200)]
    assert mean == pytest.approx(float(np.mean(samples)), rel=1e-5)
