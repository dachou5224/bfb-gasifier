"""自由空间参考轨迹的归宿与输运系数。

单格结点、过阻尼、欠阻尼、界面和离开自由空间的例子都把 v_t 放在所测气速关系里。
气速是测试给定的线性场，不是 Hamel 的鬼泡剖面。工况 1 的冒烟另标。
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from dataclasses import replace

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError,
    account_freeboard,
    ejection_velocity_samples,
    ghost_bubble_velocity,
    integrate_freeboard,
    stable_equilibria,
    terminal_settling_velocity,
    trace_until_fate,
)
from hamel.freeboard.motion import drag_groups

# 过阻尼例：小颗粒、大 Stokes 数。α 小于 σ/4 时线性化是稳定结点。
_NODE_GAS = 0.5
_NODE_MU = 2.0e-5
_NODE_DENSITY = 2000.0
_NODE_PHI = 1.0
_NODE_DIAMETER = 1.5e-4
# 欠阻尼例：较大颗粒，气速坡度大于 σ/4。
_SPIRAL_GAS = 0.4
_SPIRAL_MU = 1.5e-5
_SPIRAL_DENSITY = 1500.0
_SPIRAL_PHI = 0.8
_SPIRAL_DIAMETER = 1.0e-3


def _config(**changes) -> SolverConfig:
    return replace(SolverConfig(), **changes)


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def _damping_gap(rho_g: float, mu: float, diameter: float, rho_p: float, phi: float, alpha: float) -> float:
    """λ² + σλ + σα = 0 的判别式。正是结点，负是焦点。"""
    settling = terminal_settling_velocity(rho_g, mu, diameter, rho_p, phi)
    stokes, quadratic, _weight = drag_groups(rho_g, mu, diameter, rho_p, phi)
    sigma = stokes + 2.0 * quadratic * settling
    return sigma * (sigma - 4.0 * alpha)


def _linear_gas(settling: float, alpha: float, height_star: float):
    def gas(_cell, height, settling=settling, alpha=alpha, height_star=height_star):
        return settling + alpha * (height_star - float(height))

    return gas


def _column(count: int, length: float, rho_g: float, mu: float):
    edges = np.linspace(0.0, length, count + 1)
    return edges, np.full(count, rho_g), np.full(count, mu)


def _account(edges, gas, rho, mu, diameter, rho_p, phi, bubble_velocity, samples, **config):
    speeds, shares, direct = ejection_velocity_samples(bubble_velocity, samples)
    result = account_freeboard(
        edges, gas, rho, mu, diameter, rho_p, phi, speeds, shares, direct, _config(**config)
    )
    return result, direct


def test_traced_flight_matches_the_existing_integrator():
    edges, rho, mu = _column(3, 1.2, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)

    def gas(_cell, height, settling=settling):
        return settling + 3.0 + 0.2 * math.sin(float(height))

    config = _config()
    equilibria = stable_equilibria(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    assert equilibria == ()
    direct = integrate_freeboard(0.0, 1.2, edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, config)
    traced = trace_until_fate(
        0.0, 1.2, edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, config, equilibria
    )
    assert traced.fate == direct[0]
    assert traced.time_s == direct[1]
    assert traced.height_m == direct[2]
    assert traced.velocity_m_s == direct[3]


def test_overdamped_node_has_no_reversal_and_no_free_state():
    """单格稳定结点。当地样本单调趋近，捕获时刻为 0，自由状态不存在。"""
    alpha = 1.5
    height_star = 0.6
    assert _damping_gap(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, alpha) > 0.0
    edges, rho, mu = _column(1, 1.0, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    gas = _linear_gas(settling, alpha, height_star)
    config = _config(freeboard_delta_h_m=0.02, freeboard_delta_u_m_s=0.01)
    equilibria = stable_equilibria(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    assert len(equilibria) == 1
    assert equilibria[0].cells == (0,)
    start = 0.3
    traced = trace_until_fate(
        start,
        alpha * (height_star - start),
        edges,
        gas,
        rho,
        mu,
        _NODE_DIAMETER,
        _NODE_DENSITY,
        _NODE_PHI,
        config,
        equilibria,
    )
    assert traced.fate == "capture"
    assert traced.velocity_reversals == 0
    assert traced.captured_at_start
    assert traced.capture_time_s == 0.0
    assert traced.entered_cell == 0
    assert traced.time_s > 0.0
    assert abs(traced.height_m - height_star) <= config.freeboard_delta_h_m
    result, direct = _account(
        edges,
        gas,
        rho,
        mu,
        _NODE_DIAMETER,
        _NODE_DENSITY,
        _NODE_PHI,
        0.05,
        4,
        freeboard_delta_h_m=0.02,
        freeboard_delta_u_m_s=0.01,
    )
    assert result.capture_cells == (0,)
    assert result.cells[0].free_state == "自由状态不存在"
    assert result.cells[0].m_hat is None
    assert result.cells[0].a_f is None
    assert result.cells[0].theta == 1
    assert result.hover_inflow[0] == pytest.approx(1.0 - direct)
    assert result.ejection_hover == pytest.approx(1.0 - direct)
    assert result.return_share + result.elutriation_share + result.capture_share == pytest.approx(1.0)


def test_halving_delta_leaves_the_node_books_unchanged():
    alpha = 1.5
    edges, rho, mu = _column(1, 1.0, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    gas = _linear_gas(settling, alpha, 0.6)
    kwargs = dict(edges=edges, gas=gas, rho=rho, mu=mu, diameter=_NODE_DIAMETER, rho_p=_NODE_DENSITY, phi=_NODE_PHI)
    full, _direct = _account(**kwargs, bubble_velocity=0.05, samples=4, freeboard_delta_h_m=0.02, freeboard_delta_u_m_s=0.01)
    half, _direct_half = _account(
        **kwargs, bubble_velocity=0.05, samples=4, freeboard_delta_h_m=0.01, freeboard_delta_u_m_s=0.005
    )
    assert half.cells[0].free_state == full.cells[0].free_state == "自由状态不存在"
    assert half.return_share == full.return_share
    assert half.elutriation_share == full.elutriation_share
    assert half.capture_share == full.capture_share
    assert half.ejection_free == full.ejection_free
    assert half.ejection_hover == full.ejection_hover


def test_starting_at_equilibrium_is_captured_at_once():
    edges, rho, mu = _column(1, 1.0, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    gas = _linear_gas(settling, 1.5, 0.6)
    equilibria = stable_equilibria(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    traced = trace_until_fate(
        0.6, 0.0, edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, _config(), equilibria
    )
    assert traced.fate == "capture"
    assert traced.time_s == 0.0
    assert traced.captured_at_start
    assert traced.capture_time_s == 0.0
    assert traced.velocity_reversals == 0
    assert traced.residence_s == (0.0,)


def test_underdamped_trajectory_reenters_and_the_last_entry_is_kept():
    alpha = 8.0
    assert _damping_gap(_SPIRAL_GAS, _SPIRAL_MU, _SPIRAL_DIAMETER, _SPIRAL_DENSITY, _SPIRAL_PHI, alpha) < 0.0
    edges = np.array([0.0, 1.0, 1.25, 2.5])
    rho = np.full(3, _SPIRAL_GAS)
    mu = np.full(3, _SPIRAL_MU)
    settling = terminal_settling_velocity(_SPIRAL_GAS, _SPIRAL_MU, _SPIRAL_DIAMETER, _SPIRAL_DENSITY, _SPIRAL_PHI)
    height_star = 1.125
    gas = _linear_gas(settling, alpha, height_star)
    config = _config(freeboard_delta_h_m=1.0e-2, freeboard_delta_u_m_s=5.0e-2, freeboard_time_limit_s=20.0)
    equilibria = stable_equilibria(edges, gas, rho, mu, _SPIRAL_DIAMETER, _SPIRAL_DENSITY, _SPIRAL_PHI)
    assert equilibria[0].cells == (1,)
    traced = trace_until_fate(
        height_star, 1.5, edges, gas, rho, mu, _SPIRAL_DIAMETER, _SPIRAL_DENSITY, _SPIRAL_PHI, config, equilibria
    )
    assert traced.fate == "capture"
    assert traced.velocity_reversals >= 2
    assert traced.capture_time_s > 0.0
    assert not traced.captured_at_start
    assert traced.entered_cell == 1
    assert traced.upward_free[1] == 1.0
    assert traced.downward_hover[2] == 1.0
    assert traced.residence_s[1] > 0.0
    assert traced.residence_s[2] > 0.0
    halved = trace_until_fate(
        height_star,
        1.5,
        edges,
        gas,
        rho,
        mu,
        _SPIRAL_DIAMETER,
        _SPIRAL_DENSITY,
        _SPIRAL_PHI,
        _config(freeboard_delta_h_m=5.0e-3, freeboard_delta_u_m_s=2.5e-2, freeboard_time_limit_s=20.0),
        equilibria,
    )
    assert halved.capture_time_s == traced.capture_time_s
    assert halved.entered_cell == traced.entered_cell
    assert halved.residence_s == traced.residence_s
    assert halved.upward_free == traced.upward_free
    assert halved.downward_hover == traced.downward_hover


def test_interface_root_sends_each_side_into_the_cell_it_enters():
    alpha = 0.8
    assert _damping_gap(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, alpha) > 0.0
    edges, rho, mu = _column(4, 4.0, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    gas = _linear_gas(settling, alpha, 2.0)
    config = _config(freeboard_delta_h_m=0.05, freeboard_delta_u_m_s=0.02)
    equilibria = stable_equilibria(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    assert len(equilibria) == 1
    assert equilibria[0].height_m == pytest.approx(2.0)
    assert equilibria[0].cells == (1, 2)
    from_below = trace_until_fate(0.2, 0.3, edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, config, equilibria)
    from_above = trace_until_fate(
        3.5, alpha * (2.0 - 3.5), edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, config, equilibria
    )
    assert from_below.fate == "capture" and from_below.entered_cell == 1
    assert from_above.fate == "capture" and from_above.entered_cell == 2
    result, direct = _account(
        edges,
        gas,
        rho,
        mu,
        _NODE_DIAMETER,
        _NODE_DENSITY,
        _NODE_PHI,
        0.15,
        8,
        freeboard_delta_h_m=0.05,
        freeboard_delta_u_m_s=0.02,
    )
    assert result.capture_cells == (1, 2)
    assert result.return_share + result.elutriation_share + result.capture_share == pytest.approx(1.0)
    assert result.ejection_free + result.ejection_hover == pytest.approx(1.0 - direct)
    assert result.hover_inflow[1] == pytest.approx(result.capture_share)
    assert result.hover_inflow[2] == 0.0
    assert result.local_return + result.local_elutriation + sum(result.local_hover_inflow) == 4.0
    assert result.local_hover_inflow[1] > 0.0
    assert result.local_hover_inflow[2] > 0.0
    lower = result.cells[0]
    assert lower.free_state == "存在"
    assert lower.a_t > 0.0
    assert lower.a_f == 0.0
    assert lower.a_t * lower.m_hat == pytest.approx(result.hover_inflow[1])
    assert result.cells[1].free_state == "自由状态不存在"
    assert result.cells[1].theta == 1
    assert result.cells[2].free_state == "自由状态不存在"
    assert result.cells[2].theta == 1
    upper = result.cells[3]
    assert upper.free_state == "存在"
    assert upper.b_t * upper.m_hat == pytest.approx(1.0)
    assert upper.theta == 0


def test_monotonic_or_out_of_range_profiles_are_not_captured():
    edges, rho, mu = _column(2, 2.0, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)

    def increasing(_cell, height, settling=settling):
        return settling - 0.2 + 0.8 * float(height)

    def slower(_cell, height, settling=settling):
        return (settling - 1.5) - 0.3 * float(height)

    def faster(_cell, height, settling=settling):
        return (settling + 3.0) - 0.2 * float(height)

    for gas, expect_return in ((increasing, True), (slower, True), (faster, False)):
        assert stable_equilibria(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI) == ()
        result, direct = _account(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, 0.3, 8)
        assert result.capture_cells == ()
        assert result.capture_share == 0.0
        assert result.ejection_hover == 0.0
        assert all(cell.a_t == 0.0 and cell.b_t == 0.0 for cell in result.cells)
        assert all(cell.theta == 0 for cell in result.cells)
        if expect_return:
            assert result.elutriation_share == 0.0
            assert result.return_share == pytest.approx(1.0)
        else:
            assert result.return_share == pytest.approx(direct)
            assert result.elutriation_share == pytest.approx(1.0 - direct)


def test_reaction_free_solution_matches_the_reference_books():
    edges, rho, mu = _column(3, 1.5, _NODE_GAS, _NODE_MU)
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)

    def gas(_cell, _height, settling=settling):
        return settling + 4.0

    result, direct = _account(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, 1.0, 48)
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


def test_doubling_samples_and_halving_the_step_stays_within_1e_4():
    node_settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    edges, rho, mu = _column(3, 1.5, _NODE_GAS, _NODE_MU)

    def carried(_cell, _height, settling=node_settling):
        return settling + 4.0

    base, _direct = _account(edges, carried, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, 1.0, 48)
    refined, _direct_refined = _account(
        edges,
        carried,
        rho,
        mu,
        _NODE_DIAMETER,
        _NODE_DENSITY,
        _NODE_PHI,
        1.0,
        96,
        freeboard_max_step_s=1.0e-4,
    )
    _assert_books_close(base, refined, 1.0e-4)

    edges_h, rho_h, mu_h = _column(2, 2.0, _NODE_GAS, _NODE_MU)
    hovering = _linear_gas(node_settling, 0.8, 1.2)
    hover, _direct_hover = _account(
        edges_h,
        hovering,
        rho_h,
        mu_h,
        _NODE_DIAMETER,
        _NODE_DENSITY,
        _NODE_PHI,
        0.15,
        48,
        freeboard_delta_h_m=0.08,
        freeboard_delta_u_m_s=0.01,
    )
    hover_refined, _unused = _account(
        edges_h,
        hovering,
        rho_h,
        mu_h,
        _NODE_DIAMETER,
        _NODE_DENSITY,
        _NODE_PHI,
        0.15,
        96,
        freeboard_delta_h_m=0.08,
        freeboard_delta_u_m_s=0.01,
        freeboard_max_step_s=1.0e-4,
    )
    assert hover.cells[0].a_t > 0.0
    assert hover.cells[1].free_state == "自由状态不存在"
    _assert_books_close(hover, hover_refined, 1.0e-4)


def test_halving_delta_leaves_split_coefficients_unchanged():
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    edges, rho, mu = _column(2, 2.0, _NODE_GAS, _NODE_MU)
    gas = _linear_gas(settling, 0.8, 1.2)
    common = dict(
        edges=edges,
        gas=gas,
        rho=rho,
        mu=mu,
        diameter=_NODE_DIAMETER,
        rho_p=_NODE_DENSITY,
        phi=_NODE_PHI,
        bubble_velocity=0.15,
        samples=12,
    )
    full, _direct = _account(**common, freeboard_delta_h_m=0.08, freeboard_delta_u_m_s=0.01)
    half, _direct_half = _account(**common, freeboard_delta_h_m=0.04, freeboard_delta_u_m_s=0.005)
    assert full.cells[0].a_t > 0.0
    assert half.cells[0] == full.cells[0]
    assert half.return_share == full.return_share
    assert half.elutriation_share == full.elutriation_share
    assert half.capture_share == full.capture_share
    assert half.ejection_free == full.ejection_free
    assert half.ejection_hover == full.ejection_hover
    assert half.hover_inflow == full.hover_inflow


def test_oversized_delta_is_rejected_and_the_time_limit_is_only_an_error():
    settling = terminal_settling_velocity(_NODE_GAS, _NODE_MU, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    edges, rho, mu = _column(1, 1.0, _NODE_GAS, _NODE_MU)
    gas = _linear_gas(settling, 1.5, 0.6)
    with pytest.raises(FreeboardError, match="一半"):
        _account(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI, 0.05, 2, freeboard_delta_h_m=0.5)
    equilibria = stable_equilibria(edges, gas, rho, mu, _NODE_DIAMETER, _NODE_DENSITY, _NODE_PHI)
    with pytest.raises(FreeboardError, match="跨过邻域"):
        trace_until_fate(
            0.0,
            0.2,
            edges,
            gas,
            rho,
            mu,
            _NODE_DIAMETER,
            _NODE_DENSITY,
            _NODE_PHI,
            _config(freeboard_delta_h_m=1.0e-3, freeboard_delta_u_m_s=100.0),
            equilibria,
        )
    with pytest.raises(FreeboardError, match="归宿仍未确定"):
        trace_until_fate(
            0.0,
            0.2,
            edges,
            gas,
            rho,
            mu,
            _NODE_DIAMETER,
            _NODE_DENSITY,
            _NODE_PHI,
            _config(freeboard_time_limit_s=1.0e-3, freeboard_delta_h_m=0.02, freeboard_delta_u_m_s=0.01),
            equilibria,
        )


def test_adjacent_stable_roots_keep_the_final_capture_region():
    """两个稳定根的捕获格相邻。穿过第一根再进入第二根时，按第二根自己的捕获区记账。"""
    edges = np.array([0.0, 1.0, 2.0])
    rho = np.ones(2)
    mu = np.full(2, 2.0e-5)
    diameter, density, phi = 3.0e-4, 1000.0, 0.8
    settling = terminal_settling_velocity(1.0, 2.0e-5, diameter, density, phi)

    def gas(_cell, height, settling=settling):
        h = float(height)
        return settling - 2.0 * (h - 0.5) * (h - 1.0) * (h - 1.5)

    config = _config(
        freeboard_max_step_s=1.0e-3,
        freeboard_delta_h_m=0.01,
        freeboard_delta_u_m_s=0.01,
        freeboard_time_limit_s=30.0,
    )
    equilibria = stable_equilibria(edges, gas, rho, mu, diameter, density, phi)
    assert [item.cells for item in equilibria] == [(0,), (1,)]
    assert equilibria[0].height_m == pytest.approx(0.5)
    assert equilibria[1].height_m == pytest.approx(1.5)
    traced = trace_until_fate(0.75, 5.0, edges, gas, rho, mu, diameter, density, phi, config, equilibria)
    assert traced.fate == "capture"
    assert traced.entered_cell == 1
    assert traced.capture_time_s > 0.0
    assert not traced.captured_at_start
    assert traced.residence_s[0] > 0.0
    assert traced.residence_s[1] == 0.0
    assert traced.upward_hover[0] == 1.0
    assert traced.upward_free[0] == 0.0
    assert abs(traced.height_m - 1.5) <= config.freeboard_delta_h_m
    result = account_freeboard(
        edges, gas, rho, mu, diameter, density, phi, np.array([20.0]), np.array([1.0]), 0.0, config
    )
    assert result.hover_inflow == pytest.approx((0.0, 1.0))
    lower, upper = result.cells
    assert lower.theta == 1 and upper.theta == 1
    assert lower.free_state == "存在"
    assert lower.a_f == 0.0 and lower.b_f == 0.0 and lower.b_t == 0.0
    assert lower.a_t > 0.0
    assert lower.a_t * lower.m_hat == pytest.approx(1.0)
    assert lower.m_hat > 0.0
    assert upper.free_state == "自由状态不存在"
    assert upper.a_f is None


def test_neighborhood_velocity_can_still_leave_before_the_energy_barrier():
    """贴着格面、恢复力很弱时，邻域里的速度仍可能靠惯性离开捕获区。

    默认 80 s 不够确认，只能报错，不能记成起步即捕获。时间够长时保留这次离开，
    δ 减半不删掉这段自由停留。
    """
    edges = np.array([0.0, 1.0, 2.0])
    rho = np.ones(2)
    mu = np.full(2, 2.0e-5)
    diameter, density, phi = 0.003, 1000.0, 0.8
    settling = terminal_settling_velocity(1.0, 2.0e-5, diameter, density, phi)
    height_star = 0.9975
    gas = _linear_gas(settling, 0.01, height_star)
    equilibria = stable_equilibria(edges, gas, rho, mu, diameter, density, phi)
    assert len(equilibria) == 1
    assert equilibria[0].cells == (0,)
    assert equilibria[0].half_width_m == pytest.approx(0.00125)
    with pytest.raises(FreeboardError, match="归宿仍未确定"):
        trace_until_fate(
            height_star, 0.0099, edges, gas, rho, mu, diameter, density, phi, SolverConfig(), equilibria
        )
    config = _config(freeboard_time_limit_s=200.0)
    traced = trace_until_fate(
        height_star, 0.0099, edges, gas, rho, mu, diameter, density, phi, config, equilibria
    )
    assert traced.fate == "capture"
    assert traced.entered_cell == 0
    assert not traced.captured_at_start
    assert traced.capture_time_s > 0.0
    assert traced.residence_s[0] > 0.0
    assert traced.residence_s[1] > 0.0
    assert traced.upward_free[0] == 1.0
    assert traced.downward_hover[1] == 1.0
    halved = trace_until_fate(
        height_star,
        0.0099,
        edges,
        gas,
        rho,
        mu,
        diameter,
        density,
        phi,
        _config(freeboard_time_limit_s=400.0, freeboard_delta_h_m=5.0e-4, freeboard_delta_u_m_s=5.0e-3),
        equilibria,
    )
    assert halved.capture_time_s == traced.capture_time_s
    assert halved.residence_s == traced.residence_s
    assert halved.upward_free == traced.upward_free
    assert halved.downward_hover == traced.downward_hover
    assert halved.residence_s[0] > 0.0


def _assert_books_close(left, right, tolerance: float) -> None:
    assert _relative(left.return_share, right.return_share) < tolerance
    assert _relative(left.elutriation_share, right.elutriation_share) < tolerance
    assert _relative(left.capture_share, right.capture_share) < tolerance
    assert _relative(left.ejection_free, right.ejection_free) < tolerance
    assert _relative(left.ejection_hover, right.ejection_hover) < tolerance
    for former, latter in zip(left.cells, right.cells):
        assert former.free_state == latter.free_state
        assert former.theta == latter.theta
        if former.free_state == "自由状态不存在":
            assert latter.a_f is None
            continue
        for name in ("a_f", "a_t", "b_f", "b_t"):
            assert _relative(getattr(former, name), getattr(latter, name)) < tolerance


@pytest.mark.slow
def test_case1_smoke_agrees_with_t27_scenarios_b_and_c():
    """工况 1 全部粒级、情景 B 和 C 的积分，约 4 分钟。

    气体量用 T27 的临时估算，床料取 ADR 0024 的 ε_mf=0.46、φ_S=0.8。
    T27：这两个情景下 d_k = 0,3466 mm 在九个角点都会悬停。这里只跑这一个床料点。
    按 CHARTER §7 标为 slow。默认 pytest 不跑；单独跑用 pytest -m slow。
    """
    from tools.precheck_hover import (
        BED_TEMPERATURE_K,
        FREEBOARD_HEIGHT_M,
        GAS_CONSTANT,
        MODEL_CLASSES_MM,
        PRESSURE_PA,
        bed_representative_diameter_m,
        case1_gas_scenario_b,
        case1_gas_scenario_c,
        surface_state,
    )

    edges = np.linspace(0.0, FREEBOARD_HEIGHT_M, 9)
    diameters = np.array([0.15, *([0.25] * 8)])
    scenarios = (
        ("B", case1_gas_scenario_b()),
        ("C", case1_gas_scenario_c()),
    )
    lines = []
    for label, gas in scenarios:
        state = surface_state(0.46, 0.8, bed_representative_diameter_m(), gas)
        primary = gas["primary_mol_s"] * GAS_CONSTANT * BED_TEMPERATURE_K / PRESSURE_PA
        total = (gas["primary_mol_s"] + gas["secondary_mol_s"]) * GAS_CONSTANT * BED_TEMPERATURE_K / PRESSURE_PA
        flows = np.array([primary, *([total] * 7)])
        rho = np.full(8, gas["rho_g"])
        mu = np.full(8, gas["mu_g"])

        def velocity(cell, height, flows=flows, edges=edges, diameters=diameters, state=state):
            return ghost_bubble_velocity(cell, height, edges, diameters, flows, state["u_b"], state["d_b"])

        lines.append(f"情景 {label}  u_b={state['u_b']:.4f} m/s  d_b={state['d_b']:.4f} m  ρ_P={state['rho_p']:.1f}")
        for diameter_mm, _fraction, band in MODEL_CLASSES_MM:
            speeds, shares, direct = ejection_velocity_samples(state["u_b"], SolverConfig().freeboard_ejection_samples)
            result = account_freeboard(
                edges,
                velocity,
                rho,
                mu,
                diameter_mm * 1.0e-3,
                state["rho_p"],
                0.8,
                speeds,
                shares,
                direct,
                SolverConfig(),
            )
            lines.append(
                f"  {band} mm  d={diameter_mm:.4f}  回落={result.return_share:.4f}"
                f"  带出={result.elutriation_share:.4f}  捕获={result.capture_share:.4f}"
                f"  h*={','.join(f'{item:.3f}' for item in result.stable_heights_m) or '无'}"
                f"  捕获区={result.capture_cells or '无'}"
            )
            for index, cell in enumerate(result.cells):
                if cell.free_state == "自由状态不存在":
                    lines.append(f"    格 {index}  自由状态不存在  θ={cell.theta}")
                else:
                    lines.append(
                        f"    格 {index}  a^f={cell.a_f:.4g}  a^t={cell.a_t:.4g}"
                        f"  b^f={cell.b_f:.4g}  b^t={cell.b_t:.4g}  θ={cell.theta}"
                    )
            if abs(diameter_mm - 0.3466) < 1.0e-6:
                assert result.capture_cells
                assert result.capture_share > 0.0
            if abs(diameter_mm - 0.0506) < 1.0e-6:
                assert result.capture_share == 0.0
                assert result.elutriation_share > result.return_share
            if abs(diameter_mm - 2.4885) < 1.0e-6:
                assert result.capture_share == 0.0
                assert result.return_share > result.elutriation_share
    print("\n".join(lines))
