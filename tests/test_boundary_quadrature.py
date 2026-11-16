"""粒级内边界定位求积与极小负存量归零。解释性选择：ADR 0058。"""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.freeboard.accounting import _touch_max, boundary_weights
from hamel.solver.driver import _floor_and_recheck, _negative_stocks
from hamel.solver.gates import ResidualBlocks


def _partition(xs, signature, delta=1.0, bisections=6):
    nodes = list(xs)

    def probe(x):
        nodes.append(x)
        return len(nodes) - 1, signature(x)

    weights, counts = boundary_weights(xs, [signature(x) for x in xs], probe, delta, bisections)
    return weights, nodes, counts


def _mass_by_signature(weights, nodes, signature):
    masses = {}
    for index, value in weights.items():
        masses[signature(nodes[index])] = masses.get(signature(nodes[index]), 0.0) + value
    return masses


def test_equal_signatures_reproduce_the_equal_weight_sum():
    """签名处处相同时，每个取样点得 w、端点得 0。ADR 0058 第 1 项。"""
    xs = [0.0, 0.5, 1.5, 2.0]
    weights, nodes, counts = _partition(xs, lambda x: "A")
    assert counts == {"bisections": 0, "thirds": 0}
    assert weights == {1: 1.0, 2: 1.0}
    assert sum(weights.values()) == pytest.approx(len(xs) - 2)


def test_two_signatures_locate_the_boundary_within_the_error_bound():
    xs = [0.0, 1.0]
    weights, nodes, counts = _partition(xs, lambda x: "A" if x < 0.3 else "B")
    assert abs(weights[0] - 0.3) <= 2.0 ** -7
    assert sum(weights.values()) == pytest.approx(1.0)
    assert counts["thirds"] == 0


def test_three_signatures_and_opposite_location_errors():
    """A/C/B：两侧定位误差方向相反，C 段质量误差 ≤ w·2⁻ⁿ_b。ADR 0058 测试要求。"""
    def signature(x):
        return "A" if x < 0.093751 else ("C" if x < 0.406249 else "B")

    weights, nodes, _counts = _partition([0.0, 1.0], signature)
    masses = _mass_by_signature(weights, nodes, signature)
    assert abs(masses["A"] - 0.093751) <= 2.0 ** -7
    assert abs(masses["B"] - (1.0 - 0.406249)) <= 2.0 ** -7
    assert abs(masses["C"] - (0.406249 - 0.093751)) <= 2.0 ** -6
    assert masses["C"] == pytest.approx(0.296875)


def test_a_shifted_three_signature_block_changes_masses_within_one_step():
    def make(shift):
        return lambda x: "A" if x < 0.099999 + shift else ("C" if x < 0.4 + shift else "B")

    first_weights, first_nodes, _ = _partition([0.0, 1.0], make(0.0))
    second_weights, second_nodes, _ = _partition([0.0, 1.0], make(2.0e-6))
    first = _mass_by_signature(first_weights, first_nodes, make(0.0))
    second = _mass_by_signature(second_weights, second_nodes, make(2.0e-6))
    assert max(abs(first[key] - second[key]) for key in ("A", "B", "C")) <= 2.0 ** -6


def test_end_interval_a_c_b_assigns_each_segment_to_its_representative():
    # N_d = 1：两个端部区间，[0, 0.5] 内为 A/C/B。
    xs = [0.0, 0.5, 1.0]

    def signature(x):
        return "A" if x < 0.25 else ("C" if x < 0.5 else "B")

    weights, nodes, _counts = _partition(xs, signature)
    masses = _mass_by_signature(weights, nodes, signature)
    assert weights[0] == pytest.approx(0.2421875)  # 端点只得 A 段
    assert signature(nodes[0]) == "A"
    assert weights[1] == pytest.approx(0.5078125)  # B 段归最近取样点
    assert masses["C"] == pytest.approx(0.25)
    assert sum(weights.values()) == pytest.approx(1.0)


def test_boundary_moving_past_a_sample_point_changes_masses_within_one_step():
    xs = [0.0, 0.5, 1.0]
    left_sig = lambda x: "A" if x < 0.499999 else "B"
    right_sig = lambda x: "A" if x < 0.500001 else "B"
    left_weights, left_nodes, _ = _partition(xs, left_sig)
    right_weights, right_nodes, _ = _partition(xs, right_sig)
    left = _mass_by_signature(left_weights, left_nodes, left_sig)
    right = _mass_by_signature(right_weights, right_nodes, right_sig)
    assert max(abs(left[key] - right[key]) for key in ("A", "B")) <= 2.0 ** -6


def test_equal_signatures_match_the_equal_weight_reference():
    """签名处处相同时，逐格系数与等权求和的参考实现相同。ADR 0058 测试要求。"""
    from hamel.freeboard.accounting import account_freeboard, within_class_samples
    from hamel.freeboard.ejection import ejection_velocity_samples
    from hamel.freeboard.motion import terminal_settling_velocity

    edges = np.linspace(0.0, 2.0, 3)
    rho, mu = np.full(2, 1.2), np.full(2, 2.0e-5)
    speed = terminal_settling_velocity(1.2, 2.0e-5, 2.0e-4, 1200.0, 0.8) + 4.0
    gas = lambda _cell, _height: speed
    velocities, shares, direct = ejection_velocity_samples(speed, 16)
    samples = within_class_samples(2.0e-4, 4.0e-4, diameters=4, sphericities=2)
    config = SolverConfig()
    result = account_freeboard(
        edges, gas, rho, mu, 2.0e-4, 1200.0, 0.8, velocities, shares, direct, config, samples=samples
    )
    assert result.boundary[1] == 0
    residence, upward = np.zeros(2), np.zeros(2)
    from hamel.freeboard import trace_until_fate

    for diameter, phi, weight in samples:
        ledgers = [trace_until_fate(0.0, float(u), edges, gas, rho, mu, diameter, 1200.0, phi, config)
                   for u in velocities]
        for ledger, share in zip(ledgers, shares):
            factor = weight * float(share)
            residence += factor * np.asarray(ledger.residence_s)
            upward += factor * np.asarray(ledger.upward_free)
    for cell, book in enumerate(result.cells):
        assert book.a_f == pytest.approx(upward[cell] / residence[cell], rel=1.0e-12)
    assert result.elutriation_share == pytest.approx(float(np.sum(shares)))


def test_node_weight_is_the_maximum_over_trajectory_positions():
    touched = {}
    for value in (0.1, 0.3, 0.2):
        _touch_max(touched, 1, value)
    _touch_max(touched, 2, 0.05)
    assert touched[1] == 0.3
    assert touched[2] == 0.05


def test_invalid_solver_parameters_are_rejected():
    with pytest.raises(ValueError):
        SolverConfig(freeboard_boundary_bisections=-1)
    with pytest.raises(ValueError):
        SolverConfig(freeboard_boundary_bisections=1.5)
    with pytest.raises(ValueError):
        SolverConfig(freeboard_boundary_bisections=True)
    with pytest.raises(ValueError):
        SolverConfig(final_negative_floor_kg=-1.0)
    with pytest.raises(ValueError):
        SolverConfig(final_negative_floor_kg=math.inf)
    with pytest.raises(ValueError):
        SolverConfig(final_negative_floor_kg=True)


def _frozen(bed=1, free=1, classes=1):
    return SimpleNamespace(
        n_bed=bed, n_freeboard=free, n_class=classes,
        hover=np.zeros((free, classes), dtype=bool),
        free_absent=np.zeros((free, classes), dtype=bool),
    )


def _values(frozen, target, value):
    from hamel.solver.vector import unknown_slots

    vector = np.zeros(len(unknown_slots(frozen)))
    for slot in unknown_slots(frozen):
        if slot.name.endswith(target):
            vector[slot.index] = value
            return vector, slot
    raise AssertionError("找不到目标变量")


def _checked(residual=0.0, ash=None):
    return SimpleNamespace(
        residual=np.array([residual]), blocks=ResidualBlocks(np.array([0]), np.array([], int), np.array([], int), 1.0, 1.0, 1.0),
        ash_relative=ash,
    )


def test_negative_stock_in_the_floor_range_is_zeroed_and_rechecked():
    frozen = _frozen()
    vector, slot = _values(frozen, "char.bed.k1", -1.0e-40)
    calls = []

    def consistent(values, _round):
        calls.append(values.copy())
        return _checked()

    values, report, reason = _floor_and_recheck(vector, frozen, SolverConfig(), consistent, 3, {"solid": 1.0})
    assert reason is None and report["applied"] and report["count"] == 1
    assert report["components"][0][0] == slot.name
    assert values[slot.index] == 0.0 and len(calls) == 1
    assert report["ratios_after"]["solid"] <= 1.0e-6


def test_negative_stock_below_the_floor_is_not_zeroed():
    frozen = _frozen()
    vector, slot = _values(frozen, "char.bed.k1", -1.0e-20)
    values, report, _reason = _floor_and_recheck(vector, frozen, SolverConfig(), lambda *_: _checked(), 1, None)
    assert not report["applied"] and report["below_floor"]
    assert values[slot.index] == -1.0e-20


def test_post_floor_block_residual_above_the_gate_is_not_converged():
    frozen = _frozen()
    vector, _slot = _values(frozen, "char.bed.k1", -1.0e-40)
    _values_out, _report, reason = _floor_and_recheck(
        vector, frozen, SolverConfig(), lambda *_: _checked(residual=1.0), 1, None
    )
    assert reason == "置零后 2(b) 复核未通过"


def test_post_floor_ash_balance_above_the_gate_is_not_converged():
    frozen = _frozen()
    vector, _slot = _values(frozen, "char.bed.k1", -1.0e-40)
    _values_out, _report, reason = _floor_and_recheck(
        vector, frozen, SolverConfig(), lambda *_: _checked(ash=1.0), 1, None
    )
    assert reason == "置零后灰平衡复核未通过"


def test_no_negative_stock_is_unchanged():
    frozen = _frozen()
    vector, slot = _values(frozen, "char.bed.k1", 1.0e-3)
    values, report, reason = _floor_and_recheck(vector, frozen, SolverConfig(), lambda *_: _checked(), 1, None)
    assert report["applied"] is False and reason is None
    assert np.array_equal(values, vector) and values[slot.index] == 1.0e-3


def test_negative_stocks_only_returns_the_two_lists():
    frozen = _frozen()
    vector, _slot = _values(frozen, "char.bed.k1", -1.0e-40)
    _values_out, floored, below = _negative_stocks(vector, frozen, 1.0e-30)
    assert len(floored) == 1 and not below
