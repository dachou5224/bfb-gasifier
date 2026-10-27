"""对数未知量与可分辨步长。求解手段：ADR 0042。阻尼与收敛门的测试仍在 test_driver.py。"""
import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.solver.gates import ResidualBlocks, inner_met
from hamel.solver.newton import (
    _SQRT_EPS, _advance, _column_scale, _jacobian, _log_mask, _opening_step, damped_newton, normalized_merit,
)


def _gas_blocks(size):
    return ResidualBlocks(np.arange(size), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)


def test_positive_variables_stay_positive_at_every_damping_step_and_zero_stays_linear():
    """正值流量在 y 上走任意 2⁻ⁱ；值为 0 的流量仍按 x 相加。"""
    gas = np.array([1.0e-6])
    mask = _log_mask(gas, ["gas_outflow"])
    assert mask[0]
    downward = np.array([-40.0])
    for index in range(31):
        trial = _advance(gas, downward, 2.0 ** (-index), mask)
        assert trial[0] > 0.0
    zero = np.array([0.0])
    assert not _log_mask(zero, ["gas_outflow"])[0]
    additive = _opening_step(0.0, 4.0, False)
    assert additive == (_SQRT_EPS * 4.0, _SQRT_EPS * 4.0)

    seen = []

    def feasible(trial):
        seen.append(np.array(trial, dtype=float, copy=True))
        return bool(np.isfinite(trial).all())

    def residual(trial):
        value = float(trial[0])
        return np.array([1.0 / value - 10.0, float(trial[1])])

    config = SolverConfig()
    blocks = ResidualBlocks(np.array([0, 1]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)
    result = damped_newton(
        np.array([1.0, 0.0]), residual, np.ones(2),
        lambda values: inner_met(values, blocks, config.residual_tolerance),
        feasible, lambda values: normalized_merit(values, blocks), config,
        kinds=["gas_outflow", "gas_outflow"],
    )
    assert result.status == "CONVERGED"
    assert result.x[0] == pytest.approx(0.1)
    assert result.x[1] == pytest.approx(0.0)
    assert seen
    assert all(item[0] > 0.0 for item in seen)


def test_trace_step_grows_until_the_column_is_nonzero_and_the_cap_is_recorded():
    """√ε·|x| 上不可分辨时放大步长。到上限仍不可分辨的列记入 unresolved。"""
    trace = 1.0e-20
    point = np.array([trace, trace])
    scales = np.ones(2)
    kinds = ["gas_outflow", "char"]

    def residual(trial):
        grown = float(trial[0])
        stuck = float(trial[1])
        return np.array([1.0 + 1.0e6 * grown, 1.0 + stuck])

    base = residual(point)
    unresolved = []
    steps = np.zeros(2)
    matrix = _jacobian(
        point, base, residual, scales, lambda trial: True, blocks=_gas_blocks(2),
        log_mask=_log_mask(point, kinds), names=["grown", "stuck"], unresolved=unresolved, steps_out=steps,
    )
    assert steps[0] > _SQRT_EPS
    assert matrix[0, 0] != 0.0
    assert [item["name"] for item in unresolved] == ["stuck"]
    assert unresolved[0]["x"] == trace
    assert unresolved[0]["step"] == pytest.approx(1.0)
    assert steps[1] == pytest.approx(1.0)


def test_log_difference_matches_the_analytic_derivative():
    """y 空间前向差分对上 ∂F/∂y = x·∂F/∂x。"""
    value = 1.0e-3
    point = np.array([value])

    def residual(trial):
        return np.array([float(trial[0]) ** 0.3])

    base = residual(point)
    mask = _log_mask(point, ["char"])
    matrix = _jacobian(point, base, residual, np.ones(1), lambda trial: True, blocks=_gas_blocks(1), log_mask=mask)
    analytic = 0.3 * value ** 0.3
    assert abs(matrix[0, 0] - analytic) / abs(analytic) < 1.0e-6
    assert _opening_step(value, 1.0, True)[0] == _SQRT_EPS


def test_column_scale_is_one_on_log_variables():
    state = np.array([2.0, 0.0, 5.0])
    scales = np.array([10.0, 10.0, 10.0])
    column = _column_scale(state, scales, np.array([True, False, False]))
    assert column.tolist() == [1.0, 10.0, 5.0]
    assert _column_scale(state, scales, None).tolist() == [2.0, 10.0, 5.0]


def test_log_iteration_reaches_the_same_root_as_the_linear_iteration():
    """已在 x 空间收敛的良态例，对数迭代给出同一解。"""
    target = np.array([1.0, 2.0, 3.0])
    config = SolverConfig()
    blocks = ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 10.0, 10.0, 10.0)
    start = np.array([0.2, 0.2, 0.2])

    def run(kinds):
        return damped_newton(
            start, lambda trial: np.asarray(trial, dtype=float) - target, np.ones(3),
            lambda values: inner_met(values, blocks, config.residual_tolerance),
            lambda trial: bool(np.all(trial > 0.0)),
            lambda values: normalized_merit(values, blocks), config, kinds=kinds,
        )

    plain = run(None)
    logged = run(["gas_outflow", "char", "ash"])
    assert plain.status == logged.status == "CONVERGED"
    assert plain.x == pytest.approx([1.0, 2.0, 3.0])
    assert logged.x == pytest.approx(plain.x, rel=1.0e-6, abs=1.0e-6)
    assert inner_met(plain.residual, blocks, config.residual_tolerance)
    assert inner_met(logged.residual, blocks, config.residual_tolerance)


def test_underflowing_log_step_is_rejected_and_the_accepted_point_stays_positive():
    """极大负步长使 exp 下溢为 0 时记为坐标域不可评价，接受的对数分量仍严格为正。"""
    from hamel.solver.newton import _damp

    x = np.array([1.0])
    mask = np.array([True])
    seen = []

    def evaluate(trial):
        seen.append(float(trial[0]))
        assert float(trial[0]) > 0.0
        return np.array([float(trial[0])])

    def feasible(trial):
        return bool(np.all(trial >= 0.0))

    downward = []
    accepted = _damp(
        x, np.array([-1000.0]), 1.0, evaluate, feasible, lambda values: float(np.linalg.norm(values)),
        30, None, downward, None, mask,
    )
    assert downward[0]["status"] == "reject"
    assert downward[0]["reason"] == "对数坐标还原越界"
    assert downward[0]["accepted"] is False
    assert accepted is not None
    assert accepted[0][0] > 0.0
    assert all(value > 0.0 for value in seen)
    upward = []
    overflow = _damp(
        x, np.array([1000.0]), 1.0, evaluate, feasible, lambda values: float(np.linalg.norm(values)),
        30, None, upward, None, mask,
    )
    assert upward[0]["reason"] == "对数坐标还原越界"
    assert overflow is None


def test_linear_trace_below_machine_epsilon_uses_a_larger_step():
    """线性痕量变量在 √ε·|x| 上不可分辨，放大后列非零，且不超过旧步长。"""
    value = 1.0e-20
    scale = 1.0
    point = np.array([value])

    def residual(trial):
        return np.array([1.0 + 1.0e6 * float(trial[0])])

    unresolved = []
    steps = np.zeros(1)
    matrix = _jacobian(
        point, residual(point), residual, np.array([scale]), lambda trial: True,
        blocks=_gas_blocks(1), unresolved=unresolved, steps_out=steps,
    )
    start, cap = _opening_step(value, scale, False)
    assert start == pytest.approx(_SQRT_EPS * value)
    assert steps[0] > start
    assert steps[0] <= cap
    assert matrix[0, 0] != 0.0
    assert unresolved == []
    assert cap == pytest.approx(_SQRT_EPS * scale)
