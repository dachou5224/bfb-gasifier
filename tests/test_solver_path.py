"""ADR 0054 第 1、2 项：变量变换与三层试探规则。"""
from __future__ import annotations

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.solver.gates import ResidualBlocks, inner_met
from hamel.solver.newton import (
    IDENTITY, TrialRule, Transform, _grouped_jacobian, _jacobian, damped_newton, normalized_merit,
)

_POWER = Transform(np.zeros(0, dtype=int), np.array([0]), np.array([0.5]), np.array([1.0]))


def _rule(transform=IDENTITY, *, stock=(), gas=(), temperature=(), tau=1.0e-12, tau_n=None, names=("x",)):
    return TrialRule(
        transform, tuple(names), np.array(stock, dtype=int), np.array(gas, dtype=int),
        np.array(temperature, dtype=int), (), np.zeros(0), np.zeros(0), tau, tau_n,
    )


def _blocks(size):
    return ResidualBlocks(np.arange(size), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)


def _run(start, evaluate, rule, scales=None, config=None):
    config = SolverConfig() if config is None else config
    size = np.size(start)
    blocks = _blocks(size)
    return damped_newton(
        np.asarray(start, dtype=float), evaluate,
        np.ones(size) if scales is None else np.asarray(scales, dtype=float),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        rule.trial, lambda residual: normalized_merit(residual, blocks), config,
        transform=rule.transform,
    )


def test_sqrt_problem_takes_one_step_in_s_and_overshoots_in_x():
    """单调 √x 型标量问题：在 s = √x 上一步到位，在 x 上冲过头。"""
    evaluate = lambda x: np.array([np.sqrt(x[0]) - 1.0])
    solved = _run([np.sqrt(10.0)], evaluate, _rule(_POWER))
    assert solved.status == "CONVERGED" and solved.iterations == 1
    assert solved.x[0] == pytest.approx(1.0)
    # 在 x 上：x1 = 2√x0 − x0 = −3,68 < 0，试探点不可评价。
    plain = _run([10.0], evaluate, _rule())
    assert plain.iterations > 1


def test_negative_power_variable_is_rejected_without_evaluating_the_residual():
    """负 s 拒绝并减半，不评价残差；减半后落回 s ≥ 0 才接受。"""
    calls = {"n": 0}

    def evaluate(x):
        calls["n"] += 1
        return np.array([x[0] + 2.0])

    journal = {"iterations": [], "blocks": _blocks(1), "names": ["s"], "explain": _rule(_POWER).explain}
    config = SolverConfig(nu_max=1)
    blocks = _blocks(1)
    damped_newton(
        np.array([1.0]), evaluate, np.ones(1),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        _rule(_POWER).trial, lambda residual: normalized_merit(residual, blocks), config,
        journal=journal, transform=_POWER,
    )
    trials = journal["iterations"][0]["trials"]
    rejected = [item for item in trials if item["status"] != "ok"]
    assert rejected and all("幂变量 s < 0" in item["reason"] for item in rejected)
    # 在 s 上的牛顿步 Δs = −1,5 使 s = −0,5；减半一次后 s = 0,25 ≥ 0 才被评价并接受。
    # 起点 1 次 + 雅可比 1 列 1 次 + 接受点 1 次；被拒的负 s 试探没有评价残差。
    assert calls["n"] == 3
    assert trials[0]["accepted"] is False
    assert trials[1]["accepted"] is True and trials[1]["i"] == 1 and trials[1]["merit"] < 3.0


def test_zero_power_variable_is_evaluable_and_the_sign_is_not_folded():
    """s = 0 可评价；负 s 的逆映射带符号，不平方、不取绝对值。"""
    assert _POWER.to_x(np.array([0.0]))[0] == 0.0
    # n = 0,5：s = √Ṅ，逆映射 Ṅ = s²。
    assert _POWER.to_x(np.array([2.0]))[0] == 4.0
    # 负 s 的逆映射带符号，−s 与 +s 给出不同的 Ṅ，不平方折返。
    assert _POWER.to_x(np.array([-2.0]))[0] == -4.0
    assert _rule(_POWER).trial(np.array([0.0]))
    assert not _rule(_POWER).trial(np.array([-0.1]))


def test_zero_stock_with_a_negative_step_is_accepted_by_the_three_layer_rule():
    """零存量、牛顿分量为负：三层规则下可接受，旧的非负可行条件给 α = 0。"""
    evaluate = lambda x: np.array([x[0] + 1.0])
    rule = _rule(stock=[0])
    solved = _run([0.0], evaluate, rule)
    # 接受负存量，但收敛解的第 3 层检查拦下它。
    assert solved.status == "CONVERGED"
    detail = rule.final(solved.x)
    assert not detail["ok"] and "负存量" in detail["message"]
    old = _rule().trial
    strict = lambda trial: bool(np.all(np.asarray(trial, dtype=float) >= 0.0))
    blocks = _blocks(1)
    config = SolverConfig()
    stopped = damped_newton(
        np.array([0.0]), evaluate, np.ones(1),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        strict, lambda residual: normalized_merit(residual, blocks), config,
    )
    assert stopped.status == "DAMPING_FAILED" and stopped.iterations == 1
    assert old(np.array([0.0]))


def test_untransformed_gas_within_tau_n_does_not_block_the_full_step():
    """Ar 列在 −τ_n 内不阻挡整步；旧的可行条件在起点就判不可行。"""
    tau_n = 1.0e-12
    evaluate = lambda x: np.asarray(x, dtype=float) - 1.0
    start = np.array([0.0, -0.4 * tau_n])
    rule = _rule(gas=[1], tau_n=tau_n)
    solved = _run(start, evaluate, rule)
    assert solved.status == "CONVERGED" and solved.iterations == 1
    strict = lambda trial: bool(np.all(np.asarray(trial, dtype=float) >= 0.0))
    blocks = _blocks(2)
    config = SolverConfig()
    stopped = damped_newton(
        start, evaluate, np.ones(2),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        strict, lambda residual: normalized_merit(residual, blocks), config,
    )
    assert stopped.status == "NOT_CONVERGED" and stopped.iterations == 0


def test_grouped_jacobian_matches_columnwise_on_y():
    """分组雅可比在 y 上与逐列差分逐位相同。ADR 0046 第 5 项；ADR 0054 第 1 项不改变列的位置。"""
    transform = Transform(np.array([0, 3]), np.array([1, 4]), np.array([0.5, 0.3]), np.array([2.0, 1.5]))
    y = np.linspace(0.4, 1.2, 6)
    nodes = np.arange(6)
    layout = (nodes, np.zeros(6, dtype=int), np.zeros(6, dtype=bool))

    def evaluate(x):
        return 0.3 * x + 0.02 * x * x + 1.0e-3

    residual = evaluate(transform.to_x(y))
    grouped, stats = _grouped_jacobian(y, residual, evaluate, np.ones(6), lambda trial: True, layout, transform=transform)
    reference = _jacobian(y, residual, evaluate, np.ones(6), lambda trial: True, transform=transform)
    assert np.array_equal(grouped, reference)
    assert stats["fallbacks"] == 0


def test_scales_follow_the_transform():
    """对数列取 1，幂列取 ṅ_ref^n，其余按 ADR 0043。"""
    transform = Transform(np.array([0]), np.array([1]), np.array([0.5]), np.array([2.0]))
    assert transform.scaled(np.array([3.0, 4.0, 5.0])).tolist() == [1.0, 2.0, 5.0]
