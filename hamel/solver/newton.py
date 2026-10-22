"""内层阻尼牛顿。解 J·Δx = −F，不显式求逆。
Source: Hamel Eq.2.9, S. 20
解释性选择：ADR 0011 第 5 条。求解手段：ADR 0002；ADR 0017 第 6–8 项；ADR 0022。
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from hamel.kinetics.common import KineticsInputError
from hamel.properties.data import PropertyError
from hamel.solver.frozen import FreeboardDomainError, TrialNotEvaluable

_SQRT_EPS = math.sqrt(float(np.finfo(float).eps))
_SHRINK = (TrialNotEvaluable, KineticsInputError, PropertyError)


class SolveTimeout(TimeoutError):
    """墙钟到点。状态记 TIMEOUT，不换策略。求解器参数：ADR 0037。"""


@dataclass
class NewtonResult:
    """成功只看调用方传入的内层门。ν_max 与阻尼失败都不是成功。"""

    status: str
    reason: str
    x: np.ndarray
    residual: np.ndarray | None
    iterations: int
    merit: float | None


def damped_newton(x0, evaluate, scales, converged, feasible, merit_of, config, deadline=None) -> NewtonResult:
    """Wirsum 第 4、5 步。merit_of 是归一残差的欧氏范数。
    雅可比每次在当前点用前向差分重算。不按块合并列；合并不改变各列的值。
    """
    x = np.asarray(x0, dtype=float).copy()
    scales = np.asarray(scales, dtype=float)
    residual, status = _at(x, evaluate, feasible)
    if status == "domain":
        return NewtonResult("DOMAIN_STOP", "初值落在自由空间有效域外", x, residual, 0, None)
    if status != "ok":
        return NewtonResult("NOT_CONVERGED", "初值不可行或残差不是有限数", x, residual, 0, None)
    merit = merit_of(residual)
    if converged(residual):
        return NewtonResult("CONVERGED", "", x, residual, 0, merit)
    limit = int(config.damping_index_limit)
    for iteration in range(1, int(config.nu_max) + 1):
        _clock(deadline)
        try:
            matrix = _jacobian(x, residual, evaluate, scales, feasible)
            delta = np.linalg.solve(matrix, -residual)
            accepted = _damp(x, delta, merit, evaluate, feasible, merit_of, limit, deadline)
        except FreeboardDomainError as exc:
            return NewtonResult("DOMAIN_STOP", str(exc), x, residual, iteration, merit)
        except np.linalg.LinAlgError:
            return NewtonResult("NOT_CONVERGED", "雅可比不能求解", x, residual, iteration, merit)
        if accepted is None:
            return NewtonResult("DAMPING_FAILED", "没有使归一残差严格下降的可行步", x, residual, iteration, merit)
        x, residual, merit = accepted
        if converged(residual):
            return NewtonResult("CONVERGED", "", x, residual, iteration, merit)
    return NewtonResult("NOT_CONVERGED", "内层超过 nu_max", x, residual, int(config.nu_max), merit)


def normalized_merit(residual, blocks) -> float:
    """欧氏范数取在各块除以 ADR 0012 参考流之后。求解手段：ADR 0017 第 7 项。"""
    scaled = np.asarray(residual, dtype=float).copy()
    _divide(scaled, blocks.gas, blocks.n_dot_ref)
    _divide(scaled, blocks.solid, blocks.m_dot_ref)
    _divide(scaled, blocks.energy, blocks.q_dot_ref)
    return float(np.linalg.norm(scaled))


def _jacobian(x, residual, evaluate, scales, feasible):
    columns = np.empty((x.size, x.size))
    for index in range(x.size):
        step = _SQRT_EPS * max(abs(float(x[index])), float(scales[index]))
        if not np.isfinite(step) or step == 0.0:
            raise np.linalg.LinAlgError("差分步长不是正数")
        trial = x.copy()
        trial[index] += step
        values, status = _at(trial, evaluate, feasible)
        if status == "domain":
            raise FreeboardDomainError("差分点落在自由空间有效域外")
        if status != "ok":
            raise np.linalg.LinAlgError("差分点不可评价")
        columns[:, index] = (values - residual) / step
    return columns


def _damp(x, delta, merit, evaluate, feasible, merit_of, limit, deadline):
    for index in range(limit + 1):
        _clock(deadline)
        values, status = _at(x + delta / (2.0 ** index), evaluate, feasible)
        if status == "domain":
            raise FreeboardDomainError("试探点落在自由空间有效域外")
        if status != "ok":
            continue
        trial_merit = merit_of(values)
        if trial_merit < merit:
            return x + delta / (2.0 ** index), values, trial_merit
    return None


def _at(x, evaluate, feasible):
    if not np.isfinite(x).all() or not feasible(x):
        return None, "reject"
    try:
        values = np.asarray(evaluate(x), dtype=float)
    except FreeboardDomainError:
        return None, "domain"
    except _SHRINK:
        return None, "reject"
    if values.shape != x.shape or not np.isfinite(values).all():
        return None, "reject"
    return values, "ok"


def _divide(scaled, index, reference) -> None:
    if np.size(index):
        scaled[np.asarray(index, dtype=int)] /= reference


def _clock(deadline) -> None:
    if deadline is not None and deadline():
        raise SolveTimeout("运行超过墙钟上限")
