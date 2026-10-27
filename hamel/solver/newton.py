"""内层阻尼牛顿。残差仍在物理未知量 x 上评价。不显式求逆。
正值气体流量、焦炭与灰存量以 y = ln x 迭代，试探点为 y + Δy/2ⁱ。
求牛顿步时解 (R·J·C)·z = −R·F，Δy = C·z。求解手段：ADR 0042。
Source: Hamel Eq.2.9, S. 20
解释性选择：ADR 0011 第 5 条。求解手段：ADR 0002；ADR 0017 第 6–8 项；ADR 0022；ADR 0042。
journal 只追加迭代和试探记录，不参与接受判定。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from hamel.kinetics.common import KineticsInputError
from hamel.properties.data import PropertyError
from hamel.solver.frozen import FreeboardDomainError, TrialNotEvaluable

_EPS = float(np.finfo(float).eps)
_SQRT_EPS = math.sqrt(_EPS)
_LOG_KINDS = ("gas_outflow", "char", "ash")
_GROWTH = 100.0
_GROWTH_LIMIT = 4
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


def damped_newton(x0, evaluate, scales, converged, feasible, merit_of, config, deadline=None, journal=None, blocks=None, kinds=None) -> NewtonResult:
    """Wirsum 第 4、5 步。merit_of 是归一残差的欧氏范数。
    雅可比每次在当前点用前向差分重算。不按块合并列；合并不改变各列的值。
    正值流量与存量在 y = ln x 上差分。求解手段：ADR 0042。
    """
    x = np.asarray(x0, dtype=float).copy()
    scales = np.asarray(scales, dtype=float)
    residual, status, _reason = _at(x, evaluate, feasible)
    if status == "domain":
        return NewtonResult("DOMAIN_STOP", "初值落在自由空间有效域外", x, residual, 0, None)
    if status != "ok":
        return NewtonResult("NOT_CONVERGED", "初值不可行或残差不是有限数", x, residual, 0, None)
    merit = merit_of(residual)
    if converged(residual):
        return NewtonResult("CONVERGED", "", x, residual, 0, merit)
    limit = int(config.damping_index_limit)
    explain = None if journal is None else journal.get("explain")
    names = None if journal is None else journal.get("names")
    for iteration in range(1, int(config.nu_max) + 1):
        _clock(deadline)
        trials = [] if journal is not None else None
        matrix = None
        delta = None
        unresolved = []
        log_mask = _log_mask(x, kinds)
        try:
            matrix = _jacobian(
                x, residual, evaluate, scales, feasible, blocks=blocks, log_mask=log_mask,
                names=names, unresolved=unresolved,
            )
            delta = _newton_step(matrix, residual, x, scales, blocks, log_mask)
            accepted = _damp(x, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials, explain, log_mask)
        except FreeboardDomainError as exc:
            _note(journal, iteration, merit, residual, delta, trials, scales, unresolved)
            return NewtonResult("DOMAIN_STOP", str(exc), x, residual, iteration, merit)
        except np.linalg.LinAlgError:
            _note(journal, iteration, merit, residual, delta, trials, scales, unresolved)
            return NewtonResult("NOT_CONVERGED", "雅可比不能求解", x, residual, iteration, merit)
        _note(journal, iteration, merit, residual, delta, trials, scales, unresolved)
        if accepted is None:
            _dump(journal, x, residual, matrix, delta, scales, log_mask)
            return NewtonResult("DAMPING_FAILED", "没有使归一残差严格下降的可行步", x, residual, iteration, merit)
        x, residual, merit = accepted
        if converged(residual):
            return NewtonResult("CONVERGED", "", x, residual, iteration, merit)
    return NewtonResult("NOT_CONVERGED", "内层超过 nu_max", x, residual, int(config.nu_max), merit)


def normalized_merit(residual, blocks) -> float:
    """欧氏范数取在各块除以 ADR 0012 参考流之后。求解手段：ADR 0017 第 7 项。"""
    return float(np.linalg.norm(_normalize(residual, blocks)))


def _difference_basis(value, scale) -> float:
    """x ≠ 0 时 |x|，x = 0 时 s。求解手段：ADR 0041。"""
    magnitude = abs(float(value))
    if magnitude != 0.0:
        return magnitude
    return float(scale)


def _log_mask(x, kinds):
    """当前值为正的气体流量、焦炭与灰取对数。恰为 0 的保持线性。求解手段：ADR 0042。"""
    mask = np.zeros(np.size(x), dtype=bool)
    if kinds is None:
        return mask
    for index, kind in enumerate(kinds):
        if kind in _LOG_KINDS and float(x[index]) > 0.0:
            mask[index] = True
    return mask


def _opening_step(value, scale, logarithmic):
    """线性起始按 ADR 0041，对数起始 √ε。上限是旧步长或 1。求解手段：ADR 0042。"""
    if logarithmic:
        return _SQRT_EPS, 1.0
    span = max(abs(float(value)), float(scale))
    return _SQRT_EPS * _difference_basis(value, scale), _SQRT_EPS * span


def _advance(x, delta, factor, log_mask):
    """沿 Δy 走 factor 倍。对数还原越界时不把 0 当试探点。求解手段：ADR 0042。"""
    point = np.array(x, dtype=float, copy=True)
    step = np.asarray(delta, dtype=float) * factor
    keep = ~log_mask
    point[keep] += step[keep]
    with np.errstate(over="ignore", under="ignore"):
        point[log_mask] = np.exp(np.log(x[log_mask]) + step[log_mask])
    return point


def _column_scale(x, scales, log_mask):
    """线性变量 max(|x|, s·[x=0])，对数变量 1。求解手段：ADR 0042。"""
    column = np.empty(x.size)
    for index in range(x.size):
        if log_mask is not None and log_mask[index]:
            column[index] = 1.0
        else:
            column[index] = _difference_basis(x[index], scales[index])
    return column


def _jacobian(x, residual, evaluate, scales, feasible, *, blocks=None, log_mask=None, names=None, unresolved=None, steps_out=None):
    """前向差分。不可分辨的列步长乘 10²，最多 4 次，且不超过上限。求解手段：ADR 0042。"""
    columns = np.empty((x.size, x.size))
    level = float(np.linalg.norm(_normalize(residual, blocks)))
    ceiling = _EPS * level
    for index in range(x.size):
        logarithmic = log_mask is not None and bool(log_mask[index])
        step, cap = _opening_step(x[index], scales[index], logarithmic)
        if not np.isfinite(step) or step <= 0.0 or not np.isfinite(cap) or cap <= 0.0:
            raise np.linalg.LinAlgError("差分步长不是正数")
        chosen = None
        for attempt in range(_GROWTH_LIMIT + 1):
            trial = x.copy()
            if logarithmic:
                trial[index] = math.exp(math.log(float(x[index])) + step)
            else:
                trial[index] += step
            values, status, _reason = _at(trial, evaluate, feasible)
            if status == "domain":
                raise FreeboardDomainError("差分点落在自由空间有效域外")
            if status != "ok":
                raise np.linalg.LinAlgError("差分点不可评价")
            chosen = values - residual
            if float(np.linalg.norm(_normalize(chosen, blocks))) > ceiling or step >= cap or attempt == _GROWTH_LIMIT:
                break
            step = min(step * _GROWTH, cap)
        columns[:, index] = chosen / step
        if steps_out is not None:
            steps_out[index] = step
        if float(np.linalg.norm(_normalize(chosen, blocks))) <= ceiling:
            _remember(unresolved, names, index, x[index], step)
    return columns


def _remember(unresolved, names, index, value, step) -> None:
    """到上限仍不可分辨的列。求解手段：ADR 0042。"""
    if unresolved is None:
        return
    label = str(index) if names is None or index >= len(names) else names[index]
    unresolved.append({"index": int(index), "name": label, "x": float(value), "step": float(step)})


def _newton_step(matrix, residual, x, scales, blocks, log_mask=None):
    """解 (R·J·C)·z = −R·F，Δy = C·z。不显式求逆。求解手段：ADR 0042。"""
    column = _column_scale(x, scales, log_mask)
    if not np.isfinite(column).all() or np.any(column <= 0.0):
        raise np.linalg.LinAlgError("列缩放不是正数")
    scaled = np.asarray(matrix, dtype=float) * column
    rhs = -np.asarray(residual, dtype=float)
    row = _row_weights(residual.size, blocks)
    if row is not None:
        scaled = scaled * row[:, None]
        rhs = rhs * row
    return column * np.linalg.solve(scaled, rhs)


def _row_weights(size, blocks):
    """R 的对角元，与归一残差同一组参考流。无分块时不缩放。"""
    if blocks is None:
        return None
    weights = np.ones(size, dtype=float)
    _fill_row(weights, blocks.gas, blocks.n_dot_ref)
    _fill_row(weights, blocks.solid, blocks.m_dot_ref)
    _fill_row(weights, blocks.energy, blocks.q_dot_ref)
    return weights


def _fill_row(weights, index, reference) -> None:
    if not np.isfinite(reference) or float(reference) <= 0.0:
        raise np.linalg.LinAlgError("行缩放的参考流不是正数")
    if np.size(index):
        weights[np.asarray(index, dtype=int)] = 1.0 / float(reference)


def _damp(x, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials=None, explain=None, log_mask=None):
    active = log_mask is not None and bool(np.any(log_mask))
    for index in range(limit + 1):
        _clock(deadline)
        if active:
            point = _advance(x, delta, 1.0 / (2.0 ** index), log_mask)
            outside = not np.isfinite(point[log_mask]).all() or bool(np.any(point[log_mask] <= 0.0))
        else:
            point = x + delta / (2.0 ** index)
            outside = False
        if outside:
            values, status, reason = None, "reject", "对数坐标还原越界"
        else:
            values, status, reason = _at(point, evaluate, feasible, explain)
        trial_merit = merit_of(values) if status == "ok" else None
        if trials is not None:
            trials.append({
                "i": index, "step": float(2.0 ** (-index)), "status": status, "reason": reason,
                "merit": None if trial_merit is None else float(trial_merit), "accepted": False,
            })
        if status == "domain":
            raise FreeboardDomainError("试探点落在自由空间有效域外")
        if status != "ok":
            continue
        if trial_merit < merit:
            if trials is not None:
                trials[-1]["accepted"] = True
            return point, values, trial_merit
    return None


def _at(x, evaluate, feasible, explain=None):
    if not np.isfinite(x).all():
        return None, "reject", "x 非有限"
    if not feasible(x):
        reason = "feasible 不满足"
        if explain is not None:
            try:
                reason = explain(x)
            except Exception:
                reason = "feasible 不满足"
        return None, "reject", reason
    try:
        values = np.asarray(evaluate(x), dtype=float)
    except FreeboardDomainError:
        return None, "domain", "试探点落在自由空间有效域外"
    except _SHRINK as exc:
        return None, "reject", f"{type(exc).__name__}: {exc}"
    if values.shape != x.shape:
        return None, "reject", "残差形状与未知量不一致"
    if not np.isfinite(values).all():
        return None, "reject", "残差非有限"
    return values, "ok", ""


def _normalize(residual, blocks):
    scaled = np.asarray(residual, dtype=float).copy()
    if blocks is None:
        return scaled
    _divide(scaled, blocks.gas, blocks.n_dot_ref)
    _divide(scaled, blocks.solid, blocks.m_dot_ref)
    _divide(scaled, blocks.energy, blocks.q_dot_ref)
    return scaled


def _note(journal, iteration, merit, residual, delta, trials, scales, unresolved=None) -> None:
    if journal is None:
        return
    blocks = journal.get("blocks")
    scaled = _normalize(residual, blocks)
    names = journal.get("names") or []
    labels = _block_labels(blocks)
    order = np.argsort(-np.abs(scaled))[:10]
    largest = [{
        "index": int(item), "block": labels.get(int(item), ""),
        "name": names[int(item)] if int(item) < len(names) else str(int(item)),
        "value": float(scaled[item]),
    } for item in order]
    step = {"norm": None, "max": None}
    if delta is not None:
        scaled_step = np.asarray(delta, dtype=float) / np.asarray(scales, dtype=float)
        index = int(np.argmax(np.abs(scaled_step)))
        step = {
            "norm": float(np.linalg.norm(scaled_step)),
            "max": {
                "index": index, "block": labels.get(index, ""),
                "name": names[index] if index < len(names) else str(index),
                "value": float(scaled_step[index]), "raw": float(delta[index]),
            },
        }
    journal.setdefault("iterations", []).append({
        "iteration": iteration, "merit": float(merit),
        "blocks": {
            "gas": _block_norm(scaled, None if blocks is None else blocks.gas),
            "solid": _block_norm(scaled, None if blocks is None else blocks.solid),
            "energy": _block_norm(scaled, None if blocks is None else blocks.energy),
            "bottom": journal.get("bottom_norm"),
        },
        "largest": largest, "step": step, "trials": [] if trials is None else trials,
        "unresolved": [] if unresolved is None else unresolved,
    })


def _block_labels(blocks) -> dict:
    labels = {}
    if blocks is None:
        return labels
    for title, index in (("气体", blocks.gas), ("固体", blocks.solid), ("能量", blocks.energy)):
        for item in np.asarray(index, dtype=int).tolist():
            labels[int(item)] = title
    return labels


def _block_norm(scaled, index) -> float:
    if index is None:
        return 0.0
    chosen = np.asarray(index, dtype=int)
    if chosen.size == 0:
        return 0.0
    return float(np.linalg.norm(scaled[chosen]))


def _dump(journal, x, residual, matrix, delta, scales, log_mask=None) -> None:
    if journal is None or not journal.get("failure_path") or matrix is None or delta is None:
        return
    path = Path(journal["failure_path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    blocks = journal.get("blocks")
    payload = {"x": np.asarray(x), "F": np.asarray(residual), "J": np.asarray(matrix), "dx": np.asarray(delta), "scales": np.asarray(scales)}
    if log_mask is not None:
        payload["log_mask"] = np.asarray(log_mask, dtype=bool)
    if blocks is not None:
        payload.update(
            gas=np.asarray(blocks.gas, dtype=int), solid=np.asarray(blocks.solid, dtype=int),
            energy=np.asarray(blocks.energy, dtype=int),
            n_dot_ref=np.array(blocks.n_dot_ref), m_dot_ref=np.array(blocks.m_dot_ref), q_dot_ref=np.array(blocks.q_dot_ref),
        )
    names = journal.get("names")
    if names is not None:
        payload["names"] = np.asarray(names)
    np.savez(path, **payload)
    journal["failure_npz"] = str(path)


def _divide(scaled, index, reference) -> None:
    if np.size(index):
        scaled[np.asarray(index, dtype=int)] /= reference


def _clock(deadline) -> None:
    if deadline is not None and deadline():
        raise SolveTimeout("运行超过墙钟上限")
