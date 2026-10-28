"""内层阻尼牛顿。残差在物理未知量 x 上评价。不显式求逆。
差分步长 √ε·max(|x|, s)。求牛顿步时解 (R·J·C)·z = −R·F，Δx = C·z。
求解手段：ADR 0043。
Source: Hamel Eq.2.9, S. 20
解释性选择：ADR 0011 第 5 条。求解手段：ADR 0002；ADR 0017 第 6–8 项；ADR 0022；ADR 0043。
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
_SHRINK = (TrialNotEvaluable, KineticsInputError, PropertyError)


class SolveTimeout(TimeoutError):
    """墙钟到点。状态记 TIMEOUT，不换策略。求解器参数：ADR 0037。"""


class JacobianMismatchError(RuntimeError):
    """分组雅可比与逐列差分不逐位相同。求解手段：ADR 0046 第 5 项。"""

    def __init__(self, rows, columns):
        self.rows, self.columns = rows, columns
        pairs = ", ".join(f"({int(row)}, {int(column)})" for row, column in list(zip(rows, columns))[:8])
        super().__init__(f"分组雅可比与逐列差分不一致，共 {len(rows)} 处：{pairs}")


@dataclass
class NewtonResult:
    """成功只看调用方传入的内层门。ν_max 与阻尼失败都不是成功。"""

    status: str
    reason: str
    x: np.ndarray
    residual: np.ndarray | None
    iterations: int
    merit: float | None


def damped_newton(x0, evaluate, scales, converged, feasible, merit_of, config, deadline=None, journal=None, blocks=None, jacobian_layout=None, verify_jacobian=False) -> NewtonResult:
    """Wirsum 第 4、5 步。merit_of 是归一残差的欧氏范数。
    雅可比每次在当前点用前向差分重算。不按块合并列；合并不改变各列的值。
    求解手段：ADR 0043。
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
    for iteration in range(1, int(config.nu_max) + 1):
        _clock(deadline)
        trials = [] if journal is not None else None
        matrix = None
        delta = None
        try:
            if jacobian_layout is None:
                matrix = _jacobian(x, residual, evaluate, scales, feasible, blocks=blocks)
            else:
                matrix, stats = _grouped_jacobian(x, residual, evaluate, scales, feasible, jacobian_layout)
                if journal is not None:
                    journal.setdefault("jacobian", []).append(stats)
                if verify_jacobian and iteration == 1:
                    found = np.nonzero(matrix != _jacobian(x, residual, evaluate, scales, feasible))
                    if found[0].size:
                        raise JacobianMismatchError(*found)
            delta = _newton_step(matrix, residual, x, scales, blocks)
            accepted = _damp(x, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials, explain)
        except JacobianMismatchError as exc:
            if journal is not None:
                journal["jacobian_mismatch"] = {"rows": [int(item) for item in exc.rows], "columns": [int(item) for item in exc.columns]}
            _note(journal, iteration, merit, residual, delta, trials, scales)
            return NewtonResult("JACOBIAN_MISMATCH", str(exc), x, residual, iteration, merit)
        except FreeboardDomainError as exc:
            _note(journal, iteration, merit, residual, delta, trials, scales)
            return NewtonResult("DOMAIN_STOP", str(exc), x, residual, iteration, merit)
        except np.linalg.LinAlgError:
            _note(journal, iteration, merit, residual, delta, trials, scales)
            return NewtonResult("NOT_CONVERGED", "雅可比不能求解", x, residual, iteration, merit)
        _note(journal, iteration, merit, residual, delta, trials, scales)
        if accepted is None:
            _dump(journal, x, residual, matrix, delta, scales)
            return NewtonResult("DAMPING_FAILED", "没有使归一残差严格下降的可行步", x, residual, iteration, merit)
        x, residual, merit = accepted
        if converged(residual):
            return NewtonResult("CONVERGED", "", x, residual, iteration, merit)
    return NewtonResult("NOT_CONVERGED", "内层超过 nu_max", x, residual, int(config.nu_max), merit)


def normalized_merit(residual, blocks) -> float:
    """欧氏范数取在各块除以 ADR 0012 参考流之后。求解手段：ADR 0017 第 7 项。"""
    return float(np.linalg.norm(_normalize(residual, blocks)))


def _difference_basis(value, scale) -> float:
    """max(|x|, s)。求解手段：ADR 0043。"""
    return max(abs(float(value)), float(scale))


def _column_scale(x, scales):
    """C 的对角元 max(|x|, s)。求解手段：ADR 0043。"""
    return np.array([_difference_basis(x[index], scales[index]) for index in range(np.size(x))])


def _jacobian(x, residual, evaluate, scales, feasible, *, blocks=None, columns=None):
    """前向差分。步长 √ε·max(|x|, s)。求解手段：ADR 0043。
    columns 给定时只算这些列，供 ADR 0046 的逐列退回与非相邻列使用。
    """
    chosen = range(x.size) if columns is None else columns
    matrix = np.zeros((x.size, x.size))
    for index in chosen:
        step = _SQRT_EPS * _difference_basis(x[index], scales[index])
        if not np.isfinite(step) or step <= 0.0:
            raise np.linalg.LinAlgError("差分步长不是正数")
        trial = x.copy()
        trial[index] += step
        values, status, _reason = _at(trial, evaluate, feasible)
        if status == "domain":
            raise FreeboardDomainError("差分点落在自由空间有效域外")
        if status != "ok":
            raise np.linalg.LinAlgError("差分点不可评价")
        matrix[:, index] = (values - residual) / step
    return matrix


def _grouped_jacobian(x, residual, evaluate, scales, feasible, layout):
    """节点序号模 3 相同、节点内序号相同的列同时扰动；每行只取邻域内被扰动的那一列。
    非相邻耦合列单独逐列差分；组内试探点不可行或超域时该组退回逐列。求解手段：ADR 0046。
    """
    nodes, within, separate = layout
    matrix = np.zeros((x.size, x.size))
    rows, groups = {}, {}
    for row in range(x.size):
        rows.setdefault(int(nodes[row]), []).append(row)
    for column in range(x.size):
        if not separate[column]:
            groups.setdefault((int(nodes[column]) % 3, int(within[column])), []).append(column)
    fallbacks = 0
    for columns in groups.values():
        steps = {column: _SQRT_EPS * _difference_basis(x[column], scales[column]) for column in columns}
        if not all(np.isfinite(step) and step > 0.0 for step in steps.values()):
            raise np.linalg.LinAlgError("差分步长不是正数")
        trial = x.copy()
        for column, step in steps.items():
            trial[column] += step
        values, status, _reason = _at(trial, evaluate, feasible)
        if status != "ok":
            fallbacks += 1
            matrix[:, columns] = _jacobian(x, residual, evaluate, scales, feasible, columns=columns)[:, columns]
            continue
        for column in columns:
            node = int(nodes[column])
            chosen = [row for item in (node - 1, node, node + 1) for row in rows.get(item, ())]
            matrix[chosen, column] = (values[chosen] - residual[chosen]) / steps[column]
    separate = np.flatnonzero(separate)
    matrix[:, separate] = _jacobian(x, residual, evaluate, scales, feasible, columns=separate)[:, separate]
    return matrix, {"fallbacks": fallbacks, "groups": len(groups), "separate": int(separate.size)}


def _newton_step(matrix, residual, x, scales, blocks):
    """解 (R·J·C)·z = −R·F，Δx = C·z。不显式求逆。求解手段：ADR 0043。"""
    column = _column_scale(x, scales)
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


def _damp(x, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials=None, explain=None):
    for index in range(limit + 1):
        _clock(deadline)
        point = x + delta / (2.0 ** index)
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


def _note(journal, iteration, merit, residual, delta, trials, scales) -> None:
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


def _dump(journal, x, residual, matrix, delta, scales) -> None:
    if journal is None or not journal.get("failure_path") or matrix is None or delta is None:
        return
    path = Path(journal["failure_path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    blocks = journal.get("blocks")
    payload = {"x": np.asarray(x), "F": np.asarray(residual), "J": np.asarray(matrix), "dx": np.asarray(delta), "scales": np.asarray(scales)}
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
