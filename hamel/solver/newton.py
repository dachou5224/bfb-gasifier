"""内层阻尼牛顿。迭代在变换后的自变量 y 上进行，残差在 x = x(y) 上按原式评价。不显式求逆。
差分步长 √ε·max(|y|, s)。求牛顿步时解 (R·J·C)·z = −R·F，Δy = C·z。
求解手段：ADR 0043；ADR 0054 第 1、2 项。
Source: Hamel Eq.2.9, S. 20
解释性选择：ADR 0011 第 5 条。求解手段：ADR 0002；ADR 0017 第 6 项；ADR 0022；ADR 0046；ADR 0054。
journal 只追加迭代和试探记录，不参与接受判定。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from hamel.kinetics.common import KineticsInputError
from hamel.properties.data import PropertyError
from hamel.solver.frozen import FreeboardDomainError, TrialNotEvaluable

_EPS = float(np.finfo(float).eps)
_SQRT_EPS = math.sqrt(_EPS)
_SHRINK = (TrialNotEvaluable, KineticsInputError, PropertyError)
# 差分点上自由空间有效域被拒时先缩步再判不可评价。求解手段：ADR 0054 第 2 项第 1 层。
_DIFFERENCE_SHRINK = 5


class SolveTimeout(TimeoutError):
    """墙钟到点。状态记 TIMEOUT，不换策略。求解器参数：ADR 0037。"""


class JacobianMismatchError(RuntimeError):
    """分组雅可比与逐列差分不逐位相同。求解手段：ADR 0046 第 5 项。"""

    def __init__(self, rows, columns):
        self.rows, self.columns = rows, columns
        pairs = ", ".join(f"({int(row)}, {int(column)})" for row, column in list(zip(rows, columns))[:8])
        super().__init__(f"分组雅可比与逐列差分不一致，共 {len(rows)} 处：{pairs}")


@dataclass(frozen=True)
class Transform:
    """ADR 0054 第 1 项：y ↔ x 的一对映射。log 取 z = ln Ṅ，power 取 s = Ṅ^n；残差函数不知道变换。"""

    # Python ≥ 3.11 不接受不可哈希的 ndarray 作默认值，数组字段用 default_factory（PR #152 审查 P1）。
    log: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    power: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    order: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    power_scale: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))

    def to_x(self, y):
        """s < 0 带符号映为负 Ṅ：只让试探被拒，不平方或取绝对值折返为正。"""
        x = np.array(y, dtype=float)
        with np.errstate(over="ignore", invalid="ignore"):
            x[self.log] = np.exp(x[self.log])
            x[self.power] = np.sign(x[self.power]) * np.abs(x[self.power]) ** (1.0 / self.order)
        return x

    def from_x(self, x):
        y = np.array(x, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            y[self.log] = np.log(y[self.log])
            y[self.power] = np.sign(y[self.power]) * np.abs(y[self.power]) ** self.order
        return y

    def scaled(self, scales):
        """差分与列缩放的尺度 s_y：对数列 1，幂列 ṅ_ref^n，其余按 ADR 0043。"""
        out = np.array(scales, dtype=float)
        out[self.log], out[self.power] = 1.0, self.power_scale
        return out


IDENTITY = Transform()


def build_transform(slots, x, n_ref) -> Transform:
    """ADR 0054 第 1 项：O₂ 取 s = Ṅ^n；其余起点上 Ṅ > 0 的气体出流取 ln Ṅ；起点为 0 的不变换。"""
    values, log, power, order = np.asarray(x, dtype=float), [], [], []
    for slot in slots:
        if slot.kind != "gas_outflow":
            continue
        if slot.species == "O2":
            power.append(slot.index)
            order.append(0.5 if slot.phase == "suspension" else 0.3)
        elif values[slot.index] > 0.0:
            log.append(slot.index)
    power, order = np.array(power, dtype=int), np.array(order, dtype=float)
    return Transform(np.array(log, dtype=int), power, order, np.asarray(n_ref, dtype=float) ** order)


@dataclass(frozen=True)
class TrialRule:
    """ADR 0054 第 2 项：第 2 层（试探点）与第 3 层（收敛解）检查，内层与初值构造共用。

    τ = stock_roundoff_relative·M_tot，τ_n 为同一系数乘 ṅ_ref。
    """

    transform: Transform
    names: tuple
    stock: np.ndarray
    gas: np.ndarray
    temperature: np.ndarray
    bed_cells: tuple
    geometric: np.ndarray
    fuel: np.ndarray
    tau: float
    tau_n: float | None = None
    bottom_ash: bool = False

    def inert(self, values) -> np.ndarray:
        """各床格惰性物：几何存量 − 焦炭 − 未转化燃料。ADR 0009 第 3 条。"""
        char = [np.sum(np.asarray(values)[list(cell)]) for cell in self.bed_cells]
        return self.geometric - np.array(char) - self.fuel

    def trial(self, y) -> bool:
        return not self.explain(y)

    def explain(self, y) -> str:
        """第 2 层的第一条不满足项；全部满足时返回空串。"""
        values = np.asarray(y, dtype=float)
        if not np.isfinite(values).all():
            return "y 非有限"
        checks = (
            ("幂变量 s < 0", self.transform.power, values[self.transform.power] < 0.0),
            ("不变换的气体出流 < −τ_n", self.gas, False if self.tau_n is None else values[self.gas] < -self.tau_n),
            ("T ≤ 0", self.temperature, values[self.temperature] <= 0.0),
        )
        for label, index, bad in checks:
            if np.size(index) and np.any(bad):
                first = int(index[np.argmax(bad)])
                return f"{label}（{self.names[first]} = {values[first]:.3e}）"
        fraction = self.inert(values)[-1] / self.geometric[-1] if self.bottom_ash else 1.0
        return "" if fraction > 0.0 else f"最底床格灰分为 {fraction:.3e}，排灰总流率没有定义"

    def final(self, x, drain_ash: float | None = None) -> dict:
        """第 3 层：存量 ≥ −τ、各床格惰性物 ≥ 0、不变换的气体 ≥ −τ_n；初值构造另查排灰灰分量。"""
        values = np.asarray(x, dtype=float)
        inert = self.inert(values)
        cells = np.flatnonzero(inert < 0.0)
        within = {self.names[int(i)]: float(values[i]) for i in self.stock if -self.tau <= values[i] < 0.0}
        below = [float(value) for value in values[self.stock] if value < -self.tau]
        gas = [] if self.tau_n is None else [float(values[i]) for i in self.gas if values[i] < -self.tau_n]
        message = (
            f"收敛解出现 {len(below)} 项负存量（小于 −τ = {self.tau:.3e} kg）" if below else
            f"收敛解在床区格 {int(cells[0])}（衡算下标）的惰性物为负：{inert[cells[0]]:.3e} kg" if cells.size else
            f"收敛解有 {len(gas)} 项不变换的气体出流低于 −τ_n" if gas else
            f"收敛解的排灰灰分量为负：{drain_ash:.3e} kg/s" if self.bottom_ash and drain_ash is not None and drain_ash < 0.0 else
            "收敛解最底床格灰分为零，排灰总流率没有定义" if self.bottom_ash and inert[-1] / self.geometric[-1] <= 0.0 else ""
        )
        return {"ok": not message, "message": message, "stocks_within_tau": within, "stocks_below_tau": below[:8],
                "gas_below_tau_n": gas[:8], "inert_kg": inert.tolist(), "inert_negative_cells": cells.tolist()}


def build_rule(slots, frozen, transform, *, tau: float, tau_n: float | None = None, unknowns=None, bottom_ash: bool = False) -> TrialRule:
    """从变量名表与冻结量给出三层规则的下标。求解手段：ADR 0054 第 2 项。"""
    chosen = None if unknowns is None else set(int(item) for item in unknowns)
    bed_cells = [[] for _ in range(int(frozen.n_bed))]
    stock, gas, temperature = [], [], []
    transformed = set(transform.log.tolist()) | set(transform.power.tolist())
    for slot in slots:
        if chosen is not None and int(slot.index) not in chosen:
            continue
        if slot.kind in ("char", "ash"):
            stock.append(slot.index)
            if slot.zone == "bed":
                bed_cells[slot.balance_index].append(slot.index)
        elif slot.kind == "temperature":
            temperature.append(slot.index)
        elif slot.kind == "gas_outflow" and slot.index not in transformed:
            gas.append(slot.index)
    return TrialRule(
        transform, tuple(slot.name for slot in slots), np.array(stock, dtype=int), np.array(gas, dtype=int),
        np.array(temperature, dtype=int), tuple(tuple(cell) for cell in bed_cells),
        np.asarray(frozen.geometric_mass_kg, dtype=float), np.asarray(frozen.fuel_bed_kg.sum(axis=1), dtype=float),
        float(tau), None if tau_n is None else float(tau_n), bool(bottom_ash),
    )


@dataclass
class NewtonResult:
    """成功只看调用方传入的内层门。ν_max 与阻尼失败都不是成功。"""

    status: str
    reason: str
    x: np.ndarray
    residual: np.ndarray | None
    iterations: int
    merit: float | None


def damped_newton(x0, evaluate, scales, converged, feasible, merit_of, config, deadline=None, journal=None, blocks=None, jacobian_layout=None, verify_jacobian=False, transform=IDENTITY, final_check=None) -> NewtonResult:
    """Wirsum 第 4、5 步。merit_of 是归一残差的欧氏范数。
    雅可比每次在当前点用前向差分重算。给定布局时按 ADR 0046 分组；合并不改变各列的值。
    迭代、差分与阻尼在 y 上；evaluate 仍接收 x = x(y)。求解手段：ADR 0043；ADR 0046；ADR 0054。
    """
    y = np.asarray(transform.from_x(np.asarray(x0, dtype=float)), dtype=float)
    scales = np.asarray(scales, dtype=float)
    residual, status, _reason = _at(y, evaluate, feasible, transform)
    if status == "domain":
        return NewtonResult("DOMAIN_STOP", "初值落在自由空间有效域外", transform.to_x(y), residual, 0, None)
    if status != "ok":
        return NewtonResult("NOT_CONVERGED", "初值不可行或残差不是有限数", transform.to_x(y), residual, 0, None)
    merit = merit_of(residual)
    if converged(residual):
        return _converged(y, residual, 0, merit, transform, final_check, journal)
    limit = int(config.damping_index_limit)
    explain = None if journal is None else journal.get("explain")
    for iteration in range(1, int(config.nu_max) + 1):
        _clock(deadline)
        trials = [] if journal is not None else None
        matrix = None
        delta = None
        try:
            if jacobian_layout is None:
                matrix = _jacobian(y, residual, evaluate, scales, feasible, transform=transform, blocks=blocks)
            else:
                matrix, stats = _grouped_jacobian(y, residual, evaluate, scales, feasible, jacobian_layout, transform=transform)
                if journal is not None:
                    journal.setdefault("jacobian", []).append(stats)
                if verify_jacobian and iteration == 1:
                    found = np.nonzero(matrix != _jacobian(y, residual, evaluate, scales, feasible, transform=transform))
                    if found[0].size:
                        raise JacobianMismatchError(*found)
                    if journal is not None:
                        journal["jacobian_verified"] = True
            delta = _newton_step(matrix, residual, y, scales, blocks)
            accepted = _damp(y, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials, explain, transform)
        except JacobianMismatchError as exc:
            if journal is not None:
                journal["jacobian_mismatch"] = {"rows": [int(item) for item in exc.rows], "columns": [int(item) for item in exc.columns]}
            _note(journal, iteration, merit, residual, delta, trials, scales, False)
            return NewtonResult("JACOBIAN_MISMATCH", str(exc), transform.to_x(y), residual, iteration, merit)
        except FreeboardDomainError as exc:
            _note(journal, iteration, merit, residual, delta, trials, scales, False)
            return NewtonResult("DOMAIN_STOP", str(exc), transform.to_x(y), residual, iteration, merit)
        except np.linalg.LinAlgError as exc:
            _note(journal, iteration, merit, residual, delta, trials, scales, False)
            return NewtonResult("NOT_CONVERGED", f"雅可比不能求解：{exc}", transform.to_x(y), residual, iteration, merit)
        _note(journal, iteration, merit, residual, delta, trials, scales, True)
        if journal is not None and final_check is not None:
            # 第 3 层检查结果（逐格惰性物）随每次迭代记录。ADR 0054 第 2 项第 3 层。
            journal["iterations"][-1]["final_check"] = final_check(transform.to_x(y))
        if accepted is None:
            _dump(journal, y, residual, matrix, delta, scales, transform)
            return NewtonResult("DAMPING_FAILED", "没有使归一残差严格下降的可行步", transform.to_x(y), residual, iteration, merit)
        y, residual, merit = accepted
        if converged(residual):
            return _converged(y, residual, iteration, merit, transform, final_check, journal)
    return NewtonResult("NOT_CONVERGED", "内层超过 nu_max", transform.to_x(y), residual, int(config.nu_max), merit)


def _converged(y, residual, iteration, merit, transform, final_check, journal):
    """收敛解上做第 3 层检查；不满足即停止并报告，不截断。求解手段：ADR 0054 第 2 项。"""
    x = transform.to_x(y)
    detail = None if final_check is None else final_check(x)
    if detail is not None and journal is not None:
        journal["final_check"] = detail
    if detail is not None and not detail.get("ok", True):
        return NewtonResult("FINAL_CHECK_FAILED", detail.get("message", "收敛解不满足第 3 层检查"), x, residual, iteration, merit)
    return NewtonResult("CONVERGED", "", x, residual, iteration, merit)



def normalized_merit(residual, blocks) -> float:
    """欧氏范数取在各块除以 ADR 0012 参考流之后。求解手段：ADR 0017 第 7 项。"""
    return float(np.linalg.norm(_normalize(residual, blocks)))


def _difference_basis(value, scale) -> float:
    """max(|y|, s_y)。求解手段：ADR 0043；ADR 0054 第 1 项。"""
    return max(abs(float(value)), float(scale))


def _column_scale(y, scales):
    """C 的对角元 max(|y|, s_y)。求解手段：ADR 0043。"""
    return np.array([_difference_basis(y[index], scales[index]) for index in range(np.size(y))])


def _jacobian(y, residual, evaluate, scales, feasible, *, transform=IDENTITY, blocks=None, columns=None):
    """前向差分。步长 √ε·max(|y|, s_y)；columns 给定时只算这些列。求解手段：ADR 0043；ADR 0054。"""
    chosen = range(y.size) if columns is None else columns
    matrix = np.zeros((y.size, y.size))
    for index in chosen:
        step = _SQRT_EPS * _difference_basis(y[index], scales[index])
        if not np.isfinite(step) or step <= 0.0:
            raise np.linalg.LinAlgError("差分步长不是正数")
        for _shrink in range(_DIFFERENCE_SHRINK + 1):
            trial = y.copy()
            trial[index] += step
            values, status, _reason = _at(trial, evaluate, feasible, transform)
            # 自由空间有效域在差分点上缩步；越出速率式范围按不可评价处理（ADR 0054 第 2 项第 1 层）。
            if status != "domain":
                break
            step *= 0.5
        if status != "ok":
            raise np.linalg.LinAlgError(f"差分点不可评价：{_reason}")
        matrix[:, index] = (values - residual) / step
    return matrix


def _grouped_jacobian(y, residual, evaluate, scales, feasible, layout, *, transform=IDENTITY):
    """节点序号模 3 相同、节点内序号相同的列同时扰动，每行只取邻域内被扰动的那一列；非相邻列单独逐列，组内不可评价时退回逐列。求解手段：ADR 0046。"""
    nodes, within, separate = layout
    matrix = np.zeros((y.size, y.size))
    rows, groups = {}, {}
    for row in range(y.size):
        rows.setdefault(int(nodes[row]), []).append(row)
    for column in range(y.size):
        if not separate[column]:
            groups.setdefault((int(nodes[column]) % 3, int(within[column])), []).append(column)
    fallbacks = 0
    for columns in groups.values():
        steps = {column: _SQRT_EPS * _difference_basis(y[column], scales[column]) for column in columns}
        if not all(np.isfinite(step) and step > 0.0 for step in steps.values()):
            raise np.linalg.LinAlgError("差分步长不是正数")
        trial = y.copy()
        for column, step in steps.items():
            trial[column] += step
        values, status, _reason = _at(trial, evaluate, feasible, transform)
        if status != "ok":
            fallbacks += 1
            matrix[:, columns] = _jacobian(y, residual, evaluate, scales, feasible, transform=transform, columns=columns)[:, columns]
            continue
        for column in columns:
            node = int(nodes[column])
            chosen = [row for item in (node - 1, node, node + 1) for row in rows.get(item, ())]
            matrix[chosen, column] = (values[chosen] - residual[chosen]) / steps[column]
    separate = np.flatnonzero(separate)
    matrix[:, separate] = _jacobian(y, residual, evaluate, scales, feasible, transform=transform, columns=separate)[:, separate]
    return matrix, {"fallbacks": fallbacks, "groups": len(groups), "separate": int(separate.size)}


def _newton_step(matrix, residual, y, scales, blocks):
    """解 (R·J·C)·z = −R·F，Δy = C·z。不显式求逆。求解手段：ADR 0043。"""
    column = _column_scale(y, scales)
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


def _damp(y, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials=None, explain=None, transform=IDENTITY):
    for index in range(limit + 1):
        _clock(deadline)
        point = y + delta / (2.0 ** index)
        values, status, reason = _at(point, evaluate, feasible, transform, explain)
        trial_merit = merit_of(values) if status == "ok" else None
        if trials is not None:
            trials.append({
                "i": index, "step": float(2.0 ** (-index)), "status": status, "reason": reason,
                "merit": None if trial_merit is None else float(trial_merit), "accepted": False,
            })
        if status != "ok":
            # 超域与越出速率式范围都只在试探点上缩步；起点与收敛解超域仍停止（ADR 0054 第 2 项第 1 层）。
            continue
        if trial_merit < merit:
            if trials is not None:
                trials[-1]["accepted"] = True
            return point, values, trial_merit
    return None


def _at(y, evaluate, feasible, transform=IDENTITY, explain=None):
    if not np.isfinite(y).all():
        return None, "reject", "y 非有限"
    x = transform.to_x(y)
    if not np.isfinite(x).all():
        return None, "reject", "x = x(y) 非有限"
    if not feasible(y):
        reason = "feasible 不满足"
        if explain is not None:
            try:
                reason = explain(y)
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


def _note(journal, iteration, merit, residual, delta, trials, scales, jacobian_ok=True) -> None:
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
        "largest": largest, "step": step, "jacobian_ok": bool(jacobian_ok),
        "trials": [] if trials is None else trials,
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


def _dump(journal, y, residual, matrix, delta, scales, transform=IDENTITY) -> None:
    if journal is None or not journal.get("failure_path") or matrix is None or delta is None:
        return
    path = Path(journal["failure_path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    blocks = journal.get("blocks")
    # dx 是 y 上的牛顿步；x 是它对应的物理点。ADR 0054 第 1 项。
    payload = {"y": np.asarray(y), "x": np.asarray(transform.to_x(y)), "F": np.asarray(residual), "J": np.asarray(matrix), "dx": np.asarray(delta), "scales": np.asarray(scales)}
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
