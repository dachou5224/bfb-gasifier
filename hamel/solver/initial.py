"""初值构造：固相按焦炭块 → 床区灰分配 → 自由空间灰块的顺序线性求解，并与预计算交替。

求解手段：ADR 0045。残差取 ``evaluate_residual`` 的对应固相行，不另写方程。
冻结核式的量取自本轮预计算；气相、温度与冻结量在块内固定。
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from hamel.solver.frozen import AssemblyError, FreeboardDomainError, TrialNotEvaluable
from hamel.solver.residual import evaluate_residual
from hamel.solver.sensitivity import diameter_at, section_volume_m3
from hamel.solver.vector import ReactorState, pack_state, unknown_slots, unpack_state

_SQRT_EPS = math.sqrt(float(np.finfo(float).eps))
# 焦炭行对灰未知量的导数：按 m_tot/m_dot_ref 归一后的上限。结构上应为零。
_DERIVATIVE_TOLERANCE = 1.0e-8


class InitialStateError(AssemblyError):
    """初值构造按 ADR 0045 第 4 项停止：线性系统奇异、负存量、床区无焦炭或回代不过门。"""


class DomainRejected(InitialStateError):
    """块求解或回代的评价点落在 ADR 0026 第 5 项的有效域外。state 是拒绝处的固相（未截断）。"""

    def __init__(self, message: str, state: ReactorState, phase: str):
        super().__init__(message)
        self.state = state
        self.phase = phase


class NegativeStock(InitialStateError):
    """块解出现负存量。items 给出变量名、解值与原值，供报告逐格逐粒级列出。"""

    def __init__(self, message: str, state: ReactorState, phase: str, items: list, condition: float, residual: float):
        super().__init__(message)
        self.state = state
        self.phase = phase
        self.items = items
        self.condition = condition
        self.residual = residual


@dataclass(frozen=True)
class SolidSolve:
    """一次块顺序求解的结果与检查量。"""

    state: ReactorState
    char_residual: float
    solid_residual: float
    char_ash_derivative: float
    minimum_char_kg: float
    bed_char_kg: float


def solve_solids(state: ReactorState, frozen, m_dot_ref: float, gate: float = 1.0e-6) -> SolidSolve:
    """焦炭块（线性）→ 床区灰分配 → 自由空间灰块（线性）→ 完整固相行回代。

    先核对焦炭行对灰未知量的导数为零；不成立即说明结构假定不成立，停止上报。
    """
    slots = unknown_slots(frozen)
    char_slots = tuple(slot for slot in slots if slot.kind == "char")
    ash_slots = tuple(slot for slot in slots if slot.kind == "ash")
    x = np.asarray(pack_state(state, frozen), dtype=float)
    scale = float(np.sum(frozen.geometric_mass_kg))
    base = _values(x, frozen, "起点")
    derivative = _char_ash_derivative(x, base, char_slots, ash_slots, scale, frozen)
    if derivative * scale > _DERIVATIVE_TOLERANCE * m_dot_ref:
        raise InitialStateError(
            f"焦炭行对灰未知量的导数不为零：归一值 {derivative * scale / m_dot_ref:.3e}，结构假定不成立"
        )
    x = _solve_block(x, base, char_slots, scale, m_dot_ref, frozen, "焦炭块")
    _require_bed_ash(x, frozen)
    x = _solve_block(x, _values(x, frozen, "床区灰分配"), ash_slots, scale, m_dot_ref, frozen, "自由空间灰块")
    checked = _values(x, frozen, "回代")
    solid = np.array([checked[slot.index] for slot in slots if slot.kind in ("char", "ash")])
    residual = float(np.max(np.abs(solid)) / m_dot_ref) if solid.size else 0.0
    if residual > gate:
        where = sorted(
            ((abs(checked[slot.index]), slot.name) for slot in slots if slot.kind in ("char", "ash")),
            reverse=True,
        )[0]
        raise InitialStateError(
            f"回代不满足固相门：{residual:.3e} > {gate:.1e}，最大项 {where[1]} = {where[0]:.3e} kg/s"
        )
    solved = unpack_state(x, frozen)
    char_rows = [checked[slot.index] for slot in char_slots]
    char_residual = float(np.max(np.abs(char_rows)) / m_dot_ref) if char_rows else 0.0
    char_state = (solved.bed_char_kg, solved.freeboard_char_free_kg, solved.freeboard_char_hover_kg)
    minimum = float(min(np.min(item) for item in char_state if item.size))
    return SolidSolve(solved, char_residual, residual, derivative, minimum, float(solved.bed_char_kg.sum()))


def interface_flows(frozen, state: ReactorState) -> dict:
    """床区各界面上、下行流与自由空间回落。与预检查报告的配对相同。"""
    coeff = frozen.b_free_per_s[-1] + frozen.b_hover_per_s[-1]
    mobile = state.freeboard_char_free_kg[-1] + state.freeboard_ash_free_kg[-1]
    return {
        "up": np.asarray(frozen.k_auf_per_s * frozen.geometric_mass_kg, dtype=float),
        "down": np.asarray(frozen.m_dot_ab_kg_s, dtype=float),
        "return": np.asarray(coeff * mobile, dtype=float),
    }


def solve_prelaunch(state: ReactorState, frozen, m_dot_ref: float, gate: float = 1.0e-6, records: list | None = None) -> dict:
    """固定床区焦炭，只解自由空间的焦炭块与灰块。求解手段：ADR 0049 第 3 步。

    回代只检查自由空间的焦炭行与灰行；床区固相行此时尚未求解，只记录残差。
    records 给出时，每个子块解完就追加一条记录，便于灰块失败时仍能报告第 3 步的存量与体积比。
    """
    slots = unknown_slots(frozen)
    char_slots = tuple(slot for slot in slots if slot.kind == "char" and slot.zone == "freeboard")
    ash_slots = tuple(slot for slot in slots if slot.kind == "ash")
    rows = [slot.index for slot in slots if slot.zone == "freeboard" and slot.kind in ("char", "ash")]
    bed = [slot.index for slot in slots if slot.zone == "bed" and slot.kind in ("char", "ash")]
    x = np.asarray(pack_state(state, frozen), dtype=float)
    scale = float(np.sum(frozen.geometric_mass_kg))
    x = _solve_block(x, _values(x, frozen, "预启动焦炭块"), char_slots, scale, m_dot_ref, frozen, "自由空间焦炭块")
    if records is not None:
        records.append({"stage": "预启动自由空间焦炭块（ADR 0049 第 3 步）", **_snapshot(unpack_state(x, frozen), frozen)})
    _require_bed_ash(x, frozen)
    x = _solve_block(x, _values(x, frozen, "预启动灰块"), ash_slots, scale, m_dot_ref, frozen, "自由空间灰块")
    checked = _values(x, frozen, "预启动回代")
    residual = float(np.max(np.abs(checked[rows])) / m_dot_ref) if rows else 0.0
    bed_residual = float(np.max(np.abs(checked[bed])) / m_dot_ref) if bed else 0.0
    solved = unpack_state(x, frozen)
    if residual > gate:
        largest = max(((abs(checked[slot.index]), slot.name) for slot in slots if slot.zone == "freeboard" and slot.kind in ("char", "ash")))
        raise InitialStateError(
            f"预启动回代不满足自由空间固相门：{residual:.3e} > {gate:.1e}，最大项 {largest[1]} = {largest[0]:.3e} kg/s"
        )
    if records is not None:
        records.append({"stage": "预启动自由空间灰块（ADR 0049 第 3 步）", **_snapshot(solved, frozen)})
    return {
        "state": solved,
        "freeboard_residual": residual,
        "bed_residual": bed_residual,
        **_snapshot(solved, frozen),
    }


def _snapshot(state, frozen) -> dict:
    """自由空间各格、各状态的存量与体积比。ADR 0049 报告要求。"""
    return {
        "char_free_kg": np.asarray(state.freeboard_char_free_kg, dtype=float).tolist(),
        "char_hover_kg": np.asarray(state.freeboard_char_hover_kg, dtype=float).tolist(),
        "ash_free_kg": np.asarray(state.freeboard_ash_free_kg, dtype=float).tolist(),
        "ash_hover_kg": np.asarray(state.freeboard_ash_hover_kg, dtype=float).tolist(),
        "volume_ratio": _volume_ratio(state, frozen),
    }


def _volume_ratio(state, frozen) -> list:
    """每格自由空间固体体积与 (1−ε_mf)·V 之比。只报告，不截断。"""
    ratios = []
    for cell in range(frozen.n_freeboard):
        solid = float(
            state.freeboard_char_free_kg[cell].sum() + state.freeboard_char_hover_kg[cell].sum()
            + state.freeboard_ash_free_kg[cell].sum() + state.freeboard_ash_hover_kg[cell].sum()
        )
        height = frozen.freeboard_height_m
        geometric = frozen.n_freeboard - 1 - cell
        z0 = geometric * height
        d0 = diameter_at(z0, frozen.bed_diameter_m, frozen.freeboard_diameter_m)
        d1 = diameter_at(z0 + height, frozen.bed_diameter_m, frozen.freeboard_diameter_m)
        volume = math.pi * (0.5 * d0) ** 2 * height if d0 == d1 else section_volume_m3(z0, z0 + height, frozen.bed_diameter_m, frozen.freeboard_diameter_m)
        ratios.append(solid / frozen.rho_p_kg_m3 / (volume * (1.0 - frozen.epsilon_mf)))
    return ratios


def construct(initial, precalc, m_dot_ref: float, *, max_loops: int = 5, tolerance: float = 1.0e-2, gate: float = 1.0e-6, records=None):
    """ADR 0049 的预启动（先解自由空间块）之后，按 ADR 0045 反复“解固相块 → 重做预计算”。

    相对变化以 ADR 0012 的固定参考尺度 m_dot_ref 作分母。做满 max_loops 次时以最后一次的固相为初值。
    records 给出时，每次循环与有效域拒绝都写进该列表，便于报告保留未完成的一步。
    求解手段：ADR 0045、ADR 0049。
    """
    log = [] if records is None else records
    result = initial
    started = _prelaunch(initial, precalc, m_dot_ref, gate, log)
    if started is None:
        return initial, log, "停止：预启动阶段未完成"
    result = started
    previous = interface_flows(result.frozen, result.state)
    reached = False
    for index in range(1, int(max_loops) + 1):
        try:
            solved = solve_solids(result.state, result.frozen, m_dot_ref, gate)
        except DomainRejected as exc:
            log.append({
                "loop": index, "domain_rejected": exc.phase, "solid_residual": None,
                "note": "未完成（有效域拒绝）", "reason": str(exc),
            })
            raise
        except NegativeStock as exc:
            log.append({
                "loop": index, "negative_count": len(exc.items), "solid_residual": None,
                "note": f"停止（负存量，{exc.phase}）", "reason": str(exc),
            })
            raise
        result = precalc(solved.state)
        current = interface_flows(result.frozen, result.state)
        change = {name: float(np.max(np.abs(current[name] - previous[name])) / m_dot_ref) for name in current}
        log.append({
            "loop": index,
            "char_ash_derivative": solved.char_ash_derivative,
            "char_residual": solved.char_residual,
            "solid_residual": solved.solid_residual,
            "minimum_char_kg": solved.minimum_char_kg,
            "bed_char_kg": solved.bed_char_kg,
            "up_kg_s": current["up"].tolist(),
            "down_kg_s": current["down"].tolist(),
            "return_kg_s": current["return"].tolist(),
            "change": change,
        })
        previous = current
        if max(change.values()) < tolerance:
            reached = True
            break
    note = "" if reached else f"交替 {int(max_loops)} 次后相对变化仍未全部小于 {tolerance:g}"
    return result, log, note


def _prelaunch(initial, precalc, m_dot_ref, gate, log):
    """ADR 0049 第 3、4 步：解自由空间块，再用这一状态重做预计算（回馈已解固相）。"""
    stages: list = []
    try:
        stage = solve_prelaunch(initial.state, initial.frozen, m_dot_ref, gate, stages)
    finally:
        log.extend(stages)
    record = {name: value for name, value in stage.items() if name != "state"}
    record.update({"stage": "预启动（ADR 0049 第 3 步）", "loops": None})
    before = interface_flows(initial.frozen, initial.state)
    result = precalc(stage["state"])
    after = interface_flows(result.frozen, result.state)
    record["interface_before"] = {name: np.asarray(value, dtype=float).tolist() for name, value in before.items()}
    record["interface_after"] = {name: np.asarray(value, dtype=float).tolist() for name, value in after.items()}
    log.append(record)
    return result


def _char_ash_derivative(x, base, char_slots, ash_slots, scale, frozen) -> float:
    """焦炭行对每个灰列的前向差分，按 m_tot 归一再除以 m_dot_ref 之前的值。"""
    if not char_slots or not ash_slots:
        return 0.0
    rows = np.array([slot.index for slot in char_slots])
    values = np.array([base[index] for index in rows])
    largest = 0.0
    for slot in ash_slots:
        step = _SQRT_EPS * max(abs(float(x[slot.index])), scale)
        trial = x.copy()
        trial[slot.index] += step
        changed = _values(trial, frozen, "焦炭行对灰列")[rows]
        largest = max(largest, float(np.max(np.abs(changed - values))) / step)
    return largest


def _solve_block(x, base, columns, scale, m_dot_ref, frozen, name):
    """稠密线性求解，带 ADR 0043 的行列缩放。块内残差对这批未知量线性，差分即系数。"""
    indices = np.array([slot.index for slot in columns])
    matrix = np.zeros((indices.size, indices.size))
    steps = np.array([_SQRT_EPS * max(abs(float(x[index])), scale) for index in indices])
    values = np.array([base[index] for index in indices])
    for position, index in enumerate(indices):
        trial = x.copy()
        trial[index] += steps[position]
        matrix[:, position] = (_values(trial, frozen, name)[indices] - values) / steps[position]
    column = np.array([max(abs(float(x[index])), scale) for index in indices])
    scaled = matrix * (column / m_dot_ref)[None, :]
    condition = float(np.linalg.cond(scaled))
    try:
        delta = column * np.linalg.solve(scaled, -values / m_dot_ref)
    except np.linalg.LinAlgError as exc:
        raise InitialStateError(f"{name}的线性系统奇异，条件数 {condition:.3e}") from exc
    if not np.isfinite(delta).all():
        raise InitialStateError(f"{name}的解不是有限数")
    updated = x.copy()
    updated[indices] += delta
    negative = [
        (columns[position].name, float(updated[indices[position]]), float(x[indices[position]]))
        for position in range(indices.size) if updated[indices[position]] < 0.0
    ]
    if negative:
        smallest = min(negative, key=lambda item: item[1])
        raise NegativeStock(
            f"{name}解出负存量：{smallest[0]} = {smallest[1]:.3e} kg，共 {len(negative)} 项，条件数 {condition:.3e}",
            unpack_state(updated, frozen), name,
            [{"name": name_, "solved_kg": value, "initial_kg": before} for name_, value, before in negative],
            condition, float(np.max(np.abs(values)) / m_dot_ref),
        )
    return updated


def _require_bed_ash(x, frozen) -> None:
    """床区灰分配 A_i,k = I_i·C_i,k/ΣC_i,k 的定义条件：ΣC > 0 且 I_i ≥ 0。"""
    state = unpack_state(x, frozen)
    total = state.bed_char_kg.sum(axis=1)
    inert = frozen.geometric_mass_kg - total - frozen.fuel_bed_kg.sum(axis=1)
    for cell in range(frozen.n_bed):
        if total[cell] <= 0.0:
            raise InitialStateError(f"床区格 {cell}（衡算下标）的焦炭总量为零，灰分配没有定义")
        if inert[cell] < 0.0:
            raise InitialStateError(f"床区格 {cell}（衡算下标）的惰性物为负：{inert[cell]:.3e} kg")


def _values(x, frozen, phase: str = "") -> np.ndarray:
    """取残差。有效域拒绝换成带状态与阶段的上报：不截断存量，也不关闭检查。"""
    point = np.asarray(x, dtype=float)
    try:
        return np.asarray(evaluate_residual(point, frozen).values, dtype=float)
    except TrialNotEvaluable as exc:
        raise InitialStateError(f"初值构造的残差没有定义：{exc}") from exc
    except FreeboardDomainError as exc:
        raise DomainRejected(f"{phase or '求解'}的有效域拒绝：{exc}", unpack_state(point, frozen), phase) from exc
