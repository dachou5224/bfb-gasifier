"""初值构造：第一次预计算之后，对固相子系统做牛顿求解，再与预计算交替。

求解手段：ADR 0052（取代 ADR 0045 第 3 项的块顺序求解与 ADR 0049）。
未知量与残差取 ``unknown_slots``／``evaluate_residual`` 的固相部分，不另写方程。
每次评价残差时按当前状态重算排灰（ADR 0014 A5）与床区整体向下流（ADR 0009 第 2 项），
复用 ``precalc._drain`` 与 ``precalc.bed_downflow``。气相、温度与其余冻结量固定。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.solver.frozen import AssemblyError, FreeboardDomainError, TrialNotEvaluable
from hamel.solver.gates import arrays_allowed
from hamel.solver.newton import damped_newton
from hamel.solver.precalc import PrecalcError, _drain, bed_downflow, freeboard_return_kg_s
from hamel.solver.residual import evaluate_residual
from hamel.solver.sensitivity import diameter_at, section_volume_m3
from hamel.solver.vector import ReactorState, pack_state, unknown_slots, unpack_state


class InitialStateError(AssemblyError):
    """初值构造按 ADR 0052 第 3 项停止：阻尼失败、超过 ν_max、收敛解越界或预计算报错。

    diagnostics 带停止时的牛顿记录，供报告写明哪一步、哪一项不可行。
    """

    def __init__(self, message: str, diagnostics: dict | None = None):
        super().__init__(message)
        self.diagnostics = {} if diagnostics is None else diagnostics


class DomainRejected(InitialStateError):
    """评价点落在 ADR 0026 第 5 项的有效域外。state 是拒绝处的固相（未截断）。"""

    def __init__(self, message: str, state: ReactorState, phase: str):
        super().__init__(message)
        self.state = state
        self.phase = phase


@dataclass(frozen=True)
class SolidSolve:
    """一次固相牛顿的结果与报告量。求解手段：ADR 0052 第 1 项。"""

    state: ReactorState
    iterations: int
    merit_history: tuple[float, ...]
    solid_residual: float
    drain_ash_kg_s: float
    drain_char_kg_s: np.ndarray
    m_dot_ab_kg_s: np.ndarray
    frozen_down_kg_s: np.ndarray
    bottom_closure_kg_s: float
    bed_char_kg: float
    bed_inert_kg: np.ndarray


def interface_flows(frozen, state: ReactorState) -> dict:
    """床区各界面上、下行流与自由空间回落。与预检查报告的配对相同。"""
    coeff = frozen.b_free_per_s[-1] + frozen.b_hover_per_s[-1]
    mobile = state.freeboard_char_free_kg[-1] + state.freeboard_ash_free_kg[-1]
    return {
        "up": np.asarray(frozen.k_auf_per_s * frozen.geometric_mass_kg, dtype=float),
        "down": np.asarray(frozen.m_dot_ab_kg_s, dtype=float),
        "return": np.asarray(coeff * mobile, dtype=float),
    }


def recompute_flows(state: ReactorState, frozen, context, reaction=None):
    """按当前状态重算排灰（ADR 0014 A5）与床区整体向下流（ADR 0009 第 2 项）。

    写法复用预计算的同一函数 ``_drain`` 与 ``bed_downflow``：上行流取冻结的 k_auf·M_i；
    回落取当前自由空间最下格的自由焦炭与灰存量乘 b^f 再加冻结的燃料回落；反应汇取当前
    床区焦炭消耗；排灰取前一项。求解手段：ADR 0052 第 1 项第 2 步。
    """
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    drain_ash, drain_char = _drain(
        context.case, state, context.released, inert, frozen.geometric_mass_kg, context.books,
        context.redirected_ash_out_kg_s, context.redirected_fuel_out_kg_s,
    )
    returned = freeboard_return_kg_s(
        frozen.b_free_per_s[-1], state.freeboard_char_free_kg[-1], state.freeboard_ash_free_kg[-1],
        float(frozen.fuel_down_free_kg_s[-1].sum()),
    ) + float(context.redirected_return_kg_s)
    down, closure = bed_downflow(
        frozen.k_auf_per_s * frozen.geometric_mass_kg, context.case.fuel_kg_s, reaction, returned,
        float(drain_ash + np.sum(drain_char)),
    )
    return (
        np.asarray(inert, dtype=float), float(drain_ash), np.asarray(drain_char, dtype=float),
        np.asarray(down, dtype=float), float(closure), float(returned),
    )


def solid_feasible(state: ReactorState, frozen, context) -> bool:
    """ADR 0017 第 7 项的可行条件取固相部分，另加排灰灰分量非负。

    求解手段：ADR 0052 第 1 项第 3 步。气相与温度固定，故不重复检查。
    """
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    stocks = (
        state.bed_char_kg, state.freeboard_char_free_kg, state.freeboard_char_hover_kg,
        state.freeboard_ash_free_kg, state.freeboard_ash_hover_kg,
    )
    temperatures = np.concatenate([state.bed_temperature_k, state.freeboard_temperature_k])
    if not arrays_allowed((), stocks, temperatures, inert):
        return False
    try:
        _, drain_ash, _, _, _, _ = recompute_flows(state, frozen, context)
    except (PrecalcError, ArithmeticError):
        return False
    return bool(drain_ash >= 0.0)


def _infeasible_clause(state: ReactorState, frozen, context) -> str:
    """按 ADR 0017 第 7 项的顺序报告第一条不满足的可行条件。只用于阻尼记录。"""
    values = np.asarray(pack_state(state, frozen), dtype=float)
    for slot in unknown_slots(frozen):
        if slot.kind in ("char", "ash") and values[slot.index] < 0.0:
            return f"{slot.name} < 0"
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    if np.any(inert < 0.0):
        return f"m_inert < 0（床区格 {int(np.flatnonzero(inert < 0.0)[0])}）"
    try:
        _, drain_ash, _, _, _, _ = recompute_flows(state, frozen, context)
    except (PrecalcError, ArithmeticError) as exc:
        return f"排灰没有定义：{exc}"
    if drain_ash < 0.0:
        return f"排灰灰分量为负：{drain_ash:.3e} kg/s"
    return "可行条件不满足"


def require_converged_solution(state: ReactorState, frozen, residual, flows, m_dot_ref: float, gate: float) -> float:
    """ADR 0052 第 3 项：收敛解上任一存量为负、床格惰性物为负、排灰灰分量为负或最底床格
    灰分为零时停止并报告；同时复核固相门。返回固相行的归一最大残差。"""
    inert, drain_ash, _, _, _, _ = flows
    values = np.asarray(pack_state(state, frozen), dtype=float)
    slots = [slot for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")]
    negative = [(slot.name, float(values[slot.index])) for slot in slots if values[slot.index] < 0.0]
    if negative:
        text = ", ".join(f"{name} = {value:.3e} kg" for name, value in negative[:6])
        raise InitialStateError(f"初值构造的收敛解出现 {len(negative)} 项负存量：{text}")
    bad = [(cell, float(value)) for cell, value in enumerate(inert) if value < 0.0]
    if bad:
        raise InitialStateError(f"初值构造的收敛解在床区格 {bad[0][0]}（衡算下标）的惰性物为负：{bad[0][1]:.3e} kg")
    if drain_ash < 0.0:
        raise InitialStateError(f"初值构造的收敛解排灰灰分量为负：{drain_ash:.3e} kg/s")
    fraction = float(inert[-1] / frozen.geometric_mass_kg[-1])
    if fraction <= 0.0:
        raise InitialStateError("初值构造的收敛解最底床格灰分为零，排灰总流率没有定义")
    norm = float(np.max(np.abs(residual)) / m_dot_ref)
    if norm > gate:
        raise InitialStateError(f"初值构造的收敛解不满足固相门：{norm:.3e} > {gate:.1e}")
    return norm


def solve_solids(state: ReactorState, frozen, context, m_dot_ref: float, *, gate: float = 1.0e-6, config=None, journal: dict | None = None) -> SolidSolve:
    """以全部固相为未知量的阻尼牛顿。求解手段：ADR 0052 第 1 项第 3 步。

    残差取原残差的固相行；每次评价按当前状态重算排灰与床区向下流。差分步长与行列缩放
    按 ADR 0043（由 ``damped_newton`` 给出），阻尼、可行条件与迭代上限按 ADR 0017 第 7、8 项，
    可行条件另加排灰灰分量非负。试探点的有效域拒绝按不可行缩步；起点与收敛解超域时停止。
    """
    limits = SolverConfig() if config is None else config
    slots = unknown_slots(frozen)
    rows = np.array([slot.index for slot in slots if slot.kind in ("char", "ash")], dtype=int)
    if rows.size == 0:
        raise InitialStateError("初值构造的固相子系统没有未知量")
    base = np.asarray(pack_state(state, frozen), dtype=float)
    scale = float(np.sum(frozen.geometric_mass_kg))
    scales = np.full(rows.size, scale)

    def point(y, trial: bool = True):
        x = base.copy()
        x[rows] = y
        current = unpack_state(x, frozen)
        try:
            reaction = -np.asarray(evaluate_residual(x, frozen).char_consumption_bed_kg_s, dtype=float).sum(axis=1)
            flows = recompute_flows(current, frozen, context, reaction)
            _, drain_ash, drain_char, down, _, _ = flows
            updated = replace(frozen, drain_ash_kg_s=drain_ash, drain_char_kg_s=drain_char, m_dot_ab_kg_s=down)
            values = np.asarray(evaluate_residual(x, updated).values, dtype=float)
        except FreeboardDomainError as exc:
            if trial:
                raise TrialNotEvaluable(f"有效域拒绝：{exc}") from exc
            raise
        except (PrecalcError, ArithmeticError) as exc:
            if trial:
                raise TrialNotEvaluable(f"排灰没有定义：{exc}") from exc
            raise PrecalcError(f"初值构造的排灰没有定义：{exc}") from exc
        return values[rows], current, flows

    def feasible(y) -> bool:
        x = base.copy()
        x[rows] = y
        return solid_feasible(unpack_state(x, frozen), frozen, context)

    def converged(residual) -> bool:
        return bool(np.max(np.abs(residual)) / m_dot_ref <= gate)

    def merit(residual) -> float:
        return float(np.linalg.norm(residual)) / m_dot_ref

    def explain(y) -> str:
        x = base.copy()
        x[rows] = y
        return _infeasible_clause(unpack_state(x, frozen), frozen, context)

    log = {} if journal is None else journal
    log.setdefault("iterations", [])
    log.setdefault("names", [slots[index].name for index in rows])
    log.setdefault("explain", explain)
    try:
        point(base[rows], trial=False)
    except FreeboardDomainError as exc:
        raise DomainRejected(f"初值构造起点的有效域拒绝：{exc}", state, "起点") from exc
    except (PrecalcError, TrialNotEvaluable) as exc:
        raise InitialStateError(f"初值构造起点不可评价：{exc}") from exc
    outcome = damped_newton(base[rows], lambda y: point(y)[0], scales, converged, feasible, merit, limits, journal=log)
    if outcome.status != "CONVERGED":
        diagnostics = {
            "status": outcome.status, "reason": outcome.reason,
            "iterations": int(outcome.iterations), "nu_max": int(limits.nu_max),
            "zero_stocks": int(np.sum(base[rows] == 0.0)),
            "merit_history": [entry.get("merit") for entry in log.get("iterations", ())],
            "journal": {key: value for key, value in log.items() if key != "explain"},
        }
        raise InitialStateError(
            f"初值构造的固相牛顿未收敛：{outcome.status}（{outcome.reason}），迭代 {outcome.iterations}，"
            f"ν_max = {int(limits.nu_max)}",
            diagnostics,
        )
    x = base.copy()
    x[rows] = np.asarray(outcome.x, dtype=float)
    try:
        checked, current, flows = point(x[rows], trial=False)
    except FreeboardDomainError as exc:
        raise DomainRejected(f"初值构造收敛解的有效域拒绝：{exc}", unpack_state(x, frozen), "收敛解") from exc
    except (PrecalcError, TrialNotEvaluable) as exc:
        raise InitialStateError(f"初值构造收敛解不可评价：{exc}") from exc
    inert, drain_ash, drain_char, down, closure, _returned = flows
    residual = require_converged_solution(current, frozen, checked, flows, m_dot_ref, gate)
    history = tuple(float(entry["merit"]) for entry in log.get("iterations", ()))
    history += (float(outcome.merit) if outcome.merit is not None else residual,)
    return SolidSolve(
        current, int(outcome.iterations), history, residual, drain_ash, drain_char, down,
        np.asarray(frozen.m_dot_ab_kg_s, dtype=float), closure, float(current.bed_char_kg.sum()), inert,
    )


def construct(initial, precalc, m_dot_ref: float, *, gate: float = 1.0e-6, tolerance: float = 1.0e-2, max_loops: int = 5, config=None, records=None):
    """第一次预计算之后，反复“固相牛顿 → 按 driver.solid_feedback 回馈并重做预计算”。

    precalc(state, drain) 接收解出的固相与它按同一状态重算的排灰（灰分量，各粒级焦炭），
    后者作下一次预计算的 previous_drain_kg_s。交替判据与次数沿用 ADR 0045 第 3 项第 5 步：
    相邻两次预计算的床区界面上下行流与自由空间回落的相对变化都小于 tolerance，或做满
    max_loops 次。求解手段：ADR 0052 第 1 项第 5 步。
    """
    limits = SolverConfig() if config is None else config
    log = [] if records is None else records
    result = initial
    previous = interface_flows(result.frozen, result.state)
    reached = False
    for index in range(1, int(max_loops) + 1):
        try:
            solved = solve_solids(result.state, result.frozen, result.solid_context, m_dot_ref, gate=gate, config=limits)
        except InitialStateError as exc:
            inert = (
                result.frozen.geometric_mass_kg - result.state.bed_char_kg.sum(axis=1)
                - result.frozen.fuel_bed_kg.sum(axis=1)
            )
            log.append({
                "loop": index, "stopped": True, "reason": str(exc),
                "frozen_down_kg_s": np.asarray(result.frozen.m_dot_ab_kg_s, dtype=float).tolist(),
                "drain_ash_kg_s": float(result.frozen.drain_ash_kg_s),
                "drain_char_kg_s": np.asarray(result.frozen.drain_char_kg_s, dtype=float).tolist(),
                "bed_char_kg": float(result.state.bed_char_kg.sum()),
                "bed_inert_kg": np.asarray(inert, dtype=float).tolist(),
                **(exc.diagnostics or {}),
            })
            raise
        result = precalc(solved.state, (solved.drain_ash_kg_s, np.asarray(solved.drain_char_kg_s, dtype=float)))
        current = interface_flows(result.frozen, result.state)
        change = {name: float(np.max(np.abs(current[name] - previous[name])) / m_dot_ref) for name in current}
        log.append(_record(index, solved, result, current, change))
        previous = current
        if max(change.values()) < tolerance:
            reached = True
            break
    note = "" if reached else f"交替 {int(max_loops)} 次后相对变化仍未全部小于 {tolerance:g}"
    return result, log, note


def _record(index: int, solved: SolidSolve, result, current: dict, change: dict) -> dict:
    """一次交替的报告量。ADR 0052 报告项：存量、排灰、带出、回落、两类 Ṁ_ab、闭合差、体积比。"""
    state, frozen = result.state, result.frozen
    top = state.freeboard_ash_free_kg[0] + state.freeboard_ash_hover_kg[0]
    elutriation = float(np.sum((frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * top))
    return {
        "loop": index,
        "iterations": solved.iterations,
        "residual_history": list(solved.merit_history),
        "solid_residual": solved.solid_residual,
        "bed_char_kg": solved.bed_char_kg,
        "bed_inert_kg": np.asarray(solved.bed_inert_kg, dtype=float).tolist(),
        "drain_ash_kg_s": solved.drain_ash_kg_s,
        "drain_char_kg_s": np.asarray(solved.drain_char_kg_s, dtype=float).tolist(),
        "elutriated_ash_kg_s": elutriation,
        "return_kg_s": float(np.sum(current["return"])),
        "frozen_down_kg_s": np.asarray(solved.frozen_down_kg_s, dtype=float).tolist(),
        "down_kg_s": np.asarray(solved.m_dot_ab_kg_s, dtype=float).tolist(),
        "up_kg_s": np.asarray(current["up"], dtype=float).tolist(),
        "bottom_closure_kg_s": solved.bottom_closure_kg_s,
        "hover_char_kg": float(np.sum(state.freeboard_char_hover_kg)),
        "volume_ratio": _volume_ratio(state, frozen),
        "change": change,
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
