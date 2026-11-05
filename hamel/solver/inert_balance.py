"""惰性料衡算校正（Bild 2.2 „Abgleich der Inertmaterialbilanz“）：ADR 0055 第 2、3 项。

E(w) 在上一轮内层解上冻结气相、温度、各粒级焦炭存量与床料总量，只改床灰的细灰份额 w，
只以自由空间灰存量（两个状态、各格各粒级）为未知量、只取灰衡算行，重新求自由空间灰的稳态，
取顶部带出灰（闭合 C 的改送出口灰计入，闭合 B 的改送回床不计入）。解释性选择：ADR 0055。
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from hamel.solver.frozen import AssemblyError
from hamel.solver.newton import IDENTITY, damped_newton
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots, unpack_state

FINE_TOLERANCE = 1.0e-6
MAX_EVALUATIONS = 50
ASH_RELATIVE_TOLERANCE = 1.0e-4
_SUBSOLVE_GATE = 1.0e-6


class InertBalanceError(AssemblyError):
    """ADR 0055 第 3、5 项：端点响应越界、子求解不收敛或求根失败。报告并停止，不截断。"""

    def __init__(self, message: str, info: dict | None = None):
        super().__init__(message)
        self.info = {} if info is None else info


def _text(reason: str, info: dict) -> str:
    return (
        f"{reason}：T = {info['target']:.6e} kg/s，E₀ = {info['E0']:.6e}，E₁ = {info['E1']:.6e}，"
        f"w_f = {info['initial']:.6e}，E 评价 {info['evaluations']} 次"
    )


def fine_share(response, initial: float, target: float, *, tolerance: float = FINE_TOLERANCE, limit: int = MAX_EVALUATIONS) -> tuple:
    """端点响应加 Illinois 括号求根，解 E(w) = T。第 5 项三条停止成立时抛 ``InertBalanceError``。"""
    goal = float(target)
    if not np.isfinite(goal) or goal <= 0.0:
        raise InertBalanceError(f"飞灰目标 T 不是正数：{goal}")
    seen: dict[float, float] = {}

    def at(w: float) -> float:
        key = float(w)
        if key not in seen:
            seen[key] = float(response(key))
        return seen[key]

    e0, e1, e = at(0.0), at(1.0), at(initial)
    info = {"target": goal, "E0": e0, "E1": e1, "initial": float(initial), "E": e, "w": float(initial)}
    if abs(e - goal) / goal < tolerance:
        info.update(stop="satisfied", evaluations=len(seen))
        return float(initial), info
    coarse = (e0 - goal) / goal >= tolerance
    if coarse or (goal - e1) / goal >= tolerance:
        info.update(stop="coarse" if coarse else "fine", E=e0 if coarse else e1, w=0.0 if coarse else 1.0,
                    evaluations=len(seen))
        raise InertBalanceError(_text("粗灰带出超过目标" if coarse else "细灰不足", info), info)
    a, b, fa, fb = 0.0, 1.0, e0 - goal, e1 - goal
    for _ in range(int(limit)):
        if fb == fa:
            break
        w = (a * fb - b * fa) / (fb - fa)
        residual = at(w) - goal
        if abs(residual) / goal < tolerance:
            info.update(stop="root", w=float(w), E=at(w), evaluations=len(seen))
            return float(w), info
        if residual * fb > 0.0:
            b, fb, fa = w, residual, 0.5 * fa
        else:
            a, fa, fb = w, residual, 0.5 * fb
    info.update(stop="failed", evaluations=len(seen))
    raise InertBalanceError(_text("求根失败", info), info)


def free_ash_state(frozen, state, config, m_ref: float) -> tuple:
    """只以自由空间灰存量为未知量、只取灰衡算行求稳态（起点：``state`` 的灰存量）。

    返回 (完整未知向量, 迭代数)。子求解不收敛或有灰存量小于 −τ 时抛 ``InertBalanceError``。
    """
    rows = np.array([slot.index for slot in unknown_slots(frozen) if slot.kind == "ash"], dtype=int)
    if rows.size == 0:
        raise InertBalanceError("自由空间没有灰未知量，E(w) 没有定义")
    base = pack_state(state, frozen)
    scale = float(np.sum(frozen.geometric_mass_kg))

    def values(y):
        x = base.copy()
        x[rows] = y
        return evaluate_residual(x, frozen).values[rows]

    outcome = damped_newton(
        base[rows], values, np.full(rows.size, scale),
        lambda r: bool(np.max(np.abs(r)) / m_ref <= _SUBSOLVE_GATE),
        lambda y: bool(np.isfinite(y).all()), lambda r: float(np.linalg.norm(r) / m_ref),
        config, transform=IDENTITY,
    )
    if outcome.status != "CONVERGED":
        raise InertBalanceError(f"E(w) 的子求解未收敛：{outcome.status}（{outcome.reason}），迭代 {outcome.iterations}")
    x = base.copy()
    x[rows] = np.asarray(outcome.x, dtype=float)
    tau = float(config.stock_roundoff_relative) * scale
    if np.any(x[rows] < -tau):
        raise InertBalanceError(f"E(w) 的子求解解出负灰存量：最小 {float(np.min(x[rows])):.3e} kg < −τ = {-tau:.3e}")
    return x, int(outcome.iterations)


def top_ash(frozen, values) -> float:
    """顶部带出灰。闭合 C 的改送出口灰计入（ADR 0026 第 6 项），闭合 B 的改送回床不计入。"""
    state = unpack_state(values, frozen)
    total = float(np.sum(frozen.a_free_per_s[0] * state.freeboard_ash_free_kg[0]))
    if frozen.hover_closure == "C":
        total += float(evaluate_residual(values, frozen).redirect_ash_kg_s)
    return total


def update_fine_share(frozen, state, config, m_ref: float) -> tuple:
    """每次预计算更新 w_f：E₀、E₁ → 已满足就保留 → 查停止 → Illinois 求根。ADR 0055 第 3 项。"""
    iterations: list[int] = []
    cache: dict[float, float] = {}

    def response(w: float) -> float:
        key = float(w)
        if key not in cache:
            trial = replace(frozen, bed_ash_fine_share=key)
            x, count = free_ash_state(trial, state, config, m_ref)
            iterations.append(count)
            cache[key] = top_ash(trial, x)
        return cache[key]

    w, info = fine_share(response, float(frozen.bed_ash_fine_share), float(frozen.ash_target_kg_s))
    half = response(0.5)
    middle = 0.5 * (info["E0"] + info["E1"])
    info.update(E_half=half, linearity=None if middle == 0.0 else abs(half - middle) / abs(middle),
                subsolve_iterations=iterations)
    return float(w), info
