"""初值构造的诊断：固定初值气相与温度，固相、排灰与床区向下流联立牛顿。只报告。

不改模型。做法：
1. 运行到初值构造之前，取初值状态与预计算函数；
2. 自由空间与床区温度统一设为 T，做一次预计算，得冻结量；
3. 对全部固相未知量做牛顿迭代。每次评价残差时，按当前状态重算两项冻结量：
   - 排灰（ADR 0014 A5）：灰分量 = 燃料带入灰 − 带出灰，总流率按最底床格的灰分；
   - 床区整体向下流 Ṁ_ab（ADR 0009 第 2 项递推）：回落取当前自由空间最下格的自由存量，
     反应汇取当前床区焦炭消耗，排灰取上式；
   对照组（--frozen-down）只重算排灰，Ṁ_ab 用预计算冻结值；
4. 迭代中暂停有效域检查（只作诊断），收敛后报告床区惰性物、排灰、带出灰、体积比。

结果写入 reports/research-initial-solids/。
"""

from __future__ import annotations

import argparse
import json
import signal
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.driver as driver
import hamel.solver.initial as initial
import hamel.solver.precalc as precalc
import hamel.solver.residual as residual
from hamel.cases.config import SolverConfig
from hamel.solver.vector import pack_state, unknown_slots, unpack_state

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "research-initial-solids"


class _Captured(Exception):
    pass


class _Timeout(BaseException):
    """继承 BaseException：任何 except Exception 都不会吞掉全局期限（CHARTER §7）。"""


def _capture(context: dict) -> None:
    """运行到初值构造之前；预计算中的排灰余额为负时记 0（排灰在牛顿内按状态重算）。"""
    original_drain = precalc._drain

    def drain(case, state, released, inert, masses, books, ash_out=0.0, fuel_out=0.0):
        context["case"] = case
        try:
            return original_drain(case, state, released, inert, masses, books, ash_out, fuel_out)
        except precalc.PrecalcError as exc:
            if "带出灰" not in str(exc):
                raise
            return 0.0, np.zeros(np.shape(state.bed_char_kg)[1])

    def construct(first, function, m_ref, **kwargs):
        context.update(first=first, precalc=function, m_ref=m_ref)
        raise _Captured()

    precalc._drain = drain
    driver.construct = construct
    try:
        driver.solve_case(str(_CASE), "1", deadline=lambda: False, config=replace(SolverConfig(), outer_round_limit=1))
    except _Captured:
        pass


def _a5(x, frozen, case):
    state = unpack_state(x, frozen)
    moisture = float(case.physical.moisture_mass_fraction)
    ash_in = case.fuel_kg_s * float(case.physical.ash_mass_fraction_wf) * (1.0 - moisture)
    carried = float(np.sum(np.asarray(state.freeboard_ash_free_kg[0]) * np.asarray(frozen.a_free_per_s[0])))
    fuel_out = float(np.sum(frozen.fuel_up_free_kg_s[0])) + float(np.sum(frozen.fuel_exit_bed_kg_s))
    drain_ash = ash_in - carried - ash_in / case.fuel_kg_s * fuel_out
    bottom = float(frozen.geometric_mass_kg[-1])
    char = np.asarray(state.bed_char_kg[-1], dtype=float)
    fraction = (bottom - char.sum() - float(np.sum(frozen.fuel_bed_kg[-1]))) / bottom
    total = drain_ash / fraction if fraction > 0.0 else np.inf
    return drain_ash, total * char / bottom, {"ash_in": ash_in, "carried": carried, "fraction": fraction}


def _consistent(x, frozen, case, recompute_down: bool):
    drain_ash, drain_char, _ = _a5(x, frozen, case)
    updated = replace(frozen, drain_ash_kg_s=float(drain_ash), drain_char_kg_s=np.asarray(drain_char, dtype=float))
    if not recompute_down:
        return updated, None
    assembly = residual.evaluate_residual(x, updated)
    consumption = np.asarray(assembly.char_consumption_bed_kg_s, dtype=float).sum(axis=1)
    state = unpack_state(x, frozen)
    returned = float(np.sum(np.asarray(frozen.b_free_per_s[-1]) * (
        np.asarray(state.freeboard_char_free_kg[-1]) + np.asarray(state.freeboard_ash_free_kg[-1])
    ))) + float(np.sum(frozen.fuel_down_free_kg_s[-1]))
    up = np.asarray(frozen.k_auf_per_s) * np.asarray(frozen.geometric_mass_kg)
    down, closure = precalc.bed_downflow(up, case.fuel_kg_s, -consumption, returned, float(drain_ash + np.sum(drain_char)))
    info = {"down": np.asarray(down).tolist(), "frozen_down": np.asarray(frozen.m_dot_ab_kg_s).tolist(),
            "up": up.tolist(), "returned": returned, "consumption": float(consumption.sum()), "closure": float(closure)}
    return replace(updated, m_dot_ab_kg_s=np.asarray(down, dtype=float)), info


def _newton(x, frozen, case, m_ref, recompute_down, limit=40):
    def values(point):
        current, _ = _consistent(point, frozen, case, recompute_down)
        return np.asarray(residual.evaluate_residual(point, current).values, dtype=float)

    rows = np.array([slot.index for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")])
    scale = float(np.sum(frozen.geometric_mass_kg))
    history = []
    for _ in range(limit):
        r = values(x)
        if not np.all(np.isfinite(r[rows])):
            raise FloatingPointError("残差不是有限数（最底床格灰分降到 0 时排灰无定义）")
        norm = float(np.max(np.abs(r[rows])) / m_ref)
        history.append(norm)
        if norm < 1.0e-8:
            break
        steps = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(x[rows]), scale)
        jacobian = np.zeros((rows.size, rows.size))
        for position, index in enumerate(rows):
            trial = x.copy()
            trial[index] += steps[position]
            jacobian[:, position] = (values(trial)[rows] - r[rows]) / steps[position]
        column = np.maximum(np.abs(x[rows]), scale)
        delta = column * np.linalg.solve(jacobian * (column / m_ref)[None, :], -r[rows] / m_ref)
        damping = 1.0
        for _ in range(30):
            trial = x.copy()
            trial[rows] += damping * delta
            try:
                if np.max(np.abs(values(trial)[rows])) < np.max(np.abs(r[rows])) * (1.0 - 1.0e-4 * damping):
                    break
            except _Timeout:
                raise
            except Exception:
                pass
            damping *= 0.5
        x = trial
    return x, history


def _summary(x, frozen, case, recompute_down, temperature, history):
    state = unpack_state(x, frozen)
    bed = np.asarray(state.bed_char_kg, dtype=float)
    inert = np.asarray(frozen.geometric_mass_kg) - bed.sum(axis=1) - np.asarray(frozen.fuel_bed_kg).sum(axis=1)
    drain_ash, drain_char, info = _a5(x, frozen, case)
    _, flows = _consistent(x, frozen, case, recompute_down)
    solids = [slot.index for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")]
    return {
        "temperature_k": temperature,
        "recompute_down": recompute_down,
        "residual_history": history,
        "converged": bool(history and history[-1] < 1.0e-8),
        "bed_char_kg": float(bed.sum()),
        "bed_char_by_class_kg": bed.sum(axis=0).tolist(),
        "bed_inert_by_cell_kg": inert.tolist(),
        "minimum_solid_kg": float(np.min(x[solids])),
        "drain_ash_kg_s": float(drain_ash),
        "drain_char_kg_s": float(np.sum(drain_char)),
        "elutriated_ash_kg_s": info["carried"],
        "ash_in_kg_s": info["ash_in"],
        "bottom_ash_fraction": info["fraction"],
        "hover_char_kg": float(np.sum(state.freeboard_char_hover_kg)),
        "volume_ratio": initial._volume_ratio(state, frozen),
        "flows": flows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--temperatures", default="1073,1123,1173")
    parser.add_argument("--timeout-s", type=float, default=3600.0)
    args = parser.parse_args()
    _OUT.mkdir(parents=True, exist_ok=True)

    def expire(signum, frame):
        raise _Timeout()

    signal.signal(signal.SIGALRM, expire)
    signal.alarm(int(args.timeout_s))
    context: dict = {}
    runs = []
    status = "ERROR"
    check = residual._Context.domain
    try:
        _capture(context)
        residual._Context.domain = lambda self: None  # 只作诊断
        first, function, m_ref = context["first"], context["precalc"], context["m_ref"]
        for temperature in (float(item) for item in args.temperatures.split(",")):
            hot = replace(
                first.state,
                bed_temperature_k=np.full(np.shape(first.state.bed_temperature_k), temperature),
                freeboard_temperature_k=np.full(np.shape(first.state.freeboard_temperature_k), temperature),
            )
            result = function(hot)
            for recompute_down in (False, True):
                x = np.asarray(pack_state(result.state, result.frozen), dtype=float)
                try:
                    x, history = _newton(x, result.frozen, context["case"], m_ref, recompute_down)
                    runs.append(_summary(x, result.frozen, context["case"], recompute_down, temperature, history))
                except (_Timeout, KeyboardInterrupt):
                    raise
                except Exception as exc:  # 只作诊断：记下失败，继续下一组
                    runs.append({"temperature_k": temperature, "recompute_down": recompute_down, "converged": False,
                                 "error": f"{type(exc).__name__}: {exc}"})
                    print(json.dumps(runs[-1], ensure_ascii=False), flush=True)
                    continue
                print(json.dumps({key: runs[-1][key] for key in ("temperature_k", "recompute_down", "converged", "bed_char_kg", "bed_inert_by_cell_kg", "drain_ash_kg_s", "elutriated_ash_kg_s")}), flush=True)
        status = "OK"
    except _Timeout:
        status = "TIMEOUT"
    finally:
        signal.alarm(0)
        residual._Context.domain = check
    (_OUT / "study.json").write_text(json.dumps({"status": status, "runs": runs}, ensure_ascii=False, indent=1), encoding="utf-8")
    return 0 if status == "OK" else 1


if __name__ == "__main__":
    sys.exit(main())
