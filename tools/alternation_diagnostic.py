"""诊断：初值构造中“固相牛顿 ↔ 预计算”交替的不动点，以及其后的第 1 轮内层。只报告，不改模型。

须在 T26m（PR #141，提交 e454bde）的代码上运行：本工具依赖 ADR 0052 修订后的试探点可行条件。

与正式写法只差一处（诊断）：ADR 0014 A5 的余额为负时不停止，排灰取同式的负值；
预计算与固相牛顿的残差中都如此。收敛解上的检查（存量 ≥ −τ、床格惰性物、排灰灰分量
非负、最底床格灰分为正）照常。交替次数上限放宽到 --loops。

输出：reports/research-solver-path/alternation.json；--artifacts 给出时，内层失败点的
状态、残差与雅可比矩阵写入该目录（failure.npz、frozen.pkl，只留本地，不入库）。
"""

from __future__ import annotations

import argparse
import json
import signal
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.driver as driver
import hamel.solver.initial as initial
import hamel.solver.precalc as precalc
from hamel.cases.config import SolverConfig

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "research-solver-path" / "alternation.json"


class _Timeout(BaseException):
    """继承 BaseException：任何 except Exception 都不会吞掉全局期限（CHARTER §7）。"""


def _ash_in(case) -> float:
    return case.fuel_kg_s * float(case.physical.ash_mass_fraction_wf) * (1.0 - float(case.physical.moisture_mass_fraction))


def _top_af(books) -> list:
    return [0.0 if book.cells[-1].a_f is None else float(book.cells[-1].a_f) for book in books]


def signed_drain(case, state, released, inert, masses, books, ash_out_redirect=0.0, fuel_out_redirect=0.0):
    """A5 同式（precalc._drain），余额为负时不停止。诊断用。"""
    ash_in = _ash_in(case)
    carried = float(np.dot(np.asarray(state.freeboard_ash_free_kg[0], dtype=float), _top_af(books)))
    fuel_top = 0.0 if released.fuel_up_kg_s is None else float(released.fuel_up_kg_s[0].sum())
    n = case.n_free
    fuel_drain = float(released.undried_exit_kg_s[n:].sum() + released.dried_exit_kg_s[n:].sum())
    ash_out = carried + ash_in / case.fuel_kg_s * (fuel_top + fuel_drain + float(fuel_out_redirect)) + float(ash_out_redirect)
    fraction = float(inert[-1] / masses[-1])
    if fraction <= 0.0:
        raise precalc.PrecalcError("最底格灰分质量分数不是正数")
    total = (ash_in - ash_out) / fraction
    return ash_in - ash_out, total * np.asarray(state.bed_char_kg[-1], dtype=float) / float(masses[-1])


def run(loops: int, artifacts: str | None, limit_s: float) -> dict:
    log: dict = {"precalc": [], "solids": [], "loops": None}
    current: dict = {}
    original = precalc._drain, precalc._diameters, initial.solve_solids, initial._drain, driver.construct

    def diameters(diameter_m, fraction, phi_s):
        out = original[1](diameter_m, fraction, phi_s)
        current.update(d_bar_m=out[0], d_sv_m=out[1])
        return out

    def guarded(*args):
        try:
            return original[0](*args), False
        except precalc.PrecalcError as exc:
            if "带出灰" not in str(exc):
                raise
            return signed_drain(*args), True

    def precalc_drain(case, state, released, inert, masses, books, *rest):
        out, negative = guarded(case, state, released, inert, masses, books, *rest)
        a_f = _top_af(books)
        top = np.asarray(state.freeboard_ash_free_kg[0], dtype=float)
        log["precalc"].append(dict(
            current, ash_in_kg_s=_ash_in(case), a_f_top_per_s=a_f, top_free_ash_kg=top.tolist(),
            elutriated_ash_kg_s=float(np.dot(top, a_f)), drain_ash_kg_s=float(out[0]), a5_negative=negative,
        ))
        return out

    def solids_drain(*args):
        return guarded(*args)[0]

    def solve(state, frozen, context, m_dot_ref, **kwargs):
        try:
            solved = original[2](state, frozen, context, m_dot_ref, **kwargs)
        except Exception as exc:
            log["solids"].append({"error": f"{type(exc).__name__}: {exc}"[:500]})
            raise
        top = np.asarray(solved.state.freeboard_ash_free_kg[0], dtype=float)
        log["solids"].append({
            "iterations": solved.iterations, "bed_char_kg": solved.bed_char_kg,
            "bed_char_by_class_kg": np.asarray(solved.state.bed_char_kg).sum(axis=0).tolist(),
            "drain_ash_kg_s": solved.drain_ash_kg_s, "top_free_ash_kg": top.tolist(),
            "a_f_top_frozen_per_s": np.asarray(frozen.a_free_per_s[0]).tolist(),
            "elutriated_ash_kg_s": float(np.dot(frozen.a_free_per_s[0], top)),
            "bed_inert_min_kg": float(np.min(solved.bed_inert_kg)),
        })
        return solved

    def construct(*args, **kwargs):
        kwargs["max_loops"] = loops
        log["loops"] = kwargs.get("records")
        return original[4](*args, **kwargs)

    precalc._drain, precalc._diameters, initial.solve_solids, initial._drain, driver.construct = (
        precalc_drain, diameters, solve, solids_drain, construct,
    )

    def expire(signum, frame):
        raise _Timeout()

    previous = signal.signal(signal.SIGALRM, expire)
    signal.alarm(int(limit_s) + 1)
    started = time.monotonic()
    try:
        result = driver.solve_case(
            str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=1), artifact_dir=artifacts,
            deadline=lambda: time.monotonic() - started >= limit_s,
        )
        log["status"] = [result.status, result.reason]
        rounds = (result.diagnostics or {}).get("rounds") or []
        log["inner_round_1"] = rounds[0].get("newton") if rounds else None
    except _Timeout:
        log["status"] = ["TIMEOUT", ""]
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
        precalc._drain, precalc._diameters, initial.solve_solids, initial._drain, driver.construct = original
    log["elapsed_s"] = time.monotonic() - started
    return log


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loops", type=int, default=10)
    parser.add_argument("--artifacts", default=None)
    parser.add_argument("--limit-s", type=float, default=5400.0)
    args = parser.parse_args(argv)
    log = run(args.loops, args.artifacts, args.limit_s)
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    _OUT.write_text(json.dumps(log, ensure_ascii=False, indent=1, default=lambda o: np.asarray(o).tolist()), encoding="utf-8")
    print(log["status"], f"{log['elapsed_s']:.0f} s")


if __name__ == "__main__":
    main()
