"""诊断：初值构造中“固相牛顿 ↔ 预计算”交替的不动点，以及 ADR 0053 两项的闭合。

只报告，不改模型。本工具只运行到初值构造结束，不进入主求解的内层（T26n 卡 C）。

与正式写法只差一处（诊断）：ADR 0014 A5 的余额为负时不停止，排灰取同式的负值；
预计算与固相牛顿的残差中都如此。收敛解上的检查（存量 ≥ −τ、床格惰性物、排灰灰分量
非负、最底床格灰分为正）照常。交替次数上限放宽到 --loops。

输出：--out 指定的 JSON（默认 reports/research-solver-path/alternation.json）；--artifacts 只作
    solve_case 的目录参数（本工具不写内层失败点）。每条交替记录最底格闭合差、两项补入量、
    修正前后的 Ṁ_ab 与床区灰的边界收支差（ADR 0053 第 4 项）。
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


class _InitialDone(BaseException):
    """初值构造已结束：用它跳出 solve_case，不运行主求解的内层。"""


def signed_drain(case, state, released, inert, masses, books, ash_out_redirect=0.0, fuel_out_redirect=0.0):
    """A5 同式（``precalc._drain_terms``），只取消余额的符号拒绝。诊断用。

    直接复用 ``_drain_terms``，不复制公式：本工具与正式写法只差“余额为负时不停止”，
    最底格灰分为正的定义域检查照旧由 ``_drain_terms`` 报错。
    """
    terms = precalc._drain_terms(
        case, state, released, inert, masses, books, ash_out_redirect, fuel_out_redirect,
    )
    return float(terms["ash"]), np.asarray(terms["char"], dtype=float)


def _top_af(books) -> list:
    return [0.0 if book.cells[-1].a_f is None else float(book.cells[-1].a_f) for book in books]


def itemize(frozen, state, context, reaction) -> dict:
    """同一状态上的闭合差、ADR 0053 两项、修正前后的 Ṁ_ab 与床区灰的边界收支差。只报告。"""
    inert, drain_ash, drain_char, down, closure, returned = initial.recompute_flows(state, frozen, context, reaction)
    loss, fuel_exit = precalc.bed_release_terms(context.released, frozen.n_freeboard)
    up = np.asarray(frozen.k_auf_per_s * frozen.geometric_mass_kg, dtype=float)
    masses = np.asarray(frozen.geometric_mass_kg, dtype=float)
    drain = float(drain_ash + np.sum(drain_char))

    def flow(loss_terms, exit_kg_s):
        return precalc.bed_downflow(up, context.case.fuel_kg_s, reaction, returned, drain, loss_terms, exit_kg_s)

    before, without_loss, without_fuel = flow(None, 0.0), flow(None, fuel_exit), flow(loss, 0.0)
    share = np.asarray(inert, dtype=float) / masses
    returning = float(np.sum(np.asarray(frozen.b_free_per_s[-1]) * np.asarray(state.freeboard_ash_free_kg[-1])))
    terms = precalc._drain_terms(
        context.case, state, context.released, inert, masses, context.books,
        context.redirected_ash_out_kg_s, context.redirected_fuel_out_kg_s,
    )
    ash_gen_total = float(np.sum(np.asarray(frozen.ash_gen_bed_kg_s)) + np.sum(frozen.ash_gen_free_kg_s) + np.sum(frozen.ash_gen_hover_kg_s))
    return {
        "closure_kg_s": closure, "closure_before_kg_s": before[1],
        "closure_without_loss_kg_s": without_loss[1], "closure_without_fuel_kg_s": without_fuel[1],
        "down_kg_s": np.asarray(down, dtype=float).tolist(),
        "down_before_kg_s": np.asarray(before[0], dtype=float).tolist(),
        "loss_kg_s": loss.tolist(), "loss_total_kg_s": float(np.sum(loss)), "fuel_exit_kg_s": fuel_exit,
        "feed_kg_s": float(context.case.fuel_kg_s), "returned_kg_s": float(returned),
        "up_kg_s": up.tolist(), "reaction_kg_s": np.asarray(reaction, dtype=float).tolist(),
        "drain_ash_kg_s": float(drain_ash), "drain_char_kg_s": np.asarray(drain_char, dtype=float).tolist(),
        "drain_kg_s": drain, "drain_total_a5_kg_s": float(drain_ash / share[-1]),
        "bottom_mass_fractions": {
            "ash": float(share[-1]),
            "char": float(np.asarray(state.bed_char_kg)[-1].sum() / masses[-1]),
            "fuel": float(np.asarray(frozen.fuel_bed_kg)[-1].sum() / masses[-1]),
        },
        # 床区灰（m_inert = M − 焦炭 − 燃料，不是牛顿未知量）的边界收支差。ADR 0053 第 4 项。
        "bed_ash_gap_kg_s": float(np.sum(np.asarray(frozen.ash_gen_bed_kg_s))) + returning
        - float(up[0] * share[0]) - float(drain_ash),
        "release_mass_gap_bed_kg_s": float(np.sum(np.asarray(context.released.mass_gap_kg_s)[frozen.n_freeboard:])),
        "ash_in_kg_s": float(terms["ash_in"]), "ash_gen_total_kg_s": ash_gen_total,
        "fuel_ash_kg_s": float(terms["fuel_ash"]), "carried_ash_kg_s": float(terms["carried"]),
    }


def run(loops: int, artifacts: str | None, limit_s: float) -> dict:
    log: dict = {"precalc": [], "solids": [], "loops": None}
    current: dict = {}
    original = precalc._drain, precalc._diameters, initial.solve_solids, driver.construct

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
        loss, fuel_exit = precalc.bed_release_terms(released, case.n_free)
        log["precalc"].append(dict(
            current, a5_negative=negative, drain_ash_kg_s=float(out[0]),
            elutriated_ash_kg_s=float(np.dot(_top_af(books), np.asarray(state.freeboard_ash_free_kg[0], dtype=float))),
            loss_total_kg_s=float(np.sum(loss)), fuel_exit_kg_s=fuel_exit,
        ))
        return out

    def solve(state, frozen, context, m_dot_ref, **kwargs):
        try:
            solved = original[2](state, frozen, context, m_dot_ref, **kwargs)
        except Exception as exc:
            log["solids"].append({"error": f"{type(exc).__name__}: {exc}"[:500]})
            raise
        top = np.asarray(solved.state.freeboard_ash_free_kg[0], dtype=float)
        reaction = -np.asarray(initial.evaluate_residual(
            initial.pack_state(solved.state, frozen), frozen,
        ).char_consumption_bed_kg_s, dtype=float).sum(axis=1)
        record = {
            "iterations": solved.iterations, "bed_char_kg": solved.bed_char_kg,
            "bed_char_by_class_kg": np.asarray(solved.state.bed_char_kg).sum(axis=0).tolist(),
            "drain_ash_kg_s": solved.drain_ash_kg_s, "top_free_ash_kg": top.tolist(),
            "a_f_top_frozen_per_s": np.asarray(frozen.a_free_per_s[0]).tolist(),
            "elutriated_ash_kg_s": float(np.dot(frozen.a_free_per_s[0], top)),
            "bed_inert_min_kg": float(np.min(solved.bed_inert_kg)),
        }
        record.update(itemize(frozen, solved.state, context, reaction))
        log["solids"].append(record)
        return solved

    def construct(*args, **kwargs):
        kwargs["max_loops"] = loops
        log["loops"] = kwargs.get("records")
        out = original[3](*args, **kwargs)
        log["initial"] = {"note": out[2]}
        raise _InitialDone()

    precalc._drain, precalc._diameters, initial.solve_solids, driver.construct = (
        precalc_drain, diameters, solve, construct,
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
    except _InitialDone:
        log["status"] = ["INITIAL_ONLY", "初值构造结束；未运行主求解的内层"]
    except _Timeout:
        log["status"] = ["TIMEOUT", ""]
    except Exception as exc:
        log["status"] = [type(exc).__name__, str(exc)[:500]]
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
        precalc._drain, precalc._diameters, initial.solve_solids, driver.construct = original
    log["elapsed_s"] = time.monotonic() - started
    return log


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loops", type=int, default=10)
    parser.add_argument("--artifacts", default=None)
    parser.add_argument("--limit-s", type=float, default=5400.0)
    parser.add_argument("--out", default=str(_OUT))
    args = parser.parse_args(argv)
    log = run(args.loops, args.artifacts, args.limit_s)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(log, ensure_ascii=False, indent=1, default=lambda o: np.asarray(o).tolist()), encoding="utf-8")
    print(log["status"], f"{log['elapsed_s']:.0f} s")


if __name__ == "__main__":
    main()
