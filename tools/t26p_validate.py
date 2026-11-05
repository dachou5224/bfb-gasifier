"""T26p 验证：ADR 0055 的逐轮灰平衡（工况 1、φ = 1）。只报告，不改 ``hamel/``。

Corella 外推只在本工具里设置（``_CORELLA_HIGH_K = +inf``），不进入正式路径。两个子命令：
``extrapolated`` 从头跑外层（最多 ``--rounds`` 轮，默认 outer_round_limit）；``frozen`` 不外推，
只跑第 1 轮。逐轮记录 w_f、E、E₀、E₁、评价次数、T、Ṁ_A、子求解迭代、E(0,5) 与 (E₀+E₁)/2 的相对差、
内层状态与迭代数、各块归一残差、最底床格与最高床温、床区焦炭存量、第 5 项三条检查与 (b) 复核。全程硬超时。
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
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.solver.gates import block_ratios
from hamel.solver.newton import SolveTimeout
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "t26p-validation"
_REFS = None
class _Timeout(BaseException):
    """继承 BaseException：任何 except Exception 都不会吞掉全局期限（CHARTER §7）。"""


def _row(index, state, frozen, monitor, prec) -> dict:
    global _REFS
    if _REFS is None:
        _REFS = driver._references(assemble_case(str(_CASE), "1").physical)
    ratios = block_ratios(evaluate_residual(pack_state(state, frozen), frozen).values,
                          driver._blocks(frozen, *_REFS))
    row = {
        "round": index, "w_f": float(frozen.bed_ash_fine_share), "phi": float(frozen.fly_ash_fraction),
        "bed_char_kg": float(np.sum(state.bed_char_kg)),
        "bottom_bed_c": float(state.bed_temperature_k[-1]) - 273.15,
        "top_freeboard_c": float(state.freeboard_temperature_k[0]) - 273.15,
        "drain_ash_kg_s": float(frozen.drain_ash_kg_s), "drain_char_kg_s": float(np.sum(frozen.drain_char_kg_s)),
        "blocks": ratios, "ash_relative_open": monitor.ash_relative,
    }
    balance = (prec.reported or {}).get("ash_balance")
    if balance:
        row["ash_balance"] = {key: value.tolist() if isinstance(value, np.ndarray) else value
                              for key, value in balance.items()}
        target = float(balance["target"])
        row["checks"] = {
            "coarse": (float(balance["E0"]) - target) / target,
            "fine": (target - float(balance["E1"])) / target,
            "closed": abs(float(balance["E"]) - target) / target,
            "evaluations": balance["evaluations"],
        }
    return row


def run(extrapolate: bool, rounds: int, limit_s: float) -> dict:
    if extrapolate:
        import hamel.kinetics.homogeneous as homogeneous
        homogeneous._CORELLA_HIGH_K = float("inf")
    box: dict = {"pre": [], "rows": [], "initial_precalcs": None}
    original_pre, original_monitor = driver.precalculate, driver._monitor

    def precalculate(*args, **kwargs):
        result = original_pre(*args, **kwargs)
        box["pre"].append(result)
        return result

    def monitor(state, frozen, assembly=None):
        if box["initial_precalcs"] is None:
            box["initial_precalcs"] = len(box["pre"])
        value = original_monitor(state, frozen, assembly)
        box["rows"].append(_row(len(box["rows"]) + 1, state, frozen, value, box["pre"][-1]))
        return value

    driver.precalculate, driver._monitor = precalculate, monitor
    started = time.monotonic()

    def expire(signum, frame):
        raise _Timeout()

    signal.signal(signal.SIGALRM, expire)
    signal.setitimer(signal.ITIMER_REAL, float(limit_s))
    summary: dict = {}
    try:
        result = driver.solve_case(
            str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=int(rounds)),
            deadline=lambda: time.monotonic() - started >= limit_s,
        )
        summary.update(status=result.status, reason=result.reason, rounds=result.rounds,
                       iterations=result.iterations, balance=result.balance,
                       consistent_ash_relative=result.diagnostics.get("consistent_ash_relative"))
        summary["outer"] = [{"round": item.get("round"), "inner": item.get("inner"), "reason": item.get("reason"),
                             "ash_relative": item.get("ash_relative")}
                            for item in result.diagnostics.get("rounds") or []]
        initial = result.diagnostics.get("initial") or {}
        loops = initial.get("loops") or []
        summary["initial"] = {"loops": len(loops), "note": initial.get("note"),
                              "last": dict(loops[-1]) if loops else None}
    except (_Timeout, SolveTimeout) as exc:
        summary.update(status="TIMEOUT", reason=str(exc) or "运行超过墙钟上限")
    except Exception as exc:  # noqa: BLE001 - 运行记录要留下任何异常
        summary.update(status=type(exc).__name__, reason=str(exc)[:800])
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        driver.precalculate, driver._monitor = original_pre, original_monitor
    summary["elapsed_s"] = time.monotonic() - started
    summary["initial_precalcs"] = box["initial_precalcs"]
    return {"summary": summary, "rounds": box["rows"], "extrapolated": extrapolate}


def _plain(value):
    if isinstance(value, (np.floating, np.integer)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return str(value)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit-s", type=float, default=4.0 * 3600.0)
    parser.add_argument("--rounds", type=int, default=SolverConfig().outer_round_limit)
    parser.add_argument("--out", default=str(_OUT))
    parser.add_argument("--only", choices=("frozen", "extrapolated"))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for name, extrapolate, rounds in (("frozen", False, 1), ("extrapolated", True, args.rounds)):
        if args.only is not None and args.only != name:
            continue
        payload = run(extrapolate, rounds, args.limit_s)
        (out / f"{name}.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=1, default=_plain), encoding="utf-8",
        )
        print(name, payload["summary"]["status"], payload["summary"].get("reason"), flush=True)


if __name__ == "__main__":
    main()
