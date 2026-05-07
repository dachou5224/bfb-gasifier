#!/usr/bin/env python3
"""对比 global-NR 初始化策略：`vorabrechnung` vs `gs_warmup`。

默认基于 HTW LU shared global-NR 配置：
- config: tests.validation_case_utils.build_phase1_htw_lu_global_nr_reactor_config
- solve kwargs: tests.validation_case_utils.PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS
"""

from __future__ import annotations

import argparse
import copy
import json
import multiprocessing as mp
from pathlib import Path
from time import perf_counter
from typing import Any
import sys

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


def _run_once(*, init_strategy: str, solve_kwargs: dict[str, Any], gs_warmup_steps: int) -> dict[str, Any]:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)

    reactor = Reactor(cfg)
    kwargs = dict(solve_kwargs)
    kwargs["nr_init_strategy"] = str(init_strategy)
    if str(init_strategy) == "gs_warmup":
        kwargs["nr_gs_warmup_steps"] = int(max(gs_warmup_steps, 0))
    t0 = perf_counter()
    result, monitor = solve_with_nr_monitor(
        reactor,
        kwargs,
        check_x0=True,
    )
    wall_s = perf_counter() - t0

    lambdas = result.get("nr_accepted_lambda_history", []) or []
    n_lambda_none = sum(1 for v in lambdas if v is None)
    n_lambda_ok = sum(1 for v in lambdas if v is not None)

    return {
        "nr_init_strategy": str(init_strategy),
        "wall_s": float(wall_s),
        "monitor_wall_s": float(monitor.get("wall_time_s", 0.0) or 0.0),
        "x0_ok": bool((monitor.get("x0_sanity") or {}).get("ok", False)),
        "converged": bool(result.get("converged", False)),
        "converged_outer": bool(result.get("converged_outer", False)),
        "converged_inner_nr": bool(result.get("converged_inner_nr", False)),
        "converged_fully": bool(result.get("converged_fully", False)),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("inf"))),
        "residual": float(result.get("residual", float("inf"))),
        "nr_total_s": float(result.get("nr_total_s", 0.0) or 0.0),
        "nr_init_s_total": float(result.get("nr_init_s_total", 0.0) or 0.0),
        "nr_vorabrechnung_s": float(result.get("nr_vorabrechnung_s", 0.0) or 0.0),
        "nr_outer_iters": int(result.get("nr_outer_iters", 0) or 0),
        "nr_inner_n_iter_total": int(result.get("n_iter", 0) or 0),
        "nr_line_search_eval_total": int((result.get("nr_counts") or {}).get("line_search_evaluations", 0)),
        "nr_line_search_backtracks_total": int((result.get("nr_counts") or {}).get("line_search_backtracks", 0)),
        "nr_lambda_accept_count": int(n_lambda_ok),
        "nr_lambda_reject_count": int(n_lambda_none),
        "T_exit_K": float(result.get("T_exit", 0.0) or 0.0),
        "carbon_conv": float(result.get("carbon_conv", 0.0) or 0.0),
    }


def _run_once_worker(queue: mp.Queue, init_strategy: str, solve_kwargs: dict[str, Any], gs_warmup_steps: int) -> None:
    try:
        queue.put(_run_once(init_strategy=init_strategy, solve_kwargs=solve_kwargs, gs_warmup_steps=gs_warmup_steps))
    except Exception as exc:  # pragma: no cover - benchmark robustness path
        queue.put({"nr_init_strategy": init_strategy, "error": repr(exc)})


def _run_once_with_timeout(
    *,
    init_strategy: str,
    solve_kwargs: dict[str, Any],
    gs_warmup_steps: int,
    timeout_s: int,
) -> dict[str, Any]:
    q: mp.Queue = mp.Queue(maxsize=1)
    p = mp.Process(
        target=_run_once_worker,
        args=(q, init_strategy, solve_kwargs, gs_warmup_steps),
        daemon=True,
    )
    t0 = perf_counter()
    p.start()
    p.join(timeout=float(max(timeout_s, 1)))
    elapsed = perf_counter() - t0
    if p.is_alive():
        p.terminate()
        p.join(2.0)
        return {
            "nr_init_strategy": str(init_strategy),
            "timeout": True,
            "wall_s": float(elapsed),
        }
    if not q.empty():
        out = q.get_nowait()
        if isinstance(out, dict):
            out.setdefault("timeout", False)
            out.setdefault("wall_s", float(elapsed))
            return out
    return {
        "nr_init_strategy": str(init_strategy),
        "error": "no_result",
        "wall_s": float(elapsed),
        "timeout": False,
    }


def _mean_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    if not rows:
        return {}
    ok_rows = [r for r in rows if not r.get("timeout", False) and "error" not in r]
    if not ok_rows:
        return {
            "n_runs": len(rows),
            "ok_runs": 0,
            "timeout_count": int(sum(1 for r in rows if r.get("timeout", False))),
            "error_count": int(sum(1 for r in rows if "error" in r)),
        }
    keys_float = [
        "wall_s",
        "rms_scaled_final",
        "residual",
        "nr_total_s",
        "nr_init_s_total",
        "nr_vorabrechnung_s",
        "T_exit_K",
        "carbon_conv",
    ]
    keys_int = [
        "nr_outer_iters",
        "nr_inner_n_iter_total",
        "nr_line_search_eval_total",
        "nr_line_search_backtracks_total",
        "nr_lambda_accept_count",
        "nr_lambda_reject_count",
    ]
    out: dict[str, Any] = {
        "n_runs": len(rows),
        "ok_runs": len(ok_rows),
        "timeout_count": int(sum(1 for r in rows if r.get("timeout", False))),
        "error_count": int(sum(1 for r in rows if "error" in r)),
        "converged_count": int(sum(1 for r in ok_rows if r.get("converged"))),
        "converged_fully_count": int(sum(1 for r in ok_rows if r.get("converged_fully"))),
    }
    for k in keys_float:
        out[f"{k}_mean"] = float(sum(float(r[k]) for r in ok_rows) / len(ok_rows))
    for k in keys_int:
        out[f"{k}_mean"] = float(sum(int(r[k]) for r in ok_rows) / len(ok_rows))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark NR init strategy: vorabrechnung vs gs_warmup")
    parser.add_argument("--repeat", type=int, default=1, help="number of runs per setting")
    parser.add_argument("--max-global-iter", type=int, default=20, help="override solve max_global_iter")
    parser.add_argument("--tol-global", type=float, default=1.0, help="override solve tol_global")
    parser.add_argument("--gs-warmup-steps", type=int, default=1, help="warmup steps when strategy=gs_warmup")
    parser.add_argument("--timeout-s", type=int, default=300, help="timeout per run/strategy in seconds")
    parser.add_argument("--output-json", type=Path, default=None, help="optional path to dump benchmark JSON")
    args = parser.parse_args()

    solve_kwargs = copy.deepcopy(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    solve_kwargs["max_global_iter"] = int(max(args.max_global_iter, 1))
    solve_kwargs["tol_global"] = float(args.tol_global)

    rows_vorab: list[dict[str, Any]] = []
    rows_gs: list[dict[str, Any]] = []
    for _ in range(int(max(args.repeat, 1))):
        rows_vorab.append(
            _run_once_with_timeout(
                init_strategy="vorabrechnung",
                solve_kwargs=solve_kwargs,
                gs_warmup_steps=0,
                timeout_s=int(max(args.timeout_s, 1)),
            )
        )
        rows_gs.append(
            _run_once_with_timeout(
                init_strategy="gs_warmup",
                solve_kwargs=solve_kwargs,
                gs_warmup_steps=int(max(args.gs_warmup_steps, 0)),
                timeout_s=int(max(args.timeout_s, 1)),
            )
        )

    m_vorab = _mean_metrics(rows_vorab)
    m_gs = _mean_metrics(rows_gs)

    if "wall_s_mean" in m_vorab and "wall_s_mean" in m_gs:
        speedup = (m_gs["wall_s_mean"] - m_vorab["wall_s_mean"]) / max(m_gs["wall_s_mean"], 1e-12)
    else:
        speedup = 0.0

    print("=== NR Init Strategy Benchmark (LU) ===")
    print(f"repeat={args.repeat}, solve_kwargs={solve_kwargs}")
    print("")
    print("[A] nr_init_strategy=vorabrechnung")
    print(json.dumps(m_vorab, ensure_ascii=False, indent=2))
    print("")
    print(f"[B] nr_init_strategy=gs_warmup (steps={int(max(args.gs_warmup_steps, 0))})")
    print(json.dumps(m_gs, ensure_ascii=False, indent=2))
    print("")
    print(f"wall-time speedup (vorabrechnung vs gs_warmup): {speedup*100.0:.2f}%")

    payload = {
        "solve_kwargs": solve_kwargs,
        "gs_warmup_steps": int(max(args.gs_warmup_steps, 0)),
        "vorabrechnung_runs": rows_vorab,
        "gs_warmup_runs": rows_gs,
        "vorabrechnung_mean": m_vorab,
        "gs_warmup_mean": m_gs,
        "wall_speedup_fraction": float(speedup),
    }
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"saved: {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
