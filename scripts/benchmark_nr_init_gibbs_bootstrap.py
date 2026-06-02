#!/usr/bin/env python3
"""对比 global-NR 初始化策略（NR-only：vorabrechnung smoke benchmark）。"""

from __future__ import annotations

import argparse
import copy
import json
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


def _run_once(*, solve_kwargs: dict[str, Any]) -> dict[str, Any]:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    reactor = Reactor(cfg)
    kwargs = dict(solve_kwargs)
    kwargs["nr_init_strategy"] = "vorabrechnung"
    t0 = perf_counter()
    result, monitor = solve_with_nr_monitor(reactor, kwargs, check_x0=True)
    wall_s = perf_counter() - t0
    return {
        "nr_init_strategy": "vorabrechnung",
        "wall_s": float(wall_s),
        "x0_ok": bool((monitor.get("x0_sanity") or {}).get("ok", False)),
        "converged": bool(result.get("converged", False)),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("inf"))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Smoke benchmark: NR init via vorabrechnung")
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--max-global-iter", type=int, default=20)
    parser.add_argument("--tol-global", type=float, default=1.0)
    parser.add_argument("--output-json", type=Path, default=None)
    args = parser.parse_args()

    solve_kwargs = copy.deepcopy(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    solve_kwargs["max_global_iter"] = int(max(args.max_global_iter, 1))
    solve_kwargs["tol_global"] = float(args.tol_global)

    rows = [_run_once(solve_kwargs=solve_kwargs) for _ in range(int(max(args.repeat, 1)))]
    print(json.dumps({"runs": rows, "solve_kwargs": solve_kwargs}, ensure_ascii=False, indent=2))
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
