"""审计 LU shared NR 口径下的 wall-time 热点与 Jacobian 策略收益。

用法：
    python3 scripts/audit_global_nr_walltime.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


def run_case(strategy: str) -> dict:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    solve_kwargs = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    solve_kwargs["nr_jacobian_strategy"] = strategy

    result, monitor = solve_with_nr_monitor(
        reactor,
        {**solve_kwargs, "verbose": False},
        check_x0=True,
    )
    t_prof = result["T_profile"]
    peak_idx = max(range(len(t_prof)), key=lambda i: t_prof[i])

    return {
        "strategy": strategy,
        "wall_s": float(monitor.get("wall_time_s", float("nan"))),
        "converged": result.get("converged"),
        "nr_init_strategy": result.get("nr_init_strategy"),
        "nr_init_s_total": result.get("nr_init_s_total"),
        "nr_vorabrechnung_s": result.get("nr_vorabrechnung_s"),
        "outer_iters": result.get("nr_outer_iters"),
        "n_iter": result.get("n_iter"),
        "nr_total_s": result.get("nr_total_s"),
        "nr_gs_warmup_s": result.get("nr_gs_warmup_s"),
        "nr_inner_solve_s_total": result.get("nr_inner_solve_s_total"),
        "nr_inner_budget_total": result.get("nr_inner_budget_total"),
        "nr_inner_budget_used": result.get("nr_inner_budget_used"),
        "Texit": float(t_prof[-1]),
        "Tpeak": float(t_prof[peak_idx]),
        "peak_idx": int(peak_idx),
        "carbon_conv": float(result["carbon_conv"]),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "residual": result.get("residual"),
        "nr_timing": result.get("nr_timing") or {},
        "nr_counts": result.get("nr_counts") or {},
        "nr_lambdas": result.get("nr_accepted_lambda_history") or [],
        "nr_line_search_trials": result.get("nr_line_search_trial_counts") or [],
        "nr_clip_history": result.get("nr_clip_history") or [],
        "x0_ok": bool((monitor.get("x0_sanity") or {}).get("ok", False)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit LU global-NR wall time.")
    parser.add_argument(
        "--include-dense",
        action="store_true",
        help="Also run the very slow dense_fd Jacobian path for comparison.",
    )
    args = parser.parse_args()

    strategies = ["block_tridiag_fd"]
    if args.include_dense:
        strategies.insert(0, "dense_fd")

    rows = [run_case(strategy) for strategy in strategies]

    print("=" * 100)
    print("LU global NR wall-time audit")
    print("=" * 100)
    for row in rows:
        timing = row["nr_timing"]
        counts = row["nr_counts"]
        print(f"[{row['strategy']}]")
        print(
            f"  wall_s={row['wall_s']:.2f}s converged={row['converged']} "
            f"outer_iters={row['outer_iters']} n_iter={row['n_iter']}"
        )
        print(f"  x0_ok={row['x0_ok']}")
        print(
            f"  init={row['nr_init_strategy']} "
            f"init_total={row['nr_init_s_total']:.2f}s "
            f"vorabrechnung={row['nr_vorabrechnung_s']:.2f}s "
            f"  solver_total={row['nr_total_s']:.2f}s "
            f"gs_warmup={row['nr_gs_warmup_s']:.2f}s "
            f"inner_solve_total={row['nr_inner_solve_s_total']:.2f}s "
            f"budget={row['nr_inner_budget_used']}/{row['nr_inner_budget_total']}"
        )
        print(
            f"  Texit={row['Texit']:.1f}K Tpeak={row['Tpeak']:.1f}K@cell{row['peak_idx']} "
            f"Xc={row['carbon_conv']:.4f}"
        )
        print(
            f"  rms_scaled_final={row['rms_scaled_final']:.3e} "
            f"residual={float(row['residual']):.3e}"
        )
        print(
            f"  timing: total={timing.get('total_s', float('nan')):.2f}s "
            f"jacobian={timing.get('jacobian_build_s', float('nan')):.2f}s "
            f"linear={timing.get('linear_solve_s', float('nan')):.2f}s "
            f"line_search={timing.get('line_search_residual_s', float('nan')):.2f}s"
        )
        print(
            f"  counts: global_residual_calls={counts.get('global_residual_calls')} "
            f"jacobian_cell_residual_calls={counts.get('jacobian_cell_residual_calls')} "
            f"jacobian_nnz_last={counts.get('jacobian_nnz_last')}"
        )
        lambdas = row["nr_lambdas"]
        trials = row["nr_line_search_trials"]
        clip_hist = row["nr_clip_history"]
        avg_trials = (sum(int(v) for v in trials) / len(trials)) if trials else float("nan")
        max_clip_frac = max((float(item.get("clipped_fraction", 0.0)) for item in clip_hist), default=float("nan"))
        max_temp_clipped = max((int(item.get("temp_clipped_cells", 0)) for item in clip_hist), default=0)
        print(
            f"  damping: avg_trials={avg_trials:.2f} "
            f"accepted_lambdas_tail={lambdas[-5:]} "
            f"max_clip_frac={max_clip_frac:.3f} "
            f"max_temp_clipped_cells={max_temp_clipped}"
        )
        print()

    if len(rows) == 2:
        dense, block = rows
        print("Comparison:")
        print(f"  wall-time speedup = {dense['wall_s'] / max(block['wall_s'], 1e-9):.2f}x")
        print(
            f"  jacobian-build speedup = "
            f"{dense['nr_timing'].get('jacobian_build_s', float('nan')) / max(block['nr_timing'].get('jacobian_build_s', 1e-9), 1e-9):.2f}x"
        )
        print(
            f"  jacobian cell residual call reduction = "
            f"{dense['nr_counts'].get('jacobian_cell_residual_calls', 0)} -> "
            f"{block['nr_counts'].get('jacobian_cell_residual_calls', 0)}"
        )
    else:
        print("Tip: pass --include-dense to compare against the very slow dense_fd baseline.")
    print("=" * 100)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
