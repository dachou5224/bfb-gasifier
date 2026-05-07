"""审计 LU global NR 在不同初始化策略下的 wall-time 与分支差异。

用法：
    python3 scripts/audit_global_nr_init_strategies.py
"""

from __future__ import annotations

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


def run_case(init_strategy: str, warmup_steps: int | None = None) -> dict:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    solve_kwargs = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    solve_kwargs["nr_init_strategy"] = init_strategy
    if warmup_steps is not None:
        solve_kwargs["nr_gs_warmup_steps"] = warmup_steps
    else:
        solve_kwargs.pop("nr_gs_warmup_steps", None)

    result, monitor = solve_with_nr_monitor(
        reactor,
        {**solve_kwargs, "verbose": False},
        check_x0=True,
    )
    t_prof = result["T_profile"]
    peak_idx = max(range(len(t_prof)), key=lambda i: t_prof[i])

    return {
        "label": (
            init_strategy
            if warmup_steps is None
            else f"{init_strategy}(steps={warmup_steps})"
        ),
        "nr_init_strategy": result.get("nr_init_strategy"),
        "nr_gs_warmup_steps": result.get("nr_gs_warmup_steps"),
        "wall_s": float(monitor.get("wall_time_s", float("nan"))),
        "nr_init_s_total": result.get("nr_init_s_total"),
        "nr_vorabrechnung_s": result.get("nr_vorabrechnung_s"),
        "nr_gs_warmup_s": result.get("nr_gs_warmup_s"),
        "nr_total_s": result.get("nr_total_s"),
        "outer_iters": result.get("nr_outer_iters"),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "residual": result.get("residual"),
        "Texit": float(t_prof[-1]),
        "Tpeak": float(t_prof[peak_idx]),
        "peak_idx": int(peak_idx),
        "carbon_conv": float(result["carbon_conv"]),
        "x0_ok": bool((monitor.get("x0_sanity") or {}).get("ok", False)),
    }


def main() -> int:
    rows = [
        run_case("vorabrechnung"),
        run_case("gs_warmup", warmup_steps=1),
        run_case("gs_warmup", warmup_steps=2),
    ]

    print("=" * 100)
    print("LU global NR init-strategy audit")
    print("=" * 100)
    for row in rows:
        print(f"[{row['label']}]")
        print(
            f"  init={row['nr_init_strategy']} warmup_steps={row['nr_gs_warmup_steps']} "
            f"wall_s={row['wall_s']:.2f}s solver_total={row['nr_total_s']:.2f}s"
        )
        print(f"  x0_ok={row['x0_ok']}")
        print(
            f"  init_total={row['nr_init_s_total']:.2f}s "
            f"vorabrechnung={row['nr_vorabrechnung_s']:.2f}s "
            f"gs_warmup={row['nr_gs_warmup_s']:.2f}s"
        )
        print(
            f"  outer_iters={row['outer_iters']} n_iter={row['n_iter']} "
            f"rms_scaled_final={row['rms_scaled_final']:.3e} residual={float(row['residual']):.3e}"
        )
        print(
            f"  Texit={row['Texit']:.1f}K Tpeak={row['Tpeak']:.1f}K@cell{row['peak_idx']} "
            f"Xc={row['carbon_conv']:.4f}"
        )
        print()

    fastest = min(rows, key=lambda row: row["wall_s"])
    coolest_exit = min(rows, key=lambda row: row["Texit"])
    lowest_rms = min(rows, key=lambda row: row["rms_scaled_final"])
    print("Summary:")
    print(f"  fastest = {fastest['label']} ({fastest['wall_s']:.2f}s)")
    print(f"  coolest Texit = {coolest_exit['label']} ({coolest_exit['Texit']:.1f}K)")
    print(f"  lowest rms = {lowest_rms['label']} ({lowest_rms['rms_scaled_final']:.3e})")
    print("=" * 100)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
