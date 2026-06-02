#!/usr/bin/env python3
"""LU 工况下 freeboard 净方向审计。

目的：
1. 对比 bed exit 与 reactor exit 的净变化
2. 统计 freeboard 沿程主要反应积分
3. 明确当前 freeboard 主导的是哪条 chemistry direction
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
)


TRACK_SPECIES = ("CO", "CO2", "H2", "H2O", "CH4", "TAR1", "TAR2", "O2")
TRACK_REACTIONS_THESIS = ("R5", "R6", "R7", "R8", "R9", "R10", "R11")


def main() -> int:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    bed_wet = result["bed_exit_gas"]
    exit_wet = result["exit_gas"]
    diag_impl = result.get("freeboard_reaction_diag_impl", result["freeboard_reaction_diag"])
    diag = result.get("freeboard_reaction_diag_thesis", diag_impl)
    wet_profiles = result["freeboard_gas_profiles_wet"]

    print("=" * 156)
    print("LU freeboard direction audit")
    print("=" * 156)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"trajectory={result['freeboard_trajectory_model']} "
        f"Texit={result['reactor_exit_T']:.1f}K "
        f"carry_end={result['freeboard_carry_ratio_profile'][-1]:.3f} "
        f"recycle={result['freeboard_cyclone_recycle_candidate_char_ash_kg_s']:.4f} kg/s"
    )
    numbering = result.get("reaction_numbering", {})
    if numbering:
        print(
            f"reaction_label_authority={numbering.get('authority', 'unknown')} | "
            "table=freeboard_reaction_diag_thesis (impl view kept in *_impl)"
        )
    print("-" * 156)
    print(f"{'species':>8} {'bed_exit_wet':>14} {'reactor_exit_wet':>18} {'delta':>12}")
    print("-" * 156)
    for sp in TRACK_SPECIES:
        bed_val = float(bed_wet.get(sp, 0.0))
        exit_val = float(exit_wet.get(sp, 0.0))
        print(f"{sp:>8s} {bed_val:>14.6f} {exit_val:>18.6f} {exit_val - bed_val:>12.6f}")

    print("-" * 156)
    print(f"{'rxn':>8} {'sum':>14} {'max_cell':>14} {'sign_hint':>12}")
    print("-" * 156)
    for rxn in TRACK_REACTIONS_THESIS:
        vals = [float(v) for v in diag.get(rxn, [])]
        s = sum(vals)
        vmax = max(vals, key=abs) if vals else 0.0
        sign = "forward" if s > 0 else ("reverse" if s < 0 else "inactive")
        print(f"{rxn:>8s} {s:>14.6f} {vmax:>14.6f} {sign:>12s}")

    print("-" * 156)
    print("profile tail:")
    for sp in ("CO", "CO2", "H2", "CH4", "TAR1", "TAR2", "O2"):
        prof = wet_profiles.get(sp, [])
        if prof:
            print(
                f"  {sp}: first={prof[0]:.6f} mid={prof[len(prof)//2]:.6f} last={prof[-1]:.6f}"
            )

    print("-" * 156)
    print("interpretation:")
    print("  - 若 R11≈0 且 TAR1/TAR2 全程接近 0，则当前 freeboard 不是 tar-limited chemistry。")
    print("  - 若 R8 显著改变 CO/CO2/H2，但 source-of-truth 仍指向催化 WGSR，则默认集应保持关闭，仅作显式审计项。")
    print("  - 若 R7 是主要活跃项而 CH4 仍接近 0，则 CH4 问题仍主要来自 bed 出口进入 freeboard 前的来源。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
