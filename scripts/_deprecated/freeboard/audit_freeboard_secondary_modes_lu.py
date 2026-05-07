#!/usr/bin/env python3
"""比较 secondary inlet 的 lumped vs distributed_uniform 模式。"""

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
    apply_secondary_air_repartition,
    build_phase2_htw_lu_freeboard_reactor_config,
)


MODES = ("lumped", "distributed_uniform")
SEC_AIR_FRAC = 0.05
REFINE = 8


def main() -> int:
    print("=" * 176)
    print("LU freeboard secondary injection mode audit")
    print("=" * 176)
    for mode in MODES:
        cfg = build_phase2_htw_lu_freeboard_reactor_config()
        cfg.freeboard_secondary_local_refine = REFINE
        cfg.freeboard_secondary_injection_mode = mode
        apply_secondary_air_repartition(
            cfg,
            secondary_air_frac=SEC_AIR_FRAC,
            injection_xi=cfg.freeboard_secondary_injection_xi,
            secondary_T_K=cfg.T_inlet,
        )
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
            check_x0=True,
        )
        xi = result["freeboard_axial_xi"]
        o2 = result["freeboard_gas_profiles_wet"]["O2"]
        temp = result["freeboard_T_profile"]
        co = result["freeboard_gas_profiles_wet"]["CO"]
        h2 = result["freeboard_gas_profiles_wet"]["H2"]
        window = [
            (x, o, t, c, h)
            for x, o, t, c, h in zip(xi, o2, temp, co, h2)
            if 0.55 <= x <= 0.75
        ]
        print("-" * 176)
        print(
            f"mode={mode} refine={REFINE} sec_air_frac={SEC_AIR_FRAC:.2f} "
            f"Texit={result['reactor_exit_T']:.1f}K maxO2_window={max(o for _, o, _, _, _ in window):.5f}"
        )
        print_nr_monitor(monitor, prefix="NR monitor")
        print(f"{'xi':>8} {'O2':>8} {'T[K]':>8} {'CO':>8} {'H2':>8}")
        for x, o, t, c, h in window:
            print(f"{x:>8.3f} {o:>8.4f} {t:>8.1f} {c:>8.4f} {h:>8.4f}")
    print("-" * 176)
    print("note:")
    print("  - distributed_uniform 只是把同一总 secondary air 分摊到局部细分 slice。")
    print("  - 若它能抬出非零 O2 局部峰，则问题主要在当前 lumped injection operator。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
