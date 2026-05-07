#!/usr/bin/env python3
"""LU 工况下 freeboard 反应集合灵敏度审计。"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_validation_case_node,
)


SETS = {
    "default_no_r8": ("R5", "R6", "R7", "R10", "R11", "R12"),
    "legacy_with_r8": ("R5", "R6", "R7", "R8", "R10", "R11", "R12"),
    "no_r7": ("R5", "R6", "R8", "R10", "R11", "R12"),
    "no_r7_r8": ("R5", "R6", "R10", "R11", "R12"),
    "oxidation_only": ("R5", "R6", "R10", "R12"),
    "oxidation_plus_wgs": ("R5", "R6", "R8", "R10", "R12"),
}


def main() -> int:
    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    ref_dry = ref["exit_gas_dry_mol_frac"]
    print("=" * 156)
    print("LU freeboard reaction-set audit")
    print("=" * 156)
    print(
        f"{'set':>18} {'Texit':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'CH4':>8} "
        f"{'n_iter':>8} {'wall[s]':>8} {'|R7_th|':>10} {'|R8_th|':>10} {'|R9_th|':>10} {'|R11_th|':>10}"
    )
    print("-" * 156)

    for name, reactions in SETS.items():
        cfg = build_phase2_htw_lu_freeboard_reactor_config()
        cfg.freeboard_enabled_reactions = reactions
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
            check_x0=True,
        )
        diag_impl = result.get("freeboard_reaction_diag_impl", result["freeboard_reaction_diag"])
        diag = result.get("freeboard_reaction_diag_thesis", diag_impl)
        sum_r7 = sum(abs(v) for v in diag.get("R7", []))
        sum_r8 = sum(abs(v) for v in diag.get("R8", []))
        sum_r9 = sum(abs(v) for v in diag.get("R9", []))
        sum_r11 = sum(abs(v) for v in diag.get("R11", []))
        dry = result["exit_gas_dry"]
        print(
            f"{name:>18s} {result['reactor_exit_T']:>8.1f} {dry['CO']:>8.4f} {dry['CO2']:>8.4f} "
            f"{dry['H2']:>8.4f} {dry['CH4']:>8.4f} {result['n_iter']:>8d} {monitor['wall_time_s']:>8.2f} "
            f"{sum_r7:>10.4f} {sum_r8:>10.4f} {sum_r9:>10.4f} {sum_r11:>10.4f}"
        )
        print(
            f"{'':>18} ref={'':>4} {float(ref_dry['CO']):>8.4f} {float(ref_dry['CO2']):>8.4f} "
            f"{float(ref_dry['H2']):>8.4f} {float(ref_dry['CH4']):>8.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
