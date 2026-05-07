#!/usr/bin/env python3
"""Audit bed-to-freeboard solid-fraction transition on current shared LU NR chain."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.physics.phase_fractions import calc_bulk_solid_holdup
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_case_LU,
)


def main() -> int:
    case = load_case_LU()
    cfg = build_phase2_htw_lu_freeboard_reactor_config(case)
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    print("=" * 156)
    print("LU bed solid-transition audit")
    print("=" * 156)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"converged={result.get('converged')} n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"bed_exit_T={float(result['bed_T_profile'][-1]):.1f}K reactor_exit_T={float(result['reactor_exit_T']):.1f}K"
    )
    print("-" * 156)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'d_b':>8} {'u_b':>8} {'eps_b':>8} "
        f"{'dense_share':>11} {'eps_d_void':>11} {'eps_avg':>9} {'dense_solid':>11} {'bulk_solid':>11} {'K_bd':>8}"
    )
    print("-" * 156)

    bulk_solid = []
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        dense_share = 1.0 - float(cell.eps_b)
        eps_avg = float(cell.eps_b + dense_share * float(cell.eps_d_voidage))
        dense_solid = 1.0 - float(cell.eps_d_voidage)
        bulk = calc_bulk_solid_holdup(float(cell.eps_b), float(cell.eps_d_voidage))
        bulk_solid.append(bulk)
        print(
            f"{i:>4d} {xi:>6.2f} {cell.T:>8.1f} {cell.d_b:>8.4f} {cell.u_b:>8.4f} {cell.eps_b:>8.4f} "
            f"{dense_share:>11.4f} {cell.eps_d_voidage:>11.4f} {eps_avg:>9.4f} {dense_solid:>11.4f} {bulk:>11.4f} {cell.K_bd:>8.4f}"
        )

    trend = bulk_solid[-1] - bulk_solid[0]
    print("-" * 156)
    print(
        f"bulk-solid summary: bottom={bulk_solid[0]:.4f} mid={bulk_solid[len(bulk_solid)//2]:.4f} top={bulk_solid[-1]:.4f} delta_top-bottom={trend:+.4f}"
    )
    print("readout:")
    print("  - `dense_share = 1 - eps_b` is the bed-level fraction occupied by the suspension phase.")
    print("  - `eps_avg = eps_b + (1 - eps_b) * eps_d_voidage` is the whole-bed average voidage.")
    print("  - `dense_solid = 1 - eps_d_voidage` is the solid fraction within that suspension phase.")
    print("  - `bulk_solid = (1 - eps_b) * (1 - eps_d_voidage)` is the whole-bed solid holdup proxy.")
    print("  - a physical bed->freeboard transition should show `bulk_solid` declining toward the top, but not collapsing to ~0 through the whole bed.")
    if trend > 0.0:
        print("  - current result is still suspicious: bulk_solid increases toward the bed top instead of declining.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
