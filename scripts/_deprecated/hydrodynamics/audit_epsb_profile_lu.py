#!/usr/bin/env python3
"""Audit Eq.3.24 numerator/denominator along the LU bed on current shared NR chain."""

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
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    print("=" * 160)
    print("LU Eq.3.24 epsilon_b profile audit")
    print("=" * 160)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"converged={result.get('converged')} n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"bed_exit_T={float(result['bed_T_profile'][-1]):.1f}K"
    )
    print("-" * 160)
    print(
        f"{'cell':>4} {'xi/H':>6} {'u0':>8} {'u_d':>8} {'u_b':>8} {'d_b':>8} "
        f"{'num=u0-u_d':>12} {'den':>10} {'eps_b':>8} {'eps_avg':>9} {'bulk_solid':>11}"
    )
    print("-" * 160)

    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        num = float(cell.u0 - cell.u_d)
        den = float(cell.u_b + (2.7 - 1.0) * cell.u_d)
        eps_avg = float(cell.eps_b + (1.0 - cell.eps_b) * cell.eps_d_voidage)
        bulk = calc_bulk_solid_holdup(float(cell.eps_b), float(cell.eps_d_voidage))
        print(
            f"{i:>4d} {xi:>6.2f} {cell.u0:>8.4f} {cell.u_d:>8.4f} {cell.u_b:>8.4f} {cell.d_b:>8.4f} "
            f"{num:>12.4f} {den:>10.4f} {cell.eps_b:>8.4f} {eps_avg:>9.4f} {bulk:>11.4f}"
        )

    print("-" * 160)
    print("readout:")
    print("  - Eq.3.24 uses `eps_b = (u0 - u_d) / (u_b + (n_b - 1) * u_d)`, here with `n_b = 2.7`.")
    print("  - if `u0 - u_d` stays nearly flat but the denominator rises with height, `eps_b` will decline upward.")
    print("  - Hamel bed physics expects the whole-bed average voidage to rise upward inside the bed.")
    print("  - so an upward-declining `eps_b` / `eps_avg` trend is a direct hydrodynamics inconsistency, not a bulk-solid sign error.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
