#!/usr/bin/env python3
"""LU backsolve 分支下 lambda_b 策略隔离审计。"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.physics.bubble_dynamics import bubble_lifetime, integrate_bubble_diameter
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    u0_ref = float(reactor.cells[0].u0)
    u_mf_ref = float(reactor.cells[0].u_mf)
    n_points = cfg.n_cells + 1
    h_arr, db_cur = integrate_bubble_diameter(
        u0=u0_ref,
        u_mf=u_mf_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        n_points=n_points,
        method="hilligardt_ode",
        lambda_strategy="current",
    )
    _, db_ham = integrate_bubble_diameter(
        u0=u0_ref,
        u_mf=u_mf_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        n_points=n_points,
        method="hilligardt_ode",
        lambda_strategy="hamel_280",
    )

    print("=" * 160)
    print("LU lambda-strategy audit on backsolve_visible_epsb branch")
    print("=" * 160)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"u0_ref={u0_ref:.4f} u_mf_ref={u_mf_ref:.4f}"
    )
    print("-" * 160)
    print(f"{'cell':>4} {'xi/H':>6} {'d_b cell':>10} {'ODE_cur':>10} {'ODE_ham':>10} {'lam_cur':>10} {'lam_ham':>10}")
    print("-" * 160)

    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        idx = min(int(round((n_points - 1) * xi)), n_points - 1)
        lam_cur = bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="current")
        lam_ham = bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="hamel_280", u_mf=cell.u_mf)
        print(
            f"{i:>4d} {xi:>6.2f} {cell.d_b:>10.4f} {float(db_cur[idx]):>10.4f} {float(db_ham[idx]):>10.4f} "
            f"{lam_cur:>10.4f} {lam_ham:>10.4f}"
        )

    print("-" * 160)
    print(
        f"top d_b mean: current={float(np.mean(db_cur[-4:])):.4f} hamel_280={float(np.mean(db_ham[-4:])):.4f}"
    )
    print("readout:")
    print("  - 若固定 xi_b=0.35 时，`hamel_280` 仍显著压小 d_b，则 `lambda_b` 本身就是下一优先级 source-of-truth。")
    print("  - 若 `hamel_280` 反而更接近当前 cell/MW 量级，则问题应主要回到 `xi_b`。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
