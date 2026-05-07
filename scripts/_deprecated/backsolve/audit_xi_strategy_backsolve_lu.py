#!/usr/bin/env python3
"""LU backsolve 分支下 xi_b 策略隔离审计。"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.physics.bubble_dynamics import integrate_bubble_diameter, resolve_xi_b
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
    u_d_ref = float(reactor.cells[0].u_d)
    n_points = cfg.n_cells + 1
    _, db_fix = integrate_bubble_diameter(
        u0=u0_ref,
        u_mf=u_mf_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        u_d=u_d_ref,
        D_bed=cfg.D_bed,
        n_points=n_points,
        method="hilligardt_ode",
        lambda_strategy="hamel_280",
        xi_strategy="fixed_035",
    )
    _, db_dyn = integrate_bubble_diameter(
        u0=u0_ref,
        u_mf=u_mf_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        u_d=u_d_ref,
        D_bed=cfg.D_bed,
        n_points=n_points,
        method="hilligardt_ode",
        lambda_strategy="hamel_280",
        xi_strategy="hamel_regime",
    )

    print("=" * 164)
    print("LU xi-strategy audit on backsolve_visible_epsb branch")
    print("=" * 164)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"u0_ref={u0_ref:.4f} u_mf_ref={u_mf_ref:.4f} u_d_ref={u_d_ref:.4f}"
    )
    print("-" * 164)
    print(f"{'cell':>4} {'xi/H':>6} {'alpha':>9} {'xi_fix':>8} {'xi_dyn':>8} {'d_b cell':>10} {'ODE_fix':>10} {'ODE_dyn':>10}")
    print("-" * 164)
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        idx = min(int(round((n_points - 1) * xi)), n_points - 1)
        alpha = float(cell.u_b / max(cell.u_d, 1e-12))
        xi_fix = resolve_xi_b(alpha, strategy="fixed_035")
        xi_dyn = resolve_xi_b(alpha, strategy="hamel_regime")
        print(
            f"{i:>4d} {xi:>6.2f} {alpha:>9.2f} {xi_fix:>8.3f} {xi_dyn:>8.3f} "
            f"{cell.d_b:>10.4f} {float(db_fix[idx]):>10.4f} {float(db_dyn[idx]):>10.4f}"
        )

    print("-" * 164)
    print(f"top d_b mean: fixed_035={float(np.mean(db_fix[-4:])):.4f} hamel_regime={float(np.mean(db_dyn[-4:])):.4f}")
    print("readout:")
    print("  - 这里的 alpha 已按 Hamel 原文改为 `u_b/u_d`，不再是 `u_b/u_mf`。")
    print("  - 若在 `lambda_b=hamel_280` 固定后，`hamel_regime` 仍显著压小 d_b，则 `xi_b` 也是主偏差之一。")
    print("  - 若 fixed_035 与 hamel_regime 差别远小于 current vs hamel_280，则下一优先级仍是 `lambda_b`。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
