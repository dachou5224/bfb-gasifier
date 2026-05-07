#!/usr/bin/env python3
"""LU backsolve 分支下 d_b0 口径隔离审计。"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.physics.bubble_dynamics import bubble_diameter_ode, initial_bubble_diameter
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


def _integrate_with_db0(
    *,
    u0: float,
    u_mf: float,
    P: float,
    H_bed: float,
    d_b0: float,
    n_points: int,
    lambda_strategy: str,
    xi_strategy: str,
) -> np.ndarray:
    h_arr = np.linspace(0.0, H_bed, n_points)
    sol = solve_ivp(
        fun=lambda h, y: bubble_diameter_ode(
            h,
            y,
            u0,
            u_mf,
            P,
            lambda_strategy=lambda_strategy,
            xi_strategy=xi_strategy,
        ),
        t_span=(0.0, H_bed),
        y0=[max(float(d_b0), 1e-4)],
        t_eval=h_arr,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"bubble ODE integration failed: {sol.message}")
    return np.maximum(sol.y[0], 1e-4)


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    cell0 = reactor.cells[0]
    A_bed = np.pi / 4.0 * cfg.D_bed**2
    n_or = max(cell0.geo.N_or, 1)
    u0 = float(cell0.u0)
    u_mf = float(cell0.u_mf)
    q0_per_orifice = u0 * A_bed / n_or
    qex_per_orifice = max(u0 - u_mf, 0.0) * A_bed / n_or

    db0_current = initial_bubble_diameter(u0, u_mf, A_bed, n_or)
    db0_hamel_q0 = 1.3 * ((q0_per_orifice**2) / 9.81) ** 0.2
    db0_hamel_qex = 1.3 * ((qex_per_orifice**2) / 9.81) ** 0.2

    n_points = cfg.n_cells + 1
    db_fix_cur = _integrate_with_db0(
        u0=u0,
        u_mf=u_mf,
        P=cfg.P,
        H_bed=cfg.H_bed,
        d_b0=db0_current,
        n_points=n_points,
        lambda_strategy="hamel_280",
        xi_strategy="fixed_035",
    )
    db_fix_hamel_q0 = _integrate_with_db0(
        u0=u0,
        u_mf=u_mf,
        P=cfg.P,
        H_bed=cfg.H_bed,
        d_b0=db0_hamel_q0,
        n_points=n_points,
        lambda_strategy="hamel_280",
        xi_strategy="fixed_035",
    )
    db_fix_hamel_qex = _integrate_with_db0(
        u0=u0,
        u_mf=u_mf,
        P=cfg.P,
        H_bed=cfg.H_bed,
        d_b0=db0_hamel_qex,
        n_points=n_points,
        lambda_strategy="hamel_280",
        xi_strategy="fixed_035",
    )

    print("=" * 184)
    print("LU d_b0-strategy audit on backsolve_visible_epsb branch")
    print("=" * 184)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"u0={u0:.4f} u_mf={u_mf:.4f} A_bed={A_bed:.4f} N_or={n_or}"
    )
    print(
        f"d_b0 current(Darton)={db0_current:.5f} m | "
        f"d_b0 Hamel(q0/orifice)={db0_hamel_q0:.5f} m | "
        f"d_b0 Hamel(qex/orifice)={db0_hamel_qex:.5f} m"
    )
    print("-" * 184)
    print(f"{'cell':>4} {'xi/H':>6} {'d_b cell':>10} {'ODE_cur':>10} {'ODE_q0':>10} {'ODE_qex':>10}")
    print("-" * 184)
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        idx = min(int(round((n_points - 1) * xi)), n_points - 1)
        print(
            f"{i:>4d} {xi:>6.2f} {cell.d_b:>10.4f} {float(db_fix_cur[idx]):>10.4f} "
            f"{float(db_fix_hamel_q0[idx]):>10.4f} {float(db_fix_hamel_qex[idx]):>10.4f}"
        )
    print("-" * 184)
    print(
        f"top d_b mean: current={float(np.mean(db_fix_cur[-4:])):.4f} "
        f"hamel_q0={float(np.mean(db_fix_hamel_q0[-4:])):.4f} "
        f"hamel_qex={float(np.mean(db_fix_hamel_qex[-4:])):.4f}"
    )
    print("notes:")
    print("  - `Hamel(q0/orifice)` uses d_b0 = 1.3*(Vdot0^2/g)^0.2 with Vdot0 inferred as total superficial volumetric flow per distributor opening.")
    print("  - `Hamel(qex/orifice)` uses the same form but with excess gas (u0-u_mf) per opening; this is an audit-only interpretation, not yet source-of-truth.")
    print("readout:")
    print("  - 若两种 Hamel-style d_b0 只带来小幅变化，则当前 gap 主要不在 d_b0，而在 growth/decay 主体。")
    print("  - 若 q0/qex 解释能显著拉高 ODE top d_b，则下一步应回到 distributor/orifice 物理定义。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
