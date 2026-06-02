#!/usr/bin/env python3
"""LU backsolve 分支下 constant-bottom ODE vs piecewise-local ODE 审计。"""

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


def _piecewise_local_ode(
    reactor: Reactor,
    *,
    lambda_strategy: str,
    xi_strategy: str,
) -> np.ndarray:
    cells = reactor.cells
    cfg = reactor.config
    A_bed = np.pi / 4.0 * cfg.D_bed**2
    d = max(initial_bubble_diameter(cells[0].u0, cells[0].u_mf, A_bed, cells[0].geo.N_or), 1e-3)
    out = np.zeros(len(cells) + 1, dtype=float)
    out[0] = d
    for i, cell in enumerate(cells):
        dh = float(cell.geo.dh)
        sol = solve_ivp(
            fun=lambda h, y: bubble_diameter_ode(
                h,
                y,
                float(cell.u0),
                float(cell.u_mf),
                float(cell.P),
                lambda_strategy=lambda_strategy,
                xi_strategy=xi_strategy,
            ),
            t_span=(0.0, dh),
            y0=[max(float(d), 1e-4)],
            t_eval=[dh],
            method="RK45",
            rtol=1e-8,
            atol=1e-10,
        )
        if not sol.success:
            raise RuntimeError(f"piecewise-local bubble ODE failed at cell {i}: {sol.message}")
        d = max(float(sol.y[0][-1]), 1e-4)
        out[i + 1] = d
    return out


def _constant_bottom_ode(
    reactor: Reactor,
    *,
    lambda_strategy: str,
    xi_strategy: str,
) -> np.ndarray:
    cells = reactor.cells
    cfg = reactor.config
    A_bed = np.pi / 4.0 * cfg.D_bed**2
    d0 = max(initial_bubble_diameter(cells[0].u0, cells[0].u_mf, A_bed, cells[0].geo.N_or), 1e-3)
    h_arr = np.linspace(0.0, cfg.H_bed, len(cells) + 1)
    sol = solve_ivp(
        fun=lambda h, y: bubble_diameter_ode(
            h,
            y,
            float(cells[0].u0),
            float(cells[0].u_mf),
            float(cfg.P),
            lambda_strategy=lambda_strategy,
            xi_strategy=xi_strategy,
        ),
        t_span=(0.0, cfg.H_bed),
        y0=[d0],
        t_eval=h_arr,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"constant-bottom bubble ODE failed: {sol.message}")
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

    db_const = _constant_bottom_ode(reactor, lambda_strategy="hamel_280", xi_strategy="fixed_035")
    db_local = _piecewise_local_ode(reactor, lambda_strategy="hamel_280", xi_strategy="fixed_035")

    print("=" * 176)
    print("LU constant-bottom vs piecewise-local bubble ODE audit")
    print("=" * 176)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("-" * 176)
    print(f"{'cell':>4} {'xi/H':>6} {'u0':>8} {'u_mf':>8} {'d_b cell':>10} {'ODE_const':>10} {'ODE_local':>10}")
    print("-" * 176)
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        idx = i + 1
        print(
            f"{i:>4d} {xi:>6.2f} {cell.u0:>8.4f} {cell.u_mf:>8.4f} {cell.d_b:>10.4f} "
            f"{float(db_const[idx]):>10.4f} {float(db_local[idx]):>10.4f}"
        )
    print("-" * 176)
    print(
        f"top d_b mean: const={float(np.mean(db_const[-4:])):.4f} "
        f"local={float(np.mean(db_local[-4:])):.4f} "
        f"cell={float(np.mean([c.d_b for c in reactor.cells[-4:]])):.4f}"
    )
    print("readout:")
    print("  - 若 piecewise-local ODE 仍接近 constant-bottom ODE，则剩余 gap 基本不在 `u0/u_mf` 的沿程变化。")
    print("  - 若 piecewise-local ODE 显著抬高 d_b，则之前的常参 ODE 对比本身就低估了 bubble growth。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
