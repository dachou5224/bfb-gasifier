#!/usr/bin/env python3
"""LU 工况下 bubble-chain 的 bed hydrodynamics 审计。"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import P0_HAMEL, g
from src.core.reactor import Reactor
from src.physics.bubble_dynamics import (
    bubble_lifetime,
    bubble_rise_velocity,
    initial_bubble_diameter,
    integrate_bubble_diameter,
)
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


def _lambda_hamel(u_mf: float, P: float) -> float:
    return 280.0 * float(u_mf) / g * (float(P) / P0_HAMEL) ** (-0.7)


def _dynamic_xi(alpha_b: float) -> float:
    if alpha_b <= 1.0:
        return max(0.0, 1.0 - alpha_b**3)
    return 0.0


def _bubble_ode_variant(
    *,
    h: float,
    y: np.ndarray,
    u0: float,
    u_mf: float,
    u_d: float,
    P: float,
    psi_b: float,
    xi_mode: str,
    lambda_mode: str,
) -> np.ndarray:
    d_b = max(float(y[0]), 1e-4)
    u_b = bubble_rise_velocity(u0, u_mf, d_b, psi_b)

    if xi_mode == "fixed_035":
        xi_b = 0.35
    elif xi_mode == "dynamic":
        alpha_b = u_b / max(u_d, 1e-12)
        xi_b = _dynamic_xi(alpha_b)
    else:
        raise ValueError(f"unknown xi_mode={xi_mode!r}")

    if lambda_mode == "current":
        lam_b = bubble_lifetime(d_b, u_b, P, strategy="current")
    elif lambda_mode == "hamel":
        lam_b = _lambda_hamel(u_mf, P)
    else:
        raise ValueError(f"unknown lambda_mode={lambda_mode!r}")

    excess = max(u0 - u_mf, 1e-10)
    eps_b = np.clip(excess / max(u_b, 1e-12), 1e-6, 0.8)
    eps_b_13 = eps_b ** (1.0 / 3.0)
    coeff_6pi = (6.0 / np.pi) ** (1.0 / 3.0)
    denom = max(1.0 - xi_b * coeff_6pi * eps_b_13, 0.01)
    growth = (2.0 / (9.0 * np.pi)) * eps_b_13 / denom
    decay = d_b / (3.0 * max(lam_b * u_b, 1e-10))
    return np.array([growth - decay], dtype=float)


def _integrate_variant(
    *,
    u0: float,
    u_mf: float,
    u_d: float,
    P: float,
    H_bed: float,
    D_bed: float,
    N_or: int,
    n_points: int,
    psi_b: float,
    xi_mode: str,
    lambda_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    A_bed = np.pi / 4.0 * D_bed**2
    d_b0 = max(initial_bubble_diameter(u0, u_mf, A_bed, N_or), 1e-3)
    h_arr = np.linspace(0.0, H_bed, n_points)
    sol = solve_ivp(
        fun=lambda h, y: _bubble_ode_variant(
            h=h,
            y=y,
            u0=u0,
            u_mf=u_mf,
            u_d=u_d,
            P=P,
            psi_b=psi_b,
            xi_mode=xi_mode,
            lambda_mode=lambda_mode,
        ),
        t_span=(0.0, H_bed),
        y0=[d_b0],
        t_eval=h_arr,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"bubble ODE integration failed: {sol.message}")
    return h_arr, np.maximum(sol.y[0], 1e-4)


def main() -> int:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
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
    h_ref, db_mw = integrate_bubble_diameter(
        u0=u0_ref,
        u_mf=u_mf_ref,
        u_d=u_d_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        n_points=n_points,
        method="mori_wen",
    )
    _, db_ode_fixed_current = _integrate_variant(
        u0=u0_ref,
        u_mf=u_mf_ref,
        u_d=u_d_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        N_or=reactor.cells[0].geo.N_or,
        n_points=n_points,
        psi_b=0.76,
        xi_mode="fixed_035",
        lambda_mode="current",
    )
    _, db_ode_dyn_current = _integrate_variant(
        u0=u0_ref,
        u_mf=u_mf_ref,
        u_d=u_d_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        N_or=reactor.cells[0].geo.N_or,
        n_points=n_points,
        psi_b=0.76,
        xi_mode="dynamic",
        lambda_mode="current",
    )
    _, db_ode_dyn_hamel = _integrate_variant(
        u0=u0_ref,
        u_mf=u_mf_ref,
        u_d=u_d_ref,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        N_or=reactor.cells[0].geo.N_or,
        n_points=n_points,
        psi_b=0.76,
        xi_mode="dynamic",
        lambda_mode="hamel",
    )

    print("=" * 196)
    print("LU bubble-chain audit on backsolve_visible_epsb branch")
    print("=" * 196)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"u0_ref={u0_ref:.4f} u_mf_ref={u_mf_ref:.4f} u_d_ref={u_d_ref:.4f} P={cfg.P/1e6:.2f}MPa"
    )
    print("variants:")
    print("  - MW: current main-path Mori-Wen")
    print("  - ODE_fix_cur: xi_b=0.35 + current lambda")
    print("  - ODE_dyn_cur: xi_b follows Hamel Eq.3.34 with alpha_b=u_b/u_d + current lambda")
    print("  - ODE_dyn_hamel: xi_b follows Hamel Eq.3.34 with alpha_b=u_b/u_d + 280*u_mf/g*(P/P0)^(-0.7)")
    print("-" * 196)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'u_b':>8} {'eps_b':>8} {'alpha':>10} "
        f"{'d_b cell':>10} {'MW':>10} {'ODE_fix':>10} {'ODE_dyn':>10} {'ODE_ham':>10} "
        f"{'lam_cur':>10} {'lam_ham':>10}"
    )
    print("-" * 196)

    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        idx = min(int(round((n_points - 1) * xi)), n_points - 1)
        alpha_b = float(cell.u_b / max(cell.u_d, 1e-12))
        lam_cur = float(bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="current"))
        lam_ham = float(_lambda_hamel(cell.u_mf, cfg.P))
        print(
            f"{i:>4d} {xi:>6.2f} {cell.T:>8.1f} {cell.u_b:>8.4f} {cell.eps_b:>8.4f} {alpha_b:>10.3f} "
            f"{cell.d_b:>10.4f} {float(db_mw[idx]):>10.4f} {float(db_ode_fixed_current[idx]):>10.4f} "
            f"{float(db_ode_dyn_current[idx]):>10.4f} {float(db_ode_dyn_hamel[idx]):>10.4f} "
            f"{lam_cur:>10.4f} {lam_ham:>10.4f}"
        )

    print("-" * 196)
    print("readout:")
    print("  - 若 `ODE_dyn_hamel` 显著压小 top-bed d_b，而 `MW` 与 `ODE_fix_cur` 仍偏大，则剩余主偏差更可能在 `xi_b/lambda_b` 而不是 `u_d`。")
    print("  - 若 `MW` 与 `ODE_dyn_hamel` 量级接近，则下一步更该回到 `K_bd` 而不是继续追 `d_b`。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
