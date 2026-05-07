#!/usr/bin/env python3
"""LU backsolve 分支下 bubble ODE 各项拆解。"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import P0, g
from src.core.reactor import Reactor
from src.physics.bubble_dynamics import bubble_lifetime
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


def _lambda_hamel(u_mf: float, P: float) -> float:
    return 280.0 * float(u_mf) / g * (float(P) / P0) ** (-0.7)


def _dynamic_xi(alpha_b: float) -> float:
    if alpha_b <= 1.0:
        return max(0.0, 1.0 - alpha_b**3)
    return 0.0


def _growth_term(eps_b: float, xi_b: float) -> float:
    eps_b = float(np.clip(eps_b, 1e-8, 0.95))
    eps13 = eps_b ** (1.0 / 3.0)
    coeff = (6.0 / np.pi) ** (1.0 / 3.0)
    denom = max(1.0 - xi_b * coeff * eps13, 0.01)
    return (2.0 / (9.0 * np.pi)) * eps13 / denom


def _decay_term(d_b: float, u_b: float, lam_b: float) -> float:
    return float(d_b) / (3.0 * max(float(lam_b) * float(u_b), 1e-12))


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    print("=" * 196)
    print("LU bubble ODE term audit on backsolve_visible_epsb branch")
    print("=" * 196)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("-" * 196)
    print(
        f"{'cell':>4} {'xi/H':>6} {'d_b':>8} {'u_b':>8} {'eps_vis':>8} {'eps_exc':>8} {'alpha':>8} "
        f"{'g_fix':>9} {'g_dyn_vis':>11} {'g_dyn_exc':>11} {'dec_cur':>9} {'dec_ham':>9} {'dd_fix':>9} {'dd_dyn_ham':>12}"
    )
    print("-" * 196)

    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        alpha = float(cell.u_b / max(cell.u_mf, 1e-12))
        eps_vis = float(cell.eps_b)
        eps_exc = float(np.clip((cell.u0 - cell.u_mf) / max(cell.u_b, 1e-12), 1e-8, 0.95))
        xi_fix = 0.35
        xi_dyn = _dynamic_xi(alpha)

        g_fix = _growth_term(eps_vis, xi_fix)
        g_dyn_vis = _growth_term(eps_vis, xi_dyn)
        g_dyn_exc = _growth_term(eps_exc, xi_dyn)
        dec_cur = _decay_term(cell.d_b, cell.u_b, bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="current"))
        dec_ham = _decay_term(cell.d_b, cell.u_b, _lambda_hamel(cell.u_mf, cfg.P))
        dd_fix = g_fix - dec_cur
        dd_dyn_ham = g_dyn_vis - dec_ham

        print(
            f"{i:>4d} {xi:>6.2f} {cell.d_b:>8.4f} {cell.u_b:>8.4f} {eps_vis:>8.4f} {eps_exc:>8.4f} {alpha:>8.2f} "
            f"{g_fix:>9.4f} {g_dyn_vis:>11.4f} {g_dyn_exc:>11.4f} {dec_cur:>9.4f} {dec_ham:>9.4f} {dd_fix:>9.4f} {dd_dyn_ham:>12.4f}"
        )

    print("-" * 196)
    print("readout:")
    print("  - 若 `eps_vis` 与 `eps_exc` 接近，而 `g_dyn_vis` 与 `g_dyn_exc` 也几乎相同，则 ODE 崩塌主因不在 `eps_b` 定义。")
    print("  - 若 `g_fix >> g_dyn_vis`，说明把 `xi_b` 从 0.35 切到 fast-bubble 动态 0 会显著压缩 growth term。")
    print("  - 若 `dd_dyn_ham` 全床都接近 0 或为负，而 `dd_fix` 为正，则当前 source-of-truth 最该继续核 `xi_b/lambda_b`。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
