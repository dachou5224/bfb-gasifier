#!/usr/bin/env python3
"""LU backsolve 分支下 bubble ODE 所需 closure 反推审计。"""

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


def _growth_prefactors(eps_b: float) -> tuple[float, float]:
    eps = float(np.clip(eps_b, 1e-8, 0.95))
    eps13 = eps ** (1.0 / 3.0)
    a = (2.0 / (9.0 * np.pi)) * eps13
    b = (6.0 / np.pi) ** (1.0 / 3.0) * eps13
    return a, b


def _growth(eps_b: float, xi_b: float) -> float:
    a, b = _growth_prefactors(eps_b)
    denom = max(1.0 - xi_b * b, 0.01)
    return a / denom


def _lambda_required_for_balance(d_b: float, u_b: float, growth: float) -> float:
    return float(d_b) / max(3.0 * float(u_b) * float(growth), 1e-12)


def _xi_required_for_balance(eps_b: float, target_growth: float) -> float | None:
    a, b = _growth_prefactors(eps_b)
    if target_growth <= 0.0:
        return None
    # target_growth = a / (1 - xi*b)
    xi = (1.0 - a / target_growth) / max(b, 1e-12)
    return float(xi)


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    print("=" * 208)
    print("LU bubble required-closure audit on backsolve_visible_epsb branch")
    print("=" * 208)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("-" * 208)
    print(
        f"{'cell':>4} {'xi/H':>6} {'d_b':>8} {'u_b':>8} {'eps_b':>8} "
        f"{'g_fix':>9} {'lam_cur':>9} {'lam_ham':>9} {'lam_req':>9} "
        f"{'req/cur':>9} {'req/ham':>9} {'xi_req(cur)':>12} {'xi_req(ham)':>12}"
    )
    print("-" * 208)

    req_over_cur: list[float] = []
    req_over_ham: list[float] = []
    xi_req_cur_vals: list[float] = []
    xi_req_ham_vals: list[float] = []

    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        g_fix = _growth(float(cell.eps_b), 0.35)
        lam_cur = float(bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="current"))
        lam_ham = float(bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="hamel_280", u_mf=cell.u_mf))
        lam_req = _lambda_required_for_balance(cell.d_b, cell.u_b, g_fix)
        ratio_cur = lam_req / max(lam_cur, 1e-12)
        ratio_ham = lam_req / max(lam_ham, 1e-12)

        # Required xi if lambda is frozen and we want growth == decay
        target_g_cur = float(cell.d_b) / max(3.0 * lam_cur * cell.u_b, 1e-12)
        target_g_ham = float(cell.d_b) / max(3.0 * lam_ham * cell.u_b, 1e-12)
        xi_req_cur = _xi_required_for_balance(float(cell.eps_b), target_g_cur)
        xi_req_ham = _xi_required_for_balance(float(cell.eps_b), target_g_ham)

        req_over_cur.append(ratio_cur)
        req_over_ham.append(ratio_ham)
        if xi_req_cur is not None:
            xi_req_cur_vals.append(xi_req_cur)
        if xi_req_ham is not None:
            xi_req_ham_vals.append(xi_req_ham)

        print(
            f"{i:>4d} {xi:>6.2f} {cell.d_b:>8.4f} {cell.u_b:>8.4f} {cell.eps_b:>8.4f} "
            f"{g_fix:>9.4f} {lam_cur:>9.4f} {lam_ham:>9.4f} {lam_req:>9.4f} "
            f"{ratio_cur:>9.2f} {ratio_ham:>9.2f} "
            f"{(float('nan') if xi_req_cur is None else xi_req_cur):>12.3f} "
            f"{(float('nan') if xi_req_ham is None else xi_req_ham):>12.3f}"
        )

    print("-" * 208)
    print(
        f"median required-lambda multipliers: req/current={float(np.median(req_over_cur)):.2f} "
        f"req/hamel={float(np.median(req_over_ham)):.2f}"
    )
    print(
        f"median xi_required: for current-lambda={float(np.median(xi_req_cur_vals)):.3f} "
        f"for hamel-lambda={float(np.median(xi_req_ham_vals)):.3f}"
    )
    print("readout:")
    print("  - 若 `lam_req` 普遍是 current/hamel 的数倍，而 `xi_req` 要求超过 ~0.35 甚至接近/超过 1，则主偏差更像 `lambda_b` 太小。")
    print("  - 若 `xi_req(ham)` 已落在合理区间，而 `lam_req/ham` 接近 1，则下一步应更多回到 `xi_b`。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
