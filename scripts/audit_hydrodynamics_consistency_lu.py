#!/usr/bin/env python3
"""LU 工况 bed hydrodynamics 与 Hamel 一致性审计。"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import P0, g
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from src.physics.bubble_dynamics import (
    bubble_lifetime,
    bubble_rise_velocity,
    integrate_bubble_diameter,
)
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


def _hamel_lambda_b(u_mf: float, P: float) -> float:
    """Hamel 本地提取：lambda_b = 280*u_mf/g*(P/P0)^(-0.7)."""
    return 280.0 * float(u_mf) / g * (float(P) / P0) ** (-0.7)


def _raw_eps_d_voidage(u_d: float, u_mf: float, eps_mf: float, n_rz: float) -> float:
    ratio = max(float(u_d) / max(float(u_mf), 1e-12), 1.0)
    return float(eps_mf) * math.pow(ratio, 1.0 / max(float(n_rz), 1e-12))


def _wein_ud_potential(u_mf: float) -> float:
    return 1.45 * float(u_mf)


def main() -> int:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    reactor = Reactor(cfg)

    # 先看初始化 hydrodynamics
    for cell in reactor.cells:
        cell.calc_hydrodynamics()
    init_rows = []
    for i, cell in enumerate(reactor.cells):
        init_rows.append(
            {
                "cell": i,
                "xi": float(cell.geo.h_center / cfg.H_bed),
                "eps_b": float(cell.eps_b),
                "eps_d_void": float(cell.eps_d_voidage),
                "u0": float(cell.u0),
                "u_mf": float(cell.u_mf),
                "u_d": float(cell.u_d),
                "u_b": float(cell.u_b),
                "d_b": float(cell.d_b),
                "n_rz": float(cell.n_rz),
            }
        )

    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    # 以 solved bed-bottom 条件走一条 ODE bubble path，仅作当前主路径对照
    ref_u0 = float(reactor.cells[0].u0)
    ref_umf = float(reactor.cells[0].u_mf)
    h_arr, db_ode = integrate_bubble_diameter(
        u0=ref_u0,
        u_mf=ref_umf,
        P=cfg.P,
        H_bed=cfg.H_bed,
        D_bed=cfg.D_bed,
        n_points=reactor.config.n_cells + 1,
        method="hilligardt_ode",
        lambda_strategy="current",
        xi_strategy="fixed_035",
    )

    print("=" * 196)
    print("LU bed hydrodynamics consistency audit")
    print("=" * 196)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"P={cfg.P/1e6:.2f}MPa H_bed={cfg.H_bed:.2f}m thesis_mode={cfg.thesis_mode}"
    )
    print("implementation checkpoints:")
    print("  - thesis_mode now locks the bed path to Hilligardt/Heinbockel ODE chain (`d_b/u_b/lambda_b`).")
    print("  - Wein Eq.3.12 `u_d = 1.45*u_mf` is still the default Geldart-B closure, but low-`u0` cases now trigger a micro-bubbling guard instead of allowing `u_d >= u0`.")
    print("  - current `K_bd` path still uses Eq.3.50 mixed form; fast/slow bubble branch logic is not active in main path.")
    print("-" * 196)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'u0':>8} {'u_mf':>8} {'u_d':>8} {'u_d/u_mf':>10} "
        f"{'eps_b':>8} {'alpha=ub/umf':>12} {'alpha=ub/u_d':>14} "
        f"{'eps_d_raw':>10} {'eps_d':>8} {'solid':>8} {'n_rz':>8} "
        f"{'u_d,pot':>9} {'guard?':>7} "
        f"{'d_b MW':>9} {'d_b ODE':>9} {'lam_cur':>10} {'lam_hamel':>10}"
    )
    print("-" * 196)
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        ratio = float(cell.u_d / max(cell.u_mf, 1e-12))
        alpha_umf = float(cell.u_b / max(cell.u_mf, 1e-12))
        alpha_ud = float(cell.u_b / max(cell.u_d, 1e-12))
        eps_d_raw = _raw_eps_d_voidage(cell.u_d, cell.u_mf, cell.solid.eps_mf, cell.n_rz)
        u_d_pot = _wein_ud_potential(cell.u_mf)
        guarded = bool(cell.u_d < min(u_d_pot, cell.u0) - 1e-9)
        ode_idx = min(int(round((len(db_ode) - 1) * xi)), len(db_ode) - 1)
        lam_cur = bubble_lifetime(cell.d_b, cell.u_b, cfg.P, strategy="current")
        lam_hamel = _hamel_lambda_b(cell.u_mf, cfg.P)
        print(
            f"{i:>4d} {xi:>6.2f} {cell.T:>8.1f} {cell.u0:>8.4f} {cell.u_mf:>8.4f} {cell.u_d:>8.4f} {ratio:>10.3f} "
            f"{cell.eps_b:>8.4f} {alpha_umf:>12.3f} {alpha_ud:>14.5f} "
            f"{eps_d_raw:>10.4f} {cell.eps_d_voidage:>8.4f} {1.0-cell.eps_d_voidage:>8.4f} {cell.n_rz:>8.3f} "
            f"{u_d_pot:>9.4f} {str(guarded):>7} "
            f"{cell.d_b:>9.4f} {float(db_ode[ode_idx]):>9.4f} {lam_cur:>10.4f} {lam_hamel:>10.4f}"
        )

    print("-" * 196)
    print("top-bed initialization snapshot:")
    print(f"{'cell':>4} {'xi/H':>6} {'u0':>8} {'u_mf':>8} {'u_d':>8} {'eps_b':>8} {'eps_d':>8}")
    for row in init_rows[-4:]:
        print(
            f"{row['cell']:>4d} {row['xi']:>6.2f} {row['u0']:>8.4f} {row['u_mf']:>8.4f} "
            f"{row['u_d']:>8.4f} {row['eps_b']:>8.4f} {row['eps_d_void']:>8.4f}"
        )

    print("-" * 196)
    print("readout:")
    print("  - 若 `guard?=True` 仍大量出现，说明当前核心矛盾已转为 `u0` / `u_mf` 量级，而不再是 Eq.3.12 直接把 `u_d` 算过头。")
    print("  - 若 solved `eps_d_raw` 全床都 >> 1 且被 clip 到 0.99，而初始化时 top-bed `eps_d≈0.55` 正常，则主问题是 solved-state hydrodynamics 闭环。")
    print("  - 这里同时保留 `alpha=ub/umf` 与 Hamel 原文 `alpha=ub/u_d`，用于审计当前 specs/techspec 的旧口径偏差。")
    print("  - 若 `lam_cur` 与 `lam_hamel` 相差 1-2 个数量级，则当前 ODE bubble path 也不能直接视为与 Hamel 一致。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
