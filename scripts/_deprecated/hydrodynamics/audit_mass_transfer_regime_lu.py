#!/usr/bin/env python3
"""LU 工况 mass-transfer regime / K_bd 审计。"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import g
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from src.core.species import gas_diffusivity_correlation
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


_CLOSURES = ("current", "backsolve_visible_epsb")


def _kbd_fast_branch(*, u_mf: float, d_b: float, D_g: float, u_b: float, eps_mf: float) -> tuple[float, float, float]:
    k_bc = 4.5 * u_mf / max(d_b, 1e-12) + 5.85 * math.sqrt(max(D_g, 0.0)) * (g ** 0.25) / (max(d_b, 1e-12) ** 1.25)
    k_cd = 6.78 * math.sqrt(max(D_g, 0.0) * max(u_b, 0.0) * max(eps_mf, 0.0) / (max(d_b, 1e-12) ** 3))
    k_bd = 1.0 / max((1.0 / max(k_bc, 1e-12)) + (1.0 / max(k_cd, 1e-12)), 1e-12)
    return float(k_bc), float(k_cd), float(k_bd)


def _kbd_slow_preto(*, u_mf: float, d_b: float, chi: float = 1.0) -> float:
    return float(4.5 * chi * u_mf / max(d_b, 1e-12))


def main() -> int:
    case = load_case_LU()
    print("=" * 232)
    print("LU mass-transfer regime / K_bd audit")
    print("=" * 232)
    print(
        "compare two alpha definitions and three K_bd lenses:\n"
        "  1. current implementation Eq.3.50 mixed form\n"
        "  2. fast-bubble harmonic branch (Kunii/Levenspiel style, via chat-record extract)\n"
        "  3. Preto Eq.3.49 lens with discrete chi (slow=1, fast min=0.1)"
    )

    for closure in _CLOSURES:
        cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
        cfg.hydrodynamics_u_d_closure = closure
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
            check_x0=True,
        )

        disagreement = 0
        print("\n" + "-" * 232)
        print_nr_monitor(monitor, prefix=f"NR monitor [{closure}]")
        print(
            f"closure={closure} Texit={result['T_profile'][-1]:.1f}K n_iter={result['n_iter']} "
            f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
        )
        print("-" * 232)
        print(
            f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'eps_b':>8} {'u_d':>8} {'u_b':>8} "
            f"{'alpha_umf':>11} {'reg_umf':>8} {'alpha_ud':>11} {'reg_ud':>8} "
            f"{'Kbd_cur':>10} {'Kbc':>10} {'Kcd':>10} {'Kbd_fast':>10} {'Kbd_slow':>10}"
        )
        print("-" * 232)

        for i, cell in enumerate(reactor.cells):
            xi = float(cell.geo.h_center / cfg.H_bed)
            alpha_umf = float(cell.u_b / max(cell.u_mf, 1e-12))
            alpha_ud = float(cell.u_b / max(cell.u_d, 1e-12))
            reg_umf = "slow" if alpha_umf < 1.0 else "fast"
            reg_ud = "slow" if alpha_ud < 1.0 else "fast"
            if reg_umf != reg_ud:
                disagreement += 1

            D_g = float(gas_diffusivity_correlation(cell.T, cell.P))
            k_bc, k_cd, k_fast = _kbd_fast_branch(
                u_mf=float(cell.u_mf),
                d_b=float(cell.d_b),
                D_g=D_g,
                u_b=float(cell.u_b),
                eps_mf=float(cell.solid.eps_mf),
            )
            chi_preto = 1.0 if reg_ud == "slow" else 0.1
            k_slow = _kbd_slow_preto(u_mf=float(cell.u_mf), d_b=float(cell.d_b), chi=chi_preto)
            print(
                f"{i:>4d} {xi:>6.2f} {cell.T:>8.1f} {cell.eps_b:>8.4f} {cell.u_d:>8.4f} {cell.u_b:>8.4f} "
                f"{alpha_umf:>11.3f} {reg_umf:>8} {alpha_ud:>11.3f} {reg_ud:>8} "
                f"{cell.K_bd:>10.3f} {k_bc:>10.3f} {k_cd:>10.3f} {k_fast:>10.3f} {k_slow:>10.3f}"
            )

        print("-" * 232)
        print(f"regime disagreements (alpha_umf vs alpha_ud): {disagreement}/{len(reactor.cells)} cells")
        print(
            "readout: Preto Eq.3.49 is used here only as a literature lens. "
            "Hamel final main path remains Eq.3.50, where chi no longer appears explicitly."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
