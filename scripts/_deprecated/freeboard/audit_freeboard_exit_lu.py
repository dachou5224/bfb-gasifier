#!/usr/bin/env python3
"""LU 工况下 freeboard-aware 出口审计。"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_validation_case_node,
    CASE_LU_VALIDATION_KEY,
)


def main() -> int:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )
    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    ref_dry = ref["exit_gas_dry_mol_frac"]

    print("=" * 120)
    print("LU freeboard-aware exit audit")
    print("=" * 120)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"freeboard_active={result['freeboard_active']} "
        f"H_freeboard={cfg.H_freeboard:.2f}m n_freeboard_cells={cfg.n_freeboard_cells} "
        f"beta_A={result['freeboard_beta_A']:.4f}"
    )
    print(
        f"bed_exit_T={result['bed_T_profile'][-1]:.1f}K "
        f"reactor_exit_T={result['reactor_exit_T']:.1f}K "
        f"ref_exit_T={float(ref['exit_temperature_K']):.1f}K"
    )
    print("-" * 120)
    print(f"{'species':>8} {'bed_exit':>12} {'reactor_exit':>14} {'ref_exit':>10}")
    print("-" * 120)
    for sp in ("CO", "CO2", "H2", "CH4"):
        print(
            f"{sp:>8s} {result['bed_exit_gas_dry'].get(sp, 0.0):>12.4f} "
            f"{result['exit_gas_dry'].get(sp, 0.0):>14.4f} {float(ref_dry.get(sp, 0.0)):>10.4f}"
        )
    print("-" * 120)
    print(
        f"entrained eject={result['freeboard_entrained_eject_char_ash_kg_s']:.4f} kg/s  "
        f"freeboard exit char={result['freeboard_entrained_exit_char_kg_s']:.4f} kg/s  "
        f"ash={result['freeboard_entrained_exit_ash_kg_s']:.4f} kg/s"
    )
    print(
        f"return-to-bed char={result['freeboard_entrained_return_char_kg_s']:.4f} kg/s  "
        f"ash={result['freeboard_entrained_return_ash_kg_s']:.4f} kg/s"
    )
    print(
        f"cyclone capture char={result['freeboard_cyclone_capture_char_kg_s']:.4f} kg/s  "
        f"ash={result['freeboard_cyclone_capture_ash_kg_s']:.4f} kg/s  "
        f"recycle candidate={result['freeboard_cyclone_recycle_candidate_char_ash_kg_s']:.4f} kg/s"
    )
    if result["freeboard_carry_ratio_profile"]:
        print(
            f"u0_end={result['freeboard_u0_profile_m_s'][-1]:.4f} m/s  "
            f"u_gb_end={result['freeboard_u_gb_profile_m_s'][-1]:.4f} m/s  "
            f"u_p_end={result['freeboard_u_p_mean_profile_m_s'][-1]:.4f} m/s  "
            f"u_t_end={result['freeboard_u_t_mean_profile_m_s'][-1]:.4f} m/s  "
            f"carry_ratio_end={result['freeboard_carry_ratio_profile'][-1]:.3f}"
        )
    print("-" * 120)
    print(
        f"bed cells={len(result['bed_T_profile'])} "
        f"freeboard cells={len(result['freeboard_T_profile'])} "
        f"reactor axial points={len(result['T_profile'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
