#!/usr/bin/env python3
"""LU 工况下 bed-exit -> reactor-exit 反应预算审计。"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
DEPRECATED_CELL_AUDITS_DIR = SCRIPTS_DIR / "_deprecated" / "cell"
if str(DEPRECATED_CELL_AUDITS_DIR) not in sys.path:
    sys.path.insert(0, str(DEPRECATED_CELL_AUDITS_DIR))

from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from audit_ch4_h2o_paths_lu import _analyze_cell
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_case_LU,
)


TRACK_RXNS_THESIS = ("R5", "R6", "R7", "R8", "R9", "R10", "R11")
TRACK_SPECIES = ("CO", "CO2", "H2", "H2O", "CH4", "O2")


def main() -> int:
    case = load_case_LU()
    cfg = build_phase2_htw_lu_freeboard_reactor_config(case)
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    bed_exit_wet = result["bed_exit_gas"]
    bed_exit_dry = result["bed_exit_gas_dry"]
    exit_wet = result["exit_gas"]
    exit_dry = result["exit_gas_dry"]
    fb_diag_impl = result.get("freeboard_reaction_diag_impl", result["freeboard_reaction_diag"])
    fb_diag = result.get("freeboard_reaction_diag_thesis", fb_diag_impl)
    top_cells = []
    for i in range(max(len(reactor.cells) - 3, 0), len(reactor.cells)):
        row = _analyze_cell(reactor.cells[i])
        top_cells.append((i, row))

    print("=" * 176)
    print("LU bed-exit -> reactor-exit budget audit")
    print("=" * 176)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"secondary_agent={case['secondary_injection_agent']} "
        f"secondary_Nm3_h={case['secondary_agent_Nm3_h']:.1f} "
        f"xi_inj={case['secondary_injection_xi']:.2f} "
        f"Texit={result['reactor_exit_T']:.1f}K bed_exit_T={result['bed_T_profile'][-1]:.1f}K "
        f"n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    numbering = result.get("reaction_numbering", {})
    if numbering:
        print(
            f"reaction_label_authority={numbering.get('authority', 'unknown')} | "
            "table=freeboard_reaction_diag_thesis (impl view in *_impl)"
        )
    print("-" * 176)
    print(f"{'species':>8} {'bed_exit_wet':>12} {'reactor_exit_wet':>18} {'delta_wet':>12} {'bed_dry':>12} {'exit_dry':>12}")
    print("-" * 176)
    for sp in TRACK_SPECIES:
        print(
            f"{sp:>8s} {float(bed_exit_wet.get(sp, 0.0)):>12.5f} "
            f"{float(exit_wet.get(sp, 0.0)):>18.5f} "
            f"{float(exit_wet.get(sp, 0.0) - bed_exit_wet.get(sp, 0.0)):>12.5f} "
            f"{float(bed_exit_dry.get(sp, 0.0)):>12.5f} {float(exit_dry.get(sp, 0.0)):>12.5f}"
        )

    print("-" * 176)
    print(f"{'freeboard rxn':>14} {'sum_extent':>14} {'max_slice':>14}")
    print("-" * 176)
    for rxn in TRACK_RXNS_THESIS:
        vals = [float(v) for v in fb_diag.get(rxn, [])]
        print(f"{rxn:>14s} {sum(vals):>14.6f} {(max(vals, key=abs) if vals else 0.0):>14.6f}")

    print("-" * 176)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'H2O':>8} {'CH4':>8} "
        f"{'R4_area':>10} {'R7':>10} {'R8':>10} {'CH4_net':>10} {'H2O_net':>10}"
    )
    print("-" * 176)
    for i, row in top_cells:
        ch4_net = sum(float(v) for v in row["ch4_terms"].values())
        h2o_net = sum(float(v) for v in row["h2o_terms"].values())
        print(
            f"{i:>4d} {float(reactor.cells[i].geo.h_center/cfg.H_bed):>6.2f} {row['T']:>8.1f} "
            f"{row['wet']['CO']:>8.4f} {row['wet']['CO2']:>8.4f} {row['wet']['H2']:>8.4f} "
            f"{row['wet']['H2O']:>8.4f} {row['wet']['CH4']:>8.4f} "
            f"{row['raw']['r4_area']:>10.3f} {row['raw']['ext7']:>10.3f} {row['raw']['ext8']:>10.3f} "
            f"{ch4_net:>10.3f} {h2o_net:>10.3f}"
        )

    print("-" * 176)
    print("readout:")
    print("  - 顶床层 3 个 cell 用来判断进入 freeboard 前，CO/H2/H2O/CH4 已经被床层 chemistry 推到什么方向。")
    print("  - freeboard reaction sum 用来判断 reactor-exit 主要是延续床顶组成，还是被 freeboard 显著再分配。")
    print("  - 若 freeboard Δ 远小于 bed-top composition bias，则后续优先级应回到 bed chemistry。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
