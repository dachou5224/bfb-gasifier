#!/usr/bin/env python3
"""LU 工况顶床层 chemistry 专项审计。"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from audit_ch4_h2o_paths_lu import _analyze_cell
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


XI_MIN = 0.65


def main() -> int:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    rows = []
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / cfg.H_bed)
        if xi < XI_MIN:
            continue
        row = _analyze_cell(cell)
        rows.append((i, xi, row))

    print("=" * 188)
    print("LU top-bed chemistry audit (thesis numbering)")
    print("=" * 188)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"xi_min={XI_MIN:.2f} Texit={result['T_profile'][-1]:.1f}K "
        f"Tbed_top={result['bed_T_profile'][-1]:.1f}K n_iter={result['n_iter']} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("labels: thesis 编号优先；R3_area(i4)=Boudouard，R9_raw(i7)=methane reforming")
    print("-" * 188)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'H2O':>8} {'CH4':>8} "
        f"{'R2_area':>10} {'R3_area(i4)':>12} {'R9_raw(i7)':>12} {'R8_raw':>10} "
        f"{'CH4_net':>10} {'H2O_net':>10} {'O2b':>8} {'O2d':>8} {'H2Olim':>8}"
    )
    print("-" * 188)

    sum_r2 = sum_r4 = sum_r7 = sum_r8 = 0.0
    sum_ch4 = sum_h2o = 0.0
    for i, xi, row in rows:
        ch4_net = sum(float(v) for v in row["ch4_terms"].values())
        h2o_net = sum(float(v) for v in row["h2o_terms"].values())
        sum_r2 += float(row["raw"]["r2_area"])
        sum_r4 += float(row["raw"]["r4_area"])
        sum_r7 += float(row["raw"]["ext7"])
        sum_r8 += float(row["raw"]["ext8"])
        sum_ch4 += ch4_net
        sum_h2o += h2o_net
        print(
            f"{i:>4d} {xi:>6.2f} {row['T']:>8.1f} "
            f"{row['wet']['CO']:>8.4f} {row['wet']['CO2']:>8.4f} {row['wet']['H2']:>8.4f} {row['wet']['H2O']:>8.4f} {row['wet']['CH4']:>8.4f} "
            f"{row['raw']['r2_area']:>10.3f} {row['raw']['r4_area']:>12.3f} {row['raw']['ext7']:>12.3f} {row['raw']['ext8']:>10.3f} "
            f"{ch4_net:>10.3f} {h2o_net:>10.3f} "
            f"{row['limit_factor_o2_bubble']:>8.3f} {row['limit_factor_o2_dense']:>8.3f} {row['limit_factor_h2o']:>8.3f}"
        )

    print("-" * 188)
    print(
        f"top-bed sums: R2_area={sum_r2:.3f} R3_area(i4)={sum_r4:.3f} "
        f"R9_raw(i7)={sum_r7:.3f} R8_raw={sum_r8:.3f} CH4_net={sum_ch4:.3f} H2O_net={sum_h2o:.3f}"
    )
    print("-" * 188)
    print("readout:")
    print("  - Hamel thesis-aligned 路径下，`R9_raw(i7)` 应保持非负；若再次出现持续负值，说明模型重新引入了逆向甲烷化语义。")
    print("  - 若 `R9_raw(i7)≈0` 且 `yCH4≈0`，说明顶床层甲烷蒸汽重整已很弱，CH4 主要由更下游/更早阶段决定。")
    print("  - 若 `R2_area` 与 `H2O_net` 绝对值持续很大，而 `R8_raw` 仅弱负，则 H2O 赤字主控仍是 char-steam gasification。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
