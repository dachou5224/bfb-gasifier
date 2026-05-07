#!/usr/bin/env python3
"""shared global NR 口径下 LU 工况前几格 O2 追踪。

复用 ``audit_oxygen_reaction_trace.py`` 的单 cell O2 分解逻辑，
但固定到当前 shared ``global_nr`` 开发口径，重点看 `cell0-2`：

- bubble / dense O2 供给
- bubble / dense O2 需求
- limiter 与 shortfall
- bubble->dense 的 O2 transfer potential
- VM 可氧化 O2 当量
- CH4 净代理项
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from audit_oxygen_reaction_trace import _cell_o2_trace
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    print("=" * 180)
    print("LU shared global NR O2 trace")
    print("=" * 180)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"init={result.get('nr_init_strategy')} jacobian={result.get('nr_jacobian_strategy')} "
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("-" * 180)
    print(
        f"{'cell':>4} {'T[K]':>8} {'O2sup_b':>9} {'O2sup_d':>9} {'O2dem_b':>9} {'O2dem_d':>9} "
        f"{'lim_b':>7} {'lim_d':>7} {'short_b':>9} {'short_d':>9} {'xferPot':>9} "
        f"{'bubble%':>8} {'vmO2eq':>9} {'CH4net':>9}"
    )
    print("-" * 180)
    for i, cell in enumerate(reactor.cells[:3]):
        cell.calc_hydrodynamics()
        cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
        cell.calc_exchange()
        row = _cell_o2_trace(cell)
        print(
            f"{i:>4d} {row['T']:>8.1f} {row['O2_supply_bubble']:>9.3f} {row['O2_supply_dense']:>9.3f} "
            f"{row['bubble_o2_unlimited']:>9.3f} {row['dense_o2_unlimited']:>9.3f} "
            f"{row['limit_factor_o2_bubble']:>7.3f} {row['limit_factor_o2_dense']:>7.3f} "
            f"{row['O2_shortfall_bubble']:>9.3f} {row['O2_shortfall_dense']:>9.3f} {row['O2_transfer_potential_bd']:>9.3f} "
            f"{100.0 * row['bubble_share_of_effective_o2']:>7.1f}% {row['vm_oxidizable_o2eq']:>9.3f} {row['ch4_net_proxy']:>9.3f}"
        )
    print("-" * 180)
    print("Interpretation:")
    print("  - 若 bubble shortfall≈0 且 dense shortfall 很大，优先看 dense O2 budget / N_ex，而不是 bubble O2 limiter。")
    print("  - 若 cell0 的 vmO2eq 已不小且 CH4net 仍强负，说明问题更像底部局部氧化重叠，而不是上部 reforming。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
