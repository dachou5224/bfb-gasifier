#!/usr/bin/env python3
"""审计 freeboard secondary inlet 的 pre-mix / post-mix / post-rxn 三态。"""

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
    load_case_LU,
)


REFINE = 8
MODE = "distributed_uniform"


def main() -> int:
    case = load_case_LU()
    cfg = build_phase2_htw_lu_freeboard_reactor_config(case)
    cfg.freeboard_secondary_local_refine = REFINE
    cfg.freeboard_secondary_injection_mode = MODE

    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )
    xi = result["freeboard_secondary_observation_xi"]
    if not xi:
        print("secondary observation is empty")
        return 0

    pre_t = result["freeboard_secondary_observation_pre_mix_T"]
    mix_t = result["freeboard_secondary_observation_post_mix_pre_rxn_T"]
    post_t = result["freeboard_secondary_observation_post_rxn_T"]
    pre = result["freeboard_secondary_observation_pre_mix_wet_gas"]
    mix = result["freeboard_secondary_observation_post_mix_pre_rxn_wet_gas"]
    post = result["freeboard_secondary_observation_post_rxn_wet_gas"]

    print("=" * 176)
    print("LU freeboard secondary observation audit")
    print("=" * 176)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"secondary_agent={case['secondary_injection_agent']} "
        f"secondary_Nm3_h={case['secondary_agent_Nm3_h']:.1f} "
        f"refine={REFINE} mode={MODE} "
        f"inj_seg={result['freeboard_secondary_injection_segment']}"
    )
    print(
        f"max O2 post-mix={max(mix['O2']):.5f} "
        f"max O2 post-rxn={max(post['O2']):.5f}"
    )
    print(
        f"{'xi':>8} "
        f"{'Tpre':>8} {'Tmix':>8} {'Tpost':>8} "
        f"{'O2pre':>8} {'O2mix':>8} {'O2post':>8} "
        f"{'COpre':>8} {'COmix':>8} {'COpost':>8} "
        f"{'H2pre':>8} {'H2mix':>8} {'H2post':>8}"
    )
    for row in zip(
        xi,
        pre_t,
        mix_t,
        post_t,
        pre["O2"],
        mix["O2"],
        post["O2"],
        pre["CO"],
        mix["CO"],
        post["CO"],
        pre["H2"],
        mix["H2"],
        post["H2"],
    ):
        print(
            f"{row[0]:>8.3f} "
            f"{row[1]:>8.1f} {row[2]:>8.1f} {row[3]:>8.1f} "
            f"{row[4]:>8.4f} {row[5]:>8.4f} {row[6]:>8.4f} "
            f"{row[7]:>8.4f} {row[8]:>8.4f} {row[9]:>8.4f} "
            f"{row[10]:>8.4f} {row[11]:>8.4f} {row[12]:>8.4f}"
        )
    print("-" * 176)
    print("note:")
    print("  - post-mix/pre-rxn 用来判定文献里的 O2 spike 是否只是在当前输出采样口径下被隐藏。")
    print("  - 若 O2mix > 0 但 O2post ~ 0，则说明局部快氧化存在，只是 endpoint profile 留不住。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
