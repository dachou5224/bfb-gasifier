#!/usr/bin/env python3
"""LU 工况下 freeboard hydrodynamics 灵敏度审计。

重点检查：
1. `u_gb/u_p0` 标度
2. `beta_A` 衰减强度

输出关注：
- entrained solids 的 freeboard 末端存活量
- cyclone recycle candidate
- reactor-exit 温度与干基主气相
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_validation_case_node,
)


TRAJECTORY_MODELS = ("surrogate_exp", "force_balance", "analytical_wirsum")
U_GB_SCALES = (1.00, 1.53, 2.00)
BETA_A_SCALES = (0.25, 0.50, 1.00, 1.50)


def main() -> int:
    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    ref_dry = ref["exit_gas_dry_mol_frac"]

    print("=" * 184)
    print("LU freeboard hydrodynamics sensitivity audit")
    print("=" * 184)
    print(
        f"{'model':>14} {'u_gb':>6} {'betaA':>6} {'Texit':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'CH4':>8} "
        f"{'eject':>8} {'exit_s':>8} {'return':>8} {'recycle':>8} {'survive':>8} {'n_iter':>7} {'wall[s]':>8}"
    )
    print("-" * 184)

    for model in TRAJECTORY_MODELS:
        for u_scale in U_GB_SCALES:
            for beta_scale in BETA_A_SCALES:
                cfg = build_phase2_htw_lu_freeboard_reactor_config()
                cfg.freeboard_trajectory_model = model
                cfg.freeboard_u_gb_scale = float(u_scale)
                cfg.freeboard_beta_a_scale = float(beta_scale)
                reactor = Reactor(cfg)
                result, monitor = solve_with_nr_monitor(
                    reactor,
                    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
                    check_x0=True,
                )

                eject = float(result["freeboard_entrained_eject_char_ash_kg_s"])
                exit_s = float(result["freeboard_entrained_exit_char_kg_s"] + result["freeboard_entrained_exit_ash_kg_s"])
                returned = float(result["freeboard_entrained_return_char_kg_s"] + result["freeboard_entrained_return_ash_kg_s"])
                recycle = float(result["freeboard_cyclone_recycle_candidate_char_ash_kg_s"])
                survive = exit_s / max(eject, 1e-12)
                dry = result["exit_gas_dry"]

                print(
                    f"{model:>14} {u_scale:>6.2f} {beta_scale:>6.2f} "
                    f"{result['reactor_exit_T']:>8.1f} "
                    f"{dry['CO']:>8.4f} {dry['CO2']:>8.4f} {dry['H2']:>8.4f} {dry['CH4']:>8.4f} "
                    f"{eject:>8.4f} {exit_s:>8.4f} {returned:>8.4f} {recycle:>8.4f} {survive:>8.4f} "
                    f"{result['n_iter']:>7d} {float(monitor.get('wall_time_s', 0.0)):>8.2f}"
                )
            print("-" * 184)
        print("-" * 184)
    print(
        f"{'ref':>14} {'':>6} {'':>6} "
        f"{float(ref['exit_temperature_K']):>8.1f} "
        f"{float(ref_dry['CO']):>8.4f} {float(ref_dry['CO2']):>8.4f} "
        f"{float(ref_dry['H2']):>8.4f} {float(ref_dry['CH4']):>8.4f}"
    )
    print("-" * 184)

    print("notes:")
    print("  - `survive = entrained_exit / entrained_eject`。")
    print("  - 若 `beta_A` 缩小后 `survive` 大幅回升，而出口气变化不大，则当前主要问题是 freeboard 固体衰减过快。")
    print("  - 若 `u_gb_scale` 改动也强烈影响 `Texit`/气相，则说明 freeboard 停留时间本身在主控 reactor-exit。")
    print("  - `force_balance` 用 Haider-Levenspiel 阻力 + 重力/浮力；若其 survival 显著高于 `surrogate_exp`，说明旧指数衰减对 solids 过强。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
