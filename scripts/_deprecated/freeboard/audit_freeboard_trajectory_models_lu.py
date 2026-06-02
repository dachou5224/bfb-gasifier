#!/usr/bin/env python3
"""LU 工况下 freeboard trajectory model 对比审计。

用较小的对照集快速比较：
- 旧的 `surrogate_exp`
- 数值力平衡 `force_balance`
- 解析法 `analytical_wirsum`

关注点：
- entrained solids survival 是否仍被 `beta_A` 直接压没
- reactor-exit 是否几乎不受 trajectory model 影响
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
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
)


CASES = (
    ("surrogate_exp", 1.53, 1.00),
    ("surrogate_exp", 1.53, 0.50),
    ("surrogate_exp", 1.53, 0.25),
    ("force_balance", 1.53, 1.00),
    ("force_balance", 1.53, 0.50),
    ("force_balance", 1.53, 0.25),
    ("analytical_wirsum", 1.53, 1.00),
    ("analytical_wirsum", 1.53, 0.50),
    ("analytical_wirsum", 1.53, 0.25),
)


def main() -> int:
    print("=" * 152)
    print("LU freeboard trajectory-model audit")
    print("=" * 152)
    print(
        f"{'model':>14} {'u_gb':>6} {'betaA':>6} {'Texit':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'CH4':>8} "
        f"{'eject':>8} {'exit_s':>8} {'return':>8} {'survive':>8} {'n_iter':>7}"
    )
    print("-" * 152)

    for model, u_scale, beta_scale in CASES:
        cfg = build_phase2_htw_lu_freeboard_reactor_config()
        cfg.freeboard_trajectory_model = model
        cfg.freeboard_u_gb_scale = float(u_scale)
        cfg.freeboard_beta_a_scale = float(beta_scale)
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
            check_x0=True,
        )

        eject = float(result["freeboard_entrained_eject_char_ash_kg_s"])
        exit_s = float(result["freeboard_entrained_exit_char_kg_s"] + result["freeboard_entrained_exit_ash_kg_s"])
        returned = float(result["freeboard_entrained_return_char_kg_s"] + result["freeboard_entrained_return_ash_kg_s"])
        survive = exit_s / max(eject, 1e-12)
        dry = result["exit_gas_dry"]

        print(
            f"{model:>14} {u_scale:>6.2f} {beta_scale:>6.2f} "
            f"{result['reactor_exit_T']:>8.1f} "
            f"{dry['CO']:>8.4f} {dry['CO2']:>8.4f} {dry['H2']:>8.4f} {dry['CH4']:>8.4f} "
            f"{eject:>8.4f} {exit_s:>8.4f} {returned:>8.4f} {survive:>8.4f} {result['n_iter']:>7d}"
        )
        print(
            f"{'':>14} {'':>6} {'':>6} {'':>8} {'':>8} {'':>8} {'':>8} {'':>8} "
            f"{'':>8} {'':>8} {'':>8} {'':>8} {'wall':>4}={monitor['wall_time_s']:.2f}s "
            f"{'rms':>3}={monitor['rms_scaled_final']:.3e} {'x0_ok':>5}={bool((monitor.get('x0_sanity') or {}).get('ok', False))}"
        )

    print("-" * 152)
    print("notes:")
    print("  - `force_balance` 若把 survival 拉到远高于 `surrogate_exp`，说明旧指数衰减对 solids 过强。")
    print("  - `analytical_wirsum` 应在不劣化出口指标的前提下，提供更好的轨迹模块效率。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
