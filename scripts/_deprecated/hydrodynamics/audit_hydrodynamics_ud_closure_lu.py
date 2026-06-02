#!/usr/bin/env python3
"""LU 工况 `u_d` closure 灵敏度审计。"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from src.core.reactor import Reactor
from _nr_monitor import solve_with_nr_monitor
from audit_ch4_h2o_paths_lu import _analyze_cell
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_case_LU,
)


_CLOSURES = (
    "current",
    "umf_over_epsmf",
    "hilligardt_eq311",
    "wein_1992_eq312",
    "backsolve_visible_epsb",
)
_XI_TOP_MIN = 0.65


def _safe_float(val: object) -> float:
    try:
        out = float(val)
    except (TypeError, ValueError):
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def main() -> int:
    case = load_case_LU()

    print("=" * 212)
    print("LU u_d-closure hydrodynamics sensitivity audit")
    print("=" * 212)
    print(
        "closures: current | umf_over_epsmf | hilligardt_eq311 | wein_1992_eq312 | backsolve_visible_epsb\n"
        "target: quantify whether the current solved-state `u_d ~= u0` closure is the first-order cause of "
        "`eps_b -> 0.01 floor` and `eps_d_voidage -> 0.99 clip`, and how strongly top-bed chemistry follows it."
    )
    print("-" * 212)
    print(
        f"{'closure':<24} {'Texit[K]':>9} {'n_iter':>7} {'rms':>10} "
        f"{'x0_ok':>6} "
        f"{'eps_b top':>10} {'eps_d top':>10} {'solid top':>10} {'u_d/u_mf top':>13} "
        f"{'R2 top':>11} {'R7 top':>11} {'R8 top':>11} {'CH4 top':>11} {'H2O top':>11} "
        f"{'COdry':>8} {'CO2dry':>8} {'H2dry':>8} {'CH4dry':>8}"
    )
    print("-" * 212)

    for closure in _CLOSURES:
        cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
        cfg.hydrodynamics_u_d_closure = closure
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
            check_x0=True,
        )

        top_cells = [c for c in reactor.cells if float(c.geo.h_center / cfg.H_bed) >= _XI_TOP_MIN]
        eps_b_top = float(sum(float(c.eps_b) for c in top_cells) / max(len(top_cells), 1))
        eps_d_top = float(sum(float(c.eps_d_voidage) for c in top_cells) / max(len(top_cells), 1))
        solid_top = float(sum(1.0 - float(c.eps_d_voidage) for c in top_cells) / max(len(top_cells), 1))
        ud_ratio_top = float(
            sum(float(c.u_d / max(c.u_mf, 1e-12)) for c in top_cells) / max(len(top_cells), 1)
        )

        sum_r2 = sum_r7 = sum_r8 = sum_ch4 = sum_h2o = 0.0
        for cell in top_cells:
            row = _analyze_cell(cell)
            sum_r2 += float(row["raw"]["r2_area"])
            sum_r7 += float(row["raw"]["ext7"])
            sum_r8 += float(row["raw"]["ext8"])
            sum_ch4 += sum(float(v) for v in row["ch4_terms"].values())
            sum_h2o += sum(float(v) for v in row["h2o_terms"].values())

        dry = result.get("exit_gas_dry", {})
        print(
            f"{closure:<24} {result['T_profile'][-1]:>9.1f} {int(result['n_iter']):>7d} "
            f"{_safe_float(result.get('rms_scaled_final')):>10.3e} "
            f"{str(bool(monitor.get('x0_sanity', {}).get('ok', False))):>6} "
            f"{eps_b_top:>10.4f} {eps_d_top:>10.4f} {solid_top:>10.4f} {ud_ratio_top:>13.3f} "
            f"{sum_r2:>11.3f} {sum_r7:>11.3f} {sum_r8:>11.3f} {sum_ch4:>11.3f} {sum_h2o:>11.3f} "
            f"{_safe_float(dry.get('CO')):>8.4f} {_safe_float(dry.get('CO2')):>8.4f} "
            f"{_safe_float(dry.get('H2')):>8.4f} {_safe_float(dry.get('CH4')):>8.4f}"
        )

    print("-" * 212)
    print("readout:")
    print("  - 若替代 closure 仅通过修正 `u_d` 就能把 top-bed `eps_d_voidage` 从 ~0.99 拉回 0.45-0.60，说明当前主偏差首先来自 hydrodynamics。")
    print("  - 若 `R2/R7/R8/CH4/H2O` 随 closure 明显漂移，则 bed chemistry 结果不能在当前 `u_d` 口径下直接解读。")
    print("  - `hilligardt_eq311` 与 `wein_1992_eq312` 才是当前最接近 Hamel 原文的 `u_d` audit branches。")
    print("  - `backsolve_visible_epsb` 并不代表已经锁定原文公式；它只是把当前两条 Hamel 提取式 (`eps_b` 与 visible bubble) 组合成一个显式 closure。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
