#!/usr/bin/env python3
"""LU 工况下 cell0 R6 / O2 / VM overlap 灵敏度审计。

固定 shared ``global_nr`` 口径，比较以下旋钮对出口与 cell0 局部状态的影响：

- `r6_scale`
- `gas_inlet_dense_frac`
- 两者联动

目的：
- 区分 `CH4≈0` 更像是 `R6 kinetics` 问题，还是底部 `O2` 相分配问题
- 量化每个旋钮对 `CO/CO2/H2/CH4/Xc/Texit` 的联动影响
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


CASES = [
    ("baseline_dense030", 1.0, 0.30),
    ("dense020", 1.0, 0.20),
    ("dense040", 1.0, 0.40),
    ("r6_half_dense030", 0.50, 0.30),
    ("r6_half_dense040", 0.50, 0.40),
    ("r6_quarter_dense040", 0.25, 0.40),
]


def _run_case(name: str, r6_scale: float, dense: float) -> dict[str, float | str | int]:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    cfg.r6_scale = float(r6_scale)
    cfg.gas_inlet_dense_frac = float(dense)

    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )
    dry = result["exit_gas_dry"]
    cell0 = reactor.cells[0]
    wet0 = cell0._mole_fractions("combined")

    return {
        "name": name,
        "r6_scale": float(r6_scale),
        "dense_frac": float(dense),
        "Texit": float(result["T_profile"][-1]),
        "T0": float(cell0.T),
        "Xc": float(result["carbon_conv"]),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "n_iter": int(result.get("n_iter", 0)),
        "wall_time_s": float(monitor.get("wall_time_s", float("nan"))),
        "x0_ok": bool((monitor.get("x0_sanity") or {}).get("ok", False)),
        "CO": float(dry["CO"]),
        "CO2": float(dry["CO2"]),
        "H2": float(dry["H2"]),
        "CH4": float(dry["CH4"]),
        "cell0_O2_wet": float(wet0[idx["O2"]]),
        "cell0_CH4_wet": float(wet0[idx["CH4"]]),
        "cell0_H2O_wet": float(wet0[idx["H2O"]]),
    }


def main() -> int:
    rows = [_run_case(*case) for case in CASES]

    print("=" * 172)
    print("LU cell0 overlap sensitivity audit (shared global NR)")
    print("=" * 172)
    print(
        f"{'case':<20} {'r6':>5} {'dense':>6} {'Texit':>8} {'T0':>8} {'Xc':>8} {'rms':>10} {'iter':>5} {'wall':>7} {'x0_ok':>7} "
        f"{'CO':>8} {'CO2':>8} {'H2':>8} {'CH4':>8} {'O2@0':>8} {'CH4@0':>8} {'H2O@0':>8}"
    )
    print("-" * 172)
    for row in rows:
        print(
            f"{str(row['name']):<20} {row['r6_scale']:>5.2f} {row['dense_frac']:>6.2f} "
            f"{row['Texit']:>8.1f} {row['T0']:>8.1f} {row['Xc']:>8.4f} {row['rms']:>10.3e} {row['n_iter']:>5d} "
            f"{row['wall_time_s']:>7.2f} {str(row['x0_ok']):>7} "
            f"{row['CO']:>8.4f} {row['CO2']:>8.4f} {row['H2']:>8.4f} {row['CH4']:>8.4f} "
            f"{row['cell0_O2_wet']:>8.4f} {row['cell0_CH4_wet']:>8.4f} {row['cell0_H2O_wet']:>8.4f}"
        )

    print("-" * 172)
    print("Notes:")
    print("  - 若降低 `dense_frac` 比降低 `r6_scale` 更能改善 `CO/CO2`，说明主控旋钮更像底部 O2 相分配。")
    print("  - 若 `cell0_CH4_wet` 始终接近 0，则 CH4 问题不只是出口干基定义，而是底部源汇重叠。")
    print("  - 该脚本只做开发审计，不直接改 shared baseline。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
