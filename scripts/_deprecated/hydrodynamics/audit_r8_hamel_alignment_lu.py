#!/usr/bin/env python3
"""LU 工况下 R8 实现与 Hamel 口径对齐审计。

固定 shared ``global_nr`` 口径，定量比较以下 R8 变体：

1. 当前代码：`P_bar` + `a_R8 = 1.0`
2. 仅改催化因子：`P_bar` + `a_R8 = 0.02`
3. 仅改压力基准：`P_atm` + `a_R8 = 1.0`
4. 两者同时改：`P_atm` + `a_R8 = 0.02`

目的不是直接改主模型，而是先回答：
- 当前偏差主要来自 `atm/bar` 还是 `a_R8`
- 在 LU / 25 bar 工况下，这两类口径差异的量级各是多少
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.kinetics.arrhenius import k_standard
from src.kinetics.gas_reactions import R8_E_Rg, R8_k0, wgsr_equilibrium_constant
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


def _rate_variant(
    *,
    T: float,
    P: float,
    y_co: float,
    y_h2o: float,
    y_co2: float,
    y_h2: float,
    pressure_basis: str,
    catalyst_factor: float,
) -> float:
    T_safe = max(float(T), 300.0)
    P_safe = float(P)
    c_fac = P_safe / (Rg * T_safe)
    c_co = max(float(y_co), 0.0) * c_fac
    c_h2o = max(float(y_h2o), 0.0) * c_fac
    c_co2 = max(float(y_co2), 0.0) * c_fac
    c_h2 = max(float(y_h2), 0.0) * c_fac

    if pressure_basis == "bar":
        p_dimless = P_safe / 1e5
    elif pressure_basis == "atm":
        p_dimless = P_safe / 101325.0
    else:
        raise ValueError(f"Unsupported pressure_basis={pressure_basis!r}")

    k_std = k_standard(R8_k0, R8_E_Rg * Rg, T_safe)
    exp_corr = math.exp(max(min(-8.91 + 5.553 / T_safe, 100.0), -100.0))
    exponent = max(0.5 - p_dimless / 250.0, 0.01)
    k_fwd = float(catalyst_factor) * k_std * (p_dimless ** exponent) * exp_corr
    k_eq = wgsr_equilibrium_constant(T_safe)
    return k_fwd * (c_co * c_h2o - c_co2 * c_h2 / max(k_eq, 1e-300))


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    totals = {
        "cur_bar_a1": 0.0,
        "bar_a002": 0.0,
        "atm_a1": 0.0,
        "atm_a002": 0.0,
    }

    print("=" * 166)
    print("LU R8 Hamel-alignment audit (shared global NR)")
    print("=" * 166)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"init={result.get('nr_init_strategy')} jacobian={result.get('nr_jacobian_strategy')} "
        f"converged={result.get('converged')} n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"Texit={result['T_profile'][-1]:.1f}K"
    )
    print("-" * 166)
    print(
        f"{'cell':>4} {'xi':>5} {'T[K]':>8} {'CO':>7} {'CO2':>7} {'H2':>7} {'H2O':>7} "
        f"{'cur':>11} {'a=0.02':>11} {'atm':>11} {'atm+a=0.02':>13} "
        f"{'a002/cur':>10} {'atm/cur':>10}"
    )
    print("-" * 166)

    for i, cell in enumerate(reactor.cells):
        y_d = cell._mole_fractions("d")
        T = float(cell.T)
        P = float(cell.P)
        xi = float(cell.geo.h_center / cfg.H_bed)
        y_co = float(y_d[idx["CO"]])
        y_co2 = float(y_d[idx["CO2"]])
        y_h2 = float(y_d[idx["H2"]])
        y_h2o = float(y_d[idx["H2O"]])

        r_cur = _rate_variant(
            T=T, P=P, y_co=y_co, y_h2o=y_h2o, y_co2=y_co2, y_h2=y_h2,
            pressure_basis="bar", catalyst_factor=1.0,
        )
        r_a002 = _rate_variant(
            T=T, P=P, y_co=y_co, y_h2o=y_h2o, y_co2=y_co2, y_h2=y_h2,
            pressure_basis="bar", catalyst_factor=0.02,
        )
        r_atm = _rate_variant(
            T=T, P=P, y_co=y_co, y_h2o=y_h2o, y_co2=y_co2, y_h2=y_h2,
            pressure_basis="atm", catalyst_factor=1.0,
        )
        r_atm_a002 = _rate_variant(
            T=T, P=P, y_co=y_co, y_h2o=y_h2o, y_co2=y_co2, y_h2=y_h2,
            pressure_basis="atm", catalyst_factor=0.02,
        )

        ext_cur = r_cur * cell.V_d
        ext_a002 = r_a002 * cell.V_d
        ext_atm = r_atm * cell.V_d
        ext_atm_a002 = r_atm_a002 * cell.V_d

        totals["cur_bar_a1"] += ext_cur
        totals["bar_a002"] += ext_a002
        totals["atm_a1"] += ext_atm
        totals["atm_a002"] += ext_atm_a002

        ratio_a002 = ext_a002 / max(abs(ext_cur), 1e-300)
        ratio_atm = ext_atm / max(abs(ext_cur), 1e-300)

        print(
            f"{i:>4d} {xi:>5.2f} {T:>8.1f} {y_co:>7.3f} {y_co2:>7.3f} {y_h2:>7.3f} {y_h2o:>7.3f} "
            f"{ext_cur:>11.3e} {ext_a002:>11.3e} {ext_atm:>11.3e} {ext_atm_a002:>13.3e} "
            f"{ratio_a002:>10.3e} {ratio_atm:>10.3e}"
        )

    print("-" * 166)
    base = max(abs(totals["cur_bar_a1"]), 1e-300)
    print("Integrated dense-phase cell sources [mol/s-equivalent over V_d]:")
    print(f"  current      (bar, a=1.0 ) : {totals['cur_bar_a1']:.6e}")
    print(f"  Hamel-a only (bar, a=0.02): {totals['bar_a002']:.6e}")
    print(f"  atm only     (atm, a=1.0 ) : {totals['atm_a1']:.6e}")
    print(f"  both         (atm, a=0.02): {totals['atm_a002']:.6e}")
    print("Ratios vs current:")
    print(f"  a=0.02 / current     = {totals['bar_a002'] / base:.6e}")
    print(f"  atm / current        = {totals['atm_a1'] / base:.6e}")
    print(f"  atm+a=0.02 / current = {totals['atm_a002'] / base:.6e}")
    print("Notes:")
    print("  - 若 `atm/current` 接近 1，则 `atm/bar` 只是文档一致性问题，不是当前主偏差来源。")
    print("  - 若 `a=0.02/current` 远小于 1，则 `R8_a_R8` 才是更可能推高 WGSR 强度的主因。")
    print("  - 本脚本只比较当前解上的局部速率量级，不直接改变整炉解。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
