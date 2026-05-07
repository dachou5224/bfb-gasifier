#!/usr/bin/env python3
"""LU 工况关键反应进程审计。

固定 shared ``global_nr`` 口径，沿 reactor 轴向追踪以下关键反应。
输出默认使用 **Hamel (1999) thesis 编号**：

- `R3` Boudouard（当前实现入口 `rate_R4_effective`）
- `R8` Water-gas shift
- `R9` methane reforming（当前实现入口 `rate_R7`）

目标：

1. 审计当前 ``T`` / 组成下的动力学量级；
2. 对 `R8` 追踪 `K_eq`, `Q_p`, driving force；
3. 对 `R3` 追踪 Langmuir-Hinshelwood 分母项与 `CO2/CO` 压力代理，
   明确它当前是“动力学受限”还是“组成受限”。

说明：
- 当前代码库**没有** `R3`（Boudouard）的显式热力学 `K_eq` 实现，因此脚本只输出 `Q_proxy`
  与动力学分母项，不伪造 `R3` 的平衡判据。
- 当前 thesis `R9` 在代码中是**单向 reforming**（实现入口 `rate_R7`），因此脚本不再伪造 reversible split。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import Rg
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from src.core.species import gas_diffusivity_correlation
from src.core.species import GAS_SPECIES_INDEX as idx
from src.kinetics.arrhenius import k_standard
from src.kinetics.char_reactions import (
    R4_kb_E,
    R4_kb_k0,
    R4_kc_E,
    R4_kc_k0,
    R4_kf_E,
    R4_kf_k0,
    d_core_from_spm_char_conversion,
    rate_R4,
    rate_R4_effective,
)
from src.kinetics.gas_reactions import (
    R8_E_Rg,
    R8_a_R8,
    R8_k0,
    wgsr_equilibrium_constant,
    rate_R7,
    rate_R8,
)
from src.thermodynamics.equilibrium import (
    calc_reaction_quotient,
)
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


def _safe_div(num: float, den: float) -> float:
    return float(num / den) if abs(den) > 1e-300 else float("inf")


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    print("=" * 160)
    print("LU reaction progress audit (shared global NR)")
    print("=" * 160)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"init={result.get('nr_init_strategy')} jacobian={result.get('nr_jacobian_strategy')} "
        f"converged={result.get('converged')} outer_iters={result.get('nr_outer_iters')} "
        f"n_iter={result.get('n_iter')} rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"Texit={result['T_profile'][-1]:.1f}K"
    )
    print("-" * 160)
    print(
        f"{'cell':>4} {'xi':>5} {'T[K]':>8} "
        f"{'CO':>7} {'CO2':>7} {'H2':>7} {'H2O':>7} {'CH4':>7} "
        f"{'R3eff':>10} {'Q3proxy':>10} {'denLH':>10} "
        f"{'R9raw(i7)':>12} "
        f"{'K8':>10} {'Q8':>10} {'drv8':>9} {'r8':>10}"
    )
    print("-" * 160)

    for i, cell in enumerate(reactor.cells):
        cell.calc_hydrodynamics()
        tau = cell.geo.dh / max(cell.u_mf, 1e-3)
        cell.compute_vorabrechnung(tau)
        cell.calc_exchange()

        y_d = cell._mole_fractions("d")
        C_d = cell._concentrations("d")

        T = float(cell.T)
        xi = float(cell.geo.h_center / cfg.H_bed)

        p_co2 = max(float(C_d[idx["CO2"]] * Rg * T), 0.0)
        p_co = max(float(C_d[idx["CO"]] * Rg * T), 0.0)
        D_g = gas_diffusivity_correlation(T, cell.P)
        d_core = d_core_from_spm_char_conversion(
            float(np.clip(1.0 - np.sum(cell.m_solid[:, 0]) / max(np.sum(cell.m_solid_zu[:, 0] + cell.m_solid_in[:, 0] + cell.m_solid_rez[:, 0]), 1e-12), 0.0, 1.0 - 1e-9)),
            float(cell.solid.d_p),
        )
        r4_kin = rate_R4(T, p_co2, p_co)
        r4_eff = rate_R4_effective(T, p_co2, p_co, float(cell.solid.d_p), D_g, d_core)
        kb = k_standard(R4_kb_k0, R4_kb_E, T)
        kc = k_standard(R4_kc_k0, R4_kc_E, T)
        q4_proxy = _safe_div(p_co**2, max(p_co2, 1e-30))
        den_lh = 1.0 + kb * p_co2 + kc * p_co

        C_CH4 = max(float(C_d[idx["CH4"]]), 0.0)
        C_H2O = max(float(C_d[idx["H2O"]]), 0.0)
        C_CO = max(float(C_d[idx["CO"]]), 0.0)
        C_H2 = max(float(C_d[idx["H2"]]), 0.0)
        r9 = rate_R7(T, C_CH4, C_H2O, C_CO, C_H2)

        K8 = wgsr_equilibrium_constant(T)
        Q8 = calc_reaction_quotient("R8", y_d, cell.P, idx)
        y_co = max(float(y_d[idx["CO"]]), 0.0)
        y_h2o = max(float(y_d[idx["H2O"]]), 0.0)
        y_co2 = max(float(y_d[idx["CO2"]]), 0.0)
        y_h2 = max(float(y_d[idx["H2"]]), 0.0)
        drv8 = 1.0 - (Q8 / max(K8, 1e-300))
        P_bar = cell.P / 1e5
        k8_std = k_standard(R8_k0, R8_E_Rg * Rg, max(T, 300.0))
        exp_corr = np.exp(np.clip(-8.91 + 5.553 / max(T, 300.0), -100.0, 100.0))
        P_factor = P_bar ** max(0.5 - P_bar / 250.0, 0.01)
        k8_fwd = R8_a_R8 * k8_std * P_factor * exp_corr
        c_fac = cell.P / (Rg * max(T, 300.0))
        c_co = y_co * c_fac
        c_h2o = y_h2o * c_fac
        c_co2 = y_co2 * c_fac
        c_h2 = y_h2 * c_fac
        r8_fwd = k8_fwd * c_co * c_h2o
        r8_rev = k8_fwd * c_co2 * c_h2 / max(K8, 1e-300)
        r8 = rate_R8(
            T,
            cell.P,
            float(y_d[idx["CO"]]),
            float(y_d[idx["H2O"]]),
            float(y_d[idx["CO2"]]),
            float(y_d[idx["H2"]]),
        )

        print(
            f"{i:>4d} {xi:>5.2f} {T:>8.1f} "
            f"{float(y_d[idx['CO']]):>7.3f} {float(y_d[idx['CO2']]):>7.3f} "
            f"{float(y_d[idx['H2']]):>7.3f} {float(y_d[idx['H2O']]):>7.3f} {float(y_d[idx['CH4']]):>7.4f} "
            f"{r4_eff:>10.3e} {q4_proxy:>10.3e} {den_lh:>10.3e} "
            f"{r9:>12.3e} "
            f"{K8:>10.3e} {Q8:>10.3e} {drv8:>9.2e} {r8:>10.3e}"
        )

        print(
            f"     detail: R3[k_kin={r4_kin:.3e}, Dg={D_g:.3e}, kb*PCO2={kb*p_co2:.3e}, kc*PCO={kc*p_co:.3e}] "
            f"R9(raw={r9:.3e}) "
            f"R8[fwd={r8_fwd:.3e}, rev={r8_rev:.3e}]"
        )

    print("-" * 160)
    print("Notes:")
    print("  - R9(raw) 指 thesis R9 methane reforming，对应当前实现入口 `rate_R7`。")
    print("  - R8: drv = 1 - Qp/Keq; drv<0 表示当前局部组成推动逆向。")
    print("  - R3: 当前代码库未实现显式热力学 Keq，本脚本仅输出 Q3proxy=P_CO^2/P_CO2 与 LH 分母项。")
    print("  - 若上部 cell 持续出现 R8 正向且 H2O 大幅为负源项，应优先审 H2O deficit 而非热峰。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
