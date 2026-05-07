"""Oxygen-reaction trace audit for temperature-peak diagnosis.

Purpose
-------
For selected tuning scenarios, decompose per-cell O2-related reactions and
compare them against temperature peaks.

Focus:
- O2 supply per cell
- bubble-phase O2 demand (`R5b`, thesis `R7i6b`) and limiter participation
- suspension/char/tar O2 demand before and after O2 limiting
- VM/挥发分释放带来的可氧化组分与其 O2 当量
- all fast homogeneous oxidation reactions as one bucket
- peak-cell neighborhood analysis

Usage
-----
    .venv/bin/python scripts/audit_oxygen_reaction_trace.py
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import src.kinetics.char_reactions as cr
import src.kinetics.gas_reactions as gr
from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as IDX
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension, rate_R6, rate_R7, rate_R9, rate_R12
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10, rate_R11_bubble, rate_R11_suspension
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_refined_config, load_case_LU


@dataclass
class Scenario:
    name: str
    dense: float
    heat_loss: float
    r5s: float = 1.0
    r6s: float = 1.0
    r4s: float = 1.0


@contextmanager
def kinetics_scaling(r5_scale: float, r6_scale: float, r4_scale: float):
    b_r5 = gr.R5_suspension_k_T05
    b_r6 = gr.R6_k0
    b_r4 = cr.R4_kf_k0
    try:
        gr.R5_suspension_k_T05 = b_r5 * r5_scale
        gr.R6_k0 = b_r6 * r6_scale
        cr.R4_kf_k0 = b_r4 * r4_scale
        yield
    finally:
        gr.R5_suspension_k_T05 = b_r5
        gr.R6_k0 = b_r6
        cr.R4_kf_k0 = b_r4


def _cell_o2_trace(cell) -> dict:
    idx = IDX
    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    D_g = gr.P0 * 0.0  # placeholder to avoid lint in script-only file
    from src.core.species import gas_diffusivity_correlation

    D_g = gas_diffusivity_correlation(T, P)

    # bubble O2 consumers
    r5b = rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b
    r6b = rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]]) * cell.V_b
    r12b = rate_R12(T, C_b[idx["H2"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b if cell.enable_r12 else 0.0

    # dense/suspension O2 consumers
    r5d = rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * cell.V_d
    r6d = rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]]) * cell.V_d
    r12d = rate_R12(T, C_d[idx["H2"]], C_d[idx["O2"]], P, y_d, idx) * cell.V_d if cell.enable_r12 else 0.0
    ext7 = rate_R7(T, C_d[idx["CH4"]], C_d[idx["H2O"]], C_d[idx["CO"]], C_d[idx["H2"]]) * cell.V_d

    ext9 = 0.0
    if not cell.use_gibbs_minor:
        ext9 = rate_R9(T, C_d[idx["H2S"]], C_d[idx["O2"]]) * cell.V_d

    C_tar_b = max(C_b[idx["TAR1"]] + C_b[idx["TAR2"]], 0.0)
    C_tar_d = max(C_d[idx["TAR1"]] + C_d[idx["TAR2"]], 0.0)
    ext10b = rate_R10(T, C_tar_b, C_b[idx["O2"]], P, cell.fuel_type) * cell.V_b
    ext10d = rate_R10(T, C_tar_d, C_d[idx["O2"]], P, cell.fuel_type) * cell.V_d
    rho_cat = cell._catalyst_bulk_density()
    ext11b = rate_R11_bubble(T, C_tar_b) * cell.V_b
    ext11d = rate_R11_suspension(T, C_tar_d, rho_cat) * cell.V_d

    stoich_r10 = get_lumped_tar_stoichiometry("R10", cell.fuel_type)
    stoich_r11 = get_lumped_tar_stoichiometry("R11", cell.fuel_type)
    nu_o2_r10 = max(0.0, -float(stoich_r10.get("O2", 0.0)))
    nu_h2o_r11 = max(0.0, -float(stoich_r11.get("H2O", 0.0)))

    # VM/pyrolysis gas release currently cached by compute_vorabrechnung.
    gas_src_vm = np.maximum(cell._vm_gas_source_cache, 0.0)
    vm_co = float(gas_src_vm[idx["CO"]])
    vm_ch4 = float(gas_src_vm[idx["CH4"]])
    vm_h2 = float(gas_src_vm[idx["H2"]])
    vm_h2s = float(gas_src_vm[idx["H2S"]])
    vm_tar = float(gas_src_vm[idx["TAR1"]] + gas_src_vm[idx["TAR2"]])
    # O2-equivalent if freshly released volatiles are rapidly oxidized by the implemented oxidation set.
    vm_oxidizable_o2eq = (
        0.5 * vm_co
        + 1.5 * vm_ch4
        + 0.5 * vm_h2
        + 1.5 * vm_h2s
        + nu_o2_r10 * vm_tar
    )

    areas = cell._calc_char_surface_area_per_class()
    total_area = float(np.sum(areas))
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    alpha = 1.0
    if total_area > 0:
        X = cell._compute_char_conversion()
        dc = d_core_from_spm_char_conversion(X, cell.solid.d_p)
        r1, alpha = rate_R1(T, C_d[idx["O2"]], cell.solid.d_p, D_g, dc, cell.fuel_type)
        r2 = rate_R2(T, C_d[idx["H2O"]], cell.solid.d_p, D_g, dc)
        r3 = rate_R3(T, C_d[idx["H2"]], cell.solid.d_p, D_g, dc)
        p_co2 = max(C_d[idx["CO2"]] * Rg * T, 0.0)
        p_co = max(C_d[idx["CO"]] * Rg * T, 0.0)
        r4 = rate_R4_effective(T, p_co2, p_co, cell.solid.d_p, D_g, dc)

    o2_supply_bubble = max(
        cell.N_zu_b[idx["O2"]] + cell.N_b_in[idx["O2"]] + max(-cell.N_ex[idx["O2"]], 0.0) + cell.N_rez_b[idx["O2"]],
        1e-12,
    )
    o2_supply_dense = max(
        cell.N_zu_d[idx["O2"]] + cell.N_d_in[idx["O2"]] + max(cell.N_ex[idx["O2"]], 0.0) + cell.N_rez_d[idx["O2"]],
        1e-12,
    )

    bubble_o2_demand_unlimited = r5b + 1.5 * r6b + 0.5 * r12b + nu_o2_r10 * ext10b
    dense_o2_demand_unlimited = (
        r5d + 1.5 * r6d + 0.5 * r12d + alpha * r1 * total_area + 1.5 * ext9 + nu_o2_r10 * ext10d
    )
    o2_budget_bubble = 0.995 * o2_supply_bubble
    o2_budget_dense = 0.995 * o2_supply_dense
    limit_factor_o2_bubble = min(1.0, o2_budget_bubble / max(bubble_o2_demand_unlimited, 1e-9))
    limit_factor_o2_dense = min(1.0, o2_budget_dense / max(dense_o2_demand_unlimited, 1e-9))
    limit_factor_o2 = min(limit_factor_o2_bubble, limit_factor_o2_dense)

    bubble_o2_demand_limited = bubble_o2_demand_unlimited * limit_factor_o2_bubble
    tar_o2_demand_total = nu_o2_r10 * (ext10b + ext10d)
    dense_o2_demand_limited = dense_o2_demand_unlimited * limit_factor_o2_dense
    o2_slack_bubble = max(o2_budget_bubble - bubble_o2_demand_limited, 0.0)
    o2_slack_dense = max(o2_budget_dense - dense_o2_demand_limited, 0.0)
    o2_shortfall_bubble = max(bubble_o2_demand_unlimited - o2_budget_bubble, 0.0)
    o2_shortfall_dense = max(dense_o2_demand_unlimited - o2_budget_dense, 0.0)
    o2_transfer_potential_bd = min(o2_slack_bubble, o2_shortfall_dense)
    total_o2_demand_effective = bubble_o2_demand_limited + dense_o2_demand_limited

    # Fast homogeneous oxidation bucket: gas-phase + tar oxidation, excluding char combustion R1.
    homogeneous_o2_unlimited = bubble_o2_demand_unlimited + (r5d + 1.5 * r6d + 0.5 * r12d + 1.5 * ext9 + tar_o2_demand_total)
    homogeneous_o2_limited_effective = bubble_o2_demand_limited + (
        (r5d + 1.5 * r6d + 0.5 * r12d + 1.5 * ext9 + tar_o2_demand_total) * limit_factor_o2_dense
    )
    homogeneous_share_of_effective_o2 = homogeneous_o2_limited_effective / max(total_o2_demand_effective, 1e-12)

    # Endothermic group proxy (not O2-related but helpful for peak-cause interpretation)
    h2o_supply = max(
        cell.N_zu_d[idx["H2O"]] + cell.N_d_in[idx["H2O"]]
        + cell.N_zu_b[idx["H2O"]] + cell.N_b_in[idx["H2O"]]
        + cell.N_rez_d[idx["H2O"]] + cell.N_rez_b[idx["H2O"]],
        1e-12,
    )
    h2o_consume = r2 * total_area + nu_h2o_r11 * (ext11b + ext11d)

    N_out = np.maximum(cell.N_d + cell.N_b, 0.0)
    return {
        "T": float(T),
        "O2_supply": float(o2_supply_bubble + o2_supply_dense),
        "O2_supply_bubble": float(o2_supply_bubble),
        "O2_supply_dense": float(o2_supply_dense),
        "O2_budget_bubble": float(o2_budget_bubble),
        "O2_budget_dense": float(o2_budget_dense),
        "O2_out": float(N_out[idx["O2"]]),
        "limit_factor_o2": float(limit_factor_o2),
        "limit_factor_o2_bubble": float(limit_factor_o2_bubble),
        "limit_factor_o2_dense": float(limit_factor_o2_dense),
        "O2_slack_bubble": float(o2_slack_bubble),
        "O2_slack_dense": float(o2_slack_dense),
        "O2_shortfall_bubble": float(o2_shortfall_bubble),
        "O2_shortfall_dense": float(o2_shortfall_dense),
        "O2_transfer_potential_bd": float(o2_transfer_potential_bd),
        "bubble_o2_unlimited": float(bubble_o2_demand_unlimited),
        "homogeneous_o2_unlimited": float(homogeneous_o2_unlimited),
        "homogeneous_o2_limited_effective": float(homogeneous_o2_limited_effective),
        "homogeneous_share_of_effective_o2": float(homogeneous_share_of_effective_o2),
        "dense_o2_unlimited": float(dense_o2_demand_unlimited),
        "dense_o2_limited": float(dense_o2_demand_limited),
        "total_o2_effective": float(total_o2_demand_effective),
        "r1_char_comb_o2": float(alpha * r1 * total_area),
        "r5b_o2": float(r5b),
        "r6b_o2": float(1.5 * r6b),
        "r12b_o2": float(0.5 * r12b),
        "r5d_o2": float(r5d),
        "r6d_o2": float(1.5 * r6d),
        "r12d_o2": float(0.5 * r12d),
        "r9_o2": float(1.5 * ext9),
        "r10_o2": float(tar_o2_demand_total),
        "vm_co": vm_co,
        "vm_ch4": vm_ch4,
        "vm_h2": vm_h2,
        "vm_tar": vm_tar,
        "ch4_src_vm": vm_ch4,
        "ch4_sink_r6b": float(r6b),
        "ch4_sink_r6d": float(r6d),
        "ch4_sink_r7": float(max(ext7, 0.0)),
        "ch4_net_proxy": float(vm_ch4 - r6b - r6d - max(ext7, 0.0)),
        "vm_oxidizable_o2eq": float(vm_oxidizable_o2eq),
        "endothermic_proxy_r2": float(r2 * total_area),
        "endothermic_proxy_r4": float(r4 * total_area),
        "h2o_supply": float(h2o_supply),
        "h2o_consume_proxy": float(h2o_consume),
        "bubble_share_of_effective_o2": float(bubble_o2_demand_limited / max(total_o2_demand_effective, 1e-12)),
    }


def run_scenario(sc: Scenario):
    with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
        case = load_case_LU()
        cfg = build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)
        cfg.gas_inlet_dense_frac = sc.dense
        cfg.heat_loss_frac = sc.heat_loss
        reactor = Reactor(cfg)
        res = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    rows = []
    for cell in reactor.cells:
        cell.calc_hydrodynamics()
        cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
        cell.calc_exchange()
        rows.append(_cell_o2_trace(cell))

    T_profile = np.array(res["T_profile"], dtype=float)
    i_peak = int(np.argmax(T_profile))
    return cfg, reactor, res, rows, i_peak


def print_scenario(sc: Scenario) -> None:
    cfg, reactor, res, rows, i_peak = run_scenario(sc)
    print("\n" + "=" * 132)
    print(f"Scenario: {sc.name}")
    print(
        f"dense={sc.dense:.2f} heat_loss={sc.heat_loss:.2f} "
        f"impl_scales[R5={sc.r5s:.2f}, R6(thesis R7)={sc.r6s:.2f}, R4(thesis R3)={sc.r4s:.2f}] | "
        f"Texit={res['T_profile'][-1]:.1f}K Tpeak={res['T_profile'][i_peak]:.1f}K@cell{i_peak} Xc={res['carbon_conv']:.3f}"
    )
    print("-" * 132)
    print("labels: thesis 编号优先；R7i6=CH4 oxidation，R6i12=H2 oxidation，R3proxy=impl r4 Boudouard")
    print(
        f"{'cell':>4} {'T[K]':>8} {'O2sup':>9} {'O2out':>9} {'lim':>6} {'O2_homo':>9} {'O2_bub':>9} {'O2_den_l':>10} {'homo%':>7} {'VM_O2eq':>9} {'R1':>8} {'R5b':>8} {'R7i6b':>8} {'R6i12b':>9} {'R5d':>8} {'R7i6d':>8} {'R6i12d':>9} {'R10':>8}"
    )
    for i, r in enumerate(rows):
        peak_mark = "*" if i == i_peak else " "
        print(
            f"{i:>3d}{peak_mark} {r['T']:>8.1f} {r['O2_supply']:>9.3f} {r['O2_out']:>9.3f} {r['limit_factor_o2']:>6.3f} "
            f"{r['homogeneous_o2_limited_effective']:>9.3f} {r['bubble_o2_unlimited']:>9.3f} {r['dense_o2_limited']:>10.3f} {100*r['homogeneous_share_of_effective_o2']:>6.1f}% {r['vm_oxidizable_o2eq']:>9.3f} "
            f"{r['r1_char_comb_o2']:>8.3f} {r['r5b_o2']:>8.3f} {r['r6b_o2']:>8.3f} {r['r12b_o2']:>8.3f} {r['r5d_o2']:>8.3f} {r['r6d_o2']:>8.3f} {r['r12d_o2']:>8.3f} {r['r10_o2']:>8.3f}"
        )

    print("-" * 132)
    lo = max(0, i_peak - 1)
    hi = min(len(rows) - 1, i_peak + 1)
    print("Peak neighborhood analysis:")
    for i in range(lo, hi + 1):
        r = rows[i]
        print(
            f"  cell{i}: T={r['T']:.1f}K, O2_supply={r['O2_supply']:.3f}"
            f" [b={r['O2_supply_bubble']:.3f}, d={r['O2_supply_dense']:.3f}], "
            f"O2_budget=[b={r['O2_budget_bubble']:.3f}, d={r['O2_budget_dense']:.3f}], O2_out={r['O2_out']:.3f}, "
            f"limiter={r['limit_factor_o2']:.3f} [b={r['limit_factor_o2_bubble']:.3f}, d={r['limit_factor_o2_dense']:.3f}], "
            f"homogeneous_O2={r['homogeneous_o2_limited_effective']:.3f} ({100*r['homogeneous_share_of_effective_o2']:.1f}% of effective O2 demand), "
            f"bubble_O2={r['bubble_o2_unlimited']:.3f} ({100*r['bubble_share_of_effective_o2']:.1f}% of effective O2 demand), dense_O2_limited={r['dense_o2_limited']:.3f}, "
            f"slack=[b={r['O2_slack_bubble']:.3f}, d={r['O2_slack_dense']:.3f}] "
            f"shortfall=[b={r['O2_shortfall_bubble']:.3f}, d={r['O2_shortfall_dense']:.3f}] "
            f"transfer_bd={r['O2_transfer_potential_bd']:.3f}, "
            f"VM_O2eq={r['vm_oxidizable_o2eq']:.3f} [CO={r['vm_co']:.3f}, CH4={r['vm_ch4']:.3f}, H2={r['vm_h2']:.3f}, TAR={r['vm_tar']:.3f}], "
            f"R1={r['r1_char_comb_o2']:.3f}, R5b={r['r5b_o2']:.3f}, R7i6b={r['r6b_o2']:.3f}, R6i12b={r['r12b_o2']:.3f}, "
            f"R5d={r['r5d_o2']:.3f}, R7i6d={r['r6d_o2']:.3f}, R6i12d={r['r12d_o2']:.3f}, R10={r['r10_o2']:.3f}, "
            f"CH4[src_vm={r['ch4_src_vm']:.3f}, sink_R7i6b={r['ch4_sink_r6b']:.3f}, sink_R7i6d={r['ch4_sink_r6d']:.3f}, sink_R9i7={r['ch4_sink_r7']:.3f}, net={r['ch4_net_proxy']:.3f}], "
            f"R2proxy={r['endothermic_proxy_r2']:.3f}, R3proxy={r['endothermic_proxy_r4']:.3f}"
        )

    # Automatic reasoning hints
    peak = rows[i_peak]
    reasons: list[str] = []
    if peak["bubble_share_of_effective_o2"] > 0.35:
        reasons.append("峰值 cell 的有效 O2 需求中，bubble 相氧化占比较高；需重点审视 bubble 相氧化与传质分配。")
    if peak["homogeneous_share_of_effective_o2"] > 0.7:
        reasons.append("峰值 cell 主要由快速均相氧化（CO/CH4/H2/TAR/H2S）控制，而不是单纯 char 燃烧；应优先关注挥发分释放后立即氧化的链路。")
    if peak["vm_oxidizable_o2eq"] > 0.2 * max(peak["O2_supply"], 1e-12):
        reasons.append("峰值 cell 存在明显的 VM/挥发分可氧化当量，说明挥发分刚释放后被快速均相氧化，可能直接促成 O2 burn-out 与温峰。")
    if peak["r1_char_comb_o2"] > max(peak["r5d_o2"] + peak["r6d_o2"], 1e-9):
        reasons.append("峰值 cell 以 char combustion (R1) 为主导 O2 汇，说明温峰更像是炭燃烧驱动，而不只是气相氧化。")
    if peak["limit_factor_o2"] < 0.7:
        reasons.append("峰值 cell 已触发明显 O2 限速，说明 dense/char 路径本身就处于强氧受限状态。")
    if peak["bubble_o2_unlimited"] > 0.5 * peak["O2_supply"]:
        reasons.append("bubble 相 O2 需求已达到供氧量的显著比例，需优先审视 bubble O2 分配与 bubble-phase oxidation。")
    if peak["O2_transfer_potential_bd"] > 0.05 * max(peak["O2_supply"], 1e-12):
        reasons.append("峰值 cell 同时存在显著 bubble O2 slack 与 dense O2 shortfall，说明 O2 可能被 limiter 语义锁在 bubble 相内。")
    if peak["endothermic_proxy_r2"] + peak["endothermic_proxy_r4"] < 0.2 * max(peak["dense_o2_limited"] + peak["bubble_o2_unlimited"], 1e-12):
        reasons.append("峰值 cell 的吸热气化代理项（R2/R3）相对偏弱，放热/吸热失衡可能促成温升。")

    if reasons:
        print("Reason hints:")
        for x in reasons:
            print(f"  - {x}")


def main() -> int:
    scenarios = [
        Scenario(name="stable_window_candidate", dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50),
        Scenario(name="peaky_branch_reference", dense=0.35, heat_loss=0.12, r5s=1.00, r6s=1.00, r4s=1.00),
    ]
    for sc in scenarios:
        print_scenario(sc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
