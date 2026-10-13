"""Wesseling 1 的两区条件性边界。几何和已知物流换成这套中试炉，不读认证端点。

蒸汽取 0。全部表载空气进气泡入口。热值沿用具名核示意值，壁面热损为 0。
轨迹窗口仍是 5 s，不因自由空间变高而加长。

Source: data/validation_cases.json CASE_HTW_WESSELING_1;
docs/audit_reports/plan_ab_diagnosis_2026-09-28.md 工作包 4
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import (
    default_char_feed,
    edges_from_x,
    moist_inlet_mol_s,
)
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    default_x0_named,
    freeze_stage_q,
    outer_P_named,
    rebuild_feed_derived_inlets,
)
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_named_f_p2 import capture_named_reference
from src.solvers.wirsum_1997_pyrolysis_ledger import FuelAtomFeed, dry_mass_from_cho_atoms
from src.solvers.wirsum_1997_stiffness_full_gas import MAJOR
from src.solvers.wesseling_case_contract import (
    external_element_balance_kmol_h,
    wesseling1_conditional_boundary,
    wesseling1_external_inlet_kmol_h,
)

_KMOL_H_PER_MOL_S: float = 3.6  # [kmol/h] / [mol/s]


def build_wesseling1_two_zone_q() -> NamedStageQ:
    """按条件边界建 stage 11 的 q。反应开关与具名核 E8/P3 相同。"""
    bound = wesseling1_conditional_boundary()
    dry = bound["dry_element_kmol_h"]
    nC = float(dry["C"]) / _KMOL_H_PER_MOL_S
    feed = FuelAtomFeed(
        nC_fuel=nC,
        nH_fuel=float(dry["H"]) / _KMOL_H_PER_MOL_S,
        nO_fuel=float(dry["O"]) / _KMOL_H_PER_MOL_S,
        nN_fuel=float(dry["N"]) / _KMOL_H_PER_MOL_S,
        nS_fuel=float(dry["S"]) / _KMOL_H_PER_MOL_S,
        nC_fixed=nC * float(bound["nC_fixed_share_of_fuel_carbon"]),
        nC_vm=nC * float(bound["nC_vm_share_of_fuel_carbon"]),
        m_moist_kg_s=float(bound["moisture_kg_h"]) / 3600.0,
        m_dry_kg_s=0.0,
        h_f_dry_J_kg=float(default_char_feed().h_f_dry_J_kg),
        ash_dry_wt=float(bound["ash_dry_wt_pct"]),
        VM_daf=float(bound["VM_daf_pct"]),
    )
    feed = replace(feed, m_dry_kg_s=dry_mass_from_cho_atoms(feed))
    enabled = freeze_e8_q(Q_W=0.0).enabled_rxns
    nq = freeze_stage_q(
        11,
        Q_W=float(bound["Q_W"]),
        enabled_rxns=enabled,
        T_in=float(bound["T_in_K"]),
        feed=feed,
    )
    sp = nq.energy.species
    n_in_b = np.zeros(7, dtype=np.float64)
    n_in_b[MAJOR.index("O2")] = float(bound["O2_mol_s"])
    n_in_b[MAJOR.index("N2")] = float(bound["N2_mol_s"])
    sp2 = replace(
        sp,
        geo_fb=replace(
            sp.geo_fb,
            area=float(bound["area_m2"]),
            H_freeboard=float(bound["H_freeboard_m"]),
        ),
        geo_ex=replace(
            sp.geo_ex,
            P=float(bound["pressure_Pa"]),
            V_d=float(bound["V_d_m3"]),
            V_b=float(bound["V_b_m3"]),
        ),
        H_bed=float(bound["H_bed_m"]),
        N_in_d=np.zeros(7, dtype=np.float64),
        N_in_b=n_in_b,
    )
    nq = replace(nq, energy=replace(nq.energy, species=sp2))
    ref = replace(
        capture_named_reference(nq),
        trajectory_continue_s=15.0,
        trajectory_extra_physical_s=20.0,
        trajectory_max_work_steps=20000,
        trajectory_wall_s=5.0,
    )
    return replace(nq, reference=ref)


def _steady_inventories_kg(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """用当前冻结的抛射、回落和旋风系数解炭、灰、自由空间库存。

    灰的质量流量不直接写入库存。反应消耗和粒级迁移不在这张初值矩阵里。
    不读取认证端点。
    """
    from src.solvers.wirsum_1997_named_f import ash_feed_kg_s
    from src.solvers.wirsum_1997_named_f_e8 import cyclone_equilibrium_mass_kg

    q = outer_P_named(x, nq)
    sp = q.energy.species
    k_eject = np.maximum(np.asarray(sp.K_eject, dtype=np.float64), 1.0e-6)
    w_leave = np.asarray(sp.w_leave, dtype=np.float64)
    w_return = np.asarray(sp.w_return, dtype=np.float64)
    m_feed = np.maximum(np.asarray(sp.m_feed, dtype=np.float64), 0.0)
    k_up = float(sp.k.K_auf)
    k_down = float(sp.k.K_ab)
    denom = k_up * k_up + k_up * k_down + k_down * k_down
    recycle_factor = k_down * (k_up + k_down) / denom if denom > 0.0 else 0.0
    removal = np.maximum((1.0 - w_return) - recycle_factor * w_leave, 1.0e-3)
    bed = m_feed / k_eject / removal
    leave = w_leave * k_eject * bed
    cyclone = cyclone_equilibrium_mass_kg(float(np.sum(leave)), k_up, k_down)
    fly = k_up / max(k_up + k_down, 1.0e-30)
    ash = float(ash_feed_kg_s(nq)) / max(fly * float(k_eject[-1]), 1.0e-6)
    tau = np.maximum(np.asarray(q.hydro.tau_s, dtype=np.float64), 0.0)
    out = x.copy()
    out[14:16] = np.maximum(bed, 1.0e-12)
    out[37:39] = np.maximum(cyclone, 0.0)
    out[39] = max(ash, 1.0e-12)
    out[46:48] = np.maximum(leave * tau, 0.0)
    return out


def wesseling1_physical_x0(nq: NamedStageQ) -> npt.NDArray[np.float64]:
    """密相放热解气，气泡和自由空间放空气。库存由冻结输运闭合，不读认证端点。

    温度在 900、1100、1300 K 里取分区能量绝对值之和最小的一档。这是条件焓下的
    候选，不是测得的炉温。
    """
    from src.solvers.wirsum_1997_named_f import coupled_residual_named

    x = default_x0_named(nq)
    sp = nq.energy.species
    dense = np.asarray(sp.n_py_d, dtype=np.float64) + np.asarray(sp.n_moist, dtype=np.float64)
    if float(dense[5]) <= 1.0e-12:
        dense = dense.copy()
        dense[5] = 1.0e-12
    air = np.asarray(sp.N_in_b, dtype=np.float64)
    x[0:7] = dense
    x[7:14] = air
    x[16:23] = air
    x[23:30] = air
    x[30:37] = air
    best = x
    best_energy = float("inf")
    for temperature in (900.0, 1100.0, 1300.0):
        trial = x.copy()
        trial[-2] = float(temperature)
        trial[-1] = float(temperature)
        trial = _steady_inventories_kg(trial, nq)
        residual = coupled_residual_named(trial, outer_P_named(trial, nq))
        energy = float(np.sum(np.abs(residual[-2:])))
        if energy < best_energy:
            best = trial
            best_energy = energy
    return best


def _log_blend(named: float, target: float, s: float) -> float:
    """正量按 ``named * (target/named)^s``。s=0 是具名核，s=1 是目标边界。"""
    a = float(named)
    b = float(target)
    if a > 0.0 and b > 0.0:
        return float(np.exp((1.0 - float(s)) * np.log(a) + float(s) * np.log(b)))
    return (1.0 - float(s)) * a + float(s) * b


def q_on_wesseling1_path(s: float) -> NamedStageQ:
    """具名核到 Wesseling 1 条件边界的输入。0<s<1 不叫 Wesseling。

    正的长度、面积和摩尔流走对数插值。炭粒级仍是具名核的两档，不把缺失的
    10 个燃料粒级写进去。干料生成焓仍是示意值，不用热值填这个字段。
    """
    step = float(s)
    if step < -1.0e-15 or step > 1.0 + 1.0e-15:
        raise ValueError("path coordinate s is outside [0, 1]")
    named = replace(freeze_e8_q(Q_W=0.0), stage=11)
    if named.reference is None:
        named = replace(named, reference=capture_named_reference(named))
    if abs(step) <= 1.0e-15:
        return named
    target = build_wesseling1_two_zone_q()
    if abs(step - 1.0) <= 1.0e-15:
        return target
    nf = named.energy.feed
    tf = target.energy.feed
    nC = _log_blend(nf.nC_fuel, tf.nC_fuel, step)
    share_named = float(nf.nC_fixed) / max(float(nf.nC_fuel), 1.0e-30)
    share_target = float(tf.nC_fixed) / max(float(tf.nC_fuel), 1.0e-30)
    share = (1.0 - step) * share_named + step * share_target
    feed = replace(
        nf,
        nC_fuel=nC,
        nH_fuel=_log_blend(nf.nH_fuel, tf.nH_fuel, step),
        nO_fuel=_log_blend(nf.nO_fuel, tf.nO_fuel, step),
        nN_fuel=_log_blend(nf.nN_fuel, tf.nN_fuel, step),
        nS_fuel=_log_blend(nf.nS_fuel, tf.nS_fuel, step),
        nC_fixed=nC * share,
        nC_vm=nC * (1.0 - share),
        m_moist_kg_s=_log_blend(nf.m_moist_kg_s, tf.m_moist_kg_s, step),
        h_f_dry_J_kg=float(nf.h_f_dry_J_kg),
        ash_dry_wt=(1.0 - step) * float(nf.ash_dry_wt) + step * float(tf.ash_dry_wt),
        VM_daf=(1.0 - step) * float(nf.VM_daf) + step * float(tf.VM_daf),
    )
    feed = replace(feed, m_dry_kg_s=dry_mass_from_cho_atoms(feed))
    n_in_b = _log_blend_vec(named.energy.species.N_in_b, target.energy.species.N_in_b, step)
    n_in_d = _log_blend_vec(named.energy.species.N_in_d, target.energy.species.N_in_d, step)
    sp = named.energy.species
    tsp = target.energy.species
    sp2 = replace(
        sp,
        geo_fb=replace(
            sp.geo_fb,
            area=_log_blend(sp.geo_fb.area, tsp.geo_fb.area, step),
            H_freeboard=_log_blend(sp.geo_fb.H_freeboard, tsp.geo_fb.H_freeboard, step),
        ),
        geo_ex=replace(
            sp.geo_ex,
            P=_log_blend(sp.geo_ex.P, tsp.geo_ex.P, step),
            V_d=_log_blend(sp.geo_ex.V_d, tsp.geo_ex.V_d, step),
            V_b=_log_blend(sp.geo_ex.V_b, tsp.geo_ex.V_b, step),
        ),
        H_bed=_log_blend(sp.H_bed, tsp.H_bed, step),
        N_in_b=n_in_b,
        N_in_d=n_in_d,
        n_moist=moist_inlet_mol_s(float(feed.m_moist_kg_s)),
    )
    nq = replace(
        named,
        energy=replace(named.energy, species=sp2, feed=feed, T_in=(1.0 - step) * float(named.energy.T_in) + step * float(target.energy.T_in)),
        Q_W=0.0,
    )
    nq = rebuild_feed_derived_inlets(nq)
    continue_s = (1.0 - step) * float(named.reference.trajectory_continue_s) + step * float(
        target.reference.trajectory_continue_s
    )
    ref = replace(
        capture_named_reference(nq),
        trajectory_continue_s=float(continue_s),
        trajectory_extra_physical_s=20.0,
        trajectory_max_work_steps=20000,
        trajectory_wall_s=5.0,
    )
    return replace(nq, reference=ref)


def _log_blend_vec(named: npt.NDArray[np.float64], target: npt.NDArray[np.float64], s: float) -> npt.NDArray[np.float64]:
    a = np.asarray(named, dtype=np.float64)
    b = np.asarray(target, dtype=np.float64)
    out = np.zeros_like(a)
    for i in range(int(a.size)):
        out[i] = _log_blend(float(a[i]), float(b[i]), s)
    return out


def wesseling1_input_path() -> dict[str, Any]:
    """从具名核到 Wesseling 1 的输入比例。0 到 1 之间的点不叫 Wesseling。"""
    named = freeze_e8_q(Q_W=0.0)
    target = build_wesseling1_two_zone_q()
    n_in = named.energy.species.N_in_b
    t_in = target.energy.species.N_in_b
    rows = {
        "area_m2": (float(named.reference.area), float(target.reference.area)),
        "H_bed_m": (float(named.reference.H_bed), float(target.reference.H_bed)),
        "H_freeboard_m": (float(named.reference.H_freeboard), float(target.reference.H_freeboard)),
        "O2_mol_s": (float(n_in[5]), float(t_in[5])),
        "fuel_C_mol_s": (float(named.energy.feed.nC_fuel), float(target.energy.feed.nC_fuel)),
    }
    return {
        "s_named_kernel": 0.0,
        "s_wesseling1": 1.0,
        "intermediate_name": "not CASE_HTW_WESSELING_1",
        "quantities": {
            key: {"named": a, "target": b, "ratio": b / a if abs(a) > 0.0 else None}
            for key, (a, b) in rows.items()
        },
        "blend": "log on positive lengths, areas and molar flows; linear when an endpoint is zero",
        "intermediate_identity": "artificial homotopy; area, height and phase volumes are blended separately",
        "historical_s_0_22": "numerical root before feed water entered n_moist; not a physical certification frontier",
        "first_step_s": 0.02,
        "target_reached": False,
        "cold_start_of_target_run": False,
        "fuel_heating_value": None,
        "fuel_size_applied_to_bins": False,
    }


def wesseling1_input_reconciliation(nq: NamedStageQ | None = None) -> dict[str, Any]:
    """适配器实际使用量和合同逐项对照。粒径不一致就标出来，不改成相符。"""
    bound = wesseling1_conditional_boundary()
    q = build_wesseling1_two_zone_q() if nq is None else nq
    feed = q.energy.feed
    sp = q.energy.species
    dry = bound["dry_element_kmol_h"]
    used = {
        "C_mol_s": float(feed.nC_fuel),
        "H_mol_s": float(feed.nH_fuel),
        "O_mol_s": float(feed.nO_fuel),
        "moisture_kg_s": float(feed.m_moist_kg_s),
        "O2_mol_s": float(sp.N_in_b[MAJOR.index("O2")]),
        "N2_mol_s": float(sp.N_in_b[MAJOR.index("N2")]),
        "area_m2": float(q.reference.area),
        "H_bed_m": float(q.reference.H_bed),
        "V_d_m3": float(q.reference.V_d),
        "V_b_m3": float(q.reference.V_b),
        "Q_W": float(q.Q_W),
        "T_in_K": float(q.energy.T_in),
        "d_fine_m": float(sp.d_p[0]),
        "d_coarse_m": float(sp.d_p[1]),
    }
    contract = {
        "C_mol_s": float(dry["C"]) / _KMOL_H_PER_MOL_S,
        "H_mol_s": float(dry["H"]) / _KMOL_H_PER_MOL_S,
        "O_mol_s": float(dry["O"]) / _KMOL_H_PER_MOL_S,
        "moisture_kg_s": float(bound["moisture_kg_h"]) / 3600.0,
        "O2_mol_s": float(bound["O2_mol_s"]),
        "N2_mol_s": float(bound["N2_mol_s"]),
        "area_m2": float(bound["area_m2"]),
        "H_bed_m": float(bound["H_bed_m"]),
        "V_d_m3": float(bound["V_d_m3"]),
        "V_b_m3": float(bound["V_b_m3"]),
        "Q_W": float(bound["Q_W"]),
        "T_in_K": float(bound["T_in_K"]),
    }
    rows = []
    for key, expected in contract.items():
        rows.append({
            "name": key,
            "contract": expected,
            "adapter": used[key],
            "match": abs(used[key] - expected) <= 1.0e-6 * max(1.0, abs(expected)),
        })
    return {
        "rows": rows,
        "all_applied_inputs_match": all(row["match"] for row in rows),
        "secondary_air_Nm3_h": bound["secondary_air_Nm3_h"],
        "secondary_air_Nm3_h_historical": bound["secondary_air_Nm3_h_historical"],
        "secondary_air_xi": bound["secondary_air_xi"],
        "secondary_air_xi_digitized": bound["secondary_air_xi_digitized"],
        "secondary_air_placed": False,
        "steam_status": bound["steam_status"],
        "bed_material_status": bound["bed_material_status"],
        "recycle_mass_flow_status": bound["recycle_mass_flow_status"],
        "fuel_size_mm": bound["fuel_size_mm"],
        "adapter_char_diameter_m": (used["d_fine_m"], used["d_coarse_m"]),
        "fuel_size_applied_to_bins": False,
        "formation_enthalpy_J_kg": float(feed.h_f_dry_J_kg),
        "Q_W": float(q.Q_W),
    }


def freeboard_wet_mole_fraction(x: npt.NDArray[np.float64]) -> dict[str, float]:
    """自由空间湿基摩尔分数。两区出口取这段气体，不取旋风后。"""
    n_fb = np.asarray(x[16:23], dtype=np.float64)
    total = float(np.sum(n_fb))
    if total <= 0.0:
        return {name: 0.0 for name in MAJOR}
    return {name: float(n_fb[i] / total) for i, name in enumerate(MAJOR)}


def gas_phase_element_gap_kmol_h(x: npt.NDArray[np.float64]) -> dict[str, float]:
    """外部进料减去自由空间气体。固体和焦油不在这项里，所以它不是闭合账。"""
    n_fb = np.asarray(x[16:23], dtype=np.float64)
    carbon = n_fb[0] + n_fb[1] + n_fb[4]
    hydrogen = 2.0 * n_fb[2] + 2.0 * n_fb[3] + 4.0 * n_fb[4]
    oxygen = n_fb[0] + 2.0 * n_fb[1] + n_fb[3] + 2.0 * n_fb[5]
    nitrogen = 2.0 * n_fb[6]
    outlet = {
        "C": carbon * _KMOL_H_PER_MOL_S,
        "H": hydrogen * _KMOL_H_PER_MOL_S,
        "O": oxygen * _KMOL_H_PER_MOL_S,
        "N": nitrogen * _KMOL_H_PER_MOL_S,
    }
    inlet = wesseling1_external_inlet_kmol_h()
    inlet_use = {key: float(inlet[key]) for key in ("C", "H", "O", "N")}
    return external_element_balance_kmol_h(inlet_use, outlet)
