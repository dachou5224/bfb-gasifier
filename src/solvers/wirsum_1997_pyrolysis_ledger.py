"""阶段 D：热解映射账本。不改生产默认，不做 Startwert。

Hegermann Eq. 4.21–4.22 的 ``x_C,Koks`` / ``x_C,Teer`` 是
**kg C / kg C in Brennstoff waf**。生产 ``allocate_*`` 却把返回值
当成输入 ``nC`` 的分率；若 ``nC`` 来自 ``nC_vm``，基准不一致。

轴向释放用占据概率积分
``ṅ_{py,i,j} = ṅ_{j,∞} ∫ p_i(t) dX_j(t)``。
这是建议的离散表示，**不是** Hamel 闭合式。``p_i`` 必须具名来源，
不能把 ``H/u_mf`` 换个符号再藏进去。

Source: Hamel (1999) §2.1, Eq. 4.21–4.22, Tabelle 4.3;
docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q4
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import T_REF
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.energy_balance_defaults import (
    H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
    H_FORMATION_LIQUID_WATER_J_PER_KG,
    solid_component_hf_j_per_kg,
)
from src.core.species import (
    GAS_SPECIES,
    SOLID_H_MODE_INTEGRAL,
    TAR_SURROGATE_FORMULA,
    delta_h_char_integral,
    enthalpy_molar,
    get_tar_component_mapping,
    h_solid_mass,
)
from src.solvers.vorabrechnung.vorab_policy import VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF
from src.thermal.hegermann_yield import (
    hegermann_applicability,
    xc_koks_hegermann,
    xc_teer_hegermann,
)

CARBON_BASIS_FUEL_WAF: str = "fuel_waf_c"
CARBON_BASIS_VM_RELEASED: str = "vm_released_c"
CARBON_BASIS_VM_RESCALED: str = "vm_rescaled_from_fuel_waf"
OCCUPANCY_PLUG_FLOW_AGE: str = "plug_flow_named_cell_tau"
MAPPING_EXPR: str = "n_py_i_j = n_j_inf * integral p_i(t) dX_j(t)"
REMAINDER_DESTINATION: str = "unconverted_solid_outlet_throughput"
TAR_CONC_DILUTE_SHARED_SLOT: str = "dilute_shared_bed_exit_equal_C"
TAR_CONC_STRICT_EOS: str = "strict_phase_eos_majors_plus_tar"


@dataclass(frozen=True)
class FuelAtomFeed:
    """原料原子流 [mol-atom/s]。水分单独 [kg/s]。"""

    nC_fuel: float
    nH_fuel: float
    nO_fuel: float
    nN_fuel: float
    nS_fuel: float
    nC_fixed: float
    nC_vm: float
    m_moist_kg_s: float
    m_dry_kg_s: float
    h_f_dry_J_kg: float
    ash_dry_wt: float
    VM_daf: float


def dry_mass_from_cho_atoms(feed: FuelAtomFeed) -> float:
    """由 C/H/O 原子流和灰分分率还原干料 [kg/s]。水分另计。"""
    m_cho = (
        float(feed.nC_fuel) * ATOMIC_MASS_KG_PER_MOL["C"]
        + float(feed.nH_fuel) * ATOMIC_MASS_KG_PER_MOL["H"]
        + float(feed.nO_fuel) * ATOMIC_MASS_KG_PER_MOL["O"]
        + float(feed.nN_fuel) * ATOMIC_MASS_KG_PER_MOL["N"]
        + float(feed.nS_fuel) * ATOMIC_MASS_KG_PER_MOL["S"]
    )
    w_ash = float(feed.ash_dry_wt) / 100.0
    return m_cho / max(1.0 - w_ash, 1.0e-12)


def hegermann_fuel_waf_fractions(P_Pa: float, T_K: float) -> tuple[float, float]:
    """Eq. 4.21 / 4.22，相对燃料 waf 碳。"""
    return float(xc_koks_hegermann(P_Pa, T_K)), float(xc_teer_hegermann(P_Pa, T_K))


def fractions_on_named_carbon_pool(
    x_koks_fuel: float,
    x_teer_fuel: float,
    *,
    basis: str,
) -> tuple[float, float, str]:
    """把燃料基准分率转到调用方的 ``nC`` 池。不改生产 ``resolve_*``。"""
    name = str(basis).strip().lower()
    if name == CARBON_BASIS_FUEL_WAF:
        return float(x_koks_fuel), float(x_teer_fuel), name
    if name == CARBON_BASIS_VM_RELEASED:
        return float(x_koks_fuel), float(x_teer_fuel), "mismatch_fuel_frac_on_vm_pool"
    if name == CARBON_BASIS_VM_RESCALED:
        volatile = max(1.0 - float(x_koks_fuel), 1.0e-12)
        return 0.0, float(x_teer_fuel) / volatile, name
    raise ValueError(f"unknown carbon basis {basis!r}")


def feedstock_carbon_split(
    feed: FuelAtomFeed,
    *,
    P_Pa: float,
    T_K: float,
) -> dict[str, Any]:
    """原料 C 对 Hegermann 目标的账。

    具名核焦炭定义是 ``x_C,Koks(T,P) * nC_fuel``（燃料 waf 碳）。
    工业分析 ``nC_fixed`` / ``nC_vm`` 只作元数据，不参与产率分割。
    ``hegermann_koks_minus_fixed`` 是该双重定义的基准检查，不是第二套焦炭。

    Source: Hamel (1999) Eq. 4.21–4.22, Tabelle 4.3
    """
    xk, xt = hegermann_fuel_waf_fractions(P_Pa, T_K)
    nC_koks = xk * feed.nC_fuel
    nC_teer = xt * feed.nC_fuel
    nC_gas = max(feed.nC_fuel - nC_koks - nC_teer, 0.0)
    koks_from_vm = max(nC_koks - feed.nC_fixed, 0.0)
    extra_fixed = max(feed.nC_fixed - nC_koks, 0.0)
    appl = hegermann_applicability(P_Pa, T_K)
    return {
        "x_C_koks_fuel_waf": xk,
        "x_C_teer_fuel_waf": xt,
        "nC_koks_target": nC_koks,
        "nC_teer_target": nC_teer,
        "nC_gas_target": nC_gas,
        "nC_vol_target": nC_teer + nC_gas,
        "nC_fixed": feed.nC_fixed,
        "nC_vm": feed.nC_vm,
        "koks_from_vm": koks_from_vm,
        "extra_fixed_vs_hegermann": extra_fixed,
        "hegermann_koks_minus_fixed": nC_koks - feed.nC_fixed,
        "sum_targets": nC_koks + nC_teer + nC_gas,
        "basis": CARBON_BASIS_FUEL_WAF,
        "coke_definition": CARBON_BASIS_FUEL_WAF,
        "proximate_nC_fixed_is_metadata": True,
        "gibbs_leftover_not_returned_to_koks": True,
        "hegermann": appl,
    }


def occupancy_plug_flow(
    tau_s: npt.NDArray[np.float64],
    *,
    source: str,
) -> dict[str, Any]:
    """各格年龄区间。``source`` 必须写出 ``p_i`` 的来源，禁止匿名。"""
    tau = np.asarray(tau_s, dtype=np.float64)
    assert tau.ndim == 1 and tau.size >= 1
    if float(np.min(tau)) < 0.0:
        raise ValueError("cell tau must be non-negative [s]")
    edges = np.concatenate([[0.0], np.cumsum(tau)])
    return {
        "model": OCCUPANCY_PLUG_FLOW_AGE,
        "p_i_source": str(source),
        "tau_s": tau.copy(),
        "age_edges_s": edges,
        "mapping_expr": MAPPING_EXPR,
        "is_hamel_closed_form": False,
        "hides_H_over_umf": source == VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF,
    }


def map_release_by_occupancy(
    n_inf: float,
    X_of_age: npt.NDArray[np.float64],
    occ: dict[str, Any],
) -> npt.NDArray[np.float64]:
    """``ṅ_i = ṅ_∞ (X(t_{i+1}) − X(t_i))``，X 在格端年龄上取值。

    与 ``∫ p_i dX`` 在 plug-flow（格内 p_i=1）时等价。
    """
    X = np.asarray(X_of_age, dtype=np.float64)
    edges = np.asarray(occ["age_edges_s"], dtype=np.float64)
    assert X.shape == edges.shape
    dX = np.diff(np.clip(X, 0.0, 1.0))
    return float(n_inf) * dX


def _gas_atoms(prod: dict[str, float], *, fuel_type: str) -> dict[str, float]:
    tar_map = get_tar_component_mapping(fuel_type)
    out = {"C": 0.0, "H": 0.0, "O": 0.0, "N": 0.0, "S": 0.0}
    out["C"] += float(prod.get("CO", 0.0) + prod.get("CO2", 0.0) + prod.get("CH4", 0.0))
    out["H"] += float(
        2.0 * prod.get("H2", 0.0)
        + 2.0 * prod.get("H2O", 0.0)
        + 4.0 * prod.get("CH4", 0.0)
        + 3.0 * prod.get("NH3", 0.0)
        + 2.0 * prod.get("H2S", 0.0)
    )
    out["O"] += float(prod.get("CO", 0.0) + 2.0 * prod.get("CO2", 0.0) + prod.get("H2O", 0.0))
    out["N"] += float(prod.get("NH3", 0.0))
    out["S"] += float(prod.get("H2S", 0.0))
    for label in ("TAR1", "TAR2"):
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map[label]]
        n = float(prod.get(label, 0.0))
        out["C"] += c_tar * n
        out["H"] += h_tar * n
    return out


def _gas_enthalpy_W(prod: dict[str, float], T: float) -> float:
    """气相+焦油焓流 [W]。``enthalpy_molar`` 已含生成焓，不再加一遍 hf。"""
    h = 0.0
    for name in GAS_SPECIES:
        n = float(prod.get(name, 0.0))
        if abs(n) > 0.0:
            h += n * enthalpy_molar(name, T)
    return h


def remainder_mass_kg_s(
    nC: float,
    nH: float,
    nO: float,
    nN: float = 0.0,
    nS: float = 0.0,
) -> float:
    """未释放挥发分质量流 [kg/s]，由剩余原子加和。不含焦、灰、水分。"""
    return (
        float(nC) * ATOMIC_MASS_KG_PER_MOL["C"]
        + float(nH) * ATOMIC_MASS_KG_PER_MOL["H"]
        + float(nO) * ATOMIC_MASS_KG_PER_MOL["O"]
        + float(nN) * ATOMIC_MASS_KG_PER_MOL["N"]
        + float(nS) * ATOMIC_MASS_KG_PER_MOL["S"]
    )


def unreleased_vm_enthalpy_W(
    feed: FuelAtomFeed,
    m_unreleased_kg_s: float,
    T: float,
) -> float:
    """未转化 VM 焓流 [W]。生成焓用干料 VM 列；显热用炭 ``c_p`` 积分作代理。

    具名核不为 VM 另建库存 K。该流是稳态未转化固体出口。
    Source: Hamel (1999) Eq. 2.8 固相焓分解
    """
    hf = solid_component_hf_j_per_kg(
        h_f_dry=feed.h_f_dry_J_kg, ash_dry_wt=feed.ash_dry_wt, VM_daf=feed.VM_daf
    )
    h = float(hf[1]) + float(delta_h_char_integral(float(T), T_REF))
    return float(m_unreleased_kg_s) * h


def map_one_state(
    feed: FuelAtomFeed,
    *,
    P_Pa: float,
    T_K: float,
    X_vm: npt.NDArray[np.float64],
    tau_s: npt.NDArray[np.float64],
    occupancy_source: str,
    fuel_type: str = "coal",
) -> dict[str, Any]:
    """一次外层映射 q。源项属 q，出口 N 不在这里改写。

    焦 / 焦油分率相对燃料 waf 碳（Hegermann Eq. 4.21–4.22）。
    DAEM ``X`` 只乘挥发分碳池 ``(1-x_k) nC_fuel``，不乘焦炭。
    未释放挥发分走固体未转化出口，不为 VM 另建 K。

    Source: Hamel (1999) Eq. 4.21–4.22, Tabelle 4.3; §2.1 DAEM
    """
    from src.core.cell_pyrolysis import allocate_pyrolysis_products_hamel

    if str(fuel_type) != "coal":
        raise ValueError(
            "named Hegermann split is Rheinische Braunkohle only; "
            f"got fuel_type={fuel_type!r}"
        )
    occ = occupancy_plug_flow(tau_s, source=occupancy_source)
    split = feedstock_carbon_split(feed, P_Pa=P_Pa, T_K=T_K)
    X = np.asarray(X_vm, dtype=np.float64)
    xk = float(split["x_C_koks_fuel_waf"])
    xt = float(split["x_C_teer_fuel_waf"])
    nC_koks = float(split["nC_koks_target"])
    nC_vol = float(split["nC_vol_target"])
    f_koks_vol, f_tar_vol, _rescale = fractions_on_named_carbon_pool(
        xk, xt, basis=CARBON_BASIS_VM_RESCALED
    )
    nC_rel = map_release_by_occupancy(nC_vol, X, occ)
    nH_rel = map_release_by_occupancy(feed.nH_fuel, X, occ)
    nO_rel = map_release_by_occupancy(feed.nO_fuel, X, occ)
    nN_rel = map_release_by_occupancy(feed.nN_fuel, X, occ)
    nS_rel = map_release_by_occupancy(feed.nS_fuel, X, occ)
    cells: list[dict[str, Any]] = []
    nC_gas_sum = 0.0
    nH_gas_sum = 0.0
    nO_gas_sum = 0.0
    nN_gas_sum = 0.0
    nS_gas_sum = 0.0
    H_prod = 0.0
    T_eval = float(T_K)
    for i, nC in enumerate(nC_rel):
        prod = allocate_pyrolysis_products_hamel(
            T_K=T_K,
            P_Pa=P_Pa,
            nC=float(nC),
            nH=float(nH_rel[i]),
            nO=float(nO_rel[i]),
            nN=float(nN_rel[i]),
            nS=float(nS_rel[i]),
            fuel_type=fuel_type,
            pyrolysis_tar_carbon_frac=f_tar_vol,
            pyrolysis_koks_carbon_frac=f_koks_vol,
        )
        atoms = _gas_atoms(prod, fuel_type=fuel_type)
        nC_gas_sum += atoms["C"]
        nH_gas_sum += atoms["H"]
        nO_gas_sum += atoms["O"]
        nN_gas_sum += atoms["N"]
        nS_gas_sum += atoms["S"]
        cells.append(
            {
                "n_py_mol_s": {sp: float(prod.get(sp, 0.0)) for sp in GAS_SPECIES},
                "atoms": atoms,
                "nC_released_vm": float(nC),
            }
        )
        H_prod += _gas_enthalpy_W(prod, T_eval)
    X_end = float(np.clip(X[-1], 0.0, 1.0))
    nC_unreleased = nC_vol * (1.0 - X_end)
    nH_unreleased = float(feed.nH_fuel) * (1.0 - X_end)
    nO_unreleased = float(feed.nO_fuel) * (1.0 - X_end)
    nN_unreleased = float(feed.nN_fuel) * (1.0 - X_end)
    nS_unreleased = float(feed.nS_fuel) * (1.0 - X_end)
    m_unreleased = remainder_mass_kg_s(
        nC_unreleased, nH_unreleased, nO_unreleased, nN_unreleased, nS_unreleased
    )
    nC_accounted = nC_koks + nC_unreleased + nC_gas_sum
    nH_accounted = nH_unreleased + nH_gas_sum
    nO_accounted = nO_unreleased + nO_gas_sum
    ash_frac = float(feed.ash_dry_wt)
    if ash_frac > 1.0:
        ash_frac = ash_frac / 100.0
    m_ash = float(feed.m_dry_kg_s) * ash_frac
    m_koks = nC_koks * ATOMIC_MASS_KG_PER_MOL["C"]
    H_in = (
        feed.m_dry_kg_s * feed.h_f_dry_J_kg
        + feed.m_moist_kg_s * H_FORMATION_LIQUID_WATER_J_PER_KG
    )
    H_koks = m_koks * h_solid_mass(
        "char",
        T_eval,
        h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    H_unreleased = unreleased_vm_enthalpy_W(feed, m_unreleased, T_eval)
    H_ash = m_ash * h_solid_mass("ash", T_eval, h_f=0.0, mode=SOLID_H_MODE_INTEGRAL)
    H_moist = feed.m_moist_kg_s * H_FORMATION_LIQUID_WATER_J_PER_KG
    H_solid_remain = H_koks + H_unreleased + H_ash + H_moist
    return {
        "occupancy": occ,
        "split": split,
        "carbon_basis": CARBON_BASIS_FUEL_WAF,
        "vm_pool_basis": CARBON_BASIS_FUEL_WAF,
        "remainder_destination": REMAINDER_DESTINATION,
        "cells": cells,
        "nC_koks": nC_koks,
        "nC_released_vm_sum": float(np.sum(nC_rel)),
        "nC_gas_plus_tar": nC_gas_sum,
        "nC_unreleased_vm": nC_unreleased,
        "nH_unreleased": nH_unreleased,
        "nO_unreleased": nO_unreleased,
        "nN_unreleased": nN_unreleased,
        "nS_unreleased": nS_unreleased,
        "m_unreleased_kg_s": m_unreleased,
        "m_koks_kg_s": m_koks,
        "m_ash_kg_s": m_ash,
        "nC_accounted": nC_accounted,
        "nH_accounted": nH_accounted,
        "nO_accounted": nO_accounted,
        "nC_in": feed.nC_fuel,
        "atom_C_gap": abs(nC_accounted - feed.nC_fuel),
        "atom_H_gap": abs(nH_accounted - feed.nH_fuel),
        "atom_O_gap": abs(nO_accounted - feed.nO_fuel),
        "hegermann_koks_minus_fixed": float(split["hegermann_koks_minus_fixed"]),
        "hegermann_coke_basis_gap": abs(nC_koks - xk * feed.nC_fuel),
        "hegermann": split["hegermann"],
        "X_end": X_end,
        "H_in_W": H_in,
        "H_prod_gas_tar_W": H_prod,
        "H_koks_W": H_koks,
        "H_unreleased_W": H_unreleased,
        "H_ash_W": H_ash,
        "H_moist_W": H_moist,
        "H_solid_remain_W": H_solid_remain,
        "H_map_gap_W": H_in - H_prod - H_solid_remain,
        "enthalpy_T_K": T_eval,
        "gas_hf_counted_once": True,
        "koks_from_vm_pool": 0.0,
        "not_startwert": True,
        "not_vm_K": True,
        "proximate_nC_fixed_is_metadata": True,
    }


def two_outer_maps(
    feed: FuelAtomFeed,
    *,
    state_a: dict[str, Any],
    state_b: dict[str, Any],
) -> dict[str, Any]:
    """两次外层 ``q=P(x)``。各次原子闭合；同一公式，不重写未知 N。"""
    qa = map_one_state(feed, **state_a)
    qb = map_one_state(feed, **state_b)
    return {
        "map_a": qa,
        "map_b": qb,
        "same_formula": True,
        "rewrites_unknown_N": False,
        "default_policy_unchanged": True,
    }


def example_lu_like_feed() -> FuelAtomFeed:
    """示意原料，不是实验夹带。单位 SI。"""
    nC_fuel = 10.0
    nC_vm = 5.342
    return FuelAtomFeed(
        nC_fuel=nC_fuel,
        nH_fuel=12.0,
        nO_fuel=4.0,
        nN_fuel=0.1,
        nS_fuel=0.05,
        nC_fixed=nC_fuel - nC_vm,
        nC_vm=nC_vm,
        m_moist_kg_s=0.02,
        m_dry_kg_s=0.20,
        h_f_dry_J_kg=-3.96e6,
        ash_dry_wt=11.41,
        VM_daf=53.42,
    )
