"""热化学恒等式：纯碳 char、液态水、干料 HHV 分解。

不把「当前 x0 能量残差是否下降」当作物性正确性判据。

Source: Hamel (1999) 第 5 章标准反应焓；docs/audit_reports/secondary_diagnosis_2026-09-14.md §2
"""

from __future__ import annotations

from typing import Mapping

from src.core.constants import T_REF
from src.core.energy_balance_defaults import (
    H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
    H_FORMATION_LIQUID_WATER_J_PER_KG,
    solid_component_hf_j_per_kg,
)
from src.core.species import (
    HF_H2O_LIQUID_J_PER_MOL,
    MW_C_KG_PER_MOL,
    enthalpy_molar,
    formation_enthalpy_dry_fuel,
)

# Hamel 第 5 章 / 纯碳手册值 [J/mol C]，容限按物性拟合精度，不追求跨来源小数对齐。
HAMEL_C_PLUS_O2_TO_CO2_J_PER_MOL: float = -393.5e3
HAMEL_C_PLUS_H2O_G_TO_CO_H2_J_PER_MOL: float = 131.4e3
HAMEL_C_PLUS_CO2_TO_2CO_J_PER_MOL: float = 172.5e3


def char_molar_formation_enthalpy_j_per_mol() -> float:
    """纯碳 char 摩尔生成焓 [J/mol]。石墨基准为 0。"""
    return float(H_FORMATION_GRAPHITE_CHAR_J_PER_KG) * float(MW_C_KG_PER_MOL)


def standard_char_reaction_enthalpy_j_per_mol(name: str, *, T: float = T_REF) -> float:
    """按当前气体焓 + 纯碳 char 计算标准反应焓 [J/mol C]。

    ``C+O2→CO2``、``C+H2O(g)→CO+H2``、``C+CO2→2CO``。
    """
    hf_c = char_molar_formation_enthalpy_j_per_mol()
    t = float(T)
    if name == "C+O2→CO2":
        return float(enthalpy_molar("CO2", t) - hf_c - enthalpy_molar("O2", t))
    if name == "C+H2O(g)→CO+H2":
        return float(
            enthalpy_molar("CO", t)
            + enthalpy_molar("H2", t)
            - hf_c
            - enthalpy_molar("H2O", t)
        )
    if name == "C+CO2→2CO":
        return float(2.0 * enthalpy_molar("CO", t) - hf_c - enthalpy_molar("CO2", t))
    raise ValueError(f"Unknown char reaction {name!r}")


def dry_fuel_mass_fractions(*, ash_dry_wt: float, VM_daf: float) -> tuple[float, float, float]:
    """1 kg 干料的 char / VM / ash 质量分率 [-]。"""
    w_ash = float(ash_dry_wt) / 100.0
    w_vm = (float(VM_daf) / 100.0) * (1.0 - w_ash)
    w_char = max(1.0 - w_ash - w_vm, 0.0)
    return float(w_char), float(w_vm), float(w_ash)


def dry_fuel_component_enthalpy_j_per_kg(
    *,
    h_f_dry: float,
    ash_dry_wt: float,
    VM_daf: float,
) -> float:
    """char+VM+ash 分列后再相加的干料生成焓 [J/kg]。应等于 ``h_f_dry``。"""
    w_char, w_vm, w_ash = dry_fuel_mass_fractions(ash_dry_wt=ash_dry_wt, VM_daf=VM_daf)
    hf_char, hf_vm, _hf_moist, hf_ash = solid_component_hf_j_per_kg(
        h_f_dry=float(h_f_dry),
        ash_dry_wt=float(ash_dry_wt),
        VM_daf=float(VM_daf),
    )
    return float(w_char * hf_char + w_vm * hf_vm + w_ash * hf_ash)


def dry_fuel_hhv_identity_gap_j_per_kg(
    *,
    C: float,
    H: float,
    O: float,
    S: float,
    HHV_dry_MJ_kg: float,
    ash_dry_wt: float,
    VM_daf: float,
) -> dict[str, float]:
    """HHV 反算干料焓 vs 分列相加。gap 应为 0。"""
    h_f_dry = formation_enthalpy_dry_fuel(C, H, O, S, HHV_dry_MJ_kg)
    h_components = dry_fuel_component_enthalpy_j_per_kg(
        h_f_dry=h_f_dry,
        ash_dry_wt=ash_dry_wt,
        VM_daf=VM_daf,
    )
    return {
        "h_f_dry_J_kg": float(h_f_dry),
        "h_components_J_kg": float(h_components),
        "gap_J_kg": float(h_components - h_f_dry),
        "h_f_liquid_water_J_kg": float(H_FORMATION_LIQUID_WATER_J_PER_KG),
        "h_f_liquid_water_J_mol": float(HF_H2O_LIQUID_J_PER_MOL),
    }


def liquid_to_vapor_enthalpy_j_per_mol(*, T: float = T_REF) -> float:
    """液态水 → 水蒸气焓差 [J/mol]。应为正（汽化吸热）。"""
    return float(enthalpy_molar("H2O", float(T)) - HF_H2O_LIQUID_J_PER_MOL)


def hamel_char_reaction_targets_j_per_mol() -> Mapping[str, float]:
    return {
        "C+O2→CO2": HAMEL_C_PLUS_O2_TO_CO2_J_PER_MOL,
        "C+H2O(g)→CO+H2": HAMEL_C_PLUS_H2O_G_TO_CO_H2_J_PER_MOL,
        "C+CO2→2CO": HAMEL_C_PLUS_CO2_TO_2CO_J_PER_MOL,
    }
