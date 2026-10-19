"""莱茵褐煤的焦炭、焦油产率与热解气平衡组成。"""
from __future__ import annotations

import math
from dataclasses import dataclass

from hamel.properties.data import atomic_weight
from hamel.properties.equilibrium import EquilibriumResult, gibbs_minimization
from hamel.properties.thermo import water_free_to_ash_free
from hamel.pyrolysis.particle import PyrolysisError

# 表 4.3。两式的 ξ 不通用。
_CHAR = (0.3027, 0.1047, 0.9356, -0.08578)
_TAR = (-0.2168, 0.04148, 0.1566, 1.333)
# S. 69 的举例清单里、GLOSSARY 已有的组分。焦油不进入平衡。ADR 0015 第 9 项。
PYROLYSIS_GAS_SPECIES = ("CO", "CO2", "CH4", "C2H4", "NH3", "H2S", "H2O", "H2")

@dataclass(frozen=True)
class LigniteYields:
    """质量都按每 kg 燃料（waf）。苯、萘的量是 mol/kg waf。"""

    x_c_char: float
    x_c_tar: float
    x_char: float
    x_tar: float
    x_c_gas: float
    x_h_gas: float
    x_c: float
    x_h: float
    x_o: float
    x_n: float
    x_s: float
    benzene_mol_per_kg: float
    naphthalene_mol_per_kg: float

def carbon_retained_in_char(pressure_pa: float, temperature_k: float) -> float:
    """kg C in char / kg C in fuel (waf)。

    Source: Hamel Eq.4.21, S. 63
    Source: Hamel Tab.4.3, S. 64
    """
    pressure = float(pressure_pa)
    temperature = float(temperature_k)
    if pressure <= 0.0 or temperature <= 0.0:
        raise PyrolysisError("产率式的压力和温度必须为正")
    factor, exponent, base, power = _CHAR
    return factor * pressure**exponent * base ** (1.0e-6 * pressure) * temperature**power

def carbon_retained_in_tar(pressure_pa: float, temperature_k: float) -> float:
    """kg C in tar / kg C in fuel (waf)。

    Source: Hamel Eq.4.22, S. 64
    Source: Hamel Tab.4.3, S. 64
    """
    pressure = float(pressure_pa)
    temperature = float(temperature_k)
    if pressure <= 0.0 or temperature <= 0.0:
        raise PyrolysisError("产率式的压力和温度必须为正")
    factor, exponent, power, offset = _TAR
    return factor * pressure**exponent * temperature**power + offset

def lignite_yields(
    pressure_pa: float,
    temperature_k: float,
    carbon_mass_fraction_wf: float,
    hydrogen_mass_fraction_wf: float,
    oxygen_mass_fraction_wf: float,
    nitrogen_mass_fraction_wf: float,
    sulfur_mass_fraction_wf: float,
    ash_mass_fraction_wf: float,
) -> LigniteYields:
    """式 (4.23)、(4.24)、(4.26)、(4.27)。wf→waf 用 1/(1−灰)。不读表 7.2 的挥发分。

    Source: Hamel Eq.4.23, S. 64
    Source: Hamel Eq.4.24, S. 64
    Source: Hamel Eq.4.26, S. 69
    Source: Hamel Eq.4.27, S. 69
    """
    ash = float(ash_mass_fraction_wf)
    carbon = water_free_to_ash_free(carbon_mass_fraction_wf, ash)
    hydrogen = water_free_to_ash_free(hydrogen_mass_fraction_wf, ash)
    oxygen = water_free_to_ash_free(oxygen_mass_fraction_wf, ash)
    nitrogen = water_free_to_ash_free(nitrogen_mass_fraction_wf, ash)
    sulfur = water_free_to_ash_free(sulfur_mass_fraction_wf, ash)
    if min(carbon, hydrogen, oxygen, nitrogen, sulfur) < 0.0:
        raise PyrolysisError("元素质量分数不能为负")
    if carbon == 0.0:
        raise PyrolysisError("无水无灰基碳为 0 时，式 (4.24) 的氢碳比没有定义")
    char_carbon = carbon_retained_in_char(pressure_pa, temperature_k)
    tar_carbon = carbon_retained_in_tar(pressure_pa, temperature_k)
    # x_Koks 按燃料（waf）。焦炭在模型中是纯碳。ADR 0015 第 3 项。
    char_mass = carbon * char_carbon
    tar_mass = carbon * tar_carbon * (1.0 + hydrogen / carbon)
    gas_carbon = carbon - char_mass - tar_carbon * carbon
    gas_hydrogen = hydrogen - carbon * tar_carbon * (hydrogen / carbon)
    if char_mass < 0.0 or tar_mass < 0.0 or gas_carbon < 0.0 or gas_hydrogen < 0.0:
        raise PyrolysisError("产率给出负的焦炭、焦油或气相元素，不能分配")
    benzene, naphthalene = _tar_moles(carbon * tar_carbon, tar_mass - carbon * tar_carbon)
    return LigniteYields(
        char_carbon, tar_carbon, char_mass, tar_mass, gas_carbon, gas_hydrogen,
        carbon, hydrogen, oxygen, nitrogen, sulfur, benzene, naphthalene,
    )

def _tar_moles(carbon_kg: float, hydrogen_kg: float) -> tuple[float, float]:
    """苯与萘的比例使摩尔 H/C 等于原煤。解释性选择：ADR 0015 第 4 项。"""
    carbon_mol = carbon_kg / (atomic_weight("C") * 1.0e-3)
    hydrogen_mol = hydrogen_kg / (atomic_weight("H") * 1.0e-3)
    naphthalene = 0.5 * (carbon_mol - hydrogen_mol)
    benzene = (hydrogen_mol - 8.0 * naphthalene) / 6.0
    if benzene < 0.0 or naphthalene < 0.0 or not math.isfinite(benzene + naphthalene):
        raise PyrolysisError("原煤的摩尔 H/C 不能写成苯和萘的非负混合物")
    return benzene, naphthalene

def pyrolysis_gas(
    yields: LigniteYields, temperature_k: float, pressure_pa: float,
) -> EquilibriumResult:
    """元素池为气相 C、H 加燃料（waf）全部 O、N、S。干燥水不并入。

    解释性选择：ADR 0015 第 9 项。
    Source: Hamel Eq.4.26, S. 69
    Source: Hamel Eq.4.28, S. 70
    """
    amounts = {
        "C": yields.x_c_gas / (atomic_weight("C") * 1.0e-3),
        "H": yields.x_h_gas / (atomic_weight("H") * 1.0e-3),
        "O": yields.x_o / (atomic_weight("O") * 1.0e-3),
        "N": yields.x_n / (atomic_weight("N") * 1.0e-3),
        "S": yields.x_s / (atomic_weight("S") * 1.0e-3),
    }
    return gibbs_minimization(
        temperature_k, pressure_pa, amounts, species=list(PYROLYSIS_GAS_SPECIES),
    )
