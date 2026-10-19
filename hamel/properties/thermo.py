"""NASA 多项式、理想气体与燃料生成焓。"""
from __future__ import annotations

import math

from hamel.properties.data import (
    PropertyError,
    PropertyRangeError,
    ash_record,
    atomic_weight,
    gas_constant_j_mol_k,
    is_gas,
    molar_mass_kg_mol,
    nasa_species,
    reference_pressure_pa,
    standard_pressure_pa,
)

# ADR 0018 第 1 项。气体从 1 bar 换到 P_0；凝相不换。
_P0 = standard_pressure_pa()
_P_BAR = reference_pressure_pa()
def _interval(name: str, temperature_k: float) -> dict:
    """ADR 0021 第 8 项：H2O(L) 越出 273.15–600 K 报错，不外推。其他组分同样不外推。"""
    record = nasa_species(name)
    temperature = float(temperature_k)
    if not math.isfinite(temperature):
        raise PropertyError(f"{name} 的温度不是有限数")
    intervals = record["intervals"]
    for index, interval in enumerate(intervals):
        low = float(interval["t_low_K"])
        high = float(interval["t_high_K"])
        last = index == len(intervals) - 1
        if low <= temperature < high or (last and temperature == high):
            return interval
        if temperature == high and not last:
            return intervals[index + 1]
    low = float(intervals[0]["t_low_K"])
    high = float(intervals[-1]["t_high_K"])
    raise PropertyRangeError(f"{name} 的温度 {temperature} K 不在 {low}–{high} K 内")
def _dimensionless(name: str, temperature_k: float) -> tuple[float, float, float]:
    """NASA/TP-2002-211556 式 (1)–(3)，p. 2。返回 c_p/R、H/(RT)、S/R。"""
    temperature = float(temperature_k)
    a1, a2, a3, a4, a5, a6, a7, b1, b2 = _interval(name, temperature)["coefficients"]
    log_t = math.log(temperature)
    cp_over_r = (
        a1 / temperature**2 + a2 / temperature + a3 + a4 * temperature
        + a5 * temperature**2 + a6 * temperature**3 + a7 * temperature**4
    )
    enthalpy_over_rt = (
        -a1 / temperature**2 + a2 * log_t / temperature + a3 + a4 * temperature / 2
        + a5 * temperature**2 / 3 + a6 * temperature**3 / 4 + a7 * temperature**4 / 5
        + b1 / temperature
    )
    entropy_over_r = (
        -a1 / (2 * temperature**2) - a2 / temperature + a3 * log_t + a4 * temperature
        + a5 * temperature**2 / 2 + a6 * temperature**3 / 3 + a7 * temperature**4 / 4 + b2
    )
    return cp_over_r, enthalpy_over_rt, entropy_over_r


def heat_capacity_j_mol_k(name: str, temperature_k: float) -> float:
    """摩尔定压热容，J/(mol·K)。 Source: NASA/TP-2002-211556 Eq.(1), p. 2"""
    return _dimensionless(name, temperature_k)[0] * gas_constant_j_mol_k()


def enthalpy_j_mol(name: str, temperature_k: float) -> float:
    """绝对焓，含生成焓，J/mol。不随标准压力变化。 Source: NASA/TP-2002-211556 Eq.(2), p. 2"""
    temperature = float(temperature_k)
    return _dimensionless(name, temperature)[1] * gas_constant_j_mol_k() * temperature


def formation_enthalpy_j_mol(name: str) -> float:
    """298.15 K 的标准生成焓，J/mol。SO2 只作燃料生成焓的参照。"""
    return float(nasa_species(name)["formation_enthalpy_J_mol"])


def entropy_j_mol_k(name: str, temperature_k: float, standard_pressure_pa: float | None = None) -> float:
    """摩尔熵，J/(mol·K)。气体默认换到 P_0；传入 1e5 则停在报告的 1 bar。 Source: NASA/TP-2002-211556 Eq.(3), p. 2"""
    entropy = _dimensionless(name, temperature_k)[2] * gas_constant_j_mol_k()
    if not is_gas(name):
        return entropy
    pressure = _P0 if standard_pressure_pa is None else float(standard_pressure_pa)
    if pressure <= 0:
        raise PropertyError("标准压力必须为正")
    return entropy - gas_constant_j_mol_k() * math.log(pressure / _P_BAR)


def chemical_potential_j_mol(
    name: str, temperature_k: float, standard_pressure_pa: float | None = None,
) -> float:
    """μ^0 = h − T·s。气体的 s 默认在 P_0。ADR 0018 第 1 项。"""
    temperature = float(temperature_k)
    enthalpy = enthalpy_j_mol(name, temperature)
    entropy = entropy_j_mol_k(name, temperature, standard_pressure_pa)
    return enthalpy - temperature * entropy


def heat_capacity_volume_j_mol_k(name: str, temperature_k: float) -> float:
    """理想气体 C_v = C_p − R。ADR 0021 第 3 项。"""
    if not is_gas(name):
        raise PropertyError(f"{name} 不是理想气体，不能用 C_p − R")
    return heat_capacity_j_mol_k(name, temperature_k) - gas_constant_j_mol_k()


def density_kg_m3(temperature_k: float, pressure_pa: float, mole_fractions: dict[str, float]) -> float:
    """理想气体密度。附录 A1 的假定，S. 163。 Source: Hamel Eq.A.7, S. 164"""
    temperature = float(temperature_k)
    pressure = float(pressure_pa)
    if temperature <= 0 or pressure <= 0:
        raise PropertyError("温度和压力必须为正")
    total = 0.0
    mass = 0.0
    for name, fraction in mole_fractions.items():
        if fraction < 0:
            raise PropertyError(f"{name} 的摩尔分数为负")
        total += fraction
        mass += fraction * molar_mass_kg_mol(name)
    if total <= 0:
        raise PropertyError("摩尔分数之和必须为正")
    return pressure * (mass / total) / (gas_constant_j_mol_k() * temperature)


def volume_flow_m3_s(amount_mol_s: float, temperature_k: float, pressure_pa: float) -> float:
    """理想气体体积流，m³/s。 Source: Hamel Eq.A.7, S. 164"""
    if float(temperature_k) <= 0 or float(pressure_pa) <= 0:
        raise PropertyError("温度和压力必须为正")
    return float(amount_mol_s) * gas_constant_j_mol_k() * float(temperature_k) / float(pressure_pa)


def vaporization_enthalpy_j_mol(temperature_k: float) -> float:
    """h_v = h(H2O,g) − h(H2O,L)。液态水越界由焓函数报错。 Source: NASA/TP-2002-211556 Eq.(2), p. 2"""
    return enthalpy_j_mol("H2O", temperature_k) - enthalpy_j_mol("H2O_L", temperature_k)


def ash_specific_heat_j_per_kg_k() -> float:
    """950 J/(kg·K)，无温度区间。ADR 0018 第 6 项。"""
    return float(ash_record()["specific_heat_J_per_kg_K"])


def char_specific_enthalpy_j_per_kg(temperature_k: float) -> float:
    """焦炭按石墨 C(gr)。ADR 0018 第 5 项。"""
    return enthalpy_j_mol("C_gr", temperature_k) / molar_mass_kg_mol("C_gr")


def water_free_to_ash_free(mass_fraction_wf: float, ash_mass_fraction_wf: float) -> float:
    """x(waf) = x(wf) / (1 − x_ash,wf)。ADR 0015 第 2 项。"""
    ash = float(ash_mass_fraction_wf)
    if ash < 0 or ash >= 1:
        raise PropertyError("干基灰分必须在 [0, 1) 内")
    return float(mass_fraction_wf) / (1.0 - ash)


def fuel_formation_enthalpy_j_per_kg(
    heating_value_j_per_kg: float,
    *,
    carbon: float,
    hydrogen: float,
    nitrogen: float = 0.0,
    sulfur: float = 0.0,
    moisture: float = 0.0,
) -> float:
    """h_f,B = Σ ν_p h_f,p + H_u，单位 J/kg。质量分数与 H_u 同一基准。 产物为 CO2、H2O(g)、SO2、N2。燃料水分计入 H2O(g)。氧和灰不产生额外产物。 SO2 只用 298.15 K 生成焓。ADR 0018 第 4 项；ADR 0021 第 6 项。"""
    heating_value = float(heating_value_j_per_kg)
    if heating_value < 0:
        raise PropertyError("低位热值 H_u 按定义为正")
    def moles(mass_fraction: float, element: str, divisor: int = 1) -> float:
        mass = float(mass_fraction)
        if mass < 0:
            raise PropertyError("质量分数不能为负")
        weight = atomic_weight_kg(element)
        return mass / (divisor * weight)
    products = {
        "CO2": moles(carbon, "C"),
        "H2O": moles(hydrogen, "H", 2) + float(moisture) / molar_mass_kg_mol("H2O"),
        "SO2": moles(sulfur, "S"),
        "N2": moles(nitrogen, "N", 2),
    }
    if float(moisture) < 0:
        raise PropertyError("水分不能为负")
    enthalpy = heating_value
    for name, amount in products.items():
        enthalpy += amount * formation_enthalpy_j_mol(name)
    return enthalpy


def atomic_weight_kg(element: str) -> float:
    return atomic_weight(element) * 1e-3


def mixture_enthalpy_j(temperature_k: float, amounts_mol: dict[str, float]) -> float:
    """Σ N_j h_j，J。"""
    total = 0.0
    for name, amount in amounts_mol.items():
        if amount:
            total += float(amount) * enthalpy_j_mol(name, temperature_k)
    return total
