"""低压气体粘度、导热系数、Fuller 二元扩散与相间 D_g。"""
from __future__ import annotations

import math

from hamel.properties.data import (
    PropertyError,
    PropertyRangeError,
    _bundle,
    atoms,
    molar_mass_g_mol,
    transport_species,
)
from hamel.properties.thermo import heat_capacity_volume_j_mol_k

# 式 (9-3.8)、(9-3.9)：η 为 μP，M 为 g/mol，σ 为 Å。1 μP = 1e-7 Pa·s。
_CE_PREFACTOR = 26.69e-7
# 式 (9-4.3)，PDF 460。0.3 ≤ T* ≤ 100。
_NEUFELD = (1.16145, 0.14874, 0.52487, 0.77320, 2.16178, 2.43787)
_BAR_PA = 1.0e5
def _positive(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise PropertyError(f"{name} 必须为正")
    return number


def collision_integral(reduced_temperature: float) -> float:
    """Neufeld 等 (1972)，Poling 式 (9-4.3)。"""
    reduced = _positive("T*", reduced_temperature)
    if not 0.3 <= reduced <= 100.0:
        raise PropertyRangeError(f"T* = {reduced} 不在式 (9-4.3) 的 0.3–100 内")
    factor_a, factor_b, factor_c, factor_d, factor_e, factor_f = _NEUFELD
    return (
        factor_a * reduced ** (-factor_b)
        + factor_c * math.exp(-factor_d * reduced)
        + factor_e * math.exp(-factor_f * reduced)
    )


def lennard_jones(name: str) -> tuple[float, float]:
    """返回 σ（Å）与 ε/k（K）。萘按式 (9-4.7)、(9-4.8)。解释性选择：ADR 0021 第 1 项。"""
    record = transport_species(name)
    sigma = record["sigma_A"]
    epsilon = record["epsilon_K"]
    if sigma is not None and epsilon is not None:
        return float(sigma), float(epsilon)
    method = _bundle()["transport"]["naphthalene_lennard_jones"]
    epsilon = float(record["tc_K"]) / float(method["epsilon_factor"])
    sigma = float(method["sigma_factor"]) * float(record["vc_cm3_mol"]) ** (1.0 / 3.0)
    return sigma, epsilon


def fuller_volume(name: str) -> float:
    """表 11-1。简单分子用印出的值，其余用原子增量加芳香环。"""
    record = transport_species(name)
    if record["fuller_volume"] is not None:
        return float(record["fuller_volume"])
    increments = _bundle()["transport"]["fuller_increments"]
    total = float(record["aromatic_rings"]) * float(increments["aromatic_ring"])
    for element, count in atoms(name).items():
        if element not in increments:
            raise PropertyError(f"{name} 的元素 {element} 没有 Fuller 增量")
        total += count * float(increments[element])
    return total
def _require_reduced_pressure(name: str, pressure_pa: float) -> None:
    """ADR 0021 第 9 项：不作压力修正；P_r ≥ 2 报错。"""
    critical = float(transport_species(name)["pc_bar"]) * _BAR_PA
    if float(pressure_pa) / critical >= 2.0:
        raise PropertyRangeError(f"{name} 的对比压力不小于 2")


def viscosity_pa_s(name: str, temperature_k: float, pressure_pa: float) -> float:
    """Chapman–Enskog 式 (9-3.9) 加式 (9-4.3)。 Source: Poling Eq.(9-3.9), PDF 459"""
    temperature = _positive("temperature_k", temperature_k)
    pressure = _positive("pressure_pa", pressure_pa)
    _require_reduced_pressure(name, pressure)
    sigma, epsilon = lennard_jones(name)
    omega = collision_integral(temperature / epsilon)
    return _CE_PREFACTOR * math.sqrt(molar_mass_g_mol(name) * temperature) / (sigma**2 * omega)
def _wilke_phi(mu_i: float, mu_j: float, mass_i: float, mass_j: float) -> float:
    """式 (9-5.14)。"""
    return (1.0 + math.sqrt(mu_i / mu_j) * (mass_j / mass_i) ** 0.25) ** 2 / math.sqrt(
        8.0 * (1.0 + mass_i / mass_j)
    )


def mixture_viscosity_pa_s(
    mole_fractions: dict[str, float], temperature_k: float, pressure_pa: float,
) -> float:
    """Wilke 式 (9-5.13)。 Source: Poling Eq.(9-5.13), PDF 476"""
    names = [name for name, fraction in mole_fractions.items() if fraction > 0]
    if not names:
        raise PropertyError("混合物没有正的摩尔分数")
    viscosity = {name: viscosity_pa_s(name, temperature_k, pressure_pa) for name in names}
    mass = {name: molar_mass_g_mol(name) for name in names}
    total = 0.0
    for name in names:
        denominator = 0.0
        for other in names:
            denominator += mole_fractions[other] * _wilke_phi(
                viscosity[name], viscosity[other], mass[name], mass[other],
            )
        total += mole_fractions[name] * viscosity[name] / denominator
    return total


def chung_conductivity_w_m_k(
    viscosity_pa_s_value: float,
    molar_mass_kg_mol: float,
    cv_j_mol_k: float,
    acentric_factor: float,
    reduced_temperature: float,
    gas_constant_j_mol_k: float,
) -> float:
    """Chung 式 (10-3.14)。β 一律用 ω 的关联，极性组分不另取 0.758。解释性选择：ADR 0021 第 3 项。"""
    viscosity = _positive("viscosity", viscosity_pa_s_value)
    mass = _positive("molar_mass", molar_mass_kg_mol)
    heat_capacity = _positive("Cv", cv_j_mol_k)
    reduced = _positive("Tr", reduced_temperature)
    gas_constant = _positive("R", gas_constant_j_mol_k)
    alpha = heat_capacity / gas_constant - 1.5
    beta = 0.7862 - 0.7109 * acentric_factor + 1.3168 * acentric_factor**2
    collision_number = 2.0 + 10.5 * reduced**2
    numerator = 0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * collision_number
    denominator = 0.6366 + beta * collision_number + 1.061 * alpha * beta
    psi = 1.0 + alpha * numerator / denominator
    return 3.75 * psi * viscosity * gas_constant / mass


def thermal_conductivity_w_m_k(name: str, temperature_k: float, pressure_pa: float) -> float:
    """纯气体导热系数。粘度用 Chapman–Enskog，C_v = C_p − R。 Source: Poling Eq.(10-3.14), PDF 565"""
    from hamel.properties.thermo import gas_constant_j_mol_k
    temperature = _positive("temperature_k", temperature_k)
    record = transport_species(name)
    viscosity = viscosity_pa_s(name, temperature, pressure_pa)
    return chung_conductivity_w_m_k(
        viscosity,
        molar_mass_g_mol(name) * 1e-3,
        heat_capacity_volume_j_mol_k(name, temperature),
        float(record["omega"]),
        temperature / float(record["tc_K"]),
        gas_constant_j_mol_k(),
    )
def _translational_ratio(name_i: str, name_j: str, temperature_k: float) -> float:
    """式 (10-6.5)。Γ 用式 (10-3.9)。"""
    def gamma_times_factor(name: str) -> tuple[float, float]:
        record = transport_species(name)
        critical_temperature = float(record["tc_K"])
        gamma = 210.0 * (
            critical_temperature * molar_mass_g_mol(name) ** 3 / float(record["pc_bar"]) ** 4
        ) ** (1.0 / 6.0)
        reduced = temperature_k / critical_temperature
        factor = math.exp(0.0464 * reduced) - math.exp(-0.2412 * reduced)
        return gamma, factor
    gamma_i, factor_i = gamma_times_factor(name_i)
    gamma_j, factor_j = gamma_times_factor(name_j)
    return (gamma_j * factor_i) / (gamma_i * factor_j)


def mixture_thermal_conductivity_w_m_k(
    mole_fractions: dict[str, float], temperature_k: float, pressure_pa: float,
) -> float:
    """Wassiljewa 式 (10-6.1)，A_ij 用式 (10-6.2)、(10-6.5)，ε = 1。解释性选择：ADR 0021 第 3 项。 Source: Poling Eq.(10-6.1), PDF 583"""
    temperature = _positive("temperature_k", temperature_k)
    names = [name for name, fraction in mole_fractions.items() if fraction > 0]
    if not names:
        raise PropertyError("混合物没有正的摩尔分数")
    conductivity = {
        name: thermal_conductivity_w_m_k(name, temperature, pressure_pa) for name in names
    }
    mass = {name: molar_mass_g_mol(name) for name in names}
    total = 0.0
    for name in names:
        denominator = 0.0
        for other in names:
            ratio = _translational_ratio(name, other, temperature)
            interaction = (
                1.0 + math.sqrt(ratio) * (mass[name] / mass[other]) ** 0.25
            ) ** 2 / math.sqrt(8.0 * (1.0 + mass[name] / mass[other]))
            denominator += mole_fractions[other] * interaction
        total += mole_fractions[name] * conductivity[name] / denominator
    return total


def fuller_binary_diffusivity_m2_s(
    temperature_k: float,
    pressure_pa: float,
    molar_mass_a_g_mol: float,
    molar_mass_b_g_mol: float,
    volume_a: float,
    volume_b: float,
) -> float:
    """Fuller 式 (11-4.4)。返回 m²/s。P 在式中为 bar。 Source: Poling Eq.(11-4.4), PDF 633"""
    temperature = _positive("temperature_k", temperature_k)
    pressure_bar = _positive("pressure_pa", pressure_pa) / _BAR_PA
    mass_a = _positive("M_A", molar_mass_a_g_mol)
    mass_b = _positive("M_B", molar_mass_b_g_mol)
    reduced_mass = 2.0 / (1.0 / mass_a + 1.0 / mass_b)
    denominator = (
        pressure_bar
        * math.sqrt(reduced_mass)
        * (volume_a ** (1.0 / 3.0) + volume_b ** (1.0 / 3.0)) ** 2
    )
    cm2_s = 0.00143 * temperature**1.75 / denominator
    return cm2_s * 1.0e-4


def binary_diffusivity_m2_s(
    name_a: str, name_b: str, temperature_k: float, pressure_pa: float,
) -> float:
    """两组分的 Fuller 二元扩散系数。P_r ≥ 2 报错。不用于相间交换。"""
    pressure = _positive("pressure_pa", pressure_pa)
    _require_reduced_pressure(name_a, pressure)
    _require_reduced_pressure(name_b, pressure)
    return fuller_binary_diffusivity_m2_s(
        temperature_k, pressure,
        molar_mass_g_mol(name_a), molar_mass_g_mol(name_b),
        fuller_volume(name_a), fuller_volume(name_b),
    )


def interphase_diffusivity_m2_s(temperature_k: float, pressure_pa: float) -> float:
    """相间交换的 D_g，与组分无关。解释性选择：ADR 0021 第 4 项，取 cell 气体温度。 Source: Hamel Eq.5.20, S. 87"""
    temperature = _positive("temperature_k", temperature_k)
    pressure = _positive("pressure_pa", pressure_pa)
    return 3.13e-4 * (temperature / 1500.0) ** 1.75 * (101_300.0 / pressure)
