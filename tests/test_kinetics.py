# -*- coding: utf-8 -*-
"""第 5 章速率。期望值是按印刷式手算，不用旧代码输出。"""
from __future__ import annotations

import inspect
import math
from pathlib import Path

import pytest

from hamel.kinetics import (
    KineticsInputError,
    arrhenius_rate_constant,
    mechanism_factor,
    particle_temperature,
    r1_combustion_rate,
    r10_tar_oxidation_rate,
    r11_corella_rate,
    r11_serio_rate,
    r2_steam_gasification_rate,
    r3_boudouard_rate,
    r4_hydrogasification_rate,
    r5_co_oxidation_rate,
    r6_hydrogen_oxidation_rate,
    r7_methane_oxidation_rate,
    r8_equilibrium_constant,
    r8_water_gas_shift_rate,
    r9_methane_reforming_rate,
    reaction_stoichiometry,
)

# 测试传入的气体常数，ADR 0010 第 8 条。库内不写死这个数。
R_G = 8.314
COMPOSITION = {
    "C": {"C": 1.0},
    "O2": {"O": 2.0},
    "CO": {"C": 1.0, "O": 1.0},
    "CO2": {"C": 1.0, "O": 2.0},
    "H2": {"H": 2.0},
    "H2O": {"H": 2.0, "O": 1.0},
    "CH4": {"C": 1.0, "H": 4.0},
    "C6H6": {"C": 6.0, "H": 6.0},
    "C10H8": {"C": 10.0, "H": 8.0},
}


def _elements(coefficients):
    totals = {}
    for species, nu in coefficients.items():
        for element, count in COMPOSITION[species].items():
            totals[element] = totals.get(element, 0.0) + nu * count
    return totals


def test_arrhenius_factory_matches_hand_calculation():
    # k = 2.72 * exp(−188000 / (8.314 * 1100))
    temperature = 1100.0
    expected = 2.72 * math.exp(-188_000.0 / (R_G * temperature))
    assert arrhenius_rate_constant(
        2.72,
        temperature,
        activation_energy_j_mol=188_000.0,
        gas_constant_j_mol_k=R_G,
    ) == pytest.approx(expected)
    # 负活化能：k_b = 8.31e−7 * exp(−(−38900) / (R T))，比前因子大。
    negative = arrhenius_rate_constant(
        8.31e-7,
        temperature,
        activation_energy_j_mol=-38_900.0,
        gas_constant_j_mol_k=R_G,
    )
    assert negative == pytest.approx(8.31e-7 * math.exp(38_900.0 / (R_G * temperature)))
    assert negative > 8.31e-7


def test_arrhenius_exponent_is_only_in_the_factory():
    root = Path(__file__).resolve().parents[1] / "hamel" / "kinetics"
    for path in root.rglob("*.py"):
        if path.name == "arrhenius.py":
            continue
        assert "math.exp(-" not in path.read_text(encoding="utf-8")


def test_mechanism_factor_hand_calculation():
    # p_c = 2500 * exp(−6240 / 1000)
    # d = 2e−4 m 落在 5e−5 与 1e−3 之间：
    # Φ = [2 p_c + 2 − (p_c / 0.95) * (1e3 * d − 0.05)] / (p_c + 2)
    temperature = 1000.0
    diameter = 2.0e-4
    p_c = 2500.0 * math.exp(-6240.0 / temperature)
    expected = (2.0 * p_c + 2.0 - (p_c / 0.95) * (1.0e3 * diameter - 0.05)) / (p_c + 2.0)
    assert mechanism_factor(diameter, temperature) == pytest.approx(expected)
    # d = 1e−3 m 时中段公式给出 1，与 d > 1e−3 的分支相接。
    assert mechanism_factor(1.0e-3, temperature) == pytest.approx(1.0)
    assert mechanism_factor(2.0e-3, temperature) == 1.0


def test_r1_hand_calculation_and_zero_oxygen():
    # d_k = d_c = 2e−3 m，Φ = 1
    # k_ch = 1.22 * 1000 * exp(−10300 / 1000)
    # T_m = (1000 + 1100) / 2 = 1050
    # D = 3.13e−4 * (1050 / 1500)^1.75 * (101300 / 101300)
    # ν = μ / ρ，Re = u d / ν，Sc = μ / (ρ D)
    # Sh = 2 + 0.6 * Re^0.5 * Sc^0.33
    # k_d = Sh * D / d
    # k_1 = 1 / (1/(Φ k_d) + 1/k_ch)
    # 单颗粒进度 = π d² * k_1 * C_O2
    # N = m / (ρ_C * π d³ / 6)
    diameter = 2.0e-3
    temperature = 1000.0
    particle = 1100.0
    oxygen = 2.0
    pressure = 101_300.0
    velocity = 0.4
    viscosity = 4.0e-5
    density_gas = 0.3
    mass = 0.05
    density_char = 800.0
    phi = 1.0
    k_chemical = 1.22 * temperature * math.exp(-10_300.0 / temperature)
    film = 0.5 * (temperature + particle)
    diffusivity = 3.13e-4 * (film / 1500.0) ** 1.75 * (101_300.0 / pressure)
    kinematic = viscosity / density_gas
    reynolds = velocity * diameter / kinematic
    schmidt = viscosity / (density_gas * diffusivity)
    sherwood = 2.0 + 0.6 * reynolds**0.5 * schmidt**0.33
    k_diffusion = sherwood * diffusivity / diameter
    k_combined = 1.0 / (1.0 / (phi * k_diffusion) + 1.0 / k_chemical)
    per_particle = math.pi * diameter**2 * k_combined * oxygen
    count = mass / (density_char * math.pi * diameter**3 / 6.0)
    rate = r1_combustion_rate(
        mass,
        0.0,
        density_char,
        diameter,
        diameter,
        temperature,
        temperature,
        particle,
        temperature,
        oxygen,
        pressure,
        velocity,
        viscosity,
        density_gas,
    )
    assert rate.specific == pytest.approx(per_particle)
    assert rate.extent_mol_s == pytest.approx(per_particle * count)
    zero = r1_combustion_rate(
        mass,
        0.0,
        density_char,
        diameter,
        diameter,
        temperature,
        temperature,
        particle,
        temperature,
        0.0,
        pressure,
        velocity,
        viscosity,
        density_gas,
    )
    assert zero.specific == 0.0
    assert zero.extent_mol_s == 0.0


def test_particle_temperature_satisfies_energy_balance():
    # 式 (5.22)：R_C Δh = A α (T_P − T_∞) + A σ ε (T_P^4 − T_∞^4)
    # R_C 是印刷式的碳消耗速率，等于单颗粒进度的相反数。
    diameter = 1.0e-3
    gas = 1000.0
    enthalpy = -110.5e3
    alpha = 80.0
    emissivity = 0.85
    sigma = 5.6704e-8
    oxygen = 0.2
    pressure = 1.0e6
    velocity = 0.2
    viscosity = 4.0e-5
    density_gas = 0.4
    solved = particle_temperature(
        diameter,
        diameter,
        gas,
        gas,
        gas,
        oxygen,
        pressure,
        velocity,
        viscosity,
        density_gas,
        enthalpy,
        alpha,
        emissivity,
        sigma,
        gas,
        1.0e-6,
        40,
    )
    rate = r1_combustion_rate(
        0.01,
        0.0,
        700.0,
        diameter,
        diameter,
        gas,
        gas,
        solved,
        gas,
        oxygen,
        pressure,
        velocity,
        viscosity,
        density_gas,
    )
    area = math.pi * diameter**2
    left = -rate.specific * enthalpy
    right = area * alpha * (solved - gas) + area * sigma * emissivity * (solved**4 - gas**4)
    assert left == pytest.approx(right, abs=1.0e-6)


def test_r2_and_r4_hand_calculation():
    # R2：k = 3.42 * T * exp(−15600 / T)，单颗粒 π d² k C，再乘 N
    # R4：前因子为 3.42e−3，其余相同，浓度是 H₂
    diameter = 1.5e-3
    temperature = 1150.0
    mass = 0.02
    density = 650.0
    concentration = 3.0
    count = mass / (density * math.pi * diameter**3 / 6.0)
    area = math.pi * diameter**2
    k_r2 = 3.42 * temperature * math.exp(-15_600.0 / temperature)
    rate_r2 = r2_steam_gasification_rate(mass, 0.0, density, diameter, concentration, temperature)
    assert rate_r2.specific == pytest.approx(area * k_r2 * concentration)
    assert rate_r2.extent_mol_s == pytest.approx(rate_r2.specific * count)
    k_r4 = 3.42e-3 * temperature * math.exp(-15_600.0 / temperature)
    rate_r4 = r4_hydrogasification_rate(mass, 0.0, density, diameter, concentration, temperature)
    assert rate_r4.specific == pytest.approx(area * k_r4 * concentration)
    assert rate_r4.extent_mol_s == pytest.approx(rate_r4.specific * count)


def test_r3_hand_calculation_scales_with_mass_not_diameter():
    # k_a = 2.72 exp(−188000/(R T))，k_b、k_c 的活化能为负
    # R* = k_a P_CO2 / (1 + k_b P_CO2 + k_c P_CO)，单位 1/s
    # 进度 = R* * m / M_C。式中没有 d_k。
    temperature = 1100.0
    p_co2 = 2.0e4
    p_co = 1.0e4
    mass = 0.2
    molar_mass = 0.012
    k_a = 2.72 * math.exp(-188_000.0 / (R_G * temperature))
    k_b = 8.31e-7 * math.exp(38_900.0 / (R_G * temperature))
    k_c = 2.088e-14 * math.exp(24_500.0 / (R_G * temperature))
    specific = k_a * p_co2 / (1.0 + k_b * p_co2 + k_c * p_co)
    rate = r3_boudouard_rate(mass, molar_mass, p_co2, p_co, temperature, R_G)
    assert rate.specific == pytest.approx(specific)
    assert rate.extent_mol_s == pytest.approx(specific * mass / molar_mass)
    doubled = r3_boudouard_rate(2.0 * mass, molar_mass, p_co2, p_co, temperature, R_G)
    assert doubled.extent_mol_s == pytest.approx(2.0 * rate.extent_mol_s)
    assert doubled.specific == pytest.approx(rate.specific)
    assert "d_k_m" not in inspect.signature(r3_boudouard_rate).parameters
    assert "m_ash_kg" not in inspect.signature(r3_boudouard_rate).parameters
    assert "rho_p_kg_m3" not in inspect.signature(r3_boudouard_rate).parameters


def test_surface_rates_count_char_plus_ash_at_the_same_density():
    # ADR 0017 第 19 项：N = (m_C + m_A) / (ρ_P · π d³ / 6)。
    # 焦炭与灰各 0.01 kg、ρ_P 不变时，cell 速率是只计焦炭的两倍；单颗粒比速率不变。
    diameter = 1.0e-3
    density = 800.0
    temperature = 1100.0
    char_only = r2_steam_gasification_rate(0.01, 0.0, density, diameter, 1.0, temperature)
    with_ash = r2_steam_gasification_rate(0.01, 0.01, density, diameter, 1.0, temperature)
    assert with_ash.specific == pytest.approx(char_only.specific)
    assert with_ash.extent_mol_s == pytest.approx(2.0 * char_only.extent_mol_s)
    char_only_r4 = r4_hydrogasification_rate(0.01, 0.0, density, diameter, 1.0, temperature)
    with_ash_r4 = r4_hydrogasification_rate(0.01, 0.01, density, diameter, 1.0, temperature)
    assert with_ash_r4.extent_mol_s == pytest.approx(2.0 * char_only_r4.extent_mol_s)
    shared = (diameter, diameter, 1000.0, 1000.0, 1000.0, 1000.0, 1.0, 1.0e5, 0.2, 4.0e-5, 0.3)
    char_only_r1 = r1_combustion_rate(0.01, 0.0, density, *shared)
    with_ash_r1 = r1_combustion_rate(0.01, 0.01, density, *shared)
    assert with_ash_r1.specific == pytest.approx(char_only_r1.specific)
    assert with_ash_r1.extent_mol_s == pytest.approx(2.0 * char_only_r1.extent_mol_s)


def test_r5_hand_calculation_for_both_phases():
    # 气泡：k = 1.91e6 exp(−8056/T)，dC/dt = k C_CO C_O2^0.3 C_H2O^0.5
    # 悬浮：k = 1.8 T^0.5，不乘 exp，也不乘水。进度再乘气体体积。
    temperature = 1050.0
    volume = 0.03
    carbon_monoxide = 1.5
    oxygen = 0.4
    water = 0.8
    k_bubble = 1.91e6 * math.exp(-8056.0 / temperature)
    specific_bubble = k_bubble * carbon_monoxide * oxygen**0.3 * water**0.5
    bubble = r5_co_oxidation_rate(
        "bubble", carbon_monoxide, oxygen, temperature, volume, water
    )
    assert bubble.specific == pytest.approx(specific_bubble)
    assert bubble.extent_mol_s == pytest.approx(specific_bubble * volume)
    specific_suspension = 1.8 * temperature**0.5 * carbon_monoxide * oxygen**0.5
    dry = r5_co_oxidation_rate("suspension", carbon_monoxide, oxygen, temperature, volume)
    wet = r5_co_oxidation_rate(
        "suspension", carbon_monoxide, oxygen, temperature, volume, 9.0
    )
    assert dry.specific == pytest.approx(specific_suspension)
    assert wet.specific == pytest.approx(dry.specific)
    assert dry.extent_mol_s == pytest.approx(specific_suspension * volume)


def test_r6_and_r7_hand_calculation():
    # R6：k = 1.631e9 T^1.5 exp(−3420/T)，dC/dt = k C_H2^1.5 C_O2
    # R7：k = 3.552e11 / T * exp(−15700/T)，dC/dt = k C_CH4 C_O2
    temperature = 1000.0
    volume = 0.02
    hydrogen = 0.3
    oxygen = 0.5
    methane = 0.15
    k_r6 = 1.631e9 * temperature**1.5 * math.exp(-3420.0 / temperature)
    specific_r6 = k_r6 * hydrogen**1.5 * oxygen
    rate_r6 = r6_hydrogen_oxidation_rate(hydrogen, oxygen, temperature, volume)
    assert rate_r6.specific == pytest.approx(specific_r6)
    assert rate_r6.extent_mol_s == pytest.approx(specific_r6 * volume)
    k_r7 = 3.552e11 / temperature * math.exp(-15_700.0 / temperature)
    specific_r7 = k_r7 * methane * oxygen
    rate_r7 = r7_methane_oxidation_rate(methane, oxygen, temperature, volume)
    assert rate_r7.specific == pytest.approx(specific_r7)
    assert rate_r7.extent_mol_s == pytest.approx(specific_r7 * volume)


def test_r8_equilibrium_constant_is_synthetic_formula_only():
    # 合成数据，仅检验公式。Σνμ = 0 时 K = 1；Σνμ 增大时 K 下降。
    # Σνμ = −μ_CO − μ_H2O + μ_H2 + μ_CO2。这里只改 μ_H2。
    temperature = 1100.0
    assert r8_equilibrium_constant(0.0, 0.0, 0.0, 0.0, temperature, R_G) == pytest.approx(1.0)
    high = r8_equilibrium_constant(0.0, 0.0, R_G * temperature, 0.0, temperature, R_G)
    low = r8_equilibrium_constant(0.0, 0.0, -R_G * temperature, 0.0, temperature, R_G)
    assert high == pytest.approx(math.exp(-1.0))
    assert low == pytest.approx(math.exp(1.0))
    assert low > 1.0 > high


def test_r8_is_zero_at_equilibrium_and_changes_sign():
    # K = 1（四个化学势都是 0）。y_eq = y_CO2 y_H2 / (K y_H2O)
    # y = (0.2, 0.2, 0.2, 0.2) 时 y_CO = y_eq，式 (5.45) 的差分为 0。
    # p = 1e6 Pa = 1e6/101325 atm
    # 进度 = 0.02 * 2.77e8 * (y_CO − y_eq) * exp(−13971/T) * p^(0.5 − p/250) * exp(−8.91 + 5553/T) * m_ash
    temperature = 1100.0
    pressure = 1.0e6
    ash = 0.4
    shared = dict(
        mu_co_j_mol=0.0,
        mu_h2o_j_mol=0.0,
        mu_h2_j_mol=0.0,
        mu_co2_j_mol=0.0,
        temperature_k=temperature,
        pressure_pa=pressure,
        gas_constant_j_mol_k=R_G,
        m_ash_kg=ash,
    )
    balanced = r8_water_gas_shift_rate(0.2, 0.2, 0.2, 0.2, **shared)
    assert balanced.specific == pytest.approx(0.0, abs=1e-12)
    assert balanced.extent_mol_s == pytest.approx(0.0, abs=1e-12)
    pressure_atm = pressure / 101_325.0
    factor = (
        0.02
        * 2.77e8
        * math.exp(-13_971.0 / temperature)
        * pressure_atm ** (0.5 - pressure_atm / 250.0)
        * math.exp(-8.91 + 5553.0 / temperature)
    )
    # y_eq = 0.2 * 0.2 / (1 * 0.2) = 0.2。y_CO = 0.5 时差分为 0.3。
    rich = r8_water_gas_shift_rate(0.5, 0.2, 0.2, 0.2, **shared)
    assert rich.specific == pytest.approx(factor * 0.3)
    assert rich.extent_mol_s == pytest.approx(rich.specific * ash)
    assert rich.specific > 0.0
    lean = r8_water_gas_shift_rate(0.1, 0.2, 0.2, 0.2, **shared)
    assert lean.specific == pytest.approx(factor * (0.1 - 0.2))
    assert lean.specific < 0.0


def test_r9_hand_calculation():
    # 印刷式 R_CH4 = −6.113e−2 * exp(−137327/(R T)) * P_CH4
    # 进度比速率是这个数的相反数，再乘灰质量。P 以 Pa 计。
    temperature = 1164.0
    pressure = 8_000.0
    ash = 0.25
    specific = 6.113e-2 * math.exp(-137_327.0 / (R_G * temperature)) * pressure
    rate = r9_methane_reforming_rate(pressure, temperature, R_G, ash)
    assert rate.specific == pytest.approx(specific)
    assert rate.extent_mol_s == pytest.approx(specific * ash)
    assert rate.specific > 0.0


def test_r10_hand_calculation():
    # 芳烃：k_10 = 20700，E/R = 9650 K
    # 进度比速率 = k_10 exp(−9650/T) * T * P^0.3 * C_tar^0.5 * C_O2
    temperature = 1100.0
    pressure = 1.2e6
    tar = 0.04
    oxygen = 0.2
    volume = 0.015
    specific = (
        20_700.0
        * math.exp(-9650.0 / temperature)
        * temperature
        * pressure**0.3
        * tar**0.5
        * oxygen
    )
    rate = r10_tar_oxidation_rate(tar, oxygen, temperature, pressure, volume)
    assert rate.specific == pytest.approx(specific)
    assert rate.extent_mol_s == pytest.approx(specific * volume)


def test_tar_secondary_hand_calculation():
    # Serio 第一对：k = 5.42e4 exp(−100500/(R T))，比速率 = k C，再乘气体体积。
    # Corella 高温段（T ≥ 740+273.15）：k = 0.7 exp(−63100/(R T))，再乘催化剂质量。
    # 低温段：k = 7.5e−4 exp(−5900/(R T))。740 °C 用高温段。
    temperature = 1164.0
    tar = 0.05
    volume = 0.01
    catalyst = 1.2
    k_serio = 5.42e4 * math.exp(-100.5e3 / (R_G * temperature))
    serio = r11_serio_rate(tar, temperature, R_G, volume)
    assert serio.specific == pytest.approx(k_serio * tar)
    assert serio.extent_mol_s == pytest.approx(k_serio * tar * volume)
    k_high = 0.7 * math.exp(-63.1e3 / (R_G * temperature))
    corella = r11_corella_rate(tar, temperature, R_G, catalyst)
    assert corella.specific == pytest.approx(k_high * tar)
    assert corella.extent_mol_s == pytest.approx(k_high * tar * catalyst)
    split = 740.0 + 273.15
    at_split = r11_corella_rate(tar, split, R_G, catalyst)
    high_at_split = 0.7 * math.exp(-63.1e3 / (R_G * split)) * tar
    assert at_split.specific == pytest.approx(high_at_split)
    low = 1000.0
    k_low = 7.5e-4 * math.exp(-5.9e3 / (R_G * low))
    cool = r11_corella_rate(tar, low, R_G, catalyst)
    assert cool.specific == pytest.approx(k_low * tar)
    assert r11_corella_rate(tar, 600.0 + 273.15, R_G, catalyst).specific > 0.0
    assert r11_corella_rate(tar, 900.0 + 273.15, R_G, catalyst).specific > 0.0
    with pytest.raises(KineticsInputError):
        r11_corella_rate(tar, 600.0 + 273.15 - 1.0, R_G, catalyst)
    with pytest.raises(KineticsInputError):
        r11_corella_rate(tar, 900.0 + 273.15 + 1.0, R_G, catalyst)


def test_every_reaction_conserves_elements():
    cases = [("R2", None, None), ("R3", None, None), ("R4", None, None)]
    cases += [("R5", None, None), ("R6", None, None), ("R7", None, None)]
    cases += [("R8", None, None), ("R9", None, None)]
    for reaction, phi, tar in cases:
        totals = _elements(reaction_stoichiometry(reaction, phi, tar))
        assert totals == pytest.approx({element: 0.0 for element in totals})
    for phi in (1.0, 1.5, 2.0, mechanism_factor(2.0e-4, 1000.0)):
        totals = _elements(reaction_stoichiometry("R1", phi_c=phi))
        assert totals == pytest.approx({element: 0.0 for element in totals})
    for reaction in ("R10", "R11"):
        for tar in ("C6H6", "C10H8"):
            totals = _elements(reaction_stoichiometry(reaction, tar=tar))
            assert totals == pytest.approx({element: 0.0 for element in totals})


def test_missing_required_argument_has_no_default():
    with pytest.raises(TypeError):
        particle_temperature(1.0e-3, 1.0e-3, 1000.0, 1000.0, 1000.0, 0.1, 1.0e5, 0.1, 4.0e-5, 0.3, -1.0e5, 10.0, 0.85)
    with pytest.raises(KineticsInputError):
        r8_water_gas_shift_rate(0.2, 0.0, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 1100.0, 1.0e5, R_G, 0.1)
