"""物性与附录 A1。期望值来自 NASA 表、Poling 算例或按原文手算。"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from hamel.kinetics import r8_equilibrium_constant
from hamel.properties import (
    GAS_SPECIES,
    EquilibriumError,
    PropertyRangeError,
    adiabatic_gas_equilibrium,
    ash_specific_heat_j_per_kg_k,
    balance_deviation,
    chemical_potential_j_mol,
    chung_conductivity_w_m_k,
    collision_integral,
    enthalpy_j_mol,
    entropy_j_mol_k,
    formation_enthalpy_j_mol,
    fuel_formation_enthalpy_j_per_kg,
    fuller_binary_diffusivity_m2_s,
    fuller_volume,
    gas_constant_j_mol_k,
    gibbs_minimization,
    heat_capacity_j_mol_k,
    initial_lambdas,
    interphase_diffusivity_m2_s,
    is_converged,
    lennard_jones,
    molar_mass_kg_mol,
    thermal_conductivity_w_m_k,
    viscosity_pa_s,
    water_free_to_ash_free,
)
from hamel.properties.data import nasa_species
from hamel.properties.equilibrium import _c_values, _stoichiometry, _total_moles

_ROOT = Path(__file__).resolve().parents[1]
_SPECIES = (
    "Ar", "CH4", "CO", "CO2", "C2H4", "C6H6", "C10H8", "H2", "H2O", "H2S",
    "NH3", "N2", "O2", "SO2", "C_gr", "H2O_L",
)


def test_nasa_table_at_298() -> None:
    """生成焓与 c_p(298.15 K) 对上附录 B。"""
    for name in _SPECIES:
        table = nasa_species(name)["table_298"]
        expected_hf = table["formation_enthalpy_kJ_mol"] * 1000.0
        assert formation_enthalpy_j_mol(name) == pytest.approx(expected_hf, abs=1.0)
        assert enthalpy_j_mol(name, 298.15) == pytest.approx(expected_hf, abs=5.0)
        assert heat_capacity_j_mol_k(name, 298.15) == pytest.approx(table["cp_J_mol_K"], abs=0.02)


def test_chemical_potential_temperature_derivative_is_minus_entropy() -> None:
    gas_constant = gas_constant_j_mol_k()
    temperature = 1100.0
    step = 0.5
    for name in ("CO", "H2O", "H2", "CO2", "C_gr"):
        ahead = chemical_potential_j_mol(name, temperature + step)
        behind = chemical_potential_j_mol(name, temperature - step)
        derivative = (ahead - behind) / (2 * step)
        assert derivative == pytest.approx(-entropy_j_mol_k(name, temperature), rel=1e-4, abs=1e-3)
        assert gas_constant > 0


def test_standard_pressure_conversion() -> None:
    """气体熵与 μ 从 1 bar 换到 101300 Pa；焓不变。ADR 0018 第 1 项。"""
    temperature = 1000.0
    gas_constant = gas_constant_j_mol_k()
    shift = gas_constant * math.log(101_300.0 / 100_000.0)
    for name in GAS_SPECIES:
        entropy_bar = entropy_j_mol_k(name, temperature, standard_pressure_pa=100_000.0)
        entropy_p0 = entropy_j_mol_k(name, temperature)
        assert entropy_p0 - entropy_bar == pytest.approx(-shift, rel=1e-12)
        mu_bar = chemical_potential_j_mol(name, temperature, standard_pressure_pa=100_000.0)
        mu_p0 = chemical_potential_j_mol(name, temperature)
        assert mu_p0 - mu_bar == pytest.approx(gas_constant * temperature * math.log(101_300.0 / 100_000.0))
    liquid_bar = entropy_j_mol_k("H2O_L", 300.0, standard_pressure_pa=100_000.0)
    liquid_p0 = entropy_j_mol_k("H2O_L", 300.0, standard_pressure_pa=101_300.0)
    assert liquid_p0 == pytest.approx(liquid_bar)


def test_liquid_water_rejects_temperature_outside_its_interval() -> None:
    with pytest.raises(PropertyRangeError):
        enthalpy_j_mol("H2O_L", 250.0)
    with pytest.raises(PropertyRangeError):
        enthalpy_j_mol("H2O_L", 601.0)
    enthalpy_j_mol("H2O_L", 298.15)
    enthalpy_j_mol("H2O_L", 600.0)


def test_graphite_combustion_formation_enthalpy_is_zero() -> None:
    """石墨：H_u = −h_f,CO2 / M，产物 1 mol CO2，反算 h_f,B = 0。ADR 0018 第 4 项。"""
    molar_mass = molar_mass_kg_mol("C_gr")
    heating_value = -formation_enthalpy_j_mol("CO2") / molar_mass
    enthalpy = fuel_formation_enthalpy_j_per_kg(heating_value, carbon=1.0, hydrogen=0.0)
    assert enthalpy == pytest.approx(0.0, abs=1e-6)


def test_fuel_formation_enthalpy_hand_calculation() -> None:
    """干基、无水分。h_f,B = ν_CO2 h_CO2 + ν_H2O h_H2O + ν_SO2 h_SO2 + ν_N2 h_N2 + H_u。

    ν_CO2 = 0.80/0.012011，ν_H2O = 0.05/(2·0.001008)，
    ν_SO2 = 0.03/0.03206，ν_N2 = 0.02/(2·0.014007)，H_u = 2.5e7 J/kg。
    """
    heating_value = 25_000_000.0
    carbon, hydrogen, sulfur, nitrogen = 0.80, 0.05, 0.03, 0.02
    expected = heating_value
    expected += (carbon / 0.012011) * formation_enthalpy_j_mol("CO2")
    expected += (hydrogen / (2 * 0.001008)) * formation_enthalpy_j_mol("H2O")
    expected += (sulfur / 0.03206) * formation_enthalpy_j_mol("SO2")
    expected += (nitrogen / (2 * 0.014007)) * formation_enthalpy_j_mol("N2")
    got = fuel_formation_enthalpy_j_per_kg(
        heating_value, carbon=carbon, hydrogen=hydrogen, sulfur=sulfur, nitrogen=nitrogen,
    )
    assert got == pytest.approx(expected, rel=1e-12)
    assert water_free_to_ash_free(0.638, 0.043) == pytest.approx(0.638 / 0.957)


def test_ash_heat_capacity_is_the_recorded_constant() -> None:
    assert ash_specific_heat_j_per_kg_k() == 950.0


def test_hydrogen_pseudo_convergence_is_rejected() -> None:
    """只含 H2、b_H = 1 mol、y = 0.5 时元素平衡成立，归一化不成立。ADR 0017 第 16 项。"""
    stoichiometry = np.array([[2.0]])
    element_amount = np.array([1.0])
    mole_fraction = np.array([0.5])
    total = _total_moles(stoichiometry, element_amount, mole_fraction)
    amounts = total * mole_fraction
    deviation = balance_deviation(stoichiometry, element_amount, amounts)
    residual = abs(float(mole_fraction[0]) - 1.0)
    assert total == pytest.approx(1.0)
    assert deviation == pytest.approx(0.0, abs=1e-12)
    assert not is_converged(residual, deviation)


def test_initial_lambdas_satisfy_equation_a39() -> None:
    names = ["H2", "H2O", "CO", "CO2"]
    elements = ["H", "O", "C"]
    b_plus = np.array([4.0, 2.0, 1.0])
    stoichiometry = _stoichiometry(names, elements)
    temperature, pressure = 1200.0, 500_000.0
    c_values = _c_values(names, temperature, pressure)
    mole_fraction = np.full(len(names), 1.0 / len(names))
    total = _total_moles(stoichiometry, b_plus, mole_fraction)
    lambdas = initial_lambdas(stoichiometry, b_plus, c_values, mole_fraction, total)
    bracket = mole_fraction * (1.0 - np.log(mole_fraction))
    bracket = bracket + mole_fraction * (stoichiometry.T @ lambdas - c_values)
    left = stoichiometry @ bracket
    assert left == pytest.approx(b_plus / total, rel=1e-10, abs=1e-10)


def test_oxygen_deficient_feed_makes_no_oxygenated_gas() -> None:
    result = gibbs_minimization(
        1000.0, 101_300.0, {"H": 2.0, "O": 0.0, "C": 0.0},
        ["H2", "H2O", "CO", "CO2", "O2"],
    )
    for name in ("H2O", "CO", "CO2", "O2"):
        assert result.mole_fractions[name] == 0.0
        assert result.amounts_mol[name] == 0.0
    assert result.mole_fractions["H2"] == pytest.approx(1.0, abs=1e-8)
    oxygen = 0.0
    for name, amount in result.amounts_mol.items():
        oxygen += amount * nasa_species(name)["atoms"].get("O", 0)
    assert oxygen == 0.0


def test_equilibrium_conserves_elements() -> None:
    feed = {"C": 1.0, "O": 2.0, "H": 2.0, "N": 4.0}
    result = gibbs_minimization(1100.0, 500_000.0, feed)
    totals = {element: 0.0 for element in feed}
    for name, amount in result.amounts_mol.items():
        for element, count in nasa_species(name)["atoms"].items():
            if element in totals:
                totals[element] += amount * count
    for element, amount in feed.items():
        assert totals[element] == pytest.approx(amount, rel=1e-8, abs=1e-8)
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-8)


def test_water_gas_shift_matches_figure_5_11() -> None:
    """用 P_0 下的 μ^0 按式 (5.43) 计算，与数字化读数比较。"""
    payload = json.loads((_ROOT / "data" / "digitized" / "fig_5_11.json").read_text(encoding="utf-8"))
    gas_constant = gas_constant_j_mol_k()
    reading_error = payload["est_error"]["K_eq"]
    temperature_error = payload["est_error"]["T_K"]
    worst = 0.0
    for temperature, reading in payload["points"]:
        def equilibrium_constant(point: float) -> float:
            return r8_equilibrium_constant(
                chemical_potential_j_mol("CO", point),
                chemical_potential_j_mol("H2O", point),
                chemical_potential_j_mol("H2", point),
                chemical_potential_j_mol("CO2", point),
                point,
                gas_constant,
            )

        calculated = equilibrium_constant(temperature)
        low = equilibrium_constant(temperature - temperature_error)
        high = equilibrium_constant(temperature + temperature_error)
        band_low = min(low, high) - reading_error
        band_high = max(low, high) + reading_error
        worst = max(worst, abs(calculated - reading))
        assert band_low <= reading <= band_high
    print(f"fig 5.11 max |K_calc - K_read| = {worst:.4f}")


def test_naphthalene_fuller_volume_and_lennard_jones() -> None:
    assert fuller_volume("C10H8") == pytest.approx(140.88, abs=1e-9)
    assert fuller_volume("C6H6") == pytest.approx(90.96, abs=1e-9)
    sigma, epsilon = lennard_jones("C10H8")
    assert epsilon == pytest.approx(748.4 / 1.2593, rel=1e-12)
    assert sigma == pytest.approx(0.809 * 413.0 ** (1.0 / 3.0), rel=1e-12)


def test_poling_examples() -> None:
    """例 10-1 的 Ω_v、Chung 导热，以及例 11-3 的 Fuller 结果。"""
    assert collision_integral(1.020) == pytest.approx(1.576, abs=0.001)
    conductivity = chung_conductivity_w_m_k(
        8.72e-6, 72.151e-3, 135.8, 0.229, 373.0 / 460.39, 8.314,
    )
    assert conductivity == pytest.approx(0.0229, abs=0.00015)
    diffusivity = fuller_binary_diffusivity_m2_s(298.0, 1.0e5, 76.5, 29.0, 80.25, 19.7)
    assert diffusivity / 1e-4 == pytest.approx(0.096, abs=0.001)


def test_reduced_pressure_at_or_above_two_raises() -> None:
    with pytest.raises(PropertyRangeError):
        viscosity_pa_s("H2", 800.0, 2.0 * 13.0 * 1.0e5)
    viscosity_pa_s("H2", 800.0, 500_000.0)


def test_polar_and_nonpolar_conductivity_against_table_10_2() -> None:
    """表 10-2（PDF 569）：NH3 在 273 K 的实验值 0.0219，CO2 在 300 K 为 0.0167 W/(m·K)。"""
    ammonia = thermal_conductivity_w_m_k("NH3", 273.0, 1.0e5)
    carbon_dioxide = thermal_conductivity_w_m_k("CO2", 300.0, 1.0e5)
    ammonia_error = (ammonia - 0.0219) / 0.0219
    carbon_dioxide_error = (carbon_dioxide - 0.0167) / 0.0167
    print(
        f"NH3 273 K error {ammonia_error:.3%}; CO2 300 K error {carbon_dioxide_error:.3%}"
    )
    assert abs(ammonia_error) < 0.5
    assert abs(carbon_dioxide_error) < 0.2


def test_interphase_diffusivity_follows_equation_5_20() -> None:
    temperature, pressure = 1200.0, 500_000.0
    expected = 3.13e-4 * (temperature / 1500.0) ** 1.75 * (101_300.0 / pressure)
    assert interphase_diffusivity_m2_s(temperature, pressure) == pytest.approx(expected)


def test_adiabatic_gas_equilibrium_matches_feed_enthalpy() -> None:
    temperature = 1000.0
    pressure = 500_000.0
    feed = {"C": 1.0, "O": 2.0}
    amounts = {"CO2": 1.0}
    enthalpy = enthalpy_j_mol("CO2", temperature) * 1.0
    result = adiabatic_gas_equilibrium(enthalpy, pressure, feed, ["CO", "CO2", "O2"])
    product = sum(
        amount * enthalpy_j_mol(name, result.temperature_k)
        for name, amount in result.amounts_mol.items()
    )
    assert product == pytest.approx(enthalpy, rel=1e-4)
    assert result.temperature_k == pytest.approx(temperature, abs=15.0)


def test_equilibrium_failure_is_reported() -> None:
    with pytest.raises(EquilibriumError):
        gibbs_minimization(1000.0, 101_300.0, {"C": 1.0}, ["H2"])


def test_source_only_feed_keeps_the_given_amounts() -> None:
    result = gibbs_minimization(
        1000.0, 500_000.0, {"N": 2.0}, ["N2"], source_amounts_mol={"N2": 1.0},
    )
    assert result.amounts_mol["N2"] == pytest.approx(1.0)
    assert result.total_mol == pytest.approx(1.0)
    assert result.mole_fractions["N2"] == pytest.approx(1.0)


def test_sources_stay_beside_equilibrium_gas() -> None:
    result = gibbs_minimization(
        1000.0, 500_000.0, {"C": 1.0, "O": 2.0, "N": 2.0},
        ["CO", "CO2", "O2", "N2"], source_amounts_mol={"N2": 1.0},
    )
    assert result.amounts_mol["N2"] == pytest.approx(1.0)
    totals = {"C": 0.0, "O": 0.0, "N": 0.0}
    for name, amount in result.amounts_mol.items():
        for element, count in nasa_species(name)["atoms"].items():
            if element in totals:
                totals[element] += amount * count
    assert totals == pytest.approx({"C": 1.0, "O": 2.0, "N": 2.0}, abs=1e-8)
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-8)
