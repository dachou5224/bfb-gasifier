# -*- coding: utf-8 -*-
"""T22：干燥、热解产率、热解气与释放份额。"""
from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from hamel.cases.config import SolverConfig
from hamel.properties.data import atomic_weight, gas_constant_j_mol_k, molar_mass_kg_mol
from hamel.pyrolysis import (
    PyrolysisDomainError,
    PyrolysisError,
    ReleaseTimes,
    carbon_retained_in_char,
    carbon_retained_in_tar,
    lignite_particle_properties,
    lignite_yields,
    particle_release_times,
    pyrolysis_duration_s,
    pyrolysis_gas,
    release_distribution,
    remaining_time_integral,
)
from hamel.pyrolysis.particle import (
    conductivity_w_m_k,
    corrected_vaporization_j_kg,
    drying_energy_relative_gap,
    evaporation_temperature_k,
    heat_transfer_coefficient,
)

_CASE = json.loads(Path("data/cases/vtt_pressurized.json").read_text(encoding="utf-8"))
_FUEL = _CASE["fuels"]["rhein_braunkohle"]
_RUN = _CASE["cases"]["1"]
# 819 °C 是冒烟时指定的床温，不是模型结果。
_BED_K = 819.0 + 273.15
_PRESSURE = float(_RUN["pressure"]["value"])
_MOISTURE = float(_RUN["moisture_mass_fraction"]["value"])
_ASH = float(_FUEL["ash_mass_fraction_wf"]["value"])
_WF = {
    "C": float(_FUEL["carbon_mass_fraction_wf"]["value"]),
    "H": float(_FUEL["hydrogen_mass_fraction_wf"]["value"]),
    "O": float(_FUEL["oxygen_mass_fraction_wf"]["value"]),
    "N": float(_FUEL["nitrogen_mass_fraction_wf"]["value"]),
    "S": float(_FUEL["sulfur_mass_fraction_wf"]["value"]),
}
# ADR 0014：燃料加入温度。不是入口气体温度。
_FEED_K = 298.15
_GAS = {"N2": 0.79, "O2": 0.21}


def _yields():
    return lignite_yields(_PRESSURE, _BED_K, _WF["C"], _WF["H"], _WF["O"], _WF["N"], _WF["S"], _ASH)


def _hand_carbon_fractions():
    pressure, temperature = _PRESSURE, _BED_K
    char = 0.3027 * pressure**0.1047 * 0.9356 ** (1.0e-6 * pressure) * temperature ** (-0.08578)
    tar = -0.2168 * pressure**0.04148 * temperature**0.1566 + 1.333
    return char, tar


def test_each_yield_equation_matches_a_hand_evaluation():
    char_carbon, tar_carbon = _hand_carbon_fractions()
    assert carbon_retained_in_char(_PRESSURE, _BED_K) == pytest.approx(char_carbon, rel=1e-12)
    assert carbon_retained_in_tar(_PRESSURE, _BED_K) == pytest.approx(tar_carbon, rel=1e-12)
    carbon = _WF["C"] / (1.0 - _ASH)
    hydrogen = _WF["H"] / (1.0 - _ASH)
    yields = _yields()
    assert yields.x_c == pytest.approx(carbon)
    assert yields.x_h == pytest.approx(hydrogen)
    assert yields.x_char == pytest.approx(carbon * char_carbon, rel=1e-12)
    assert yields.x_tar == pytest.approx(carbon * tar_carbon * (1.0 + hydrogen / carbon), rel=1e-12)
    assert yields.x_c_gas == pytest.approx(carbon - yields.x_char - tar_carbon * carbon, rel=1e-12)
    assert yields.x_h_gas == pytest.approx(hydrogen - carbon * tar_carbon * (hydrogen / carbon), rel=1e-12)
    carbon_mol = yields.x_c * tar_carbon / (atomic_weight("C") * 1.0e-3)
    hydrogen_mol = (yields.x_tar - yields.x_c * tar_carbon) / (atomic_weight("H") * 1.0e-3)
    naphthalene = 0.5 * (carbon_mol - hydrogen_mol)
    benzene = (hydrogen_mol - 8.0 * naphthalene) / 6.0
    assert yields.benzene_mol_per_kg == pytest.approx(benzene, rel=1e-12)
    assert yields.naphthalene_mol_per_kg == pytest.approx(naphthalene, rel=1e-12)
    ratio = (6.0 * benzene + 8.0 * naphthalene) / (6.0 * benzene + 10.0 * naphthalene)
    assert ratio == pytest.approx(hydrogen_mol / carbon_mol, rel=1e-12)


def test_tar_outside_benzene_naphthalene_range_raises():
    with pytest.raises(PyrolysisError):
        lignite_yields(_PRESSURE, _BED_K, 0.20, 0.50, 0.20, 0.01, 0.01, 0.08)


def test_pyrolysis_gas_conserves_the_element_pool_without_drying_water():
    yields = _yields()
    equilibrium = pyrolysis_gas(yields, _BED_K, _PRESSURE)
    assert sum(equilibrium.mole_fractions.values()) == pytest.approx(1.0, abs=1e-8)
    expected = {
        "C": yields.x_c_gas / (atomic_weight("C") * 1.0e-3),
        "H": yields.x_h_gas / (atomic_weight("H") * 1.0e-3),
        "O": yields.x_o / (atomic_weight("O") * 1.0e-3),
        "N": yields.x_n / (atomic_weight("N") * 1.0e-3),
        "S": yields.x_s / (atomic_weight("S") * 1.0e-3),
    }
    from hamel.properties.data import atoms

    got = {name: 0.0 for name in expected}
    for species, amount in equilibrium.amounts_mol.items():
        for element, count in atoms(species).items():
            got[element] += amount * count
    for name, target in expected.items():
        assert got[name] == pytest.approx(target, rel=1e-8, abs=1e-10)


def test_exponential_decay_integral_plus_tail_is_the_time_constant():
    tau = 2.5
    grid = np.linspace(0.0, 40.0, 8001)
    series = np.exp(-grid / tau)
    assert remaining_time_integral(grid, series, tail=float(series[-1]) * tau) == pytest.approx(tau, rel=1e-6)


def test_isothermal_daem_matches_an_independent_quadrature():
    material = lignite_particle_properties()
    gas_constant = gas_constant_j_mol_k()
    mean, sigma, frequency = material.activation_mean_j_mol, material.activation_sigma_j_mol, material.frequency_s

    def integrand(energy: float) -> float:
        density = math.exp(-0.5 * ((energy - mean) / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))
        return density / (frequency * math.exp(-energy / (gas_constant * _BED_K)))

    reference, _estimated = quad(integrand, mean - 12.0 * sigma, mean + 12.0 * sigma, epsabs=1.0e-12, limit=400)
    closed = math.exp(mean / (gas_constant * _BED_K) + sigma**2 / (2.0 * gas_constant**2 * _BED_K**2)) / frequency
    assert reference == pytest.approx(closed, rel=1e-8)
    times = particle_release_times(1.0e-3, _BED_K, _PRESSURE, 0.0, _BED_K, {"N2": 1.0}, 0.0, 1.0, material, SolverConfig())
    assert times.tau_dry_s == 0.0
    assert times.t_v_s == pytest.approx(reference, rel=1.0e-3)
    assert times.tail_error_bound_s < 1.0e-3 * times.t_v_s


def test_volatile_time_shorter_than_water_time_raises():
    with pytest.raises(PyrolysisDomainError):
        pyrolysis_duration_s(2.0, 1.5)


def test_zero_moisture_and_zero_volatiles_do_not_divide():
    material, config = lignite_particle_properties(), SolverConfig()
    dry = particle_release_times(1.0e-4, _BED_K, _PRESSURE, 0.0, _BED_K, {"N2": 1.0}, 0.0, 1.0, material, config)
    assert dry.tau_dry_s == 0.0 and dry.t_v_s > 0.0
    wet = particle_release_times(1.0e-4, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, {"N2": 1.0}, 0.0, 0.0, material, config)
    assert wet.t_v_s == 0.0 and wet.tau_pyr_s == 0.0 and wet.tau_dry_s > 0.0
    assert particle_release_times(1.0e-4, _BED_K, _PRESSURE, 0.0, _BED_K, {"N2": 1.0}, 0.0, 0.0, material, config) == ReleaseTimes(0.0, 0.0, 0.0, 0.0, 0.0)


def test_drying_energy_gap_shrinks_when_the_radial_grid_is_refined():
    material, diameter = lignite_particle_properties(), 3.466e-4
    coarse = drying_energy_relative_gap(diameter, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, _GAS, 0.5, material, SolverConfig())
    fine = drying_energy_relative_gap(
        diameter, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, _GAS, 0.5, material,
        SolverConfig(pyrolysis_radial_intervals=48),
    )
    assert coarse < 5.0e-3
    assert fine < coarse


def test_small_biot_drying_time_matches_the_lumped_shrinking_core():
    material = replace(lignite_particle_properties(), diffusivity_m2_s=1.0e-4)
    diameter = 3.466e-4
    evaporation = evaporation_temperature_k(_PRESSURE)
    corrected, dry_basis = corrected_vaporization_j_kg(_MOISTURE, evaporation, _FEED_K, material.heat_capacity_j_kg_k)
    alpha, biot = heat_transfer_coefficient(
        diameter, _BED_K, evaporation, _PRESSURE, _GAS, 0.0, conductivity_w_m_k(material),
    )
    assert biot < 1.0e-3
    duration = 0.5 * diameter * material.density_kg_m3 * dry_basis * corrected / (3.0 * alpha * (_BED_K - evaporation))
    times = particle_release_times(diameter, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, _GAS, 0.0, 1.0, material, SolverConfig())
    assert times.t_w_s == pytest.approx(0.5 * duration, rel=2.0e-3)


def _relative(left: float, right: float) -> float:
    return abs(left - right) / right


def test_time_radial_and_activation_refinement_converges_release_times():
    material = lignite_particle_properties()
    diameter = float(_RUN["fuel_size_class_diameters"]["value"][-1])
    base = particle_release_times(diameter, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, _GAS, 0.5, 1.0, material, SolverConfig())
    variants = (
        SolverConfig(pyrolysis_radial_intervals=48),
        SolverConfig(pyrolysis_time_rtol=1.0e-7),
        SolverConfig(pyrolysis_activation_nodes=24),
    )
    for config in variants:
        other = particle_release_times(diameter, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, _GAS, 0.5, 1.0, material, config)
        assert _relative(base.t_w_s, other.t_w_s) < 1.0e-3
        assert _relative(base.t_v_s, other.t_v_s) < 1.0e-3
        assert _relative(base.tau_pyr_s, other.tau_pyr_s) < 1.0e-3


def test_single_cell_release_matches_the_feed_and_closes():
    yields = _yields()
    feed = np.array([0.01, 0.02])
    zeros = np.zeros((1, 2))
    result = release_distribution(
        feed, 0.05, [0.0], [0.10], zeros, zeros, zeros, np.array([2.0, 3.0]), np.array([4.0, 5.0]),
        _MOISTURE, _ASH, yields, _BED_K, _PRESSURE,
    )
    waf = (1.0 - _MOISTURE) * (1.0 - _ASH)
    assert result.water_kg_s[0] == pytest.approx(_MOISTURE * float(np.sum(feed)))
    assert result.char_kg_s[0] == pytest.approx(waf * yields.x_char * feed)
    assert result.ash_kg_s[0] == pytest.approx((1.0 - _MOISTURE) * _ASH * feed)
    assert result.undried_kg[0] == pytest.approx(feed * np.array([2.0, 3.0]))
    assert np.max(np.abs(result.undried_exit_kg_s)) == 0.0
    assert np.max(np.abs(result.mass_gap_kg_s)) < 1.0e-12
    for gap in result.element_gap_kg_s.values():
        assert np.max(np.abs(gap)) < 1.0e-12
    instant = release_distribution(
        feed, 0.05, [0.0], [0.10], zeros, zeros, zeros, np.zeros(2), np.zeros(2),
        _MOISTURE, _ASH, yields, _BED_K, _PRESSURE,
    )
    assert np.all(instant.undried_kg == 0.0)
    assert np.all(instant.dried_kg == 0.0)
    assert instant.water_kg_s[0] == pytest.approx(_MOISTURE * float(np.sum(feed)))


def test_feed_on_a_shared_face_enters_the_upper_cell():
    from hamel.pyrolysis import feed_cell_index

    bottom = [0.30, 0.15, 0.0]
    top = [0.45, 0.30, 0.15]
    assert feed_cell_index(0.15, bottom, top) == 1
    assert feed_cell_index(0.12, bottom, top) == 2


def _inlet_velocity_and_gas():
    """一次风、蒸汽、吹扫氮在指定床温下的表观速度。二次风不计入。这不是流动模型的输出。"""
    air = float(_RUN["primary_air_amount_flow"]["value"])
    purge = float(_RUN["nitrogen_purge_amount_flow"]["value"])
    steam = float(_RUN["steam_mass_flow"]["value"]) / molar_mass_kg_mol("H2O")
    oxygen, nitrogen, water = 0.21 * air, 0.79 * air + purge, steam
    total = oxygen + nitrogen + water
    area = math.pi * (0.15 / 2.0) ** 2
    velocity = total * gas_constant_j_mol_k() * _BED_K / (_PRESSURE * area)
    return velocity, {"O2": oxygen / total, "N2": nitrogen / total, "H2O": water / total}


def test_case1_size_classes_and_three_cell_release_close():
    material = lignite_particle_properties()
    yields = _yields()
    velocity, gas = _inlet_velocity_and_gas()
    diameters = [float(value) for value in _RUN["fuel_size_class_diameters"]["value"]]
    fractions = np.array(_RUN["fuel_size_class_mass_fractions"]["value"], dtype=float)
    feed = float(_RUN["fuel_mass_flow"]["value"]) * fractions
    dry = np.empty(len(diameters))
    pyro = np.empty(len(diameters))
    for index, diameter in enumerate(diameters):
        times = particle_release_times(
            diameter, _BED_K, _PRESSURE, _MOISTURE, _FEED_K, gas, velocity, 1.0, material, SolverConfig(),
        )
        dry[index], pyro[index] = times.tau_dry_s, times.tau_pyr_s
        assert times.t_v_s >= times.t_w_s
    bottom = np.array([0.30, 0.15, 0.0])
    top = np.array([0.45, 0.30, 0.15])
    upward = np.full((3, len(diameters)), 0.20)
    downward = np.full((3, len(diameters)), 0.05)
    leaving = np.zeros((3, len(diameters)))
    leaving[0, :] = 0.01
    # 三 cell 示例不是装置的 16 格。0.12 m 是工况的给料高度，落在本示例最下层。
    result = release_distribution(
        feed, 0.12, bottom, top, upward, downward, leaving, dry, pyro,
        _MOISTURE, _ASH, yields, _BED_K, _PRESSURE,
    )
    assert result.water_kg_s.shape == (3,)
    assert np.max(np.abs(result.mass_gap_kg_s)) < 1.0e-9
    for gap in result.element_gap_kg_s.values():
        assert np.max(np.abs(gap)) < 1.0e-9
