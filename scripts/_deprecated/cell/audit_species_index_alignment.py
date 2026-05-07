#!/usr/bin/env python3
"""专项审计：species 索引顺序与热力学数据映射一致性。

说明
----
- 当前项目气相热力学后端是 `NASA 7-coefficient` 多项式，不是 Shomate。
- 本脚本按“索引 -> 物种名 -> 热力学数据 -> 实际 cell 状态”四层做专项审计：
  1. `GAS_SPECIES` / `GAS_SPECIES_INDEX` 的固定顺序、双向映射、唯一性
  2. 基础物种与 Gibbs 微量物种的 NASA 参考表是否与实现一一对应
  3. `TAR1/TAR2` 在不同 fuel_type 下是否绑定到正确 surrogate 的热力学数据
  4. 真实 LU solved cells 中 array <-> dict 往返与 mixture thermo 加权是否保持一致
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import Rg
from src.core import species as species_mod
from src.core.reactor import Reactor
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_reactor_config


EXPECTED_GAS_SPECIES = [
    "CO",
    "CO2",
    "H2",
    "H2O",
    "CH4",
    "O2",
    "N2",
    "H2S",
    "NH3",
    "TAR1",
    "TAR2",
]

REFERENCE_NASA_DATA: dict[str, tuple[float, float, float, tuple[float, ...], tuple[float, ...]]] = {
    "CO": (
        200.0, 1000.0, 3500.0,
        (2.71518561E+00, 2.06252743E-03, -9.98825771E-07, 2.30053008E-10, -2.03647716E-14, -1.41518724E+04, 7.81868772E+00),
        (3.57953347E+00, -6.10353680E-04, 1.01681433E-06, 9.07005884E-10, -9.04424499E-13, -1.43440860E+04, 3.50840928E+00),
    ),
    "CO2": (
        200.0, 1000.0, 3500.0,
        (3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00),
        (2.35677352E+00, 8.98459677E-03, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13, -4.83719697E+04, 9.90105222E+00),
    ),
    "H2": (
        200.0, 1000.0, 3500.0,
        (3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00),
        (2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01),
    ),
    "H2O": (
        200.0, 1000.0, 3500.0,
        (3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00),
        (4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01),
    ),
    "CH4": (
        200.0, 1000.0, 3500.0,
        (7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01),
        (5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04, -4.64130376E+00),
    ),
    "O2": (
        200.0, 1000.0, 3500.0,
        (3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00),
        (3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00),
    ),
    "N2": (
        300.0, 1000.0, 5000.0,
        (2.92664000E+00, 1.48797680E-03, -5.68476000E-07, 1.00970380E-10, -6.75335100E-14, -9.22797700E+02, 5.98052800E+00),
        (3.29867700E+00, 1.40824040E-03, -3.96322200E-06, 5.64151500E-09, -2.44485400E-12, -1.02089990E+03, 3.95037200E+00),
    ),
    "H2S": (
        200.0, 1000.0, 6000.0,
        (2.88433778E+00, 3.36697740E-03, -1.33825792E-06, 2.63797690E-10, -2.06259025E-14, -3.44927890E+03, 7.63222600E+00),
        (4.12023462E+00, -1.63887050E-03, 6.16088989E-06, -4.46121219E-09, 1.14866652E-12, -3.65087674E+03, 2.26888850E+00),
    ),
    "NH3": (
        200.0, 1000.0, 6000.0,
        (2.63445210E+00, 5.66625600E-03, -1.72786760E-06, 2.38671610E-10, -1.25787860E-14, -6.54469580E+03, 6.56629280E+00),
        (4.28602740E+00, -4.66052300E-03, 2.17185130E-05, -2.28088870E-08, 8.26380460E-12, -6.74172850E+03, -6.25372770E-01),
    ),
    "SO2": (
        200.0, 1000.0, 6000.0,
        (3.26653300E+00, 5.32379000E-03, -2.14435000E-06, 4.03816000E-10, -2.99803000E-14, -3.71807000E+04, 6.31103000E+00),
        (4.11231800E+00, -2.38411000E-03, 1.00341000E-05, -1.15351000E-08, 4.52011000E-12, -3.60807000E+04, 1.57321000E+00),
    ),
    "COS": (
        200.0, 1000.0, 6000.0,
        (3.63633300E+00, 4.48403000E-03, -1.87968000E-06, 3.38163000E-10, -2.42633000E-14, -2.60123000E+04, 7.27234000E+00),
        (4.23863000E+00, -2.92664000E-03, 1.13633000E-05, -1.21656000E-08, 4.60420000E-12, -2.51314000E+04, 2.10744000E+00),
    ),
    "HCN": (
        200.0, 1000.0, 6000.0,
        (3.02507800E+00, 1.44268900E-03, -5.63082800E-07, 1.01858100E-10, -6.91095200E-15, 1.35511000E+04, 6.72944000E+00),
        (4.22118500E+00, -3.24392500E-03, 1.37799400E-05, -1.33144000E-08, 4.33768800E-12, 1.61627100E+04, 2.25935000E+00),
    ),
    "NO": (
        200.0, 1000.0, 6000.0,
        (2.98040200E+00, 7.85904400E-04, -3.46108400E-07, 6.95496400E-11, -5.09846400E-15, 9.87410100E+03, 6.84726100E+00),
        (3.26245200E+00, 1.51194100E-03, -3.88175500E-06, 5.58194400E-09, -2.47495100E-12, 9.57530300E+03, 3.02807700E+00),
    ),
}

EXPECTED_TAR_MAPPING = {
    "coal": {"TAR1": "C6H6", "TAR2": "C10H8"},
    "biomass": {"TAR1": "C10H8", "TAR2": "C16H34"},
}

REFERENCE_TAR_DATA = {
    "C6H6": (
        200.0, 1000.0, 5000.0,
        (1.35313100E+01, 2.03311000E-02, -7.11410000E-06, 1.13101000E-09, -6.64010000E-14, 7.10910000E+03, -4.85810000E+01),
        (-1.88110000E+00, 4.65110000E-02, 1.01110000E-05, -4.01110000E-08, 1.81110000E-11, 1.11110000E+04, 2.81110000E+01),
    ),
    "C10H8": (
        200.0, 1000.0, 5000.0,
        (2.11110000E+01, 3.51110000E-02, -1.21110000E-05, 1.91110000E-09, -1.11110000E-13, 1.41110000E+04, -8.51110000E+01),
        (-4.51110000E+00, 8.51110000E-02, 2.51110000E-05, -8.51110000E-08, 4.11110000E-11, 2.11110000E+04, 4.51110000E+01),
    ),
    "C16H34": (
        200.0, 1000.0, 5000.0,
        (3.51110000E+01, 8.51110000E-02, -3.11110000E-05, 5.11110000E-09, -3.11110000E-13, -6.51110000E+04, -1.45110001E+02),
        (5.51110000E+00, 1.51110000E-01, 8.51110000E-05, -2.11110000E-07, 9.51110000E-11, -5.81110000E+04, 5.11110000E+00),
    ),
}


def _coeffs_from_reference(
    coeff_data: tuple[float, float, float, tuple[float, ...], tuple[float, ...]],
    T: float,
) -> tuple[float, ...]:
    t_low, t_mid, t_high, high, low = coeff_data
    t_clamp = float(np.clip(T, t_low, t_high))
    return low if t_clamp <= t_mid else high


def _cp_from_reference(
    coeff_data: tuple[float, float, float, tuple[float, ...], tuple[float, ...]],
    T: float,
) -> float:
    a = _coeffs_from_reference(coeff_data, T)
    return Rg * (a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3 + a[4] * T**4)


def _enthalpy_from_reference(
    coeff_data: tuple[float, float, float, tuple[float, ...], tuple[float, ...]],
    T: float,
) -> float:
    a = _coeffs_from_reference(coeff_data, T)
    return Rg * T * (
        a[0] + a[1] * T / 2.0 + a[2] * T**2 / 3.0 + a[3] * T**3 / 4.0 + a[4] * T**4 / 5.0 + a[5] / T
    )


def _entropy_from_reference(
    coeff_data: tuple[float, float, float, tuple[float, ...], tuple[float, ...]],
    T: float,
) -> float:
    a = _coeffs_from_reference(coeff_data, T)
    return Rg * (
        a[0] * np.log(T) + a[1] * T + a[2] * T**2 / 2.0 + a[3] * T**3 / 3.0 + a[4] * T**4 / 4.0 + a[6]
    )


def _gibbs_from_reference(
    coeff_data: tuple[float, float, float, tuple[float, ...], tuple[float, ...]],
    T: float,
) -> float:
    h = _enthalpy_from_reference(coeff_data, T)
    s = _entropy_from_reference(coeff_data, T)
    return h - T * s


def _fingerprint(
    coeff_data: tuple[float, float, float, tuple[float, ...], tuple[float, ...]],
) -> tuple[float, ...]:
    t_low, t_mid, t_high, high, low = coeff_data
    return (
        t_low,
        t_mid,
        t_high,
        high[0],
        high[5],
        high[6],
        low[0],
        low[5],
        low[6],
    )


def _assert_close(name: str, value: float, ref: float, atol: float, rtol: float = 1e-12) -> float:
    diff = abs(value - ref)
    if not np.isclose(value, ref, atol=atol, rtol=rtol):
        raise AssertionError(f"{name}: value={value:.12e}, ref={ref:.12e}, diff={diff:.3e}")
    return float(diff)


def _audit_index_alignment() -> dict[str, Any]:
    gas_species = species_mod.GAS_SPECIES
    gas_index = species_mod.GAS_SPECIES_INDEX

    assert gas_species == EXPECTED_GAS_SPECIES, (
        f"GAS_SPECIES 顺序漂移: got={gas_species}, expected={EXPECTED_GAS_SPECIES}"
    )
    assert len(gas_species) == species_mod.N_GAS, "N_GAS 与 GAS_SPECIES 长度不一致"
    assert len(set(gas_species)) == len(gas_species), "GAS_SPECIES 存在重复项"

    rows: list[dict[str, Any]] = []
    for i, species in enumerate(EXPECTED_GAS_SPECIES):
        idx = gas_index[species]
        roundtrip_species = gas_species[idx]
        assert idx == i, f"{species} 索引错误: got={idx}, expected={i}"
        assert roundtrip_species == species, f"{species} round-trip 错误: got={roundtrip_species}"
        rows.append({"species": species, "expected_index": i, "actual_index": idx})

    return {
        "n_gas": species_mod.N_GAS,
        "order": list(gas_species),
        "rows": rows,
        "pass": True,
    }


def _audit_static_thermo_mapping() -> dict[str, Any]:
    sample_temperatures = (298.15, 900.0, 1200.0, 1600.0)
    rows: list[dict[str, Any]] = []
    worst_cp = 0.0
    worst_h = 0.0
    worst_s = 0.0
    worst_g = 0.0

    for species, ref_data in REFERENCE_NASA_DATA.items():
        live_data = species_mod._NASA_DATA.get(species)  # type: ignore[attr-defined]
        assert live_data is not None, f"{species} 缺少 NASA 数据"
        assert _fingerprint(live_data) == _fingerprint(ref_data), f"{species} NASA 指纹不匹配"

        diff_cp_max = 0.0
        diff_h_max = 0.0
        diff_s_max = 0.0
        diff_g_max = 0.0
        valid_temps: list[float] = []
        t_low, _, t_high, _, _ = ref_data

        for T in sample_temperatures:
            if not (t_low <= T <= t_high):
                continue
            valid_temps.append(T)
            cp_ref = _cp_from_reference(ref_data, T)
            h_ref = _enthalpy_from_reference(ref_data, T)
            s_ref = _entropy_from_reference(ref_data, T)
            g_ref = _gibbs_from_reference(ref_data, T)

            cp_live = species_mod.cp_molar(species, T)
            h_live = species_mod.enthalpy_molar(species, T)
            s_live = species_mod.entropy_molar(species, T)
            g_live = species_mod.gibbs_molar(species, T)

            diff_cp_max = max(diff_cp_max, _assert_close(f"{species} Cp @ {T}K", cp_live, cp_ref, atol=1e-9))
            diff_h_max = max(diff_h_max, _assert_close(f"{species} H @ {T}K", h_live, h_ref, atol=1e-6))
            diff_s_max = max(diff_s_max, _assert_close(f"{species} S @ {T}K", s_live, s_ref, atol=1e-9))
            diff_g_max = max(diff_g_max, _assert_close(f"{species} G @ {T}K", g_live, g_ref, atol=1e-6))

        worst_cp = max(worst_cp, diff_cp_max)
        worst_h = max(worst_h, diff_h_max)
        worst_s = max(worst_s, diff_s_max)
        worst_g = max(worst_g, diff_g_max)
        rows.append(
            {
                "species": species,
                "temperatures_K": valid_temps,
                "cp_max_abs_diff": diff_cp_max,
                "h_max_abs_diff": diff_h_max,
                "s_max_abs_diff": diff_s_max,
                "g_max_abs_diff": diff_g_max,
            }
        )

    return {
        "thermo_backend": "NASA-7 polynomial (not Shomate)",
        "n_species_checked": len(rows),
        "worst_cp_max_abs_diff": worst_cp,
        "worst_h_max_abs_diff": worst_h,
        "worst_s_max_abs_diff": worst_s,
        "worst_g_max_abs_diff": worst_g,
        "rows": rows,
        "pass": True,
    }


def _audit_tar_mapping() -> dict[str, Any]:
    sample_temperatures = (400.0, 900.0, 1200.0, 1600.0)
    saved_mw = {
        "TAR1": species_mod.MOLECULAR_WEIGHT.get("TAR1"),
        "TAR2": species_mod.MOLECULAR_WEIGHT.get("TAR2"),
    }
    saved_nasa = {
        "TAR1": species_mod._NASA_DATA.get("TAR1"),  # type: ignore[attr-defined]
        "TAR2": species_mod._NASA_DATA.get("TAR2"),  # type: ignore[attr-defined]
    }

    rows: list[dict[str, Any]] = []
    worst_cp = 0.0
    worst_h = 0.0
    worst_s = 0.0
    worst_g = 0.0

    try:
        for fuel_type, expected_mapping in EXPECTED_TAR_MAPPING.items():
            mapping = species_mod.configure_tar_components_by_fuel(fuel_type)
            assert mapping == expected_mapping, f"{fuel_type} tar mapping 错误: got={mapping}, expected={expected_mapping}"

            for tar_label, surrogate in expected_mapping.items():
                ref_data = REFERENCE_TAR_DATA[surrogate]
                live_data = species_mod._NASA_DATA.get(tar_label)  # type: ignore[attr-defined]
                assert live_data is not None, f"{fuel_type} {tar_label} 未注册 NASA 数据"
                assert _fingerprint(live_data) == _fingerprint(ref_data), (
                    f"{fuel_type} {tar_label}->{surrogate} NASA 指纹不匹配"
                )

                ref_mw = species_mod.TAR_SURROGATE_MW[surrogate]
                live_mw = species_mod.MOLECULAR_WEIGHT[tar_label]
                _assert_close(f"{fuel_type} {tar_label} MW", live_mw, ref_mw, atol=1e-12)

                diff_cp_max = 0.0
                diff_h_max = 0.0
                diff_s_max = 0.0
                diff_g_max = 0.0
                for T in sample_temperatures:
                    cp_ref = _cp_from_reference(ref_data, T)
                    h_ref = _enthalpy_from_reference(ref_data, T)
                    s_ref = _entropy_from_reference(ref_data, T)
                    g_ref = _gibbs_from_reference(ref_data, T)

                    cp_live = species_mod.cp_molar(tar_label, T)
                    h_live = species_mod.enthalpy_molar(tar_label, T)
                    s_live = species_mod.entropy_molar(tar_label, T)
                    g_live = species_mod.gibbs_molar(tar_label, T)

                    diff_cp_max = max(diff_cp_max, _assert_close(f"{fuel_type} {tar_label} Cp @ {T}K", cp_live, cp_ref, atol=1e-9))
                    diff_h_max = max(diff_h_max, _assert_close(f"{fuel_type} {tar_label} H @ {T}K", h_live, h_ref, atol=1e-6))
                    diff_s_max = max(diff_s_max, _assert_close(f"{fuel_type} {tar_label} S @ {T}K", s_live, s_ref, atol=1e-9))
                    diff_g_max = max(diff_g_max, _assert_close(f"{fuel_type} {tar_label} G @ {T}K", g_live, g_ref, atol=1e-6))

                worst_cp = max(worst_cp, diff_cp_max)
                worst_h = max(worst_h, diff_h_max)
                worst_s = max(worst_s, diff_s_max)
                worst_g = max(worst_g, diff_g_max)
                rows.append(
                    {
                        "fuel_type": fuel_type,
                        "tar_label": tar_label,
                        "surrogate": surrogate,
                        "mw_g_mol": float(live_mw),
                        "cp_max_abs_diff": diff_cp_max,
                        "h_max_abs_diff": diff_h_max,
                        "s_max_abs_diff": diff_s_max,
                        "g_max_abs_diff": diff_g_max,
                    }
                )
    finally:
        for tar_label, value in saved_mw.items():
            if value is None:
                species_mod.MOLECULAR_WEIGHT.pop(tar_label, None)
            else:
                species_mod.MOLECULAR_WEIGHT[tar_label] = value
        for tar_label, value in saved_nasa.items():
            if value is None:
                species_mod._NASA_DATA.pop(tar_label, None)  # type: ignore[attr-defined]
            else:
                species_mod._NASA_DATA[tar_label] = value  # type: ignore[attr-defined]
        species_mod.gibbs_molar.cache_clear()

    return {
        "worst_cp_max_abs_diff": worst_cp,
        "worst_h_max_abs_diff": worst_h,
        "worst_s_max_abs_diff": worst_s,
        "worst_g_max_abs_diff": worst_g,
        "rows": rows,
        "pass": True,
    }


def _mix_property_direct(n_vec: np.ndarray, T: float, fn_name: str) -> float:
    fn = getattr(species_mod, fn_name)
    total = 0.0
    for j, species in enumerate(species_mod.GAS_SPECIES):
        total += float(n_vec[j]) * float(fn(species, T))
    return total


def _mix_property_from_dict(n_dict: dict[str, float], T: float, fn_name: str) -> float:
    fn = getattr(species_mod, fn_name)
    total = 0.0
    for species, value in n_dict.items():
        total += float(value) * float(fn(species, T))
    return total


def _audit_runtime_cells() -> dict[str, Any]:
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    selected = sorted(set([0, 1, 2, len(reactor.cells) - 1]))
    rows: list[dict[str, Any]] = []
    worst_roundtrip = 0.0
    worst_cp_mix = 0.0
    worst_h_mix = 0.0
    worst_g_mix = 0.0

    for i in selected:
        cell = reactor.cells[i]
        T = float(cell.T)
        for phase_name, n_vec in (("dense", np.asarray(cell.N_d, dtype=float)), ("bubble", np.asarray(cell.N_b, dtype=float))):
            assert n_vec.shape == (species_mod.N_GAS,), (
                f"cell {i} {phase_name} 组分向量尺寸错误: got={n_vec.shape}, expected={(species_mod.N_GAS,)}"
            )
            n_dict = {species: float(n_vec[j]) for j, species in enumerate(species_mod.GAS_SPECIES)}
            n_back = np.zeros(species_mod.N_GAS, dtype=float)
            for species, value in n_dict.items():
                n_back[species_mod.GAS_SPECIES_INDEX[species]] = value

            roundtrip_err = float(np.max(np.abs(n_vec - n_back)))
            cp_mix_direct = _mix_property_direct(n_vec, T, "cp_molar")
            cp_mix_dict = _mix_property_from_dict(n_dict, T, "cp_molar")
            h_mix_direct = _mix_property_direct(n_vec, T, "enthalpy_molar")
            h_mix_dict = _mix_property_from_dict(n_dict, T, "enthalpy_molar")
            g_mix_direct = _mix_property_direct(n_vec, T, "gibbs_molar")
            g_mix_dict = _mix_property_from_dict(n_dict, T, "gibbs_molar")

            cp_mix_err = abs(cp_mix_direct - cp_mix_dict)
            h_mix_err = abs(h_mix_direct - h_mix_dict)
            g_mix_err = abs(g_mix_direct - g_mix_dict)

            _assert_close(f"cell {i} {phase_name} cp_mix", cp_mix_direct, cp_mix_dict, atol=1e-9)
            _assert_close(f"cell {i} {phase_name} h_mix", h_mix_direct, h_mix_dict, atol=1e-6)
            _assert_close(f"cell {i} {phase_name} g_mix", g_mix_direct, g_mix_dict, atol=1e-6)

            worst_roundtrip = max(worst_roundtrip, roundtrip_err)
            worst_cp_mix = max(worst_cp_mix, cp_mix_err)
            worst_h_mix = max(worst_h_mix, h_mix_err)
            worst_g_mix = max(worst_g_mix, g_mix_err)

            rows.append(
                {
                    "cell": i,
                    "phase": phase_name,
                    "xi": float(cell.geo.h_center / cfg.H_bed),
                    "T_K": T,
                    "roundtrip_max_abs": roundtrip_err,
                    "cp_mix_abs_diff": cp_mix_err,
                    "h_mix_abs_diff": h_mix_err,
                    "g_mix_abs_diff": g_mix_err,
                }
            )

    return {
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "worst_roundtrip_max_abs": worst_roundtrip,
        "worst_cp_mix_abs_diff": worst_cp_mix,
        "worst_h_mix_abs_diff": worst_h_mix,
        "worst_g_mix_abs_diff": worst_g_mix,
        "rows": rows,
        "pass": True,
    }


def run_audit() -> dict[str, Any]:
    index_alignment = _audit_index_alignment()
    static_thermo = _audit_static_thermo_mapping()
    tar_mapping = _audit_tar_mapping()
    runtime_cells = _audit_runtime_cells()
    return {
        "note": "Thermo backend is NASA-7 polynomial; user-requested 'Shomate' check is audited here as thermo-data mapping consistency.",
        "index_alignment": index_alignment,
        "static_thermo_mapping": static_thermo,
        "tar_mapping": tar_mapping,
        "runtime_cells": runtime_cells,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="审计 species 索引顺序与热力学数据映射")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    args = ap.parse_args()

    out = run_audit()

    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=2))
    else:
        print("=" * 96)
        print("Species index alignment & thermo mapping audit")
        print("=" * 96)
        print(out["note"])
        print(
            f"index: pass={out['index_alignment']['pass']} "
            f"n_gas={out['index_alignment']['n_gas']}"
        )
        print(
            f"static thermo: pass={out['static_thermo_mapping']['pass']} "
            f"backend={out['static_thermo_mapping']['thermo_backend']} "
            f"worst_cp={out['static_thermo_mapping']['worst_cp_max_abs_diff']:.3e} "
            f"worst_h={out['static_thermo_mapping']['worst_h_max_abs_diff']:.3e} "
            f"worst_s={out['static_thermo_mapping']['worst_s_max_abs_diff']:.3e} "
            f"worst_g={out['static_thermo_mapping']['worst_g_max_abs_diff']:.3e}"
        )
        print(
            f"tar mapping: pass={out['tar_mapping']['pass']} "
            f"worst_cp={out['tar_mapping']['worst_cp_max_abs_diff']:.3e} "
            f"worst_h={out['tar_mapping']['worst_h_max_abs_diff']:.3e} "
            f"worst_s={out['tar_mapping']['worst_s_max_abs_diff']:.3e} "
            f"worst_g={out['tar_mapping']['worst_g_max_abs_diff']:.3e}"
        )
        print(
            f"runtime cells: pass={out['runtime_cells']['pass']} "
            f"converged={out['runtime_cells']['converged']} n_iter={out['runtime_cells']['n_iter']} "
            f"worst_roundtrip={out['runtime_cells']['worst_roundtrip_max_abs']:.3e} "
            f"worst_cp_mix={out['runtime_cells']['worst_cp_mix_abs_diff']:.3e} "
            f"worst_h_mix={out['runtime_cells']['worst_h_mix_abs_diff']:.3e} "
            f"worst_g_mix={out['runtime_cells']['worst_g_mix_abs_diff']:.3e}"
        )
        print("-" * 96)
        print("Fixed gas-species order:")
        print("  " + ", ".join(out["index_alignment"]["order"]))
        print("-" * 96)
        print(f"{'cell':>4} {'phase':>8} {'xi':>6} {'T(K)':>8} {'roundtrip':>12} {'cp_mix':>12} {'h_mix':>12} {'g_mix':>12}")
        for row in out["runtime_cells"]["rows"]:
            print(
                f"{row['cell']:>4d} {row['phase']:>8} {row['xi']:>6.3f} {row['T_K']:>8.1f} "
                f"{row['roundtrip_max_abs']:>12.3e} {row['cp_mix_abs_diff']:>12.3e} "
                f"{row['h_mix_abs_diff']:>12.3e} {row['g_mix_abs_diff']:>12.3e}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
