from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_pyrolysis import allocate_pyrolysis_products_elemental, calc_drying_pyrolysis_sources
from src.core.species import TAR_SURROGATE_FORMULA, calc_tar_surrogate_fractions, get_tar_component_mapping


def test_allocate_pyrolysis_products_elemental_preserves_element_budgets() -> None:
    prod = allocate_pyrolysis_products_elemental(
        nC=10.0,
        nH=12.0,
        nO=4.0,
        nN=0.2,
        nS=0.1,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=0.2,
    )

    tar_map = get_tar_component_mapping("coal")
    nC_out = prod["CO"] + prod["CO2"] + prod["CH4"]
    nH_out = 2.0 * prod["H2"] + 2.0 * prod["H2O"] + 4.0 * prod["CH4"] + 3.0 * prod["NH3"] + 2.0 * prod["H2S"]
    for tar_label in ("TAR1", "TAR2"):
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map[tar_label]]
        nC_out += c_tar * prod[tar_label]
        nH_out += h_tar * prod[tar_label]
    nO_out = prod["CO"] + 2.0 * prod["CO2"] + prod["H2O"]
    assert nC_out <= 10.0 + 1e-9
    assert nH_out <= 12.0 + 1e-9
    assert nO_out <= 4.0 + 1e-9
    assert prod["NH3"] == 0.2
    assert prod["H2S"] == 0.1

    tar_frac = calc_tar_surrogate_fractions("coal")
    assert prod["TAR1"] > 0.0
    assert prod["TAR2"] > 0.0
    ratio = prod["TAR1"] / max(prod["TAR1"] + prod["TAR2"], 1e-12)
    assert np.isclose(ratio, tar_frac[tar_map["TAR1"]], atol=1e-9)


def test_calc_drying_pyrolysis_sources_returns_consistent_signs() -> None:
    bundle = calc_drying_pyrolysis_sources(
        tau=5.0,
        T=1173.15,
        P=2.5e6,
        d_p=1.0e-3,
        moisture_wt=16.9,
        ash_dry_wt=11.41,
        C_dry=61.5,
        H_dry=4.1,
        O_dry=21.8,
        nitrogen_fraction=0.6,
        sulfur_fraction=0.01,
        sulfur_volatile_frac=0.5,
        pyrolysis_tar_carbon_frac=0.2,
        fuel_type="coal",
        m_vm_in=0.05,
        m_moist_in=0.02,
        solid_shape=(1, 4),
        vm_index=1,
        moisture_index=2,
    )

    assert bundle.gas_source.shape == (11,)
    assert bundle.solid_sink.shape == (1, 4)
    assert bundle.gas_source.sum() > 0.0
    assert bundle.solid_sink[0, 1] <= 0.0
    assert bundle.solid_sink[0, 2] <= 0.0
    assert 0.0 <= bundle.x_dry <= 1.0
    assert 0.0 <= bundle.x_vm <= 1.0


def test_calc_drying_pyrolysis_sources_preserves_total_sink_across_size_classes() -> None:
    bundle = calc_drying_pyrolysis_sources(
        tau=5.0,
        T=1173.15,
        P=2.5e6,
        d_p=1.0e-3,
        moisture_wt=16.9,
        ash_dry_wt=11.41,
        C_dry=61.5,
        H_dry=4.1,
        O_dry=21.8,
        nitrogen_fraction=0.6,
        sulfur_fraction=0.01,
        sulfur_volatile_frac=0.5,
        pyrolysis_tar_carbon_frac=0.2,
        fuel_type="coal",
        m_vm_in=0.05,
        m_moist_in=0.02,
        solid_shape=(4, 4),
        m_vm_in_classes=np.array([0.02, 0.01, 0.015, 0.005], dtype=float),
        m_moist_in_classes=np.array([0.01, 0.005, 0.003, 0.002], dtype=float),
        char_index=0,
        vm_index=1,
        moisture_index=2,
    )

    vm_release = -float(np.sum(bundle.solid_sink[:, 1]))
    moist_release = -float(np.sum(bundle.solid_sink[:, 2]))
    assert vm_release <= 0.05 + 1e-9
    assert moist_release <= 0.02 + 1e-9
    assert bundle.solid_sink[0, 1] < bundle.solid_sink[-1, 1]
    assert bundle.solid_sink[0, 2] < bundle.solid_sink[-1, 2]


def test_calc_drying_pyrolysis_sources_short_circuits_when_no_vm_or_moisture(monkeypatch: pytest.MonkeyPatch) -> None:
    def _should_not_run(*args, **kwargs):
        raise AssertionError("solve_drying_CN should not run when there is no VM/moisture feed")

    monkeypatch.setattr("src.core.cell_pyrolysis.solve_drying_CN", _should_not_run)

    bundle = calc_drying_pyrolysis_sources(
        tau=5.0,
        T=1173.15,
        P=2.5e6,
        d_p=1.0e-3,
        moisture_wt=16.9,
        ash_dry_wt=11.41,
        C_dry=61.5,
        H_dry=4.1,
        O_dry=21.8,
        nitrogen_fraction=0.6,
        sulfur_fraction=0.01,
        sulfur_volatile_frac=0.5,
        pyrolysis_tar_carbon_frac=0.2,
        fuel_type="coal",
        m_vm_in=0.0,
        m_moist_in=0.0,
        solid_shape=(1, 4),
        vm_index=1,
        moisture_index=2,
    )

    assert np.allclose(bundle.gas_source, 0.0)
    assert np.allclose(bundle.solid_sink, 0.0)
    assert bundle.x_dry == 0.0
    assert bundle.x_vm == 0.0
