from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_pyrolysis import (
    _allocate_pyrolysis_tar_and_minors,
    allocate_pyrolysis_products_elemental,
    allocate_pyrolysis_products_hamel,
    calc_drying_pyrolysis_sources,
)
from tests.validation_case_utils import load_case_LU
from src.core.species import TAR_SURROGATE_FORMULA, calc_tar_surrogate_fractions, get_tar_component_mapping


def test_allocate_pyrolysis_products_elemental_preserves_element_budgets() -> None:
    # 无 C₂H₄ 时偏碳池 An=b 不可行；扣 Koks 后与生产 Hegermann 路径同形
    prod = allocate_pyrolysis_products_elemental(
        nC=10.0,
        nH=12.0,
        nO=4.0,
        nN=0.2,
        nS=0.1,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=0.2,
        pyrolysis_koks_carbon_frac=0.55,
        T_K=1100.0,
        P_Pa=2.5e6,
    )

    tar_map = get_tar_component_mapping("coal")
    nC_out = prod["CO"] + prod["CO2"] + prod["CH4"]
    nH_out = (
        2.0 * prod["H2"]
        + 2.0 * prod["H2O"]
        + 4.0 * prod["CH4"]
        + 3.0 * prod["NH3"]
        + 2.0 * prod["H2S"]
    )
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
        pyrolysis_tar_yield_mode="hegermann",
    )

    from src.core.species import N_GAS

    assert bundle.gas_source.shape == (N_GAS,)
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
        # fixed + 仅扣 Teer、Koks=0 时气相 C 过多，An=b 不可行；Hamel 路径用 Hegermann 同时扣 Koks/Teer
        pyrolysis_tar_yield_mode="hegermann",
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


def _major_gas_element_rates(prod: dict[str, float]) -> dict[str, float]:
    from src.core.elemental_ledger import gas_element_molar_rates
    from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX

    rates = np.zeros(len(GAS_SPECIES))
    for sp in ("CO", "CO2", "CH4", "H2", "H2O"):
        rates[GAS_SPECIES_INDEX[sp]] = prod.get(sp, 0.0)
    return gas_element_molar_rates(rates, fuel_type="coal")


def test_apply_pyrolysis_co_reduce_h2_inject_reduces_co_without_ch4_touch():
    from src.core.cell_pyrolysis import apply_pyrolysis_co_reduce_h2_inject

    prod = {"CO": 6.0, "CH4": 4.0, "H2": 1.0, "H2O": 2.0, "CO2": 0.0}
    out = apply_pyrolysis_co_reduce_h2_inject(prod, 0.5, h2_per_trimmed_co=2.0)
    assert out["CO"] == pytest.approx(3.0)
    assert out["CH4"] == pytest.approx(4.0)
    assert out["H2"] == pytest.approx(7.0)
    assert out["H2O"] == pytest.approx(5.0)


def test_apply_pyrolysis_ch4_to_h2_co_split_preserves_ch_and_o():
    from src.core.cell_pyrolysis import apply_pyrolysis_ch4_to_h2_co_split

    prod = {"CO": 6.0, "CH4": 4.0, "H2": 0.0, "CO2": 0.0, "H2O": 1.0}
    out = apply_pyrolysis_ch4_to_h2_co_split(prod, 0.5)
    assert out["CH4"] == pytest.approx(2.0)
    assert out["H2"] == pytest.approx(4.0)
    assert out["CO"] == pytest.approx(8.0)
    # C: 6+4=10 vs 8+2=10; H: 16 vs 4+8=12... CH4 4*4=16 H, H2 4*2=8, CH4 left 2*4=8, total 16 ok


def test_allocate_pyrolysis_hamel_gibbs_lu_respects_tar_gas_pool() -> None:
    """LU：Koks+Teer 扣完后，§4.3 Gibbs 用完气相 C/H/O（An=b）。"""
    from tests.validation_case_utils import load_case_LU
    from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
    from src.thermal.hegermann_yield import xc_koks_hegermann, xc_teer_hegermann

    case = load_case_LU()
    fuel_kg_s = case["fuel_feed"] / 3600.0
    ash_frac = case["ash_dry_wt"] / 100.0
    dry = fuel_kg_s * (1.0 - case["moisture_wt"] / 100.0)
    daf = dry * (1.0 - ash_frac)
    vm = daf * case["VM_daf"] / 100.0
    to_daf = 1.0 / max(1.0 - ash_frac, 1e-9)

    nC = vm * case["C_dry"] / 100.0 * to_daf / ATOMIC_MASS_KG_PER_MOL["C"]
    nH = vm * case["H_dry"] / 100.0 * to_daf / ATOMIC_MASS_KG_PER_MOL["H"]
    nO = vm * case["O_dry"] / 100.0 * to_daf / ATOMIC_MASS_KG_PER_MOL["O"]
    T_K = 1173.0
    xk = xc_koks_hegermann(case["P"], T_K)
    xt = xc_teer_hegermann(case["P"], T_K)

    _, c_gas, h_gas, o_gas = _allocate_pyrolysis_tar_and_minors(
        nC=nC,
        nH=nH,
        nO=nO,
        nN=0.0,
        nS=0.0,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=xt,
        pyrolysis_koks_carbon_frac=xk,
        target_tar_hc_ratio=None,
    )

    prod = allocate_pyrolysis_products_hamel(
        T_K=T_K,
        P_Pa=case["P"],
        nC=nC,
        nH=nH,
        nO=nO,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=xt,
        pyrolysis_koks_carbon_frac=xk,
    )
    assert prod["CO"] + prod["CH4"] > 0.0

    out = _major_gas_element_rates(prod)
    assert out["C"] == pytest.approx(c_gas, abs=1e-6)
    assert out["H"] == pytest.approx(h_gas, abs=1e-6)
    assert out["O"] == pytest.approx(o_gas, abs=1e-6)


def test_allocate_pyrolysis_hamel_gibbs_element_closure_after_tar() -> None:
    """§4.3 Gibbs 须用完 Teer/Koks 分出后的 C/H/O 气相池。"""
    nC, nH, nO = 10.0, 12.0, 4.0
    f_koks = 0.55
    _, c_gas, h_gas, o_gas = _allocate_pyrolysis_tar_and_minors(
        nC=nC,
        nH=nH,
        nO=nO,
        nN=0.0,
        nS=0.0,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=0.2,
        pyrolysis_koks_carbon_frac=f_koks,
        target_tar_hc_ratio=None,
    )
    prod = allocate_pyrolysis_products_hamel(
        T_K=1100.0,
        P_Pa=2.5e6,
        nC=nC,
        nH=nH,
        nO=nO,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=0.2,
        pyrolysis_koks_carbon_frac=f_koks,
    )

    out = _major_gas_element_rates(prod)
    assert out["C"] == pytest.approx(c_gas, abs=1e-6)
    assert out["H"] == pytest.approx(h_gas, abs=1e-6)
    assert out["O"] == pytest.approx(o_gas, abs=1e-6)
