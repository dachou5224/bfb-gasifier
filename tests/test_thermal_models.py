"""干燥/热解子模型测试（Chapter 4 对应方程）。"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def test_corrected_evaporation_enthalpy_larger_than_latent():
    from src.thermal.drying import corrected_evaporation_enthalpy, H_EVAP

    hvp = corrected_evaporation_enthalpy(w0_tr=0.2, T0=300.0, T_e=373.15)
    assert hvp > H_EVAP


def test_nusselt_particle_positive():
    from src.thermal.drying import nusselt_particle

    nu = nusselt_particle(Re=50.0, Pr=0.7)
    assert nu > 2.0


def test_saturation_temperature_water_increases_with_pressure():
    from src.thermal.drying import saturation_temperature_water

    t_atm = saturation_temperature_water(101_325.0)
    t_25bar = saturation_temperature_water(2.5e6)
    assert 370.0 <= t_atm <= 376.0
    assert 490.0 <= t_25bar <= 510.0
    assert t_25bar > t_atm


def test_daem_default_coal_params_match_chapter4():
    from src.thermal.devolatilization import DAEM_PARAMS

    coal = DAEM_PARAMS["brown_coal"]
    assert np.isclose(coal.E0, 192_000.0)
    assert np.isclose(coal.sigma, 40_000.0)
    assert np.isclose(coal.A, 1.67e13)


def test_daem_conversion_radial_bounds():
    from src.thermal.devolatilization import daem_conversion_radial

    r = np.linspace(0.0, 1.0e-3, 8)
    t = np.linspace(0.0, 50.0, 120)
    # 构造温度场：表面更热，中心更冷
    T_rt = np.vstack(
        [
            700.0 + 300.0 * (ri / r[-1]) + 120.0 * (t / t[-1])
            for ri in r
        ]
    )
    x = daem_conversion_radial(T_rt, t, r)
    assert 0.0 <= x <= 1.0


def test_drying_returns_radial_history():
    from src.thermal.drying import solve_drying_CN

    out = solve_drying_CN(
        d_p=1.0e-3,
        T_bed=1100.0,
        t_total=1.0,
        Nr=8,
        Nt=20,
        return_history=True,
    )
    assert "T_history_rt" in out and "r_nodes" in out
    assert out["T_history_rt"].shape == (9, 21)
    assert len(out["r_nodes"]) == 9


def test_pyrolysis_allocator_element_conservation():
    from src.core.cell import Cell, SolidProps
    from src.core.species import get_tar_component_mapping, TAR_SURROGATE_FORMULA

    c = Cell(solid=SolidProps())
    c.fuel_type = "coal"
    prod = c._allocate_pyrolysis_products_elemental(nC=10.0, nH=12.0, nO=4.0)

    tar_map = get_tar_component_mapping("coal")
    nC_out = prod["CO"] + prod["CO2"] + prod["CH4"]
    nH_out = 2.0 * prod["H2"] + 2.0 * prod["H2O"] + 4.0 * prod["CH4"]
    for tar_label in ("TAR1", "TAR2"):
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map[tar_label]]
        nC_out += c_tar * prod[tar_label]
        nH_out += h_tar * prod[tar_label]
    nO_out = (
        prod["CO"] + 2.0 * prod["CO2"] + prod["H2O"]
    )

    assert nC_out <= 10.0 + 1e-9
    assert nH_out <= 12.0 + 1e-9
    assert nO_out <= 4.0 + 1e-9
    assert all(v >= 0.0 for v in prod.values())


def test_pyrolysis_allocator_co2_h2o_fallback_closure():
    from src.core.cell_pyrolysis import allocate_pyrolysis_products_elemental

    prod = allocate_pyrolysis_products_elemental(
        nC=1.0,
        nH=2.0,
        nO=4.0,  # 氧富集，要求 fallback 闭合
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=0.0,
    )
    assert prod["CO2"] > 0.0 or prod["H2O"] > 0.0
