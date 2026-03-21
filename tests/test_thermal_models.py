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
    c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map["TAR1"]]

    nC_out = (
        prod["CO"] + prod["CO2"] + prod["CH4"] + c_tar * prod["TAR1"]
    )
    nH_out = (
        2.0 * prod["H2"] + 2.0 * prod["H2O"] + 4.0 * prod["CH4"] + h_tar * prod["TAR1"]
    )
    nO_out = (
        prod["CO"] + 2.0 * prod["CO2"] + prod["H2O"]
    )

    assert nC_out <= 10.0 + 1e-9
    assert nH_out <= 12.0 + 1e-9
    assert nO_out <= 4.0 + 1e-9
    assert all(v >= 0.0 for v in prod.values())
