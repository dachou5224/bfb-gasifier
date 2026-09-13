"""P5 续：热解碳分率接线。偏碳池在 koks=0 时 An=b 不可行。"""

from __future__ import annotations

import pytest

from src.core.cell_pyrolysis import allocate_pyrolysis_products_hamel
from src.core.pyrolysis_carbon_split import resolve_pyrolysis_carbon_fractions
from src.core.species import TAR_SURROGATE_FORMULA, get_tar_component_mapping


def _gas_tar_atoms(prod: dict[str, float], fuel_type: str = "coal") -> tuple[float, float, float]:
    tar_map = get_tar_component_mapping(fuel_type)
    n_c = prod["CO"] + prod["CO2"] + prod["CH4"]
    n_h = 2.0 * prod["H2"] + 2.0 * prod["H2O"] + 4.0 * prod["CH4"]
    for label in ("TAR1", "TAR2"):
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map[label]]
        n_c += c_tar * prod[label]
        n_h += h_tar * prod[label]
    n_o = prod["CO"] + 2.0 * prod["CO2"] + prod["H2O"]
    return n_c, n_h, n_o


def test_fixed_zero_koks_makes_sanity_pool_infeasible() -> None:
    f_koks, f_tar = resolve_pyrolysis_carbon_fractions(
        fuel_type="coal",
        yield_mode="fixed",
        P_Pa=2.5e6,
        T_K=1200.0,
        pyrolysis_tar_carbon_frac=0.2,
    )
    assert f_koks == 0.0
    assert f_tar == 0.2
    with pytest.raises(RuntimeError, match="An=b"):
        allocate_pyrolysis_products_hamel(
            T_K=1200.0,
            P_Pa=2.5e6,
            nC=10.0,
            nH=12.0,
            nO=4.0,
            fuel_type="coal",
            pyrolysis_tar_carbon_frac=f_tar,
            pyrolysis_koks_carbon_frac=f_koks,
        )


def test_hegermann_split_makes_sanity_pool_feasible_without_overspend() -> None:
    f_koks, f_tar = resolve_pyrolysis_carbon_fractions(
        fuel_type="coal",
        yield_mode="hegermann",
        P_Pa=2.5e6,
        T_K=1200.0,
        pyrolysis_tar_carbon_frac=0.2,
    )
    assert f_koks > 0.5
    assert 0.0 < f_tar < 0.3
    prod = allocate_pyrolysis_products_hamel(
        T_K=1200.0,
        P_Pa=2.5e6,
        nC=10.0,
        nH=12.0,
        nO=4.0,
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=f_tar,
        pyrolysis_koks_carbon_frac=f_koks,
    )
    n_c, n_h, n_o = _gas_tar_atoms(prod)
    assert n_c <= 10.0 + 1e-9
    assert n_h <= 12.0 + 1e-9
    assert n_o <= 4.0 + 1e-9
    assert n_o == pytest.approx(4.0, abs=1e-9)
    assert n_h == pytest.approx(12.0, abs=1e-9)
    # Koks 碳留在固相，气相+焦油碳必须明显少于输入
    assert n_c < 10.0 * (1.0 - 0.5)


def test_resolve_rejects_unknown_mode() -> None:
    with pytest.raises(ValueError, match="pyrolysis_tar_yield_mode"):
        resolve_pyrolysis_carbon_fractions(
            fuel_type="coal",
            yield_mode="dump_surplus_c",
            P_Pa=2.5e6,
            T_K=1200.0,
        )
