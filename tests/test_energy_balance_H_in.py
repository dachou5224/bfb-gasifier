from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_balances import calc_energy_balance_H_in, calc_energy_balance_residual
from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_INLET_FRACTION,
    cell_inlet_fraction_heat_loss_W,
)
from src.core.species import configure_tar_components_by_fuel, GAS_SPECIES_INDEX


def test_calc_energy_balance_H_in_matches_residual_without_loss() -> None:
    configure_tar_components_by_fuel("coal")
    gas = np.zeros(11, dtype=np.float64)
    gas[4] = 2.0
    solid = np.array([[0.4, 0.2, 0.1, 0.3]], dtype=np.float64)
    kw = dict(
        N_b_in=gas * 0.3,
        N_d_in=gas * 0.7,
        T_in_gas=1100.0,
        N_zu_b=np.zeros_like(gas),
        N_zu_d=np.zeros_like(gas),
        T_zu_gas=1100.0,
        N_rez_b=np.zeros_like(gas),
        N_rez_d=np.zeros_like(gas),
        T_rez_gas=1100.0,
        m_solid_rez=np.zeros_like(solid),
        T_rez_solid=1100.0,
        m_solid_in=np.zeros_like(solid),
        T_in_solid=1100.0,
        m_solid_zu=solid,
        T_zu_solid=1100.0,
        N_b=gas * 0.3,
        N_d=gas * 0.7,
        m_solid=solid,
        T=1100.0,
        heat_loss_frac=0.0,
        ash_dry_wt=11.41,
        VM_daf=53.42,
        h_f_dry=-5.0e5,
        h_cache={},
    )
    H_in = calc_energy_balance_H_in(
        N_b_in=kw["N_b_in"],
        N_d_in=kw["N_d_in"],
        T_in_gas=kw["T_in_gas"],
        N_zu_b=kw["N_zu_b"],
        N_zu_d=kw["N_zu_d"],
        T_zu_gas=kw["T_zu_gas"],
        N_rez_b=kw["N_rez_b"],
        N_rez_d=kw["N_rez_d"],
        T_rez_gas=kw["T_rez_gas"],
        m_solid_rez=kw["m_solid_rez"],
        T_rez_solid=kw["T_rez_solid"],
        m_solid_in=kw["m_solid_in"],
        T_in_solid=kw["T_in_solid"],
        m_solid_zu=kw["m_solid_zu"],
        T_zu_solid=kw["T_zu_solid"],
        ash_dry_wt=kw["ash_dry_wt"],
        VM_daf=kw["VM_daf"],
        h_f_dry=kw["h_f_dry"],
        h_cache=kw["h_cache"],
    )
    res0 = calc_energy_balance_residual(**kw, heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION)
    assert res0 == 0.0
    res_frac = calc_energy_balance_residual(**{**kw, "heat_loss_frac": 0.1})
    assert res_frac == pytest.approx(-0.1 * max(H_in, 0.0))


class _FakeCell:
    heat_loss_frac = 0.1
    solid_enthalpy_mode = "hamel_eq28_component"
    T_in_gas = T_zu_gas = T_rez_gas = T_in_solid = T_zu_solid = T_rez_solid = 1100.0
    T_solid_auf_in = T_solid_ab_in = 1100.0
    _h_cache = {}

    def __init__(self) -> None:
        self.N_b_in = np.zeros(11, dtype=np.float64)
        self.N_d_in = np.zeros(11, dtype=np.float64)
        self.N_zu_b = np.zeros(11, dtype=np.float64)
        self.N_zu_d = np.zeros(11, dtype=np.float64)
        self.N_rez_b = np.zeros(11, dtype=np.float64)
        self.N_rez_d = np.zeros(11, dtype=np.float64)
        self.m_solid_rez = np.zeros((1, 4), dtype=np.float64)
        self.m_solid_in = np.zeros((1, 4), dtype=np.float64)
        self.m_solid_zu = np.zeros((1, 4), dtype=np.float64)
        self.m_solid_auf_in = np.zeros((1, 4), dtype=np.float64)
        self.m_solid_ab_in = np.zeros((1, 4), dtype=np.float64)

    class solid:
        ash_dry_wt = 11.41
        VM_daf = 53.42
        h_f_dry = -5.0e5


def test_cell_inlet_fraction_heat_loss_zero_when_H_in_non_positive() -> None:
    configure_tar_components_by_fuel("coal")
    cell = _FakeCell()
    assert cell_inlet_fraction_heat_loss_W(cell) == 0.0


def test_cell_inlet_fraction_heat_loss_scales_with_positive_H_in() -> None:
    configure_tar_components_by_fuel("coal")
    cell = _FakeCell()
    cell.N_b_in[GAS_SPECIES_INDEX["H2"]] = 50.0
    q1 = cell_inlet_fraction_heat_loss_W(cell)
    cell.N_b_in[GAS_SPECIES_INDEX["H2"]] = 100.0
    q2 = cell_inlet_fraction_heat_loss_W(cell)
    assert q1 > 0.0 and q2 > q1


def test_inlet_fraction_skips_heat_loss_when_H_in_non_positive() -> None:
    configure_tar_components_by_fuel("coal")
    gas = np.zeros(11, dtype=np.float64)
    # 液态水分生成焓为负；纯碳 char 不再把 H_in 做成负值。
    solid = np.array([[0.0, 0.0, 1.0, 0.0]], dtype=np.float64)
    kw = dict(
        N_b_in=gas,
        N_d_in=gas,
        T_in_gas=1100.0,
        N_zu_b=np.zeros_like(gas),
        N_zu_d=np.zeros_like(gas),
        T_zu_gas=1100.0,
        N_rez_b=np.zeros_like(gas),
        N_rez_d=np.zeros_like(gas),
        T_rez_gas=1100.0,
        m_solid_rez=np.zeros_like(solid),
        T_rez_solid=1100.0,
        m_solid_in=np.zeros_like(solid),
        T_in_solid=1100.0,
        m_solid_zu=solid,
        T_zu_solid=1100.0,
        N_b=gas,
        N_d=gas,
        m_solid=solid,
        T=1100.0,
        ash_dry_wt=11.41,
        VM_daf=53.42,
        h_f_dry=-5.0e5,
        h_cache={},
        heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION,
    )
    H_in = calc_energy_balance_H_in(
        N_b_in=kw["N_b_in"],
        N_d_in=kw["N_d_in"],
        T_in_gas=kw["T_in_gas"],
        N_zu_b=kw["N_zu_b"],
        N_zu_d=kw["N_zu_d"],
        T_zu_gas=kw["T_zu_gas"],
        N_rez_b=kw["N_rez_b"],
        N_rez_d=kw["N_rez_d"],
        T_rez_gas=kw["T_rez_gas"],
        m_solid_rez=kw["m_solid_rez"],
        T_rez_solid=kw["T_rez_solid"],
        m_solid_in=kw["m_solid_in"],
        T_in_solid=kw["T_in_solid"],
        m_solid_zu=kw["m_solid_zu"],
        T_zu_solid=kw["T_zu_solid"],
        ash_dry_wt=kw["ash_dry_wt"],
        VM_daf=kw["VM_daf"],
        h_f_dry=kw["h_f_dry"],
        h_cache=kw["h_cache"],
    )
    assert H_in < 0.0
    res0 = calc_energy_balance_residual(**{**kw, "heat_loss_frac": 0.0})
    res1 = calc_energy_balance_residual(**{**kw, "heat_loss_frac": 0.1})
    assert res0 == res1
