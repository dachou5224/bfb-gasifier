"""P0.2：关键 ReactorConfig knob 必须接到 Cell / residual，不能是死旋钮。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import S_CHAR
from src.core.connectivity.gas_inlet import (
    resolve_bottom_gas_inlet_dense_fraction,
    set_bottom_cell_feeds,
)
from src.core.freeboard.trajectory_common import (
    _haider_high_re_drag_constant,
    _hamel_exact_high_re_drag_constant,
    _wirsum_abc,
    _wirsum_abc_exact_audit,
)
from src.core.hamel_parity_defaults import HAMEL_BED_SOLID_STATE_MODEL
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.global_nr_solver import global_residual, pack_reactor
from tests.global_nr_fixtures import call_generate_initial_x0_phase1, phase1_axial_t_profile
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config

RATE_SCALE_ATTRS = ("r2_scale", "r4_scale", "r5_scale", "r7_scale", "r8_scale")


def _two_cell_reactor(**overrides) -> Reactor:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    for key, value in overrides.items():
        setattr(cfg, key, value)
    reactor = Reactor(cfg)
    call_generate_initial_x0_phase1(reactor, cfg, phase1_axial_t_profile(cfg))
    return reactor


@pytest.mark.parametrize("attr", RATE_SCALE_ATTRS)
def test_rate_scale_copies_from_reactor_config_to_cells(attr: str) -> None:
    reactor = _two_cell_reactor(**{attr: 0.37})
    for cell in reactor.cells:
        assert getattr(cell, attr) == pytest.approx(0.37)


@pytest.mark.parametrize("attr", RATE_SCALE_ATTRS)
def test_rate_scale_changes_reaction_source_when_reactants_present(attr: str) -> None:
    reactor = _two_cell_reactor(**{attr: 1.0})
    cell = reactor.cells[0]
    cell.T = 1200.0
    cell.m_solid[:, S_CHAR] = np.maximum(cell.m_solid[:, S_CHAR], 0.25)
    for sp in ("H2O", "CO2", "CO", "CH4", "H2", "O2"):
        j = GAS_SPECIES_INDEX[sp]
        cell.N_d[j] = max(float(cell.N_d[j]), 1.0)
        cell.N_b[j] = max(float(cell.N_b[j]), 0.2)
    cell.calc_hydrodynamics()
    setattr(cell, attr, 1.0)
    cell.calc_reactions()
    r_on = np.concatenate([cell.R_gas_d, cell.R_gas_b, cell.R_solid.ravel()])
    setattr(cell, attr, 0.0)
    cell.calc_reactions()
    r_off = np.concatenate([cell.R_gas_d, cell.R_gas_b, cell.R_solid.ravel()])
    assert not np.allclose(r_on, r_off, rtol=0.0, atol=0.0)


def test_extent_limiters_thesis_flag_copies_onto_cells() -> None:
    off = _two_cell_reactor(extent_limiters_enabled_thesis=False)
    on = _two_cell_reactor(extent_limiters_enabled_thesis=True)
    assert all(cell.extent_limiters_enabled is False for cell in off.cells)
    assert all(cell.extent_limiters_enabled is True for cell in on.cells)


def test_extent_limiters_change_reaction_rates_or_residual() -> None:
    reactor = _two_cell_reactor(extent_limiters_enabled_thesis=False)
    cell = reactor.cells[0]
    x = pack_reactor(reactor.cells)
    f_off = global_residual(x.copy(), reactor.cells, reactor._apply_all_bc_for_nr)
    cell.calc_reactions()
    r_off = np.concatenate([cell.R_gas_d, cell.R_gas_b, cell.R_solid.ravel()])
    for item in reactor.cells:
        item.extent_limiters_enabled = True
    f_on = global_residual(x.copy(), reactor.cells, reactor._apply_all_bc_for_nr)
    cell.calc_reactions()
    r_on = np.concatenate([cell.R_gas_d, cell.R_gas_b, cell.R_solid.ravel()])
    assert (not np.allclose(r_off, r_on)) or (not np.allclose(f_off, f_on))


def test_gas_inlet_split_strategy_fixed_vs_hydro_changes_dense_fraction() -> None:
    reactor = _two_cell_reactor(gas_inlet_split_strategy="fixed", gas_inlet_dense_frac=0.30)
    bot = reactor.cells[0]
    bot.u0 = 0.40
    bot.u_d = 0.20
    bot.eps_b = 0.25
    cfg_fixed = reactor.config
    frac_fixed = resolve_bottom_gas_inlet_dense_fraction(bot, cfg_fixed)
    cfg_fixed.gas_inlet_split_strategy = "precalc_hydrodynamic_flux"
    bot._vorab_bottom_gas_inlet_dense_frac = np.nan
    frac_hydro = resolve_bottom_gas_inlet_dense_fraction(bot, cfg_fixed)
    assert frac_fixed == pytest.approx(0.30)
    assert frac_hydro == pytest.approx(0.20 * 0.75 / 0.40)
    assert frac_hydro != pytest.approx(frac_fixed)


def test_gas_inlet_split_strategy_reaches_bottom_feeds() -> None:
    reactor = _two_cell_reactor(gas_inlet_split_strategy="fixed", gas_inlet_dense_frac=0.10)
    set_bottom_cell_feeds(reactor.cells, reactor.config)
    n2_idx = GAS_SPECIES_INDEX["N2"]
    bot = reactor.cells[0]
    total_n2 = float(bot.N_zu_d[n2_idx] + bot.N_zu_b[n2_idx])
    assert total_n2 > 0.0
    assert bot.N_zu_d[n2_idx] == pytest.approx(0.10 * total_n2)


def test_thesis_mode_forces_holdup_transport_on_cells() -> None:
    reactor = _two_cell_reactor(thesis_mode=True, bed_solid_state_model="legacy_stream")
    assert all(cell.solid_state_model == HAMEL_BED_SOLID_STATE_MODEL for cell in reactor.cells)


def test_freeboard_coeff_models_change_wirsum_abc() -> None:
    kwargs = dict(
        u_g=8.0,
        u_p=12.0,
        rho_g=0.8,
        rho_p=1400.0,
        d_p=2.0e-3,
        mu_g=4.0e-5,
        phi_s=0.75,
    )
    z_haider = _haider_high_re_drag_constant(0.75)
    z_hamel = _hamel_exact_high_re_drag_constant(0.75)
    assert z_haider != pytest.approx(z_hamel)
    a0, b0, c0 = _wirsum_abc(**kwargs)
    a1, b1, c1, *_ = _wirsum_abc_exact_audit(**kwargs)
    assert (a0, b0, c0) != pytest.approx((a1, b1, c1))


def test_holdup_k0_anchor_copies_from_reactor_config() -> None:
    reactor = _two_cell_reactor(holdup_k0_anchor="none")
    assert all(cell.holdup_k0_anchor == "none" for cell in reactor.cells)


def test_holdup_k0_anchor_none_changes_cell_residual_when_K_zero() -> None:
    from src.core.cell import Cell

    cell = Cell()
    cell.solid_state_model = "holdup_transport"
    cell.K_solid_auf[:] = 0.0
    cell.K_solid_ab[:] = 0.0
    cell.m_solid[0, S_CHAR] = 0.25
    cell.m_solid_auf_in[0, S_CHAR] = 1.0
    cell.holdup_k0_anchor = "engineering"
    f_eng = cell.calc_solid_balance()
    cell.holdup_k0_anchor = "none"
    f_none = cell.calc_solid_balance()
    assert f_eng[0, S_CHAR] == pytest.approx(0.75)
    assert f_none[0, S_CHAR] == pytest.approx(1.0)


def test_holdup_k0_anchor_rejects_unknown_value() -> None:
    from src.core.reactor import ReactorConfig

    with pytest.raises(ValueError, match="holdup_k0_anchor"):
        ReactorConfig(holdup_k0_anchor="clip")
