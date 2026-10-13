from __future__ import annotations

import numpy as np
import pytest

from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    HEAT_LOSS_MODE_INLET_FRACTION,
    apply_freeboard_substep_gas_enthalpy_heat_loss,
    compound_segment_heat_loss_fraction,
)
from src.core.freeboard_segment import simulate_freeboard
from src.core.reactor import Reactor, ReactorConfig, resolve_freeboard_heat_loss_bundle
from src.core.species import GAS_SPECIES_INDEX, N_GAS, configure_tar_components_by_fuel


def _minimal_freeboard_kwargs(**overrides) -> dict:
    configure_tar_components_by_fuel("coal")
    idx = GAS_SPECIES_INDEX
    n_in = np.zeros(N_GAS, dtype=np.float64)
    n_in[idx["CO"]] = 1.0
    n_in[idx["H2"]] = 1.0
    n_in[idx["H2O"]] = 1.0
    n_in[idx["N2"]] = 3.0
    base = dict(
        N_in=n_in,
        T_in=1100.0,
        P=2.5e6,
        D_bed=0.6,
        H_freeboard=2.0,
        n_cells=2,
        u_b_bed_top=1.2,
        d_b_bed_top=0.12,
        eps_b_bed_top=0.3,
        eps_d_void_bed_top=0.5,
        rho_solid_bed_top=1700.0,
        d_p_classes_bed_top=np.array([7e-4]),
        m_char_classes_bed_top=np.array([0.02]),
        m_ash_classes_bed_top=np.array([0.01]),
        fuel_type="coal",
        trajectory_model="analytical_wirsum",
        enabled_reactions=("R5", "R6", "R7", "R8", "R12"),
        heat_loss_frac=0.0,
    )
    base.update(overrides)
    return base


def test_compound_segment_heat_loss_fraction_matches_legacy_formula() -> None:
    seg = compound_segment_heat_loss_fraction(0.1, 10)
    legacy = 1.0 - (1.0 - 0.1) ** (1.0 / 10.0)
    assert seg == pytest.approx(legacy)


def test_freeboard_substep_skips_inlet_fraction_when_H_non_positive() -> None:
    h0 = -1.0e5
    h1 = apply_freeboard_substep_gas_enthalpy_heat_loss(
        h0,
        heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION,
        heat_loss_frac_total=0.1,
        n_cells=2,
    )
    assert h1 == h0


def test_explicit_freeboard_substep_deducts_fixed_watts() -> None:
    h0 = 1.0e6
    h1 = apply_freeboard_substep_gas_enthalpy_heat_loss(
        h0,
        heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        heat_loss_frac_total=0.0,
        n_cells=2,
        n_substeps=12,
        n_slices=1,
        Q_wall_W=240.0,
        Q_heat_exchanger_W=0.0,
    )
    assert h1 == pytest.approx(h0 - 10.0)


def test_simulate_freeboard_explicit_heat_loss_lowers_exit_temperature() -> None:
    base = _minimal_freeboard_kwargs(heat_loss_frac=0.0)
    out_none = simulate_freeboard(**base)
    out_cooled = simulate_freeboard(
        **_minimal_freeboard_kwargs(
            heat_loss_frac=0.0,
            heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
            Q_wall_W=2.0e4,
            Q_heat_exchanger_W=5.0e3,
        )
    )
    assert out_cooled["exit_T"] < out_none["exit_T"]


def test_resolve_freeboard_heat_loss_bundle_splits_explicit_q_by_height() -> None:
    cfg = ReactorConfig(
        n_cells=4,
        H_bed=6.0,
        H_freeboard=3.0,
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        Q_wall_W_total=900.0,
        Q_heat_exchanger_W_total=300.0,
    )
    bundle = resolve_freeboard_heat_loss_bundle(cfg)
    assert bundle.heat_loss_mode == HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT
    assert bundle.Q_wall_W == pytest.approx(300.0)
    assert bundle.Q_heat_exchanger_W == pytest.approx(100.0)


def test_reactor_assigns_freeboard_cell_explicit_q() -> None:
    cfg = ReactorConfig(
        n_cells=2,
        n_freeboard_cells=2,
        H_bed=3.0,
        H_freeboard=3.0,
        thesis_mode=True,
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        Q_wall_W_total=200.0,
        Q_heat_exchanger_W_total=100.0,
    )
    reactor = Reactor(cfg)
    assert len(reactor.freeboard_cells) == 2
    assert reactor.freeboard_cells[0].Q_wall_W == pytest.approx(50.0)
    assert reactor.freeboard_cells[0].Q_heat_exchanger_W == pytest.approx(25.0)
