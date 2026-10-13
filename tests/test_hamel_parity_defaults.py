"""Hamel P1 默认：bed holdup_transport + freeboard exact_hamel。"""

from __future__ import annotations

import pytest

from src.core.hamel_parity_defaults import (
    HAMEL_BED_SOLID_STATE_MODEL,
    HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE,
    HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL,
    HAMEL_GIBBS_INITIAL_GUESS_MODE,
    HAMEL_GIBBS_TOL_MODE,
    HAMEL_SIZE_MIGRATION_INVENTORY_MODE,
    LEGACY_BED_SOLID_STATE_MODEL,
    LEGACY_FEED_SOLID_SIZE_DISTRIBUTION_MODE,
    LEGACY_FREEBOARD_TRAJECTORY_COEFF_MODEL,
    LEGACY_GIBBS_INITIAL_GUESS_MODE,
    LEGACY_GIBBS_TOL_MODE,
    LEGACY_SIZE_MIGRATION_INVENTORY_MODE,
    PREPROJ_SOLID_HOLDUP_REACTION_RATE_MULTIPLIER,
    PREPROJ_PHASE_SPLIT_REACTION_RATE_MULTIPLIER,
    resolve_feed_solid_mass_fractions,
)
from src.core.connectivity.gas_inlet import set_bottom_cell_feeds
from src.core.reactor import Reactor, ReactorConfig
from src.thermodynamics.gibbs_hamel_reduced import ReducedHamelGibbsSolver
from src.thermodynamics.gibbs_minimizer import GibbsMinimizer


def test_gibbs_solver_defaults_use_hamel_a40_and_holub_tol() -> None:
    class _DummySpeciesDB:
        def get_mu0(self, species: str, T: float) -> float:
            _ = T
            return {"H2": 0.0, "H2O": -220_000.0, "O2": 0.0}[species]

        def get_atom_count(self, species: str, element: str) -> int:
            atoms = {
                "H2": {"H": 2, "O": 0},
                "H2O": {"H": 2, "O": 1},
                "O2": {"H": 0, "O": 2},
            }
            return int(atoms.get(species, {}).get(element, 0))

    minimizer = GibbsMinimizer(_DummySpeciesDB())
    _, min_diag = minimizer.solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        return_diag=True,
    )
    reduced = ReducedHamelGibbsSolver(_DummySpeciesDB())
    _, red_diag = reduced.solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        return_diag=True,
    )
    for diag in (min_diag, red_diag):
        assert diag["initial_guess_mode"] == HAMEL_GIBBS_INITIAL_GUESS_MODE
        assert diag["tol_mode"] == HAMEL_GIBBS_TOL_MODE
        assert diag["converged"] is True


def test_legacy_gibbs_modes_remain_opt_in() -> None:
    class _DummySpeciesDB:
        def get_mu0(self, species: str, T: float) -> float:
            _ = T
            return {"H2": 0.0, "H2O": -220_000.0, "O2": 0.0}[species]

        def get_atom_count(self, species: str, element: str) -> int:
            atoms = {
                "H2": {"H": 2, "O": 0},
                "H2O": {"H": 2, "O": 1},
                "O2": {"H": 0, "O": 2},
            }
            return int(atoms.get(species, {}).get(element, 0))

    _, diag = GibbsMinimizer(_DummySpeciesDB()).solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        tol_mode=LEGACY_GIBBS_TOL_MODE,
        initial_guess_mode=LEGACY_GIBBS_INITIAL_GUESS_MODE,
        return_diag=True,
    )
    assert diag["initial_guess_mode"] == LEGACY_GIBBS_INITIAL_GUESS_MODE
    assert diag["tol_mode"] == LEGACY_GIBBS_TOL_MODE


def test_reactor_config_defaults_to_hamel_solid_and_freeboard_coeffs() -> None:
    from src.core.energy_balance_defaults import HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT

    cfg = ReactorConfig()
    assert cfg.bed_solid_state_model == HAMEL_BED_SOLID_STATE_MODEL
    assert cfg.freeboard_trajectory_coeff_model == HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL
    assert cfg.size_migration_inventory_mode == HAMEL_SIZE_MIGRATION_INVENTORY_MODE
    assert cfg.energy_balance_heat_loss_mode == HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT


def test_reactor_bed_cells_use_configured_solid_state_model() -> None:
    reactor = Reactor(ReactorConfig(n_cells=3))
    assert all(cell.solid_state_model == HAMEL_BED_SOLID_STATE_MODEL for cell in reactor.cells)
    assert all(
        cell.size_migration_inventory_mode == HAMEL_SIZE_MIGRATION_INVENTORY_MODE
        for cell in reactor.cells
    )


def test_legacy_solid_and_freeboard_coeffs_remain_opt_in() -> None:
    from src.core.energy_balance_defaults import HEAT_LOSS_MODE_INLET_FRACTION

    cfg = ReactorConfig(
        n_cells=2,
        bed_solid_state_model=LEGACY_BED_SOLID_STATE_MODEL,
        freeboard_trajectory_coeff_model=LEGACY_FREEBOARD_TRAJECTORY_COEFF_MODEL,
        size_migration_inventory_mode=LEGACY_SIZE_MIGRATION_INVENTORY_MODE,
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION,
    )
    reactor = Reactor(cfg)
    assert all(cell.solid_state_model == LEGACY_BED_SOLID_STATE_MODEL for cell in reactor.cells)
    assert all(
        cell.size_migration_inventory_mode == LEGACY_SIZE_MIGRATION_INVENTORY_MODE
        for cell in reactor.cells
    )


def test_thesis_mode_forces_holdup_even_if_legacy_configured() -> None:
    cfg = ReactorConfig(
        n_cells=2,
        thesis_mode=True,
        bed_solid_state_model=LEGACY_BED_SOLID_STATE_MODEL,
    )
    reactor = Reactor(cfg)
    assert all(cell.solid_state_model == HAMEL_BED_SOLID_STATE_MODEL for cell in reactor.cells)


def test_feed_solid_mass_fractions_default_to_coarsest_class() -> None:
    frac = resolve_feed_solid_mass_fractions(n_size_classes=10, mode=HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE)
    assert float(frac.sum()) == pytest.approx(1.0)
    assert float(frac[-1]) == pytest.approx(1.0)
    assert float(frac[0]) == pytest.approx(0.0)


def test_set_bottom_cell_feeds_uses_mass_fractions_not_uniform_nk() -> None:
    cfg = ReactorConfig(n_cells=2, n_age_classes=10, fuel_feed=1.0)
    reactor = Reactor(cfg)
    set_bottom_cell_feeds(reactor.cells, cfg)
    bot = reactor.cells[0]
    assert float(bot.m_solid_zu[-1, 0]) > 0.0
    assert float(bot.m_solid_zu[0, 0]) == pytest.approx(0.0)


def test_reactor_config_init_solid_preprojection_default_off() -> None:
    cfg = ReactorConfig()
    assert cfg.feed_solid_size_distribution_mode == HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE
    assert cfg.nr_init_bottom_zone_solid_energy_preprojection_thesis is False
    assert cfg.nr_init_bottom_primary_exchange_preprojection_thesis is False
    assert cfg.nr_init_bottom_major_local_balance_preprojection_thesis is False
    assert cfg.nr_init_bottom_major_phase_split_preprojection_thesis is False
    assert cfg.nr_init_bottom_joint_preprojection_thesis is False
    assert cfg.preproj_solid_holdup_reaction_rate_multiplier_thesis == pytest.approx(
        PREPROJ_SOLID_HOLDUP_REACTION_RATE_MULTIPLIER
    )
    assert cfg.preproj_gas_reaction_rate_multiplier_thesis == pytest.approx(0.0)
    assert cfg.preproj_phase_split_reaction_rate_multiplier_thesis == pytest.approx(
        PREPROJ_PHASE_SPLIT_REACTION_RATE_MULTIPLIER
    )


def test_legacy_uniform_feed_distribution_opt_in() -> None:
    import numpy as np

    frac = resolve_feed_solid_mass_fractions(
        n_size_classes=4,
        mode=LEGACY_FEED_SOLID_SIZE_DISTRIBUTION_MODE,
    )
    assert np.allclose(frac, 0.25)
