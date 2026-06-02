from types import SimpleNamespace

import numpy as np
import pytest

from src.core.cell import Cell
from src.core.connectivity import (
    align_bottom_primary_gas_state_to_inlet_split,
    preproject_bottom_major_gas_state_to_local_balance,
    preproject_bottom_primary_gas_state_to_exchange_closure,
    preproject_bottom_total_gas_and_temperature_to_energy_closure,
    resolve_bottom_gas_inlet_dense_fraction,
    set_bottom_cell_feeds,
    snapshot_bottom_gas_inlet_split_from_vorabrechnung,
)
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


def _cfg(*, dense_frac=0.2, strategy="fixed"):
    return SimpleNamespace(
        gas_inlet_dense_frac=dense_frac,
        gas_inlet_split_strategy=strategy,
        O2_feed=10.0,
        H2O_feed=5.0,
        N2_feed=20.0,
        moisture_wt=10.0,
        ash_dry_wt=5.0,
        VM_daf=50.0,
        fuel_feed=1.0,
        n_age_classes=1,
        T_inlet=900.0,
    )


def test_bottom_gas_split_legacy_fixed_uses_configured_fraction():
    cell = Cell()
    frac = resolve_bottom_gas_inlet_dense_fraction(cell, _cfg(dense_frac=0.2, strategy="fixed"))

    assert frac == pytest.approx(0.2)


def test_bottom_gas_split_can_come_from_vorabrechnung_hydrodynamic_flux():
    cell = Cell()
    cell.u0 = 1.25
    cell.u_d = 0.50
    cell.eps_b = 0.40

    frac = resolve_bottom_gas_inlet_dense_fraction(
        cell,
        _cfg(dense_frac=0.2, strategy="precalc_hydrodynamic_flux"),
    )

    assert frac == pytest.approx(0.50 * (1.0 - 0.40) / 1.25)


def test_bottom_cell_feeds_use_hydrodynamic_split_for_dense_and_bubble_rows():
    cell = Cell()
    cell.u0 = 1.25
    cell.u_d = 0.50
    cell.eps_b = 0.40
    cfg = _cfg(dense_frac=0.2, strategy="precalc_hydrodynamic_flux")

    set_bottom_cell_feeds([cell], cfg)

    idx = GAS_SPECIES_INDEX
    dense_frac = 0.50 * (1.0 - 0.40) / 1.25
    assert cell.N_zu_d[idx["O2"]] == pytest.approx(cfg.O2_feed * dense_frac)
    assert cell.N_zu_b[idx["O2"]] == pytest.approx(cfg.O2_feed * (1.0 - dense_frac))


def test_bottom_gas_split_hydrodynamic_strategy_falls_back_before_precalc_exists():
    cell = Cell()
    frac = resolve_bottom_gas_inlet_dense_fraction(
        cell,
        _cfg(dense_frac=0.2, strategy="precalc_hydrodynamic_flux"),
    )

    assert frac == pytest.approx(0.2)


def test_bottom_gas_split_hydrodynamic_strategy_uses_frozen_vorabrechnung_snapshot():
    cell = Cell()
    cfg = _cfg(dense_frac=0.2, strategy="precalc_hydrodynamic_flux")
    cell.u0 = 1.25
    cell.u_d = 0.50
    cell.eps_b = 0.40
    snapshot_bottom_gas_inlet_split_from_vorabrechnung([cell], cfg)

    cell.u0 = 0.75
    cell.u_d = 0.45
    cell.eps_b = 0.20

    assert resolve_bottom_gas_inlet_dense_fraction(cell, cfg) == pytest.approx(0.50 * 0.60 / 1.25)


def test_bottom_primary_gas_state_alignment_preserves_total_and_uses_frozen_split():
    cell = Cell()
    cfg = _cfg(dense_frac=0.2, strategy="precalc_hydrodynamic_flux")
    cell._vorab_bottom_gas_inlet_dense_frac = 0.25
    idx = GAS_SPECIES_INDEX
    for sp, total in {"O2": 8.0, "H2O": 12.0, "N2": 20.0}.items():
        j = idx[sp]
        cell.N_d[j] = 0.5 * total
        cell.N_b[j] = 0.5 * total
    cell.N_d[idx["CO"]] = 3.0
    cell.N_b[idx["CO"]] = 7.0

    align_bottom_primary_gas_state_to_inlet_split([cell], cfg)

    for sp, total in {"O2": 8.0, "H2O": 12.0, "N2": 20.0}.items():
        j = idx[sp]
        assert cell.N_d[j] + cell.N_b[j] == pytest.approx(total)
        assert cell.N_d[j] / total == pytest.approx(0.25)
    assert cell.N_d[idx["CO"]] == pytest.approx(3.0)
    assert cell.N_b[idx["CO"]] == pytest.approx(7.0)


def test_bottom_primary_exchange_preprojection_reduces_bed0_primary_phase_stiffness():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX

    product_species = ("CO", "CO2", "H2", "CH4")
    product_before = {
        sp: (float(cell.N_d[idx[sp]]), float(cell.N_b[idx[sp]]))
        for sp in product_species
    }
    primary_totals = {
        sp: float(cell.N_d[idx[sp]] + cell.N_b[idx[sp]])
        for sp in ("O2", "H2O", "N2")
    }

    align_bottom_primary_gas_state_to_inlet_split([cell], cfg)
    before = np.asarray(cell.residuals(), dtype=float)
    before_primary_norm = np.linalg.norm(before[[idx["O2"], idx["H2O"], idx["N2"]]], ord=2)
    before_n2_pair = max(abs(float(before[idx["N2"]])), abs(float(before[idx["N2"] + N_GAS])))

    changed = preproject_bottom_primary_gas_state_to_exchange_closure([cell], cfg)

    after = np.asarray(cell.residuals(), dtype=float)
    after_primary_norm = np.linalg.norm(after[[idx["O2"], idx["H2O"], idx["N2"]]], ord=2)
    after_n2_pair = max(abs(float(after[idx["N2"]])), abs(float(after[idx["N2"] + N_GAS])))

    assert changed is True
    assert after_primary_norm < 0.75 * before_primary_norm
    assert after_n2_pair < 0.25 * before_n2_pair
    for sp, total in primary_totals.items():
        j = idx[sp]
        assert float(cell.N_d[j] + cell.N_b[j]) == pytest.approx(total)
    for sp, (dense_before, bubble_before) in product_before.items():
        j = idx[sp]
        assert float(cell.N_d[j]) == pytest.approx(dense_before)
        assert float(cell.N_b[j]) == pytest.approx(bubble_before)


def test_bottom_major_gas_preprojection_reduces_bed0_product_residual_without_moving_primary_totals():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX

    preproject_bottom_primary_gas_state_to_exchange_closure([cell], cfg)
    before = np.asarray(cell.residuals(), dtype=float)
    before_gas_max = max(
        float(np.max(np.abs(before[:N_GAS]))),
        float(np.max(np.abs(before[N_GAS : 2 * N_GAS]))),
    )
    primary_totals = {
        sp: float(cell.N_d[idx[sp]] + cell.N_b[idx[sp]])
        for sp in ("O2", "H2O", "N2")
    }

    changed = preproject_bottom_major_gas_state_to_local_balance([cell], cfg)

    after = np.asarray(cell.residuals(), dtype=float)
    after_gas_max = max(
        float(np.max(np.abs(after[:N_GAS]))),
        float(np.max(np.abs(after[N_GAS : 2 * N_GAS]))),
    )

    assert changed is True or after_gas_max <= before_gas_max
    if changed:
        assert after_gas_max < 0.75 * before_gas_max
    for sp, total in primary_totals.items():
        j = idx[sp]
        assert float(cell.N_d[j] + cell.N_b[j]) == pytest.approx(total)
    assert float(cell.N_d[idx["CO"]] + cell.N_b[idx["CO"]]) > 1e-6
    assert float(cell.N_d[idx["H2"]] + cell.N_b[idx["H2"]]) > 1e-6


def test_bottom_joint_gas_temperature_preprojection_reduces_energy_and_total_gas_residual():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    cell = reactor.cells[0]
    diag = getattr(cell, "_bottom_joint_preprojection_diag", None)

    assert diag is not None
    assert diag["final_energy_abs_W"] < 0.75 * diag["initial_energy_abs_W"]
    assert diag["final_combined_gas_max_mol_s"] < 0.75 * diag["initial_combined_gas_max_mol_s"]
    assert diag["final_gas_max_mol_s"] <= diag["initial_gas_max_mol_s"]
    assert float(cell.T) >= float(cell._nr_temperature_min_K) - 1e-9
    assert float(cell.T) <= float(cell._nr_temperature_max_K) + 1e-9

    idx = GAS_SPECIES_INDEX
    # Direct drying/pyrolysis gas remains suspension/dense sourced; this
    # preprojection may seed dense outlet support, but must not seed TAR/NH3 into
    # the bubble outlet as a hidden source.
    for sp in ("NH3", "TAR1", "TAR2"):
        j = idx[sp]
        assert float(cell.N_b[j]) == pytest.approx(0.0, abs=1e-12)


def test_bed0_ch4_oxidation_is_limited_by_available_ch4_supply():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    reactor._apply_all_bc_for_nr()
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX

    # A tiny positive CH4 outlet should not allow the R6/R7 consumption terms to
    # consume tens of mol/s of methane merely because O2 is abundant.
    cell.N_d[idx["CH4"]] = 8.0e-4
    cell.N_b[idx["CH4"]] = 2.0e-4
    cell._thermo_cache_valid = False
    cell._hydro_cache_valid = False

    cell.residuals()

    dense_ch4_supply = (
        float(cell.N_zu_d[idx["CH4"]])
        + float(cell.N_d_in[idx["CH4"]])
        + max(float(cell.N_ex[idx["CH4"]]), 0.0)
        + float(cell.N_rez_d[idx["CH4"]])
        + float(cell._vm_gas_source_cache[idx["CH4"]])
    )
    bubble_ch4_supply = (
        float(cell.N_zu_b[idx["CH4"]])
        + float(cell.N_b_in[idx["CH4"]])
        + max(-float(cell.N_ex[idx["CH4"]]), 0.0)
        + float(cell.N_rez_b[idx["CH4"]])
    )
    dense_consumption = max(float(cell._vm_gas_source_cache[idx["CH4"]]) - float(cell.R_gas_d[idx["CH4"]]), 0.0)
    bubble_consumption = max(-float(cell.R_gas_b[idx["CH4"]]), 0.0)

    assert dense_consumption <= 1.01 * max(dense_ch4_supply, 1e-6)
    assert bubble_consumption <= 1.01 * max(bubble_ch4_supply, 1e-6)
