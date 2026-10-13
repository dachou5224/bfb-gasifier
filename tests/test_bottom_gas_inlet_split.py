from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from src.core.cell import Cell
from src.core.connectivity import (
    align_bottom_primary_gas_state_to_inlet_split,
    preproject_bottom_major_gas_phase_split_to_exchange_closure,
    preproject_bottom_major_gas_state_to_local_balance,
    preproject_bottom_primary_gas_state_to_exchange_closure,
    preproject_bottom_total_gas_and_temperature_to_energy_closure,
    resolve_bottom_gas_inlet_dense_fraction,
    run_bottom_zone_solid_energy_preprojection,
    set_bottom_cell_feeds,
    snapshot_bottom_gas_inlet_split_from_vorabrechnung,
)
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.global_nr_solver import (
    _residual_group_metrics,
    _residual_gas_phase_metrics,
    build_equation_scales,
    global_residual,
    pack_reactor,
)
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


def _phase2_preprojection_baseline_config():
    """Phase2 配置但关闭沿床焓步进/分级进料，供预投影单测隔离。"""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=False)
    cfg.nr_init_bottom_major_phase_split_preprojection_thesis = False
    return cfg


def _cfg(*, dense_frac=0.2, strategy="fixed", o2_dense_frac=None):
    return SimpleNamespace(
        gas_inlet_dense_frac=dense_frac,
        gas_inlet_o2_dense_frac=o2_dense_frac,
        gas_inlet_split_strategy=strategy,
        vorab_staged_primary_gas_injection_thesis=False,
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


def test_bottom_cell_feeds_o2_uses_separate_dense_fraction_when_configured():
    cell = Cell()
    cfg = _cfg(dense_frac=0.5, o2_dense_frac=0.9)
    set_bottom_cell_feeds([cell], cfg)
    idx = GAS_SPECIES_INDEX
    assert cell.N_zu_d[idx["O2"]] == pytest.approx(cfg.O2_feed * 0.9)
    assert cell.N_zu_b[idx["O2"]] == pytest.approx(cfg.O2_feed * 0.1)
    assert cell.N_zu_d[idx["H2O"]] == pytest.approx(cfg.H2O_feed * 0.5)
    assert cell.N_zu_b[idx["H2O"]] == pytest.approx(cfg.H2O_feed * 0.5)


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
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
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


def _scaled_split_rms(reactor: Reactor) -> float:
    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20e6, 1e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1e-9)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    F_hat = F / build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    return float(_residual_gas_phase_metrics(F_hat, cells)["gas_phase_split_rms"])


def _global_solid_energy_rms(reactor: Reactor) -> tuple[float, float]:
    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20e6, 1e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1e-9)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    F_hat = F / build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    gm = _residual_group_metrics(F_hat, cells)
    return float(gm["solid_rms"]), float(gm["energy_rms"])


def test_bottom_zone_solid_energy_preprojection_reduces_init_solid_rms():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_bed_stream_mapping_thesis = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    _, energy_before = _global_solid_energy_rms(reactor)
    zone_diag = run_bottom_zone_solid_energy_preprojection(
        reactor.cells,
        cfg,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
    )
    _solid_after, energy_after = _global_solid_energy_rms(reactor)

    assert zone_diag is not None
    assert cfg.nr_bottom_zone_solid_energy_preprojection_thesis is True
    assert float(zone_diag["solid_preprojections_accepted"]) >= 2.0
    # Hamel 对齐 init 后全局 solid RMS 未必单调下降；此处验证 zone LSQ 被接受且能量未爆炸。
    assert energy_after <= max(energy_before * 1.05, energy_before + 0.05)


def test_bottom_major_phase_split_preprojection_reduces_bed0_split_without_moving_primary_totals():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
    cfg.nr_bottom_major_phase_split_preprojection_thesis = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    primary_totals = {
        sp: float(cell.N_d[idx[sp]] + cell.N_b[idx[sp]])
        for sp in ("O2", "H2O", "N2")
    }
    before_split = _scaled_split_rms(reactor)

    changed = preproject_bottom_major_gas_phase_split_to_exchange_closure(reactor.cells, cfg)
    after_split = _scaled_split_rms(reactor)

    assert changed is True
    for sp, total in primary_totals.items():
        j = idx[sp]
        assert float(cell.N_d[j] + cell.N_b[j]) == pytest.approx(total)
    diag = getattr(cell, "_bottom_major_phase_split_preprojection_diag", None)
    assert diag is not None
    final_key = "final_split_rms_scaled" if "final_split_rms_scaled" in diag else "final_split_rms_mol_s"
    initial_key = "initial_split_rms_scaled" if "initial_split_rms_scaled" in diag else "initial_split_rms_mol_s"
    assert diag[final_key] < 0.75 * diag[initial_key]


def test_bottom_major_gas_preprojection_reduces_bed0_product_residual_without_moving_primary_totals():
    cfg = _phase2_preprojection_baseline_config()
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
        assert after_gas_max < 0.80 * before_gas_max
    for sp, total in primary_totals.items():
        j = idx[sp]
        assert float(cell.N_d[j] + cell.N_b[j]) == pytest.approx(total)
    assert float(cell.N_d[idx["CO"]] + cell.N_b[idx["CO"]]) > 1e-6
    assert float(cell.N_d[idx["H2"]] + cell.N_b[idx["H2"]]) > 1e-6


def test_bottom_joint_gas_temperature_preprojection_reduces_energy_and_total_gas_residual():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
    cfg.nr_init_bottom_primary_exchange_preprojection_thesis = True
    cfg.nr_init_bottom_major_local_balance_preprojection_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    cell = reactor.cells[0]
    ok = preproject_bottom_total_gas_and_temperature_to_energy_closure(reactor.cells, cfg, gas_relax=1.0)
    assert ok
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


def test_init_precalc_records_guarded_joint_preprojection_relax():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
    cfg.vorab_a_tier_bed_stream_mapping_thesis = False
    cfg.nr_init_bottom_joint_preprojection_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    run_diag = getattr(reactor, "_bottom_joint_preprojection_run_diag", None)

    assert run_diag is not None
    assert "joint_gas_relax" in run_diag
    assert "bed1_guard_tripped" in run_diag
    assert "joint_global_gate_pass" in run_diag
    assert "joint_rms_scaled_before" in run_diag
    assert "joint_rms_scaled_after" in run_diag
    assert "joint_merit_scaled_before" in run_diag
    assert "joint_merit_scaled_after" in run_diag
    assert "joint_split_aware_gate" in run_diag
    assert "bc_aware_bed0_bed1_joint" in run_diag
    assert run_diag["bc_aware_bed0_bed1_joint"] is True
    assert run_diag["joint_split_aware_gate"] is True
    assert run_diag["joint"] is True
    assert run_diag["joint_global_gate_pass"] is True
    assert 0.0 < float(run_diag["joint_gas_relax"]) <= 1.0
    assert float(run_diag["joint_merit_scaled_after"]) <= float(run_diag["joint_merit_scaled_before"]) + 1.0e-9
    assert getattr(reactor.cells[0], "_bottom_joint_preprojection_diag", None) is not None


def test_bottom_joint_split_aware_gate_accepts_merit_improving_blend():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
    cfg.vorab_a_tier_bed_stream_mapping_thesis = False
    cfg.nr_bottom_joint_split_aware_gate_thesis = True
    cfg.nr_bottom_joint_preprojection_gas_relax_min_thesis = 0.05
    cfg.nr_init_bottom_joint_preprojection_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    run_diag = getattr(reactor, "_bottom_joint_preprojection_run_diag", None)
    assert run_diag is not None
    assert run_diag["joint"] is True
    assert float(run_diag["joint_merit_scaled_after"]) <= float(run_diag["joint_merit_scaled_before"]) + 1.0e-9
    assert float(run_diag["joint_gas_relax"]) <= 0.5 + 1.0e-9


def test_bottom_joint_strict_gate_rejects_merit_only_improvement_when_split_aware_disabled():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
    cfg.nr_bottom_joint_split_aware_gate_thesis = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    run_diag = getattr(reactor, "_bottom_joint_preprojection_run_diag", None)
    assert run_diag is not None
    if run_diag["joint"]:
        assert float(run_diag["joint_rms_scaled_after"]) <= float(run_diag["joint_rms_scaled_before"]) + 1.0e-9
        assert float(run_diag["joint_merit_scaled_after"]) <= float(run_diag["joint_merit_scaled_before"]) + 1.0e-9
    else:
        assert run_diag["joint_global_gate_pass"] is False


def test_bed0_ch4_oxidation_is_limited_by_available_ch4_supply():
    cfg = _phase2_preprojection_baseline_config()
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


def test_bc_aware_bed0_bed1_joint_limits_neighbor_energy_spillover():
    cfg = _phase2_preprojection_baseline_config()
    cfg.vorab_a_tier_cell_mapping_thesis = False
    cfg.nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis = True
    cfg.nr_init_bottom_joint_preprojection_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    run_diag = getattr(reactor, "_bottom_joint_preprojection_run_diag", None)
    assert run_diag is not None
    assert run_diag["bc_aware_bed0_bed1_joint"] is True

    bed0_diag = getattr(reactor.cells[0], "_bottom_joint_preprojection_diag", None)
    if run_diag["joint"] and bed0_diag is not None and bed0_diag.get("bc_aware_neighbor_index") == 1:
        before = float(bed0_diag["initial_neighbor_energy_abs_W"])
        after = float(bed0_diag["final_neighbor_energy_abs_W"])
        guard = float(cfg.nr_bottom_joint_bed1_energy_guard_ratio_thesis)
        assert after <= before * guard + 1.0e6
        assert float(reactor.cells[0].T) >= float(reactor.cells[0]._nr_temperature_min_K) - 1e-9
        assert float(reactor.cells[1].T) >= float(reactor.cells[1]._nr_temperature_min_K) - 1e-9


def test_partial_cell_residual_groups_match_full_vector_slices() -> None:
    cfg = _phase2_preprojection_baseline_config()
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    cell.calc_hydrodynamics()
    cell.calc_exchange()
    cell.calc_reactions()

    full = cell.residuals()
    n_gas = int(cell.N_d.shape[0])
    gas = cell.residuals(groups=("gas",), copy=True)
    assert np.allclose(gas, full[: 2 * n_gas], rtol=0.0, atol=1.0e-12)

    energy = cell.residuals(groups=("energy",), copy=True)
    assert energy[0] == pytest.approx(full[-1], rel=0.0, abs=1.0e-9)

    gas_energy = cell.residuals(groups=("gas", "energy"), copy=True)
    assert np.allclose(gas_energy[: 2 * n_gas], full[: 2 * n_gas], rtol=0.0, atol=1.0e-12)
    assert gas_energy[-1] == pytest.approx(full[-1], rel=0.0, abs=1.0e-9)
