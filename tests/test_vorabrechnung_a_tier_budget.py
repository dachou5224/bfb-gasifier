"""A 档 Vorabrechnung 预算链：tau_bed、累积 O2/干燥/热解语义。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import S_VM
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.vorabrechnung import (
    apply_vorabrechnung_bed_cell_mapping,
    apply_vorabrechnung_bed_stream_mapping,
    build_vorabrechnung_cell_budgets,
    derive_nr_temperature_seed_from_vorab_budgets,
    estimate_axial_T_profile,
    generate_initial_x0,
    vorabrechnung_tau_for_cell,
)
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import (
    build_phase1_htw_lu_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
)


def _prepare_reactor_for_budget(cfg, *, n_cells: int = 4) -> Reactor:
    cfg.n_cells = n_cells
    reactor = Reactor(cfg)
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    for cell in reactor.cells:
        cell.calc_hydrodynamics()
    return reactor


def test_build_vorabrechnung_cell_budgets_uses_tau_bed_not_per_cell_dh_over_umf():
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = _prepare_reactor_for_budget(cfg, n_cells=4)
    tau_bed = vorabrechnung_tau_for_cell(reactor.cells[0])

    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )

    assert len(budgets) == len(reactor.cells)
    assert budgets[0].tau_cumulative_s < budgets[-1].tau_cumulative_s
    top = budgets[-1]
    assert abs(top.tau_cumulative_s - top.axial_xi * tau_bed) <= 1e-6
    per_cell_tau = reactor.cells[0].geo.dh / max(reactor.cells[0].u_mf, 1e-3)
    assert top.tau_cumulative_s > 2.0 * per_cell_tau


def test_build_vorabrechnung_cell_budgets_o2_depletes_with_cumulative_vm():
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = _prepare_reactor_for_budget(cfg, n_cells=5)

    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )

    o2_bottom = budgets[0].o2_remaining_mol_s
    o2_top = budgets[-1].o2_remaining_mol_s
    assert o2_bottom >= o2_top
    assert budgets[-1].nC_vm_mol_s >= budgets[0].nC_vm_mol_s
    assert budgets[-1].x_vm_cumulative >= budgets[0].x_vm_cumulative


def test_generate_initial_x0_with_cell_budgets_uses_budget_o2_not_exp_heuristic():
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = _prepare_reactor_for_budget(cfg, n_cells=3)
    T_profile, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )

    generate_initial_x0(
        cells=reactor.cells,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        T_profile=T_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
        cell_budgets=budgets,
    )

    idx = GAS_SPECIES_INDEX
    o2_top = float(reactor.cells[-1].N_d[idx["O2"]] + reactor.cells[-1].N_b[idx["O2"]])
    o2_bottom = float(reactor.cells[0].N_d[idx["O2"]] + reactor.cells[0].N_b[idx["O2"]])
    assert o2_top <= o2_bottom + 1e-9
    for cell in reactor.cells:
        assert np.all(np.isfinite(cell.N_d))
        assert np.all(cell.N_d >= 0.0)


def test_incremental_pyrolysis_major_mol_s_nonnegative():
    from src.solvers.vorabrechnung import VorabrechnungCellBudget, _incremental_pyrolysis_major_mol_s

    prev = VorabrechnungCellBudget(
        cell_index=0,
        axial_xi=0.1,
        tau_cumulative_s=1.0,
        t_budget_K=900.0,
        x_dry_cumulative=0.2,
        x_vm_cumulative=0.1,
        o2_remaining_mol_s=10.0,
        nC_vm_mol_s=0.5,
        nH_vm_mol_s=0.4,
        nO_vm_mol_s=0.1,
    )
    curr = VorabrechnungCellBudget(
        cell_index=1,
        axial_xi=0.3,
        tau_cumulative_s=2.0,
        t_budget_K=950.0,
        x_dry_cumulative=0.4,
        x_vm_cumulative=0.25,
        o2_remaining_mol_s=8.0,
        nC_vm_mol_s=1.0,
        nH_vm_mol_s=0.8,
        nO_vm_mol_s=0.2,
    )
    inc = _incremental_pyrolysis_major_mol_s(curr, prev, T_K=950.0, P_Pa=2.5e6)
    assert inc["CO"] > 0.0
    assert inc["O2"] == 0.0
    assert inc["N2"] == 0.0


def test_apply_vorabrechnung_bed_stream_mapping_anchors_bed0_primary_to_zu():
    from src.core.connectivity.gas_inlet import resolve_bottom_o2_dense_fraction

    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 3
    cfg.gas_inlet_split_strategy = "precalc_hydrodynamic_flux"
    reactor = _prepare_reactor_for_budget(cfg, n_cells=3)
    reactor._snapshot_bottom_gas_inlet_split_from_vorabrechnung()
    T_profile, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    generate_initial_x0(
        reactor.cells,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        T_profile=T_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
        cell_budgets=budgets,
    )
    idx = GAS_SPECIES_INDEX
    c0 = reactor.cells[0]
    o2_before = float(c0.N_d[idx["O2"]] + c0.N_b[idx["O2"]])
    assert o2_before < float(cfg.O2_feed) * 0.5

    diag = apply_vorabrechnung_bed_stream_mapping(reactor.cells, budgets, cfg)
    assert diag["bed_stream_cells_mapped"] == pytest.approx(1.0)
    o2_dense = float(resolve_bottom_o2_dense_fraction(c0, cfg))
    for sp in ("O2", "H2O", "N2"):
        j = idx[sp]
        total = float(max(c0.N_zu_d[j] + c0.N_zu_b[j], 0.0))
        assert float(c0.N_d[j] + c0.N_b[j]) == pytest.approx(total, rel=1e-9, abs=1e-12)
    j_o2 = idx["O2"]
    total_o2 = float(max(c0.N_zu_d[j_o2] + c0.N_zu_b[j_o2], 0.0))
    assert float(c0.N_d[j_o2]) == pytest.approx(total_o2 * o2_dense, rel=1e-6, abs=1e-9)
    j_h2o = idx["H2O"]
    assert float(c0.N_d[j_h2o]) != pytest.approx(float(c0.N_zu_d[j_h2o]), rel=1e-3)


def test_apply_vorabrechnung_bed_cell_mapping_respects_inlet_split_on_bed0():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 3
    cfg.gas_inlet_split_strategy = "precalc_hydrodynamic_flux"
    reactor = _prepare_reactor_for_budget(cfg, n_cells=3)
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
    reactor._snapshot_bottom_gas_inlet_split_from_vorabrechnung()
    T_profile, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    generate_initial_x0(
        reactor.cells,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        T_profile=T_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
        cell_budgets=budgets,
    )
    diag = apply_vorabrechnung_bed_cell_mapping(reactor.cells, budgets, cfg)
    assert diag["bed_cells_mapped"] == pytest.approx(3.0)
    c0 = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    from src.core.connectivity import resolve_bottom_gas_inlet_dense_fraction

    dense_frac = resolve_bottom_gas_inlet_dense_fraction(c0, cfg)
    nd = float(c0.N_d[idx["N2"]])
    nb = float(c0.N_b[idx["N2"]])
    assert nd / max(nd + nb, 1e-12) == pytest.approx(dense_frac, rel=1e-3, abs=1e-3)


def test_seed_major_product_holdup_from_vorab_gas_balance_reduces_bed0_product_residuals():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_a_tier_bed_stream_upstream_closure_thesis = False
    cfg.nr_bottom_joint_preprojection_thesis = False
    cfg.nr_bottom_major_phase_split_preprojection_thesis = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    diag = getattr(reactor, "_vorab_a_tier_product_holdup_seed_diag", None)
    assert diag is not None
    assert float(diag["bed_cells_targeted"]) >= 1.0
    assert float(diag["closure_total_mol_s"]) >= 0.0


def test_apply_vorabrechnung_bed_solid_reactive_profile_uses_budget_vm_dry():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 3
    reactor = _prepare_reactor_for_budget(cfg, n_cells=3)
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    from src.core.connectivity import set_bottom_cell_feeds
    from src.solvers.vorabrechnung import apply_vorabrechnung_bed_solid_reactive_profile_from_budget

    set_bottom_cell_feeds(reactor.cells, cfg)
    diag = apply_vorabrechnung_bed_solid_reactive_profile_from_budget(
        reactor.cells,
        budgets,
        cfg,
    )
    assert diag["bed_solid_reactive_cells_profiled"] >= 1.0
    c0 = reactor.cells[0]
    vm_total = float(np.sum(c0.m_solid[:, S_VM]))
    assert vm_total >= 0.0


def test_apply_vorabrechnung_bed_axial_primary_profile_sets_o2_from_budget():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 4
    cfg.vorab_a_tier_bed_axial_primary_profile_thesis = True
    cfg.vorab_a_tier_bed_axial_primary_blend_thesis = 1.0
    reactor = _prepare_reactor_for_budget(cfg, n_cells=4)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    from src.core.connectivity import set_bottom_cell_feeds, propagate_upstream
    from src.solvers.vorabrechnung import (
        apply_vorabrechnung_bed_axial_primary_profile_from_budget,
        apply_vorabrechnung_bed_stream_mapping,
    )

    set_bottom_cell_feeds(reactor.cells, cfg)
    apply_vorabrechnung_bed_stream_mapping(reactor.cells, budgets, cfg)
    for i in range(1, len(reactor.cells)):
        propagate_upstream(reactor.cells, cfg, i)
    bed1 = reactor.cells[1]
    idx = GAS_SPECIES_INDEX
    o2_before = float(bed1.N_d[idx["O2"]] + bed1.N_b[idx["O2"]])
    diag = apply_vorabrechnung_bed_axial_primary_profile_from_budget(reactor.cells, budgets, cfg)
    assert diag["bed_axial_primary_cells_profiled"] >= 1.0
    o2_after = float(bed1.N_d[idx["O2"]] + bed1.N_b[idx["O2"]])
    assert abs(o2_after - float(budgets[1].o2_remaining_mol_s)) <= 1e-6
    assert o2_after <= o2_before + 1e-9 or budgets[1].o2_remaining_mol_s < budgets[0].o2_remaining_mol_s


def test_stream_mapping_clears_inflated_n2_holdup_on_non_injection_beds():
    """分级末端后床格不得保留抬高的 N2 x0；H2O 在 thesis 旗标下同样透传（§5.27–5.28）。"""
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 4
    cfg.vorab_a_tier_bed_axial_primary_profile_thesis = True
    cfg.vorab_a_tier_post_staged_h2o_passthrough_thesis = True
    reactor = _prepare_reactor_for_budget(cfg, n_cells=4)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    from src.core.connectivity import set_bottom_cell_feeds, propagate_upstream
    from src.solvers.vorabrechnung import apply_vorabrechnung_bed_stream_mapping

    set_bottom_cell_feeds(reactor.cells, cfg)
    idx = GAS_SPECIES_INDEX
    j_n2 = idx["N2"]
    j_h2o = idx["H2O"]
    for bi in (2, 3):
        reactor.cells[bi].N_d[j_n2] = 80.0
        reactor.cells[bi].N_b[j_n2] = 20.0
        reactor.cells[bi].N_d[j_h2o] = 40.0
        reactor.cells[bi].N_b[j_h2o] = 10.0
        reactor.cells[bi]._thermo_cache_valid = False

    apply_vorabrechnung_bed_stream_mapping(reactor.cells, budgets, cfg)
    for i in range(1, len(reactor.cells)):
        propagate_upstream(reactor.cells, cfg, i)

    for bi in (1, 2, 3):
        cell = reactor.cells[bi]
        for j in (j_n2, j_h2o):
            n_hold = float(cell.N_d[j] + cell.N_b[j])
            n_in = float(cell.N_d_in[j] + cell.N_b_in[j])
            n_zu = float(cell.N_zu_d[j] + cell.N_zu_b[j])
            assert n_hold == pytest.approx(n_in + n_zu, rel=1e-6, abs=1e-9)
        assert float(cell.N_d[j_n2] + cell.N_b[j_n2]) < 90.0
        assert float(cell.N_d[j_h2o] + cell.N_b[j_h2o]) < 45.0


def test_apply_vorabrechnung_bed_axial_primary_profile_blends_o2_partially():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 4
    cfg.vorab_a_tier_bed_axial_primary_profile_thesis = True
    cfg.vorab_a_tier_bed_axial_primary_blend_thesis = 0.35
    reactor = _prepare_reactor_for_budget(cfg, n_cells=4)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    from src.core.connectivity import set_bottom_cell_feeds
    from src.solvers.vorabrechnung import (
        apply_vorabrechnung_bed_axial_primary_profile_from_budget,
        apply_vorabrechnung_bed_stream_mapping,
    )

    set_bottom_cell_feeds(reactor.cells, cfg)
    apply_vorabrechnung_bed_stream_mapping(reactor.cells, budgets, cfg)
    bed1 = reactor.cells[1]
    idx = GAS_SPECIES_INDEX
    o2_before = float(bed1.N_d[idx["O2"]] + bed1.N_b[idx["O2"]])
    o2_target = float(budgets[1].o2_remaining_mol_s)
    diag = apply_vorabrechnung_bed_axial_primary_profile_from_budget(reactor.cells, budgets, cfg)
    o2_after = float(bed1.N_d[idx["O2"]] + bed1.N_b[idx["O2"]])
    expected = o2_before + 0.35 * (o2_target - o2_before)
    assert diag["bed_axial_primary_blend_factor"] == 0.35
    assert abs(o2_after - expected) <= 1e-6
    if abs(o2_target - o2_before) > 1e-9:
        assert abs(o2_after - o2_target) > 1e-6
        assert abs(o2_after - o2_before) > 1e-9


def test_init_vorab_cell_abgleich_bridge_runs_on_phase2_case():
    from src.solvers.vorab_cell_abgleich import assess_vorabrechnung_cell_alignment

    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_a_tier_bed_stream_upstream_closure_thesis = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    abgleich_diag = getattr(reactor, "_vorab_a_tier_init_abgleich_diag", None)
    assert abgleich_diag is not None
    assert "aligned" in abgleich_diag
    assert "final" in abgleich_diag
    assert float(abgleich_diag["bridge_passes"]) >= 0.0
    final = abgleich_diag["final"]
    assert float(final["bed0_primary_holdup_rel_max"]) <= 0.06


def test_bed01_phase_split_respects_max_bed_index():
    from src.core.connectivity import set_bottom_cell_feeds
    from src.solvers.vorabrechnung import (
        apply_vorabrechnung_bed01_phase_split_after_vorab_closure,
        apply_vorabrechnung_bed_stream_mapping,
    )

    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 4
    cfg.vorab_a_tier_bed01_phase_split_after_vorab_closure_thesis = True
    cfg.vorab_a_tier_bed01_phase_split_max_bed_index_thesis = 0
    reactor = _prepare_reactor_for_budget(cfg, n_cells=4)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    set_bottom_cell_feeds(reactor.cells, cfg)
    apply_vorabrechnung_bed_stream_mapping(reactor.cells, budgets, cfg)
    diag = apply_vorabrechnung_bed01_phase_split_after_vorab_closure(reactor.cells, cfg)
    # max_bed_index=0 → 只碰 bed0；accepted 可为 0（stream mapping 已交换闭合）
    assert diag["bed01_phase_split_cells"] == 1.0
    assert diag["bed01_phase_split_accepted"] in (0.0, 1.0)


def test_assess_vorabrechnung_cell_alignment_detects_bed0_primary_drift():
    from src.solvers.vorab_cell_abgleich import assess_vorabrechnung_cell_alignment

    cfg = build_phase1_htw_lu_reactor_config()
    reactor = _prepare_reactor_for_budget(cfg, n_cells=3)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    c0 = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    c0.N_d[idx["O2"]] = 1.0e-6
    c0.N_b[idx["O2"]] = 1.0e-6
    metrics = assess_vorabrechnung_cell_alignment(reactor.cells, budgets, cfg)
    assert metrics["aligned"] is False
    assert float(metrics["bed0_primary_holdup_rel_max"]) > 0.5


def test_vorab_a_tier_cell_mapping_init_path_skips_bottom_joint():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_a_tier_cell_mapping_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    mapping_diag = getattr(reactor, "_vorab_a_tier_cell_mapping_diag", None)
    assert mapping_diag is not None
    assert mapping_diag.get("final_primary_inlet_split_align") is True
    joint_diag = getattr(reactor, "_bottom_joint_preprojection_run_diag", None)
    assert joint_diag is not None
    assert joint_diag.get("skipped_for_vorab_cell_mapping") is True


def test_assess_vorabrechnung_cell_alignment_includes_solid_reactive_metric():
    from src.solvers.vorab_cell_abgleich import assess_vorabrechnung_cell_alignment

    cfg = build_phase1_htw_lu_reactor_config()
    cfg.vorab_a_tier_bed_solid_reactive_profile_thesis = True
    reactor = _prepare_reactor_for_budget(cfg, n_cells=3)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
    )
    from src.solvers.vorabrechnung import apply_vorabrechnung_bed_solid_reactive_profile_from_budget

    apply_vorabrechnung_bed_solid_reactive_profile_from_budget(
        reactor.cells,
        budgets,
        cfg,
    )
    metrics = assess_vorabrechnung_cell_alignment(reactor.cells, budgets, cfg)
    assert "solid_reactive_budget_rel_max" in metrics
    assert float(metrics["solid_reactive_budget_rel_max"]) <= 0.05


def test_init_precalc_stores_frozen_vorab_cell_budgets():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_a_tier_bed_stream_upstream_closure_thesis = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    budgets = getattr(reactor, "_vorab_a_tier_cell_budgets", None)
    assert budgets is not None
    assert len(budgets) >= reactor.config.n_cells


def test_march_vorabrechnung_bed_temperature_starts_near_inlet():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 6
    cfg.vorab_bed_temperature_march_thesis = True
    cfg.vorab_staged_primary_gas_injection_thesis = False
    reactor = _prepare_reactor_for_budget(cfg, n_cells=6)
    T_profile, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
        cfg=cfg,
    )
    legacy_bottom = float(
        estimate_axial_T_profile(
            n_cells=6,
            T_inlet=cfg.T_inlet,
            O2_feed=cfg.O2_feed,
            H2O_feed=cfg.H2O_feed,
            N2_feed=cfg.N2_feed,
            fuel_feed_kg_s=cfg.fuel_feed,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            moisture_wt=cfg.moisture_wt,
            P=cfg.P,
        )[0]
    )
    assert float(T_profile[0]) <= legacy_bottom - 100.0
    assert float(T_profile[0]) <= float(cfg.T_inlet) + 250.0
    assert float(budgets[0].t_budget_K) == pytest.approx(float(T_profile[0]))
    if len(T_profile) >= 3:
        assert float(T_profile[2]) >= float(T_profile[0]) - 1.0


def test_staged_primary_gas_feeds_conserves_primary_moles():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 4
    cfg.vorab_staged_primary_gas_injection_thesis = True
    cfg.vorab_staged_primary_gas_bed_cells_thesis = 3
    reactor = Reactor(cfg)
    from src.core.connectivity import distribute_staged_primary_gas_to_bed_cells

    diag = distribute_staged_primary_gas_to_bed_cells(reactor.cells, cfg)
    assert diag["staged_primary_gas_bed_cells"] == pytest.approx(3.0)
    idx = GAS_SPECIES_INDEX
    total_o2 = sum(
        float(reactor.cells[i].N_zu_d[idx["O2"]] + reactor.cells[i].N_zu_b[idx["O2"]])
        for i in range(3)
    )
    total_h2o = sum(
        float(reactor.cells[i].N_zu_d[idx["H2O"]] + reactor.cells[i].N_zu_b[idx["H2O"]])
        for i in range(3)
    )
    assert total_o2 == pytest.approx(float(cfg.O2_feed))
    assert total_h2o == pytest.approx(float(cfg.H2O_feed))
    assert float(reactor.cells[0].N_zu_d[idx["O2"]] + reactor.cells[0].N_zu_b[idx["O2"]]) > float(
        reactor.cells[2].N_zu_d[idx["O2"]] + reactor.cells[2].N_zu_b[idx["O2"]]
    )


def test_staged_primary_gas_custom_layer_weights():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 4
    cfg.vorab_staged_primary_gas_injection_thesis = True
    cfg.vorab_staged_primary_gas_bed_cells_thesis = 3
    cfg.vorab_staged_primary_gas_layer_weights_thesis = (0.04, 0.78, 0.17)
    reactor = Reactor(cfg)
    from src.core.connectivity import distribute_staged_primary_gas_to_bed_cells

    diag = distribute_staged_primary_gas_to_bed_cells(reactor.cells, cfg)
    assert diag["staged_primary_gas_layer_weights"] == pytest.approx(
        [0.04 / 0.99, 0.78 / 0.99, 0.17 / 0.99]
    )
    idx = GAS_SPECIES_INDEX
    o2 = [
        float(reactor.cells[i].N_zu_d[idx["O2"]] + reactor.cells[i].N_zu_b[idx["O2"]])
        for i in range(3)
    ]
    assert sum(o2) == pytest.approx(float(cfg.O2_feed))
    assert o2[0] == pytest.approx(0.04 / 0.99 * float(cfg.O2_feed), rel=1e-6)
    assert o2[1] == pytest.approx(0.78 / 0.99 * float(cfg.O2_feed), rel=1e-6)
    assert o2[2] == pytest.approx(0.17 / 0.99 * float(cfg.O2_feed), rel=1e-6)


def test_march_vorabrechnung_widens_nr_temperature_fence_upper_bound():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_bed_temperature_march_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    bed0 = reactor.cells[0]
    marched = getattr(reactor, "_vorab_marched_t_budget_K", None)
    assert marched is not None
    # Cold-budget march keeps bed0 budget near inlet; NR x₀ may cold-anchor at
    # seed_bottom, but the fence upper bound must still open to the march anchor.
    assert float(marched[0]) <= float(cfg.T_inlet) + 250.0
    assert float(getattr(bed0, "_nr_temperature_max_K", 0.0)) >= float(
        getattr(cfg, "vorab_nr_temperature_fence_upper_anchor_K_thesis", 1400.0)
    )
    assert float(bed0.T) == pytest.approx(
        float(cfg.vorab_nr_temperature_seed_bottom_K_thesis), abs=1.0
    )
    assert float(bed0._nr_temperature_max_K) > float(bed0.T) + 500.0


def test_build_nr_temperature_seed_adiabatic_anchor_matches_legacy_clip():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 6
    cfg.vorab_bed_temperature_march_thesis = True
    cfg.vorab_nr_temperature_seed_mode_thesis = "adiabatic_anchor"
    reactor = _prepare_reactor_for_budget(cfg, n_cells=6)
    T_marched, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
        cfg=cfg,
    )
    from src.solvers.vorabrechnung import build_nr_temperature_seed_for_march_mode

    T_nr = build_nr_temperature_seed_for_march_mode(
        reactor.cells,
        budgets,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
        cfg=cfg,
    )
    T_legacy, _ = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
        T_bottom_cap_K=float(cfg.vorab_bottom_temperature_cap_K_thesis),
        cfg=None,
    )
    assert float(T_marched[0]) <= float(cfg.T_inlet) + 250.0
    assert float(T_nr[0]) == pytest.approx(float(T_legacy[0]), rel=1e-6)
    assert float(T_nr[0]) >= 1000.0


def test_derive_nr_temperature_seed_decoupled_from_marched_budget():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 6
    cfg.vorab_bed_temperature_march_thesis = True
    cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
    cfg.vorab_staged_primary_gas_injection_thesis = False
    reactor = _prepare_reactor_for_budget(cfg, n_cells=6)
    T_marched, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
        cfg=cfg,
    )
    T_nr = derive_nr_temperature_seed_from_vorab_budgets(
        reactor.cells,
        budgets,
        O2_feed=cfg.O2_feed,
        cfg=cfg,
    )
    assert float(T_marched[0]) <= float(cfg.T_inlet) + 250.0
    assert float(budgets[0].t_budget_K) <= float(cfg.T_inlet) + 250.0
    assert float(T_nr[0]) >= float(cfg.vorab_nr_temperature_seed_bottom_K_thesis) - 1.0
    assert float(T_nr[-1]) >= float(T_nr[0])
    assert float(T_nr[0]) > float(T_marched[0]) + 200.0


def test_cold_anchor_pins_bed0_only_upper_beds_follow_o2_progress():
    """bed0 cold-anchor must not replace upper-bed O₂/ξ progress with index leap."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 6
    cfg.vorab_bed_temperature_march_thesis = True
    cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
    cfg.vorab_nr_temperature_seed_bed0_cold_anchor_thesis = True
    cfg.vorab_staged_primary_gas_injection_thesis = False
    reactor = _prepare_reactor_for_budget(cfg, n_cells=6)
    _, budgets = build_vorabrechnung_cell_budgets(
        reactor.cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        fuel_type=cfg.fuel_type,
        P=cfg.P,
        cfg=cfg,
    )
    T_anchor = derive_nr_temperature_seed_from_vorab_budgets(
        reactor.cells, budgets, O2_feed=cfg.O2_feed, cfg=cfg,
    )
    cfg.vorab_nr_temperature_seed_bed0_cold_anchor_thesis = False
    T_free = derive_nr_temperature_seed_from_vorab_budgets(
        reactor.cells, budgets, O2_feed=cfg.O2_feed, cfg=cfg,
    )
    T_bottom = float(cfg.vorab_nr_temperature_seed_bottom_K_thesis)
    assert float(T_anchor[0]) == pytest.approx(T_bottom, abs=1.0)
    # Upper beds match the O₂/ξ profile (same as no-anchor), not an index leap.
    assert float(T_anchor[1]) == pytest.approx(float(T_free[1]), rel=1e-6)
    assert float(T_anchor[2]) == pytest.approx(float(T_free[2]), rel=1e-6)
