"""Vorab VM 脱挥发区分区释放。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import S_VM
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.vorab_cell_abgleich import assess_vorabrechnung_cell_alignment
from src.solvers.vorabrechnung import (
    _incremental_pyrolysis_major_mol_s,
    _major_product_init_gas_source_mol_s,
    _zero_vorab_budget,
    build_vorabrechnung_cell_budgets,
    seed_pyrolysis_gas_dense_holdup_for_nr_x0,
)
from src.thermal.devolatilization import (
    vm_devolatilization_zone_cumulative_fraction,
    vm_devolatilization_zone_increment_fraction,
)
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


def _prepare_reactor(cfg, *, n_cells: int = 6) -> Reactor:
    cfg.n_cells = n_cells
    reactor = Reactor(cfg)
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    for cell in reactor.cells:
        cell.calc_hydrodynamics()
    return reactor


def test_vm_zone_fractions_sum_to_total():
    x_total = 0.95
    n_zone = 3
    cum = [vm_devolatilization_zone_cumulative_fraction(j, n_zone, x_total) for j in range(n_zone)]
    assert cum[0] < cum[1] < cum[2]
    assert abs(cum[-1] - x_total) <= 1e-12
    inc = [
        vm_devolatilization_zone_increment_fraction(j, n_zone, x_total)
        for j in range(n_zone)
    ]
    assert all(0.0 < v <= 1.0 for v in inc)


def test_vm_zone_front_loaded_puts_more_in_bed0_than_linear():
    x_total = 1.0
    n_zone = 3
    lin0 = vm_devolatilization_zone_cumulative_fraction(0, n_zone, x_total, profile="linear")
    front0 = vm_devolatilization_zone_cumulative_fraction(
        0, n_zone, x_total, profile="front_loaded", front_power=0.5
    )
    assert front0 > lin0
    # ζ=1/3 → √(1/3)≈0.577
    assert abs(front0 - (1.0 / 3.0) ** 0.5) <= 1e-12
    front_end = vm_devolatilization_zone_cumulative_fraction(
        n_zone - 1, n_zone, x_total, profile="front_loaded"
    )
    assert abs(front_end - x_total) <= 1e-12
    # 绝对释放量：bed0 前重
    d_lin = [
        vm_devolatilization_zone_cumulative_fraction(j, n_zone, x_total, profile="linear")
        - (
            0.0
            if j == 0
            else vm_devolatilization_zone_cumulative_fraction(j - 1, n_zone, x_total, profile="linear")
        )
        for j in range(n_zone)
    ]
    d_front = [
        vm_devolatilization_zone_cumulative_fraction(j, n_zone, x_total, profile="front_loaded")
        - (
            0.0
            if j == 0
            else vm_devolatilization_zone_cumulative_fraction(
                j - 1, n_zone, x_total, profile="front_loaded"
            )
        )
        for j in range(n_zone)
    ]
    assert d_front[0] > d_lin[0]
    assert abs(sum(d_front) - x_total) <= 1e-12


def test_build_vorab_budgets_spreads_vm_over_three_bed_cells():
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=False)
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 3
    reactor = _prepare_reactor(cfg, n_cells=6)

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

    x0 = float(budgets[0].x_vm_cumulative)
    x1 = float(budgets[1].x_vm_cumulative)
    x2 = float(budgets[2].x_vm_cumulative)
    x_top = float(budgets[2].x_vm_cumulative)
    assert x0 < 0.99
    assert x1 > x0
    assert abs(x2 - x_top) <= 1e-9
    assert abs(x2 - 1.0) <= 0.05

    prev = _zero_vorab_budget()
    inc_co = []
    for i in range(3):
        inc = _incremental_pyrolysis_major_mol_s(
            budgets[i],
            prev,
            T_K=float(budgets[i].t_budget_K),
            P_Pa=float(cfg.P),
            pyrolysis_tar_carbon_frac=float(reactor.cells[0].solid.pyrolysis_tar_carbon_frac),
        )
        inc_co.append(float(inc["CO"]))
        prev = budgets[i]

    assert inc_co[0] > 0.0
    assert inc_co[1] > 0.0
    assert inc_co[2] > 0.0
    # linear VM 原子增量近似均分；Hegermann x_C,Koks/Teer 随 T 变化，CO 不必严格等分
    total_co = sum(inc_co)
    assert total_co > 0.0
    shares = [c / total_co for c in inc_co]
    assert all(s > 0.15 for s in shares)
    assert max(shares) < 0.55


def test_vm_zone_solid_reactive_holdup_matches_axial_budget():
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=False)
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 3
    cfg.n_cells = 6
    reactor = _prepare_reactor(cfg, n_cells=6)

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
    from src.solvers.vorabrechnung import (
        _expected_vm_holdup_kg_vm_zone,
        _fuel_vm_mass_kg_from_bed0,
        align_vm_devolatilization_zone_vorab_state,
    )

    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    align_vm_devolatilization_zone_vorab_state(
        reactor.cells,
        budgets,
        cfg,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        reseed_holdup_chain_fn=reactor._seed_initialized_holdup_chain_for_nr,
        refresh_zone_pyrolysis_sources=True,
    )
    total_vm = _fuel_vm_mass_kg_from_bed0(reactor.cells)
    bed_indices = [i for i, c in enumerate(reactor.cells) if c.cell_type == "bed"]
    for rank in range(2):
        bi = bed_indices[rank]
        prev_x = 0.0 if rank == 0 else float(budgets[bed_indices[rank - 1]].x_vm_cumulative)
        exp = _expected_vm_holdup_kg_vm_zone(
            cell_rank=rank,
            total_vm_kg=total_vm,
            prev_x_vm_cumulative=prev_x,
            x_vm_cumulative=float(budgets[bi].x_vm_cumulative),
        )
        act = float(np.sum(reactor.cells[bi].m_solid[:, S_VM]))
        assert abs(act - exp) <= max(0.05 * max(exp, 1e-9), 1e-5)

    from tests.validation_case_utils import build_phase1_htw_lu_reactor_config

    cfg = build_phase1_htw_lu_reactor_config()
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 0
    reactor = _prepare_reactor(cfg, n_cells=4)

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

    prev = _zero_vorab_budget()
    inc0 = _incremental_pyrolysis_major_mol_s(
        budgets[0],
        prev,
        T_K=float(budgets[0].t_budget_K),
        P_Pa=float(cfg.P),
        pyrolysis_tar_carbon_frac=float(reactor.cells[0].solid.pyrolysis_tar_carbon_frac),
    )
    inc1 = _incremental_pyrolysis_major_mol_s(
        budgets[1],
        budgets[0],
        T_K=float(budgets[1].t_budget_K),
        P_Pa=float(cfg.P),
        pyrolysis_tar_carbon_frac=float(reactor.cells[0].solid.pyrolysis_tar_carbon_frac),
    )
    assert float(inc0["CO"]) > 0.0
    assert float(inc0["CO"]) > 10.0 * max(float(inc1["CO"]), 1.0e-9)


@pytest.mark.slow
def test_vm_zone_init_abgleich_product_holdup_uses_vorab_source_not_kinetics():
    """VM 分区 + init abgleich：产物 Check2 用 Vorab 固定源，避免 R_gas<0 时 rel 爆炸。"""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 3
    reactor = Reactor(cfg)

    real_seed = seed_pyrolysis_gas_dense_holdup_for_nr_x0
    pre_metrics: dict[str, float] = {}

    def capture_before_pyro(*args, **kwargs):
        budgets = reactor._vorab_a_tier_cell_budgets
        m = assess_vorabrechnung_cell_alignment(
            reactor.cells,
            budgets,
            cfg,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
        )
        pre_metrics["product_rel"] = float(m["product_source_holdup_rel_max"])
        c0 = reactor.cells[0]
        j = int(GAS_SPECIES_INDEX["CH4"])
        c0.residuals()
        pre_metrics["ch4_vm"] = _major_product_init_gas_source_mol_s(c0, j)
        pre_metrics["ch4_rgas"] = float(c0.R_gas_d[j] + c0.R_gas_b[j])
        return real_seed(*args, **kwargs)

    import src.solvers.vorabrechnung as vorab_mod

    vorab_mod.seed_pyrolysis_gas_dense_holdup_for_nr_x0 = capture_before_pyro
    try:
        run_init_and_precalc_for_global_nr(
            reactor,
            init_strategy="vorabrechnung",
            gs_warmup_steps=None,
        )
    finally:
        vorab_mod.seed_pyrolysis_gas_dense_holdup_for_nr_x0 = real_seed

    assert pre_metrics["ch4_vm"] > 0.0
    assert pre_metrics["product_rel"] < 10.0
    if pre_metrics["ch4_rgas"] < 0.0:
        assert pre_metrics["ch4_vm"] > 0.0

    abgleich = getattr(reactor, "_vorab_a_tier_init_abgleich_diag", None)
    assert abgleich is not None
    assert float(abgleich["history"][0]["product_source_holdup_rel_max"]) < 10.0
    assert float(abgleich["final"]["product_source_holdup_rel_max"]) < 1.0
    assert float(abgleich["final"]["solid_reactive_budget_rel_max"]) < 0.05
    # Phase2 cold-anchor 故意保留大 ΔT；本用例只验收产物/固相 Vorab 闭合，不要求满 aligned
    assert float(abgleich["final"]["product_source_holdup_rel_max"]) <= float(
        abgleich["final"]["tol_product_rel"]
    )
    assert float(abgleich["final"]["solid_reactive_budget_rel_max"]) <= float(
        abgleich["final"]["tol_solid_reactive_rel"]
    )
