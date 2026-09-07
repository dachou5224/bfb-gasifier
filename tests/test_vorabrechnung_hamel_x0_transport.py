"""Hamel-aligned transport x₀ 与 §4.3 热解 Gibbs 组成测试。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.vorabrechnung import (
    VorabrechnungCellBudget,
    _incremental_pyrolysis_major_mol_s,
    _resolve_product_bed_cell_indices,
    build_vorabrechnung_cell_budgets,
    generate_initial_x0,
    solve_pyrolysis_gas_gibbs_composition,
)
from tests.test_vorabrechnung_a_tier_budget import _prepare_reactor_for_budget
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


def test_resolve_product_bed_cell_limit_zero_means_all_beds():
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.n_cells = 6
    reactor = Reactor(cfg)
    indices = _resolve_product_bed_cell_indices(reactor.cells, 0)
    bed_count = sum(1 for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed")
    assert len(indices) == bed_count
    assert indices == list(range(bed_count))


def test_pyrolysis_gibbs_composition_finite_and_nonnegative():
    # 无 C₂H₄：需 H 足够（半原子口径 nH），否则 An=b 不可行
    comp = solve_pyrolysis_gas_gibbs_composition(
        T_K=1088.0,
        P_Pa=2.5e6,
        nC_mol=0.05,
        nH_mol=0.12,
        nO_mol=0.02,
        solver_mode="hamel_reduced",
    )
    for sp in ("CO", "CO2", "H2", "CH4", "H2O"):
        assert np.isfinite(comp[sp])
        assert comp[sp] >= 0.0
    assert comp["CO"] + comp["CO2"] + comp["CH4"] > 0.0


def test_pyrolysis_gibbs_at_local_T_without_temperature_floor():
    """Hamel §4.3：在给定局部 T 求解，无温度地板/升档；全量 An=b。"""
    # 可行池（已含足够 H）：O→CO，剩余 C/H 可被 CH4/H2 吃完
    comp = solve_pyrolysis_gas_gibbs_composition(
        T_K=700.0,
        P_Pa=2.5e6,
        nC_mol=0.01,
        nH_mol=0.02,
        nO_mol=0.005,
        solver_mode="hamel_reduced",
    )
    for sp in ("CO", "CO2", "H2", "CH4", "H2O"):
        assert np.isfinite(comp[sp])
        assert comp[sp] >= 0.0
    assert (
        comp["CO"] + comp["CO2"] + comp["CH4"] + comp["H2"] + comp["H2O"]
        > 0.0
    )


def test_incremental_pyrolysis_gibbs_differs_from_legacy_fixed_ratio():
    from src.solvers.vorabrechnung import _legacy_fixed_ratio_pyrolysis_increment_not_hamel

    prev = VorabrechnungCellBudget(
        cell_index=0,
        axial_xi=0.4,
        tau_cumulative_s=2.0,
        t_budget_K=1000.0,
        x_dry_cumulative=0.5,
        x_vm_cumulative=0.3,
        o2_remaining_mol_s=5.0,
        nC_vm_mol_s=0.4,
        nH_vm_mol_s=0.55,
        nO_vm_mol_s=0.08,
    )
    curr = VorabrechnungCellBudget(
        cell_index=1,
        axial_xi=0.55,
        tau_cumulative_s=3.0,
        t_budget_K=1088.0,
        x_dry_cumulative=0.7,
        x_vm_cumulative=0.45,
        o2_remaining_mol_s=3.0,
        nC_vm_mol_s=0.55,
        nH_vm_mol_s=0.75,
        nO_vm_mol_s=0.10,
    )
    dC = curr.nC_vm_mol_s - prev.nC_vm_mol_s
    dH = curr.nH_vm_mol_s - prev.nH_vm_mol_s
    dO = curr.nO_vm_mol_s - prev.nO_vm_mol_s
    heur = _legacy_fixed_ratio_pyrolysis_increment_not_hamel(dC=dC, dH=dH, dO=dO)
    gibbs = _incremental_pyrolysis_major_mol_s(
        curr,
        prev,
        T_K=1088.0,
        P_Pa=2.5e6,
        pyrolysis_gibbs_solver_mode="hamel_reduced",
    )
    assert heur["CO"] > 0.0
    assert gibbs["CO"] > 0.0
    ratio = gibbs["CO"] / max(heur["CO"], 1e-12)
    assert 0.05 < ratio < 20.0


def test_transport_x0_bed5_co_not_collapsed_vs_per_cell_gibbs():
    """transport x₀ 应沿 inflow 保持 CO，避免 ~1105 K Gibbs 分岔将 bed5 CO 清零。"""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.n_cells = 6
    cfg.n_freeboard_cells = 0
    idx = GAS_SPECIES_INDEX

    def _co_after_x0(x0_mode: str) -> tuple[float, float]:
        reactor = _prepare_reactor_for_budget(cfg, n_cells=6)
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
            cfg=cfg,
        )
        # 使用 Vorab 冷预算 t_budget_K（~1105 K 分岔场景），非 NR adiabatic T_seed
        T_est = np.array([float(b.t_budget_K) for b in budgets], dtype=np.float64)
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
            T_profile=T_est,
            fuel_type=cfg.fuel_type,
            use_hamel_major_gibbs_x0=(x0_mode == "major_gibbs_per_cell"),
            major_gibbs_solver_mode="hamel_reduced",
            cell_budgets=budgets,
            x0_holdup_mode=x0_mode,
            cfg=cfg,
            use_pyrolysis_gibbs_composition=True,
            pyrolysis_gibbs_solver_mode="hamel_reduced",
        )
        co2 = float(reactor.cells[2].N_d[idx["CO"]] + reactor.cells[2].N_b[idx["CO"]])
        co4 = float(reactor.cells[4].N_d[idx["CO"]] + reactor.cells[4].N_b[idx["CO"]])
        return co2, co4

    co2_gibbs, co4_gibbs = _co_after_x0("major_gibbs_per_cell")
    co2_transport, co4_transport = _co_after_x0("transport_from_vorab")

    assert co2_gibbs > 5.0
    assert co4_gibbs < 0.01 * co2_gibbs
    assert co4_transport > 0.5 * co2_transport
    assert co4_transport > 1.0


def test_fast_oxidation_stoich_consumes_h2_before_co_and_ch4():
    """Startwert R12→R5→R6：有限 O₂ 时优先消 H₂。"""
    from src.core.species import N_GAS
    from src.solvers.vorabrechnung.vorab_x0 import _apply_fast_oxidation_stoich_to_phase_holdup

    idx = GAS_SPECIES_INDEX
    N = np.zeros(N_GAS)
    N[idx["H2"]] = 2.0
    N[idx["O2"]] = 0.6  # 仅够 R12（需 1.0 才吃完 2 mol H₂）
    N[idx["CO"]] = 1.0
    N[idx["CH4"]] = 1.0
    diag = _apply_fast_oxidation_stoich_to_phase_holdup(N)
    assert diag["xi12"] == 1.2  # 0.6 O₂ → 1.2 H₂
    assert N[idx["H2"]] == 0.8
    assert N[idx["O2"]] < 1e-12
    assert diag["xi5"] == 0.0
    assert diag["xi6"] == 0.0
    assert diag["xi10"] == 0.0


def test_fast_oxidation_stoich_closes_tar_o2_via_r10():
    """Startwert 含 R10：消除 TAR+O₂ 共存；O₂ 充足时再吃掉产氢。"""
    from src.core.species import N_GAS
    from src.solvers.vorabrechnung.vorab_x0 import _apply_fast_oxidation_stoich_to_phase_holdup

    idx = GAS_SPECIES_INDEX
    N_lim = np.zeros(N_GAS)
    N_lim[idx["TAR1"]] = 0.4
    N_lim[idx["TAR2"]] = 0.6
    N_lim[idx["O2"]] = 10.0
    diag_lim = _apply_fast_oxidation_stoich_to_phase_holdup(N_lim, fuel_type="coal", n_passes=3)
    assert diag_lim["xi10"] > 0.0
    # O₂ 不足吃完 R10 产氢时，至少不应再与 O₂ 共存
    assert float(N_lim[idx["O2"]]) < 1e-9 or float(N_lim[idx["TAR1"]] + N_lim[idx["TAR2"]]) < 1e-9

    N = np.zeros(N_GAS)
    N[idx["TAR1"]] = 0.4
    N[idx["TAR2"]] = 0.6
    N[idx["O2"]] = 40.0
    diag = _apply_fast_oxidation_stoich_to_phase_holdup(N, fuel_type="coal", n_passes=3)
    assert diag["xi10"] > 0.0
    assert float(N[idx["TAR1"]] + N[idx["TAR2"]]) < 1e-9
    assert float(N[idx["H2"]]) < 1e-9
    assert float(N[idx["O2"]]) > 0.0  # 仍有剩余氧化剂


def test_project_bed_holdup_fast_oxidation_removes_h2_o2_coexistence():
    from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 1.0
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 2.0
    diag = project_bed_holdup_fast_oxidation_closure(reactor.cells, cfg, bed_cell_limit=1)
    assert diag["xi12_total"] >= 2.0 - 1e-9
    assert float(cell.N_d[idx["H2"]] + cell.N_b[idx["H2"]]) < 1e-12


def test_per_phase_fastox_keeps_bubble_o2_and_dense_h2():
    """分相投影不得把气泡氧拿去烧悬浮相氢。"""
    from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["H2"]] = 2.0
    cell.N_b[idx["O2"]] = 4.0
    cell.N_d[idx["N2"]] = 5.0
    cell.N_b[idx["N2"]] = 5.0
    diag = project_bed_holdup_fast_oxidation_closure(
        reactor.cells, cfg, bed_cell_limit=1, phase_update="per_phase"
    )
    assert diag["phase_update_per_phase"] == 1.0
    assert float(cell.N_d[idx["H2"]]) == pytest.approx(2.0)
    assert float(cell.N_b[idx["O2"]]) == pytest.approx(4.0)
    assert float(cell.N_d[idx["O2"]]) == pytest.approx(0.0)
    assert float(cell.N_b[idx["H2"]]) == pytest.approx(0.0)


def test_per_phase_fastox_burns_only_intra_phase_overlap():
    from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["H2"]] = 2.0
    cell.N_d[idx["O2"]] = 0.5  # 同相：2 H₂ + 0.5 O₂ → 烧掉 1 H₂，剩 1 H₂
    cell.N_b[idx["O2"]] = 4.0
    diag = project_bed_holdup_fast_oxidation_closure(
        reactor.cells, cfg, bed_cell_limit=1, phase_update="per_phase"
    )
    assert diag["xi12_total"] == pytest.approx(1.0)
    assert float(cell.N_d[idx["H2"]]) == pytest.approx(1.0)
    assert float(cell.N_d[idx["O2"]]) < 1e-12
    assert float(cell.N_b[idx["O2"]]) == pytest.approx(4.0)


def test_char_r1_stoich_consumes_leftover_o2_to_co2():
    """残氧 + 炭：α=1 时 C+O₂→CO₂。"""
    from src.core.species import N_GAS
    from src.kinetics.char_reactions import R4_M_C
    from src.solvers.vorabrechnung.vorab_x0 import _apply_char_r1_stoich_to_holdup

    idx = GAS_SPECIES_INDEX
    N = np.zeros(N_GAS)
    N[idx["O2"]] = 2.0
    m_char = np.array([1.0])  # [kg] ≫ 2 mol C
    diag = _apply_char_r1_stoich_to_holdup(N, m_char, alpha=1.0)
    assert diag["xi_c"] == pytest.approx(2.0)
    assert float(N[idx["O2"]]) < 1e-12
    assert float(N[idx["CO2"]]) == pytest.approx(2.0)
    assert float(N[idx["CO"]]) == pytest.approx(0.0)
    assert float(m_char[0]) == pytest.approx(1.0 - 2.0 * R4_M_C)


def test_project_fastox_eats_char_when_gas_fuel_gone():
    """气相燃料已空、仍有 O₂ 时，R1 投影须动炭并降氧。"""
    from src.core.cell import S_CHAR
    from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["O2"]] = 3.0
    cell.N_b[idx["N2"]] = 10.0
    cell.m_solid[:, S_CHAR] = 0.5
    o2_before = 3.0
    char_before = float(np.sum(cell.m_solid[:, S_CHAR]))
    diag = project_bed_holdup_fast_oxidation_closure(reactor.cells, cfg, bed_cell_limit=1)
    o2_after = float(cell.N_d[idx["O2"]] + cell.N_b[idx["O2"]])
    char_after = float(np.sum(cell.m_solid[:, S_CHAR]))
    assert float(diag["xi_r1_total"]) > 0.0
    assert o2_after < o2_before - 0.1
    assert char_after < char_before
    assert float(cell.N_d[idx["H2"]] + cell.N_b[idx["H2"]]) < 1e-12


def test_project_fastox_r1_flag_off_leaves_o2_when_gas_fuel_gone():
    from src.core.cell import S_CHAR
    from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["O2"]] = 3.0
    cell.N_b[idx["N2"]] = 10.0
    cell.m_solid[:, S_CHAR] = 0.5
    diag = project_bed_holdup_fast_oxidation_closure(reactor.cells, cfg, bed_cell_limit=1)
    assert float(diag.get("xi_r1_total", 0.0)) == 0.0
    assert float(cell.N_d[idx["O2"]] + cell.N_b[idx["O2"]]) == pytest.approx(3.0)


def test_bed0_eq24_two_phase_startwert_flag_off_is_noop():
    from src.solvers.vorabrechnung.vorab_x0 import apply_bed0_eq24_two_phase_startwert

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = False
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["O2"]] = 1.0
    cell.N_b[idx["O2"]] = 1.0
    diag = apply_bed0_eq24_two_phase_startwert(reactor.cells, cfg)
    assert diag["enabled"] == 0.0
    assert float(cell.N_d[idx["O2"]]) == pytest.approx(1.0)
    assert float(cell.N_b[idx["O2"]]) == pytest.approx(1.0)


def test_bed0_eq24_two_phase_puts_o2_in_bubble_and_syngas_in_dense():
    """O₂ 与合成气分相：同相不共存，总量按 zu+R 闭合。"""
    from src.solvers.vorabrechnung.vorab_x0 import apply_bed0_eq24_two_phase_startwert

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_zu_d[idx["O2"]] = 4.0
    cell.N_zu_b[idx["O2"]] = 4.0
    cell.N_d[idx["O2"]] = 1.0
    cell.N_b[idx["O2"]] = 1.0
    cell.N_d[idx["N2"]] = 5.0
    cell.N_b[idx["N2"]] = 5.0
    cell.N_zu_d[idx["N2"]] = 5.0
    cell.N_zu_b[idx["N2"]] = 5.0
    cell._vm_gas_source_cache[idx["H2"]] = 2.0
    cell._vm_gas_source_cache[idx["CO"]] = 1.5
    cell._vm_cache_valid = True
    cell.T = 1150.0
    cell.P = 2.5e6
    cell.calc_hydrodynamics()
    diag = apply_bed0_eq24_two_phase_startwert(reactor.cells, cfg, n_passes=1)
    assert diag["enabled"] == 1.0
    n_o2_d = float(cell.N_d[idx["O2"]])
    n_o2_b = float(cell.N_b[idx["O2"]])
    n_h2_d = float(cell.N_d[idx["H2"]])
    n_h2_b = float(cell.N_b[idx["H2"]])
    assert n_o2_b > n_o2_d
    assert n_h2_d > n_h2_b
    y_d = float(cell._mole_fractions("d")[idx["O2"]])
    y_b = float(cell._mole_fractions("b")[idx["O2"]])
    assert y_b > y_d + 1e-6
    overlap = n_o2_d * n_h2_d
    assert overlap < 1e-6


def test_bed0_r1_oxidizer_seed_moves_bubble_o2_and_leaves_dense_leftover():
    """stoich+seed 从气泡转入悬浮相；按相 fastox 后留下 seed 残氧。"""
    from src.solvers.vorabrechnung.vorab_x0 import apply_bed0_r1_oxidizer_seed

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    cfg.vorab_bed0_r1_oxidizer_seed_mol_s_thesis = 1.0
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["H2"]] = 2.0  # [mol/s]
    cell.N_d[idx["CO"]] = 2.0  # [mol/s]
    cell.N_b[idx["O2"]] = 8.0  # [mol/s]
    cell.N_d[idx["N2"]] = 5.0
    cell.N_b[idx["N2"]] = 5.0
    diag = apply_bed0_r1_oxidizer_seed(reactor.cells, cfg)
    assert diag["enabled"] == 1.0
    assert diag["moved_o2"] == pytest.approx(3.0)
    assert float(cell.N_d[idx["O2"]]) == pytest.approx(1.0, abs=1e-9)
    assert float(cell.N_b[idx["O2"]]) == pytest.approx(5.0, abs=1e-9)
    assert float(cell.N_d[idx["H2"]]) < 1e-12
    assert float(cell.N_d[idx["CO"]]) < 1e-12

    cfg.vorab_bed0_r1_oxidizer_seed_mol_s_thesis = 0.0
    cell.N_d[idx["O2"]] = 0.0
    diag_off = apply_bed0_r1_oxidizer_seed(reactor.cells, cfg)
    assert diag_off["enabled"] == 0.0


def test_bed0_inert_exchange_split_flag_off_is_noop():
    from src.solvers.vorabrechnung.vorab_x0 import apply_bed0_inert_exchange_split_startwert

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = False
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["N2"]] = 3.0
    cell.N_b[idx["N2"]] = 9.0
    diag = apply_bed0_inert_exchange_split_startwert(reactor.cells, cfg)
    assert diag["enabled"] == 0.0
    assert float(cell.N_d[idx["N2"]]) == pytest.approx(3.0)
    assert float(cell.N_b[idx["N2"]]) == pytest.approx(9.0)


def test_bed0_inert_exchange_split_closes_n2_and_keeps_r1_seed():
    """惰性重分贴 N₂ Eq.2.2；不挪悬浮相残氧。"""
    from src.solvers.vorabrechnung.vorab_x0 import apply_bed0_inert_exchange_split_startwert

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_zu_d[idx["N2"]] = 10.0
    cell.N_zu_b[idx["N2"]] = 10.0
    cell.N_d[idx["N2"]] = 21.5
    cell.N_b[idx["N2"]] = 20.8
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 3.0
    cell.N_d[idx["H2"]] = 0.0
    cell.N_b[idx["H2"]] = 0.0
    cell.T = 1150.0  # [K]
    cell.P = 2.5e6  # [Pa]
    cell.calc_hydrodynamics()
    cell.residuals()
    n2_res_before = float(
        cell.N_zu_d[idx["N2"]]
        + cell.N_d_in[idx["N2"]]
        + cell.N_rez_d[idx["N2"]]
        + cell.R_gas_d[idx["N2"]]
        - cell.N_d[idx["N2"]]
        + cell.N_ex[idx["N2"]]
    )
    diag = apply_bed0_inert_exchange_split_startwert(reactor.cells, cfg, n_passes=3)
    assert diag["enabled"] == 1.0
    assert float(cell.N_d[idx["O2"]]) == pytest.approx(2.0, abs=1e-12)
    assert float(cell.N_b[idx["O2"]]) == pytest.approx(3.0, abs=1e-12)
    n2_res_after = float(
        cell.N_zu_d[idx["N2"]]
        + cell.N_d_in[idx["N2"]]
        + cell.N_rez_d[idx["N2"]]
        + cell.R_gas_d[idx["N2"]]
        - cell.N_d[idx["N2"]]
        + cell.N_ex[idx["N2"]]
    )
    assert abs(n2_res_after) < abs(n2_res_before)
    assert abs(n2_res_after) < 8.0


def test_init_realigns_upper_holdup_after_fastox_startwert():
    """fastox / R1 种子后 bed1 的 char/ash 须贴 Eq.2.6 support（局部 snap，handoff §5.55）。"""
    from src.core.cell import S_ASH, S_CHAR
    from src.core.connectivity.solid_transport import _transport_inflow_support_holdup
    from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = True
    cfg.vorab_bed0_r1_oxidizer_seed_mol_s_thesis = 2.0
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(reactor, init_strategy="vorabrechnung", gs_warmup_steps=None)
    cell = reactor.cells[1]
    assert str(getattr(cell, "cell_type", "bed")) == "bed"
    support = _transport_inflow_support_holdup(cell)
    for comp in (S_CHAR, S_ASH):
        hold_c = float(np.sum(np.maximum(cell.m_solid[:, comp], 0.0)))
        sup_c = float(np.sum(support[:, comp]))
        if sup_c < 1e-6:
            continue
        assert hold_c == pytest.approx(sup_c, rel=2e-3, abs=1e-3)
    bed0 = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    bed0.residuals()
    n2_res_d = float(
        bed0.N_zu_d[idx["N2"]]
        + bed0.N_d_in[idx["N2"]]
        + bed0.N_rez_d[idx["N2"]]
        + bed0.R_gas_d[idx["N2"]]
        - bed0.N_d[idx["N2"]]
        + bed0.N_ex[idx["N2"]]
    )
    assert abs(n2_res_d) < 8.0
    assert float(bed0.N_d[idx["O2"]]) > 0.5


def test_post_fastox_h2o_passthrough_flag_off_is_noop():
    from src.solvers.vorabrechnung.vorab_cell_mapping import (
        apply_post_fastox_post_staged_h2o_passthrough,
    )

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_a_tier_post_staged_h2o_passthrough_thesis = False
    reactor = Reactor(cfg)
    diag = apply_post_fastox_post_staged_h2o_passthrough(reactor.cells, cfg)
    assert diag["enabled"] == 0.0
    assert diag["cells_updated"] == 0.0


def test_post_fastox_h2o_passthrough_closes_bed3_gap_after_full_init():
    """fastox 后 H2O 透传把 bed3 holdup 贴到 inflow（handoff §5.62）。"""
    from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = True
    cfg.vorab_bed0_r1_oxidizer_seed_mol_s_thesis = 2.0
    cfg.vorab_a_tier_post_staged_h2o_passthrough_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(reactor, init_strategy="vorabrechnung", gs_warmup_steps=None)
    idx = GAS_SPECIES_INDEX
    j = idx["H2O"]
    bed3 = reactor.cells[3]
    hold = float(bed3.N_d[j] + bed3.N_b[j])
    inf = float(bed3.N_d_in[j] + bed3.N_b_in[j])
    zu = float(bed3.N_zu_d[j] + bed3.N_zu_b[j])
    assert abs(hold - inf - zu) < 0.05
    diag = getattr(reactor, "_vorab_post_fastox_h2o_passthrough_diag", None)
    assert diag is not None
    assert diag["enabled"] == 1.0
    assert diag["cells_updated"] >= 1.0


def test_post_fastox_syngas_passthrough_flag_off_is_noop():
    from src.solvers.vorabrechnung.vorab_cell_mapping import (
        apply_post_fastox_post_staged_syngas_passthrough,
    )

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_a_tier_post_staged_syngas_passthrough_thesis = False
    reactor = Reactor(cfg)
    diag = apply_post_fastox_post_staged_syngas_passthrough(reactor.cells, cfg)
    assert diag["enabled"] == 0.0
    assert diag["cells_updated"] == 0.0


def test_post_fastox_syngas_passthrough_closes_bed3_product_gap_after_full_init():
    """fastox 后产物透传把 bed3 CO/CO2/CH4/TAR 贴到 inflow（handoff §5.63）。"""
    from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis = False
    cfg.vorab_bed0_eq24_two_phase_startwert_thesis = True
    cfg.vorab_bed0_r1_oxidizer_seed_mol_s_thesis = 2.0
    cfg.vorab_a_tier_post_staged_h2o_passthrough_thesis = False
    cfg.vorab_a_tier_post_staged_syngas_passthrough_thesis = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(reactor, init_strategy="vorabrechnung", gs_warmup_steps=None)
    idx = GAS_SPECIES_INDEX
    bed3 = reactor.cells[3]
    for sp in ("CO", "CO2", "CH4", "TAR2"):
        j = idx[sp]
        hold = float(bed3.N_d[j] + bed3.N_b[j])
        inf = float(bed3.N_d_in[j] + bed3.N_b_in[j])
        zu = float(bed3.N_zu_d[j] + bed3.N_zu_b[j])
        assert abs(hold - inf - zu) < 0.05, sp
    # 不改 H2O：CORE 过填仍在。
    j_h2o = idx["H2O"]
    h2o_hold = float(bed3.N_d[j_h2o] + bed3.N_b[j_h2o])
    h2o_in = float(bed3.N_d_in[j_h2o] + bed3.N_b_in[j_h2o])
    assert h2o_hold - h2o_in > 1.0
    # 不冲 bed0 R1 种子。
    assert float(reactor.cells[0].N_d[idx["O2"]]) > 0.5
    diag = getattr(reactor, "_vorab_post_fastox_syngas_passthrough_diag", None)
    assert diag is not None
    assert diag["enabled"] == 1.0
    assert diag["cells_updated"] >= 1.0

