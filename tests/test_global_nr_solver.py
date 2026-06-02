from __future__ import annotations

import math
import numpy as np
import pytest
import scipy.sparse as sp

from src.core.connectivity import cell_total_solid_holdup, propagate_explicit_freeboard_chain
from src.core.cell import Cell, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_kinetics import build_reaction_sources
from src.core.freeboard_bridge import (
    _spm_shrink_proxy_from_reaction_partition,
    effective_bed_top_entrained_d_p_classes,
)
from src.core.freeboard_segment import (
    _advance_particle_samples_analytical_wirsum,
    _build_size_class_bundles,
    _build_velocity_samples,
)
from src.core.reactor import Reactor, ReactorConfig, _cell_solid_outflow_component, _resolve_axial_heat_loss_distribution, _resolve_nr_init_strategy
from src.core.species import GAS_SPECIES, gas_density_ideal, gas_viscosity_power_law
from src.solvers.result_builder import build_exit_summary, finalize_global_nr_result
from src.solvers.global_nr_solver import (
    _clip_dx,
    _next_lambda_seed,
    _project_cell_residual_to_main_nr,
    _residual_gas_phase_metrics,
    _residual_group_norms,
    build_equation_scales,
    build_jacobian_fd,
    cell_offsets,
    global_residual,
    n_var,
    pack_cell,
    pack_reactor,
    unpack_reactor,
    solve_global_nr,
)
from src.solvers.nr_indexing import MAIN_NR_GAS_COUNT, gas_mode_for_nr, n_gas_var, n_solid_var
from src.solvers.structured_jacobian import build_jacobian_structure, build_solver_graph
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from src.physics.freeboard import calc_beta_a, calc_u_gb
from src.solvers.vorabrechnung import (
    estimate_axial_T_profile,
    generate_initial_x0,
    _major_gibbs_seed_is_credible,
    _major_gibbs_diag_allows_warmstart,
    _major_elements_from_feeds,
    _solve_major_gibbs_seed,
)
from src.thermal.devolatilization import devolatilization_rate_for_cell
from tests.validation_case_utils import (
    build_phase1_htw_lu_global_nr_reactor_config,
    build_phase1_htw_lu_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
    build_phase2_htw_lu_local_refined_damped_config,
    PHASE1_HTW_LU_SOLVE_KWARGS,
)


def _build_initialized_lu_reactor(n_cells: int = 3) -> Reactor:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = n_cells
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
    )
    return reactor


def test_residual_group_norms_follow_cell_by_cell_layout():
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    cell = reactor.cells[0]
    nv = n_var(cell)
    ng = 2 * MAIN_NR_GAS_COUNT
    ns = n_solid_var(cell)
    F_hat = np.zeros(nv, dtype=np.float64)
    F_hat[0] = 3.0
    F_hat[ng] = 4.0
    F_hat[nv - 1] = 5.0

    norms = _residual_group_norms(F_hat, [cell])

    assert norms["gas"] == pytest.approx(3.0 / math.sqrt(ng))
    assert norms["solid"] == pytest.approx(4.0 / math.sqrt(ns))
    assert norms["energy"] == pytest.approx(5.0)


def test_residual_gas_phase_metrics_separate_combined_and_split_modes():
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    cell = reactor.cells[0]
    nv = n_var(cell)
    F_hat = np.zeros(nv, dtype=np.float64)
    # Dense and bubble CO residuals cancel in the total gas balance but remain
    # a large phase-split residual.
    F_hat[0] = 3.0
    F_hat[MAIN_NR_GAS_COUNT] = -3.0

    metrics = _residual_gas_phase_metrics(F_hat, [cell])

    assert metrics["gas_combined_max_abs"] == pytest.approx(0.0)
    assert metrics["gas_phase_split_max_abs"] == pytest.approx(6.0)


def test_degenerate_gas_phase_metrics_ignore_bubble_anchor_as_split():
    cell = Cell()
    cell.cell_type = "freeboard"
    nv = n_var(cell)
    F_hat = np.zeros(nv, dtype=np.float64)
    F_hat[0] = 3.0
    F_hat[MAIN_NR_GAS_COUNT] = -3.0

    metrics = _residual_gas_phase_metrics(F_hat, [cell])

    assert metrics["gas_combined_max_abs"] == pytest.approx(3.0)
    assert metrics["gas_phase_split_max_abs"] == pytest.approx(0.0)


def test_freeboard_residual_projection_uses_single_total_gas_without_bubble_anchor():
    cell = Cell()
    cell.cell_type = "freeboard"
    cell.N_b[:] = 0.0
    cell.N_b[0] = 0.25
    full = np.zeros(2 * len(cell.N_d) + cell.m_solid.size + 1, dtype=np.float64)
    full[0] = 2.0
    full[len(cell.N_d)] = -0.5
    cell.residuals = lambda: full

    projected = _project_cell_residual_to_main_nr(cell)

    assert projected[0] == pytest.approx(1.5)
    assert projected.shape == (MAIN_NR_GAS_COUNT + 1,)
    assert projected[MAIN_NR_GAS_COUNT] == pytest.approx(0.0)


def test_clip_dx_allows_positive_solid_step_when_local_source_exists():
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    cell = reactor.cells[0]
    cell.solid_state_model = "holdup_transport"
    cell.m_solid.fill(0.0)
    cell.m_solid_zu.fill(0.0)
    cell.m_solid_rez.fill(0.0)
    cell.m_solid_in.fill(0.0)
    cell.m_solid_auf_in.fill(0.0)
    cell.m_solid_ab_in.fill(0.0)
    cell.R_solid.fill(0.0)
    cell.R_solid[0, S_CHAR] = 0.1
    dx = np.zeros(n_var(cell), dtype=np.float64)
    solid0 = 2 * MAIN_NR_GAS_COUNT
    dx[solid0] = 0.05

    clipped = _clip_dx(dx, [cell], ref_gas_mol_s=1.0, ref_solid_kg_s=1.0, zero_empty_solid=True)

    assert clipped[solid0] > 0.001


def test_clip_dx_limits_negative_solid_step_to_existing_inventory():
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    cell = reactor.cells[0]
    cell.solid_state_model = "holdup_transport"
    cell.m_solid.fill(0.0)
    cell.m_solid[0, S_CHAR] = 0.01
    cell.m_solid[0, S_ASH] = 0.005
    dx = np.zeros(n_var(cell), dtype=np.float64)
    solid0 = 2 * MAIN_NR_GAS_COUNT
    dx[solid0] = -10.0
    dx[solid0 + 1] = -10.0

    clipped = _clip_dx(dx, [cell], ref_gas_mol_s=1.0, ref_solid_kg_s=1.0)

    assert clipped[solid0] == pytest.approx(-0.008)
    assert clipped[solid0 + 1] == pytest.approx(-0.004)


def test_clip_dx_limits_temperature_step_when_requested():
    cell = Cell()
    dx = np.zeros(n_var(cell), dtype=np.float64)
    dx[-1] = 12000.0

    clipped = _clip_dx(dx, [cell], ref_gas_mol_s=1.0, ref_solid_kg_s=1.0, t_step_limit_K=600.0)

    assert clipped[-1] == pytest.approx(600.0)


def test_unpack_reactor_respects_vorabrechnung_temperature_fence():
    cell = Cell()
    cell._nr_temperature_min_K = 1100.0
    cell._nr_temperature_max_K = 1500.0
    x = np.zeros(n_var(cell), dtype=np.float64)
    x[-1] = 700.0

    unpack_reactor(x, [cell])

    assert cell.T == pytest.approx(1100.0)

    x[-1] = 1800.0
    unpack_reactor(x, [cell])

    assert cell.T == pytest.approx(1500.0)


def test_clip_dx_does_not_treat_side_element_last_solid_dof_as_temperature():
    cell = Cell()
    cell.cell_type = "cyclone"
    cell.solid_state_model = "holdup_transport"
    cell.m_solid.fill(10.0)
    dx = np.zeros(n_var(cell), dtype=np.float64)
    dx[-1] = 42.0

    clipped = _clip_dx(dx, [cell], ref_gas_mol_s=1.0, ref_solid_kg_s=10000.0, t_step_limit_K=1.0)

    assert clipped[-1] == pytest.approx(42.0)


def _build_initialized_thesis_freeboard_reactor(
    *,
    n_bed: int = 3,
    n_freeboard: int = 2,
) -> Reactor:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = int(n_bed)
    cfg.n_freeboard_cells = int(n_freeboard)
    cfg.major_gibbs_solver_mode = "hamel_reduced"
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    return reactor


def _phase2_first_cell_major_elements() -> tuple[dict[str, float], float, float]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    for i, cell in enumerate(reactor.cells):
        cell.T = float(t_profile[i])
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    for cell in reactor.cells:
        cell.calc_hydrodynamics()

    cell0 = reactor.cells[0]
    moisture_frac = cfg.moisture_wt / 100.0
    ash_frac = cfg.ash_dry_wt / 100.0
    vm_daf_frac = cfg.VM_daf / 100.0
    dry_feed = cfg.fuel_feed * (1.0 - moisture_frac)
    daf_feed = dry_feed * (1.0 - ash_frac)

    m_c = 12.011e-3
    m_h = 1.00794e-3
    m_o = 15.999e-3
    c_daf = (cfg.C_dry / 100.0) / max(1.0 - ash_frac, 1e-12)
    h_daf = (cfg.H_dry / 100.0) / max(1.0 - ash_frac, 1e-12)
    o_daf = (cfg.O_dry / 100.0) / max(1.0 - ash_frac, 1e-12)

    frac_height = 0.5 / max(cfg.n_cells, 1)
    o2_consumed_frac = min(1.0 - np.exp(-5.0 * frac_height), 1.0)
    o2_remaining = cfg.O2_feed * max(1.0 - o2_consumed_frac, 0.0)
    tau_est = cell0.geo.dh / max(0.05, cell0.u_mf if cell0.u_mf > 0 else 0.05)
    x_vm, _ = devolatilization_rate_for_cell(
        T_bed=float(t_profile[0]),
        tau_cell=tau_est,
        T_init=max(float(t_profile[0]) * 0.3, 400.0),
        VM_daf=vm_daf_frac,
        fuel_type="brown_coal",
    )
    x_vm_cum = min(x_vm * (frac_height + 0.1), 1.0)
    m_vm_released = daf_feed * vm_daf_frac * x_vm_cum
    n_c_vm = m_vm_released * c_daf / m_c
    n_h_vm = m_vm_released * h_daf / m_h / 2.0
    n_o_vm = m_vm_released * o_daf / m_o / 2.0
    elems = _major_elements_from_feeds(
        nC_vm=float(n_c_vm),
        nH_vm=float(n_h_vm),
        nO_vm=float(n_o_vm),
        H2O_feed=float(cfg.H2O_feed),
        N2_feed=float(cfg.N2_feed),
        o2_remaining=float(o2_remaining),
    )
    return elems, float(t_profile[0]), float(cell0.P)


def test_generate_initial_x0_major_gibbs_seed_is_finite_and_nonnegative():
    reactor = _build_initialized_lu_reactor(n_cells=2)
    for cell in reactor.cells:
        assert np.all(np.isfinite(cell.N_d))
        assert np.all(np.isfinite(cell.N_b))
        assert np.all(cell.N_d >= 0.0)
        assert np.all(cell.N_b >= 0.0)
        assert float(np.sum(cell.N_d + cell.N_b)) > 0.0


def test_generate_initial_x0_falls_back_when_major_gibbs_seed_raises(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    monkeypatch.setattr(
        "src.solvers.vorabrechnung._solve_major_gibbs_seed",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=True,
    )

    for cell in reactor.cells:
        assert np.all(np.isfinite(cell.N_d))
        assert np.all(np.isfinite(cell.N_b))
        assert float(np.sum(cell.N_d + cell.N_b)) > 0.0


def test_generate_initial_x0_falls_back_when_major_gibbs_seed_is_not_credible(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )

    monkeypatch.setattr(
        "src.solvers.vorabrechnung._solve_major_gibbs_seed",
        lambda **kwargs: (
            {
                "CO2": 1e24,
                "CO": 0.0,
                "CH4": 0.0,
                "H2": 0.0,
                "H2O": 1e24,
                "O2": 1e24,
                "N2": 1e24,
            },
            {"converged": False, "lambda": np.zeros(4), "ln_N": 10.0},
        ),
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=True,
    )

    for cell in reactor.cells:
        assert np.all(np.isfinite(cell.N_d))
        assert np.all(np.isfinite(cell.N_b))
        assert float(np.sum(cell.N_d + cell.N_b)) < 1e6


def test_generate_initial_x0_strict_mode_raises_when_major_gibbs_seed_fails(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    monkeypatch.setattr(
        "src.solvers.vorabrechnung._solve_major_gibbs_seed",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )

    with pytest.raises(RuntimeError, match="boom"):
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
            T_profile=t_profile,
            fuel_type=cfg.fuel_type,
            use_hamel_major_gibbs_x0=True,
            strict_hamel_major_gibbs_x0=True,
        )


def test_generate_initial_x0_strict_mode_rejects_noncredible_finite_seed(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    monkeypatch.setattr(
        "src.solvers.vorabrechnung._solve_major_gibbs_seed",
        lambda **kwargs: (
            {
                "CO2": 1e20,
                "CO": 0.0,
                "CH4": 0.0,
                "H2": 0.0,
                "H2O": 1e20,
                "O2": 1e20,
                "N2": 1e20,
            },
            {"converged": False, "final_residual": 1e20, "lambda": np.ones(4), "ln_N": 5.0},
        ),
    )

    with pytest.raises(RuntimeError, match="non-credible seed in strict mode"):
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
            T_profile=t_profile,
            fuel_type=cfg.fuel_type,
            use_hamel_major_gibbs_x0=True,
            strict_hamel_major_gibbs_x0=True,
        )


def test_major_gibbs_hamel_reduced_seed_is_credible_on_phase2_first_cell():
    elems, t0, p0 = _phase2_first_cell_major_elements()
    result, diag = _solve_major_gibbs_seed(
        T=t0,
        P=p0,
        elements=elems,
        lambda0=None,
        ln_N0=None,
        solver_mode="hamel_reduced",
    )
    assert diag.get("solver_mode") == "hamel_reduced"
    assert _major_gibbs_seed_is_credible(guess=result, diag=diag, elements=elems) is True
    assert float(sum(result.values())) > 0.0


def test_major_gibbs_shadow_compare_reports_both_solver_diagnostics():
    elems, t0, p0 = _phase2_first_cell_major_elements()
    result, diag = _solve_major_gibbs_seed(
        T=t0,
        P=p0,
        elements=elems,
        lambda0=None,
        ln_N0=None,
        solver_mode="shadow_compare",
    )
    assert diag.get("solver_mode") == "shadow_compare"
    assert "shadow_compare" in diag
    assert "augmented" in diag["shadow_compare"]
    assert "hamel_reduced" in diag["shadow_compare"]
    assert diag.get("selected_solver") == "hamel_reduced"
    assert float(sum(result.values())) > 0.0


def test_thesis_major_gibbs_mode_defaults_to_hamel_reduced():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.major_gibbs_solver_mode = "augmented"
    reactor = Reactor(cfg)
    assert reactor.config.major_gibbs_solver_mode == "hamel_reduced"


def test_init_precalc_uses_hamel_reduced_major_gibbs_by_default_in_thesis(monkeypatch):
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 2
    cfg.n_freeboard_cells = 1
    cfg.major_gibbs_solver_mode = "augmented"
    reactor = Reactor(cfg)
    captured_modes: list[str] = []

    def _fake_solver(**kwargs):
        captured_modes.append(str(kwargs.get("solver_mode")))
        return (
            {"CO2": 0.8, "CO": 0.2, "CH4": 0.0, "H2": 0.3, "H2O": 0.5, "O2": 0.1, "N2": 1.0},
            {"converged": True, "final_residual": 1e-12, "lambda": np.ones(4), "ln_N": 0.0},
        )

    monkeypatch.setattr("src.solvers.vorabrechnung._solve_major_gibbs_seed", _fake_solver)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    assert captured_modes
    assert set(captured_modes) == {"hamel_reduced"}


def test_generate_initial_x0_warmstart_not_propagated_from_high_residual_cell(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    captured: list[tuple[np.ndarray | None, float | None]] = []

    def _fake_solver(**kwargs):
        captured.append((kwargs.get("lambda0"), kwargs.get("ln_N0")))
        return (
            {"CO2": 0.5, "CO": 0.5, "CH4": 0.2, "H2": 0.5, "H2O": 0.5, "O2": 0.2, "N2": 0.8},
            {
                "converged": True,
                "final_residual": 1e-3,  # above warm-start threshold
                "lambda": np.ones(4),
                "ln_N": 0.5,
            },
        )

    monkeypatch.setattr("src.solvers.vorabrechnung._solve_major_gibbs_seed", _fake_solver)
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=True,
    )
    assert len(captured) == 2
    assert captured[0][0] is None and captured[0][1] is None
    assert captured[1][0] is None and captured[1][1] is None


def test_generate_initial_x0_warmstart_propagated_only_from_low_residual_cell(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    captured: list[tuple[np.ndarray | None, float | None]] = []

    def _fake_solver(**kwargs):
        captured.append((kwargs.get("lambda0"), kwargs.get("ln_N0")))
        return (
            {"CO2": 0.5, "CO": 0.5, "CH4": 0.2, "H2": 0.5, "H2O": 0.5, "O2": 0.2, "N2": 0.8},
            {
                "converged": True,
                "final_residual": 1e-12,  # below warm-start threshold
                "lambda": np.ones(4),
                "ln_N": 0.5,
            },
        )

    monkeypatch.setattr("src.solvers.vorabrechnung._solve_major_gibbs_seed", _fake_solver)
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=True,
    )
    assert len(captured) == 2
    assert captured[0][0] is None and captured[0][1] is None
    assert isinstance(captured[1][0], np.ndarray)
    assert captured[1][0] is not None and captured[1][0].shape == (4,)
    assert captured[1][1] == pytest.approx(0.5)


def test_major_gibbs_diag_allows_warmstart_requires_converged_and_low_residual():
    assert _major_gibbs_diag_allows_warmstart({"converged": True, "final_residual": 1e-12}) is True
    assert _major_gibbs_diag_allows_warmstart({"converged": True, "final_residual": 1e-3}) is False
    assert _major_gibbs_diag_allows_warmstart({"converged": False, "final_residual": 1e-12}) is False
    assert _major_gibbs_diag_allows_warmstart({"converged": True}) is False


def test_major_gibbs_diag_warmstart_rejects_shadow_compare_when_reduced_not_selected():
    assert (
        _major_gibbs_diag_allows_warmstart(
            {
                "converged": True,
                "final_residual": 1e-12,
                "solver_mode": "shadow_compare",
                "selected_solver": "augmented",
            }
        )
        is False
    )
    assert (
        _major_gibbs_diag_allows_warmstart(
            {
                "converged": True,
                "final_residual": 1e-12,
                "solver_mode": "shadow_compare",
                "selected_solver": "hamel_reduced",
            }
        )
        is True
    )


def test_major_gibbs_diag_warmstart_rejects_relaxed_or_large_closure():
    assert (
        _major_gibbs_diag_allows_warmstart(
            {
                "converged": True,
                "final_residual": 1e-12,
                "converged_by": "closure_relaxed",
                "closure_rel": 1e-8,
            }
        )
        is False
    )
    assert (
        _major_gibbs_diag_allows_warmstart(
            {
                "converged": True,
                "final_residual": 1e-12,
                "closure_rel": 1e-3,
            }
        )
        is False
    )
    assert (
        _major_gibbs_diag_allows_warmstart(
            {
                "converged": True,
                "final_residual": 1e-12,
                "closure_rel": 1e-9,
            }
        )
        is True
    )


def test_generate_initial_x0_uses_holdup_seed_for_holdup_transport_cells():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 1
    reactor = Reactor(cfg)
    reactor.cells[0].solid_state_model = "holdup_transport"
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)

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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
    )

    expected_total = cell_total_solid_holdup(reactor.cells[0])
    active_seed = float(np.sum(np.maximum(reactor.cells[0].m_solid_zu + reactor.cells[0].m_solid_in + reactor.cells[0].m_solid_rez, 0.0)))
    assert float(np.sum(reactor.cells[0].m_solid)) > 0.0
    assert float(np.sum(reactor.cells[0].m_solid)) < expected_total
    assert float(np.sum(reactor.cells[0].m_solid)) <= active_seed + 1e-9


def test_generate_initial_x0_holdup_seed_depletes_reactive_components_axially():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 3
    reactor = Reactor(cfg)
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)

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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
    )

    vm0 = float(np.sum(reactor.cells[0].m_solid[:, S_VM]))
    vm2 = float(np.sum(reactor.cells[-1].m_solid[:, S_VM]))
    moist0 = float(np.sum(reactor.cells[0].m_solid[:, S_MOISTURE]))
    moist2 = float(np.sum(reactor.cells[-1].m_solid[:, S_MOISTURE]))
    char0 = float(np.sum(reactor.cells[0].m_solid[:, S_CHAR]))
    char2 = float(np.sum(reactor.cells[-1].m_solid[:, S_CHAR]))

    assert vm2 < vm0
    assert moist2 <= moist0
    assert char2 <= char0
    for cell in reactor.cells:
        assert float(np.sum(cell.m_solid)) <= cell_total_solid_holdup(cell) + 1e-9


def test_generate_initial_x0_holdup_seed_uses_active_solid_inflows_not_total_bed_inventory():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)

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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
    )

    top = reactor.cells[-1]
    top_total = float(np.sum(top.m_solid))
    top_active_in = float(np.sum(np.maximum(top.m_solid_zu + top.m_solid_in + top.m_solid_rez, 0.0)))
    assert top_total <= top_active_in + 1e-9
    assert top_total < cell_total_solid_holdup(top)


def test_generate_initial_x0_holdup_seed_uses_upstream_solid_temperature(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
    reactor.cells[0].T_in_solid = 293.15
    reactor.cells[1].T_in_solid = 1180.0
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )

    seen_t_init: list[float] = []
    orig_dev = devolatilization_rate_for_cell

    def _recording_devolatilization_rate_for_cell(*, T_bed, tau_cell, T_init, VM_daf, fuel_type):
        seen_t_init.append(float(T_init))
        return orig_dev(T_bed=T_bed, tau_cell=tau_cell, T_init=T_init, VM_daf=VM_daf, fuel_type=fuel_type)

    monkeypatch.setattr(
        "src.solvers.vorabrechnung.devolatilization_rate_for_cell",
        _recording_devolatilization_rate_for_cell,
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
    )

    assert seen_t_init[:2] == pytest.approx([293.15, 1180.0], rel=1e-9)


def test_generate_initial_x0_uses_hydrodynamics_based_phase_split():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 1
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    cell.eps_b = 0.12
    cell.eps_d_voidage = 0.48
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
    )

    total_b = float(np.sum(cell.N_b))
    total_d = float(np.sum(cell.N_d))
    bubble_share = total_b / max(total_b + total_d, 1e-12)
    assert bubble_share == pytest.approx(0.2, rel=1e-3)


def test_generate_initial_x0_aligns_bottom_feed_species_with_frozen_inlet_split():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 1
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    cell.eps_b = 0.12
    cell.eps_d_voidage = 0.48
    cell._vorab_bottom_gas_inlet_dense_frac = 0.25
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
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
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=False,
    )

    for sp in ("O2", "H2O", "N2"):
        j = GAS_SPECIES.index(sp)
        total = float(cell.N_d[j] + cell.N_b[j])
        assert cell.N_d[j] / total == pytest.approx(0.25)


def test_run_init_and_precalc_ignores_stale_gas_inventory_in_macro_hydrodynamics_seed():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 3

    reactor_fresh = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor_fresh,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )

    reactor_stale = Reactor(build_phase1_htw_lu_reactor_config())
    reactor_stale.config.n_cells = 3
    reactor_stale = Reactor(reactor_stale.config)
    for i, cell in enumerate(reactor_stale.cells):
        cell.N_d[:] = 25.0 + float(i)
        cell.N_b[:] = 5.0 + 0.5 * float(i)
        cell.T = 900.0 + 10.0 * float(i)

    run_init_and_precalc_for_global_nr(
        reactor_stale,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )

    for fresh, stale in zip(reactor_fresh.cells, reactor_stale.cells):
        np.testing.assert_allclose(stale.N_d, fresh.N_d, rtol=0.0, atol=1e-10)
        np.testing.assert_allclose(stale.N_b, fresh.N_b, rtol=0.0, atol=1e-10)
        assert stale.eps_b == pytest.approx(fresh.eps_b, abs=1e-12)
        assert stale.u0 == pytest.approx(fresh.u0, abs=1e-12)

    total_b = float(np.sum(reactor_fresh.cells[0].N_b))
    total_d = float(np.sum(reactor_fresh.cells[0].N_d))
    bubble_share = total_b / max(total_b + total_d, 1e-12)
    assert bubble_share > 0.30


def test_init_precalc_seeds_freeboard_holdup_from_bed_top_entrainment():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    bed_top = reactor.cells[-1]

    assert float(np.sum(np.maximum(bed_top._solid_upflow_rates()[:, S_CHAR], 0.0))) > 0.0
    assert reactor._last_explicit_freeboard_closure is not None
    assert float(np.sum(reactor._last_explicit_freeboard_closure["profiles"]["solid_holdup_char_classes_kg"][0])) > 0.0

    fb0 = reactor.freeboard_cells[0]
    assert float(np.sum(np.maximum(fb0.m_solid[:, S_CHAR], 0.0))) > 0.0
    assert float(np.sum(np.maximum(fb0.m_solid_in[:, S_CHAR], 0.0))) > 0.0
    assert fb0.T_in_solid == pytest.approx(bed_top.T)
    assert float(np.sum(np.maximum(fb0.m_solid_auf_in[:, S_CHAR], 0.0))) == 0.0
    assert float(np.sum(np.maximum(fb0.m_solid_ab_in[:, S_CHAR], 0.0))) == 0.0
    assert float(np.sum(np.maximum(reactor.freeboard_cells[1].m_solid[:, S_CHAR], 0.0))) >= 0.0
    assert float(max(reactor._last_explicit_freeboard_closure.get("entrained_return_char_kg_s", 0.0), 0.0)) > 0.0


def test_build_exit_summary_exposes_freeboard_solid_holdup_profiles_for_exact_hamel():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor.config.freeboard_trajectory_coeff_model = "exact_hamel"
    fb = reactor._refresh_explicit_freeboard_transport_from_closure()

    assert fb is not None

    result = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )

    assert result["freeboard_active"] is True
    assert result["freeboard_trajectory_coeff_model"] == "exact_hamel"
    assert len(result["freeboard_solid_holdup_char_profile_kg"]) == 2
    assert len(result["freeboard_solid_holdup_ash_profile_kg"]) == 2
    assert len(result["freeboard_char_reaction_source_profile_mol_s"]) == 2
    assert result["freeboard_solid_holdup_char_profile_kg"][0] > 0.0
    assert result["freeboard_solid_holdup_char_profile_kg"][1] >= 0.0
    assert result["freeboard_solid_holdup_ash_profile_kg"][0] > 0.0
    assert result["freeboard_solid_holdup_ash_profile_kg"][1] >= 0.0
    assert abs(result["freeboard_char_reaction_source_profile_mol_s"][0]) > 0.0
    assert abs(result["freeboard_char_reaction_source_profile_mol_s"][1]) >= 0.0
    for name in ("R1", "R2", "R3", "R4"):
        assert name in result["freeboard_reaction_diag_impl"]
        assert name in result["freeboard_explicit_reaction_diag_impl"]
        assert len(result["freeboard_explicit_reaction_diag_impl"][name]) == 2
    np.testing.assert_allclose(
        result["freeboard_solid_holdup_char_profile_kg"],
        [float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) for cell in reactor.freeboard_cells],
    )
    np.testing.assert_allclose(
        result["freeboard_solid_holdup_ash_profile_kg"],
        [float(np.sum(np.maximum(cell.m_solid[:, S_ASH], 0.0))) for cell in reactor.freeboard_cells],
    )


def test_build_exit_summary_exposes_bottom_recycle_and_bed_top_solid_enthalpy_diagnostics():
    reactor = Reactor(ReactorConfig(n_cells=2, thesis_mode=True, n_freeboard_cells=0, n_age_classes=1))
    bot, top = reactor.cells
    bot.m_solid_rez[0, S_CHAR] = 0.09
    bot.m_solid_rez[0, S_ASH] = 0.02
    bot.T_rez_solid = 1125.0
    bot.m_solid_zu[0, S_CHAR] = 0.3
    bot.T_zu_solid = 293.15
    top.m_solid_ab_in[0, S_CHAR] = 0.12
    top.T_solid_ab_in = 1180.0

    result = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )

    assert result["bed0_recycle_solid_char_kg_s"] == pytest.approx(0.09)
    assert result["bed0_recycle_solid_ash_kg_s"] == pytest.approx(0.02)
    assert result["bed0_recycle_solid_T_K"] == pytest.approx(1125.0)
    assert result["bed0_recycle_solid_enthalpy_W"] == pytest.approx(
        bot._calc_solid_enthalpy_flow(bot.m_solid_rez, bot.T_rez_solid)
    )
    assert result["bed0_fresh_solid_enthalpy_W"] == pytest.approx(
        bot._calc_solid_enthalpy_flow(bot.m_solid_zu, bot.T_zu_solid)
    )
    assert result["bed_top_downflow_solid_enthalpy_W"] == pytest.approx(
        top._calc_solid_enthalpy_flow(top.m_solid_ab_in, top.T_solid_ab_in)
    )


def test_finalize_result_reports_freeboard_closure_gap_without_mutating_nr_state():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    bed_top = reactor.cells[-1]
    old_top_up = float(reactor._last_explicit_freeboard_closure["bed_top_up_char_kg_s"])
    bed_top.K_solid_auf[:, S_CHAR] *= 2.0
    current_top_up = float(np.sum(np.maximum(bed_top._solid_upflow_rates()[:, S_CHAR], 0.0)))
    assert current_top_up > old_top_up

    fb0 = reactor.freeboard_cells[0]
    fb0.T = 987.0
    fb0.N_d[:] = np.linspace(0.1, 0.2, len(GAS_SPECIES))
    t_before = float(fb0.T)
    n_d_before = np.array(fb0.N_d, copy=True)
    bed_top_ab_in_before = np.array(bed_top.m_solid_ab_in, copy=True)

    result = finalize_global_nr_result(
        reactor,
        {
            "converged": False,
            "converged_outer": False,
            "converged_inner_nr": False,
            "rms_scaled_final": 1.0,
            "n_iter": 1,
        },
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )

    assert result["freeboard_result_closure_refreshed"] is False
    assert result["freeboard_result_boundary_refreshed"] is False
    assert result["freeboard_pre_result_refresh_bed_top_up_char_kg_s"] == pytest.approx(old_top_up)
    assert result["freeboard_pre_result_refresh_current_bed_top_up_char_kg_s"] == pytest.approx(current_top_up)
    assert result["freeboard_pre_result_refresh_gap_current_bed_top_kg_s"] == pytest.approx(current_top_up - old_top_up)
    assert result["freeboard_current_bed_top_up_char_kg_s"] == pytest.approx(current_top_up)
    assert result["freeboard_bed_top_up_char_kg_s"] == pytest.approx(old_top_up)
    assert result["freeboard_post_result_refresh_gap_current_bed_top_kg_s"] == pytest.approx(current_top_up - old_top_up)
    assert result["freeboard_result_closure_refresh_iters"] == 0
    np.testing.assert_allclose(bed_top.m_solid_ab_in, bed_top_ab_in_before)
    assert reactor.freeboard_cells[0].T == pytest.approx(t_before)
    np.testing.assert_allclose(reactor.freeboard_cells[0].N_d, n_d_before)


def test_exact_hamel_changes_freeboard_holdup_distribution_vs_stable_initialized_case():
    stable = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    stable.config.freeboard_trajectory_coeff_model = "stable_mixed_drag_split"
    stable._refresh_explicit_freeboard_transport_from_closure()

    exact = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    exact.config.freeboard_trajectory_coeff_model = "exact_hamel"
    exact._refresh_explicit_freeboard_transport_from_closure()

    stable_hold = [float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) for cell in stable.freeboard_cells]
    exact_hold = [float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) for cell in exact.freeboard_cells]

    assert stable_hold[0] > 0.0 and stable_hold[1] >= 0.0
    assert exact_hold[0] > 0.0 and exact_hold[1] >= 0.0
    assert abs(exact_hold[0] - stable_hold[0]) > 1e-6
    assert abs(exact_hold[1] - stable_hold[1]) >= 0.0


def test_exact_hamel_freeboard_top_char_holdup_activates_char_reaction_source():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor._refresh_explicit_freeboard_transport_from_closure()
    fb_top = next(
        cell
        for cell in reactor.freeboard_cells
        if float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) > 0.0
    )

    assert float(np.sum(np.maximum(fb_top.m_solid[:, S_CHAR], 0.0))) > 0.0

    fb_top.calc_exchange()
    thermo = fb_top._get_local_thermo_bundle()
    bundle = build_reaction_sources(
        T=fb_top.T,
        P=fb_top.P,
        fuel_type=fb_top.fuel_type,
        V_b=fb_top.V_b,
        V_d=fb_top.V_d,
        C_b=thermo["C_b"],
        C_d=thermo["C_d"],
        y_b=thermo["y_b"],
        y_d=thermo["y_d"],
        gas_src_vm=fb_top._vm_gas_source_cache,
        solid_sink_vm=fb_top._vm_solid_sink_cache,
        areas=fb_top._calc_char_surface_area_per_class(),
        solid_d_p=fb_top.solid.d_p,
        D_g=float(thermo["D_g"]),
        char_conversion=fb_top._compute_char_conversion(),
        rho_cat=fb_top._catalyst_bulk_density(),
        enable_r12=fb_top.enable_r12,
        use_gibbs_minor=fb_top.use_gibbs_minor,
        gibbs_minor_sources=None,
        r4_scale=fb_top.r4_scale,
        r5_scale=fb_top.r5_scale,
        r6_scale=fb_top.r6_scale,
        r7_scale=fb_top.r7_scale,
        rate_multiplier=1.0,
        N_zu_d=fb_top.N_zu_d,
        N_d_in=fb_top.N_d_in,
        N_zu_b=fb_top.N_zu_b,
        N_b_in=fb_top.N_b_in,
        N_rez_d=fb_top.N_rez_d,
        N_rez_b=fb_top.N_rez_b,
        N_ex=fb_top.N_ex,
        solid_shape=fb_top.R_solid.shape,
        char_index=S_CHAR,
    )

    assert abs(float(bundle.net_molar_gas_source_char)) > 0.0


def test_build_exit_summary_can_separate_freeboard_closure_and_explicit_char_hetero_layers():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor.config.freeboard_closure_enable_char_hetero = True
    reactor.config.freeboard_explicit_enable_char_hetero = False
    reactor._refresh_explicit_freeboard_transport_from_closure()

    result = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )
    assert any(abs(v) > 0.0 for v in result["freeboard_reaction_diag_impl"]["R2"])
    assert all(abs(v) == 0.0 for v in result["freeboard_explicit_reaction_diag_impl"]["R1"])
    assert all(abs(v) == 0.0 for v in result["freeboard_explicit_reaction_diag_impl"]["R2"])

    reactor2 = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor2.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor2.config.freeboard_closure_enable_char_hetero = False
    reactor2.config.freeboard_explicit_enable_char_hetero = True
    reactor2._refresh_explicit_freeboard_transport_from_closure()

    result2 = build_exit_summary(
        reactor2,
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )
    assert all(abs(v) == 0.0 for v in result2["freeboard_reaction_diag_impl"]["R1"])
    assert all(abs(v) == 0.0 for v in result2["freeboard_reaction_diag_impl"]["R2"])
    assert any(abs(v) > 0.0 for v in result2["freeboard_explicit_reaction_diag_impl"]["R2"])


def test_explicit_freeboard_char_hetero_depletes_projected_char_holdup() -> None:
    reactor_off = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor_off.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor_off.config.freeboard_closure_enable_char_hetero = False
    reactor_off.config.freeboard_explicit_enable_char_hetero = False
    reactor_off._refresh_explicit_freeboard_transport_from_closure()

    reactor_on = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor_on.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor_on.config.freeboard_closure_enable_char_hetero = False
    reactor_on.config.freeboard_explicit_enable_char_hetero = True
    reactor_on._refresh_explicit_freeboard_transport_from_closure()

    hold_off = [float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0))) for c in reactor_off.freeboard_cells]
    hold_on = [float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0))) for c in reactor_on.freeboard_cells]
    sink_on = [
        float(np.sum(np.maximum(np.asarray(c.freeboard_explicit_char_sink_applied_kg, dtype=np.float64), 0.0)))
        for c in reactor_on.freeboard_cells
    ]

    assert any(v > 0.0 for v in sink_on)
    assert any(h_on < h_off for h_on, h_off in zip(hold_on, hold_off))


def test_build_exit_summary_exposes_freeboard_explicit_char_sink_profile() -> None:
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor.config.freeboard_closure_enable_char_hetero = False
    reactor.config.freeboard_explicit_enable_char_hetero = True
    reactor._refresh_explicit_freeboard_transport_from_closure()

    result = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )

    assert "freeboard_explicit_char_sink_applied_profile_kg" in result
    assert len(result["freeboard_explicit_char_sink_applied_profile_kg"]) == 2
    assert any(v > 0.0 for v in result["freeboard_explicit_char_sink_applied_profile_kg"])


def test_build_exit_summary_exposes_freeboard_closure_hold_up_before_after_audit() -> None:
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    reactor.config.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor.config.freeboard_closure_enable_char_hetero = False
    reactor.config.freeboard_explicit_enable_char_hetero = True
    reactor._refresh_explicit_freeboard_transport_from_closure()

    result = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
        cell_solid_outflow_component_fn=_cell_solid_outflow_component,
    )

    before = result["freeboard_solid_holdup_char_before_profile_kg"]
    after = result["freeboard_solid_holdup_char_profile_kg"]
    sink = result["freeboard_explicit_char_sink_applied_profile_kg"]
    delta = result["freeboard_transport_char_delta_profile_kg"]
    tau_eq_before = result["freeboard_equiv_char_residence_time_before_profile_s"]
    tau_eq_after = result["freeboard_equiv_char_residence_time_profile_s"]
    hold_time_mean = result["freeboard_hold_time_mean_profile_s"]
    hold_time_max = result["freeboard_hold_time_max_profile_s"]
    hold_tau_mean = result["freeboard_hold_time_to_tau_mean_ratio_profile"]
    hold_tau_max = result["freeboard_hold_time_to_tau_max_ratio_profile"]
    bed_char = result["reactor_bed_char_inventory_kg"]
    visible_char = result["reactor_visible_char_inventory_kg"]
    hold_before_bed_ratio = result["freeboard_char_holdup_before_to_bed_inventory_ratio_profile"]
    hold_after_bed_ratio = result["freeboard_char_holdup_to_bed_inventory_ratio_profile"]
    hold_before_ratio = result["freeboard_char_holdup_before_to_visible_inventory_ratio_profile"]
    hold_after_ratio = result["freeboard_char_holdup_to_visible_inventory_ratio_profile"]
    d_p_input = result["freeboard_bed_top_d_p_input_m"]
    d_p_eff = result["freeboard_bed_top_d_p_eff_m"]
    coeff_diag = result["freeboard_trajectory_coeff_diag"]
    x_top = result["freeboard_bed_top_char_conversion"]
    x_top_used = result["freeboard_bed_top_char_conversion_used"]
    x_top_local = result["freeboard_bed_top_char_conversion_local"]
    x_top_total = result["freeboard_bed_top_char_conversion_proxy_total"]
    x_top_r1_share = result["freeboard_bed_top_combustion_share_proxy"]
    x_top_comb = result["freeboard_bed_top_combustion_conversion_proxy"]
    x_top_gas = result["freeboard_bed_top_gasification_age_proxy"]
    age_quad = result["freeboard_bed_top_age_quadrature_diag"]
    top_r1_sink = result["freeboard_bed_top_r1_char_consumption_kg_s"]
    top_hetero_sink = result["freeboard_bed_top_hetero_char_consumption_kg_s"]
    top_up_char = result["freeboard_bed_top_up_char_kg_s"]

    assert len(before) == len(after) == len(sink) == len(delta) == len(tau_eq_before) == len(tau_eq_after) == len(hold_time_mean) == len(hold_time_max) == len(hold_tau_mean) == len(hold_tau_max) == len(hold_before_bed_ratio) == len(hold_after_bed_ratio) == len(hold_before_ratio) == len(hold_after_ratio) == 2
    assert any(v > 0.0 for v in sink)
    assert any(a < b for a, b in zip(after, before))
    for d, a, b in zip(delta, after, before):
        assert abs(d - (a - b)) < 1e-9
    assert all(tb >= ta for tb, ta in zip(tau_eq_before, tau_eq_after))
    assert all(hmax >= hmean >= 0.0 for hmean, hmax in zip(hold_time_mean, hold_time_max))
    assert bed_char > 0.0
    assert visible_char > 0.0
    assert all(rmax >= rmean >= 0.0 for rmean, rmax in zip(hold_tau_mean, hold_tau_max))
    assert all(r >= 0.0 for r in hold_before_bed_ratio)
    assert all(r >= 0.0 for r in hold_after_bed_ratio)
    assert all(r >= 0.0 for r in hold_before_ratio)
    assert all(r >= 0.0 for r in hold_after_ratio)
    assert len(d_p_input) == len(d_p_eff) == 10
    assert isinstance(coeff_diag, dict)
    assert result["freeboard_bed_top_d_p_eff_min_m"] == pytest.approx(min(d_p_eff))
    assert result["freeboard_bed_top_d_p_eff_max_m"] == pytest.approx(max(d_p_eff))
    assert result["freeboard_required_dp_scale_for_smallest_eff_to_grav"] >= 0.0
    assert result["freeboard_required_spm_conversion_for_smallest_input_to_grav"] >= 0.0
    assert x_top >= 0.0
    assert x_top_used >= 0.0
    assert x_top_local >= 0.0
    assert x_top_total >= 0.0
    assert 0.0 <= x_top_r1_share <= 1.0
    assert x_top_comb >= 0.0
    assert x_top_gas >= 0.0
    assert x_top >= x_top_comb
    assert isinstance(age_quad, dict)
    assert int(age_quad.get("expanded_class_count", 0)) >= len(d_p_input)
    assert top_r1_sink >= 0.0
    assert top_hetero_sink >= 0.0
    assert top_up_char >= 0.0
    assert all(d_eff <= d_in for d_eff, d_in in zip(d_p_eff, d_p_input))


def test_freeboard_cells_exclude_solid_holdup_from_global_nr_unknown_vector():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    bed = reactor.cells[0]
    fb0 = reactor.freeboard_cells[0]

    assert n_solid_var(bed) > 0
    assert n_solid_var(fb0) == 0
    assert gas_mode_for_nr(fb0) == "single"
    assert n_gas_var(fb0) == MAIN_NR_GAS_COUNT
    assert fb0.solid_state_model == "freeboard_closure"
    assert n_var(fb0) == MAIN_NR_GAS_COUNT + 1

    x_fb = pack_cell(fb0)
    assert x_fb.shape == (MAIN_NR_GAS_COUNT + 1,)
    assert float(np.sum(np.maximum(fb0.m_solid[:, S_CHAR], 0.0))) > 0.0
    np.testing.assert_allclose(fb0.calc_solid_balance(), 0.0)


def test_phase2_solver_vector_excludes_freeboard_solid_dofs():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    solver_cells = reactor._solver_cells_for_nr()
    offsets = cell_offsets(solver_cells)

    fb0_idx = len(reactor.cells)
    fb1_idx = fb0_idx + 1
    x = pack_reactor(solver_cells)

    for idx in (fb0_idx, fb1_idx):
        start = offsets[idx]
        stop = offsets[idx + 1]
        assert stop - start == MAIN_NR_GAS_COUNT + 1
    assert x.shape[0] == offsets[-1]


def test_phase2_solver_vector_uses_cell_type_specific_layouts():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    solver_cells = reactor._solver_cells_for_nr()

    expected = {
        "bed": ("two_phase", 2 * MAIN_NR_GAS_COUNT, solver_cells[0].solid.n_size_classes * 2),
        "freeboard": ("single", MAIN_NR_GAS_COUNT, 0),
        "cyclone": ("none", 0, solver_cells[0].solid.n_size_classes * 2),
        "return_leg": ("none", 0, solver_cells[0].solid.n_size_classes * 2),
    }

    for cell in solver_cells:
        kind = getattr(cell, "cell_type", "bed")
        gas_mode, ng, ns = expected[kind]
        assert gas_mode_for_nr(cell) == gas_mode
        assert n_gas_var(cell) == ng
        assert n_solid_var(cell) == ns
        assert len(_project_cell_residual_to_main_nr(cell)) == n_var(cell)


def test_bed_to_freeboard_interface_routes_total_gas_only():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    top = reactor.cells[-1]
    fb0 = reactor.freeboard_cells[0]

    total = np.linspace(1.0, 2.0, top.N_d.size)
    top.N_d[:] = 0.25 * total
    top.N_b[:] = 0.75 * total
    propagate_explicit_freeboard_chain(reactor.cells, reactor.freeboard_cells)
    inlet_a = fb0.N_d_in.copy()

    top.N_d[:] = 0.90 * total
    top.N_b[:] = 0.10 * total
    propagate_explicit_freeboard_chain(reactor.cells, reactor.freeboard_cells)

    np.testing.assert_allclose(fb0.N_d_in, inlet_a)
    np.testing.assert_allclose(fb0.N_d_in, total)
    np.testing.assert_allclose(fb0.N_b_in, 0.0)


def test_vorabrechnung_refresh_updates_freeboard_ghost_bubble_bridge():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    top = reactor.cells[-1]
    fb0 = reactor.freeboard_cells[0]

    old_bridge_u = float(fb0.freeboard_u_bed_top)
    old_bridge_d = float(fb0.freeboard_d_b_bed_top)
    top.N_b *= 2.0
    top.N_d *= 0.7

    reactor._refresh_vorabrechnung_sources_for_nr(force=False)

    assert float(top.u_b) != pytest.approx(old_bridge_u)
    assert float(top.d_b) != pytest.approx(old_bridge_d)
    assert float(fb0.freeboard_u_bed_top) == pytest.approx(float(top.u_b))
    assert float(fb0.freeboard_d_b_bed_top) == pytest.approx(float(top.d_b))
    assert float(fb0._vorab_hydro_cache["u_d"]) == pytest.approx(float(fb0.u_d))


def test_single_gas_freeboard_residual_depends_on_total_gas_not_phase_split():
    cell = Cell()
    cell.cell_type = "freeboard"
    full_a = np.zeros(2 * len(cell.N_d) + cell.m_solid.size + 1, dtype=np.float64)
    full_b = np.zeros_like(full_a)
    total = np.linspace(0.1, 0.7, len(cell.N_d))
    full_a[: len(cell.N_d)] = 0.2 * total
    full_a[len(cell.N_d) : 2 * len(cell.N_d)] = 0.8 * total
    full_b[: len(cell.N_d)] = 0.9 * total
    full_b[len(cell.N_d) : 2 * len(cell.N_d)] = 0.1 * total

    cell.residuals = lambda: full_a
    projected_a = _project_cell_residual_to_main_nr(cell)
    cell.residuals = lambda: full_b
    projected_b = _project_cell_residual_to_main_nr(cell)

    np.testing.assert_allclose(projected_a, projected_b)


def test_freeboard_reaction_and_energy_bridge_ignore_upstream_phase_split():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    top = reactor.cells[-1]
    fb0 = reactor.freeboard_cells[0]
    total = np.linspace(1.0, 2.0, top.N_d.size)

    def evaluate_with_split(dense_fraction: float) -> dict[str, np.ndarray | float]:
        top.N_d[:] = float(dense_fraction) * total
        top.N_b[:] = (1.0 - float(dense_fraction)) * total
        propagate_explicit_freeboard_chain(reactor.cells, reactor.freeboard_cells)
        full = fb0.residuals()
        return {
            "N_d_in": fb0.N_d_in.copy(),
            "N_b_in": fb0.N_b_in.copy(),
            "R_gas_d": fb0.R_gas_d.copy(),
            "R_gas_b": fb0.R_gas_b.copy(),
            "full": full.copy(),
            "projected": _project_cell_residual_to_main_nr(fb0),
            "energy": float(full[-1]),
        }

    low_dense = evaluate_with_split(0.2)
    high_dense = evaluate_with_split(0.9)

    np.testing.assert_allclose(low_dense["N_d_in"], high_dense["N_d_in"], atol=1e-14)
    np.testing.assert_allclose(low_dense["N_b_in"], 0.0, atol=1e-14)
    np.testing.assert_allclose(high_dense["N_b_in"], 0.0, atol=1e-14)
    np.testing.assert_allclose(low_dense["R_gas_d"], high_dense["R_gas_d"], atol=1e-14)
    np.testing.assert_allclose(low_dense["R_gas_b"], 0.0, atol=1e-14)
    np.testing.assert_allclose(high_dense["R_gas_b"], 0.0, atol=1e-14)
    np.testing.assert_allclose(low_dense["full"], high_dense["full"], atol=1e-14)
    np.testing.assert_allclose(low_dense["projected"], high_dense["projected"], atol=1e-14)
    assert float(low_dense["energy"]) == pytest.approx(float(high_dense["energy"]), abs=1e-14)


def test_freeboard_energy_bridge_receives_bed_top_entrained_solid_enthalpy():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    bed_top = reactor.cells[-1]
    fb0 = reactor.freeboard_cells[0]
    propagate_explicit_freeboard_chain(reactor.cells, reactor.freeboard_cells)

    expected = bed_top._solid_upflow_rates()
    np.testing.assert_allclose(fb0.m_solid_in[:, S_CHAR], expected[:, S_CHAR])
    np.testing.assert_allclose(fb0.m_solid_in[:, S_ASH], expected[:, S_ASH])
    np.testing.assert_allclose(fb0.m_solid_in[:, S_VM], 0.0)
    np.testing.assert_allclose(fb0.m_solid_in[:, S_MOISTURE], 0.0)
    assert fb0.T_in_solid == pytest.approx(bed_top.T)

    energy_with_solid_bridge = float(fb0.calc_energy_balance())
    saved = fb0.m_solid_in.copy()
    fb0.m_solid_in.fill(0.0)
    try:
        energy_without_solid_bridge = float(fb0.calc_energy_balance())
    finally:
        fb0.m_solid_in[:, :] = saved
    assert energy_with_solid_bridge != pytest.approx(energy_without_solid_bridge)


def test_bed_top_energy_bridge_receives_freeboard_downflow_enthalpy():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    bed_top = reactor.cells[-1]
    fb0 = reactor.freeboard_cells[0]

    fb0.K_solid_ab[:, S_CHAR] = 0.0
    fb0.K_solid_ab[:, S_ASH] = 0.0
    fb0.K_solid_ab[0, S_CHAR] = 0.25
    fb0.m_solid[0, S_CHAR] = 1.2
    fb0.T = 975.0
    reactor._apply_all_bc_for_nr()

    expected_downflow = np.maximum(fb0._solid_downflow_rates(), 0.0)
    np.testing.assert_allclose(bed_top.m_solid_ab_in[:, S_CHAR], expected_downflow[:, S_CHAR])
    np.testing.assert_allclose(bed_top.m_solid_ab_in[:, S_VM], 0.0)
    np.testing.assert_allclose(bed_top.m_solid_ab_in[:, S_MOISTURE], 0.0)
    assert bed_top.T_solid_ab_in == pytest.approx(fb0.T)

    energy_with_downflow = float(bed_top.calc_energy_balance())
    saved_ab = bed_top.m_solid_ab_in.copy()
    saved_t = float(bed_top.T_solid_ab_in)
    bed_top.m_solid_ab_in.fill(0.0)
    try:
        energy_without_downflow = float(bed_top.calc_energy_balance())
    finally:
        bed_top.m_solid_ab_in[:, :] = saved_ab
        bed_top.T_solid_ab_in = saved_t
    assert energy_with_downflow != pytest.approx(energy_without_downflow)


def test_wirsum_trajectory_no_longer_depends_on_internal_subsegments():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    top = reactor.cells[-1]
    reactor._apply_all_bc_for_nr()

    area = math.pi * (float(reactor.config.D_bed) ** 2) / 4.0
    dh = float(reactor.config.H_freeboard) / max(int(reactor.config.n_freeboard_cells), 1)
    n_tot = max(float(np.sum(np.maximum(top.N_b + top.N_d, 0.0))), 1e-12)
    y_map = {sp: float(max((top.N_b + top.N_d)[j], 0.0) / n_tot) for j, sp in enumerate(GAS_SPECIES)}
    rho_g = gas_density_ideal(reactor.config.P, top.T, y_map)
    mu_g = gas_viscosity_power_law(top.T, 1.8e-5)
    d_p_class, char_frac_class, m_dot_class, _ = _build_size_class_bundles(
        d_p_classes=np.asarray(top.solid.d_p_classes),
        m_char_classes=np.asarray(top.m_solid[:, S_CHAR]),
        m_ash_classes=np.asarray(top.m_solid[:, 3]),
        rho_p=float(top.solid.rho_s),
        eps_b=float(top.eps_b),
        eps_d_void=float(top.eps_d_voidage),
        u_b=float(top.u_b),
        d_b=max(float(top.d_b), 1e-9),
        area=area,
    )
    vel, weights = _build_velocity_samples(1.53 * top.u_b, 0.6, 5)
    u_p_prev = []
    m_prev = []
    d_p_samples = []
    char_frac_samples = []
    for k, m_dot in enumerate(m_dot_class):
        for j, v in enumerate(vel):
            u_p_prev.append(float(v))
            m_prev.append(float(m_dot * weights[j]))
            d_p_samples.append(float(top.solid.d_p_classes[k]))
            char_frac_samples.append(float(char_frac_class[k]))

    ug0 = calc_u_gb(top.u_b, scale=1.0)
    beta = calc_beta_a(top.d_b)
    ug1 = float(top.u0) + (ug0 - float(top.u0)) * np.exp(-beta * dh)

    baseline = None
    for n_sub in (1, 2, 3, 6, 12, 24):
        up, mout, hold, ret, slower_u, slower_z, diag = _advance_particle_samples_analytical_wirsum(
            dh=dh,
            u_g_prev=ug0,
            u_g_next=ug1,
            rho_g=rho_g,
            rho_p=float(top.solid.rho_s),
            d_p=np.asarray(d_p_samples, dtype=np.float64),
            phi_s=float(top.solid.phi_s),
            mu_g=mu_g,
            u_p_samples_prev=np.asarray(u_p_prev, dtype=np.float64),
            m_dot_samples_prev=np.asarray(m_prev, dtype=np.float64),
            char_frac_samples=np.asarray(char_frac_samples, dtype=np.float64),
            n_subsegments=n_sub,
        )
        summary = (
            float(np.sum(mout)),
            float(np.sum(ret)),
            float(np.sum(hold)),
            float(np.max(up)),
            float(np.sum(slower_u)),
            float(np.sum(slower_z)),
        )
        if baseline is None:
            baseline = summary
            assert summary[0] + summary[1] > 0.0
            assert summary[2] > 0.0
            assert diag["sign_switch"] > 0
            assert diag["fallbacks"] >= 0
        else:
            np.testing.assert_allclose(summary, baseline, rtol=0.0, atol=1e-10)


def test_freeboard_bed_top_ejection_is_conservative_with_supplied_upflow():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    top = reactor.cells[-1]
    reactor._apply_all_bc_for_nr()

    area = math.pi * (float(reactor.config.D_bed) ** 2) / 4.0
    m_char = np.maximum(np.asarray(top._solid_upflow_rates()[:, S_CHAR], dtype=np.float64), 0.0)
    m_ash = np.maximum(np.asarray(top._solid_upflow_rates()[:, S_ASH], dtype=np.float64), 0.0)
    _, _, m_dot_class, m_dot_eject = _build_size_class_bundles(
        d_p_classes=np.asarray(top.solid.d_p_classes),
        m_char_classes=m_char,
        m_ash_classes=m_ash,
        rho_p=float(top.solid.rho_s),
        eps_b=float(top.eps_b),
        eps_d_void=float(top.eps_d_voidage),
        u_b=float(top.u_b),
        d_b=max(float(top.d_b), 1e-9),
        area=area,
    )

    expected = float(np.sum(m_char + m_ash))
    assert m_dot_eject == pytest.approx(expected, rel=1e-12, abs=1e-12)
    assert float(np.sum(m_dot_class)) == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_freeboard_bed_top_entrained_diameter_shrinks_with_char_conversion(monkeypatch: pytest.MonkeyPatch):
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    top = reactor.cells[-1]
    raw = np.array(top.solid.d_p_classes, dtype=np.float64, copy=True)
    monkeypatch.setattr(top, "_compute_char_conversion", lambda: 0.875)
    d_eff, x = effective_bed_top_entrained_d_p_classes(reactor)
    expected = raw * (1.0 - 0.875) ** (1.0 / 3.0)
    np.testing.assert_allclose(d_eff, expected, rtol=0.0, atol=1e-12)
    assert float(x) == pytest.approx(0.875, abs=1e-12)


def test_bed_top_size_age_proxy_amplifies_combustion_shrink_after_gasification():
    x_spm, x_comb, x_gas = _spm_shrink_proxy_from_reaction_partition(
        total_char_conversion_proxy=0.60,
        combustion_share_proxy=0.25,
    )

    assert x_comb == pytest.approx(0.15)
    assert x_gas == pytest.approx(0.45)
    assert x_spm == pytest.approx(0.15 / 0.55)
    assert x_spm > x_comb


def test_phase1_global_nr_pack_excludes_vm_and_moisture_dead_dofs():
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    assert np.isfinite(float(result["rms_scaled_final"]))
    assert n_solid_var(reactor.cells[0]) == reactor.cells[0].solid.n_size_classes * 2

    x = pack_reactor(reactor.cells)
    F0 = global_residual(x, reactor.cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        reactor.cells,
        ref_gas_mol_s=float(reactor.config.O2_feed + reactor.config.H2O_feed + reactor.config.N2_feed),
        ref_solid_kg_s=reactor.config.fuel_feed,
        ref_energy_W=reactor.config.fuel_feed * 20e6,
    )
    J, meta = build_jacobian_fd(
        x,
        F0,
        reactor.cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_structured",
    )

    assert meta["zero_rows"] == 0
    assert meta["zero_cols"] == 0


def test_phase1_global_nr_line_search_accepts_non_boundary_step():
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    accepted = result.get("nr_accepted_lambda_history", [])
    trial_counts = result.get("nr_line_search_trial_counts", [])
    assert any(lam is not None for lam in accepted)
    assert all(int(n) >= 1 for n in trial_counts)
    assert float(result["T_profile"][-1]) > 300.0


def test_block_tridiag_fd_matches_dense_fd_on_initialized_lu_state():
    reactor_dense = _build_initialized_lu_reactor(n_cells=3)
    cfg = reactor_dense.config
    x_dense = pack_reactor(reactor_dense.cells)
    F0_dense = global_residual(x_dense, reactor_dense.cells, reactor_dense._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        reactor_dense.cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    J_dense, meta_dense = build_jacobian_fd(
        x_dense,
        F0_dense,
        reactor_dense.cells,
        reactor_dense._apply_all_bc_for_nr,
        eq_scale,
        strategy="dense_fd",
    )
    reactor_block = _build_initialized_lu_reactor(n_cells=3)
    x_block = pack_reactor(reactor_block.cells)
    F0_block = global_residual(x_block, reactor_block.cells, reactor_block._apply_all_bc_for_nr)
    J_block, meta_block = build_jacobian_fd(
        x_block,
        F0_block,
        reactor_block.cells,
        reactor_block._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_fd",
    )

    diff = (J_dense - J_block).toarray()
    assert np.max(np.abs(diff)) < 1e-8
    assert meta_block["residual_cell_calls"] <= meta_dense["residual_cell_calls"]
    assert meta_dense["zero_cols"] >= 0 and meta_dense["zero_rows"] >= 0
    assert meta_block["zero_cols"] >= 0 and meta_block["zero_rows"] >= 0


def test_structured_bed_recycle_side_block_matches_dense_fd_on_initialized_lu_state():
    reactor_dense = _build_initialized_lu_reactor(n_cells=3)
    cfg = reactor_dense.config
    x_dense = pack_reactor(reactor_dense.cells)
    F0_dense = global_residual(x_dense, reactor_dense.cells, reactor_dense._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        reactor_dense.cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    structure = build_jacobian_structure(build_solver_graph(reactor_dense.cells), reactor_dense.cells)
    assert (0, len(reactor_dense.cells) - 1) in structure.side_block_pairs

    J_dense, meta_dense = build_jacobian_fd(
        x_dense,
        F0_dense,
        reactor_dense.cells,
        reactor_dense._apply_all_bc_for_nr,
        eq_scale,
        strategy="dense_fd",
    )
    reactor_structured = _build_initialized_lu_reactor(n_cells=3)
    x_structured = pack_reactor(reactor_structured.cells)
    F0_structured = global_residual(
        x_structured,
        reactor_structured.cells,
        reactor_structured._apply_all_bc_for_nr,
    )
    J_structured, meta_structured = build_jacobian_fd(
        x_structured,
        F0_structured,
        reactor_structured.cells,
        reactor_structured._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_structured",
        apply_local_bc_fn=reactor_structured._apply_local_bc_for_nr,
        affected_residual_cells_fn=reactor_structured._affected_nr_residual_cells,
    )

    diff = (J_dense - J_structured).toarray()
    assert np.max(np.abs(diff)) < 1e-8
    assert meta_structured["residual_cell_calls"] <= meta_dense["residual_cell_calls"]
    assert meta_structured["jacobian_structure"]["side_element_count"] >= 1


def test_block_tridiag_fd_local_bc_matches_full_bc_on_initialized_lu_state():
    reactor_full = _build_initialized_lu_reactor(n_cells=3)
    cfg = reactor_full.config
    x_full = pack_reactor(reactor_full.cells)
    F0_full = global_residual(x_full, reactor_full.cells, reactor_full._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        reactor_full.cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    J_full, _ = build_jacobian_fd(
        x_full,
        F0_full,
        reactor_full.cells,
        reactor_full._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_fd",
    )
    reactor_local = _build_initialized_lu_reactor(n_cells=3)
    x_local = pack_reactor(reactor_local.cells)
    F0_local = global_residual(x_local, reactor_local.cells, reactor_local._apply_all_bc_for_nr)
    J_local, _ = build_jacobian_fd(
        x_local,
        F0_local,
        reactor_local.cells,
        reactor_local._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_fd",
        apply_local_bc_fn=reactor_local._apply_local_bc_for_nr,
    )

    diff = (J_full - J_local).toarray()
    assert np.max(np.abs(diff)) < 1e-8


def test_block_tridiag_fd_forwards_local_bc_callback():
    reactor = _build_initialized_lu_reactor(n_cells=3)
    cfg = reactor.config
    x = pack_reactor(reactor.cells)
    F0 = global_residual(x, reactor.cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        reactor.cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    touched: list[int] = []

    def _track_local_bc(cell_idx: int) -> None:
        touched.append(int(cell_idx))
        reactor._apply_local_bc_for_nr(cell_idx)

    build_jacobian_fd(
        x,
        F0,
        reactor.cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_fd",
        apply_local_bc_fn=_track_local_bc,
    )

    assert touched
    assert set(touched).issubset({0, 1, 2})


def test_band_plus_side_elements_structured_explicit_path_ignores_local_callbacks():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 2
    cfg.n_freeboard_cells = 1
    cfg.explicit_freeboard_solver_graph_enabled = True
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    cfg = reactor.config
    solver_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(solver_cells)
    F0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        solver_cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    touched: list[int] = []
    affected_calls: list[int] = []

    def _track_local_bc(cell_idx: int) -> None:
        touched.append(int(cell_idx))
        reactor._apply_local_bc_for_nr(int(cell_idx))

    def _track_affected(cell_idx: int) -> tuple[int, ...]:
        affected_calls.append(int(cell_idx))
        return reactor._affected_nr_residual_cells(int(cell_idx))

    J_struct_with_cb, meta_struct = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="band_plus_side_elements_structured",
        apply_local_bc_fn=_track_local_bc,
        affected_residual_cells_fn=_track_affected,
    )
    J_struct_no_cb, _ = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="band_plus_side_elements_structured",
    )

    diff = (J_struct_no_cb - J_struct_with_cb).toarray()
    assert np.max(np.abs(diff)) < 1e-10
    assert touched == []
    assert affected_calls == []
    assert meta_struct["jacobian_structure"]["side_element_count"] > 0
    assert meta_struct["jacobian_structure"]["structure_validation_ok"] is True


def test_band_plus_side_elements_structured_closure_owned_path_ignores_local_callbacks():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    cfg.explicit_freeboard_solver_graph_enabled = False
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )

    solver_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(solver_cells)
    F0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        solver_cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    touched: list[int] = []
    affected_calls: list[int] = []

    def _track_local_bc(cell_idx: int) -> None:
        touched.append(int(cell_idx))
        reactor._apply_local_bc_for_nr(int(cell_idx))

    def _track_affected(cell_idx: int) -> tuple[int, ...]:
        affected_calls.append(int(cell_idx))
        return reactor._affected_nr_residual_cells(int(cell_idx))

    J_struct_with_cb, meta_struct = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="band_plus_side_elements_structured",
        apply_local_bc_fn=_track_local_bc,
        affected_residual_cells_fn=_track_affected,
    )
    J_struct_no_cb, _ = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="band_plus_side_elements_structured",
    )

    diff = (J_struct_no_cb - J_struct_with_cb).toarray()
    assert np.max(np.abs(diff)) < 1e-10
    assert touched == []
    assert affected_calls == []
    assert meta_struct["jacobian_structure"]["side_element_count"] > 0
    assert meta_struct["jacobian_structure"]["structure_validation_ok"] is True


def test_band_plus_side_elements_structured_matches_dense_fd_on_thesis_freeboard_graph():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    cfg = reactor.config
    solver_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(solver_cells)
    F0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        solver_cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    J_dense, meta_dense = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="dense_fd",
    )
    J_block, meta_block = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="band_plus_side_elements_structured",
    )

    # In closure-owned freeboard mode, dense_fd reflects a fully refreshed global
    # BC path, while structured Jacobian enforces the Hamel-style side-element
    # sparsity policy for NR linearization.
    assert J_block.shape == J_dense.shape
    assert np.all(np.isfinite(J_block.data))
    assert meta_block["residual_cell_calls"] <= meta_dense["residual_cell_calls"]
    assert meta_block["jacobian_structure"]["structure_validation_ok"] is True
    assert meta_block["jacobian_structure"]["side_element_count"] > 0


def test_thesis_freeboard_temperature_perturbation_changes_bed0_residual():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    solver_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(solver_cells)
    F0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    offsets = cell_offsets(solver_cells)
    bed0_slice = slice(offsets[0], offsets[1])
    base = F0[bed0_slice].copy()

    n_bed = len(reactor.cells)
    n_freeboard = len(reactor.freeboard_cells)
    idx_fb_last = n_bed + n_freeboard - 1

    def _max_bed0_delta_after_temperature_perturb(cell_idx: int, dT: float = 1.0) -> float:
        xp = x.copy()
        temp_var_idx = offsets[cell_idx + 1] - 1
        xp[temp_var_idx] += dT
        unpack_reactor(xp, solver_cells)
        reactor._apply_all_bc_for_nr()
        Fp = global_residual(xp, solver_cells, reactor._apply_all_bc_for_nr)
        delta = Fp[bed0_slice] - base
        unpack_reactor(x, solver_cells)
        reactor._apply_all_bc_for_nr()
        return float(np.max(np.abs(delta)))

    assert _max_bed0_delta_after_temperature_perturb(idx_fb_last) > 1e-8


def test_init_precalc_seeds_upper_bed_holdup_chain_for_thesis_mode():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)

    char_holdups = [
        float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0)))
        for cell in reactor.cells
    ]
    ash_holdups = [
        float(np.sum(np.maximum(cell.m_solid[:, S_ASH], 0.0)))
        for cell in reactor.cells
    ]

    assert char_holdups[0] > 0.0
    assert char_holdups[1] > 0.0
    assert char_holdups[2] > 0.0
    assert ash_holdups[1] > 0.0
    assert ash_holdups[2] > 0.0


def test_resolve_nr_init_strategy_prefers_vorabrechnung_by_default():
    assert _resolve_nr_init_strategy(None, None) == "vorabrechnung"
    assert _resolve_nr_init_strategy(None, 0) == "vorabrechnung"
    assert _resolve_nr_init_strategy(None, 1) == "vorabrechnung"
    # NR-only：无显式 init_strategy 时恒为 vorabrechnung；gs_warmup_steps / allow_legacy_gs 不参与解析
    assert _resolve_nr_init_strategy(None, 1, allow_legacy_gs=True) == "vorabrechnung"
    assert _resolve_nr_init_strategy("paper_vorab", None) == "vorabrechnung"


def test_default_nr_jacobian_strategy_is_block_for_bed_chain():
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    assert reactor._default_nr_jacobian_strategy() == "block_tridiag_structured"


def test_default_nr_jacobian_strategy_is_side_elements_sparse_for_thesis_side_blocks():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    assert reactor._default_nr_jacobian_strategy() == "band_plus_side_elements_structured"


def test_side_elements_sparse_fd_alias_matches_dense_fd_on_thesis_graph():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    cfg = reactor.config
    solver_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(solver_cells)
    F0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        solver_cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    J_dense, _ = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="dense_fd",
    )
    J_side, _ = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="side_elements_sparse_fd",
    )
    diff = (J_dense - J_side).toarray()
    assert np.max(np.abs(diff)) < 1e-8


def test_build_solver_graph_and_structure_report_wraparound_side_blocks():
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    solver_cells = reactor._solver_cells_for_nr()
    graph = build_solver_graph(solver_cells)
    structure = build_jacobian_structure(graph, solver_cells)

    assert graph.n_bed == 3
    assert graph.n_freeboard == 2
    assert graph.has_cyclone is True
    assert graph.has_return_leg is True
    assert graph.tail_start == 2
    assert graph.head_stop == 3
    assert (0, 3) in structure.side_block_pairs
    assert (0, 6) in structure.side_block_pairs
    assert (6, 4) in structure.side_block_pairs
    assert (0, 0) in structure.band_block_pairs
    assert (1, 0) in structure.band_block_pairs


def test_phase2_result_reports_structured_jacobian_metadata():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
    )

    assert result["nr_jacobian_strategy"] == "band_plus_side_elements_structured"
    assert result["nr_linear_solver_backend"] == "structured_direct"
    assert result["nr_linear_solver_backend_last"] in {"structured_direct", "sparse_direct_fallback"}
    assert isinstance(result["nr_jacobian_structure"], dict)
    assert result["nr_structure_validation_ok"] is True
    assert result["nr_band_block_count"] >= 1
    assert result["nr_side_element_count"] >= 0
    assert result["nr_schur_size"] >= 0
    assert result["nr_layout_audit"]
    layout_by_type = {entry["cell_type"]: entry for entry in result["nr_layout_audit"]}
    assert layout_by_type["freeboard"]["gas_mode"] == "single"
    assert layout_by_type["freeboard"]["gas_unknown_count"] == MAIN_NR_GAS_COUNT
    assert layout_by_type["freeboard"]["solid_unknown_count"] == 0
    assert layout_by_type["return_leg"]["gas_mode"] == "none"
    assert layout_by_type["return_leg"]["gas_unknown_count"] == 0
    assert all(entry["unknown_count"] == entry["residual_count"] for entry in result["nr_layout_audit"])
    assert "band_lu_s" in result["nr_timing"]
    assert "side_update_s" in result["nr_timing"]


def test_next_lambda_seed_adapts_from_last_successful_damping():
    assert abs(_next_lambda_seed(0.5, 0.25, 2, 0.5) - 0.3) < 1e-12
    assert abs(_next_lambda_seed(0.5, 0.125, 1, 0.3) - 0.25) < 1e-12
    assert abs(_next_lambda_seed(0.5, 0.5, 1, 0.9) - 0.5) < 1e-12


def test_build_equation_scales_use_packed_solid_dof_reference_for_holdup_transport():
    cell = Cell()
    cell.solid_state_model = "holdup_transport"
    cell.m_solid[0, S_CHAR] = 0.2
    cell.m_solid_zu[0, S_ASH] = 0.3
    cell.R_solid[0, S_CHAR] = 0.1
    cell.R_solid[0, S_ASH] = -0.05

    scale = build_equation_scales([cell], ref_gas_mol_s=10.0, ref_solid_kg_s=0.1, ref_energy_W=100.0)

    gas_width = n_var(cell) - n_solid_var(cell) - 1
    solid_slice = scale[gas_width : gas_width + n_solid_var(cell)]
    np.testing.assert_allclose(solid_slice, [0.3, 0.35])


def test_clip_dx_uses_cell_local_solid_reference_for_holdup_transport():
    cell = Cell()
    cell.solid_state_model = "holdup_transport"
    cell.m_solid[0, S_VM] = 2.0
    cell.m_solid_auf_in[0, S_VM] = 8.0
    gas_width = n_var(cell) - n_solid_var(cell) - 1
    dx = np.zeros(n_var(cell), dtype=np.float64)
    dx[gas_width : gas_width + n_solid_var(cell)] = 100.0

    clipped = _clip_dx(dx, [cell], ref_gas_mol_s=10.0, ref_solid_kg_s=0.1)

    np.testing.assert_allclose(clipped[gas_width : gas_width + n_solid_var(cell)], 2.0)


def test_solve_global_nr_reuses_backtracked_lambda_on_next_iter(monkeypatch):
    cell = Cell()
    x0 = pack_reactor([cell])
    n_total = len(x0)

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_equation_scales",
        lambda *args, **kwargs: np.ones(n_total),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_jacobian_fd",
        lambda *args, **kwargs: (sp.eye(n_total, format="csr"), {"residual_cell_calls": 0, "nnz": n_total}),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._clip_dx",
        lambda dx, cells, ref_gas_mol_s, ref_solid_kg_s, **kwargs: dx,
    )

    state: dict[str, object] = {"accepted_base": None, "line_search_calls": 0, "current_resid": 1.0}
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.pack_reactor",
        lambda cells: np.array(state["accepted_base"] if state["accepted_base"] is not None else x0, copy=True),
    )

    def _fake_global_residual(x, cells, apply_bc_fn):
        x = np.array(x, copy=True)
        if state["accepted_base"] is None:
            state["accepted_base"] = x.copy()
            return np.ones_like(x)

        base = np.array(state["accepted_base"], copy=True)
        lam = abs(float(x[-1] - base[-1])) / float(state["current_resid"])
        call_idx = int(state["line_search_calls"]) + 1
        state["line_search_calls"] = call_idx

        if call_idx == 1:
            assert abs(lam - 0.5) < 1e-12
            return np.full_like(x, 2.0)
        if call_idx == 2:
            assert abs(lam - 0.25) < 1e-12
            state["accepted_base"] = x.copy()
            state["current_resid"] = 0.8
            return np.full_like(x, 0.8)
        if call_idx == 3:
            # 第二轮应从上一轮接受步长自适应继承，而不是回到 0.5。
            assert abs(lam - 0.2) < 1e-12
            state["accepted_base"] = x.copy()
            state["current_resid"] = 0.7
            return np.full_like(x, 0.7)
        return np.full_like(x, 0.7)

    monkeypatch.setattr("src.solvers.global_nr_solver.global_residual", _fake_global_residual)

    result = solve_global_nr(
        cells=[cell],
        apply_bc_fn=lambda: None,
        max_iter=2,
        tol_rms=0.1,
        lambda_init=0.5,
        jacobian_lag_steps=2,
    )

    assert result["accepted_lambda_history"][:2] == [0.25, 0.2]
    assert result["line_search_trial_counts"][:2] == [2, 1]
    assert result["counts"]["line_search_backtracks"] == 1


def test_solve_global_nr_records_line_search_failure_diagnostics(monkeypatch):
    cell = Cell()
    x0 = pack_reactor([cell])
    n_total = len(x0)

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_equation_scales",
        lambda *args, **kwargs: np.ones(n_total),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_jacobian_fd",
        lambda *args, **kwargs: (
            sp.eye(n_total, format="csr"),
            {"residual_cell_calls": 0, "nnz": n_total, "zero_cols": 0, "zero_rows": 0},
        ),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._clip_dx",
        lambda dx, cells, ref_gas_mol_s, ref_solid_kg_s, **kwargs: dx,
    )

    state: dict[str, np.ndarray | None] = {"base": None}

    def _always_worse_residual(x, cells, apply_bc_fn):
        x_arr = np.array(x, copy=True)
        if state["base"] is None:
            state["base"] = x_arr.copy()
            return np.ones_like(x_arr)
        if np.allclose(x_arr, state["base"], atol=1e-15, rtol=0.0):
            return np.ones_like(x_arr)
        return np.full_like(x_arr, 2.0)

    monkeypatch.setattr("src.solvers.global_nr_solver.global_residual", _always_worse_residual)

    result = solve_global_nr(
        cells=[cell],
        apply_bc_fn=lambda: None,
        max_iter=1,
        tol_rms=0.1,
        lambda_init=0.5,
        n_damp_halvings=4,
    )

    assert result["accepted_lambda_history"] == [None]
    assert result["line_search_trial_counts"] == [4]
    assert result["counts"]["line_search_failures"] == 1
    assert result["counts"]["jacobian_zero_cols_last"] == 0
    assert result["counts"]["jacobian_zero_rows_last"] == 0
    diag = result["clip_history"][0]
    assert diag["line_search_failed"] is True
    assert diag["best_trial_lambda"] == pytest.approx(0.5)
    assert diag["best_trial_rms_scaled"] == pytest.approx(2.0)


def test_solve_global_nr_retries_once_after_line_search_failure(monkeypatch):
    cell = Cell()
    x0 = pack_reactor([cell])
    n_total = len(x0)

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_equation_scales",
        lambda *args, **kwargs: np.ones(n_total),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_jacobian_fd",
        lambda *args, **kwargs: (
            sp.eye(n_total, format="csr"),
            {"residual_cell_calls": 0, "nnz": n_total, "zero_cols": 0, "zero_rows": 0},
        ),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._clip_dx",
        lambda dx, cells, ref_gas_mol_s, ref_solid_kg_s, **kwargs: dx,
    )

    state: dict[str, object] = {"base": None, "trial_calls": 0}

    def _fail_then_improve(x, cells, apply_bc_fn):
        x_arr = np.array(x, copy=True)
        if state["base"] is None:
            state["base"] = x_arr.copy()
            return np.ones_like(x_arr)
        if np.allclose(x_arr, state["base"], atol=1e-15, rtol=0.0):
            return np.ones_like(x_arr)
        state["trial_calls"] = int(state["trial_calls"]) + 1
        if int(state["trial_calls"]) <= 4:
            return np.full_like(x_arr, 2.0)
        return np.full_like(x_arr, 0.5)

    monkeypatch.setattr("src.solvers.global_nr_solver.global_residual", _fail_then_improve)

    result = solve_global_nr(
        cells=[cell],
        apply_bc_fn=lambda: None,
        max_iter=2,
        tol_rms=0.1,
        lambda_init=0.5,
        n_damp_halvings=4,
    )

    assert result["accepted_lambda_history"] == [None, 0.5]
    assert result["line_search_trial_counts"] == [4, 1]
    assert result["counts"]["line_search_failures"] == 1
    assert result["counts"]["line_search_retries"] == 1


def test_apply_all_bc_for_nr_projects_orphan_vm_and_moisture():
    reactor = _build_initialized_lu_reactor(n_cells=3)
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)

    upper = reactor.cells[1]
    upper.m_solid[:, S_VM] = 0.05
    upper.m_solid[:, S_MOISTURE] = 0.03
    assert float(np.sum(upper.m_solid_in[:, S_VM] + upper.m_solid_zu[:, S_VM] + upper.m_solid_rez[:, S_VM])) == 0.0
    assert float(np.sum(upper.m_solid_in[:, S_MOISTURE] + upper.m_solid_zu[:, S_MOISTURE] + upper.m_solid_rez[:, S_MOISTURE])) == 0.0

    reactor._apply_all_bc_for_nr()

    assert np.allclose(reactor.cells[1].m_solid[:, S_VM], 0.0)
    assert np.allclose(reactor.cells[1].m_solid[:, S_MOISTURE], 0.0)


def test_reactor_propagates_ud_closure_to_cells():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    cfg.hydrodynamics_bubble_diameter_model = "hilligardt_ode"
    cfg.hydrodynamics_psi_b_strategy = "technical_distributor"
    cfg.hydrodynamics_lambda_strategy = "hamel_280"
    cfg.hydrodynamics_xi_strategy = "hamel_regime"
    cfg.hydrodynamics_bubble_velocity_strategy = "heinbockel_eq343"
    cfg.hydrodynamics_bubble_ode_strategy = "heinbockel_eq341"
    reactor = Reactor(cfg)

    assert {cell.u_d_closure for cell in reactor.cells} == {"backsolve_visible_epsb"}
    assert {cell.bubble_diameter_model for cell in reactor.cells} == {"hilligardt_ode"}
    assert {cell.psi_b_strategy for cell in reactor.cells} == {"technical_distributor"}
    assert {cell.lambda_strategy for cell in reactor.cells} == {"hamel_280"}
    assert {cell.xi_strategy for cell in reactor.cells} == {"hamel_regime"}
    assert {cell.bubble_velocity_strategy for cell in reactor.cells} == {"heinbockel_eq343"}
    assert {cell.bubble_ode_strategy for cell in reactor.cells} == {"heinbockel_eq341"}


def test_reactor_defaults_psi_b_to_wein_eq315():
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)

    assert cfg.hydrodynamics_psi_b_strategy == "wein_1992"
    assert {cell.psi_b_strategy for cell in reactor.cells} == {"wein_1992"}


def test_reactor_defaults_minor_species_to_gibbs_equilibrium():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)

    assert cfg.use_gibbs_minor is True
    assert {cell.use_gibbs_minor for cell in reactor.cells} == {True}


def test_thesis_mode_locks_hydrodynamics_chain_to_hamel_defaults():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    # Intentionally set conflicting closures to verify thesis-mode override.
    cfg.hydrodynamics_u_d_closure = "current"
    cfg.hydrodynamics_bubble_diameter_model = "mori_wen"
    cfg.hydrodynamics_psi_b_strategy = "technical_distributor"
    cfg.hydrodynamics_lambda_strategy = "current"
    cfg.hydrodynamics_xi_strategy = "current"
    cfg.hydrodynamics_bubble_velocity_strategy = "hilligardt_eq313"
    cfg.hydrodynamics_bubble_ode_strategy = "hilligardt_eq333"

    reactor = Reactor(cfg)

    assert {cell.u_d_closure for cell in reactor.cells} == {"wein_1992_eq312"}
    assert {cell.bubble_diameter_model for cell in reactor.cells} == {"hilligardt_ode"}
    assert {cell.psi_b_strategy for cell in reactor.cells} == {"wein_1992"}
    assert {cell.lambda_strategy for cell in reactor.cells} == {"hamel_280"}
    assert {cell.xi_strategy for cell in reactor.cells} == {"hamel_regime"}
    assert {cell.bubble_velocity_strategy for cell in reactor.cells} == {"heinbockel_eq343"}
    assert {cell.bubble_ode_strategy for cell in reactor.cells} == {"heinbockel_eq341"}
    assert "R8" in cfg.freeboard_enabled_reactions


def test_reactor_solve_defaults_to_global_nr(monkeypatch):
    reactor = Reactor(build_phase1_htw_lu_reactor_config())

    def _fake_global_nr(**kwargs):
        return {"solver_path": "global_nr", **kwargs}

    def _fake_gauss_seidel(**kwargs):
        raise AssertionError("gauss_seidel should not be used by default")

    monkeypatch.setattr(reactor, "_solve_global_nr", _fake_global_nr)
    monkeypatch.setattr(reactor, "_solve_gauss_seidel", _fake_gauss_seidel)

    result = reactor.solve(max_global_iter=7, tol_global=1e-3)

    assert result["solver_path"] == "global_nr"
    assert result["max_iter"] == 7
    assert result["tol"] == 1e-3


def test_reactor_rejects_non_nr_solver_even_in_thesis_mode():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    reactor = Reactor(cfg)

    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=5, tol_global=1e-2, solver="gauss_seidel")


def test_gauss_seidel_solver_is_disabled_unless_legacy_flag_enabled():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.allow_legacy_gs = False  # 与 Phase1 共享配置默认（为 GS 基线开放 legacy）区分
    reactor = Reactor(cfg)

    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=1, tol_global=1e-3, solver="gauss_seidel")


def test_gauss_seidel_solver_still_rejected_even_with_legacy_flag():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.allow_legacy_gs = True
    reactor = Reactor(cfg)

    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=2, tol_global=1e-2, solver="gauss_seidel")


def test_global_nr_outer_loop_forces_vorabrechnung_refresh(monkeypatch):
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)

    cell = reactor.cells[0]
    cell.u_mf = 1.0
    cell._vm_cache_valid = True
    seen_valid_before_refresh: list[bool] = []
    freeze_flags_during_inner: list[bool] = []
    inner_caps: list[int] = []

    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    def _hyd():
        cell.u_mf = 1.0
        cell._vorab_hydro_cache_valid = True

    monkeypatch.setattr(cell, "calc_hydrodynamics", _hyd)

    def _vorab(tau):
        seen_valid_before_refresh.append(bool(cell._vm_cache_valid))
        cell._vm_cache_valid = True
        cell._vorab_hydro_cache_valid = True

    monkeypatch.setattr(cell, "compute_vorabrechnung", _vorab)
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: (
            freeze_flags_during_inner.append(bool(cell._freeze_vorabrechnung_inner_nr)),
            inner_caps.append(int(kwargs["max_iter"])),
            {
                "converged": False,
                "n_iter": 2,
                "n_newton_iters_attempted": int(kwargs["max_iter"]),
                "rms_scaled_final": 0.2,
                "norm_history": [],
                "accepted_lambda_history": [0.25],
                "line_search_trial_counts": [2],
                "clip_history": [{"iter": 1, "accepted_lambda": 0.25, "line_search_trials": 2}],
            },
        )[2],
    )

    result = reactor._solve_global_nr(max_iter=8, tol=1.0, verbose=False)

    assert seen_valid_before_refresh == [False, False]
    assert freeze_flags_during_inner == [True, True]
    assert inner_caps == [6, 2]
    assert cell._freeze_vorabrechnung_inner_nr is False
    assert result["nr_vorabrechnung_policy"] == "outer_refresh_fixed_inner_sources"
    assert result["nr_accepted_lambda_history"] == [0.25, 0.25]
    assert result["nr_line_search_trial_counts"] == [2, 2]
    assert result["nr_clip_history"][0]["accepted_lambda"] == 0.25
    assert result["nr_inner_budget_total"] == 8
    assert result["nr_inner_budget_used"] == 8
    assert len(result["nr_vorabrechnung_signatures"]) == 2
    assert all(sig.startswith("vorab:v1:") for sig in result["nr_vorabrechnung_signatures"])
    assert result["nr_outer_history"][0]["vorabrechnung_signature"].startswith("vorab:v1:")


def test_global_nr_initialization_prepares_hydrodynamics_freeze_cache(monkeypatch):
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1
    reactor.config.nr_outer_iter_max = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    init_cache_called = {"ok": False}

    def _init_cache():
        init_cache_called["ok"] = True
        for c in reactor._solver_cells_for_nr():
            c._vorab_hydro_cache_valid = True

    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", _init_cache)

    def _fake_solve_global_nr(**kwargs):
        assert init_cache_called["ok"] is True
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 1e-6,
            "norm_history": [1e-6],
            "accepted_lambda_history": [1.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    monkeypatch.setattr("src.solvers.global_nr_solver.solve_global_nr", _fake_solve_global_nr)

    result = reactor._solve_global_nr(max_iter=3, tol=1.0, verbose=False)
    assert init_cache_called["ok"] is True
    assert bool(result["converged"]) is True


def test_thesis_inner_nr_uses_configurable_damping_knobs(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.nr_damping_halvings_thesis = 20
    cfg.nr_lambda_min_thesis = 1.0 / 32768.0
    cfg.nr_prefer_full_step_thesis = False
    cfg.nr_step_model_thesis = "lm"
    cfg.nr_lm_mu0_thesis = 2.0e-4
    cfg.nr_lm_mu_growth_thesis = 5.0
    cfg.nr_nonmonotone_enabled_thesis = True
    cfg.nr_nonmonotone_window_thesis = 7
    cfg.nr_nonmonotone_relax_thesis = 1.03
    cfg.nr_line_search_max_trials_thesis = 4
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1
    reactor.config.nr_outer_iter_max = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    for c in reactor._solver_cells_for_nr():
        c._vorab_hydro_cache_valid = True
    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", lambda: None)
    monkeypatch.setattr(reactor, "_refresh_vorabrechnung_sources_for_nr", lambda force: None)

    captured: dict[str, object] = {}

    def _fake_solve_global_nr(**kwargs):
        captured["n_damp_halvings"] = kwargs["n_damp_halvings"]
        captured["lambda_min"] = kwargs["lambda_min"]
        captured["prefer_full_step"] = kwargs["prefer_full_step"]
        captured["step_model"] = kwargs["step_model"]
        captured["lm_mu0"] = kwargs["lm_mu0"]
        captured["lm_mu_growth"] = kwargs["lm_mu_growth"]
        captured["nonmonotone_enabled"] = kwargs["nonmonotone_enabled"]
        captured["nonmonotone_window"] = kwargs["nonmonotone_window"]
        captured["nonmonotone_relax"] = kwargs["nonmonotone_relax"]
        captured["line_search_max_trials"] = kwargs["line_search_max_trials"]
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 1e-6,
            "norm_history": [1e-6],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [2],
            "clip_history": [],
        }

    monkeypatch.setattr("src.solvers.global_nr_solver.solve_global_nr", _fake_solve_global_nr)

    _ = reactor._solve_global_nr(max_iter=3, tol=1.0, verbose=False)

    assert captured["n_damp_halvings"] == 20
    assert float(captured["lambda_min"]) == pytest.approx(1.0 / 32768.0)
    assert bool(captured["prefer_full_step"]) is False
    assert str(captured["step_model"]) == "lm"
    assert float(captured["lm_mu0"]) == pytest.approx(2.0e-4)
    assert float(captured["lm_mu_growth"]) == pytest.approx(5.0)
    assert bool(captured["nonmonotone_enabled"]) is True
    assert int(captured["nonmonotone_window"]) == 7
    assert float(captured["nonmonotone_relax"]) == pytest.approx(1.03)
    assert int(captured["line_search_max_trials"]) == 4


def test_thesis_reaction_continuation_runs_warmup_then_full_stage(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.nr_reaction_continuation_stages_thesis = (0.25, 1.0)
    cfg.nr_reaction_continuation_height_m_thesis = cfg.H_bed
    cfg.nr_reaction_continuation_inner_iter_cap_thesis = 1
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1
    reactor.config.nr_outer_iter_max = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)
    for c in reactor._solver_cells_for_nr():
        c._vorab_hydro_cache_valid = True
    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", lambda: None)
    monkeypatch.setattr(reactor, "_refresh_vorabrechnung_sources_for_nr", lambda force: None)

    seen: list[float] = []

    def _fake_solve_global_nr(**kwargs):
        seen.append(float(kwargs["cells"][0].nr_reaction_rate_multiplier))
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 1e-6,
            "max_abs_scaled_final": 1e-6,
            "norm_history": [1e-6],
            "accepted_lambda_history": [1.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
            "timing": {},
            "counts": {},
        }

    monkeypatch.setattr("src.solvers.global_nr_solver.solve_global_nr", _fake_solve_global_nr)
    result = reactor._solve_global_nr(max_iter=3, tol=1.0, verbose=False)

    assert seen == [0.25, 1.0]
    assert reactor.cells[0].nr_reaction_rate_multiplier == pytest.approx(1.0)
    assert result["nr_reaction_continuation_active"] is True
    assert result["nr_reaction_continuation_stages"] == [0.25, 1.0]
    assert result["nr_reaction_continuation_target_cells"] == 1


def test_thesis_inner_nr_rejects_unknown_step_model():
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.nr_step_model_thesis = "bad_mode"
    reactor = Reactor(cfg)
    with pytest.raises(ValueError, match="nr_step_model_thesis"):
        _ = reactor.solve(
            solver="global_nr",
            max_global_iter=1,
            tol_global=1.0,
            nr_init_strategy="vorabrechnung",
            nr_jacobian_strategy="dense_fd",
        )


def test_thesis_inner_nr_passes_ptc_knobs(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.nr_step_model_thesis = "ptc"
    cfg.nr_ptc_alpha0_thesis = 0.5
    cfg.nr_ptc_alpha_growth_thesis = 3.0
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1
    reactor.config.nr_outer_iter_max = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)
    for c in reactor._solver_cells_for_nr():
        c._vorab_hydro_cache_valid = True
    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", lambda: None)
    monkeypatch.setattr(reactor, "_refresh_vorabrechnung_sources_for_nr", lambda force: None)

    captured: dict[str, object] = {}

    def _fake_solve_global_nr(**kwargs):
        captured["step_model"] = kwargs["step_model"]
        captured["ptc_alpha0"] = kwargs["ptc_alpha0"]
        captured["ptc_alpha_growth"] = kwargs["ptc_alpha_growth"]
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 1e-6,
            "norm_history": [1e-6],
            "accepted_lambda_history": [1.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    monkeypatch.setattr("src.solvers.global_nr_solver.solve_global_nr", _fake_solve_global_nr)
    _ = reactor._solve_global_nr(max_iter=3, tol=1.0, verbose=False)

    assert str(captured["step_model"]) == "ptc"
    assert float(captured["ptc_alpha0"]) == pytest.approx(0.5)
    assert float(captured["ptc_alpha_growth"]) == pytest.approx(3.0)


def test_thesis_inner_nr_passes_equilibration_knobs(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.nr_step_model_thesis = "equil_newton"
    cfg.nr_equil_iters_thesis = 4
    cfg.nr_equil_scale_clip_thesis = 1.0e2
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1
    reactor.config.nr_outer_iter_max = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)
    for c in reactor._solver_cells_for_nr():
        c._vorab_hydro_cache_valid = True
    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", lambda: None)
    monkeypatch.setattr(reactor, "_refresh_vorabrechnung_sources_for_nr", lambda force: None)

    captured: dict[str, object] = {}

    def _fake_solve_global_nr(**kwargs):
        captured["step_model"] = kwargs["step_model"]
        captured["equil_iters"] = kwargs["equil_iters"]
        captured["equil_scale_clip"] = kwargs["equil_scale_clip"]
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 1e-6,
            "norm_history": [1e-6],
            "accepted_lambda_history": [1.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    monkeypatch.setattr("src.solvers.global_nr_solver.solve_global_nr", _fake_solve_global_nr)
    _ = reactor._solve_global_nr(max_iter=3, tol=1.0, verbose=False)

    assert str(captured["step_model"]) == "equil_newton"
    assert int(captured["equil_iters"]) == 4
    assert float(captured["equil_scale_clip"]) == pytest.approx(1.0e2)


def test_thesis_mode_rejects_gs_warmup_init_strategy(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    cell = reactor.cells[0]
    cell.u_mf = 1.0
    monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: None)
    monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: setattr(cell, "_vorab_hydro_cache_valid", True))
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 0.2,
            "norm_history": [],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
        },
    )

    with pytest.raises(ValueError, match="gs_warmup_steps"):
        reactor._solve_global_nr(
            max_iter=6,
            tol=1.0,
            verbose=False,
            init_strategy="gs_warmup",
            gs_warmup_steps=2,
        )


def test_thesis_mode_single_shot_vorabrechnung_sources_not_recomputed_each_outer(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.thesis_vorab_sources_single_shot = True
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    bed_cells = list(reactor.cells)
    calls = {"vorab": 0}

    for cell in reactor._solver_cells_for_nr():
        cell.u_mf = 1.0
        monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: None)
    for cell in bed_cells:
        monkeypatch.setattr(
            cell,
            "compute_vorabrechnung",
            lambda tau, _calls=calls: _calls.__setitem__("vorab", _calls["vorab"] + 1),
        )

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: {
            "converged": False,
            "n_iter": 2,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 0.2,
            "norm_history": [],
            "accepted_lambda_history": [0.25],
            "line_search_trial_counts": [2],
            "clip_history": [],
        },
    )

    result = reactor._solve_global_nr(max_iter=8, tol=1.0, verbose=False)

    assert calls["vorab"] == len(bed_cells)
    assert result["nr_vorabrechnung_policy"] == "single_shot_sources_outer_refresh_hydrodynamics_fixed_inner"


def test_thesis_outer_refresh_does_not_reseed_side_blocks(monkeypatch):
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 2
    cfg.n_freeboard_cells = 1
    cfg.thesis_mode = True
    reactor = Reactor(cfg)

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0, 1010.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)
    monkeypatch.setattr(reactor, "_refresh_explicit_freeboard_transport_from_closure", lambda **_: None)

    side_seed_calls = {"n": 0}

    def _track_side_seed() -> None:
        side_seed_calls["n"] += 1

    monkeypatch.setattr(reactor, "_initialize_explicit_side_block_states", _track_side_seed)

    def _fake_init_hyd_freeze() -> None:
        for cell in reactor._solver_cells_for_nr():
            cell._vorab_hydro_cache_valid = True

    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", _fake_init_hyd_freeze)
    monkeypatch.setattr(reactor, "_initialize_thesis_single_shot_vorab_sources", lambda: None)

    for cell in reactor._solver_cells_for_nr():
        monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: None)
        monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: None)

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 0.2,
            "norm_history": [0.2],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
        },
    )

    _ = reactor._solve_global_nr(max_iter=8, tol=1.0, verbose=False)

    # Seed once in init/precalc; outer refresh must not overwrite side-block states.
    assert side_seed_calls["n"] == 1


def test_thesis_outer_refresh_skips_freeboard_closure_resync(monkeypatch):
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 2
    cfg.n_freeboard_cells = 1
    cfg.thesis_mode = True
    reactor = Reactor(cfg)

    assert reactor._use_explicit_freeboard_solver_graph() is True
    assert len(reactor.freeboard_cells) == 1
    reactor.freeboard_cells[0].m_solid.fill(1.2345)

    fb_calls = {"n": 0}

    def _track_freeboard_sync():
        fb_calls["n"] += 1
        reactor.freeboard_cells[0].m_solid.fill(999.0)
        return None

    monkeypatch.setattr(reactor, "_refresh_explicit_freeboard_transport_from_closure", _track_freeboard_sync)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    refresh_calls: list[dict[str, int | bool]] = []

    def _track_refresh(cells, *, force=False, refresh_sources=True):
        refresh_calls.append(
            {
                "n_cells": len(cells),
                "force": bool(force),
                "refresh_sources": bool(refresh_sources),
            }
        )

    monkeypatch.setattr("src.solvers.vorabrechnung.refresh_vorabrechnung_for_cells", _track_refresh)

    reactor._refresh_vorabrechnung_sources_for_nr(True)

    assert fb_calls["n"] == 0
    assert np.allclose(reactor.freeboard_cells[0].m_solid, 1.2345)
    assert refresh_calls == [
        {
            "n_cells": len(reactor.cells),
            "force": False,
            "refresh_sources": False,
        },
        {
            "n_cells": len(reactor.freeboard_cells) + len(reactor.side_cells),
            "force": False,
            "refresh_sources": False,
        },
    ]


def test_thesis_single_shot_expands_outer_slots_for_full_budget(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    cfg.thesis_vorab_sources_single_shot = True
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    for cell in reactor._solver_cells_for_nr():
        cell.u_mf = 1.0
        monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: None)
        monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: None)

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": int(kwargs["max_iter"]),
            "rms_scaled_final": 0.2,
            "norm_history": [0.2],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
        },
    )

    result = reactor._solve_global_nr(max_iter=20, tol=1.0, verbose=False)
    assert result["nr_outer_max"] == 20


def test_non_thesis_mode_rejects_gs_warmup_under_nr_only_policy(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = False
    cfg.allow_legacy_gs = False
    reactor = Reactor(cfg)

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0] * cfg.n_cells))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-3,
            "norm_history": [1e-3],
            "accepted_lambda_history": [1.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
        },
    )

    with pytest.raises(ValueError, match="gs_warmup_steps"):
        reactor._solve_global_nr(max_iter=3, tol=1.0, verbose=False, init_strategy="gs_warmup", gs_warmup_steps=2)


def test_thesis_mode_never_falls_back_to_gs_when_inner_nr_stalls(monkeypatch):
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.thesis_mode = True
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    cell = reactor.cells[0]
    cell.u_mf = 1.0
    monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: None)
    monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: setattr(cell, "_vorab_hydro_cache_valid", True))

    monkeypatch.setattr(
        reactor,
        "_solve_gauss_seidel",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("no GS fallback allowed in thesis_mode")),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.solve_global_nr",
        lambda **kwargs: {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 0.5,
            "norm_history": [0.5],
            "accepted_lambda_history": [None],
            "line_search_trial_counts": [8],
            "clip_history": [],
        },
    )

    result = reactor._solve_global_nr(
        max_iter=4,
        tol=1.0,
        verbose=False,
        init_strategy="vorabrechnung",
        gs_warmup_steps=0,
    )

    assert result["nr_init_strategy"] == "vorabrechnung"
    assert result["nr_gs_warmup_steps"] == 0
    assert result["converged"] is False


def test_global_nr_stops_when_inner_stalls_without_progress(monkeypatch):
    reactor = Reactor(build_phase1_htw_lu_reactor_config())
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1

    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)

    cell = reactor.cells[0]
    cell.u_mf = 1.0

    def _hyd():
        cell.u_mf = 1.0
        cell._vorab_hydro_cache_valid = True

    monkeypatch.setattr(cell, "calc_hydrodynamics", _hyd)
    monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: setattr(cell, "_vorab_hydro_cache_valid", True))

    calls = {"count": 0}

    def _fake_solve_global_nr(**kwargs):
        calls["count"] += 1
        if calls["count"] == 1:
            cell.T += 20.0
            return {
                "converged": False,
                "n_iter": 4,
                "n_newton_iters_attempted": 4,
                "rms_scaled_final": 0.02,
                "norm_history": [0.03, 0.02],
                "accepted_lambda_history": [0.5, 0.25],
                "line_search_trial_counts": [1, 2],
                "clip_history": [],
            }
        if calls["count"] == 2:
            cell.T += 0.1
            return {
                "converged": False,
                "n_iter": 1,
                "n_newton_iters_attempted": 1,
                "rms_scaled_final": 0.0205,
                "norm_history": [0.0205],
                "accepted_lambda_history": [None],
                "line_search_trial_counts": [8],
                "clip_history": [],
            }
        raise AssertionError("outer loop should stop after no-progress stall is detected")

    monkeypatch.setattr("src.solvers.global_nr_solver.solve_global_nr", _fake_solve_global_nr)

    result = reactor._solve_global_nr(max_iter=12, tol=1.0, verbose=False)

    assert calls["count"] == 2
    assert result["nr_outer_iters"] == 2
    assert result["n_iter"] == 3  # 2 norm points from iter-1 + 1 norm point from iter-2
    assert result["nr_inner_budget_used"] == 5
    assert result["nr_outer_history"][-1]["inner_converged"] is False
    assert result["nr_outer_history"][-1]["outer_aligned"] is False


def test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_major_gibbs_x0 = False
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=4,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_active"] is True
    assert result["freeboard_closure_mode"] == "explicit_cells_with_external_closure"
    assert len(reactor.freeboard_cells) == cfg.n_freeboard_cells
    assert all(cell.cell_type == "freeboard" for cell in reactor.freeboard_cells)
    assert result["freeboard_trajectory_coeff_model"] == cfg.freeboard_trajectory_coeff_model
    assert result["side_block_active"] is True
    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None
    assert len(result["bed_T_profile"]) == cfg.n_cells
    assert len(result["bed_u0_profile_m_s"]) == cfg.n_cells
    assert len(result["bed_eps_b_profile"]) == cfg.n_cells
    assert len(result["bed_eps_d_void_profile"]) == cfg.n_cells
    assert len(result["bed_bulk_solid_fraction_profile"]) == cfg.n_cells
    assert len(result["freeboard_T_profile"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_eps_b_profile"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_eps_d_void_profile"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_solid_holdup_char_profile_kg"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_solid_holdup_ash_profile_kg"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_char_reaction_source_profile_mol_s"]) == cfg.n_freeboard_cells
    assert len(result["T_profile"]) == cfg.n_cells + cfg.n_freeboard_cells
    assert "bed_exit_gas_dry" in result
    assert "exit_gas_dry" in result
    assert result["reactor_exit_T"] == result["cyclone_T"]
    assert result["freeboard_beta_A"] is not None
    assert result["freeboard_trajectory_model"] == cfg.freeboard_trajectory_model
    assert result["freeboard_trajectory_solver"] is not None
    assert isinstance(result["freeboard_trajectory_diag"], dict)
    assert result["freeboard_reaction_diag_labeling"] == "impl_ids"
    assert "R11" in result["freeboard_reaction_diag_impl"]
    assert "R11" in result["freeboard_reaction_diag_thesis"]
    assert set(result["freeboard_explicit_reaction_diag_impl"].keys()) == {"R1", "R2", "R3", "R4"}
    assert len(result["freeboard_explicit_reaction_diag_impl"]["R1"]) == cfg.n_freeboard_cells
    assert result["reaction_numbering"]["authority"] == "hamel_1999_thesis"
    assert result["reaction_numbering"]["impl_to_thesis"]["R12"] == "R6"
    topo = result["thesis_connectivity_topology"]
    assert topo is not None
    assert topo["mode"] == "hamel_minimal_connectivity"
    assert topo["freeboard_active"] is True
    assert topo["recycle_active"] is True
    block_kinds = {block["kind"] for block in topo["blocks"]}
    assert {"bed_cell", "freeboard_cell", "interface_block", "cyclone", "return_leg"}.issubset(block_kinds)
    block_map = {block["block_id"]: block for block in topo["blocks"]}
    assert block_map["cyclone_block"]["solver_coupling"] == "explicit_state"
    assert block_map["return_leg_block"]["solver_coupling"] == "explicit_state"
    edge_kinds = {edge["kind"] for edge in topo["edges"]}
    assert {"bubble_rupture_handoff", "cyclone_inlet", "external_zirkulation", "wake_solid_upflow"}.issubset(edge_kinds)
    internal_edges = [edge for edge in topo["edges"] if edge["kind"] in {"wake_solid_upflow", "internal_zirkulation"}]
    assert internal_edges
    assert all(edge["phase_scope"] == "suspension_phase_only" for edge in internal_edges)
    assert all(edge["from_block"].startswith("bed_cell_") and edge["to_block"].startswith("bed_cell_") for edge in internal_edges)
    assert result["freeboard_secondary_injection_applied"] is True
    assert result["freeboard_entrained_eject_char_ash_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_char_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_ash_kg_s"] >= 0.0
    assert result["freeboard_entrained_return_char_kg_s"] >= 0.0
    assert result["freeboard_entrained_return_ash_kg_s"] >= 0.0
    assert result["freeboard_cyclone_capture_char_kg_s"] >= 0.0
    assert result["freeboard_cyclone_capture_ash_kg_s"] >= 0.0
    assert result["freeboard_cyclone_recycle_candidate_char_ash_kg_s"] >= 0.0
    assert len(result["freeboard_entrained_char_kg_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_entrained_ash_kg_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_entrained_return_char_profile_kg_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_entrained_return_ash_profile_kg_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_u0_profile_m_s"]) == cfg.n_freeboard_cells
    np.testing.assert_allclose(result["freeboard_T_profile"], [cell.T for cell in reactor.freeboard_cells])
    assert "cyclone_T" in result
    assert "return_leg_T" in result
    assert len(result["freeboard_u_gb_profile_m_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_u_p_mean_profile_m_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_u_t_mean_profile_m_s"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_carry_ratio_profile"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_entrained_catalyst_density_kg_m3"]) == cfg.n_freeboard_cells
    np.testing.assert_allclose(
        result["freeboard_solid_holdup_char_profile_kg"],
        [float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) for cell in reactor.freeboard_cells],
    )
    np.testing.assert_allclose(
        result["freeboard_solid_holdup_ash_profile_kg"],
        [float(np.sum(np.maximum(cell.m_solid[:, S_ASH], 0.0))) for cell in reactor.freeboard_cells],
    )
    assert any(abs(v) > 0.0 for v in result["freeboard_char_reaction_source_profile_mol_s"])
    np.testing.assert_allclose(
        result["freeboard_entrained_char_kg_s"],
        [float(np.sum(np.maximum(cell._solid_upflow_rates()[:, S_CHAR], 0.0))) for cell in reactor.freeboard_cells],
    )
    np.testing.assert_allclose(
        result["freeboard_entrained_return_char_profile_kg_s"],
        [float(np.sum(np.maximum(cell._solid_downflow_rates()[:, S_CHAR], 0.0))) for cell in reactor.freeboard_cells],
    )
    np.testing.assert_allclose(
        result["freeboard_entrained_exit_char_kg_s"],
        float(np.sum(np.maximum(reactor.freeboard_cells[-1]._solid_upflow_rates()[:, S_CHAR], 0.0))),
    )
    np.testing.assert_allclose(
        result["freeboard_entrained_return_char_kg_s"],
        float(np.sum(np.maximum(reactor.freeboard_cells[0]._solid_downflow_rates()[:, S_CHAR], 0.0))),
    )
    np.testing.assert_allclose(
        result["freeboard_cyclone_capture_char_kg_s"],
        float(np.sum(np.maximum(reactor.cyclone_cell._solid_downflow_rates()[:, S_CHAR], 0.0))),
    )


def test_thesis_connectivity_topology_is_reported_for_bed_only_case():
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    cfg.vorab_major_gibbs_x0 = False
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=3,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy="block_tridiag_structured",
    )

    topo = result["thesis_connectivity_topology"]
    assert topo is not None
    assert topo["freeboard_active"] is False
    assert result["side_block_active"] is False
    block_ids = {block["block_id"] for block in topo["blocks"]}
    assert {"bed_cell_0", "cyclone_block", "return_leg_block", "system_exit"}.issubset(block_ids)
    block_map = {block["block_id"]: block for block in topo["blocks"]}
    assert block_map["cyclone_block"]["solver_coupling"] == "post_freeboard_separation"
    assert block_map["return_leg_block"]["solver_coupling"] == "boundary_condition"
    edge_kinds = {edge["kind"] for edge in topo["edges"]}
    assert {"external_zirkulation", "internal_zirkulation", "wake_solid_upflow"}.issubset(edge_kinds)
    backmix_edges = [edge for edge in topo["edges"] if edge["kind"] == "internal_zirkulation"]
    assert backmix_edges
    assert all(edge["solver_coupling"] == "explicit_state" for edge in backmix_edges)
    assert all(edge["mechanism"] == "mass_continuity_backmixing" for edge in backmix_edges)


def test_thesis_strict_major_gibbs_seed_no_longer_raises_on_phase2_init():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_major_gibbs_x0 = True
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )
    assert result.get("nr_init_strategy") == "vorabrechnung"


def test_reaction_numbering_mapping_helper_aggregates_impl_to_thesis():
    from src.core.reaction_numbering import map_impl_diag_to_thesis

    diag_impl = {"R11b": 1.0, "R11d": 2.0, "R12": 3.0, "R9": 4.0}
    out = map_impl_diag_to_thesis(diag_impl)

    assert out["R11"] == 3.0
    assert out["R6"] == 3.0
    assert out["_impl_only:R9"] == 4.0


def test_freeboard_analytical_wirsum_trajectory_mode_smoke():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = "analytical_wirsum"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_active"] is True
    assert result["freeboard_trajectory_model"] == "analytical_wirsum"
    assert result["freeboard_trajectory_solver"] == "wirsum_analytical"
    assert isinstance(result["freeboard_trajectory_diag"], dict)
    assert result["freeboard_entrained_eject_char_ash_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_char_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_ash_kg_s"] >= 0.0


def test_freeboard_analytical_wirsum_exact_hamel_coeff_mode_smoke():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = "analytical_wirsum"
    cfg.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_active"] is True
    assert result["freeboard_trajectory_model"] == "analytical_wirsum"
    assert result["freeboard_trajectory_coeff_model"] == "exact_hamel"
    assert result["freeboard_trajectory_solver"] == "wirsum_analytical_exact_hamel_coeffs"
    assert isinstance(result["freeboard_trajectory_diag"], dict)
    assert result["freeboard_entrained_eject_char_ash_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_char_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_ash_kg_s"] >= 0.0
    assert len(result["freeboard_solid_holdup_char_profile_kg"]) == cfg.n_freeboard_cells
    assert len(result["freeboard_solid_holdup_ash_profile_kg"]) == cfg.n_freeboard_cells
    assert result["freeboard_solid_holdup_char_profile_kg"][0] > 0.0
    assert sum(result["freeboard_solid_holdup_char_profile_kg"]) > 0.0
    assert sum(result["freeboard_solid_holdup_ash_profile_kg"]) > 0.0


def test_freeboard_global_projected_hamel_trajectory_mode_smoke():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = "global_projected_hamel"
    cfg.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_active"] is True
    assert result["freeboard_trajectory_model"] == "global_projected_hamel"
    assert result["freeboard_trajectory_solver"] == "global_projected_hamel"
    assert isinstance(result["freeboard_trajectory_diag"], dict)
    assert result["freeboard_trajectory_diag"].get("projection_crossings", 0) >= 0
    assert result["freeboard_entrained_eject_char_ash_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_char_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_ash_kg_s"] >= 0.0


def test_freeboard_force_balance_trajectory_mode_smoke():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = "force_balance"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_active"] is True
    assert result["freeboard_trajectory_model"] == "force_balance"
    assert result["freeboard_entrained_eject_char_ash_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_char_kg_s"] >= 0.0
    assert result["freeboard_entrained_exit_ash_kg_s"] >= 0.0


def test_global_projected_hamel_initialized_projection_reaches_at_least_as_high_as_analytical():
    analytical = _build_initialized_thesis_freeboard_reactor(n_bed=5, n_freeboard=8)
    analytical.config.freeboard_trajectory_model = "analytical_wirsum"
    analytical.config.freeboard_trajectory_coeff_model = "exact_hamel"
    analytical._refresh_explicit_freeboard_transport_from_closure()

    projected = _build_initialized_thesis_freeboard_reactor(n_bed=5, n_freeboard=8)
    projected.config.freeboard_trajectory_model = "global_projected_hamel"
    projected.config.freeboard_trajectory_coeff_model = "exact_hamel"
    projected._refresh_explicit_freeboard_transport_from_closure()

    analytical_up = [
        float(np.sum(np.maximum(cell._solid_upflow_rates()[:, S_CHAR], 0.0)))
        for cell in analytical.freeboard_cells
    ]
    projected_up = [
        float(np.sum(np.maximum(cell._solid_upflow_rates()[:, S_CHAR], 0.0)))
        for cell in projected.freeboard_cells
    ]
    analytical_last = max((i for i, v in enumerate(analytical_up) if v > 1e-12), default=-1)
    projected_last = max((i for i, v in enumerate(projected_up) if v > 1e-12), default=-1)
    projected_diag = (projected._last_explicit_freeboard_closure or {}).get("trajectory_diag", {})

    assert projected_last >= analytical_last
    assert projected_diag.get("fallback_no_event", 0) == 0
    assert projected_diag.get("fallback_delta_neg_time", 0) >= 0
    assert projected_diag.get("projection_crossings", 0) >= 0


def test_global_projected_hamel_reactor_smoke_activates_cyclone_inlet():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 5
    cfg.n_freeboard_cells = 8
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = "global_projected_hamel"
    cfg.freeboard_trajectory_coeff_model = "exact_hamel"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    top_up = float(np.sum(np.maximum(reactor.freeboard_cells[-1]._solid_upflow_rates()[:, S_CHAR], 0.0)))
    cyclone_in = float(np.sum(np.maximum(reactor.cyclone_cell.m_solid_auf_in[:, S_CHAR], 0.0)))

    assert result["freeboard_active"] is True
    assert result["freeboard_trajectory_model"] == "global_projected_hamel"
    assert result["freeboard_trajectory_diag"].get("projection_crossings", 0) >= 0
    assert result["freeboard_trajectory_diag"].get("fallback_delta_neg_time", 0) >= 0
    assert top_up >= 0.0
    assert cyclone_in >= 0.0


def test_freeboard_secondary_injection_xi_uses_global_reactor_coordinate():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_secondary_injection_xi = 0.60
    cfg.freeboard_secondary_O2_mol_s = 0.05
    cfg.freeboard_secondary_N2_mol_s = 0.05 * 3.76
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=3,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_secondary_injection_applied"] is True
    assert result["freeboard_secondary_injection_segment"] == 3


def test_phase2_local_refined_damped_builder_sets_bed_profile_and_damping():
    cfg = build_phase2_htw_lu_local_refined_damped_config()
    assert cfg.bed_dh_profile is not None
    assert len(cfg.bed_dh_profile) == cfg.n_cells
    assert sum(cfg.bed_dh_profile) == pytest.approx(cfg.H_bed)
    assert min(cfg.bed_dh_profile) < max(cfg.bed_dh_profile)
    assert cfg.nr_damping_halvings_thesis == 18
    assert cfg.nr_lambda_min_thesis == pytest.approx(1.0 / 16384.0)


def test_freeboard_secondary_local_refine_expands_profile_near_injection():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_secondary_injection_xi = 0.60
    cfg.freeboard_secondary_O2_mol_s = 0.05
    cfg.freeboard_secondary_N2_mol_s = 0.05 * 3.76
    cfg.freeboard_secondary_local_refine = 4
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=3,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_secondary_injection_applied"] is True
    assert result["freeboard_secondary_local_refine"] == 4
    assert len(result["freeboard_T_profile"]) == cfg.n_freeboard_cells + 3
    assert len(result["freeboard_axial_xi"]) == cfg.n_freeboard_cells + 3


def test_freeboard_secondary_distributed_mode_reports_metadata():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_secondary_injection_xi = 0.60
    cfg.freeboard_secondary_O2_mol_s = 0.05
    cfg.freeboard_secondary_N2_mol_s = 0.05 * 3.76
    cfg.freeboard_secondary_local_refine = 4
    cfg.freeboard_secondary_injection_mode = "distributed_uniform"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=3,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_secondary_injection_applied"] is True
    assert result["freeboard_secondary_injection_mode"] == "distributed_uniform"
    assert result["freeboard_secondary_local_refine"] == 4


def test_freeboard_secondary_observation_exposes_post_mix_pre_rxn_state():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_secondary_injection_xi = 0.60
    cfg.freeboard_secondary_O2_mol_s = 0.05
    cfg.freeboard_secondary_N2_mol_s = 0.05 * 3.76
    cfg.freeboard_secondary_local_refine = 4
    cfg.freeboard_secondary_injection_mode = "distributed_uniform"
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=3,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )

    assert result["freeboard_secondary_injection_applied"] is True
    assert len(result["freeboard_secondary_observation_xi"]) == 4
    assert len(result["freeboard_secondary_observation_post_mix_pre_rxn_T"]) == 4
    assert max(result["freeboard_secondary_observation_post_mix_pre_rxn_wet_gas"]["O2"]) > 0.0
    assert max(result["freeboard_secondary_observation_post_mix_pre_rxn_wet_gas"]["O2"]) >= max(
        result["freeboard_secondary_observation_post_rxn_wet_gas"]["O2"]
    )
