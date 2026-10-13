"""global NR 测试共享 builder / harness（OPT-T01）。

从 ``test_global_nr_solver.py`` 提取，供参数化用例与 Jacobian 测试复用。
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from src.core.reactor import Reactor
from src.solvers.global_nr_solver import (
    build_equation_scales,
    build_jacobian_fd,
    global_residual,
    pack_reactor,
)
from src.solvers.vorabrechnung import (
    _major_elements_from_feeds,
    estimate_axial_T_profile,
    generate_initial_x0,
)
from src.thermal.devolatilization import devolatilization_rate_for_cell
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import (
    build_phase1_htw_lu_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
    build_vtt_minimal_loop_nr_acceptance_config,
)


def phase1_axial_t_profile(cfg) -> np.ndarray:
    return estimate_axial_T_profile(
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


def build_phase1_two_cell_reactor() -> tuple[Reactor, object, np.ndarray]:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    t_profile = phase1_axial_t_profile(cfg)
    return Reactor(cfg), cfg, t_profile


def call_generate_initial_x0_phase1(
    reactor: Reactor,
    cfg,
    t_profile: np.ndarray,
    *,
    use_hamel_major_gibbs_x0: bool = True,
    strict_hamel_major_gibbs_x0: bool = False,
) -> None:
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
        use_hamel_major_gibbs_x0=use_hamel_major_gibbs_x0,
        strict_hamel_major_gibbs_x0=strict_hamel_major_gibbs_x0,
    )


def major_gibbs_seed_raise(**kwargs: object) -> tuple[dict[str, float], dict[str, object]]:
    raise RuntimeError("boom")


def major_gibbs_seed_noncredible_large(**kwargs: object) -> tuple[dict[str, float], dict[str, object]]:
    return (
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
    )


def major_gibbs_seed_noncredible_finite(**kwargs: object) -> tuple[dict[str, float], dict[str, object]]:
    return (
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
    )


MAJOR_GIBBS_SEED_CASES: tuple[tuple[str, Callable[..., tuple], bool, str | None, float | None], ...] = (
    ("fallback_raise", major_gibbs_seed_raise, False, None, None),
    ("fallback_noncredible", major_gibbs_seed_noncredible_large, False, None, 1e6),
    ("strict_raise", major_gibbs_seed_raise, True, "boom", None),
    ("strict_noncredible", major_gibbs_seed_noncredible_finite, True, "non-credible seed in strict mode", None),
)


def build_initialized_lu_reactor(n_cells: int = 3) -> Reactor:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = n_cells
    reactor = Reactor(cfg)
    t_profile = phase1_axial_t_profile(cfg)
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


def build_initialized_thesis_freeboard_reactor(
    *,
    n_bed: int = 3,
    n_freeboard: int = 2,
) -> Reactor:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = int(n_bed)
    cfg.n_freeboard_cells = int(n_freeboard)
    cfg.vorab_major_gibbs_x0 = False
    cfg.major_gibbs_solver_mode = "hamel_reduced"
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    return reactor


def build_initialized_vtt_minimal_loop_reactor() -> Reactor:
    reactor = Reactor(build_vtt_minimal_loop_nr_acceptance_config())
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    return reactor


def phase2_first_cell_major_elements() -> tuple[dict[str, float], float, float]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    t_profile = phase1_axial_t_profile(cfg)
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


def setup_holdup_transport_reactor(n_cells: int) -> tuple[Reactor, object, np.ndarray]:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = n_cells
    reactor = Reactor(cfg)
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
    t_profile = phase1_axial_t_profile(cfg)
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    return reactor, cfg, t_profile


def run_holdup_generate_initial_x0(reactor: Reactor, cfg, t_profile: np.ndarray) -> None:
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


def assert_band_plus_side_elements_ignores_local_callbacks(
    *,
    n_bed: int,
    n_freeboard: int,
    explicit_freeboard_solver_graph_enabled: bool,
) -> None:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = int(n_bed)
    cfg.n_freeboard_cells = int(n_freeboard)
    cfg.explicit_freeboard_solver_graph_enabled = bool(explicit_freeboard_solver_graph_enabled)
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
    assert np.max(np.abs(diff)) < 1e-9
    assert touched == []
    assert affected_calls == []
    assert meta_struct["jacobian_structure"]["side_element_count"] > 0
    assert meta_struct["jacobian_structure"]["structure_validation_ok"] is True


def run_freeboard_trajectory_smoke(
    *,
    trajectory_model: str,
    coeff_model: str | None = None,
) -> dict[str, Any]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 2
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = trajectory_model
    if coeff_model is not None:
        cfg.freeboard_trajectory_coeff_model = coeff_model
    reactor = Reactor(cfg)
    return reactor.solve(
        solver="global_nr",
        max_global_iter=2,
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )


FREEBOARD_TRAJECTORY_SMOKE_CASES: tuple[tuple[str, str | None, str | None, bool], ...] = (
    ("analytical_wirsum", None, "wirsum_analytical_exact_hamel_coeffs", False),
    ("analytical_wirsum", "exact_hamel", "wirsum_analytical_exact_hamel_coeffs", True),
    ("global_projected_hamel", "exact_hamel", "global_projected_hamel", False),
    ("force_balance", None, None, False),
)


def patch_thesis_single_cell_outer_harness(monkeypatch, reactor: Reactor) -> None:
    monkeypatch.setattr("src.solvers.vorabrechnung.estimate_axial_T_profile", lambda **_: np.array([1000.0]))
    monkeypatch.setattr("src.solvers.vorabrechnung.generate_initial_x0", lambda **_: None)
    monkeypatch.setattr(reactor, "_apply_all_bc_for_nr", lambda: None)
    for c in reactor._solver_cells_for_nr():
        c._vorab_hydro_cache_valid = True
    monkeypatch.setattr(reactor, "_initialize_nr_hydrodynamics_freeze_cache", lambda: None)
    monkeypatch.setattr(reactor, "_refresh_vorabrechnung_sources_for_nr", lambda force: None)


def run_thesis_single_cell_inner_nr(
    monkeypatch,
    cfg,
    *,
    capture_keys: tuple[str, ...],
) -> dict[str, object]:
    reactor = Reactor(cfg)
    reactor.cells = reactor.cells[:1]
    reactor.config.n_cells = 1
    reactor.config.nr_outer_iter_max = 1
    patch_thesis_single_cell_outer_harness(monkeypatch, reactor)

    captured: dict[str, object] = {}

    def _fake_solve_global_nr(**kwargs):
        for key in capture_keys:
            captured[key] = kwargs[key]
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
    return captured


THESIS_INNER_NR_KNOB_CASES: tuple[tuple[str, dict[str, object], dict[str, object]], ...] = (
    (
        "damping_lm",
        {
            "nr_damping_halvings_thesis": 20,
            "nr_lambda_min_thesis": 1.0 / 32768.0,
            "nr_prefer_full_step_thesis": False,
            "nr_step_model_thesis": "lm",
            "nr_lm_mu0_thesis": 2.0e-4,
            "nr_lm_mu_growth_thesis": 5.0,
            "nr_nonmonotone_enabled_thesis": True,
            "nr_nonmonotone_window_thesis": 7,
            "nr_nonmonotone_relax_thesis": 1.03,
            "nr_line_search_max_trials_thesis": 4,
        },
        {
            "n_damp_halvings": 20,
            "lambda_min": 1.0 / 32768.0,
            "prefer_full_step": False,
            "step_model": "lm",
            "lm_mu0": 2.0e-4,
            "lm_mu_growth": 5.0,
            "nonmonotone_enabled": True,
            "nonmonotone_window": 7,
            "nonmonotone_relax": 1.03,
            "line_search_max_trials": 4,
        },
    ),
    (
        "ptc",
        {
            "nr_step_model_thesis": "ptc",
            "nr_ptc_alpha0_thesis": 0.5,
            "nr_ptc_alpha_growth_thesis": 3.0,
        },
        {
            "step_model": "ptc",
            "ptc_alpha0": 0.5,
            "ptc_alpha_growth": 3.0,
        },
    ),
    (
        "equil_newton",
        {
            "nr_step_model_thesis": "equil_newton",
            "nr_equil_iters_thesis": 4,
            "nr_equil_scale_clip_thesis": 1.0e2,
        },
        {
            "step_model": "equil_newton",
            "equil_iters": 4,
            "equil_scale_clip": 1.0e2,
        },
    ),
)
