"""流程图：Zellenmodell / Newton–Raphson（Check1 内层）。

封装 ``solve_global_nr`` 调用与 Vorabrechnung 冻结开关生命周期。
"""

from __future__ import annotations

import functools
from typing import TYPE_CHECKING, Callable

import numpy as np

from src.solvers.convergence import InnerConvergence

if TYPE_CHECKING:
    from src.core.cell import Cell
    from src.core.reactor import Reactor, ReactorConfig


def _solve_inner_nr(
    reactor: "Reactor",
    solver_cells: list["Cell"],
    cfg: "ReactorConfig",
    tol_rms: float,
    inner_check1: InnerConvergence,
    verbose: bool,
    inner_iter_cap: int,
    lambda_seed_outer: float,
    j_mode: str,
    j_lag: int,
) -> dict:
    # Outer refresh snapshots hydro before Pre-Inner. Pre-Inner may mutate N/T/holdup
    # and gate on live hydro; resnapshot here so Inner solves F(x; h_entry).
    from src.solvers.vorabrechnung import prepare_inner_nr_hydrodynamics_freeze

    prepare_inner_nr_hydrodynamics_freeze(solver_cells, cfg=cfg)
    for cell in solver_cells:
        if not bool(getattr(cell, "_vorab_hydro_cache_valid", False)):
            raise RuntimeError("Vorabrechnung hydrodynamics freeze cache is not initialized for inner NR.")
    for cell in solver_cells:
        cell.enable_inner_nr_vorabrechnung_freeze()
    try:
        # 运行时从模块取符号，便于测试 monkeypatch ``src.solvers.global_nr_solver.solve_global_nr``
        from src.solvers import global_nr_solver as _gnr

        use_structured_direct = _gnr.structured_jacobian_uses_direct_linear_solver(j_mode)
        if (
            _gnr.is_band_plus_side_elements_jacobian_strategy(j_mode)
            and not bool(reactor._use_explicit_freeboard_solver_graph())
        ):
            use_structured_direct = False

        from src.solvers.wirsum_damping import wirsum_damped_newton_params_from_config

        wirsum = wirsum_damped_newton_params_from_config(cfg)
        thesis_halvings = int(wirsum.n_damp_halvings)
        thesis_lambda_min = float(wirsum.lambda_min)
        thesis_prefer_full_step = bool(wirsum.prefer_full_step)
        thesis_lambda_init = float(wirsum.lambda_init)
        thesis_step_model = str(getattr(cfg, "nr_step_model_thesis", "newton")).strip().lower()
        if thesis_step_model not in {"newton", "lm", "ptc", "equil_newton"}:
            raise ValueError(
                f"Unsupported nr_step_model_thesis={thesis_step_model!r}; expected 'newton', 'lm', 'ptc', or 'equil_newton'"
            )
        thesis_lm_mu0 = float(max(getattr(cfg, "nr_lm_mu0_thesis", 1.0e-4), 1.0e-12))
        thesis_lm_mu_growth = float(max(getattr(cfg, "nr_lm_mu_growth_thesis", 10.0), 1.0))
        thesis_ptc_alpha0 = float(max(getattr(cfg, "nr_ptc_alpha0_thesis", 1.0), 1.0e-12))
        thesis_ptc_alpha_growth = float(max(getattr(cfg, "nr_ptc_alpha_growth_thesis", 2.0), 1.0))
        thesis_equil_iters = int(max(getattr(cfg, "nr_equil_iters_thesis", 3), 0))
        thesis_equil_scale_clip = float(max(getattr(cfg, "nr_equil_scale_clip_thesis", 1.0e3), 1.0))
        thesis_ls_max_trials = int(max(getattr(cfg, "nr_line_search_max_trials_thesis", 0), 0))
        thesis_nonmono_enabled = bool(getattr(cfg, "nr_nonmonotone_enabled_thesis", False))
        thesis_nonmono_window = int(max(getattr(cfg, "nr_nonmonotone_window_thesis", 5), 1))
        thesis_nonmono_relax = float(max(getattr(cfg, "nr_nonmonotone_relax_thesis", 1.02), 1.0))
        thesis_enforce_gas_split = bool(
            getattr(cfg, "nr_enforce_gas_phase_split_convergence_thesis", True)
        )
        thesis_line_search_split_merit = bool(
            getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)
        )
        thesis_line_search_energy_merit = bool(
            getattr(cfg, "nr_line_search_energy_merit_thesis", False)
        )
        # Default False: Hamel outer_fixed hydro semantics (see ReactorConfig).
        thesis_refresh_hydro_on_accept = bool(
            getattr(cfg, "nr_refresh_hydrodynamics_on_accepted_step_thesis", False)
        )

        split_cap = float(getattr(cfg, "nr_inner_split_dominant_t_step_cap_K_thesis", 150.0))
        if split_cap <= 0.0:
            split_cap_val = None
        else:
            split_cap_val = split_cap

        hysteresis_on = bool(getattr(cfg, "nr_inner_t_cap_split_hysteresis_thesis", True))
        hold_split = bool(getattr(reactor, "_nr_tcap_hold_split_dominant", False)) if hysteresis_on else False
        release_ratio = float(
            getattr(cfg, "nr_inner_split_dominant_t_cap_release_ratio_thesis", 0.70)
        )

        base_dT_budget = float(getattr(cfg, "nr_inner_abs_dT_budget_K_thesis", 0.0) or 0.0)
        cont_dT_budget = float(
            getattr(cfg, "nr_inner_abs_dT_budget_check1_continuation_K_thesis", 0.0) or 0.0
        )
        if bool(getattr(reactor, "_nr_outer_refresh_skipped", False)) and cont_dT_budget > 0.0:
            if base_dT_budget > 0.0:
                abs_dT_budget = min(base_dT_budget, cont_dT_budget)
            else:
                abs_dT_budget = cont_dT_budget
        else:
            abs_dT_budget = base_dT_budget

        result = _gnr.solve_global_nr(
            cells=solver_cells,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
            affected_residual_cells_fn=reactor._affected_nr_residual_cells,
            connectivity_graph=reactor.connectivity_graph,
            ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
            ref_solid_kg_s=cfg.fuel_feed,
            ref_energy_W=cfg.fuel_feed * 20e6,
            max_iter=inner_iter_cap,
            tol_rms=tol_rms,
            inner_convergence=inner_check1,
            lambda_init=float(
                np.clip(
                    lambda_seed_outer if thesis_prefer_full_step else thesis_lambda_init,
                    thesis_lambda_min,
                    1.0,
                )
            ),
            n_damp_halvings=(thesis_halvings if bool(cfg.thesis_mode) else 12),
            lambda_min=(thesis_lambda_min if bool(cfg.thesis_mode) else 1.0 / 1024.0),
            jacobian_strategy=j_mode,
            linear_solver_backend=(
                "structured_direct" if use_structured_direct else "sparse_direct_fallback"
            ),
            jacobian_lag_steps=j_lag,
            allow_gd_fallback=not bool(cfg.thesis_mode),
            prefer_full_step=(thesis_prefer_full_step if bool(cfg.thesis_mode) else False),
            step_model=(thesis_step_model if bool(cfg.thesis_mode) else "newton"),
            lm_mu0=(thesis_lm_mu0 if bool(cfg.thesis_mode) else 1.0e-4),
            lm_mu_growth=(thesis_lm_mu_growth if bool(cfg.thesis_mode) else 10.0),
            ptc_alpha0=(thesis_ptc_alpha0 if bool(cfg.thesis_mode) else 1.0),
            ptc_alpha_growth=(thesis_ptc_alpha_growth if bool(cfg.thesis_mode) else 2.0),
            equil_iters=(thesis_equil_iters if bool(cfg.thesis_mode) else 0),
            equil_scale_clip=(thesis_equil_scale_clip if bool(cfg.thesis_mode) else 1.0e3),
            line_search_max_trials=(thesis_ls_max_trials if bool(cfg.thesis_mode) else 0),
            nonmonotone_enabled=(thesis_nonmono_enabled if bool(cfg.thesis_mode) else False),
            nonmonotone_window=(thesis_nonmono_window if bool(cfg.thesis_mode) else 5),
            nonmonotone_relax=(thesis_nonmono_relax if bool(cfg.thesis_mode) else 1.0),
            enforce_gas_phase_split_convergence=(
                thesis_enforce_gas_split if bool(cfg.thesis_mode) else False
            ),
            tol_gas_phase_split_max=tol_rms,
            line_search_use_gas_phase_split_merit=(
                thesis_line_search_split_merit if bool(cfg.thesis_mode) else False
            ),
            line_search_use_energy_merit=(
                thesis_line_search_energy_merit if bool(cfg.thesis_mode) else False
            ),
            refresh_hydrodynamics_on_accepted_step=(
                thesis_refresh_hydro_on_accept if bool(cfg.thesis_mode) else False
            ),
            line_search_fast_oxidation_projection=(
                bool(getattr(cfg, "nr_line_search_fast_oxidation_projection_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            line_search_per_phase_fastox=(
                bool(getattr(cfg, "nr_line_search_per_phase_fastox_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            line_search_per_phase_fastox_bed_limit=int(
                getattr(cfg, "nr_line_search_per_phase_fastox_bed_limit_thesis", 0)
            ),
            clip_same_phase_oxidizer_fuel_step=(
                bool(getattr(cfg, "nr_clip_same_phase_oxidizer_fuel_step_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            allow_oxidizer_into_fuel_step=(
                bool(getattr(cfg, "nr_clip_allow_oxidizer_into_fuel_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            line_search_syngas_collapse_guard=(
                bool(getattr(cfg, "nr_line_search_syngas_collapse_guard_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            line_search_syngas_collapse_ratio=float(
                getattr(cfg, "nr_line_search_syngas_collapse_ratio_thesis", 0.05)
            ),
            line_search_syngas_floor_ratio=float(
                getattr(cfg, "nr_line_search_syngas_floor_ratio_thesis", 0.0)
            ),
            line_search_refresh_aware_accept=(
                bool(getattr(cfg, "nr_line_search_refresh_aware_accept_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            line_search_split_constrained_step=(
                bool(getattr(cfg, "nr_line_search_split_constrained_step_thesis", False))
                if bool(cfg.thesis_mode)
                else False
            ),
            clip_solid_ftb_relative_floor=bool(
                getattr(cfg, "nr_clip_solid_ftb_relative_floor_thesis", False)
            ),
            inner_t_step_default_cap_K=float(
                getattr(cfg, "nr_inner_t_step_default_cap_K_thesis", 600.0)
            ),
            inner_split_dominant_t_step_cap_K=split_cap_val,
            inner_split_dominant_merit_ratio=float(
                getattr(cfg, "nr_inner_split_dominant_merit_ratio_thesis", 0.85)
            ),
            inner_t_cap_hold_split_dominant=bool(hold_split),
            inner_t_cap_split_release_ratio=float(release_ratio),
            inner_abs_dT_budget_K=(float(abs_dT_budget) if abs_dT_budget > 0.0 else None),
            record_inner_group_history=bool(
                getattr(cfg, "nr_inner_record_group_history_thesis", False)
            ),
            verbose=verbose,
        )
        if hysteresis_on:
            reactor._nr_tcap_hold_split_dominant = bool(
                result.get("inner_t_cap_hold_split_dominant", False)
            )
        return result
    finally:
        for cell in solver_cells:
            cell.disable_inner_nr_vorabrechnung_freeze()


def build_global_nr_inner_solve_fn(
    reactor: "Reactor",
    solver_cells: list["Cell"],
    cfg: "ReactorConfig",
    tol_rms: float,
    inner_check1: InnerConvergence,
    verbose: bool,
) -> Callable[[int, float, str, int], dict]:
    """构造 ``run_global_nr_outer_abgleich(..., solve_inner_fn=...)`` 所需的内层求解闭包。"""
    reactor._nr_tcap_hold_split_dominant = False
    return functools.partial(
        _solve_inner_nr,
        reactor,
        solver_cells,
        cfg,
        tol_rms,
        inner_check1,
        verbose,
    )
