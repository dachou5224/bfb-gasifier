"""流程图：Zellenmodell / Newton–Raphson（Check1 内层）。

封装 ``solve_global_nr`` 调用与 Vorabrechnung 冻结开关生命周期。
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from src.solvers.convergence import InnerConvergence

if TYPE_CHECKING:
    from src.core.cell import Cell
    from src.core.reactor import Reactor, ReactorConfig


def build_global_nr_inner_solve_fn(
    reactor: "Reactor",
    solver_cells: list["Cell"],
    cfg: "ReactorConfig",
    tol_rms: float,
    inner_check1: InnerConvergence,
    verbose: bool,
) -> Callable[[int, float, str, int], dict]:
    """构造 ``run_global_nr_outer_abgleich(..., solve_inner_fn=...)`` 所需的内层求解闭包。"""

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        for cell in solver_cells:
            if not bool(getattr(cell, "_vorab_hydro_cache_valid", False)):
                raise RuntimeError("Vorabrechnung hydrodynamics freeze cache is not initialized for inner NR.")
        for cell in solver_cells:
            cell.enable_inner_nr_vorabrechnung_freeze()
        try:
            # 运行时从模块取符号，便于测试 monkeypatch ``src.solvers.global_nr_solver.solve_global_nr``
            from src.solvers import global_nr_solver as _gnr

            use_structured_direct = str(j_mode) in {"block_tridiag_structured", "band_plus_side_elements_structured"}
            if (
                str(j_mode) == "band_plus_side_elements_structured"
                and not bool(reactor._use_explicit_freeboard_solver_graph())
            ):
                # Closure-owned freeboard path: keep Hamel side-element Jacobian
                # structure, but use sparse-direct linear solve for robustness.
                use_structured_direct = False

            thesis_halvings = max(1, int(getattr(cfg, "nr_damping_halvings_thesis", 14)))
            thesis_lambda_min = float(max(getattr(cfg, "nr_lambda_min_thesis", 1.0 / 4096.0), 1e-12))
            thesis_prefer_full_step = bool(getattr(cfg, "nr_prefer_full_step_thesis", True))
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

            return _gnr.solve_global_nr(
                cells=solver_cells,
                apply_bc_fn=reactor._apply_all_bc_for_nr,
                apply_local_bc_fn=reactor._apply_local_bc_for_nr,
                affected_residual_cells_fn=reactor._affected_nr_residual_cells,
                ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
                ref_solid_kg_s=cfg.fuel_feed,
                ref_energy_W=cfg.fuel_feed * 20e6,
                max_iter=inner_iter_cap,
                tol_rms=tol_rms,
                inner_convergence=inner_check1,
                lambda_init=float(np.clip(lambda_seed_outer, 1.0 / 1024.0, 1.0)),
                # Thesis path keeps Hamel-style damped Newton (full step first, then
                # binary damping), but allows deeper halving before declaring failure.
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
                verbose=verbose,
            )
        finally:
            for cell in solver_cells:
                cell.disable_inner_nr_vorabrechnung_freeze()

    return _solve_inner
