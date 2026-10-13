"""Result assembly helpers for Reactor NR orchestration.

Keep `src/core/reactor.py` focused on flow orchestration (Config -> Vorabrechnung -> Cell Model -> NR),
while this module builds reporting payloads from solved states.
"""

from __future__ import annotations

import numpy as np

from src.core.cell import S_ASH, S_CHAR


def compute_carbon_conversion(reactor, cell_solid_outflow_component_fn) -> float:
    bot = reactor.cells[0]
    top = reactor.cells[-1]
    m_char_in = float(
        np.sum(
            np.maximum(
                bot.m_solid_zu[:, S_CHAR] + bot.m_solid_in[:, S_CHAR],
                0.0,
            )
        )
    )
    if reactor._use_explicit_side_block_cells() and reactor.cyclone_cell is not None:
        m_char_out = float(cell_solid_outflow_component_fn(reactor.cyclone_cell, S_CHAR, direction="up"))
    else:
        top_up_char = float(cell_solid_outflow_component_fn(top, S_CHAR, direction="up"))
        recycle_char = float(np.sum(np.maximum(bot.m_solid_rez[:, S_CHAR], 0.0)))
        m_char_out = max(top_up_char - recycle_char, 0.0)
    carbon_conv = 1.0 - (m_char_out / max(m_char_in, 1e-12))
    return float(np.clip(carbon_conv, 0.0, 1.0))


def refresh_explicit_freeboard_closure_for_result(reactor) -> dict:
    """Report explicit freeboard closure freshness without mutating NR state."""
    if not reactor._use_explicit_freeboard_solver_graph():
        return {
            "freeboard_result_closure_refreshed": False,
            "freeboard_result_closure_refresh_reason": "not_explicit_freeboard_graph",
        }
    if not getattr(reactor, "freeboard_cells", None):
        return {
            "freeboard_result_closure_refreshed": False,
            "freeboard_result_closure_refresh_reason": "no_freeboard_cells",
        }

    current_bed_top_initial = float(np.sum(np.maximum(reactor.cells[-1]._solid_upflow_rates()[:, S_CHAR], 0.0)))
    previous_fb = getattr(reactor, "_last_explicit_freeboard_closure", None)
    previous_bed_top = (
        float(previous_fb.get("bed_top_up_char_kg_s", current_bed_top_initial))
        if isinstance(previous_fb, dict)
        else current_bed_top_initial
    )
    # Do not refresh/sync the closure here. The NR solver has just produced a
    # self-consistent state vector; mutating explicit freeboard/side solid states
    # during result assembly invalidates the final residual metrics and can turn
    # side-element solid DOFs into apparent post-solve residual spikes.
    final_bed_top = current_bed_top_initial
    closure_gap = current_bed_top_initial - previous_bed_top
    return {
        "freeboard_result_closure_refreshed": False,
        "freeboard_result_boundary_refreshed": False,
        "freeboard_result_closure_refresh_reason": "diagnostic_only_no_state_mutation",
        "freeboard_pre_result_refresh_bed_top_up_char_kg_s": previous_bed_top,
        "freeboard_pre_result_refresh_current_bed_top_up_char_kg_s": current_bed_top_initial,
        "freeboard_current_bed_top_up_char_kg_s": final_bed_top,
        "freeboard_pre_result_refresh_gap_current_bed_top_kg_s": closure_gap,
        "freeboard_post_result_refresh_bed_top_up_char_kg_s": previous_bed_top,
        "freeboard_post_result_refresh_gap_current_bed_top_kg_s": closure_gap,
        "freeboard_result_closure_refresh_iters": 0,
    }


def build_exit_summary(
    reactor,
    resolve_axial_heat_loss_distribution_fn,
    cell_solid_outflow_component_fn,
) -> dict:
    from src.solvers.result_exit_summary import build_exit_summary_segments

    return build_exit_summary_segments(
        reactor,
        resolve_axial_heat_loss_distribution_fn=resolve_axial_heat_loss_distribution_fn,
        cell_solid_outflow_component_fn=cell_solid_outflow_component_fn,
    )


def finalize_global_nr_result(
    reactor,
    nr_result: dict,
    resolve_axial_heat_loss_distribution_fn,
    cell_solid_outflow_component_fn,
) -> dict:
    freeboard_result_refresh = refresh_explicit_freeboard_closure_for_result(reactor)
    carbon_conv = compute_carbon_conversion(reactor, cell_solid_outflow_component_fn)
    exit_summary = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=resolve_axial_heat_loss_distribution_fn,
        cell_solid_outflow_component_fn=cell_solid_outflow_component_fn,
    )
    return {
        "converged": nr_result["converged"],
        "converged_outer": nr_result.get("converged_outer"),
        "converged_inner_nr": nr_result.get("converged_inner_nr"),
        "converged_fully": nr_result.get("converged_outer") and nr_result.get("converged_inner_nr"),
        "rms_scaled_final": nr_result.get("rms_scaled_final"),
        "rms_scaled_gas_final": nr_result.get("rms_scaled_gas_final"),
        "rms_scaled_solid_final": nr_result.get("rms_scaled_solid_final"),
        "rms_scaled_energy_final": nr_result.get("rms_scaled_energy_final"),
        "rms_scaled_gas_combined_final": nr_result.get("rms_scaled_gas_combined_final"),
        "rms_scaled_gas_phase_split_final": nr_result.get("rms_scaled_gas_phase_split_final"),
        "rms_scaled_component_max_final": nr_result.get("rms_scaled_component_max_final"),
        "max_abs_scaled_gas_final": nr_result.get("max_abs_scaled_gas_final"),
        "max_abs_scaled_solid_final": nr_result.get("max_abs_scaled_solid_final"),
        "max_abs_scaled_energy_final": nr_result.get("max_abs_scaled_energy_final"),
        "max_abs_scaled_gas_combined_final": nr_result.get("max_abs_scaled_gas_combined_final"),
        "max_abs_scaled_gas_phase_split_final": nr_result.get("max_abs_scaled_gas_phase_split_final"),
        "max_abs_scaled_final": nr_result.get("max_abs_scaled_final"),
        "n_iter": nr_result["n_iter"],
        "residual": nr_result.get("residual", 0.0),
        "norm_history": nr_result.get("norm_history", []),
        "nr_init_strategy": nr_result.get("nr_init_strategy"),
        "nr_init_s_total": nr_result.get("nr_init_s_total"),
        "nr_vorabrechnung_s": nr_result.get("nr_vorabrechnung_s"),
        "nr_vorabrechnung_policy": nr_result.get("nr_vorabrechnung_policy"),
        "nr_gs_warmup_steps": nr_result.get("nr_gs_warmup_steps", 0),
        "nr_outer_max": nr_result.get("nr_outer_max"),
        "nr_outer_iters": nr_result.get("nr_outer_iters"),
        "nr_outer_history": nr_result.get("nr_outer_history"),
        "nr_outer_hot_restart": nr_result.get("nr_outer_hot_restart"),
        "nr_vorabrechnung_signatures": nr_result.get("nr_vorabrechnung_signatures", []),
        "nr_inner_tol_rms": nr_result.get("nr_inner_tol_rms"),
        "nr_jacobian_strategy": nr_result.get("nr_jacobian_strategy", nr_result.get("jacobian_strategy")),
        "nr_layout_audit": nr_result.get("nr_layout_audit", []),
        "nr_jacobian_structure": nr_result.get("jacobian_structure"),
        "nr_linear_solver_backend": nr_result.get("nr_linear_solver_backend", nr_result.get("linear_solver_backend")),
        "nr_linear_solver_backend_last": nr_result.get("linear_solver_backend_last"),
        "nr_side_element_count": (nr_result.get("jacobian_structure") or {}).get("side_element_count"),
        "nr_band_block_count": (nr_result.get("jacobian_structure") or {}).get("band_block_count"),
        "nr_structure_validation_ok": (nr_result.get("jacobian_structure") or {}).get("structure_validation_ok"),
        "nr_linear_fallback_used": bool((nr_result.get("counts") or {}).get("structured_fallbacks", 0)),
        "nr_schur_size": (nr_result.get("counts") or {}).get("nr_schur_size_last"),
        "nr_jacobian_lag_steps": nr_result.get("nr_jacobian_lag_steps", 1),
        "nr_reaction_continuation_active": nr_result.get("nr_reaction_continuation_active", False),
        "nr_reaction_continuation_stages": nr_result.get("nr_reaction_continuation_stages", [1.0]),
        "nr_reaction_continuation_target_cells": nr_result.get("nr_reaction_continuation_target_cells", 0),
        "nr_gs_warmup_s": nr_result.get("nr_gs_warmup_s", 0.0),
        "nr_outer_refresh_s_total": nr_result.get("nr_outer_refresh_s_total"),
        "nr_inner_solve_s_total": nr_result.get("nr_inner_solve_s_total"),
        "nr_preinner_s_total": nr_result.get("nr_preinner_s_total"),
        "nr_inner_budget_total": nr_result.get("nr_inner_budget_total"),
        "nr_inner_budget_used": nr_result.get("nr_inner_budget_used"),
        "nr_total_s": nr_result.get("nr_total_s"),
        "nr_timing": nr_result.get("timing"),
        "nr_counts": nr_result.get("counts"),
        "nr_accepted_lambda_history": nr_result.get("accepted_lambda_history", []),
        "nr_line_search_trial_counts": nr_result.get("line_search_trial_counts", []),
        "nr_clip_history": nr_result.get("clip_history", []),
        "carbon_conv": carbon_conv,
        **exit_summary,
        **freeboard_result_refresh,
    }
