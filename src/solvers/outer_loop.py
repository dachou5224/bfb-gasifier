"""Outer-loop orchestration for thesis-aligned global NR solve.

This module hosts the Abgleich loop that alternates:
1) Vorabrechnung refresh (hydrodynamics + source snapshots)
2) inner global NR solve with frozen per-cell Vorabrechnung caches

Keeping this logic outside Reactor helps preserve Hamel-style module boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Callable, List

import numpy as np

from src.core.cell import Cell
from src.solvers.convergence import OuterConvergence, OuterMeritEarlyStop


def resolve_outer_composite_merit(
    nr_result: dict,
    gas_metrics: dict[str, float],
    *,
    early_stop: OuterMeritEarlyStop,
) -> dict[str, float]:
    """Composite outer merit aligned with thesis inner line-search merit."""
    rms = float(nr_result.get("rms_scaled_final", np.inf))
    split = float(nr_result.get("rms_scaled_gas_phase_split_final", 0.0))
    if split <= 0.0 and gas_metrics.get("gas_phase_split_rms") is not None:
        split = float(gas_metrics["gas_phase_split_rms"])
    energy = float(nr_result.get("rms_scaled_energy_final", 0.0))
    merit = rms
    if early_stop.use_gas_phase_split_merit:
        merit = max(merit, split)
    if early_stop.use_energy_merit:
        merit = max(merit, energy)
    return {
        "rms_scaled": rms,
        "gas_phase_split_rms": split,
        "energy_rms": energy,
        "merit_composite": float(merit),
    }


@dataclass
class OuterLoopResult:
    converged_outer: bool
    outer_iters: int
    outer_history: list[dict]
    all_nr_history: list[float]
    all_nr_lambdas: list[float | None]
    all_nr_line_search_trials: list[int]
    all_nr_clip_history: list[dict]
    agg_nr_timing: dict[str, float]
    agg_nr_counts: dict[str, int]
    last_nr_result: dict
    outer_refresh_s_total: float
    inner_solve_s_total: float
    preinner_s_total: float
    used_inner_budget: int
    merit_early_stop: bool = False
    merit_early_stop_reason: str = ""


def run_global_nr_outer_abgleich(
    *,
    cells: List[Cell],
    outer_max: int,
    total_inner_budget: int,
    tol: float,
    jacobian_mode: str,
    jacobian_lag: int,
    tol_rms: float,
    refresh_fn: Callable[[bool], None],
    snapshot_signature_fn: Callable[[], str],
    solve_inner_fn: Callable[[int, float, str, int], dict],
    pre_inner_fn: Callable[..., dict[str, bool | float] | None] | None = None,
    post_inner_holdup_fn: Callable[[], bool] | None = None,
    pre_inner_on_continuation: bool = False,
    outer_convergence: OuterConvergence | None = None,
    gas_metrics_fn: Callable[[], dict[str, float]] | None = None,
    check2_bridge_fn: Callable[[], dict[str, bool | float]] | None = None,
    strict_check1_before_refresh: bool = False,
    check1_energy_regress_refresh_ratio: float = 0.0,
    outer_context_fn: Callable[..., None] | None = None,
    merit_early_stop: OuterMeritEarlyStop | None = None,
    hot_restart_after_y_co: bool = False,
    hot_restart_y_co_threshold: float = 0.14,
    hot_restart_bonus_inner_budget: int = 0,
    hot_restart_bonus_outer_slots: int = 0,
    y_co_top_fn: Callable[[], float] | None = None,
    hot_restart_clear_runtime_fn: Callable[[], None] | None = None,
) -> OuterLoopResult:
    """Execute the outer Abgleich loop.

    Parameters
    ----------
    solve_inner_fn
        Signature: ``(inner_iter_cap, lambda_seed_outer, jacobian_mode, jacobian_lag) -> nr_result``.
    check1_energy_regress_refresh_ratio
        If >1 and energy RMS rises by this factor vs the previous outer, re-allow
        Vorabrechnung refresh on the next outer even while Check1 is open.
    outer_context_fn
        Optional ``(*, refresh_skipped: bool)`` hook so the reactor can tighten
        inner |ΔT| budgets during Check1 continuation.
    hot_restart_after_y_co
        After bed-top ``y_CO`` enters bandwidth, reset outer gate/seed/merit state
        once on the next outer (skip_precalc-equivalent; does not re-run Vorab init).
    """
    all_nr_history: list[float] = []
    all_nr_lambdas: list[float | None] = []
    all_nr_line_search_trials: list[int] = []
    all_nr_clip_history: list[dict] = []
    agg_nr_timing: dict[str, float] = {}
    agg_nr_counts: dict[str, int] = {}
    outer_history: list[dict] = []

    outer_converged = False
    outer_iters = 0
    last_nr_result: dict = {}
    outer_refresh_s_total = 0.0
    inner_solve_s_total = 0.0
    preinner_s_total = 0.0
    used_inner_budget = 0
    last_outer_rms = np.inf
    lambda_seed_outer = 0.5
    refresh_allowed = True
    # Check2 tolerance is an independent physical alignment criterion.
    # It must not be coupled to solve()'s tol_global (often O(1) for residual scaling),
    # otherwise outer alignment becomes trivially true.
    outer_checker = outer_convergence or OuterConvergence()
    early_stop_cfg = merit_early_stop or OuterMeritEarlyStop(enabled=False)
    best_composite_merit = np.inf
    regression_streak = 0
    last_composite_merit = np.inf
    last_energy_rms = np.inf
    merit_early_stop_triggered = False
    merit_early_stop_reason = ""
    energy_regress_ratio = float(max(check1_energy_regress_refresh_ratio, 0.0))
    outer_max_eff = int(max(outer_max, 1))
    total_inner_budget_eff = int(max(total_inner_budget, 1))
    hot_restart_pending = False
    hot_restart_done = False
    y_co_thr = float(hot_restart_y_co_threshold)

    _outer = 0
    while _outer < outer_max_eff:
        hot_restarted_this_outer = False
        if hot_restart_pending and not hot_restart_done:
            hot_restarted_this_outer = True
            hot_restart_done = True
            hot_restart_pending = False
            # 等价 skip_precalc 二次 outer：复位门禁/种子/merit，不重跑 Vorab init
            refresh_allowed = True
            lambda_seed_outer = 0.5
            last_outer_rms = np.inf
            best_composite_merit = np.inf
            regression_streak = 0
            last_composite_merit = np.inf
            last_energy_rms = np.inf
            bonus_inner = int(max(hot_restart_bonus_inner_budget, 0))
            bonus_outer = int(max(hot_restart_bonus_outer_slots, 0))
            if bonus_inner > 0:
                total_inner_budget_eff += bonus_inner
            if bonus_outer > 0:
                outer_max_eff += bonus_outer
            if hot_restart_clear_runtime_fn is not None:
                hot_restart_clear_runtime_fn()

        remaining_inner_budget = total_inner_budget_eff - used_inner_budget
        if remaining_inner_budget <= 0:
            break

        outer_iters = _outer + 1
        T_before_refresh = np.array([c.T for c in cells], dtype=np.float64)
        pre_vorab_signature = snapshot_signature_fn()
        refresh_skipped_for_check1_gate = bool(strict_check1_before_refresh and not refresh_allowed)
        if outer_context_fn is not None:
            outer_context_fn(refresh_skipped=bool(refresh_skipped_for_check1_gate))

        if refresh_skipped_for_check1_gate:
            T_after_refresh = T_before_refresh
            vorab_signature = pre_vorab_signature
            dT_refresh = 0.0
            T_refresh_ref = float(np.max(np.abs(T_before_refresh))) if T_before_refresh.size else 1.0
            signature_changed_by_refresh = False
            refresh_seed_reset_threshold = max(10.0, 0.01 * max(T_refresh_ref, 1.0))
            seed_reset_by_refresh = False
            lambda_seed_for_inner = float(lambda_seed_outer)
        else:
            refresh_started = perf_counter()
            refresh_fn(True)
            outer_refresh_s_total += perf_counter() - refresh_started
            T_after_refresh = np.array([c.T for c in cells], dtype=np.float64)
            vorab_signature = snapshot_signature_fn()
            dT_refresh = float(np.max(np.abs(T_after_refresh - T_before_refresh))) if T_before_refresh.size else 0.0
            T_refresh_ref = float(np.max(np.abs(T_before_refresh))) if T_before_refresh.size else 1.0
            signature_changed_by_refresh = bool(pre_vorab_signature != vorab_signature)
            # 若 refresh 显著改变状态（典型表现：上一 outer 收敛后被 refresh 拉离），
            # 则不沿用上一 outer 的小 lambda_seed，避免把新线性化方向“饿死”在极小步长。
            refresh_seed_reset_threshold = max(10.0, 0.01 * max(T_refresh_ref, 1.0))
            seed_reset_by_refresh = bool(
                dT_refresh > refresh_seed_reset_threshold
                or (signature_changed_by_refresh and float(lambda_seed_outer) <= 0.1)
            )
            lambda_seed_for_inner = 0.5 if seed_reset_by_refresh else float(lambda_seed_outer)

        if refresh_skipped_for_check1_gate:
            # Bild 2.2 semantics: while Check1 fails, stay in inner Newton loop
            # and consume remaining budget without re-running Vorabrechnung.
            inner_iter_cap = int(remaining_inner_budget)
        else:
            remaining_outer_slots = max(outer_max_eff - _outer, 1)
            reserve_per_future_outer = 2 if total_inner_budget_eff >= 4 else 1
            inner_iter_cap = max(
                1,
                remaining_inner_budget - reserve_per_future_outer * max(remaining_outer_slots - 1, 0),
            )
            inner_iter_cap = min(int(inner_iter_cap), int(remaining_inner_budget))

        pre_inner_diag: dict[str, bool | float] | None = None
        if pre_inner_fn is not None and (
            not refresh_skipped_for_check1_gate or bool(pre_inner_on_continuation)
        ):
            preinner_started = perf_counter()
            try:
                pre_inner_diag = pre_inner_fn(
                    check1_continuation=bool(refresh_skipped_for_check1_gate)
                )
            except TypeError:
                # Backward-compatible call sites that ignore Check1 continuation.
                pre_inner_diag = pre_inner_fn()
            preinner_s_total += perf_counter() - preinner_started

        inner_started = perf_counter()
        nr_result = solve_inner_fn(
            inner_iter_cap,
            float(np.clip(lambda_seed_for_inner, 1.0 / 1024.0, 1.0)),
            jacobian_mode,
            jacobian_lag,
        )
        inner_solve_s_total += perf_counter() - inner_started

        postinner_holdup_reconciled = False
        if post_inner_holdup_fn is not None:
            postinner_holdup_reconciled = bool(post_inner_holdup_fn())

        last_nr_result = nr_result
        used_inner_budget += int(nr_result.get("n_newton_iters_attempted", nr_result.get("n_iter", 0)))
        all_nr_history.extend(nr_result.get("norm_history", []))
        all_nr_lambdas.extend(nr_result.get("accepted_lambda_history", []))
        all_nr_line_search_trials.extend(nr_result.get("line_search_trial_counts", []))
        for item in nr_result.get("clip_history", []):
            entry = dict(item)
            entry["outer_iter"] = outer_iters
            all_nr_clip_history.append(entry)
        for key, val in (nr_result.get("timing") or {}).items():
            agg_nr_timing[key] = float(agg_nr_timing.get(key, 0.0)) + float(val)
        for key, val in (nr_result.get("counts") or {}).items():
            if key == "jacobian_nnz_last":
                agg_nr_counts[key] = int(val)
            else:
                agg_nr_counts[key] = int(agg_nr_counts.get(key, 0)) + int(val)

        curr_outer_rms = float(nr_result.get("rms_scaled_final", np.inf))
        dT_outer = np.max(np.abs(np.array([c.T for c in cells], dtype=np.float64) - T_before_refresh))
        T_ref_scale = float(np.max(np.abs(T_before_refresh))) if T_before_refresh.size else 1.0
        inner_converged = bool(nr_result.get("converged"))
        gas_metrics: dict[str, float] = {}
        if gas_metrics_fn is not None and (
            outer_checker.enable_gas_phase_split_gate
            or outer_checker.enable_bed0_primary_gate
        ):
            gas_metrics = dict(gas_metrics_fn())
        merit_metrics = resolve_outer_composite_merit(
            nr_result,
            gas_metrics,
            early_stop=early_stop_cfg,
        )
        curr_composite_merit = float(merit_metrics["merit_composite"])
        if curr_composite_merit < best_composite_merit:
            best_composite_merit = curr_composite_merit
        raw_outer_aligned = bool(
            outer_checker.is_aligned(
                float(dT_outer),
                T_ref_scale,
                gas_phase_split_rms=gas_metrics.get("gas_phase_split_rms"),
                bed0_primary_rel_drift=gas_metrics.get("bed0_primary_rel_drift"),
            )
        )
        # Hamel flow semantics: evaluate Check2 only after Check1 passes.
        outer_matched = bool(inner_converged and raw_outer_aligned)
        outer_tol_effective = max(float(outer_checker.abs_tol_T), float(outer_checker.rel_tol) * max(T_ref_scale, 1.0))
        dt_check_pass = float(dT_outer) <= float(outer_tol_effective)
        gas_check_pass = True
        if outer_checker.enable_gas_phase_split_gate and gas_metrics.get("gas_phase_split_rms") is not None:
            gas_check_pass = float(gas_metrics["gas_phase_split_rms"]) <= float(
                outer_checker.max_gas_phase_split_rms
            )
        bed0_check_pass = True
        if outer_checker.enable_bed0_primary_gate and gas_metrics.get("bed0_primary_rel_drift") is not None:
            bed0_check_pass = float(gas_metrics["bed0_primary_rel_drift"]) <= float(
                outer_checker.max_bed0_primary_rel_drift
            )
        check2_bridge_diag: dict[str, bool | float] | None = None
        # Narrow bridge: only when Vorab gas/budget gates fail, not ΔT-only mismatch.
        if (
            inner_converged
            and not raw_outer_aligned
            and check2_bridge_fn is not None
            and (not gas_check_pass or not bed0_check_pass)
        ):
            check2_bridge_diag = dict(check2_bridge_fn())

        y_co_top_post = float("nan")
        if y_co_top_fn is not None:
            try:
                y_co_top_post = float(y_co_top_fn())
            except Exception:
                y_co_top_post = float("nan")
        elif pre_inner_diag is not None and pre_inner_diag.get("y_co_top") is not None:
            try:
                y_co_top_post = float(pre_inner_diag.get("y_co_top", float("nan")))
            except Exception:
                y_co_top_post = float("nan")

        hot_restart_armed_now = False
        if (
            bool(hot_restart_after_y_co)
            and (not hot_restart_done)
            and (not hot_restart_pending)
            and np.isfinite(y_co_top_post)
            and float(y_co_top_post) >= y_co_thr
        ):
            hot_restart_pending = True
            hot_restart_armed_now = True
            # 保证至少再跑 1 个 outer 执行热重启（bonus 在执行时再叠加）
            if (_outer + 1) >= outer_max_eff:
                outer_max_eff = _outer + 2

        outer_history.append(
            {
                "outer_iter": outer_iters,
                "vorabrechnung_signature_before_refresh": pre_vorab_signature,
                "vorabrechnung_signature": vorab_signature,
                "refresh_skipped_for_check1_gate": bool(refresh_skipped_for_check1_gate),
                "dT_refresh": float(dT_refresh),
                "signature_changed_by_refresh": bool(signature_changed_by_refresh),
                "refresh_seed_reset_threshold": float(refresh_seed_reset_threshold),
                "lambda_seed_reset_by_refresh": bool(seed_reset_by_refresh),
                "lambda_seed_used": float(np.clip(lambda_seed_for_inner, 1.0 / 1024.0, 1.0)),
                "inner_iter_cap": int(inner_iter_cap),
                "inner_budget_remaining_after": int(max(total_inner_budget_eff - used_inner_budget, 0)),
                "outer_hot_restarted": bool(hot_restarted_this_outer),
                "outer_hot_restart_armed": bool(hot_restart_armed_now),
                "outer_hot_restart_pending": bool(hot_restart_pending),
                "outer_hot_restart_done": bool(hot_restart_done),
                "dT_outer": float(dT_outer),
                "inner_n_iter": int(nr_result.get("n_iter", 0)),
                "inner_n_newton_iters_attempted": int(
                    nr_result.get("n_newton_iters_attempted", nr_result.get("n_iter", 0))
                ),
                "inner_rms_scaled_final": curr_outer_rms,
                "outer_energy_rms_final": float(merit_metrics["energy_rms"]),
                "outer_merit_composite": curr_composite_merit,
                "outer_merit_best_so_far": float(best_composite_merit),
                "inner_converged": inner_converged,
                "preinner_phase_split_enabled": bool((pre_inner_diag or {}).get("enabled", False)),
                "preinner_split_rms_before": float((pre_inner_diag or {}).get("split_rms_before", 0.0)),
                "preinner_split_rms_after": float((pre_inner_diag or {}).get("split_rms_after", 0.0)),
                "preinner_merit_before": float((pre_inner_diag or {}).get("merit_before", 0.0)),
                "preinner_merit_after": float((pre_inner_diag or {}).get("merit_after", 0.0)),
                "preinner_bed0_accepted": bool((pre_inner_diag or {}).get("bed0_accepted", False)),
                "preinner_bed1_accepted": bool((pre_inner_diag or {}).get("bed1_split_accepted", False)),
                "preinner_bed2_accepted": bool((pre_inner_diag or {}).get("bed2_accepted", False)),
                "preinner_global_gate_pass": bool((pre_inner_diag or {}).get("global_gate_pass", False)),
                "preinner_split_only_fallback": bool(
                    (pre_inner_diag or {}).get("split_only_fallback", False)
                ),
                "preinner_solid_rms_before": float((pre_inner_diag or {}).get("solid_rms_before", 0.0)),
                "preinner_solid_rms_after": float((pre_inner_diag or {}).get("solid_rms_after", 0.0)),
                "preinner_energy_rms_before": float((pre_inner_diag or {}).get("energy_rms_before", 0.0)),
                "preinner_energy_rms_after": float((pre_inner_diag or {}).get("energy_rms_after", 0.0)),
                "preinner_merit_driven": bool((pre_inner_diag or {}).get("merit_driven", False)),
                "preinner_skipped_all_tiers": bool((pre_inner_diag or {}).get("skipped_all_tiers", False)),
                "preinner_holdup_reconciled": bool((pre_inner_diag or {}).get("holdup_reconciled", False)),
                "preinner_holdup_flux_mismatch_before": float(
                    (pre_inner_diag or {}).get("holdup_flux_mismatch_before", 0.0)
                ),
                "preinner_holdup_skipped_check1_continuation": bool(
                    (pre_inner_diag or {}).get("holdup_skipped_check1_continuation", False)
                ),
                "preinner_run_holdup_reconcile": bool(
                    ((pre_inner_diag or {}).get("tier_plan") or {}).get("run_holdup_reconcile", False)
                ),
                "preinner_split_shock": bool(
                    ((pre_inner_diag or {}).get("tier_plan") or {}).get("split_shock", False)
                ),
                "postinner_holdup_reconciled": bool(postinner_holdup_reconciled),
                "preinner_bed1_energy_t_accepted": bool(
                    (pre_inner_diag or {}).get("bed1_energy_t_accepted", False)
                ),
                "upper_dense_absorb_enabled": bool(
                    (pre_inner_diag or {}).get("upper_dense_absorb_enabled", False)
                ),
                "upper_dense_absorb_accepted": bool(
                    (pre_inner_diag or {}).get("upper_dense_absorb_accepted", False)
                ),
                "upper_dense_absorb_beds": list(
                    (pre_inner_diag or {}).get("upper_dense_absorb_beds") or []
                ),
                "upper_dense_absorb_added_CO_mol_s": float(
                    (pre_inner_diag or {}).get("upper_dense_absorb_added_CO_mol_s", 0.0)
                ),
                "upper_dense_absorb_added_H2_mol_s": float(
                    (pre_inner_diag or {}).get("upper_dense_absorb_added_H2_mol_s", 0.0)
                ),
                "upper_co_local2x2_enabled": bool(
                    (pre_inner_diag or {}).get("upper_co_local2x2_enabled", False)
                ),
                "upper_co_local2x2_accepted": bool(
                    (pre_inner_diag or {}).get("upper_co_local2x2_accepted", False)
                ),
                "upper_co_local2x2_beds": list(
                    (pre_inner_diag or {}).get("upper_co_local2x2_beds") or []
                ),
                "upper_co_bridge_frozen": bool(
                    (pre_inner_diag or {}).get("upper_co_bridge_frozen", False)
                ),
                "y_co_top": float((pre_inner_diag or {}).get("y_co_top", float("nan"))),
                "y_co_top_post": float(y_co_top_post) if np.isfinite(y_co_top_post) else float("nan"),
                "trace_o2_rate_gate_armed": bool(
                    (pre_inner_diag or {}).get("trace_o2_rate_gate_armed", False)
                ),
                "trace_o2_rate_gate_newly_armed": bool(
                    (pre_inner_diag or {}).get("trace_o2_rate_gate_newly_armed", False)
                ),
                "trace_o2_rate_gate_y_co_top": float(
                    (pre_inner_diag or {}).get("trace_o2_rate_gate_y_co_top", float("nan"))
                ),
                "bed_top_energy_passthrough_armed": bool(
                    (pre_inner_diag or {}).get("bed_top_energy_passthrough_armed", False)
                ),
                "bed_top_energy_passthrough_newly_armed": bool(
                    (pre_inner_diag or {}).get("bed_top_energy_passthrough_newly_armed", False)
                ),
                "bed_top_energy_passthrough_y_co_top": float(
                    (pre_inner_diag or {}).get("bed_top_energy_passthrough_y_co_top", float("nan"))
                ),
                "inner_group_history": list(nr_result.get("inner_group_history") or []),
                "outer_aligned": outer_matched,
                "outer_aligned_raw": raw_outer_aligned,
                "outer_dT_tol_effective": float(outer_tol_effective),
                "gas_phase_split_rms": float(gas_metrics.get("gas_phase_split_rms", 0.0)),
                "gas_phase_split_max_abs": float(gas_metrics.get("gas_phase_split_max_abs", 0.0)),
                "bed0_primary_rel_drift": float(gas_metrics.get("bed0_primary_rel_drift", 0.0)),
                "outer_gas_split_gate_enabled": bool(outer_checker.enable_gas_phase_split_gate),
                "outer_gas_split_tol": float(outer_checker.max_gas_phase_split_rms),
                "outer_bed0_primary_gate_enabled": bool(outer_checker.enable_bed0_primary_gate),
                "outer_bed0_primary_tol": float(outer_checker.max_bed0_primary_rel_drift),
                "check2_bridge_enabled": bool((check2_bridge_diag or {}).get("enabled", False)),
                "check2_bridge_passes": float((check2_bridge_diag or {}).get("bridge_passes", 0.0)),
                "check2_bridge_aligned_after": bool((check2_bridge_diag or {}).get("aligned", False)),
                "outer_dt_check_pass": bool(dt_check_pass),
                "outer_gas_check_pass": bool(gas_check_pass),
                "outer_bed0_primary_check_pass": bool(bed0_check_pass),
            }
        )
        if outer_matched and inner_converged:
            outer_converged = True
            break

        inner_lambdas = nr_result.get("accepted_lambda_history") or []
        last_inner_lambda = inner_lambdas[-1] if inner_lambdas else None
        for lam in reversed(inner_lambdas):
            if lam is not None:
                lambda_seed_outer = float(np.clip(float(lam), 1.0 / 1024.0, 1.0))
                break

        inner_accepted = any(
            lam is not None for lam in (nr_result.get("accepted_lambda_history") or [])
        )
        no_inner_progress = int(nr_result.get("n_iter", 0)) <= 1
        rms_stalled = (
            np.isfinite(curr_outer_rms)
            and np.isfinite(last_outer_rms)
            and curr_outer_rms >= 0.999 * last_outer_rms
        )
        merit_stalled = (
            np.isfinite(curr_composite_merit)
            and np.isfinite(last_composite_merit)
            and curr_composite_merit
            >= (1.0 - float(early_stop_cfg.plateau_rel_tol)) * last_composite_merit
        )
        if strict_check1_before_refresh:
            energy_rms_now = float(merit_metrics.get("energy_rms", curr_outer_rms))
            energy_regressed = (
                energy_regress_ratio > 1.0
                and np.isfinite(energy_rms_now)
                and np.isfinite(last_energy_rms)
                and last_energy_rms < 1.0e90
                and energy_rms_now > energy_regress_ratio * max(float(last_energy_rms), 1.0e-12)
            )
            if inner_converged:
                refresh_allowed = True
            elif refresh_skipped_for_check1_gate and rms_stalled and (
                no_inner_progress or not inner_accepted
            ):
                # Check1 continuation 无线搜索接受步：下一轮允许 Vorab refresh 重建线性化。
                refresh_allowed = True
            elif energy_regressed:
                # Large energy regression under frozen Vorab: rebuild hydro/K before more ΔT.
                refresh_allowed = True
            else:
                refresh_allowed = False
            last_energy_rms = energy_rms_now
        if early_stop_cfg.enabled and outer_iters >= int(early_stop_cfg.min_outer_iters):
            if curr_composite_merit > best_composite_merit * (
                1.0 + float(early_stop_cfg.regression_rel_tol)
            ):
                regression_streak += 1
            else:
                regression_streak = 0
            if regression_streak >= int(early_stop_cfg.regression_patience):
                if hot_restart_pending and not hot_restart_done:
                    pass
                else:
                    merit_early_stop_triggered = True
                    merit_early_stop_reason = "composite_merit_regression"
                    break
        plateau_stalled = merit_stalled if early_stop_cfg.enabled else rms_stalled
        if no_inner_progress and plateau_stalled and inner_accepted:
            # 已武装热重启时勿 plateau 早停，否则 skip_precalc 等价续算永远跑不到
            if hot_restart_pending and not hot_restart_done:
                pass
            else:
                if early_stop_cfg.enabled:
                    merit_early_stop_triggered = True
                    merit_early_stop_reason = "composite_merit_plateau"
                break
        last_outer_rms = curr_outer_rms
        last_composite_merit = curr_composite_merit
        _outer += 1

    return OuterLoopResult(
        converged_outer=outer_converged,
        outer_iters=outer_iters,
        outer_history=outer_history,
        all_nr_history=all_nr_history,
        all_nr_lambdas=all_nr_lambdas,
        all_nr_line_search_trials=all_nr_line_search_trials,
        all_nr_clip_history=all_nr_clip_history,
        agg_nr_timing=agg_nr_timing,
        agg_nr_counts=agg_nr_counts,
        last_nr_result=last_nr_result,
        outer_refresh_s_total=outer_refresh_s_total,
        inner_solve_s_total=inner_solve_s_total,
        preinner_s_total=preinner_s_total,
        used_inner_budget=used_inner_budget,
        merit_early_stop=bool(merit_early_stop_triggered),
        merit_early_stop_reason=str(merit_early_stop_reason),
    )
