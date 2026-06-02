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
from src.solvers.convergence import OuterConvergence


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
    used_inner_budget: int


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
    outer_convergence: OuterConvergence | None = None,
    strict_check1_before_refresh: bool = False,
) -> OuterLoopResult:
    """Execute the outer Abgleich loop.

    Parameters
    ----------
    solve_inner_fn
        Signature: ``(inner_iter_cap, lambda_seed_outer, jacobian_mode, jacobian_lag) -> nr_result``.
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
    used_inner_budget = 0
    last_outer_rms = np.inf
    lambda_seed_outer = 0.5
    refresh_allowed = True
    # Check2 tolerance is an independent physical alignment criterion.
    # It must not be coupled to solve()'s tol_global (often O(1) for residual scaling),
    # otherwise outer alignment becomes trivially true.
    outer_checker = outer_convergence or OuterConvergence()

    for _outer in range(outer_max):
        remaining_inner_budget = total_inner_budget - used_inner_budget
        if remaining_inner_budget <= 0:
            break

        outer_iters = _outer + 1
        T_before_refresh = np.array([c.T for c in cells], dtype=np.float64)
        pre_vorab_signature = snapshot_signature_fn()
        refresh_skipped_for_check1_gate = bool(strict_check1_before_refresh and not refresh_allowed)

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
            remaining_outer_slots = max(outer_max - _outer, 1)
            reserve_per_future_outer = 2 if total_inner_budget >= 4 else 1
            inner_iter_cap = max(
                1,
                remaining_inner_budget - reserve_per_future_outer * max(remaining_outer_slots - 1, 0),
            )
            inner_iter_cap = min(int(inner_iter_cap), int(remaining_inner_budget))

        inner_started = perf_counter()
        nr_result = solve_inner_fn(
            inner_iter_cap,
            float(np.clip(lambda_seed_for_inner, 1.0 / 1024.0, 1.0)),
            jacobian_mode,
            jacobian_lag,
        )
        inner_solve_s_total += perf_counter() - inner_started

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
        raw_outer_aligned = bool(outer_checker.is_aligned(float(dT_outer), T_ref_scale))
        # Hamel flow semantics: evaluate Check2 only after Check1 passes.
        outer_matched = bool(inner_converged and raw_outer_aligned)
        outer_tol_effective = max(float(outer_checker.abs_tol_T), float(outer_checker.rel_tol) * max(T_ref_scale, 1.0))
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
                "inner_budget_remaining_after": int(max(total_inner_budget - used_inner_budget, 0)),
                "dT_outer": float(dT_outer),
                "inner_n_iter": int(nr_result.get("n_iter", 0)),
                "inner_n_newton_iters_attempted": int(
                    nr_result.get("n_newton_iters_attempted", nr_result.get("n_iter", 0))
                ),
                "inner_rms_scaled_final": curr_outer_rms,
                "inner_converged": inner_converged,
                "outer_aligned": outer_matched,
                "outer_aligned_raw": raw_outer_aligned,
                "outer_dT_tol_effective": float(outer_tol_effective),
            }
        )
        if strict_check1_before_refresh:
            refresh_allowed = bool(inner_converged)
        if outer_matched and inner_converged:
            outer_converged = True
            break

        inner_lambdas = nr_result.get("accepted_lambda_history") or []
        last_inner_lambda = inner_lambdas[-1] if inner_lambdas else None
        for lam in reversed(inner_lambdas):
            if lam is not None:
                lambda_seed_outer = float(np.clip(float(lam), 1.0 / 1024.0, 1.0))
                break

        no_inner_progress = int(nr_result.get("n_iter", 0)) <= 1
        rms_stalled = (
            np.isfinite(curr_outer_rms)
            and np.isfinite(last_outer_rms)
            and curr_outer_rms >= 0.999 * last_outer_rms
        )
        if no_inner_progress and rms_stalled:
            break
        last_outer_rms = curr_outer_rms

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
        used_inner_budget=used_inner_budget,
    )
