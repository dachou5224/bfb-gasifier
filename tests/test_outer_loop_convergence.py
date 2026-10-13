from __future__ import annotations

from dataclasses import dataclass

from src.solvers.convergence import InnerConvergence, OuterConvergence, OuterMeritEarlyStop
from src.solvers.outer_loop import resolve_outer_composite_merit, run_global_nr_outer_abgleich


@dataclass
class _FakeCell:
    T: float


def test_inner_convergence_rtol_zero_matches_absolute_atol():
    ic = InnerConvergence(atol=0.01, rtol=0.0)
    assert ic.is_converged(0.009, 1e6) is True
    assert ic.is_converged(0.011, 0.0) is False


def test_outer_convergence_alignment_threshold():
    checker = OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5)
    assert checker.is_aligned(dT_outer_max=0.2, T_ref_scale=1000.0) is True
    assert checker.is_aligned(dT_outer_max=2.0, T_ref_scale=1000.0) is False


def test_outer_convergence_gas_split_gate():
    checker = OuterConvergence(
        rel_tol=1e-3,
        abs_tol_T=0.5,
        enable_gas_phase_split_gate=True,
        max_gas_phase_split_rms=0.05,
    )
    assert checker.is_aligned(
        dT_outer_max=0.2,
        T_ref_scale=1000.0,
        gas_phase_split_rms=0.03,
    ) is True
    assert checker.is_aligned(
        dT_outer_max=0.2,
        T_ref_scale=1000.0,
        gas_phase_split_rms=0.08,
    ) is False
    assert checker.is_aligned(
        dT_outer_max=2.0,
        T_ref_scale=1000.0,
        gas_phase_split_rms=0.01,
    ) is False


def test_outer_loop_history_records_gas_alignment_fields():
    cells = [_FakeCell(T=1000.0)]
    gas_calls = {"n": 0}

    def _gas_metrics() -> dict[str, float]:
        gas_calls["n"] += 1
        return {
            "gas_phase_split_rms": 0.02,
            "gas_phase_split_max_abs": 0.04,
            "bed0_primary_rel_drift": 0.01,
        }

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        cells[0].T = 1000.1
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-3,
            "norm_history": [1e-3],
            "accepted_lambda_history": [lambda_seed_outer],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    result = run_global_nr_outer_abgleich(
        cells=cells,
        outer_max=2,
        total_inner_budget=4,
        tol=1e-6,
        jacobian_mode="exact",
        jacobian_lag=1,
        tol_rms=0.01,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(
            rel_tol=1e-3,
            abs_tol_T=0.5,
            enable_gas_phase_split_gate=True,
            max_gas_phase_split_rms=0.05,
        ),
        gas_metrics_fn=_gas_metrics,
    )
    assert result.converged_outer is True
    assert gas_calls["n"] == 1
    last = result.outer_history[-1]
    assert last["gas_phase_split_rms"] == 0.02
    assert last["outer_gas_split_gate_enabled"] is True


def test_outer_loop_runs_check2_bridge_when_gas_gate_fails():
    cells = [_FakeCell(T=1000.0)]
    bridge_calls = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _gas_metrics() -> dict[str, float]:
        return {"gas_phase_split_rms": 0.08, "bed0_primary_rel_drift": 0.01}

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        cells[0].T = 1000.1
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-3,
            "norm_history": [1e-3],
            "accepted_lambda_history": [lambda_seed_outer],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    def _check2_bridge() -> dict[str, bool | float]:
        bridge_calls["n"] += 1
        return {"enabled": True, "aligned": False, "bridge_passes": 1.0}

    run_global_nr_outer_abgleich(
        cells=cells,
        outer_max=1,
        total_inner_budget=2,
        tol=1e-6,
        jacobian_mode="exact",
        jacobian_lag=1,
        tol_rms=0.01,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(
            rel_tol=1e-3,
            abs_tol_T=0.5,
            enable_gas_phase_split_gate=True,
            max_gas_phase_split_rms=0.05,
        ),
        gas_metrics_fn=_gas_metrics,
        check2_bridge_fn=_check2_bridge,
    )
    assert bridge_calls["n"] == 1


def test_outer_loop_skips_check2_bridge_on_dt_only_mismatch():
    cells = [_FakeCell(T=1000.0)]
    bridge_calls = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _gas_metrics() -> dict[str, float]:
        return {"gas_phase_split_rms": 0.02, "bed0_primary_rel_drift": 0.01}

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        cells[0].T = 1002.0
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-3,
            "norm_history": [1e-3],
            "accepted_lambda_history": [lambda_seed_outer],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    def _check2_bridge() -> dict[str, bool | float]:
        bridge_calls["n"] += 1
        return {"enabled": True, "aligned": False, "bridge_passes": 1.0}

    run_global_nr_outer_abgleich(
        cells=cells,
        outer_max=1,
        total_inner_budget=2,
        tol=1e-6,
        jacobian_mode="exact",
        jacobian_lag=1,
        tol_rms=0.01,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(
            rel_tol=1e-3,
            abs_tol_T=0.5,
            enable_gas_phase_split_gate=True,
            max_gas_phase_split_rms=0.05,
        ),
        gas_metrics_fn=_gas_metrics,
        check2_bridge_fn=_check2_bridge,
    )
    assert bridge_calls["n"] == 0


def test_outer_loop_runs_check2_bridge_when_inner_converged_but_misaligned():
    cells = [_FakeCell(T=1000.0)]
    bridge_calls = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _gas_metrics() -> dict[str, float]:
        return {"gas_phase_split_rms": 0.08, "bed0_primary_rel_drift": 0.01}

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        cells[0].T = 1002.0
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-3,
            "norm_history": [1e-3],
            "accepted_lambda_history": [lambda_seed_outer],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    def _check2_bridge() -> dict[str, bool | float]:
        bridge_calls["n"] += 1
        return {"enabled": True, "aligned": False, "bridge_passes": 1.0}

    result = run_global_nr_outer_abgleich(
        cells=cells,
        outer_max=1,
        total_inner_budget=2,
        tol=1e-6,
        jacobian_mode="exact",
        jacobian_lag=1,
        tol_rms=0.01,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(
            rel_tol=1e-3,
            abs_tol_T=0.5,
            enable_gas_phase_split_gate=True,
            max_gas_phase_split_rms=0.05,
        ),
        gas_metrics_fn=_gas_metrics,
        check2_bridge_fn=_check2_bridge,
    )
    assert bridge_calls["n"] == 1
    assert result.outer_history[-1]["check2_bridge_passes"] == 1.0


def test_outer_loop_history_records_alignment_fields():
    cells = [_FakeCell(T=1000.0)]
    step = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        step["n"] += 1
        if step["n"] == 1:
            cells[0].T = 1002.0
            return {
                "converged": False,
                "n_iter": 1,
                "n_newton_iters_attempted": min(inner_iter_cap, 1),
                "rms_scaled_final": 1e-2,
                "norm_history": [1e-2],
                "accepted_lambda_history": [lambda_seed_outer],
                "line_search_trial_counts": [1],
                "clip_history": [],
            }
        cells[0].T = 1000.1
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": min(inner_iter_cap, 1),
            "rms_scaled_final": 1e-4,
            "norm_history": [1e-4],
            "accepted_lambda_history": [lambda_seed_outer],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=3,
        total_inner_budget=4,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
    )

    assert out.outer_iters >= 2
    assert any(item["outer_aligned"] is False for item in out.outer_history)
    assert any(item["outer_aligned"] is True for item in out.outer_history)
    assert all("outer_dT_tol_effective" in item for item in out.outer_history)
    assert all("dT_refresh" in item for item in out.outer_history)
    assert all("lambda_seed_used" in item for item in out.outer_history)


def test_outer_loop_resets_lambda_seed_after_large_refresh_jump():
    cells = [_FakeCell(T=1000.0)]
    seen_lambda: list[float] = []

    def _refresh(_force: bool) -> None:
        cells[0].T = 1055.0

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        seen_lambda.append(float(lambda_seed_outer))
        if len(seen_lambda) == 1:
            # Keep Check2 raw alignment false so outer 2 still executes refresh.
            cells[0].T = 1200.0
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": min(inner_iter_cap, 1),
            "rms_scaled_final": 1e-2,
            "norm_history": [1e-2],
            "accepted_lambda_history": [1.0 / 1024.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=1,
        total_inner_budget=1,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
    )

    assert seen_lambda == [0.5]
    assert out.outer_history[0]["lambda_seed_reset_by_refresh"] is True


def test_outer_loop_resets_lambda_seed_after_signature_change_with_tiny_seed():
    cells = [_FakeCell(T=1000.0)]
    seen_lambda: list[float] = []
    sig_state = {"v": 0}

    def _refresh(_force: bool) -> None:
        sig_state["v"] += 1

    def _signature() -> str:
        return f"sig:v{sig_state['v']}"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        seen_lambda.append(float(lambda_seed_outer))
        if len(seen_lambda) == 1:
            # Keep raw Check2 misaligned in outer-1, so outer-2 still executes refresh.
            cells[0].T = 1200.0
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": min(inner_iter_cap, 1),
            "rms_scaled_final": 1e-2,
            "norm_history": [1e-2],
            "accepted_lambda_history": [1.0 / 1024.0],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=2,
        total_inner_budget=2,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
    )

    # First outer starts from default 0.5; second outer would inherit tiny seed,
    # but refresh signature change should reset it back to 0.5.
    assert seen_lambda == [0.5, 0.5]
    assert out.outer_history[1]["signature_changed_by_refresh"] is True
    assert out.outer_history[1]["lambda_seed_reset_by_refresh"] is True


def test_strict_check1_gate_skips_refresh_until_inner_converges():
    cells = [_FakeCell(T=1000.0)]
    refresh_calls = {"n": 0}
    seen_caps: list[int] = []
    step = {"n": 0}

    def _refresh(_force: bool) -> None:
        refresh_calls["n"] += 1

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        step["n"] += 1
        seen_caps.append(int(inner_iter_cap))
        if step["n"] == 1:
            # first inner call: not converged
            cells[0].T = 1200.0
            return {
                "converged": False,
                "n_iter": 1,
                "n_newton_iters_attempted": min(inner_iter_cap, 1),
                "rms_scaled_final": 1e-2,
                "norm_history": [1e-2],
                "accepted_lambda_history": [0.25],
                "line_search_trial_counts": [1],
                "clip_history": [],
            }
        # second inner call: converged and aligned
        cells[0].T = 1000.0
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": min(inner_iter_cap, 1),
            "rms_scaled_final": 1e-4,
            "norm_history": [1e-4],
            "accepted_lambda_history": [0.25],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=3,
        total_inner_budget=6,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
        strict_check1_before_refresh=True,
    )

    assert refresh_calls["n"] >= 2
    assert seen_caps[0] == 2
    assert seen_caps[1] == 5
    assert out.outer_history[1]["refresh_skipped_for_check1_gate"] is True


def test_outer_loop_does_not_stop_early_when_inner_rejects_all_line_search_steps():
    cells = [_FakeCell(T=1000.0)]
    outer_calls = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        outer_calls["n"] += 1
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 2,
            "rms_scaled_final": 0.02,
            "norm_history": [0.02],
            "accepted_lambda_history": [None],
            "line_search_trial_counts": [12],
            "clip_history": [{"line_search_failed": True}],
            "counts": {"line_search_failures": 1},
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=5,
        total_inner_budget=10,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
    )

    assert outer_calls["n"] == 5
    assert out.outer_iters == 5


def test_outer_loop_runs_preinner_hook_after_refresh_not_when_refresh_skipped():
    cells = [_FakeCell(T=1000.0)]
    preinner_calls: list[bool] = []

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _pre_inner() -> dict[str, bool | float]:
        preinner_calls.append(True)
        return {
            "enabled": True,
            "split_rms_before": 0.02,
            "split_rms_after": 0.01,
            "merit_before": 0.02,
            "merit_after": 0.01,
            "bed0_accepted": True,
            "bed1_accepted": False,
        }

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 0.02,
            "norm_history": [0.02],
            "accepted_lambda_history": [None],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=3,
        total_inner_budget=6,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        pre_inner_fn=_pre_inner,
        pre_inner_on_continuation=False,
        strict_check1_before_refresh=True,
    )

    assert len(preinner_calls) == 2
    assert out.outer_history[0]["preinner_phase_split_enabled"] is True
    assert out.outer_history[0]["preinner_bed0_accepted"] is True
    assert out.outer_history[1]["refresh_skipped_for_check1_gate"] is True
    assert out.outer_history[1]["preinner_phase_split_enabled"] is False
    assert out.outer_history[2]["refresh_skipped_for_check1_gate"] is False
    assert out.outer_history[2]["preinner_phase_split_enabled"] is True


def test_outer_loop_runs_preinner_on_check1_continuation_when_enabled():
    cells = [_FakeCell(T=1000.0)]
    preinner_calls = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _pre_inner() -> dict[str, bool | float]:
        preinner_calls["n"] += 1
        return {"enabled": True, "bed0_accepted": True}

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 0.02,
            "norm_history": [0.02],
            "accepted_lambda_history": [None],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=3,
        total_inner_budget=6,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        pre_inner_fn=_pre_inner,
        pre_inner_on_continuation=True,
        strict_check1_before_refresh=True,
    )

    assert preinner_calls["n"] == 3
    assert out.outer_history[1]["preinner_phase_split_enabled"] is True


def test_strict_check1_reallows_refresh_when_continuation_stalls():
    """When Check1 continuation line-search stalls, next outer must Vorab-refresh again."""
    cells = [_FakeCell(T=1000.0)]
    refresh_calls = {"n": 0}
    step = {"n": 0}

    def _refresh(_force: bool) -> None:
        refresh_calls["n"] += 1

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        step["n"] += 1
        if step["n"] == 1:
            cells[0].T = 1200.0
            return {
                "converged": False,
                "n_iter": 2,
                "n_newton_iters_attempted": 2,
                "rms_scaled_final": 0.5,
                "norm_history": [0.8, 0.5],
                "accepted_lambda_history": [0.5, 0.25],
                "line_search_trial_counts": [1, 1],
                "clip_history": [],
            }
        if step["n"] == 2:
            return {
                "converged": False,
                "n_iter": 1,
                "n_newton_iters_attempted": 2,
                "rms_scaled_final": 0.5,
                "norm_history": [0.5],
                "accepted_lambda_history": [None],
                "line_search_trial_counts": [4],
                "clip_history": [{"line_search_failed": True}],
            }
        cells[0].T = 1000.0
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-4,
            "norm_history": [1e-4],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=4,
        total_inner_budget=8,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
        strict_check1_before_refresh=True,
    )

    assert refresh_calls["n"] >= 3
    assert out.converged_outer is True
    assert out.outer_history[1]["refresh_skipped_for_check1_gate"] is True
    assert out.outer_history[2]["refresh_skipped_for_check1_gate"] is False


def test_strict_check1_reallows_refresh_on_energy_regression():
    """Energy RMS spike under Check1 continuation must re-open Vorab refresh."""
    cells = [_FakeCell(T=1000.0)]
    refresh_calls = {"n": 0}
    contexts: list[bool] = []
    step = {"n": 0}

    def _refresh(_force: bool) -> None:
        refresh_calls["n"] += 1

    def _signature() -> str:
        return f"sig:{refresh_calls['n']}"

    def _context(*, refresh_skipped: bool = False) -> None:
        contexts.append(bool(refresh_skipped))

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        step["n"] += 1
        if step["n"] == 1:
            cells[0].T = 1100.0
            return {
                "converged": False,
                "n_iter": 2,
                "n_newton_iters_attempted": 2,
                "rms_scaled_final": 0.40,
                "rms_scaled_energy_final": 0.40,
                "rms_scaled_gas_phase_split_final": 0.05,
                "norm_history": [0.5, 0.40],
                "accepted_lambda_history": [0.5],
                "line_search_trial_counts": [1],
                "clip_history": [],
            }
        if step["n"] == 2:
            cells[0].T = 1300.0
            return {
                "converged": False,
                "n_iter": 2,
                "n_newton_iters_attempted": 2,
                "rms_scaled_final": 0.90,
                "rms_scaled_energy_final": 0.90,
                "rms_scaled_gas_phase_split_final": 0.06,
                "norm_history": [0.40, 0.90],
                "accepted_lambda_history": [0.25],
                "line_search_trial_counts": [1],
                "clip_history": [],
            }
        cells[0].T = 1000.0
        return {
            "converged": True,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 1e-4,
            "rms_scaled_energy_final": 1e-4,
            "rms_scaled_gas_phase_split_final": 1e-4,
            "norm_history": [1e-4],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=4,
        total_inner_budget=8,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        outer_convergence=OuterConvergence(rel_tol=1e-3, abs_tol_T=0.5),
        strict_check1_before_refresh=True,
        check1_energy_regress_refresh_ratio=1.25,
        outer_context_fn=_context,
        merit_early_stop=OuterMeritEarlyStop(
            enabled=False,
            use_gas_phase_split_merit=True,
            use_energy_merit=True,
        ),
    )

    assert out.outer_history[1]["refresh_skipped_for_check1_gate"] is True
    # After O2 energy 0.40 → 0.90 (>1.25×), O3 must refresh again.
    assert out.outer_history[2]["refresh_skipped_for_check1_gate"] is False
    assert refresh_calls["n"] >= 3
    assert True in contexts and False in contexts


def test_resolve_outer_composite_merit_uses_split_and_energy():
    early = OuterMeritEarlyStop(
        enabled=True,
        use_gas_phase_split_merit=True,
        use_energy_merit=True,
    )
    metrics = resolve_outer_composite_merit(
        {
            "rms_scaled_final": 0.03,
            "rms_scaled_gas_phase_split_final": 0.05,
            "rms_scaled_energy_final": 0.04,
        },
        {},
        early_stop=early,
    )
    assert metrics["merit_composite"] == 0.05


def test_outer_loop_merit_early_stop_on_composite_regression():
    cells = [_FakeCell(T=1000.0)]
    merits = [0.04, 0.038, 0.07, 0.08]

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        merit = merits[min(len(merits) - 1, _solve_inner.calls)]
        _solve_inner.calls += 1
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": merit * 0.5,
            "rms_scaled_gas_phase_split_final": merit * 0.8,
            "rms_scaled_energy_final": merit,
            "norm_history": [merit],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
            "counts": {},
        }

    _solve_inner.calls = 0

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=10,
        total_inner_budget=20,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        merit_early_stop=OuterMeritEarlyStop(
            enabled=True,
            min_outer_iters=2,
            regression_rel_tol=0.25,
            regression_patience=1,
        ),
    )

    assert out.merit_early_stop is True
    assert out.merit_early_stop_reason == "composite_merit_regression"
    assert out.outer_iters == 3


def test_outer_loop_merit_plateau_stops_on_composite_stall():
    cells = [_FakeCell(T=1000.0)]
    sequence = [
        (0.03, 0.03, 0.03),
        (0.029, 0.029, 0.029),
    ]

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        idx = min(_solve_inner.calls, len(sequence) - 1)
        _solve_inner.calls += 1
        rms, split, energy = sequence[idx]
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": rms,
            "rms_scaled_gas_phase_split_final": split,
            "rms_scaled_energy_final": energy,
            "norm_history": [max(rms, split, energy)],
            "accepted_lambda_history": [0.5],
            "line_search_trial_counts": [1],
            "clip_history": [],
            "counts": {},
        }

    _solve_inner.calls = 0

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=10,
        total_inner_budget=20,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        merit_early_stop=OuterMeritEarlyStop(
            enabled=True,
            min_outer_iters=1,
            plateau_rel_tol=0.002,
        ),
    )

    assert out.merit_early_stop is True
    assert out.merit_early_stop_reason == "composite_merit_plateau"
    assert out.outer_iters >= 2


def test_outer_loop_hot_restart_resets_check1_gate_after_y_co_band():
    """CO 带宽后下一 outer 强制允许 refresh（等价 skip_precalc 门禁复位）。"""
    cells = [_FakeCell(T=1000.0)]
    # 按 post-outer 调用次序：outer1 低、outer2 进带宽 → outer3 热重启
    y_seq = [0.10, 0.15, 0.15, 0.15]
    y_calls = {"n": 0}
    clear_calls = {"n": 0}

    def _refresh(_force: bool) -> None:
        return

    def _signature() -> str:
        return "sig:v1"

    def _y_co() -> float:
        idx = min(y_calls["n"], len(y_seq) - 1)
        y_calls["n"] += 1
        return float(y_seq[idx])

    def _clear() -> None:
        clear_calls["n"] += 1

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        return {
            "converged": False,
            "n_iter": 1,
            "n_newton_iters_attempted": 1,
            "rms_scaled_final": 0.05,
            "norm_history": [0.05],
            "accepted_lambda_history": [None],
            "line_search_trial_counts": [1],
            "clip_history": [],
        }

    out = run_global_nr_outer_abgleich(
        cells=cells,  # type: ignore[arg-type]
        outer_max=4,
        total_inner_budget=8,
        tol=1e-3,
        jacobian_mode="block_tridiag_structured",
        jacobian_lag=1,
        tol_rms=1e-2,
        refresh_fn=_refresh,
        snapshot_signature_fn=_signature,
        solve_inner_fn=_solve_inner,
        strict_check1_before_refresh=True,
        hot_restart_after_y_co=True,
        hot_restart_y_co_threshold=0.14,
        hot_restart_bonus_inner_budget=4,
        hot_restart_bonus_outer_slots=2,
        y_co_top_fn=_y_co,
        hot_restart_clear_runtime_fn=_clear,
    )

    armed = [h for h in out.outer_history if h.get("outer_hot_restart_armed")]
    restarted = [h for h in out.outer_history if h.get("outer_hot_restarted")]
    assert len(armed) == 1
    assert armed[0]["outer_iter"] == 2
    assert len(restarted) == 1
    assert restarted[0]["outer_iter"] == 3
    assert restarted[0]["refresh_skipped_for_check1_gate"] is False
    assert clear_calls["n"] == 1
    # outer1 refresh, outer2 Check1-skip, outer3 hot-restart refresh
    assert out.outer_history[0]["refresh_skipped_for_check1_gate"] is False
    assert out.outer_history[1]["refresh_skipped_for_check1_gate"] is True
    assert out.outer_history[2]["refresh_skipped_for_check1_gate"] is False


def test_outer_hot_restart_disabled_by_default():
    from src.core.reactor import ReactorConfig

    cfg = ReactorConfig()
    assert cfg.nr_outer_hot_restart_after_y_co_thesis is False
    assert cfg.nr_snap_trace_o2_holdup_thesis is False
    assert cfg.nr_outer_freeze_upper_co_bridge_after_y_co_thesis is False