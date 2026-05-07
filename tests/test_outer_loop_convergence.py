from __future__ import annotations

from dataclasses import dataclass

from src.solvers.convergence import InnerConvergence, OuterConvergence
from src.solvers.outer_loop import run_global_nr_outer_abgleich


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
