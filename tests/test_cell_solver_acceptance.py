from __future__ import annotations

from src.solvers.cell_solver import _prefer_candidate_metrics


def test_prefer_candidate_metrics_rejects_scaled_only_improvement():
    assert not _prefer_candidate_metrics(
        incumbent_rms=0.20,
        incumbent_res=10.0,
        candidate_rms=0.15,
        candidate_res=20.0,
        residual_tol_scaled=0.10,
    )


def test_prefer_candidate_metrics_accepts_joint_improvement():
    assert _prefer_candidate_metrics(
        incumbent_rms=0.20,
        incumbent_res=10.0,
        candidate_rms=0.15,
        candidate_res=9.0,
        residual_tol_scaled=0.10,
    )
