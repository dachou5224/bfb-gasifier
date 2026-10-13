"""审查 2026-09-15 第 7 章：导数边界、线性步、推进模式、耦合与端点认证。"""

from __future__ import annotations

import numpy as np
import pytest

from src.solvers.wirsum_1997_named_f import default_x0_named
from src.solvers.wirsum_1997_named_f_continuation import (
    certify_named_full_rate_endpoint,
    named_in_domain,
    named_linear_step,
    stage_passed,
)
from src.solvers.wirsum_1997_named_f_e4 import freeze_e4_q, residual_named_e4
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q
from src.solvers.wirsum_1997_named_f_r10 import (
    freeze_e6_r10_q,
    residual_e6_plus_scaled_r10,
)
from src.solvers.wirsum_1997_named_f_r12 import _r12_r_and_dC, freeze_r12_only_q
from src.solvers.wirsum_1997_named_f_r12_r5 import (
    _r5_bubble_r_and_dC,
    _r5_suspension_r_and_dC,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ


def test_r12_zero_oxygen_keeps_finite_right_derivative() -> None:
    C = np.zeros(N_MAJ, dtype=np.float64)
    C[IDX["H2"]] = 1.0
    r, dr = _r12_r_and_dC(1100.0, C)
    assert r == pytest.approx(0.0)
    analytic = float(dr[IDX["O2"]])
    assert analytic == pytest.approx(2.65633e12, rel=1.0e-4)
    diffs = []
    for h in (1.0e-8, 1.0e-10, 1.0e-12):
        Cp = C.copy()
        Cp[IDX["O2"]] = h
        rp, _ = _r12_r_and_dC(1100.0, Cp)
        diffs.append(float((rp - r) / h))
    assert diffs[-1] == pytest.approx(analytic, rel=1.0e-3)


def test_r12_zero_h2_has_zero_h2_derivative() -> None:
    C = np.zeros(N_MAJ, dtype=np.float64)
    C[IDX["O2"]] = 1.0
    r, dr = _r12_r_and_dC(1100.0, C)
    assert r == pytest.approx(0.0)
    assert float(dr[IDX["H2"]]) == pytest.approx(0.0)
    assert float(dr[IDX["O2"]]) == pytest.approx(0.0)


def test_r5_fractional_oxygen_has_no_finite_right_derivative() -> None:
    C = np.zeros(N_MAJ, dtype=np.float64)
    C[IDX["CO"]] = 1.0
    C[IDX["H2O"]] = 1.0
    r_b, dr_b = _r5_bubble_r_and_dC(1100.0, C)
    r_s, dr_s = _r5_suspension_r_and_dC(1100.0, C)
    assert r_b == pytest.approx(0.0)
    assert r_s == pytest.approx(0.0)
    assert float(dr_b[IDX["O2"]]) == pytest.approx(0.0)
    assert float(dr_s[IDX["O2"]]) == pytest.approx(0.0)
    C_h2 = np.zeros(N_MAJ, dtype=np.float64)
    C_h2[IDX["H2"]] = 1.0
    _, dr12 = _r12_r_and_dC(1100.0, C_h2)
    assert float(dr12[IDX["O2"]]) != 0.0


def test_r5_zero_co_keeps_finite_right_derivative() -> None:
    C = np.zeros(N_MAJ, dtype=np.float64)
    C[IDX["O2"]] = 1.0
    C[IDX["H2O"]] = 1.0
    r, dr = _r5_bubble_r_and_dC(1100.0, C)
    assert r == pytest.approx(0.0)
    assert float(dr[IDX["CO"]]) > 0.0


def test_r10_residual_does_not_clip_negative_state() -> None:
    nq6 = freeze_e6_q(Q_W=0.0)
    nq_all = freeze_e6_r10_q(Q_W=0.0)
    x = default_x0_named(nq6)
    x[38] = -1.0e-5
    assert named_in_domain(x) is False
    f_neg = residual_e6_plus_scaled_r10(x, nq6, nq_all, r10_scale=1.0)
    x_clip = x.copy()
    x_clip[:-1] = np.maximum(x_clip[:-1], 0.0)
    f_clip = residual_e6_plus_scaled_r10(x_clip, nq6, nq_all, r10_scale=1.0)
    assert not np.allclose(f_neg, f_clip)


def test_named_linear_step_reports_raw_vs_modified_jacobian() -> None:
    n = 40
    x = np.ones(n, dtype=np.float64)
    x[-1] = 800.0
    F = np.ones(n, dtype=np.float64)
    J_full = np.eye(n, dtype=np.float64)
    _, _, _, diag_full = named_linear_step(J_full, F, x, extra_diag_species=("H2",))
    assert diag_full["used_trace_diag"] is False
    assert diag_full["linear_residual_raw_l2"] < 1.0e-10
    assert diag_full["directional_phi_dot"] == pytest.approx(
        diag_full["minus_f_l2_sq"], rel=1.0e-10, abs=1.0e-10,
    )

    J_def = np.eye(n, dtype=np.float64)
    i_h2 = IDX["H2"]
    J_def[i_h2, i_h2] = 0.0
    _, _, _, diag_mod = named_linear_step(J_def, F, x, extra_diag_species=("H2",))
    assert diag_mod["used_trace_diag"] is True
    assert diag_mod["equilibrated_rank_original"] < n
    assert diag_mod["equilibrated_rank_used"] == n
    assert diag_mod["linear_residual_modified_l2"] < 1.0e-8
    assert diag_mod["linear_residual_raw_l2"] > diag_mod["linear_residual_modified_l2"]


def test_stage_passed_distinguishes_advance_modes() -> None:
    x = np.ones(8, dtype=np.float64)
    x[-1] = 800.0
    stop_ok = {"residual_l2": 1.0e-8, "converged_newton": False}
    stop_big = {"residual_l2": 1.0, "converged_newton": False}
    assert stage_passed(stop_ok, x, advance_mode="root_correcting", accepted_alphas=[])
    assert not stage_passed(stop_big, x, advance_mode="root_correcting", accepted_alphas=[0.5])
    assert stage_passed(
        stop_big, x, advance_mode="heuristic_search", accepted_alphas=[0.5],
    )
    assert not stage_passed(
        stop_big, x, advance_mode="heuristic_search", accepted_alphas=[],
    )
    x_neg = x.copy()
    x_neg[0] = -1.0
    assert not stage_passed(
        stop_ok, x_neg, advance_mode="root_correcting", accepted_alphas=[1.0],
    )


def test_e4_composite_p_differs_from_frozen_q() -> None:
    nq4 = freeze_e4_q()
    x = default_x0_named(nq4)
    x[-1] = 1312.0
    f_p = residual_named_e4(x, nq4, coupling="composite_P")
    f_q = residual_named_e4(x, nq4, coupling="frozen_q")
    assert float(np.linalg.norm(f_p - f_q)) > 1.0


def test_certify_rejects_out_of_domain_without_clip() -> None:
    nq = freeze_r12_only_q()
    x = default_x0_named(nq)
    x[0] = -1.0e-5
    cert = certify_named_full_rate_endpoint(x, nq, coupling="frozen_q")
    assert cert["in_domain"] is False
    assert cert["clip_applied"] is False
    assert cert["endpoint_certified"] is False
    assert cert["rate_scale"] == 1.0
