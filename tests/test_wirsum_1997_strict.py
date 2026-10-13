"""P2：Wirsum 1997 strict Newton（Eq. 2.21–2.23）与四组 A/B。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.solid_balance_reduced import (
    N_CELLS,
    canonical_open_chain_K,
    holdup_jacobian,
    holdup_residual,
)
from src.solvers.wirsum_1997_strict import (
    NR_SOLVER_MODE_WIRSUM_1997_STRICT,
    apply_wirsum_1997_strict_config,
    build_reduced_holdup_ab_payload,
    eval_merit,
    resolve_nr_solver_mode,
    run_wirsum_strict_ab,
    wirsum_1997_strict_newton,
)


def _feed_cell0(rate: float = 1.0) -> np.ndarray:
    m = np.zeros(N_CELLS, dtype=np.float64)
    m[0] = rate  # [kg/s]
    return m


def test_resolve_nr_solver_mode_rejects_unknown() -> None:
    assert resolve_nr_solver_mode("Engineering") == "engineering"
    with pytest.raises(ValueError, match="nr_solver_mode"):
        resolve_nr_solver_mode("clip")


def test_linear_holdup_converges_in_one_full_newton_step() -> None:
    K_up, K_down = canonical_open_chain_K()
    m_feed = _feed_cell0(1.0)
    r_dot = np.zeros(N_CELLS, dtype=np.float64)

    def residual(M: np.ndarray) -> np.ndarray:
        return holdup_residual(
            M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
        )

    def jacobian(M: np.ndarray) -> np.ndarray:
        return holdup_jacobian(
            M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
        )

    x0 = np.array([0.1, 0.1, 0.1], dtype=np.float64)  # [kg]
    result = wirsum_1997_strict_newton(residual, jacobian, x0)
    assert result.converged
    assert result.iterations == 1
    assert result.accepted_alphas == [1.0]
    assert result.n_jac_builds == 1
    dx = np.linalg.solve(jacobian(x0), -residual(x0))
    np.testing.assert_allclose(result.x, x0 + dx, atol=1e-12)
    np.testing.assert_allclose(result.x, np.array([25.0 / 6.0, 10.0 / 3.0, 2.0]), atol=1e-12)
    np.testing.assert_allclose(result.F, 0.0, atol=1e-12)


def test_trial_in_domain_skips_negative_residual_eval() -> None:
    evals: list[np.ndarray] = []

    def residual(x: np.ndarray) -> np.ndarray:
        evals.append(np.asarray(x, dtype=np.float64).copy())
        assert np.all(x >= 0.0)
        return x + 1.0  # [与 x 同形] 根在负域

    def jacobian(x: np.ndarray) -> np.ndarray:
        return np.eye(x.size)

    def in_domain(x: np.ndarray) -> bool:
        return bool(np.all(np.asarray(x) >= 0.0))

    x0 = np.array([0.25], dtype=np.float64)
    result = wirsum_1997_strict_newton(
        residual, jacobian, x0, trial_in_domain=in_domain, max_iter=8
    )
    assert all(float(v[0]) >= 0.0 for v in evals)
    assert not result.converged
    assert result.failure_reason in {"line_search", "max_iter"}


def test_K_zero_holdup_reports_singular_jacobian() -> None:
    z = np.zeros(N_CELLS, dtype=np.float64)
    m_feed = _feed_cell0(1.0)

    def residual(M: np.ndarray) -> np.ndarray:
        return holdup_residual(
            M, m_feed=m_feed, K_up=z, K_down=z, r_dot=z, k0_anchor="none"
        )

    def jacobian(M: np.ndarray) -> np.ndarray:
        return holdup_jacobian(
            M, m_feed=m_feed, K_up=z, K_down=z, r_dot=z, k0_anchor="none"
        )

    result = wirsum_1997_strict_newton(residual, jacobian, np.ones(N_CELLS))
    assert not result.converged
    assert result.failure_reason == "singular_jacobian"
    assert result.rank_last == 0
    assert result.n_zero_rows_last == N_CELLS


def test_arctan_needs_dyadic_backtracking() -> None:
    def residual(x: np.ndarray) -> np.ndarray:
        return np.arctan(x)

    def jacobian(x: np.ndarray) -> np.ndarray:
        return np.array([[1.0 / (1.0 + float(x[0]) ** 2)]], dtype=np.float64)

    x0 = np.array([10.0], dtype=np.float64)
    result = wirsum_1997_strict_newton(residual, jacobian, x0, max_iter=40)
    assert result.converged
    assert min(result.accepted_alphas) < 1.0
    np.testing.assert_allclose(result.x, 0.0, atol=1e-8)
    assert result.residual_l2_history[-1] < result.residual_l2_history[0]


def test_raw_l2_merit_is_euclidean_norm() -> None:
    F = np.array([3.0, 4.0], dtype=np.float64)
    assert eval_merit(F, kind="raw_l2", scale=None) == pytest.approx(5.0)
    scale = np.array([1.0, 2.0], dtype=np.float64)
    # RMS of [3, 2]
    assert eval_merit(F, kind="scaled_rms", scale=scale) == pytest.approx(
        np.sqrt((9.0 + 4.0) / 2.0)
    )


def _mm_holdup_fns():
    K_up, K_down = canonical_open_chain_K()
    m_feed = _feed_cell0(1.0)
    r_dot = np.zeros(N_CELLS, dtype=np.float64)
    k_sat = np.array([0.50, 0.20, 0.10], dtype=np.float64)  # [kg/s]
    M_ref = np.ones(N_CELLS, dtype=np.float64)  # [kg]

    def residual(M: np.ndarray) -> np.ndarray:
        F = holdup_residual(
            M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
        )
        return F - k_sat * M / (1.0 + M / M_ref)

    def jacobian(M: np.ndarray) -> np.ndarray:
        J = holdup_jacobian(
            M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
        ).copy()
        J[np.diag_indices(N_CELLS)] -= k_sat * M_ref / (M_ref + M) ** 2
        return J

    return residual, jacobian


def test_nonlinear_holdup_ab_fresh_beats_or_matches_frozen() -> None:
    residual, jacobian = _mm_holdup_fns()
    x0 = np.array([0.05, 0.05, 0.05], dtype=np.float64)  # [kg]
    scale = np.ones(N_CELLS, dtype=np.float64)  # [kg/s]
    table = run_wirsum_strict_ab(
        residual, jacobian, x0, scale=scale, max_iter=40, tol_merit=1.0e-10
    )
    a = table["A_raw_fresh"]
    c = table["C_raw_frozen"]
    assert a["converged"]
    assert a["n_jac_builds"] == a["iterations"]
    assert c["n_jac_builds"] == 1
    if c["converged"]:
        assert int(c["iterations"]) >= int(a["iterations"])
    np.testing.assert_allclose(a["F"], 0.0, atol=1e-8)


def test_row_weighted_holdup_raw_and_scaled_merit_histories_differ() -> None:
    residual, jacobian = _mm_holdup_fns()
    weights = np.array([1.0, 1.0e3, 1.0e-3], dtype=np.float64)

    def residual_w(M: np.ndarray) -> np.ndarray:
        return weights * residual(M)

    def jacobian_w(M: np.ndarray) -> np.ndarray:
        return weights[:, None] * jacobian(M)

    x0 = np.array([0.05, 0.05, 0.05], dtype=np.float64)  # [kg]
    table = run_wirsum_strict_ab(
        residual_w, jacobian_w, x0, scale=weights, max_iter=40, tol_merit=1.0e-10
    )
    assert table["A_raw_fresh"]["converged"]
    assert table["B_scaled_fresh"]["converged"]
    assert table["C_raw_frozen"]["n_jac_builds"] == 1
    assert table["D_scaled_frozen"]["n_jac_builds"] == 1
    assert table["A_raw_fresh"]["merit_history"] != table["B_scaled_fresh"]["merit_history"]


def test_reactor_config_strict_mode_disables_engineering_bridges() -> None:
    from src.core.reactor import ReactorConfig

    cfg = ReactorConfig(
        nr_solver_mode=NR_SOLVER_MODE_WIRSUM_1997_STRICT,
        holdup_k0_anchor="engineering",
        nr_line_search_energy_merit_thesis=True,
        nr_line_search_fast_oxidation_projection_thesis=True,
        nr_line_search_syngas_collapse_guard_thesis=True,
    )
    assert cfg.nr_solver_mode == NR_SOLVER_MODE_WIRSUM_1997_STRICT
    assert cfg.holdup_k0_anchor == "none"
    assert cfg.nr_line_search_energy_merit_thesis is False
    assert cfg.nr_line_search_fast_oxidation_projection_thesis is False
    assert cfg.nr_line_search_syngas_collapse_guard_thesis is False
    assert cfg.nr_step_model_thesis == "newton"
    apply_wirsum_1997_strict_config(cfg)
    assert cfg.nr_nonmonotone_enabled_thesis is False


def test_reduced_holdup_ab_verdicts() -> None:
    payload = build_reduced_holdup_ab_payload()
    linear = payload["linear_K_positive"]
    assert all(case["converged"] and case["iterations"] == 1 for case in linear.values())
    mm = payload["nonlinear_mm"]
    assert mm["A_raw_fresh"]["converged"] and mm["A_raw_fresh"]["iterations"] == 5
    assert mm["B_scaled_fresh"]["iterations"] == mm["A_raw_fresh"]["iterations"]
    assert not mm["C_raw_frozen"]["converged"]
    assert mm["C_raw_frozen"]["failure_reason"] == "max_iter"
    weighted = payload["nonlinear_mm_row_weighted"]
    assert weighted["A_raw_fresh"]["converged"]
    assert min(weighted["A_raw_fresh"]["accepted_alphas"]) < 1.0
    assert weighted["B_scaled_fresh"]["accepted_alphas"] == [1.0] * weighted["B_scaled_fresh"]["iterations"]
    assert weighted["C_raw_frozen"]["failure_reason"] == "line_search"
    k0 = payload["K_zero_no_constitutive_outflow"]
    assert k0["failure_reason"] == "singular_jacobian"
    assert k0["rank_last"] == 0
