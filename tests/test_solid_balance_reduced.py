"""P1.1：三格固相流平衡 vs 库存平衡的约化对比。

不启动 Reactor。对照 Model A（m_out [kg/s]）与 Model B（M [kg], m=K M）。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_balances import calc_solid_balance_residual
from src.core.solid_balance_reduced import (
    N_CELLS,
    alpha_from_K,
    canonical_open_chain_K,
    holdup_jacobian,
    holdup_residual,
    jacobian_diagnostics,
    stream_jacobian,
    stream_residual,
)


def _zeros() -> np.ndarray:
    return np.zeros(N_CELLS, dtype=np.float64)


def _feed_cell0(rate: float = 1.0) -> np.ndarray:
    m = _zeros()
    m[0] = rate  # [kg/s]
    return m


def _fd_jacobian(fn, x: np.ndarray, h: float = 1.0e-7) -> np.ndarray:
    n = x.size
    J = np.zeros((n, n), dtype=np.float64)
    f0 = fn(x)
    for j in range(n):
        xp = x.copy()
        xp[j] += h
        J[:, j] = (fn(xp) - f0) / h
    return J


def test_model_a_stream_balance_is_full_rank_and_unique() -> None:
    K_up, K_down = canonical_open_chain_K()
    alpha = alpha_from_K(K_up, K_down)
    J = stream_jacobian(alpha_up=alpha)
    diag = jacobian_diagnostics(J)
    assert diag.rank == N_CELLS
    assert diag.nullity == 0
    assert diag.n_zero_rows == 0
    assert diag.n_zero_cols == 0
    assert np.isfinite(diag.cond)

    m_feed = _feed_cell0(1.0)
    r_dot = _zeros()
    m_out = np.linalg.solve(J, -(m_feed + r_dot))
    F = stream_residual(m_out, m_feed=m_feed, alpha_up=alpha, r_dot=r_dot)
    np.testing.assert_allclose(F, 0.0, atol=1e-12)
    assert m_out[0] > 0.0
    # 床顶出流带走全部给料（无反应）
    m_up, m_down = alpha * m_out, (1.0 - alpha) * m_out
    assert m_down[0] == pytest.approx(0.0)
    np.testing.assert_allclose(m_up[-1], 1.0, atol=1e-12)


def test_model_b_holdup_with_K_positive_matches_model_a_streams() -> None:
    K_up, K_down = canonical_open_chain_K()
    alpha = alpha_from_K(K_up, K_down)
    m_feed = _feed_cell0(1.0)
    r_dot = _zeros()
    M = np.array([2.0, 1.5, 0.8], dtype=np.float64)  # [kg]
    F_b = holdup_residual(
        M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    m_out = (K_up + K_down) * M
    F_a = stream_residual(m_out, m_feed=m_feed, alpha_up=alpha, r_dot=r_dot)
    np.testing.assert_allclose(F_a, F_b, atol=1e-14)

    J_b = holdup_jacobian(
        M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    diag = jacobian_diagnostics(J_b)
    assert diag.rank == N_CELLS
    assert diag.n_zero_rows == 0
    M_star = np.linalg.solve(J_b, -m_feed)
    assert np.all(M_star > 0.0)
    F_star = holdup_residual(
        M_star, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    np.testing.assert_allclose(F_star, 0.0, atol=1e-12)


def test_model_b_K_zero_no_source_has_zero_jacobian_and_any_M_is_steady() -> None:
    K_up, K_down = _zeros(), _zeros()
    m_feed, r_dot = _zeros(), _zeros()
    M_a = np.array([0.0, 1.0, 4.0], dtype=np.float64)
    M_b = np.array([3.0, 0.2, 0.0], dtype=np.float64)
    F_a = holdup_residual(
        M_a, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    F_b = holdup_residual(
        M_b, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    np.testing.assert_allclose(F_a, 0.0)
    np.testing.assert_allclose(F_b, 0.0)
    J = holdup_jacobian(
        M_a, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    diag = jacobian_diagnostics(J)
    assert diag.rank == 0
    assert diag.nullity == N_CELLS
    assert diag.n_zero_rows == N_CELLS
    assert diag.n_zero_cols == N_CELLS


def test_model_b_K_zero_with_feed_has_no_steady_state_and_residual_ignores_M() -> None:
    K_up, K_down = _zeros(), _zeros()
    m_feed = _feed_cell0(1.0)
    r_dot = _zeros()
    M_a = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    M_b = np.array([10.0, 5.0, 1.0], dtype=np.float64)
    F_a = holdup_residual(
        M_a, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    F_b = holdup_residual(
        M_b, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    np.testing.assert_allclose(F_a, F_b)
    np.testing.assert_allclose(F_a, m_feed)
    J = holdup_jacobian(
        M_a, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    diag = jacobian_diagnostics(J)
    assert diag.rank == 0
    assert diag.n_zero_rows == N_CELLS
    # 无本构出流时给料无法离开，不存在有限库存稳态
    assert np.linalg.norm(F_a) > 0.5


def test_engineering_anchor_mixes_kg_and_kg_s_and_artificially_fills_rank() -> None:
    K_up, K_down = _zeros(), _zeros()
    m_feed = _feed_cell0(1.0)
    r_dot = _zeros()
    M = np.array([0.25, 0.0, 0.0], dtype=np.float64)  # [kg]
    F_eng = holdup_residual(
        M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="engineering"
    )
    F_phys = holdup_residual(
        M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
    )
    # 物理残差 [kg/s]；工程锚 F = source - M，0.75 无一致单位
    np.testing.assert_allclose(F_phys, np.array([1.0, 0.0, 0.0]))
    np.testing.assert_allclose(F_eng, np.array([0.75, 0.0, 0.0]))
    J_eng = holdup_jacobian(
        M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="engineering"
    )
    diag = jacobian_diagnostics(J_eng)
    assert diag.rank == N_CELLS
    assert diag.n_zero_rows == 0
    np.testing.assert_allclose(np.diag(J_eng), -1.0)


def test_first_order_reaction_supplies_legitimate_time_scale_when_K_zero() -> None:
    K_up, K_down = _zeros(), _zeros()
    m_feed = _feed_cell0(1.0)
    r_dot = _zeros()
    k_rxn = np.full(N_CELLS, 0.2, dtype=np.float64)  # [1/s]
    J = holdup_jacobian(
        _zeros(),
        m_feed=m_feed,
        K_up=K_up,
        K_down=K_down,
        r_dot=r_dot,
        k_rxn=k_rxn,
        k0_anchor="none",
    )
    diag = jacobian_diagnostics(J)
    assert diag.rank == N_CELLS
    np.testing.assert_allclose(np.diag(J), -k_rxn)
    M_star = np.linalg.solve(J, -m_feed)
    np.testing.assert_allclose(M_star, np.array([5.0, 0.0, 0.0]))
    F = holdup_residual(
        M_star,
        m_feed=m_feed,
        K_up=K_up,
        K_down=K_down,
        r_dot=r_dot,
        k_rxn=k_rxn,
        k0_anchor="none",
    )
    np.testing.assert_allclose(F, 0.0, atol=1e-12)


def test_analytic_jacobians_match_finite_difference() -> None:
    K_up, K_down = canonical_open_chain_K()
    alpha = alpha_from_K(K_up, K_down)
    m_feed = _feed_cell0(1.0)
    r_dot = np.array([-0.1, 0.0, 0.0], dtype=np.float64)  # [kg/s]
    m_out = np.array([0.9, 0.6, 0.4], dtype=np.float64)
    J_a = stream_jacobian(alpha_up=alpha)
    J_a_fd = _fd_jacobian(
        lambda x: stream_residual(x, m_feed=m_feed, alpha_up=alpha, r_dot=r_dot),
        m_out,
    )
    np.testing.assert_allclose(J_a, J_a_fd, atol=1e-6)

    M = np.array([1.1, 0.7, 0.5], dtype=np.float64)
    k_rxn = np.array([0.05, 0.0, 0.02], dtype=np.float64)
    J_b = holdup_jacobian(
        M,
        m_feed=m_feed,
        K_up=K_up,
        K_down=K_down,
        r_dot=r_dot,
        k_rxn=k_rxn,
        k0_anchor="none",
    )
    J_b_fd = _fd_jacobian(
        lambda x: holdup_residual(
            x,
            m_feed=m_feed,
            K_up=K_up,
            K_down=K_down,
            r_dot=r_dot,
            k_rxn=k_rxn,
            k0_anchor="none",
        ),
        M,
    )
    np.testing.assert_allclose(J_b, J_b_fd, atol=1e-6)


def test_production_holdup_guardrail_matches_engineering_anchor_on_one_bucket() -> None:
    """生产 ``source-M`` 与约化 engineering 锚在单桶 K=0 时同构。"""
    shape = (1, 4)
    m_solid = np.zeros(shape, dtype=np.float64)
    m_solid[0, 0] = 0.25  # [kg] char
    m_auf_in = np.zeros(shape, dtype=np.float64)
    m_auf_in[0, 0] = 1.0  # [kg/s]
    res = calc_solid_balance_residual(
        m_solid_zu=np.zeros(shape, dtype=np.float64),
        m_solid_rez=np.zeros(shape, dtype=np.float64),
        m_solid_in=np.zeros(shape, dtype=np.float64),
        m_solid_auf_in=m_auf_in,
        m_solid_ab_in=np.zeros(shape, dtype=np.float64),
        R_solid=np.zeros(shape, dtype=np.float64),
        size_migration=np.zeros(shape, dtype=np.float64),
        m_solid=m_solid,
        K_solid_auf=np.zeros(shape, dtype=np.float64),
        K_solid_ab=np.zeros(shape, dtype=np.float64),
        solid_state_model="holdup_transport",
    )
    np.testing.assert_allclose(res[0, 0], 0.75)
    F_eng = holdup_residual(
        np.array([0.25, 0.0, 0.0], dtype=np.float64),
        m_feed=np.array([1.0, 0.0, 0.0], dtype=np.float64),
        K_up=_zeros(),
        K_down=_zeros(),
        r_dot=_zeros(),
        k0_anchor="engineering",
    )
    assert F_eng[0] == pytest.approx(0.75)
