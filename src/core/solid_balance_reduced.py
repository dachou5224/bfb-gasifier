"""三格、单物种、单粒径固相稳态质量平衡的约化问题。

本模块不依赖 ``Reactor`` / ``Cell``。它只比较两种数学提法：

- Model A（流平衡）：未知量是格间质量流 ``m_out`` [kg/s]。
- Model B（库存平衡）：未知量是格内存量 ``M`` [kg]，出流由 ``m = K M`` 闭合。

当 ``K_up = K_down = 0`` 时，Model B 的稳态方程不再约束 ``M``。
生产代码里的 ``source - M`` 把 [kg/s] 与 [kg] 相减，等价于偷偷引入
``τ = 1 s`` 的出流本构；这不是 Wirsum/Hamel 的质量守恒。

Source: Wirsum (1997) Kap. 2.3; Hamel (1999) Eq. (2-6)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

N_CELLS: int = 3
# 判定零行/零列：相对矩阵尺度的绝对阈值
_ZERO_ABS_TOL: float = 1.0e-15
# 与 cell_balances._apply_holdup_transport_numerical_guardrails 相同
_K0_ATOL: float = 1.0e-14  # [1/s]
_SRC_ATOL: float = 1.0e-14  # [kg/s]


@dataclass(frozen=True)
class JacobianDiagnostics:
    """稠密雅可比的秩与零空间诊断。"""

    rank: int
    n: int
    nullity: int
    n_zero_rows: int
    n_zero_cols: int
    cond: float
    singular_values: npt.NDArray[np.float64]


def jacobian_diagnostics(J: npt.NDArray[np.float64]) -> JacobianDiagnostics:
    """对稠密 ``J`` 计算秩、零行/列数与条件数。"""
    assert J.ndim == 2 and J.shape[0] == J.shape[1], "J 必须是方阵"
    n = int(J.shape[0])
    scale = float(np.max(np.abs(J))) if J.size else 0.0
    row_zero = np.max(np.abs(J), axis=1) <= max(_ZERO_ABS_TOL, scale * _ZERO_ABS_TOL)
    col_zero = np.max(np.abs(J), axis=0) <= max(_ZERO_ABS_TOL, scale * _ZERO_ABS_TOL)
    if scale == 0.0:
        s = np.zeros(n, dtype=np.float64)
        rank = 0
        cond = float("inf")
    else:
        s = np.linalg.svd(J, compute_uv=False)
        tol = float(n) * np.finfo(np.float64).eps * float(s[0])
        rank = int(np.sum(s > tol))
        cond = float(s[0] / s[-1]) if s[-1] > 0.0 else float("inf")
    return JacobianDiagnostics(
        rank=rank,
        n=n,
        nullity=n - rank,
        n_zero_rows=int(np.sum(row_zero)),
        n_zero_cols=int(np.sum(col_zero)),
        cond=cond,
        singular_values=np.asarray(s, dtype=np.float64),
    )


def _split_outflow(
    m_out: npt.NDArray[np.float64],
    alpha_up: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """由总出流与上向分流比得到 ``m_up`` / ``m_down`` [kg/s]。"""
    assert m_out.shape == (N_CELLS,)
    assert alpha_up.shape == (N_CELLS,)
    m_up = alpha_up * m_out
    m_down = (1.0 - alpha_up) * m_out
    return m_up, m_down


def stream_residual(
    m_out: npt.NDArray[np.float64],
    *,
    m_feed: npt.NDArray[np.float64],
    alpha_up: npt.NDArray[np.float64],
    r_dot: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Model A：流平衡残差 [kg/s]。

    格 ``i``::

        F_i = m_feed_i + m_up_{i-1} + m_down_{i+1} - m_out_i + r_dot_i

    边界：床底无下向入流，床顶无下向来流（``m_down`` 出床顶视为离开系统）。
    ``r_dot`` 为已知反应汇 [kg/s]，不依赖 ``m_out``。
    """
    m_out = np.asarray(m_out, dtype=np.float64)
    m_feed = np.asarray(m_feed, dtype=np.float64)
    alpha_up = np.asarray(alpha_up, dtype=np.float64)
    r_dot = np.asarray(r_dot, dtype=np.float64)
    assert m_out.shape == (N_CELLS,)
    assert m_feed.shape == (N_CELLS,)
    assert alpha_up.shape == (N_CELLS,)
    assert r_dot.shape == (N_CELLS,)
    m_up, m_down = _split_outflow(m_out, alpha_up)
    F = np.zeros(N_CELLS, dtype=np.float64)
    for i in range(N_CELLS):
        inflow = m_feed[i]
        if i > 0:
            inflow += m_up[i - 1]
        if i + 1 < N_CELLS:
            inflow += m_down[i + 1]
        F[i] = inflow - m_out[i] + r_dot[i]
    return F


def stream_jacobian(
    *,
    alpha_up: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Model A 对 ``m_out`` 的解析雅可比 [无量纲： (kg/s)/(kg/s)]。"""
    alpha_up = np.asarray(alpha_up, dtype=np.float64)
    assert alpha_up.shape == (N_CELLS,)
    J = np.zeros((N_CELLS, N_CELLS), dtype=np.float64)
    for i in range(N_CELLS):
        J[i, i] -= 1.0
        if i > 0:
            J[i, i - 1] += alpha_up[i - 1]
        if i + 1 < N_CELLS:
            J[i, i + 1] += 1.0 - alpha_up[i + 1]
    return J


def holdup_residual(
    M: npt.NDArray[np.float64],
    *,
    m_feed: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    K_down: npt.NDArray[np.float64],
    r_dot: npt.NDArray[np.float64],
    k_rxn: npt.NDArray[np.float64] | None = None,
    k0_anchor: str = "none",
) -> npt.NDArray[np.float64]:
    """Model B：库存平衡残差。

    物理残差（``k0_anchor='none'``）单位为 [kg/s]::

        F_i = m_feed_i + K_up_{i-1} M_{i-1} + K_down_{i+1} M_{i+1}
              - (K_up_i + K_down_i) M_i + r_dot_i - k_rxn_i M_i

    ``k0_anchor='engineering'`` 复现 ``cell_balances._apply_holdup_transport_numerical_guardrails``：
    当 ``K_up+K_down=0`` 时用 ``F = net_source - M`` 或 ``F = -M``，
    把 [kg/s] 与 [kg] 相减。
    """
    M = np.asarray(M, dtype=np.float64)
    m_feed = np.asarray(m_feed, dtype=np.float64)
    K_up = np.asarray(K_up, dtype=np.float64)
    K_down = np.asarray(K_down, dtype=np.float64)
    r_dot = np.asarray(r_dot, dtype=np.float64)
    assert M.shape == (N_CELLS,)
    assert m_feed.shape == (N_CELLS,)
    assert K_up.shape == (N_CELLS,)
    assert K_down.shape == (N_CELLS,)
    assert r_dot.shape == (N_CELLS,)
    if k_rxn is None:
        k_rxn = np.zeros(N_CELLS, dtype=np.float64)
    else:
        k_rxn = np.asarray(k_rxn, dtype=np.float64)
        assert k_rxn.shape == (N_CELLS,)
    if k0_anchor not in {"none", "engineering"}:
        raise ValueError(f"未知 k0_anchor={k0_anchor!r}")

    F = np.zeros(N_CELLS, dtype=np.float64)
    for i in range(N_CELLS):
        inflow = m_feed[i]
        if i > 0:
            inflow += K_up[i - 1] * M[i - 1]
        if i + 1 < N_CELLS:
            inflow += K_down[i + 1] * M[i + 1]
        K_sum = K_up[i] + K_down[i]
        F[i] = inflow - K_sum * M[i] + r_dot[i] - k_rxn[i] * M[i]

    if k0_anchor == "none":
        return F

    # 工程锚：与 src/core/cell_balances.py 同构
    K_sum = K_up + K_down
    disconnected = (K_sum <= _K0_ATOL) & (np.abs(F) <= _SRC_ATOL)
    source_no_transport = (K_sum <= _K0_ATOL) & (np.abs(F) > _SRC_ATOL)
    F = np.where(disconnected, -M, F)
    F = np.where(source_no_transport, F - M, F)
    return F


def holdup_jacobian(
    M: npt.NDArray[np.float64],
    *,
    m_feed: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    K_down: npt.NDArray[np.float64],
    r_dot: npt.NDArray[np.float64],
    k_rxn: npt.NDArray[np.float64] | None = None,
    k0_anchor: str = "none",
) -> npt.NDArray[np.float64]:
    """Model B 对 ``M`` 的解析雅可比。

    ``k0_anchor='none'`` 时单位为 [1/s]。
    ``engineering`` 锚点行的对角元为 ``-1``，单位变为 [1] 或混用。
    """
    M = np.asarray(M, dtype=np.float64)
    if k_rxn is None:
        k_rxn = np.zeros(N_CELLS, dtype=np.float64)
    else:
        k_rxn = np.asarray(k_rxn, dtype=np.float64)
    K_up = np.asarray(K_up, dtype=np.float64)
    K_down = np.asarray(K_down, dtype=np.float64)
    J = np.zeros((N_CELLS, N_CELLS), dtype=np.float64)
    for i in range(N_CELLS):
        J[i, i] -= (K_up[i] + K_down[i] + k_rxn[i])
        if i > 0:
            J[i, i - 1] += K_up[i - 1]
        if i + 1 < N_CELLS:
            J[i, i + 1] += K_down[i + 1]

    if k0_anchor == "none":
        return J

    F = holdup_residual(
        M,
        m_feed=m_feed,
        K_up=K_up,
        K_down=K_down,
        r_dot=r_dot,
        k_rxn=k_rxn,
        k0_anchor="none",
    )
    K_sum = K_up + K_down
    disconnected = (K_sum <= _K0_ATOL) & (np.abs(F) <= _SRC_ATOL)
    source_no_transport = (K_sum <= _K0_ATOL) & (np.abs(F) > _SRC_ATOL)
    for i in range(N_CELLS):
        if disconnected[i]:
            J[i, :] = 0.0
            J[i, i] = -1.0
        elif source_no_transport[i]:
            J[i, i] = -1.0
    return J


def canonical_open_chain_K() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """开链传输：床底只向上，床顶离开系统，中间双向。

    返回 ``(K_up, K_down)``，单位 [1/s]。
    """
    K_up = np.array([0.4, 0.3, 0.5], dtype=np.float64)  # [1/s]
    K_down = np.array([0.0, 0.2, 0.0], dtype=np.float64)  # [1/s]
    return K_up, K_down


def alpha_from_K(
    K_up: npt.NDArray[np.float64],
    K_down: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """``alpha_up = K_up / (K_up+K_down)``；``K_sum=0`` 时取 1（约定出流向上离开）。"""
    K_sum = K_up + K_down
    alpha = np.ones(N_CELLS, dtype=np.float64)
    mask = K_sum > 0.0
    alpha[mask] = K_up[mask] / K_sum[mask]
    return alpha
