"""单 cell 守恒方程求解器（scipy.optimize.fsolve）。

给定上游条件（入口摩尔流率、进料），求解当前 cell 的出口状态
（N_b, N_d, m_solid, T），使得所有守恒方程残差为零。

Source: docs/CLAUDE.md Phase 5.2
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import fsolve

from src.core.cell import Cell
from src.core.species import N_GAS

# 有限差分步长（相对 + 绝对）
_EPS = np.sqrt(np.finfo(float).eps)


def _pack_state(cell: Cell) -> np.ndarray:
    """将 cell 的自由变量打包为一维向量。

    [N_d(N_GAS), N_b(N_GAS), m_solid(nk), T(1)]
    """
    return np.concatenate([
        cell.N_d,
        cell.N_b,
        cell.m_solid,
        np.array([cell.T]),
    ])


def _unpack_state(x: np.ndarray, cell: Cell) -> None:
    """将一维向量解包回 cell 状态（in-place，避免新建数组）。"""
    nk = cell.solid.n_size_classes
    np.maximum(x[:N_GAS], 0.0, out=cell.N_d)
    np.maximum(x[N_GAS:2 * N_GAS], 0.0, out=cell.N_b)
    np.maximum(x[2 * N_GAS:2 * N_GAS + nk], 0.0, out=cell.m_solid)
    cell.T = max(float(x[2 * N_GAS + nk]), 300.0)


def _residual_wrapper(x: np.ndarray, cell: Cell) -> np.ndarray:
    """fsolve 的目标函数：将 x 解包后计算残差。"""
    _unpack_state(x, cell)
    return cell.residuals()


def _jacobian_wrapper(x: np.ndarray, cell: Cell) -> np.ndarray:
    """数值 Jacobian：d(residual)/d(x)，供 fsolve 使用以替代内部有限差分。"""
    n = len(x)
    res0 = _residual_wrapper(x, cell)
    eps = _EPS * (1.0 + np.abs(x))
    eps = np.maximum(eps, 1e-12)
    J = np.zeros((n, n))
    for j in range(n):
        x_plus = x.copy()
        x_plus[j] += eps[j]
        res_plus = _residual_wrapper(x_plus, cell)
        J[:, j] = (res_plus - res0) / eps[j]
    _unpack_state(x, cell)  # 恢复 cell 状态
    return J


def solve_cell(
    cell: Cell,
    max_iter: int = 500,
    tol: float = 1e-8,
) -> dict:
    """求解单个 cell 的守恒方程。

    Parameters
    ----------
    cell     : 已设定入口条件和进料的 Cell 对象
    max_iter : fsolve 最大迭代次数
    tol      : 收敛容差

    Returns
    -------
    dict with keys:
      'converged' : bool
      'residual'  : float (最大残差绝对值)
      'info'      : fsolve info dict

    Source: docs/CLAUDE.md Phase 5.2
    """
    x0 = _pack_state(cell)

    # 确保初始猜测合理
    x0 = np.where(np.isnan(x0), 0.0, x0)
    x0 = np.where(np.isinf(x0), 0.0, x0)

    sol, info, ier, msg = fsolve(
        _residual_wrapper,
        x0,
        args=(cell,),
        fprime=_jacobian_wrapper,
        full_output=True,
        maxfev=max_iter * (len(x0) + 1),
        xtol=tol,
    )

    _unpack_state(sol, cell)
    final_res = cell.residuals()
    max_res = float(np.max(np.abs(final_res)))

    return {
        "converged": ier == 1 or max_res < tol * 100,
        "residual": max_res,
        "info": info,
        "message": msg,
    }
