"""单 cell 守恒方程求解器（scipy.optimize.fsolve）。

给定上游条件（入口摩尔流率、进料），求解当前 cell 的出口状态
（N_b, N_d, m_solid, T），使得所有守恒方程残差为零。

Source: docs/CLAUDE.md Phase 5.2
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import fsolve

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import N_GAS

# 有限差分步长（相对 + 绝对）
_EPS = np.sqrt(np.finfo(float).eps)


def _pack_state(cell: Cell) -> np.ndarray:
    """将 cell 的自由变量打包为一维向量。

    [N_d(N_GAS), N_b(N_GAS), m_solid(nk*4), T(1)]
    """
    return np.concatenate([
        cell.N_d,
        cell.N_b,
        cell.m_solid.flatten(),
        np.array([cell.T]),
    ])


def _unpack_state(x: np.ndarray, cell: Cell) -> None:
    """将一维向量解包回 cell 状态（in-place，避免新建数组）。"""
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    np.maximum(x[:N_GAS], 0.0, out=cell.N_d)
    np.maximum(x[N_GAS:2 * N_GAS], 0.0, out=cell.N_b)
    m_flat = np.maximum(x[2 * N_GAS:2 * N_GAS + nv_sol], 0.0)
    cell.m_solid[:] = m_flat.reshape((nk, N_SOLID_COMP))
    cell.T = max(float(x[2 * N_GAS + nv_sol]), 300.0)



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


from scipy.optimize import least_squares

def solve_cell(
    cell: Cell,
    max_iter: int = 20,
    tol: float = 1e-6,
    verbose: bool = False,
) -> dict:
    """使用 least_squares 稳健地求解单个 cell 的守恒方程。"""
    x0 = _pack_state(cell)
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    n_vars = len(x0)

    # 设定边界：[N_d, N_b, m_solid, T]
    lower_bounds = np.zeros(n_vars)
    lower_bounds[-1] = 300.0  # T_min
    upper_bounds = np.full(n_vars, np.inf)
    upper_bounds[-1] = 3000.0  # T_max

    # 静态缩放因子（基于进料和能量量级）
    diag = np.ones(n_vars)
    diag[:2 * N_GAS] = 10.0    # 气相流量量级 ~10
    diag[2 * N_GAS : 2 * N_GAS + nv_sol] = 0.5  # 固相量级 ~0.5
    diag[-1] = 1000.0          # 温度量级 ~1000

    res = least_squares(
        _residual_wrapper,
        x0,
        args=(cell,),
        jac=_jacobian_wrapper,
        bounds=(lower_bounds, upper_bounds),
        x_scale=diag,
        ftol=1e-4,  # 放宽门限以允许更大步进
        xtol=1e-4,
        max_nfev=max_iter * n_vars,
        diff_step=0.01,  # 强制较大的有限差分步长
        verbose=2 if verbose else 0,
    )

    _unpack_state(res.x, cell)
    final_res = cell.residuals()
    max_res = float(np.max(np.abs(final_res)))

    return {
        "converged": res.success,
        "residual": max_res,
        "status": res.status,
        "message": res.message,
    }
