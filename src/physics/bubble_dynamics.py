"""气泡动力学：气泡上升速度、气泡直径沿高度增长、慢泡/快泡判别。

提供两套气泡直径模型：
  1) Darton (1977) / Mori-Wen (1975) 经验公式（默认，广泛验证）
  2) Hilligardt (1986) ODE（需要参数 xi_b / lambda_b 标定，保留为备用）

Source: specs/02_hydrodynamics.md §2; Hamel (1999) Eq.4.4-4.6;
        Hilligardt (1986); Darton et al. (1977); Mori & Wen (1975)
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.integrate import solve_ivp

from src.core.constants import g, P0


# -----------------------------------------------------------------------
# 气泡初始直径
# -----------------------------------------------------------------------

def initial_bubble_diameter(
    u0: float,
    u_mf: float,
    A_bed: float = 0.283,
    N_or: int = 100,
) -> float:
    """气泡初始直径 d_b0 [m]（Darton 1977）。

    d_b0 = 0.54/g^0.2 * (u0-u_mf)^0.4 * (4*sqrt(A_bed/N_or))^0.8

    Parameters
    ----------
    u0    : 表观气速 [m/s]
    u_mf  : 最小流化速度 [m/s]
    A_bed : 床层截面积 [m^2] (TODO: 待用户确认，默认 D=0.6m)
    N_or  : 分布板孔数 (TODO: 待用户确认，默认 100)

    Source: Darton et al. (1977); Hamel (1999)
    """
    excess = max(u0 - u_mf, 1e-6)
    h0 = 4.0 * np.sqrt(A_bed / max(N_or, 1))
    return max(0.54 / g**0.2 * excess**0.4 * h0**0.8, 1e-3)


# -----------------------------------------------------------------------
# 最大稳定气泡直径
# -----------------------------------------------------------------------

def max_bubble_diameter(
    u0: float,
    u_mf: float,
    D_bed: float = 0.6,
) -> float:
    """Mori-Wen (1975) 最大稳定气泡直径 d_bm [m]。

    d_bm = 0.652 * [A_bed * (u0 - u_mf)]^0.4

    Parameters
    ----------
    u0    : 表观气速 [m/s]
    u_mf  : 最小流化速度 [m/s]
    D_bed : 床层直径 [m]

    Source: Mori & Wen (1975)
    """
    excess = max(u0 - u_mf, 1e-6)
    A_bed = np.pi / 4.0 * D_bed**2
    return 0.652 * (A_bed * excess)**0.4


# -----------------------------------------------------------------------
# 气泡直径沿高度分布
# -----------------------------------------------------------------------

def darton_bubble_diameter(
    h: float,
    u0: float,
    u_mf: float,
    A_bed: float = 0.283,
    N_or: int = 100,
) -> float:
    """Darton (1977) 气泡直径 d_b(h) [m]。

    d_b = 0.54/g^0.2 * (u0-u_mf)^0.4 * (h + 4*sqrt(A_bed/N_or))^0.8

    Source: Darton et al. (1977)
    """
    excess = max(u0 - u_mf, 1e-6)
    h_eff = h + 4.0 * np.sqrt(A_bed / max(N_or, 1))
    return 0.54 / g**0.2 * excess**0.4 * h_eff**0.8


def mori_wen_bubble_diameter(
    h: float,
    u0: float,
    u_mf: float,
    D_bed: float = 0.6,
    d_b0: float | None = None,
    A_bed: float = 0.283,
    N_or: int = 100,
) -> float:
    """Mori-Wen (1975) 气泡直径 d_b(h) [m]。

    d_b(h) = d_bm - (d_bm - d_b0) * exp(-0.3*h/D_bed)

    Source: Mori & Wen (1975)
    """
    d_bm = max_bubble_diameter(u0, u_mf, D_bed)
    if d_b0 is None:
        d_b0 = initial_bubble_diameter(u0, u_mf, A_bed, N_or)
    return d_bm - (d_bm - d_b0) * np.exp(-0.3 * h / D_bed)


# -----------------------------------------------------------------------
# 单气泡上升速度
# -----------------------------------------------------------------------

def single_bubble_velocity(d_b: float) -> float:
    """单气泡（孤立）上升速度 u_b,single [m/s]。

    Davies & Taylor (1950):
      u_b,single = 0.711 * sqrt(g * d_b)

    Source: Hamel (1999) Eq.4.5
    """
    return 0.711 * np.sqrt(g * max(d_b, 0.0))


def bubble_rise_velocity(
    u0: float,
    u_mf: float,
    d_b: float,
    psi_b: float = 0.76,
) -> float:
    """气泡上升速度 u_b [m/s]（Hilligardt 模型）。

    u_b = psi_b * (u0 - u_mf) + u_b,single(d_b)

    Parameters
    ----------
    u0    : 表观气速 [m/s]
    u_mf  : 最小流化速度 [m/s]
    d_b   : 气泡直径 [m]
    psi_b : 气泡相互作用因子 [-]，工业分布板默认 0.76

    Source: specs/02_hydrodynamics.md §2, Eq.4.4; Hilligardt (1986)
    """
    u_bs = single_bubble_velocity(d_b)
    return psi_b * max(u0 - u_mf, 0.0) + u_bs


# -----------------------------------------------------------------------
# 慢泡 / 快泡 判别
# -----------------------------------------------------------------------

def classify_bubble_regime(u_b: float, u_mf: float) -> str:
    """根据 alpha_b = u_b / u_mf 判别慢泡或快泡状态。

    alpha_b < 1 => 慢泡 (slow)
    alpha_b >= 1 => 快泡 (fast)

    Source: specs/02_hydrodynamics.md §2; Hamel (1999) Eq.1.2
    """
    if u_mf <= 0:
        return "fast"
    alpha_b = u_b / u_mf
    return "slow" if alpha_b < 1.0 else "fast"


# -----------------------------------------------------------------------
# 气泡平均寿命（含加压修正）
# -----------------------------------------------------------------------

def bubble_lifetime(
    d_b: float,
    u_b: float,
    P: float,
    P0_ref: float = P0,
) -> float:
    """气泡平均寿命 lambda_b [s]（含加压修正）。

    Hilligardt (1986) 简化形式：
      lambda_b = d_b / (0.5 * u_b) * (P / P0)^(-0.2)

    Source: Hilligardt (1986); Hamel (1999) Eq.4.6 参数
    """
    if u_b <= 0:
        return 1e10
    return (d_b / (0.5 * u_b)) * (P / P0_ref)**(-0.2)


# -----------------------------------------------------------------------
# Hilligardt ODE (备用，参数需标定)
# -----------------------------------------------------------------------

def bubble_diameter_ode(
    h: float,
    y: np.ndarray,
    u0: float,
    u_mf: float,
    P: float,
    psi_b: float = 0.76,
    xi_b: float = 0.35,
) -> np.ndarray:
    """气泡直径沿高度变化的微分方程 dd_b/dh。

    dd_b/dh = [2/(9*pi) * eps_b^(1/3) / (1 - xi_b*(6/pi)^(1/3)*eps_b^(1/3))]
              - d_b / (3 * lambda_b * u_b)

    NOTE: 此 ODE 的 lambda_b 参数在高压（>1 MPa）下需要额外标定，
    建议优先使用 mori_wen_bubble_diameter() 或 darton_bubble_diameter()。

    Source: specs/02_hydrodynamics.md §2, Eq.4.6; Hamel (1999)
    """
    d_b = max(y[0], 1e-4)

    u_b = bubble_rise_velocity(u0, u_mf, d_b, psi_b)
    excess = max(u0 - u_mf, 1e-10)
    eps_b = np.clip(excess / u_b, 1e-6, 0.8)

    eps_b_13 = eps_b**(1.0 / 3.0)
    coeff_6pi = (6.0 / np.pi)**(1.0 / 3.0)
    denom = max(1.0 - xi_b * coeff_6pi * eps_b_13, 0.01)
    growth = (2.0 / (9.0 * np.pi)) * eps_b_13 / denom

    lam_b = bubble_lifetime(d_b, u_b, P)
    decay = d_b / (3.0 * max(lam_b * u_b, 1e-10))

    dd_dh = growth - decay
    if y[0] < 1e-4 and dd_dh < 0:
        dd_dh = 0.0
    return np.array([dd_dh])


# -----------------------------------------------------------------------
# 气泡直径沿高度积分（默认用 Mori-Wen 解析）
# -----------------------------------------------------------------------

def integrate_bubble_diameter(
    u0: float,
    u_mf: float,
    P: float,
    H_bed: float,
    D_bed: float = 0.6,
    d_b0: float | None = None,
    n_points: int = 100,
    method: str = "mori_wen",
) -> Tuple[np.ndarray, np.ndarray]:
    """沿轴向计算气泡直径 d_b(h)。

    Parameters
    ----------
    u0      : 表观气速 [m/s]
    u_mf    : 最小流化速度 [m/s]
    P       : 压力 [Pa]（仅 ODE 方法使用）
    H_bed   : 床层高度 [m]
    D_bed   : 床层直径 [m]
    d_b0    : 初始气泡直径 [m]，None 则用 Darton 公式
    n_points : 输出点数
    method  : "mori_wen"（默认）或 "darton" 或 "hilligardt_ode"

    Returns
    -------
    h_arr  : 高度数组 [m]
    db_arr : 气泡直径数组 [m]

    Source: specs/02_hydrodynamics.md §2
    """
    A_bed = np.pi / 4.0 * D_bed**2
    if d_b0 is None:
        d_b0 = initial_bubble_diameter(u0, u_mf, A_bed)

    h_arr = np.linspace(0, H_bed, n_points)

    if method == "mori_wen":
        db_arr = np.array([
            mori_wen_bubble_diameter(h, u0, u_mf, D_bed, d_b0, A_bed)
            for h in h_arr
        ])
    elif method == "darton":
        db_arr = np.array([
            darton_bubble_diameter(h, u0, u_mf, A_bed)
            for h in h_arr
        ])
    elif method == "hilligardt_ode":
        d_b0 = max(d_b0, 1e-3)
        sol = solve_ivp(
            fun=lambda h, y: bubble_diameter_ode(h, y, u0, u_mf, P),
            t_span=(0, H_bed),
            y0=[d_b0],
            t_eval=h_arr,
            method="RK45",
            rtol=1e-8,
            atol=1e-10,
        )
        if not sol.success:
            raise RuntimeError(f"气泡直径 ODE 积分失败: {sol.message}")
        db_arr = np.maximum(sol.y[0], 1e-4)
    else:
        raise ValueError(f"未知方法: {method}")

    return h_arr, db_arr
