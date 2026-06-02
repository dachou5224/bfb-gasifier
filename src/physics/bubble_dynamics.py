"""气泡动力学：气泡上升速度、气泡直径沿高度增长、慢泡/快泡判别。

提供两套气泡直径模型：
  1) Darton (1977) / Mori-Wen (1975) 经验公式（默认，广泛验证）
  2) Hilligardt (1986) ODE（需要参数 ``xi_b / lambda_b`` 标定）

Source of truth:
- ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``
- Hamel (1999) Eq. 3.14, 3.15, 3.23, 3.41-3.44
- Hilligardt (1986); Heinbockel (1995); Darton et al. (1977); Mori & Wen (1975)
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.integrate import solve_ivp

from src.core.constants import P0_HAMEL, g, n_B_coalescence, n_b


def bubble_interaction_factor(
    u_mf: float | None = None,
    *,
    strategy: str = "technical_distributor",
) -> float:
    """Bubble interaction factor ``psi_b``.

    Supported strategies
    --------------------
    - ``technical_distributor``:
      Hamel (1999) Eq.3.14, p.28 for technical gas distributors
        psi_b = 0.76
    - ``porous_plate``:
      Hamel (1999) Eq.3.14, p.28 for porous plates
        psi_b = 0.67
    - ``wein_1992``:
      Hamel (1999) Eq.3.15, p.28
        psi_b = 0.17 * u_mf^(-0.33)
    """
    name = str(strategy).strip().lower()
    if name == "technical_distributor":
        return 0.76
    if name == "porous_plate":
        return 0.67
    if name == "wein_1992":
        if u_mf is None or u_mf <= 0.0:
            raise ValueError("bubble_interaction_factor(strategy='wein_1992') requires positive u_mf")
        return 0.17 * float(u_mf) ** (-0.33)
    raise ValueError(
        "Unsupported psi_b strategy={!r}; expected 'technical_distributor', "
        "'porous_plate', or 'wein_1992'".format(strategy)
    )


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
    exp_arg = np.clip(-0.3 * h / max(D_bed, 1e-12), -200.0, 200.0)
    return d_bm - (d_bm - d_b0) * np.exp(exp_arg)


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
    *,
    P: float | None = None,
    strategy: str = "hilligardt_eq313",
) -> float:
    """气泡上升速度 u_b [m/s]（Hilligardt 模型）。

    Supported strategies
    --------------------
    - ``hilligardt_eq313``:
        u_b = psi_b * (u0 - u_mf) + u_b,i
    - ``heinbockel_eq343``:
        Hamel (1999) Eq.3.43:
        u_b = psi_b * (u0 - u_mf) * ((P/P0)^0.2 - 1) + v_b * u_b,l
    - ``heinbockel_eq343_legacy_2p14_p07``:
        previously used aggressive pressure form kept only for audit contrast
    - ``heinbockel_eq343_legacy_plus``:
        previous project audit form kept for regression contrast

    Parameters
    ----------
    u0    : 表观气速 [m/s]
    u_mf  : 最小流化速度 [m/s]
    d_b   : 气泡直径 [m]
    psi_b : 气泡相互作用因子 [-]，技术气体分布板默认 0.76；
            原文还给出 porous plate = 0.67、Wein(1992) 关联式

    Source: Hamel (1999) Eq.3.13-3.15, p.28; Eq.3.43 p.31; Hilligardt (1986); Wein (1992); Heinbockel (1995)
    """
    name = str(strategy).strip().lower()
    u_bs = single_bubble_velocity(d_b)
    excess = max(u0 - u_mf, 0.0)
    if name == "hilligardt_eq313":
        return psi_b * excess + u_bs
    if name == "heinbockel_eq343":
        if P is None or P <= 0.0:
            raise ValueError("bubble_rise_velocity(strategy='heinbockel_eq343') requires positive P")
        pressure_term = (P / P0_HAMEL) ** 0.2 - 1.0
        return psi_b * excess * pressure_term + u_bs
    if name == "heinbockel_eq343_legacy_plus":
        if P is None or P <= 0.0:
            raise ValueError("bubble_rise_velocity(strategy='heinbockel_eq343_legacy_plus') requires positive P")
        return psi_b * excess * ((P / P0_HAMEL) ** 0.2 + 1.0) + u_bs
    if name == "heinbockel_eq343_legacy_2p14_p07":
        if P is None or P <= 0.0:
            raise ValueError(
                "bubble_rise_velocity(strategy='heinbockel_eq343_legacy_2p14_p07') requires positive P"
            )
        return psi_b * (2.14 * (P / P0_HAMEL) ** 0.7 * excess + u_bs)
    raise ValueError(
        f"Unsupported bubble velocity strategy={strategy!r}; "
        "expected 'hilligardt_eq313', 'heinbockel_eq343', "
        "'heinbockel_eq343_legacy_plus', or 'heinbockel_eq343_legacy_2p14_p07'"
    )


# -----------------------------------------------------------------------
# 慢泡 / 快泡 判别
# -----------------------------------------------------------------------

def classify_bubble_regime(u_b: float, u_d: float) -> str:
    """根据 alpha_b = u_b / u_d 判别慢泡或快泡状态。

    alpha_b < 1 => 慢泡 (slow)
    alpha_b >= 1 => 快泡 (fast)

    Source: Hamel (1999) Chapter 3.1.2, p.26-27, Eq.3.34 context
    """
    if u_d <= 0:
        return "fast"
    alpha_b = u_b / u_d
    return "slow" if alpha_b < 1.0 else "fast"


# -----------------------------------------------------------------------
# 气泡平均寿命（含加压修正）
# -----------------------------------------------------------------------

def bubble_lifetime(
    d_b: float,
    u_b: float,
    P: float,
    P0_ref: float = P0_HAMEL,
    *,
    strategy: str = "hamel_280",
    u_mf: float | None = None,
) -> float:
    """气泡平均寿命 lambda_b [s]（含加压修正）。

    Supported strategies
    --------------------
    - ``hamel_280``:
      Hamel (1999) Eq.3.35 / Eq.3.42, 引 Heinbockel (1995)：
        lambda_b = 280 * u_mf / g * (P / P0)^(-0.7)
    - ``current``:
      旧项目保留的 legacy 口径（非 Hamel 原文默认）：
        lambda_b = d_b / (0.5 * u_b) * (P / P0)^(-0.2)

    Source: Hamel (1999) Eq.3.35 p.30; Eq.3.42 p.32; Heinbockel (1995)
    """
    name = str(strategy).strip().lower()
    if name == "hamel_280":
        if u_mf is None:
            raise ValueError("bubble_lifetime(strategy='hamel_280') requires u_mf")
        return 280.0 * max(float(u_mf), 0.0) / g * (P / P0_ref) ** (-0.7)
    if u_b <= 0 and name == "current":
        return 1e10
    if name == "current":
        return (d_b / (0.5 * u_b)) * (P / P0_ref) ** (-0.2)
    raise ValueError(
        f"Unsupported lambda strategy={strategy!r}; expected 'current' or 'hamel_280'"
    )


# -----------------------------------------------------------------------
# Hilligardt ODE (备用，参数需标定)
# -----------------------------------------------------------------------

def bubble_diameter_ode(
    h: float,
    y: np.ndarray,
    u0: float,
    u_mf: float,
    P: float,
    u_d: float | None = None,
    psi_b: float = 0.76,
    xi_b: float = 0.35,
    lambda_strategy: str = "hamel_280",
    xi_strategy: str = "hamel_regime",
    velocity_strategy: str = "hilligardt_eq313",
    ode_strategy: str = "hilligardt_eq333",
    d_b_max: float | None = None,
) -> np.ndarray:
    """气泡直径沿高度变化的微分方程 dd_b/dh。

    两种 ODE 策略共享相同的 growth - decay 结构：

    ``hilligardt_eq333`` (常压 Hilligardt 1986)：
        growth = 2/(9π) * ε_b^(1/3) / (1 - ξ_b*(6/π)^(1/3)*ε_b^(1/3))
        decay  = d_b / (3 * λ_b * u_b)
        dd_b/dh = growth - decay

    ``heinbockel_eq341`` (高压 Heinbockel 1995, via Hamel 1999)：
        dd_b/dh = ((2/(9π))^(1/3) * ε_b^(1/3) * (P/P0)^(P0/P))
                  / (1 - ε_b * (P/P0)^(1/3) * ε_b^(1/3))
                  * d_b / (3 * λ_b * u_b)

    Eq.3.42 provides λ_b = 280*u_mf/g*(P/P0)^(-0.7). At high pressure (P >> P0),
    λ_b is very short, making the decay term dominant and capping d_b growth.
    This is physically correct: high pressure suppresses bubble coalescence.

    NOTE: `lambda_strategy='hamel_280'` 是两个分支共用的默认 Eq.3.42 实现。
    IMPORTANT: `chi` 属于 Preto 质量传递模型，不属于本 ODE 的 growth term。

    ``d_b_max`` (slug flow cap)：
        当 d_b ≥ d_b_max 时，强制 dd_b/dh = min(dd_dh, 0) 防止进入段塞流区域。
        Hamel 两相泡流模型在 d_b > 0.6·D_bed（段塞流判据）时物理失效。
        建议取 d_b_max = 0.6 * D_bed（Broadhurst & Becker 1975 段塞流判据）。

    Source: Hamel (1999) Eq.3.33/3.34/3.35/3.41/3.42/3.43/3.44
    """
    d_b = max(y[0], 1e-4)

    u_b = bubble_rise_velocity(u0, u_mf, d_b, psi_b, P=P, strategy=velocity_strategy)
    excess = max(u0 - u_mf, 1e-10)
    u_d_eff = float(max(u_d if u_d is not None else u_mf, 1e-12))
    eps_b = np.clip(excess / u_b, 1e-6, 0.95)

    name = str(ode_strategy).strip().lower()
    if name == "hilligardt_eq333":
        alpha_b = u_b / u_d_eff
        xi_b = resolve_xi_b(alpha_b, strategy=xi_strategy, default_xi=xi_b)

        eps_b_13 = eps_b**(1.0 / 3.0)
        coeff_6pi = (6.0 / np.pi)**(1.0 / 3.0)
        denom = max(1.0 - xi_b * coeff_6pi * eps_b_13, 0.01)
        growth = (2.0 / (9.0 * np.pi)) * eps_b_13 / denom

        lam_b = bubble_lifetime(d_b, u_b, P, strategy=lambda_strategy, u_mf=u_mf)
        decay = d_b / (3.0 * max(lam_b * u_b, 1e-10))
        dd_dh = growth - decay
    elif name == "heinbockel_eq341":
        # Hamel (1999) Eq.3.41 as transcribed in docs/hamel_submodels/03.
        eps_b_13 = eps_b ** (1.0 / 3.0)
        pressure_ratio = max(P / P0_HAMEL, 1e-12)
        denom = 1.0 - eps_b * (pressure_ratio ** (1.0 / 3.0)) * eps_b_13
        denom = float(np.sign(denom) * max(abs(denom), 1e-3))
        lam_b = bubble_lifetime(d_b, u_b, P, strategy=lambda_strategy, u_mf=u_mf)
        dd_dh = (
            ((2.0 / (9.0 * np.pi)) ** (1.0 / 3.0))
            * eps_b_13
            * (pressure_ratio ** (P0_HAMEL / max(P, 1e-12)))
            / denom
            * d_b
            / (3.0 * max(lam_b * u_b, 1e-10))
        )
    elif name == "heinbockel_eq341_legacy_growth_decay":
        # Previous project audit branch retained for regression contrast.
        u_br = n_b * u_d_eff * (P0_HAMEL / P) ** 0.15
        growth = (
            (1.0 / 3.0)
            * (psi_b * max(u0 - u_mf, 0.0) + u_br)
            / max(n_B_coalescence * u_b, 1e-12)
            * (P / P0_HAMEL) ** 0.4
            * (eps_b ** (1.0 / 3.0))
            / max(d_b ** (1.0 / 3.0), 1e-12)
        )
        lam_b = bubble_lifetime(d_b, u_b, P, strategy=lambda_strategy, u_mf=u_mf)
        decay = d_b / (3.0 * max(lam_b * u_b, 1e-10))
        dd_dh = growth - decay
    else:
        raise ValueError(
            f"Unsupported bubble ODE strategy={ode_strategy!r}; "
            "expected 'hilligardt_eq333', 'heinbockel_eq341', "
            "or 'heinbockel_eq341_legacy_growth_decay'"
        )

    if y[0] < 1e-4 and dd_dh < 0:
        dd_dh = 0.0

    # Slug flow cap: when d_b reaches d_b_max, suppress further growth
    # (Broadhurst & Becker 1975: slug flow onset at d_b > 0.6*D_bed)
    if d_b_max is not None and d_b >= d_b_max and dd_dh > 0:
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
    u_d: float | None = None,
    D_bed: float = 0.6,
    d_b0: float | None = None,
    n_points: int = 100,
    method: str = "mori_wen",
    lambda_strategy: str = "hamel_280",
    xi_strategy: str = "hamel_regime",
    velocity_strategy: str = "hilligardt_eq313",
    ode_strategy: str = "hilligardt_eq333",
    psi_b: float = 0.76,
) -> Tuple[np.ndarray, np.ndarray]:
    """沿轴向计算气泡直径 d_b(h)。

    Parameters
    ----------
    u0      : 表观气速 [m/s]
    u_mf    : 最小流化速度 [m/s]
    P       : 压力 [Pa]（仅 ODE 方法使用）
    H_bed   : 床层高度 [m]
    u_d     : Suspensionsphase 内实际气速 [m/s]；仅 ODE/Hamel-regime 使用
    D_bed   : 床层直径 [m]
    d_b0    : 初始气泡直径 [m]，None 则用 Darton 公式
    n_points : 输出点数
    method  : "mori_wen"（默认）或 "darton" 或 "hilligardt_ode"
    lambda_strategy : ODE 模式下的 `lambda_b` 口径（默认 `hamel_280`）
    xi_strategy : ODE 模式下的 `xi_b` 口径（默认 `hamel_regime`）
    velocity_strategy : ODE 模式下的 `u_b` 口径（默认 `hilligardt_eq313`）
    ode_strategy : ODE 主体形式（默认 `hilligardt_eq333`）
    psi_b   : 气泡相互作用因子 [-]，须与 cell 层一致（默认 0.76）

    Returns
    -------
    h_arr  : 高度数组 [m]
    db_arr : 气泡直径数组 [m]

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``
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
        # Slug flow prevention cap: d_b ≤ 0.6*D_bed (Broadhurst & Becker 1975)
        d_b_max = 0.6 * D_bed
        sol = solve_ivp(
            fun=lambda h, y: bubble_diameter_ode(
                h,
                y,
                u0,
                u_mf,
                P,
                u_d=u_d,
                psi_b=psi_b,
                lambda_strategy=lambda_strategy,
                xi_strategy=xi_strategy,
                velocity_strategy=velocity_strategy,
                ode_strategy=ode_strategy,
                d_b_max=d_b_max,
            ),
            t_span=(0, H_bed),
            y0=[d_b0],
            t_eval=h_arr,
            method="RK45",
            rtol=1e-8,
            atol=1e-10,
        )
        if not sol.success:
            raise RuntimeError(f"气泡直径 ODE 积分失败: {sol.message}")
        db_arr = np.minimum(np.maximum(sol.y[0], 1e-4), d_b_max)
    else:
        raise ValueError(f"未知方法: {method}")

    return h_arr, db_arr


def resolve_xi_b(
    alpha_b: float,
    *,
    strategy: str = "fixed_035",
    default_xi: float = 0.35,
) -> float:
    """Resolve Hilligardt ODE shape factor `xi_b`.

    Strategies
    ----------
    - ``fixed_035``:
      keep the current project default / tuning placeholder
    - ``hamel_regime``:
      use Hamel (1999) Chapter 3.1.2, Eq.3.34 with
        alpha_b = u_b / u_d
        xi_b = 1 - alpha_b^3  for 0 < alpha_b < 1
        xi_b = 0              for alpha_b > 1
    """
    name = str(strategy).strip().lower()
    alpha_b = float(max(alpha_b, 0.0))
    if name == "fixed_035":
        return float(default_xi)
    if name == "hamel_regime":
        if alpha_b <= 1.0:
            return float(max(0.0, 1.0 - alpha_b**3))
        return 0.0
    raise ValueError(
        f"Unsupported xi strategy={strategy!r}; expected 'fixed_035' or 'hamel_regime'"
    )
