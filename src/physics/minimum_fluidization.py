"""最小流化速度 u_mf 计算。

通过 Archimedes 数和 Ergun 方程求解 Re_mf，进而计算 u_mf。

Source: specs/02_hydrodynamics.md §1; Hamel (1999) Eq.4.2; Ergun (1952)
"""

from __future__ import annotations

import numpy as np

from src.core.constants import g


def calc_archimedes(
    rho_g: float,
    rho_s: float,
    d_p: float,
    mu_g: float,
) -> float:
    """计算 Archimedes 数 Ar。

    Ar = rho_g * (rho_s - rho_g) * g * d_p^3 / mu_g^2

    Parameters
    ----------
    rho_g : 气体密度 [kg/m³]
    rho_s : 固体颗粒密度 [kg/m³]
    d_p   : 颗粒直径 [m]
    mu_g  : 气体动力粘度 [Pa·s]

    Returns
    -------
    Ar : Archimedes 数 [-]

    Source: Hamel (1999); specs/02_hydrodynamics.md §1
    """
    assert d_p > 0 and mu_g > 0 and rho_g > 0, "输入参数必须为正"
    return rho_g * (rho_s - rho_g) * g * d_p**3 / mu_g**2


def calc_re_mf(
    Ar: float,
    eps_mf: float = 0.45,
    phi_s: float = 0.86,
) -> float:
    """求解 Ergun 方程得到 Re_mf。

    Ergun 方程：
      Ar = [150*(1 - eps_mf) / (phi_s * eps_mf^3)] * Re_mf
         + [1.75 / (phi_s * eps_mf^3)] * Re_mf^2

    这是关于 Re_mf 的二次方程 a*x^2 + b*x - Ar = 0，取正根。

    Parameters
    ----------
    Ar     : Archimedes 数 [-]
    eps_mf : 最小流化空隙率 [-]，默认 0.45
    phi_s  : 颗粒球形度 [-]，默认 0.86

    Returns
    -------
    Re_mf : 最小流化 Reynolds 数 [-]

    Source: specs/02_hydrodynamics.md §1, Eq.4.2; Ergun (1952)
    """
    assert Ar >= 0, f"Ar must be non-negative, got {Ar}"
    e3 = eps_mf**3
    b_coeff = 150.0 * (1.0 - eps_mf) / (phi_s * e3)
    a_coeff = 1.75 / (phi_s * e3)
    discriminant = b_coeff**2 + 4.0 * a_coeff * Ar
    return (-b_coeff + np.sqrt(discriminant)) / (2.0 * a_coeff)


def calc_u_mf(
    Re_mf: float,
    rho_g: float,
    mu_g: float,
    d_p: float,
) -> float:
    """由 Re_mf 计算最小流化速度 u_mf。

    u_mf = Re_mf * mu_g / (rho_g * d_p)

    Parameters
    ----------
    Re_mf : 最小流化 Reynolds 数 [-]
    rho_g : 气体密度 [kg/m³]
    mu_g  : 气体动力粘度 [Pa·s]
    d_p   : 颗粒直径 [m]

    Returns
    -------
    u_mf : 最小流化速度 [m/s]

    Source: specs/02_hydrodynamics.md §1
    """
    return Re_mf * mu_g / (rho_g * d_p)


def compute_u_mf(
    rho_g: float,
    rho_s: float,
    d_p: float,
    mu_g: float,
    eps_mf: float = 0.45,
    phi_s: float = 0.86,
) -> float:
    """一步计算最小流化速度 u_mf（便捷函数）。

    Parameters
    ----------
    rho_g  : 气体密度 [kg/m³]
    rho_s  : 固体颗粒密度 [kg/m³]
    d_p    : 颗粒直径 [m]
    mu_g   : 气体动力粘度 [Pa·s]
    eps_mf : 最小流化空隙率 [-]
    phi_s  : 颗粒球形度 [-]

    Returns
    -------
    u_mf : 最小流化速度 [m/s]

    Source: specs/02_hydrodynamics.md §1
    """
    Ar = calc_archimedes(rho_g, rho_s, d_p, mu_g)
    Re_mf = calc_re_mf(Ar, eps_mf, phi_s)
    return calc_u_mf(Re_mf, rho_g, mu_g, d_p)
