"""自由板区颗粒轨迹与夹带：beta_A、u_gb、C_D_haider。

Source: specs/02_hydrodynamics.md §5; docs/CLAUDE.md Phase 2.4
"""

from __future__ import annotations

import math


def calc_u_gb(u_b: float, scale: float = 1.53) -> float:
    """估算自由板区颗粒特征速度 u_gb。

    默认采用 techspec 建议均值：u_p0 = 1.53 * u_b。

    Source: docs/CLAUDE.md §当前已知缺口(4)
    """
    assert u_b >= 0.0, "u_b 必须非负 [m/s]"
    assert scale > 0.0, "scale 必须大于 0"
    return scale * u_b


def calc_cd_haider(Re_p: float, phi_s: float) -> float:
    """Haider-Levenspiel 非球形颗粒阻力系数。

    C_D = 24/Re * (1 + A*Re^B) + C/(1 + D/Re)

    Source: Haider & Levenspiel (1989); docs/techspec.md §4.4（自由板区）
    """
    assert Re_p > 0.0, "Re_p 必须大于 0"
    assert 0.0 < phi_s <= 1.0, "phi_s 必须在 (0,1]"

    A = math.exp(2.3288 - 6.4581 * phi_s + 2.4486 * (phi_s**2))
    B = 0.0964 + 0.5565 * phi_s
    C = math.exp(4.9050 - 13.8944 * phi_s + 18.4222 * (phi_s**2) - 10.2599 * (phi_s**3))
    D = math.exp(1.4681 + 12.2584 * phi_s - 20.7322 * (phi_s**2) + 15.8855 * (phi_s**3))
    return (24.0 / Re_p) * (1.0 + A * (Re_p**B)) + C / (1.0 + D / Re_p)


def calc_beta_a(C_d: float, rho_g: float, rho_p: float, d_p: float) -> float:
    """估算自由板区轴向衰减系数 beta_A。

    采用单颗粒阻力平衡近似：
    beta_A = 3 * C_D * rho_g / (4 * rho_p * d_p)

    Source: docs/techspec.md §4（自由板区，工程近似）
    """
    assert C_d > 0.0, "C_d 必须大于 0"
    assert rho_g > 0.0 and rho_p > 0.0, "rho_g/rho_p 必须大于 0 [kg/m^3]"
    assert d_p > 0.0, "d_p 必须大于 0 [m]"
    return 3.0 * C_d * rho_g / (4.0 * rho_p * d_p)
