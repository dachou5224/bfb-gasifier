"""自由板区颗粒轨迹与夹带：beta_A、u_gb、C_D_haider。

Source: specs/02_hydrodynamics.md §5; docs/CLAUDE.md Phase 2.4
"""

from __future__ import annotations

import math

from src.core.constants import g


def calc_u_gb(u_b: float, scale: float = 1.0) -> float:
    """估算自由板区颗粒特征速度 u_gb。

    Hamel Eq. 3.76：ghost-bubble 初始速度满足 ``u_gb,0 = u_b,ws``。
    因此默认 ``scale=1.0``。颗粒起始速度 ``u_p,0`` 与其不同，
    见 ``src/core/freeboard_segment.py::_build_velocity_samples``。

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


def calc_cd_hamel_eq384(Re_p: float, phi_s: float) -> float:
    """Hamel dissertation Eq. 3.84 drag coefficient transcription.

    c_W = 24/Re * [1 + (8.1716 * exp(-4.0655 * phi_s)) * Re^(0.0964 + 0.5565 * phi_s)]
          + (73.69 * Re * exp(-5.0784 * phi_s)) / (Re + 5.378 * exp(6.2122 * phi_s))
    """
    assert Re_p > 0.0, "Re_p 必须大于 0"
    assert 0.0 < phi_s <= 1.0, "phi_s 必须在 (0,1]"

    phi = float(phi_s)
    re = float(Re_p)
    a = 8.1716 * math.exp(-4.0655 * phi)
    b = 0.0964 + 0.5565 * phi
    num = 73.69 * re * math.exp(-5.0784 * phi)
    den = re + 5.378 * math.exp(6.2122 * phi)
    return (24.0 / re) * (1.0 + a * (re**b)) + num / max(den, 1e-12)


def calc_beta_a(d_b_ws: float) -> float:
    """按 Hamel Gl. 3.75 估算自由板 ghost-bubble 衰减系数 beta_A。

    beta_A = 1 / (3 * d_b,ws)

    其中 ``d_b_ws`` 为气泡离开床层表面时的直径。
    """
    assert d_b_ws > 0.0, "d_b_ws 必须大于 0 [m]"
    return 1.0 / (3.0 * d_b_ws)


def calc_terminal_velocity_haider(
    *,
    rho_g: float,
    rho_p: float,
    d_p: float,
    phi_s: float,
    mu_g: float,
    g_acc: float = g,
    max_iter: int = 80,
) -> float:
    """用 Haider-Levenspiel 阻力关系求自由板颗粒终端沉降速度。

    以重力/浮力与阻力平衡求解：
        0 = g * (1 - rho_g/rho_p) - 3 * C_D * rho_g * u_t^2 / (4 * rho_p * d_p)

    返回正标量 ``u_t`` [m/s]，表示相对气体的终端沉降速度量级。
    """
    assert rho_g > 0.0, "rho_g 必须大于 0 [kg/m3]"
    assert rho_p > rho_g, "rho_p 必须大于 rho_g [kg/m3]"
    assert d_p > 0.0, "d_p 必须大于 0 [m]"
    assert 0.0 < phi_s <= 1.0, "phi_s 必须在 (0,1]"
    assert mu_g > 0.0, "mu_g 必须大于 0 [Pa·s]"
    assert g_acc > 0.0, "g_acc 必须大于 0 [m/s2]"

    ut = max(1e-4, (rho_p - rho_g) * g_acc * d_p * d_p / max(18.0 * mu_g, 1e-12))
    for _ in range(max_iter):
        re_p = max(rho_g * ut * d_p / mu_g, 1e-12)
        cd = calc_cd_haider(re_p, phi_s)
        ut_new = math.sqrt(
            max(4.0 * d_p * (rho_p - rho_g) * g_acc / (3.0 * cd * rho_g), 0.0)
        )
        if abs(ut_new - ut) <= 1e-8 * max(ut_new, 1.0):
            return float(ut_new)
        ut = 0.5 * ut + 0.5 * ut_new
    return float(ut)
