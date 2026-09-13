"""Wirsum (1997) §3.3–3.4 自由板内核：模型名、抛射、Roof 3.59、解析轨迹。

本模块是 thesis-exact / engineering 的对照核，不接入生产 ``simulate_freeboard``。
现有 ``analytical_wirsum`` 是工程混合路径的别名，不是 Wirsum 原文。

Source: Wirsum (1997) §3.3–3.4 Eq. 3.59; Hamel (1999) Eq. 3.53（Roof 同构）,
        Eq. 3.73（Martens 1.53 u_b，工程/Hamel 扩展）
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt

from src.core.constants import g

FREEBOARD_MODEL_WIRSUM_1997_EXACT: str = "wirsum_1997_exact"
FREEBOARD_MODEL_ENGINEERING: str = "engineering_freeboard"
FREEBOARD_MODEL_HAMEL_EXTENDED: str = "hamel_extended"
FREEBOARD_MODEL_LEGACY_ALIAS: str = "analytical_wirsum"

WIRSUM_1997_U_P0_OVER_U_B: float = 1.2  # [-] Wirsum §3.4 平均抛射速度
MARTENS_U_P0_OVER_U_B: float = 1.53  # [-] Hamel Eq. 3.73 / Martens (1984)
WIRSUM_1997_VELOCITY_SIGMA: float = 0.60  # [m/s] 竖直抛射速度高斯标准差
WIRSUM_1997_Z_DRAG_MIN: float = 0.43  # [-] 高 Re 阻力极限下限
WIRSUM_1997_Z_DRAG_MAX: float = 1.0  # [-] 高 Re 阻力极限上限
WIRSUM_1997_Z_DRAG_DEFAULT: float = 0.44  # [-] 球形颗粒常用值，落在 [0.43, 1.0]

# Pemberton & Davidson wake（Hamel Eq. 3.54）；仅 engineering / hamel_extended
WAKE_ZETA_W: float = 0.40  # [-]
WAKE_F_W: float = 0.25  # [-]

EjectionKind = Literal["roof_eq359", "roof_plus_wake", "roof_plus_wake_renormalize"]
DragKind = Literal["wirsum_1997_cw", "hamel_z_exp", "haider_mixed"]
ParticleFate = Literal["reaches_next_cell", "returns", "undecided"]


@dataclass(frozen=True)
class FreeboardModelSpec:
    """自由板闭合的具名规格。生产路径尚未读取本对象。"""

    name: str
    u_p0_over_u_b: float  # [-]
    velocity_sigma: float  # [m/s]
    ejection: EjectionKind
    drag: DragKind


_SPECS: dict[str, FreeboardModelSpec] = {
    FREEBOARD_MODEL_WIRSUM_1997_EXACT: FreeboardModelSpec(
        name=FREEBOARD_MODEL_WIRSUM_1997_EXACT,
        u_p0_over_u_b=WIRSUM_1997_U_P0_OVER_U_B,
        velocity_sigma=WIRSUM_1997_VELOCITY_SIGMA,
        ejection="roof_eq359",
        drag="wirsum_1997_cw",
    ),
    FREEBOARD_MODEL_HAMEL_EXTENDED: FreeboardModelSpec(
        name=FREEBOARD_MODEL_HAMEL_EXTENDED,
        u_p0_over_u_b=MARTENS_U_P0_OVER_U_B,
        velocity_sigma=WIRSUM_1997_VELOCITY_SIGMA,
        ejection="roof_plus_wake",
        drag="hamel_z_exp",
    ),
    FREEBOARD_MODEL_ENGINEERING: FreeboardModelSpec(
        name=FREEBOARD_MODEL_ENGINEERING,
        u_p0_over_u_b=MARTENS_U_P0_OVER_U_B,
        velocity_sigma=WIRSUM_1997_VELOCITY_SIGMA,
        ejection="roof_plus_wake_renormalize",
        drag="haider_mixed",
    ),
}


def resolve_freeboard_model(configured: str) -> FreeboardModelSpec:
    """解析自由板模型名。``analytical_wirsum`` 映射到 engineering，不是 exact。"""
    name = str(configured).strip().lower()
    if name == FREEBOARD_MODEL_LEGACY_ALIAS:
        return _SPECS[FREEBOARD_MODEL_ENGINEERING]
    if name not in _SPECS:
        raise ValueError(
            f"Unknown freeboard_model={configured!r}; expected "
            f"{FREEBOARD_MODEL_WIRSUM_1997_EXACT!r}, "
            f"{FREEBOARD_MODEL_ENGINEERING!r}, "
            f"{FREEBOARD_MODEL_HAMEL_EXTENDED!r}, or legacy "
            f"{FREEBOARD_MODEL_LEGACY_ALIAS!r}"
        )
    return _SPECS[name]


def calc_u_p0_wirsum_1997(u_b: float) -> float:
    """Wirsum §3.4 平均抛射速度 ``u_{p,0} = 1.2 u_b`` [m/s]。

    Source: Wirsum (1997) §3.4
    """
    assert u_b >= 0.0, "u_b 必须非负 [m/s]"
    return WIRSUM_1997_U_P0_OVER_U_B * float(u_b)


def calc_u_p0_martens(u_b: float) -> float:
    """Martens / Hamel Eq. 3.73：``u_{p,0} = 1.53 u_b`` [m/s]，工程/Hamel 扩展。

    Source: Hamel (1999) Eq. 3.73; Martens (1984)
    """
    assert u_b >= 0.0, "u_b 必须非负 [m/s]"
    return MARTENS_U_P0_OVER_U_B * float(u_b)


def build_launch_velocity_samples(
    u_p0: float,
    sigma: float,
    n_bins: int,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """离散高斯竖直抛射速度与质量分数。

    ``u_i = max(u_p0 + σ ξ_i, 0)``，``ξ`` 均匀取 ``[-2, 2]``，
    权重 ``∝ exp(-ξ²/2)`` 并归一化到 1。``σ`` 是绝对标准差 [m/s]。

    Source: Wirsum (1997) §3.4；与现有 ``_build_velocity_samples`` 同构，
            但不把下限抬到 ``1e-6``（那是工程地板）。
    """
    assert n_bins >= 1, "n_bins 必须 >= 1"
    assert sigma >= 0.0, "sigma 必须非负 [m/s]"
    u0 = float(u_p0)
    if n_bins == 1:
        return np.array([max(u0, 0.0)], dtype=np.float64), np.array([1.0], dtype=np.float64)
    xi = np.linspace(-2.0, 2.0, n_bins)
    velocities = np.maximum(u0 + float(sigma) * xi, 0.0)
    weights = np.exp(-0.5 * xi * xi)
    weights = weights / np.sum(weights)
    return velocities.astype(np.float64), weights.astype(np.float64)


def calc_roof_mass_flux_eq359(
    *,
    d_p: float,
    rho_p: float,
    eps_d: float,
    eps_b: float,
    u_b: float,
    d_b: float,
) -> float:
    """Wirsum Eq. 3.59 / Hamel Eq. 3.53 屋顶抛射通量 [kg/(m² s)]。

    ``F_{0,r} = d_p · ρ_p · (1 − ε_d) · (3 ε_b u_b / d_b)``

    Source: Wirsum (1997) Eq. 3.59; Hamel (1999) Eq. 3.53
    """
    assert d_p > 0.0 and rho_p > 0.0 and d_b > 0.0
    assert u_b >= 0.0 and 0.0 <= eps_b < 1.0 and 0.0 <= eps_d < 1.0
    return (
        float(d_p)
        * float(rho_p)
        * max(1.0 - float(eps_d), 0.0)
        * (3.0 * float(eps_b) * float(u_b) / float(d_b))
    )


def calc_wake_mass_flux_engineering(
    *,
    rho_p: float,
    eps_d: float,
    eps_b: float,
    u_b: float,
    zeta_w: float = WAKE_ZETA_W,
    f_w: float = WAKE_F_W,
) -> float:
    """Wake 抛射通量 [kg/(m² s)]。Hamel Eq. 3.54，不属于 ``wirsum_1997_exact``。

    ``F_{0,w} = ζ_w · f_w · ε_b · u_b · ρ_p · (1 − ε_d)``

    Source: Hamel (1999) Eq. 3.54; Pemberton & Davidson (1986)
    """
    assert rho_p > 0.0 and u_b >= 0.0
    return (
        float(zeta_w)
        * float(f_w)
        * float(eps_b)
        * float(u_b)
        * float(rho_p)
        * max(1.0 - float(eps_d), 0.0)
    )


def eject_size_class_mass_flows(
    *,
    d_p_classes: npt.NDArray[np.float64],
    feed_classes: npt.NDArray[np.float64],
    rho_p: float,
    eps_d: float,
    eps_b: float,
    u_b: float,
    d_b: float,
    area: float,
    mode: str,
) -> npt.NDArray[np.float64]:
    """按粒径级给出床面抛射质量流 [kg/s]。

    ``wirsum_1997_exact``：只用 Roof Eq. 3.59，不 renormalize。
    ``hamel_extended``：Roof + Wake，不 renormalize。
    ``engineering_conservative`` / ``engineering_freeboard``：Roof + Wake，
    再把总和缩放到 ``sum(feed_classes)``。

    ``feed_classes`` 只提供组成分数；exact 模式不把它当成必须守住的总量。
    """
    spec_name = str(mode).strip().lower()
    if spec_name == "engineering_conservative":
        spec_name = FREEBOARD_MODEL_ENGINEERING
    if spec_name == FREEBOARD_MODEL_LEGACY_ALIAS:
        spec_name = FREEBOARD_MODEL_ENGINEERING
    spec = resolve_freeboard_model(spec_name)

    d_p = np.maximum(np.asarray(d_p_classes, dtype=np.float64), 1e-9)
    feed = np.maximum(np.asarray(feed_classes, dtype=np.float64), 0.0)
    assert d_p.shape == feed.shape
    assert area > 0.0
    total_feed = float(np.sum(feed))
    if feed.size == 0 or total_feed <= 0.0:
        return np.zeros_like(feed)

    solid_frac = feed / total_feed
    d_p_mean = float(np.sum(d_p * solid_frac))  # [m]
    roof = calc_roof_mass_flux_eq359(
        d_p=d_p_mean,
        rho_p=rho_p,
        eps_d=eps_d,
        eps_b=eps_b,
        u_b=u_b,
        d_b=d_b,
    )
    wake = 0.0
    if spec.ejection != "roof_eq359":
        wake = calc_wake_mass_flux_engineering(
            rho_p=rho_p, eps_d=eps_d, eps_b=eps_b, u_b=u_b
        )
    m_dot = (roof + wake) * area * solid_frac  # [kg/s]
    if spec.ejection == "roof_plus_wake_renormalize":
        total_raw = float(np.sum(m_dot))
        if total_raw > 0.0:
            m_dot = m_dot * (total_feed / total_raw)
    return m_dot.astype(np.float64)


def wirsum_1997_drag_z(z: float | None = None) -> float:
    """高 Re 阻力极限 ``z`` [-]。默认 0.44，论文范围约 0.43–1.0。"""
    value = WIRSUM_1997_Z_DRAG_DEFAULT if z is None else float(z)
    if value <= 0.0:
        raise ValueError(f"Wirsum drag z must be positive, got {value}")
    return value


def wirsum_1997_abc(
    *,
    u_g: float,
    u_p: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float | None = None,
) -> tuple[float, float, float, str]:
    """把 ``c_w = 24/Re + z`` 展开成 ``du_p/dt = A u_p² + B u_p + C``。

    更快 / 更慢分支由 ``u_p ? u_g`` 决定；``Δ = 4AC − B²`` 之后才选解析族。
    ``z`` 是常数，不由球形度关联式计算（那是 Hamel / Haider 扩展）。

    Source: Wirsum (1997) §3.4
    """
    assert d_p > 0.0 and rho_p > 0.0 and rho_g > 0.0 and mu_g > 0.0
    z_drag = wirsum_1997_drag_z(z)
    g_eff = g * (1.0 - float(rho_g) / float(rho_p))  # [m/s²]
    k_z = 0.75 * z_drag * float(rho_g) / (float(rho_p) * float(d_p))  # [1/m]
    k_mu = 18.0 * float(mu_g) / (float(rho_p) * float(d_p) * float(d_p))  # [1/s]
    u_gas = float(u_g)
    if float(u_p) < u_gas:
        branch = "slower"
        a = k_z
        b = -2.0 * k_z * u_gas - k_mu
        c = k_z * u_gas * u_gas + k_mu * u_gas - g_eff
    else:
        branch = "faster"
        a = -k_z
        b = 2.0 * k_z * u_gas - k_mu
        c = -k_z * u_gas * u_gas + k_mu * u_gas - g_eff
    return float(a), float(b), float(c), branch


def wirsum_1997_delta(a: float, b: float, c: float) -> float:
    """``Δ = 4AC − B²``。"""
    return 4.0 * float(a) * float(c) - float(b) * float(b)


def particle_acceleration_cw(
    u_p: float,
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float | None = None,
) -> float:
    """竖直加速度 [m/s²]：``c_w = 24/Re + z``，正方向向上。"""
    z_drag = wirsum_1997_drag_z(z)
    u_rel = float(u_p) - float(u_g)  # [m/s]
    re_p = float(rho_g) * abs(u_rel) * float(d_p) / float(mu_g)
    re_eff = max(re_p, 1e-16)
    c_w = 24.0 / re_eff + z_drag
    g_eff = g * (1.0 - float(rho_g) / float(rho_p))
    drag = 0.75 * c_w * float(rho_g) * abs(u_rel) * u_rel / (float(rho_p) * float(d_p))
    return -g_eff - drag


def quadratic_ode_state(
    u0: float,
    a: float,
    b: float,
    c: float,
    dt: float,
) -> tuple[float, float, float]:
    """解析步进 ``du/dt = A u² + B u + C``，返回 ``(u(t), Δz, Δ)``。

    ``Δz`` 是 ``u`` 在 ``[0, dt]`` 上的精确积分 [m]。

    Source: Wirsum (1997) §3.4；``Δ = 4AC − B²`` 选分支
    """
    delta = wirsum_1997_delta(a, b, c)
    if dt <= 0.0:
        return float(u0), 0.0, float(delta)

    if abs(a) < 1e-14:
        if abs(b) < 1e-14:
            u1 = float(u0 + c * dt)
            z1 = float(u0 * dt + 0.5 * c * dt * dt)
            return u1, z1, float(delta)
        exp_bt = math.exp(float(np.clip(b * dt, -80.0, 80.0)))
        u1 = float(u0 * exp_bt + (c / b) * (exp_bt - 1.0))
        z1 = float((u0 + c / b) * (exp_bt - 1.0) / b - (c / b) * dt)
        return u1, z1, float(delta)

    if delta < -1e-14:
        k = math.sqrt(-delta)
        r1 = (-b - k) / (2.0 * a)
        r2 = (-b + k) / (2.0 * a)
        den0 = u0 - r2
        if abs(den0) < 1e-14:
            den0 = 1e-14 if den0 >= 0.0 else -1e-14
        q0 = (u0 - r1) / den0
        q = q0 * math.exp(float(np.clip(-k * dt, -80.0, 80.0)))
        den = 1.0 - q
        if abs(den) < 1e-14:
            den = 1e-14 if den >= 0.0 else -1e-14
        u1 = float((r1 - q * r2) / den)
        den_q0 = 1.0 - q0
        den_q = 1.0 - q
        if abs(den_q0) < 1e-14:
            den_q0 = 1e-14 if den_q0 >= 0.0 else -1e-14
        if abs(den_q) < 1e-14:
            den_q = 1e-14 if den_q >= 0.0 else -1e-14
        z1 = float(r1 * dt + ((r1 - r2) / k) * math.log(abs(den_q / den_q0)))
        return u1, z1, float(delta)

    if abs(delta) <= 1e-14:
        r = -b / (2.0 * a)
        y0 = u0 - r
        den = 1.0 - a * y0 * dt
        if abs(den) < 1e-14:
            den = 1e-14 if den >= 0.0 else -1e-14
        u1 = float(r + y0 / den)
        z1 = float(r * dt - math.log(abs(den)) / a)
        return u1, z1, float(delta)

    omega = math.sqrt(delta)
    theta0 = math.atan((2.0 * a * u0 + b) / omega)
    theta = theta0 + 0.5 * omega * dt
    half_pi = 0.5 * math.pi
    while theta > half_pi:
        theta -= math.pi
    while theta < -half_pi:
        theta += math.pi
    theta = float(np.clip(theta, -half_pi + 1e-6, half_pi - 1e-6))
    u1 = float((-b + omega * math.tan(theta)) / (2.0 * a))
    z1 = float(
        (-b / (2.0 * a)) * dt
        + (1.0 / a) * math.log(abs(math.cos(theta0) / max(abs(math.cos(theta)), 1e-14)))
    )
    return u1, z1, float(delta)


def integrate_quadratic_analytical(
    u0: float,
    a: float,
    b: float,
    c: float,
    t_end: float,
    n_steps: int,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """固定 ``A,B,C`` 的解析轨迹 ``(t, u, z)``。"""
    assert n_steps >= 1 and t_end >= 0.0
    t = np.linspace(0.0, float(t_end), n_steps + 1)
    u = np.empty_like(t)
    z = np.empty_like(t)
    u[0] = float(u0)
    z[0] = 0.0
    for i in range(n_steps):
        dt = float(t[i + 1] - t[i])
        u[i + 1], dz, _ = quadratic_ode_state(float(u[i]), a, b, c, dt)
        z[i + 1] = z[i] + dz
    return t, u, z


def integrate_quadratic_rk4(
    u0: float,
    a: float,
    b: float,
    c: float,
    t_end: float,
    n_steps: int,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """同一 ``A,B,C`` 的 RK4 对照；高度用梯形积分。"""
    assert n_steps >= 1 and t_end >= 0.0
    dt = float(t_end) / float(n_steps)
    t = np.linspace(0.0, float(t_end), n_steps + 1)
    u = np.empty_like(t)
    z = np.empty_like(t)
    u[0] = float(u0)
    z[0] = 0.0

    def rhs(speed: float) -> float:
        return a * speed * speed + b * speed + c

    for i in range(n_steps):
        ui = float(u[i])
        k1 = rhs(ui)
        k2 = rhs(ui + 0.5 * dt * k1)
        k3 = rhs(ui + 0.5 * dt * k2)
        k4 = rhs(ui + dt * k3)
        u[i + 1] = ui + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        z[i + 1] = z[i] + 0.5 * dt * (ui + float(u[i + 1]))
    return t, u, z


def classify_particle_fate(
    u_p0: float,
    a: float,
    b: float,
    c: float,
    cell_height: float,
    *,
    dt: float = 1.0e-3,
    t_max: float = 5.0,
) -> ParticleFate:
    """固定 ``A,B,C`` 下颗粒到达下一格顶还是在格内回落。"""
    assert cell_height > 0.0 and dt > 0.0
    u = float(u_p0)
    z = 0.0
    t = 0.0
    while t < t_max:
        u_next, dz, _ = quadratic_ode_state(u, a, b, c, dt)
        z_next = z + dz
        if z_next >= cell_height:
            return "reaches_next_cell"
        if u_next <= 0.0 and z_next < cell_height:
            return "returns"
        u = u_next
        z = z_next
        t += dt
    return "undecided"
