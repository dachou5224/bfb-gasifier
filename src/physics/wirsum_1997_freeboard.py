"""Wirsum (1997) §3.3–3.4 自由板内核：模型名、抛射、Roof 3.59、解析轨迹。

本模块是 thesis-exact / engineering 的对照核，不接入生产 ``simulate_freeboard``。
现有 ``analytical_wirsum`` 是工程混合路径的别名，不是 Wirsum 原文。

Source: Wirsum (1997) §3.3–3.4 Eq. 3.59; Hamel (1999) Eq. 3.53（Roof 同构）,
        Eq. 3.73（Martens 1.53 u_b，工程/Hamel 扩展）
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Callable, Literal

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
HAMEL_PHI_S_MEAN_DEFAULT: float = 0.8  # [-] Eq. 3.86 缺省均值，未解析形状分布
HAMEL_PHI_S_MIN: float = 0.3  # [-]
HAMEL_PHI_S_MAX: float = 1.0  # [-]
HAMEL_EQ386_Z_PREFACTOR: float = 73.69  # [-] Eq. 3.84 高 Re 极限
HAMEL_EQ386_Z_EXP: float = 5.0784  # [-]
GAUSSIAN_SAMPLE_MODE: str = "truncated_pm_2sigma_uniform_xi"  # 工程离散，不是 Eq. 3.72 积分算法

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
    权重 ``∝ exp(-ξ²/2)`` 并归一化到 1。``σ`` 是绝对标准差 [m/s]，
    **不是** ``0.6 ū_{p,0}``。正文称 Varianz，公式按标准差用。

    本离散标 ``truncated_pm_2sigma_uniform_xi``，不是 Hamel Eq. 3.72
    的规范密度积分。负速度截零后再归一化。
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


def zeta_w_f_w_hamel_365(u0_over_umf: float) -> float:
    """Hamel 印刷页 40 式 (3.65) 的无量纲乘积 ``ζ_w f_w``。

    ``u_0/u_mf < 10`` 时为 0.1。
    ``u_0/u_mf ≥ 10`` 时为 ``0.5 (1 - exp(9.5 - u_0/u_mf))``。
    第二支的渐近上限是 0.5，文中 Milioli 与 Foster 的估计，不是第三段。
    现行 ``hamel_extended`` 的 0.4×0.25 只与第一支数值相同。

    Source: Hamel (1999) Gl. (3.65)，PDF 第 30 页左栏
    """
    ratio = float(u0_over_umf)
    if ratio < 0.0:
        raise ValueError(f"u0/umf must be >= 0, got {ratio}")
    if ratio < 10.0:
        return 0.1
    exponent = float(np.clip(9.5 - ratio, -100.0, 100.0))
    return float(0.5 * (1.0 - np.exp(exponent)))


def f0_hamel_367(
    *,
    d_bar: float,
    rho_p: float,
    eps_d: float,
    eps_b: float,
    u_b: float,
    d_b: float,
    zeta_f: float,
    S: float,
) -> float:
    """Hamel Eq. 3.67 总抛射通量 [kg/(m²·s)]。

    ``S`` 是式中 ``Σ_{k=2}^{n_B} (k-1) P''_k``。``S=0`` 只剩 roof，``S=1`` 只剩 wake。
    ``zeta_f`` 用 ``zeta_w_f_w_hamel_365``。这里不从库存反推，也不把两项无条件相加。
    """
    if not 0.0 <= float(S) <= 1.0:
        raise ValueError(f"coalescence weight S must be in [0, 1], got {S}")
    roof = calc_roof_mass_flux_eq359(
        d_p=d_bar, rho_p=rho_p, eps_d=eps_d, eps_b=eps_b, u_b=u_b, d_b=d_b,
    )
    wake = calc_wake_mass_flux_engineering(
        rho_p=rho_p, eps_d=eps_d, eps_b=eps_b, u_b=u_b,
        zeta_w=float(zeta_f), f_w=1.0,
    )
    return float((1.0 - float(S)) * roof + float(S) * wake)


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


def hamel_eq386_asymptotic_z(phi_s: float = HAMEL_PHI_S_MEAN_DEFAULT) -> float:
    """Hamel Eq. 3.86 高 Re 极限 ``z(φ)=73.69 exp(-5.0784 φ)`` [-]。

    原页没有唯一 ``Re_∞``；本式是 Eq. 3.84 的严格高 Re 极限，标
    ``hamel_eq386_asymptotic``。缺省 ``φ=0.8`` 是均值，未解析形状分布。
    """
    phi = float(phi_s)
    if not (HAMEL_PHI_S_MIN <= phi <= HAMEL_PHI_S_MAX):
        raise ValueError(
            f"phi_s must be in [{HAMEL_PHI_S_MIN}, {HAMEL_PHI_S_MAX}], got {phi_s}"
        )
    return HAMEL_EQ386_Z_PREFACTOR * math.exp(-HAMEL_EQ386_Z_EXP * phi)


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


def drag_equilibrium_speed_m_s(
    *,
    u_g: float,
    u_p: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
) -> float:
    """阻力平衡 ``du/dt = 0`` 的稳定速度 [m/s]。

    只用来估计还要积多少步，不代替轨迹，也不写成停留时间。
    Source: Wirsum (1997) §3.4
    """
    a, b, c, _branch = wirsum_1997_abc(
        u_g=u_g, u_p=u_p, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
    )
    if abs(float(a)) <= 1.0e-14:
        if abs(float(b)) <= 1.0e-14:
            return float(u_p)
        return float(-c / b)
    disc = float(b) * float(b) - 4.0 * float(a) * float(c)
    if disc < 0.0:
        return float(u_p)
    root = math.sqrt(disc)
    roots = ((-float(b) + root) / (2.0 * float(a)), (-float(b) - root) / (2.0 * float(a)))
    stable = [float(r) for r in roots if 2.0 * float(a) * float(r) + float(b) < 0.0]
    if not stable:
        return float(u_p)
    return float(min(stable, key=lambda r: abs(r - float(u_p))))


def drag_work_steps(remain_m: float, u_m_s: float, u_eq_m_s: float, dt_s: float) -> int:
    """剩余高度除以较慢的正速度，再除以步长。结果是步数，不是停留时间。"""
    speeds = [abs(float(u_m_s))]
    if float(u_eq_m_s) > 1.0e-6:
        speeds.append(abs(float(u_eq_m_s)))
    speed = max(min(speeds), 1.0e-3)  # [m/s]
    if float(remain_m) <= 0.0 or float(dt_s) <= 0.0:
        return 0
    return int(math.ceil(float(remain_m) / speed / float(dt_s)))


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


@dataclass(frozen=True)
class _VelocityPiece:
    """一段常系数 ``du/dt = A u² + B u + C`` 的正向解析式。

    正切段只覆盖到达第一处竖直渐近线之前，不做周期折返。
    双实根段用对数，长时间下指数下溢时速度收到稳定根。
    """

    kind: str
    a: float
    b: float
    c: float
    u0: float
    t_pole: float
    omega: float = 0.0
    theta0: float = 0.0
    k: float = 0.0
    r1: float = 0.0
    r2: float = 0.0
    q0: float = 0.0

    def speed(self, t: float) -> float:
        """这段上的速度 [m/s]。``t`` 不超过 ``t_pole``。"""
        tt = float(t)
        if tt <= 0.0:
            return float(self.u0)
        if self.kind == "tan":
            theta = self.theta0 + 0.5 * self.omega * tt
            return float((-self.b + self.omega * math.tan(theta)) / (2.0 * self.a))
        if self.kind == "roots":
            if self.q0 == 0.0:
                return float(self.r1)
            q = self._q(tt)
            den = 1.0 - q
            if abs(den) < 1.0e-15:
                return float(self.r1)
            return float((self.r1 - q * self.r2) / den)
        if self.kind == "linear":
            if abs(self.b) < 1.0e-14:
                return float(self.u0 + self.c * tt)
            ebt = math.exp(self.b * tt)
            return float(self.u0 * ebt + (self.c / self.b) * (ebt - 1.0))
        y0 = self.u0 - self.r1
        den = 1.0 - self.a * y0 * tt
        if abs(den) < 1.0e-15:
            return float(self.r1)
        return float(self.r1 + y0 / den)

    def dz(self, t: float) -> float:
        """从这段起点积分到 ``t`` 的高度增量 [m]。"""
        tt = float(t)
        if tt <= 0.0:
            return 0.0
        if self.kind == "tan":
            theta = self.theta0 + 0.5 * self.omega * tt
            cosine = math.cos(theta)
            if abs(cosine) < 1.0e-15:
                cosine = 1.0e-15 if cosine >= 0.0 else -1.0e-15
            return float(
                (-self.b / (2.0 * self.a)) * tt
                + (1.0 / self.a) * math.log(abs(math.cos(self.theta0) / cosine))
            )
        if self.kind == "roots":
            if self.q0 == 0.0:
                return float(self.r1 * tt)
            q = self._q(tt)
            alpha = (self.r1 - self.r2) / self.k
            return float(self.r1 * tt + alpha * math.log(abs((1.0 - q) / (1.0 - self.q0))))
        if self.kind == "linear":
            if abs(self.b) < 1.0e-14:
                return float(self.u0 * tt + 0.5 * self.c * tt * tt)
            ebm = math.expm1(self.b * tt)
            return float((self.u0 + self.c / self.b) * ebm / self.b - (self.c / self.b) * tt)
        y0 = self.u0 - self.r1
        den = 1.0 - self.a * y0 * tt
        if abs(den) < 1.0e-15:
            den = 1.0e-15 if den >= 0.0 else -1.0e-15
        return float(self.r1 * tt - math.log(abs(den)) / self.a)

    def time_to_speed(self, u_target: float) -> float | None:
        """第一次到达 ``u_target`` 的时间 [s]。不在这段正向轨道上则为 None。"""
        ut = float(u_target)
        if abs(ut - self.u0) <= 1.0e-15:
            return 0.0
        if self.kind == "tan":
            arg = (2.0 * self.a * ut + self.b) / self.omega
            if not math.isfinite(arg):
                return None
            theta_t = math.atan(arg)
            t = 2.0 * (theta_t - self.theta0) / self.omega
            if t <= 1.0e-15 or t >= self.t_pole - 1.0e-12:
                return None
            return float(t)
        if self.kind == "roots":
            if self.q0 == 0.0:
                return None
            if (ut - self.u0) * (self.r1 - self.u0) < 0.0 and abs(ut - self.r1) > 1.0e-12:
                return None
            den_t = ut - self.r2
            if abs(den_t) < 1.0e-15 or abs(self.q0) < 1.0e-18:
                return None
            qt = (ut - self.r1) / den_t
            if qt == 0.0 or self.q0 * qt <= 0.0:
                return None
            t = math.log(self.q0 / qt) / self.k
            if t <= 1.0e-15 or t >= self.t_pole - 1.0e-12:
                return None
            return float(t)
        if self.kind == "linear":
            if abs(self.b) < 1.0e-14:
                if abs(self.c) < 1.0e-14:
                    return None
                t = (ut - self.u0) / self.c
                return float(t) if t > 1.0e-15 else None
            base = self.u0 + self.c / self.b
            if abs(base) < 1.0e-15:
                return None
            ratio = (ut + self.c / self.b) / base
            if ratio <= 0.0:
                return None
            t = math.log(ratio) / self.b
            return float(t) if t > 1.0e-15 else None
        y0 = self.u0 - self.r1
        if abs(y0) < 1.0e-15 or abs(self.a) < 1.0e-15:
            return None
        den = y0 / (ut - self.r1)
        t = (1.0 - den) / (self.a * y0) if abs(ut - self.r1) > 1.0e-15 else None
        if t is None or t <= 1.0e-15 or t >= self.t_pole - 1.0e-12:
            return None
        return float(t)

    def _q(self, t: float) -> float:
        if self.q0 == 0.0:
            return 0.0
        log_abs = math.log(abs(self.q0)) - self.k * float(t)
        if log_abs < -80.0:
            return 0.0
        if log_abs > 40.0:
            return math.copysign(math.exp(40.0), self.q0)
        return math.copysign(math.exp(log_abs), self.q0)


def _velocity_piece(u0: float, a: float, b: float, c: float) -> _VelocityPiece:
    """由当前速度所在的 ABC 分支生成一段解析式。"""
    aa = float(a)
    bb = float(b)
    cc = float(c)
    u = float(u0)
    delta = wirsum_1997_delta(aa, bb, cc)
    if abs(aa) < 1.0e-14:
        return _VelocityPiece("linear", aa, bb, cc, u, t_pole=math.inf)
    if delta > 1.0e-12:
        omega = math.sqrt(delta)
        theta0 = math.atan((2.0 * aa * u + bb) / omega)
        t_pole = (math.pi - 2.0 * theta0) / omega
        if t_pole <= 1.0e-15:
            t_pole = 1.0e-15
        return _VelocityPiece(
            "tan", aa, bb, cc, u, t_pole=float(t_pole), omega=float(omega), theta0=float(theta0),
        )
    if delta < -1.0e-12:
        k = math.sqrt(-delta)
        r1 = (-bb - k) / (2.0 * aa)
        r2 = (-bb + k) / (2.0 * aa)
        den0 = u - r2
        if abs(u - r1) <= 1.0e-14:
            q0 = 0.0
        else:
            if abs(den0) < 1.0e-15:
                den0 = 1.0e-15 if den0 >= 0.0 else -1.0e-15
            q0 = (u - r1) / den0
        t_pole = math.inf
        if q0 > 1.0:
            t_pole = math.log(q0) / k
        return _VelocityPiece(
            "roots", aa, bb, cc, u, t_pole=float(t_pole), k=float(k), r1=float(r1), r2=float(r2), q0=float(q0),
        )
    r = -bb / (2.0 * aa)
    y0 = u - r
    t_pole = math.inf
    if abs(aa * y0) > 1.0e-15 and aa * y0 > 0.0:
        t_pole = 1.0 / (aa * y0)
    return _VelocityPiece("repeated", aa, bb, cc, u, t_pole=float(t_pole), r1=float(r))


def _abc_piece(
    u: float,
    *,
    side: str,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
) -> _VelocityPiece:
    """``side`` 固定分支。``u=u_g`` 时不能靠速度大小选择。"""
    probe = float(u_g) - 1.0 if side == "slower" else float(u_g) + 1.0
    a, b, c, _branch = wirsum_1997_abc(
        u_g=u_g, u_p=probe, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
    )
    return _velocity_piece(float(u), a, b, c)


def _bisect_height(
    z_of: Callable[[float], float],
    t_lo: float,
    t_hi: float,
    z_target: float,
    *,
    rising: bool,
) -> float:
    lo = float(t_lo)
    hi = float(t_hi)
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        zm = float(z_of(mid))
        if (rising and zm < z_target) or ((not rising) and zm > z_target):
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def _height_on_ray(
    z_of: Callable[[float], float],
    speed_of: Callable[[float], float],
    t0: float,
    z_target: float,
    *,
    rising: bool,
) -> float | None:
    """从 ``t0`` 起沿单调高度找 ``z_target``。到不了就返回 None。"""
    lo = float(t0)
    z_lo = float(z_of(lo))
    if rising and z_lo >= z_target:
        return lo if lo > 1.0e-12 else None
    if (not rising) and z_lo <= z_target:
        return lo if lo > 1.0e-12 else None
    hi = lo + 1.0
    prev = z_lo
    for _ in range(60):
        if hi > 1.0e15:
            return None
        zh = float(z_of(hi))
        if rising and zh >= z_target:
            return _bisect_height(z_of, lo, hi, z_target, rising=True)
        if (not rising) and zh <= z_target:
            return _bisect_height(z_of, lo, hi, z_target, rising=False)
        if abs(zh - prev) < 1.0e-12 and abs(float(speed_of(hi))) < 1.0e-9:
            return None
        prev = zh
        lo = hi
        hi = t0 + (hi - t0) * 2.0
    return None


def constant_property_boundary_event(
    u_p0: float,
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
    cell_height: float,
    height0: float = 0.0,
    t_cap: float | None = None,
) -> tuple[ParticleFate, float, float, float, str]:
    """常物性单元上的首次边界事件。

    返回 ``(命运, t [s], z [m], u [m/s], 原因)``。
    ``t_cap`` 只截断报告；截断不是物理边界。平衡速度为零且高度停在格内时，
    原因是 ``asymptotic_residence``，不把无穷时间写成停留时间。

    Source: Hamel (1999) Eq. 3.82–3.92；印刷页 49，先定平衡区边界上的速度
    """
    assert cell_height > 0.0
    ug = float(u_g)
    height = float(height0)
    u = float(u_p0)
    props = dict(u_g=ug, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z)
    if height <= 0.0 and u <= 0.0:
        return "returns", 0.0, 0.0, u, "bed_return"
    side = "slower" if u < ug else "faster"
    if abs(u - ug) <= 1.0e-15:
        side = "slower"
    t = 0.0
    for _switch in range(4):
        piece = _abc_piece(u, side=side, **props)
        t_turn = piece.time_to_speed(0.0)
        t_exit = piece.time_to_speed(ug) if abs(u - ug) > 1.0e-12 else None
        if t_exit is not None and t_exit <= 1.0e-12:
            t_exit = None
        knots = [0.0]
        if t_turn is not None and 1.0e-12 < t_turn < piece.t_pole:
            knots.append(float(t_turn))
        if t_exit is not None and t_exit < piece.t_pole:
            knots.append(float(t_exit))
        knots = sorted(set(knots))
        spans: list[tuple[float, float | None]] = []
        for left, right in zip(knots, knots[1:]):
            spans.append((left, right))
        if t_exit is None:
            spans.append((knots[-1], None))

        def z_abs(tt: float, base: float = height, item: _VelocityPiece = piece) -> float:
            return base + item.dz(tt)

        hit_t: float | None = None
        hit_kind = ""
        for t0, t1 in spans:
            rising = piece.speed(t0 + 1.0e-8) > 0.0
            if t1 is None:
                target = float(cell_height) if rising else 0.0
                th = _height_on_ray(z_abs, piece.speed, t0, target, rising=rising)
            else:
                z0 = z_abs(t0)
                z1 = z_abs(t1)
                th = None
                if rising and z0 < cell_height <= z1:
                    th = _bisect_height(z_abs, t0, t1, float(cell_height), rising=True)
                elif (not rising) and z1 <= 0.0 < z0:
                    th = _bisect_height(z_abs, t0, t1, 0.0, rising=False)
            if th is None or th <= 1.0e-12:
                continue
            if hit_t is None or th < hit_t:
                hit_t = float(th)
                hit_kind = "upper_boundary" if rising else "bed_return"
        if hit_t is not None:
            t_event = t + hit_t
            if t_cap is not None and t_event > float(t_cap):
                break
            z_end = float(cell_height) if hit_kind == "upper_boundary" else 0.0
            fate: ParticleFate = "reaches_next_cell" if hit_kind == "upper_boundary" else "returns"
            return fate, float(t_event), z_end, float(piece.speed(hit_t)), hit_kind
        if t_exit is None:
            if piece.kind == "roots" and abs(piece.r1) <= 1.0e-9:
                t_hold = max(knots[-1], 20.0 / piece.k if piece.k > 0.0 else 0.0)
                if math.isfinite(piece.t_pole):
                    t_hold = min(t_hold, 0.5 * piece.t_pole)
                z_hold = height + piece.dz(t_hold)
                return (
                    "undecided",
                    float(t + t_hold),
                    float(min(max(z_hold, 0.0), cell_height)),
                    float(piece.speed(t_hold)),
                    "asymptotic_residence",
                )
            return "undecided", float(t), float(height), float(u), "unresolved_event"
        if t_cap is not None and t + float(t_exit) > float(t_cap):
            break
        height = z_abs(float(t_exit))
        u = ug
        t += float(t_exit)
        side = "slower" if side == "faster" else "faster"
    if t_cap is None:
        return "undecided", float(t), float(height), float(u), "switch_limit"
    piece = _abc_piece(u, side=side, **props)
    remain = max(float(t_cap) - t, 0.0)
    if remain >= piece.t_pole:
        remain = max(piece.t_pole * 0.5, 0.0)
    return (
        "undecided",
        float(t_cap),
        float(height + piece.dz(remain)),
        float(piece.speed(remain)),
        "window",
    )


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
    """固定 ``A,B,C`` 对照：``u≤0`` 即记回落。主线命运用 ``classify_particle_fate_branching``。"""
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


def _step_branching(
    u: float,
    dt: float,
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
) -> tuple[float, float]:
    """跨 ``u_p=u_g`` 换分支后积分 ``(u, Δz)``。"""
    kwargs = dict(u_g=u_g, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z)
    a, b, c, _ = wirsum_1997_abc(u_p=u, **kwargs)
    u1, dz, _ = quadratic_ode_state(u, a, b, c, dt)
    if (u - u_g) * (u1 - u_g) >= 0.0:
        return u1, dz
    lo, hi = 0.0, float(dt)
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        um, _, _ = quadratic_ode_state(u, a, b, c, mid)
        if (u - u_g) * (um - u_g) > 0.0:
            lo = mid
        else:
            hi = mid
    _u_c, dz1, _ = quadratic_ode_state(u, a, b, c, hi)
    side = -1.0 if u > u_g else 1.0
    u_c = float(u_g) + side * 1.0e-12
    remain = max(float(dt) - hi, 0.0)
    a2, b2, c2, _ = wirsum_1997_abc(u_p=u_c, **kwargs)
    u2, dz2, _ = quadratic_ode_state(u_c, a2, b2, c2, remain)
    return u2, dz1 + dz2


def classify_particle_fate_branching(
    u_p0: float,
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
    cell_height: float,
    dt: float = 1.0e-3,
    t_max: float = 5.0,
) -> ParticleFate:
    """回床事件是轨迹回到床面 ``z≤0``。``u_p=0`` 只是转折点。超时 ``undecided``。

    气速本段冻结；跨越 ``u_p=u_g`` 必须换 ABC 分支。
    """
    assert cell_height > 0.0 and dt > 0.0
    u = float(u_p0)
    height = 0.0
    t = 0.0
    left_bed = False
    if u <= 0.0:
        return "returns"
    while t < t_max:
        u, dz = _step_branching(
            u, dt, u_g=u_g, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z
        )
        height += dz
        if height >= cell_height:
            return "reaches_next_cell"
        if height > 1.0e-12:
            left_bed = True
        if left_bed and height <= 0.0:
            return "returns"
        t += dt
    return "undecided"


def integrate_particle_fate_with_time(
    u_p0: float,
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
    cell_height: float,
    dt: float = 1.0e-3,
    t_max: float = 5.0,
    u_g_of_h: Callable[[float], float] | None = None,
    height0: float = 0.0,
    resolve_events: bool = False,
) -> tuple[ParticleFate, float, float, float]:
    """轨迹积分，返回 ``(fate, t_event [s], z_end [m], u_end [m/s])``。

    坐标原点在这一格下边界。``height0`` 允许从格顶进入。
    气速不随高度变化时，用常物性解析边界事件。``resolve_events`` 为假时，
    ``t_max`` 仍是报告窗口：窗口末的状态和返回时间属于同一时刻。
    气速随高度变化时仍逐步积分，最后一步不超过剩余窗口，越过边界的那一步
    再对时间求根。
    Source: Hamel (1999) Eq. 3.82–3.92
    """
    assert cell_height > 0.0 and dt > 0.0
    if u_g_of_h is None:
        fate, t_ev, z_end, u_end, _reason = constant_property_boundary_event(
            float(u_p0),
            u_g=u_g,
            rho_g=rho_g,
            rho_p=rho_p,
            d_p=d_p,
            mu_g=mu_g,
            z=z,
            cell_height=cell_height,
            height0=height0,
            t_cap=None if resolve_events else float(t_max),
        )
        return fate, float(t_ev), float(z_end), float(u_end)
    u = float(u_p0)
    height = float(height0)
    t = 0.0
    left_bed = height > 1.0e-12
    if height <= 0.0 and u <= 0.0:
        return "returns", 0.0, 0.0, u

    def _crossed(dt_use: float) -> tuple[float, float, float]:
        u_gas = float(u_g_of_h(height))
        u_new, dz = _step_branching(
            u, dt_use, u_g=u_gas, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
        )
        return float(u_new), float(height + dz), float(dt_use)

    while t < float(t_max) - 1.0e-15:
        dt_use = min(float(dt), float(t_max) - t)
        _u_new, height_new, _dt_done = _crossed(dt_use)
        top = height_new >= cell_height
        bed = (left_bed or height > 1.0e-12) and height_new <= 0.0
        if top or bed:
            lo = 0.0
            hi = dt_use
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                _um, zm, _ = _crossed(mid)
                if top:
                    if zm < cell_height:
                        lo = mid
                    else:
                        hi = mid
                elif zm > 0.0:
                    lo = mid
                else:
                    hi = mid
            u_hit, z_hit, _ = _crossed(hi)
            t_hit = t + hi
            if top:
                return "reaches_next_cell", float(t_hit), float(z_hit), float(u_hit)
            return "returns", float(t_hit), float(z_hit), float(u_hit)
        u = _u_new
        height = height_new
        t += dt_use
        if height > 1.0e-12:
            left_bed = True
    return "undecided", float(t), float(height), float(u)


def classify_refined_pair(
    coarse: tuple[str, float, float, float],
    fine: tuple[str, float, float, float],
    *,
    tol_t: float = 0.02,
    tol_t_abs_s: float = 0.0,
    tol_z_m: float = 1.0e-3,
    tol_u_m_s: float = 1.0e-3,
) -> str:
    """区分窗口耗尽和步长失配。两者都未到事件且末态接近时，不是 dt 失配。

    短窗上的 2% 可以小于一个积分步。绝对容差用来放过这种步长本身的差别。
    """
    fate_c, t_c, z_c, u_c = coarse
    fate_f, t_f, z_f, u_f = fine
    if fate_c == fate_f == "undecided":
        close = abs(float(z_c) - float(z_f)) <= float(tol_z_m) and abs(float(u_c) - float(u_f)) <= float(tol_u_m_s)
        return "both_window_exhausted" if close else "state_disagreement"
    if fate_c != fate_f:
        return "fate_disagreement"
    scale = max(abs(float(t_c)), abs(float(t_f)), 1.0e-6)
    allowed = max(float(tol_t) * scale, float(tol_t_abs_s))
    if abs(float(t_c) - float(t_f)) > allowed:
        return "event_time_disagreement"
    return "accepted"


def accept_refined_event(
    fate_coarse: str,
    t_coarse_s: float,
    fate_fine: str,
    t_fine_s: float,
    *,
    tol: float = 0.02,
) -> tuple[str, float] | None:
    """两次步长给出同一事件，且时刻相对差不超过 ``tol``，才接受。

    截断时间不是停留时间。步长不一致时保持未判定。
    """
    if fate_coarse != fate_fine or fate_coarse == "undecided":
        return None
    scale = max(abs(float(t_coarse_s)), abs(float(t_fine_s)), 1.0e-6)
    if abs(float(t_coarse_s) - float(t_fine_s)) > float(tol) * scale:
        return None
    return str(fate_fine), float(t_fine_s)


def continue_undecided_event(
    u_p: float,
    height: float,
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
    cell_height: float,
    dt: float,
    extra_s: float,
    u_g_of_h: Callable[[float], float] | None = None,
    tol: float = 0.02,
    pair: dict | None = None,
) -> tuple[ParticleFate, float, float, float, bool, str]:
    """从当前高度和速度再积分 ``extra_s``。

    返回 ``(命运, 附加时间 [s], 格内高度 [m], 速度 [m/s], 是否接受)``。
    附加时间是事件时刻，不是窗口长度。步长不一致时高度和速度仍是细步终态。
    """
    kwargs = dict(
        u_g=u_g, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
        cell_height=cell_height, t_max=float(extra_s), u_g_of_h=u_g_of_h, height0=height,
    )
    coarse = integrate_particle_fate_with_time(u_p, dt=float(dt), **kwargs)
    fine = integrate_particle_fate_with_time(u_p, dt=0.5 * float(dt), **kwargs)
    kind = classify_refined_pair(
        (coarse[0], coarse[1], coarse[2], coarse[3]),
        (fine[0], fine[1], fine[2], fine[3]),
        tol_t=tol,
        tol_t_abs_s=4.0 * float(dt),
    )
    if pair is not None:
        pair["coarse"] = {
            "fate": coarse[0], "t_s": float(coarse[1]), "height_m": float(coarse[2]), "u_m_s": float(coarse[3]),
        }
        pair["fine"] = {
            "fate": fine[0], "t_s": float(fine[1]), "height_m": float(fine[2]), "u_m_s": float(fine[3]),
        }
        pair["kind"] = kind
    if kind != "accepted":
        return "undecided", float(fine[1]), float(fine[2]), float(fine[3]), False, kind
    return fine[0], float(fine[1]), float(fine[2]), float(fine[3]), True, kind


def integrate_particle_through_cells(
    u_p0: float,
    cell_heights: list[float],
    *,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    z: float,
    dt: float = 1.0e-3,
    t_max: float = 5.0,
    u_g_of_h: Callable[[float], float] | None = None,
    continue_budget_s: float = 0.0,
    trace: list | None = None,
    extra_physical_s: float = 0.0,
    max_work_steps: int = 0,
    wall_s: float = 0.0,
    runtime_cap_s: float | None = None,
    resolve_events: bool = False,
) -> tuple[ParticleFate, float, int]:
    """沿多格积分。回床只发生在第 0 格下边界，出炉只发生在最顶格上边界。

    上一格回落时从下一格顶进入，不把该格的向下流直接加进床。
    气速函数的自变量是从床面算起的高度 [m]。
    ``extra_physical_s`` 是窗口之后还能飞的物理时间 [s]。
    ``max_work_steps`` 是粗步数上限，不随把秒数除以 dt 自动放大。
    ``wall_s`` 是这段续积分的墙钟 [s]。三者分开记录。都为 0 时不再续。
    ``runtime_cap_s`` 只是 ``extra_physical_s`` 的旧名。
    """
    assert cell_heights and all(h > 0.0 for h in cell_heights)
    bases = [0.0]
    for h in cell_heights[:-1]:
        bases.append(bases[-1] + float(h))
    i = 0
    height = 0.0
    u = float(u_p0)
    t = 0.0
    steps = 0
    if runtime_cap_s is not None:
        extra_physical_s = float(runtime_cap_s)
    window_left = float(t_max)
    extra_left = float(continue_budget_s)
    note: dict = {
        "pair": None,
        "work_steps": None,
        "used_physical_s": 0.0,
        "used_work_steps": 0,
        "wall_elapsed_s": 0.0,
    }

    def _stop(fate: str, reason: str) -> tuple[ParticleFate, float, int]:
        if trace is not None:
            row = {
                "fate": fate,
                "reason": reason,
                "t_s": float(t),
                "cell": int(i),
                "height_m": float(height),
                "u_m_s": float(u),
            }
            if note["pair"] is not None:
                row["coarse"] = note["pair"]["coarse"]
                row["fine"] = note["pair"]["fine"]
            if note["work_steps"] is not None:
                row["drag_work_steps"] = int(note["work_steps"])
            row["extra_physical_s"] = float(extra_physical_s)
            row["max_work_steps"] = int(max_work_steps)
            row["wall_s"] = float(wall_s)
            row["used_physical_s"] = float(note["used_physical_s"])
            row["used_work_steps"] = int(note["used_work_steps"])
            row["wall_elapsed_s"] = float(note["wall_elapsed_s"])
            trace.append(row)
        return fate, float(t), i

    def _extend_rising() -> str | None:
        """窗口用尽后，步长一致且仍上升时按运行预算再积分。

        返回事件命运；没有事件时返回 None，停止原因在 ``note['reason']``。
        """
        nonlocal t, height, u
        budgets_off = (
            float(extra_physical_s) <= 0.0
            and int(max_work_steps) <= 0
            and float(wall_s) <= 0.0
        )
        if budgets_off or float(u) <= 1.0e-6:
            note["reason"] = "both_window_exhausted"
            return None
        t_wall0 = time.monotonic()
        used_physical = 0.0
        used_work = 0
        while float(u) > 1.0e-6:
            u_gas = float(u_g if gas is None else gas(float(height)))
            u_eq = drag_equilibrium_speed_m_s(
                u_g=u_gas, u_p=float(u), rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
            )
            remain = max(float(cell_heights[i]) - float(height), 0.0)
            work = drag_work_steps(remain, float(u), u_eq, float(dt))
            note["work_steps"] = int(work)
            if float(wall_s) > 0.0 and time.monotonic() - t_wall0 >= float(wall_s):
                note["reason"] = "wall_clock_incomplete"
                break
            if float(extra_physical_s) > 0.0 and used_physical >= float(extra_physical_s) - 1.0e-15:
                note["reason"] = "physical_budget_incomplete"
                break
            if int(max_work_steps) > 0 and used_work >= int(max_work_steps):
                note["reason"] = "work_budget_incomplete"
                break
            remain_phys = (
                float(extra_physical_s) - used_physical if float(extra_physical_s) > 0.0 else 1.0e9
            )
            remain_work = (
                int(max_work_steps) - used_work if int(max_work_steps) > 0 else 10**9
            )
            slice_steps = min(200, int(remain_work), max(1, work if work > 0 else int(remain_work)))
            slice_s = min(float(slice_steps) * float(dt), remain_phys)
            if slice_s <= 1.0e-15:
                note["reason"] = "physical_budget_incomplete"
                break
            pair: dict = {}
            t_before = float(t)
            fate2, t_add, z2, u2, ok, kind2 = continue_undecided_event(
                float(u), float(height),
                u_g=u_g, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
                cell_height=float(cell_heights[i]), dt=dt,
                extra_s=float(slice_s), u_g_of_h=gas, pair=pair,
            )
            used_physical += float(t_add)
            used_work += max(1, int(math.ceil(float(t_add) / float(dt))))
            note["used_physical_s"] = used_physical
            note["used_work_steps"] = used_work
            note["wall_elapsed_s"] = time.monotonic() - t_wall0
            t += float(t_add)
            height = float(z2)
            u = float(u2)
            if pair:
                for side in ("coarse", "fine"):
                    pair[side]["t_s"] = t_before + float(pair[side]["t_s"])
                note["pair"] = pair
            if ok:
                return str(fate2)
            if kind2 != "both_window_exhausted":
                note["reason"] = kind2
                return None
            if float(extra_physical_s) > 0.0 and used_physical >= float(extra_physical_s) - 1.0e-15:
                note["reason"] = "physical_budget_incomplete"
                break
            if int(max_work_steps) > 0 and used_work >= int(max_work_steps):
                note["reason"] = "work_budget_incomplete"
                break
            if float(u) <= 1.0e-6:
                note["reason"] = "both_window_exhausted"
                break
        else:
            note["reason"] = "both_window_exhausted"
        return None

    if u <= 0.0:
        return _stop("returns", "nonpositive_launch")
    if resolve_events and u_g_of_h is None:
        visits = 0
        while visits < 8:
            visits += 1
            fate_e, dt_ev, z_end, u_end, reason = constant_property_boundary_event(
                float(u),
                u_g=u_g,
                rho_g=rho_g,
                rho_p=rho_p,
                d_p=d_p,
                mu_g=mu_g,
                z=z,
                cell_height=float(cell_heights[i]),
                height0=float(height),
                t_cap=None,
            )
            t += float(dt_ev)
            u = float(u_end)
            height = float(z_end)
            if fate_e == "undecided":
                return _stop("undecided", str(reason))
            if fate_e == "reaches_next_cell":
                if i + 1 >= len(cell_heights):
                    height = float(cell_heights[i])
                    return _stop("reaches_next_cell", "furnace_exit")
                i += 1
                height = 0.0
                continue
            if fate_e == "returns" and i == 0:
                height = 0.0
                return _stop("returns", "bed_return")
            if fate_e == "returns":
                i -= 1
                height = float(cell_heights[i])
                continue
            return _stop("undecided", "unclassified")
        return _stop("undecided", "switch_limit")
    while steps < 100000:
        h_abs = bases[i] + max(height, 0.0)
        gas = u_g_of_h if u_g_of_h is None else (lambda h, b=bases[i]: u_g_of_h(b + max(h, 0.0)))
        steps += 1
        extended = False
        if window_left > 1.0e-15:
            slice_s = window_left
            from_window = True
        elif extra_left > 1.0e-15:
            slice_s = extra_left
            from_window = False
        else:
            event = _extend_rising()
            if event is None:
                return _stop("undecided", str(note.get("reason", "window_exhausted")))
            fate = event
            z_end = float(height)
            u_end = float(u)
            extended = True
        if not extended:
            h0 = float(height)
            u0 = float(u)
            fate, dt_ev, z_end, u_end = integrate_particle_fate_with_time(
                u0,
                u_g=u_g,
                rho_g=rho_g,
                rho_p=rho_p,
                d_p=d_p,
                mu_g=mu_g,
                z=z,
                cell_height=float(cell_heights[i]),
                dt=dt,
                t_max=float(slice_s),
                u_g_of_h=gas,
                height0=h0,
            )
            if from_window:
                window_left = max(0.0, window_left - float(dt_ev))
            else:
                extra_left = max(0.0, extra_left - float(dt_ev))
            t += float(dt_ev)
            u = float(u_end)
            height = float(z_end)
            if fate == "undecided":
                if not from_window or extra_left <= 1.0e-15:
                    height = float(z_end)
                    event = _extend_rising()
                    if event is None:
                        return _stop("undecided", str(note.get("reason", "both_window_exhausted")))
                    fate = event
                    z_end = float(height)
                else:
                    pair: dict = {}
                    fate2, t_add, z2, u2, ok, kind = continue_undecided_event(
                        u, float(z_end),
                        u_g=u_g, rho_g=rho_g, rho_p=rho_p, d_p=d_p, mu_g=mu_g, z=z,
                        cell_height=float(cell_heights[i]), dt=dt,
                        extra_s=float(extra_left), u_g_of_h=gas, pair=pair,
                    )
                    spent = float(t_add)
                    extra_left = max(0.0, extra_left - spent)
                    t_before = float(t)
                    t += spent
                    height = float(z2)
                    u = float(u2)
                    if pair:
                        for side in ("coarse", "fine"):
                            pair[side]["t_s"] = t_before + float(pair[side]["t_s"])
                        note["pair"] = pair
                    if not ok:
                        if kind == "both_window_exhausted":
                            event = _extend_rising()
                            if event is None:
                                return _stop("undecided", str(note.get("reason", kind)))
                            fate = event
                            z_end = float(height)
                        else:
                            return _stop("undecided", kind)
                    else:
                        z_end = float(z2)
                        fate = fate2
        if fate == "reaches_next_cell":
            over = max(float(z_end) - float(cell_heights[i]), 0.0)
            if i + 1 >= len(cell_heights):
                return _stop("reaches_next_cell", "furnace_exit")
            i += 1
            height = over
            continue
        if fate == "returns" and i == 0:
            return _stop("returns", "bed_return")
        if fate == "returns":
            i -= 1
            height = float(cell_heights[i])
            continue
        height = float(z_end)
        return _stop("undecided", "unclassified")
    height = float(height)
    return _stop("undecided", "step_cap")


def u_gb_hamel_374(
    h: float,
    u_superficial: float,
    u_bubble: float,
    d_bubble: float,
) -> float:
    """床面气泡破裂后的气体速度 [m/s]。Hamel (1999) Eq. 3.74–3.76。

    ``u_0`` 取表观气速，``u_gb,0`` 取床面气泡速度，``β_A = 1/(3 d_b)``。
    ``h=0`` 时等于气泡速度，高度增大后趋近表观气速。
    """
    beta = 1.0 / max(3.0 * float(d_bubble), 1.0e-12)  # [1/m]
    decay = np.exp(-beta * max(float(h), 0.0))
    return float(u_superficial) + (float(u_bubble) - float(u_superficial)) * float(decay)


def drag_z_for_spec(spec: FreeboardModelSpec, *, phi_s: float = HAMEL_PHI_S_MEAN_DEFAULT) -> float:
    """按模型规格取 z。``hamel_z_exp`` 走 Eq. 3.86 极限，不是 Wirsum 默认 0.44。"""
    if spec.drag == "hamel_z_exp":
        return hamel_eq386_asymptotic_z(phi_s)
    return wirsum_1997_drag_z()
