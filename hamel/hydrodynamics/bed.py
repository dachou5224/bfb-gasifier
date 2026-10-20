"""床区流动（§3.1）与床内固体输运系数。

采用式：α_b = u_b·ε_d/u_d；式 (3.43) 的压力系数为 ψ_b·(P/P_0)^0,2，P_0 = 101300 Pa。
气泡直径用定步长 RK4。式 (3.24) 的分母在加压时用式 (3.44) 的 u_b,r，常压与印刷式相同。
g = 9,81 m/s²，与 ADR 0014 的床料质量换算相同。Geldart 类别由调用方传入，原文没有印出分类边界。
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
G = 9.81
P_0 = 101_300.0
N_B = 2.7
F_W = 0.25
_FINE_PARTICLE_M = 0.2e-3
class HydrodynamicsError(ValueError):
    """适用条件不满足，或沿床高的积分不能在原方程下继续。"""
def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or isinstance(value, np.ndarray):
        raise HydrodynamicsError(f"{name} 必须是有限标量")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise HydrodynamicsError(f"{name} 必须是有限标量") from exc
    if not np.isfinite(number):
        raise HydrodynamicsError(f"{name} 不是有限值")
    return number
def _positive(name: str, value: object) -> float:
    number = _finite(name, value)
    if number <= 0.0:
        raise HydrodynamicsError(f"{name} 必须大于 0")
    return number
def archimedes_number(diameter_m, nu_g, rho_p, rho_g) -> float:
    """Ar = g d³ (ρ_P − ρ_g) / (ν_g² ρ_g)。
    Source: Hamel Eq.3.2, S. 21
    Source: Hamel Eq.3.32, S. 28
    """
    diameter = _positive("diameter_m", diameter_m)
    kinematic = _positive("nu_g", nu_g)
    density_p = _finite("rho_p", rho_p)
    density_g = _positive("rho_g", rho_g)
    if density_p <= density_g:
        raise HydrodynamicsError("ρ_P 必须大于 ρ_g")
    return G * diameter**3 * (density_p - density_g) / (kinematic**2 * density_g)
def ergun_archimedes(epsilon_mf, phi_s, reynolds_mf) -> float:
    """由 Re_mf 按 Ergun 式写出 Ar。
    Source: Hamel Eq.3.1, S. 21
    """
    voidage = _positive("epsilon_mf", epsilon_mf)
    sphericity = _positive("phi_s", phi_s)
    reynolds = _finite("reynolds_mf", reynolds_mf)
    if not 0.0 < voidage < 1.0:
        raise HydrodynamicsError("ε_mf 必须在 (0, 1) 内")
    viscous = 150.0 * (1.0 - voidage) / (sphericity**2 * voidage**3) * reynolds
    inertial = 1.75 / (sphericity * voidage**3) * reynolds**2
    return viscous + inertial
def sphericity_diameter_product(u_mf_ambient, epsilon_mf, mu_g, rho_g, rho_p) -> float:
    """环境条件下由测得的 u_mf,0 反求 φ_S·d̄_P。物性是环境条件的物性。
    Source: Hamel Eq.3.3, S. 22
    """
    velocity = _positive("u_mf_ambient", u_mf_ambient)
    voidage = _positive("epsilon_mf", epsilon_mf)
    viscosity = _positive("mu_g", mu_g)
    density_g = _positive("rho_g", rho_g)
    density_p = _finite("rho_p", rho_p)
    if not 0.0 < voidage < 1.0:
        raise HydrodynamicsError("ε_mf 必须在 (0, 1) 内")
    if density_p <= density_g:
        raise HydrodynamicsError("ρ_P 必须大于 ρ_g")
    buoyancy = voidage**3 * (density_p - density_g) * G
    linear = 1.75 * density_g * velocity**2 / (2.0 * buoyancy)
    under = 150.0 * (1.0 - voidage) * viscosity * velocity / buoyancy
    return linear + (under + linear**2) ** 0.5
def minimum_fluidization_velocity(epsilon_mf, nu_g, phi_s, diameter_m, rho_p, rho_g) -> float:
    """操作条件下的最小流化速度。物性是操作条件的物性。
    Source: Hamel Eq.3.4, S. 22
    """
    voidage = _positive("epsilon_mf", epsilon_mf)
    kinematic = _positive("nu_g", nu_g)
    sphericity = _positive("phi_s", phi_s)
    diameter = _positive("diameter_m", diameter_m)
    density_p = _finite("rho_p", rho_p)
    density_g = _positive("rho_g", rho_g)
    if not 0.0 < voidage < 1.0:
        raise HydrodynamicsError("ε_mf 必须在 (0, 1) 内")
    if density_p <= density_g:
        raise HydrodynamicsError("ρ_P 必须大于 ρ_g")
    specific = sphericity * diameter
    inside = (
        3.1e-4
        * voidage**3
        * G
        * specific**3
        * (density_p - density_g)
        / ((1.0 - voidage) ** 2 * kinematic**2 * density_g)
    )
    if inside <= -1.0:
        raise HydrodynamicsError("式 (3.4) 的根号内不是正数")
    return 42.9 * (1.0 - voidage) * kinematic / specific * ((1.0 + inside) ** 0.5 - 1.0)
def suspension_superficial_velocity(u_mf, u_0, distributor, geldart, three_dimensional) -> float:
    """Hilligardt 的悬浮相空塔速度。只接受技术布风器、三维床、Geldart B 或 D。
    Source: Hamel Eq.3.11, S. 25
    """
    minimum = _positive("u_mf", u_mf)
    superficial = _finite("u_0", u_0)
    _require_technical_three_dimensional(distributor, three_dimensional)
    if geldart not in ("B", "D"):
        raise HydrodynamicsError("式 (3.11) 只对 Geldart B 或 D 有采用式")
    if superficial <= minimum:
        raise HydrodynamicsError("u_0 必须大于 u_mf")
    return minimum + (superficial - minimum) / 3.0
def throughflow_velocity(u_d, pressure) -> float:
    """u_b,r = 2,7 u_d (P/P_0)^(−0,15)。P = P_0 时就是式 (3.23)。
    Source: Hamel Eq.3.44, S. 32
    Source: Hamel Eq.3.23, S. 27
    """
    suspension = _positive("u_d", u_d)
    absolute = _positive("pressure", pressure)
    return N_B * suspension * (absolute / P_0) ** (-0.15)
def visible_bubble_fraction(u_0, u_d, u_b, u_br) -> float:
    """由式 (3.7)–(3.10) 解出的可见气泡份额。
    分母是 u_b + u_b,r − u_d。常压下 u_b,r = n_b u_d，与印刷式 (3.24) 相同。
    Source: Hamel Eq.3.24, S. 27
    Source: Hamel Eq.3.7, S. 24
    """
    superficial = _finite("u_0", u_0)
    suspension = _finite("u_d", u_d)
    bubble = _finite("u_b", u_b)
    through = _finite("u_br", u_br)
    denominator = bubble + through - suspension
    if denominator == 0.0 or not np.isfinite(denominator):
        raise HydrodynamicsError("气泡份额的分母为零")
    return (superficial - suspension) / denominator
def richardson_zaki_exponent(reynolds_s) -> float:
    """Hilligardt 采用的 Richardson–Zaki 指数。
    Source: Hamel Eq.3.26, S. 27
    """
    reynolds = _finite("reynolds_s", reynolds_s)
    if reynolds <= 0.0:
        raise HydrodynamicsError("Re_S 必须大于 0")
    if reynolds <= 0.2:
        return 4.65
    if reynolds <= 1.0:
        return 4.4 * reynolds ** (-0.03)
    if reynolds <= 500.0:
        return 4.4 * reynolds ** (-0.1)
    return 2.4
def suspension_voidage(epsilon_mf, u_d, u_mf, n_rz) -> float:
    """ε_d = ε_mf (u_d/u_mf)^(1/n_RZ)。
    Source: Hamel Eq.3.25, S. 27
    """
    voidage = _positive("epsilon_mf", epsilon_mf)
    suspension = _positive("u_d", u_d)
    minimum = _positive("u_mf", u_mf)
    exponent = _positive("n_rz", n_rz)
    if not 0.0 < voidage < 1.0:
        raise HydrodynamicsError("ε_mf 必须在 (0, 1) 内")
    return voidage * (suspension / minimum) ** (1.0 / exponent)
def terminal_velocity(diameter_m, nu_g, rho_p, rho_g) -> float:
    """单颗粒沉降速度。粒径是所算粒级的 d_P。
    Source: Hamel Eq.3.31, S. 28
    """
    diameter = _positive("diameter_m", diameter_m)
    kinematic = _positive("nu_g", nu_g)
    archimedes = archimedes_number(diameter, kinematic, rho_p, rho_g)
    factor = (1.0 + (1.0 / 9.0) * archimedes**0.5) ** 0.5 - 1.0
    return 18.0 * (kinematic / diameter) * factor**2
def settling_reynolds(terminal_m_s, diameter_m, rho_g, mu_g) -> float:
    """Re = u_S d_P ρ_g / μ_g。式 (3.31) 由它与式 (3.30) 解出，故即 Re_S。
    Source: Hamel Eq.3.28, S. 28
    """
    velocity = _positive("terminal_m_s", terminal_m_s)
    diameter = _positive("diameter_m", diameter_m)
    density = _positive("rho_g", rho_g)
    viscosity = _positive("mu_g", mu_g)
    return velocity * diameter * density / viscosity
def suspension_dynamic_viscosity(diameter_m) -> float:
    """d_P < 0,2 mm 用式 (3.21)，其余用式 (3.22)。粒径以 m 代入。
    Source: Hamel Eq.3.21, S. 27
    Source: Hamel Eq.3.22, S. 27
    """
    diameter = _positive("diameter_m", diameter_m)
    if diameter < _FINE_PARTICLE_M:
        return 232.0 * diameter ** (2.0 / 3.0)
    return 2.5e5 * diameter**1.5
def suspension_kinematic_viscosity(mu_d, rho_g, epsilon_d) -> float:
    """由 μ_d = ν_d ρ_g (1 − ε_d) 解出 ν_d。
    Source: Hamel Eq.3.20, S. 27
    """
    dynamic = _positive("mu_d", mu_d)
    density = _positive("rho_g", rho_g)
    voidage = _finite("epsilon_d", epsilon_d)
    solids = 1.0 - voidage
    if solids <= 0.0:
        raise HydrodynamicsError("ε_d 必须小于 1，否则悬浮相运动粘度没有定义")
    return dynamic / (density * solids)
def bubble_drag_coefficient(reynolds_b) -> float:
    """c_b = 16/Re_b + 2,64。
    Source: Hamel Eq.3.18, S. 26
    """
    reynolds = _positive("reynolds_b", reynolds_b)
    return 16.0 / reynolds + 2.64
def isolated_bubble_velocity(diameter_m, drag_coefficient) -> float:
    """单气泡速度。
    Source: Hamel Eq.3.17, S. 26
    """
    diameter = _positive("diameter_m", diameter_m)
    drag = _positive("drag_coefficient", drag_coefficient)
    return ((1.0 / drag) * (4.0 / 3.0) * diameter * G) ** 0.5
def psi_b(distributor, geldart) -> float:
    """可见气泡流份额。只对 Geldart B 给出数值。
    Source: Hamel Eq.3.14, S. 26
    """
    if geldart != "B":
        raise HydrodynamicsError("式 (3.14) 只对 Geldart B 有定义")
    if distributor == "technical":
        return 0.76
    if distributor == "porous":
        return 0.67
    raise HydrodynamicsError("布风器必须是 technical 或 porous")
def upsilon_b(vessel_diameter_m, geldart) -> float:
    """气泡群修正。d_B 是容器直径。只对 Geldart B 给出数值。
    Source: Hamel Eq.3.16, S. 26
    """
    if geldart != "B":
        raise HydrodynamicsError("式 (3.16) 只对 Geldart B 有定义")
    diameter = _positive("vessel_diameter_m", vessel_diameter_m)
    if diameter < 0.1:
        return 0.64
    if diameter <= 1.0:
        return 2.0 * diameter**0.5
    return 2.0
@dataclass(frozen=True)
class HydroSearch:
    """床区两个标量子问题的搜索参数。求解手段：ADR 0020。不设默认值。"""
    velocity_floor: float
    velocity_margin: float
    rise_bisection_steps: int
    rise_residual_tolerance: float
    density_start_factor: float
    density_growth_factor: float
    density_scan_limit: int
    density_cap: float
    density_bisection_steps: int
    density_mass_tolerance: float
def _count(name: str, value: object) -> int:
    number = _finite(name, value)
    if number < 1.0 or int(number) != number:
        raise HydrodynamicsError(f"{name} 必须是正整数")
    return int(number)
def _use_search(search: HydroSearch) -> HydroSearch:
    if not isinstance(search, HydroSearch):
        raise HydrodynamicsError("搜索参数必须由调用方作为 HydroSearch 传入")
    _positive("velocity_floor", search.velocity_floor)
    _positive("velocity_margin", search.velocity_margin)
    _count("rise_bisection_steps", search.rise_bisection_steps)
    _positive("rise_residual_tolerance", search.rise_residual_tolerance)
    _positive("density_start_factor", search.density_start_factor)
    growth = _finite("density_growth_factor", search.density_growth_factor)
    if growth <= 1.0:
        raise HydrodynamicsError("density_growth_factor 必须大于 1")
    _count("density_scan_limit", search.density_scan_limit)
    _positive("density_cap", search.density_cap)
    _count("density_bisection_steps", search.density_bisection_steps)
    _positive("density_mass_tolerance", search.density_mass_tolerance)
    return search
def bubble_rise_velocity(excess_term, upsilon, diameter_m, nu_d, search: HydroSearch) -> float:
    """解 u_b = excess_term + υ_b u_b,i(u_b)。区间、轮数和残差容差由调用方传入。
    求解手段：ADR 0020。采用式的 excess_term 由 local_bed_state 写成 ψ_b(u_0−u_mf)(P/P_0)^0,2。
    Source: Hamel Eq.3.13, S. 26
    Source: Hamel Eq.3.19, S. 27
    """
    used = _use_search(search)
    excess = _finite("excess_term", excess_term)
    swarm = _positive("upsilon", upsilon)
    diameter = _positive("diameter_m", diameter_m)
    kinematic = _positive("nu_d", nu_d)
    isolated_scale = ((4.0 / 3.0) * diameter * G / 2.64) ** 0.5
    lo = used.velocity_floor
    hi = abs(excess) + swarm * isolated_scale + used.velocity_margin
    def residual(velocity: float) -> float:
        reynolds = velocity * diameter / kinematic
        isolated = isolated_bubble_velocity(diameter, bubble_drag_coefficient(reynolds))
        return excess + swarm * isolated - velocity
    low = residual(lo)
    high = residual(hi)
    if not np.isfinite(low) or not np.isfinite(high) or low * high > 0.0:
        raise HydrodynamicsError(
            "在登记的气泡速度搜索区间内没有找到变号。这是数值搜索失败，不能据此断言原式没有根"
        )
    for _ in range(used.rise_bisection_steps):
        mid = 0.5 * (lo + hi)
        value = residual(mid)
        if low * value <= 0.0:
            hi = mid
            high = value
        else:
            lo = mid
            low = value
    root = 0.5 * (lo + hi)
    gap = abs(residual(root))
    if gap > used.rise_residual_tolerance * max(abs(root), used.velocity_floor):
        raise HydrodynamicsError(
            "气泡速度二分后原式残差仍大于登记容差。这是数值搜索失败，不能据此断言原式没有根"
        )
    return root
def alpha_b(u_b, u_d, epsilon_d) -> float:
    """α_b = u_b ε_d / u_d。解释性选择：ADR 0011 第 7 条。
    Source: Hamel Eq.3.6, S. 23
    """
    bubble = _finite("u_b", u_b)
    suspension = _positive("u_d", u_d)
    voidage = _finite("epsilon_d", epsilon_d)
    return bubble * voidage / suspension
def coalescence_parameter(alpha) -> float:
    """ξ_b。α_b < 0 时原文没有定义。
    Source: Hamel Eq.3.34, S. 29
    """
    value = _finite("alpha", alpha)
    if value < 0.0:
        raise HydrodynamicsError("α_b < 0，式 (3.34) 没有定义")
    if value <= 1.0:
        return 1.0 - value**3
    return 0.0
def bubble_lifetime(u_mf, pressure) -> float:
    """λ_b = 280 (u_mf/g) (P/P_0)^(−0,7)。P = P_0 时回到式 (3.35)。
    Source: Hamel Eq.3.42, S. 31
    Source: Hamel Eq.3.35, S. 29
    """
    minimum = _positive("u_mf", u_mf)
    absolute = _positive("pressure", pressure)
    return 280.0 * minimum / G * (absolute / P_0) ** (-0.7)
def bubble_growth_rate(diameter_m, epsilon_b, xi_b, lifetime_s, u_b, pressure) -> float:
    """加压后的 dd_b/dh。P = P_0 时第一项回到式 (3.33)。
    Source: Hamel Eq.3.41, S. 31
    Source: Hamel Eq.3.33, S. 28
    """
    diameter = _positive("diameter_m", diameter_m)
    fraction = _finite("epsilon_b", epsilon_b)
    coalescence = _finite("xi_b", xi_b)
    lifetime = _positive("lifetime_s", lifetime_s)
    velocity = _positive("u_b", u_b)
    absolute = _positive("pressure", pressure)
    if fraction <= 0.0:
        raise HydrodynamicsError("ε_b 必须大于 0")
    pressure_factor = (absolute / P_0) ** (P_0 / absolute)
    denominator = 1.0 - coalescence * (6.0 / np.pi) ** (1.0 / 3.0) * fraction ** (1.0 / 3.0)
    if denominator <= 0.0:
        raise HydrodynamicsError("气泡长大式的分母不是正数")
    coalescence_term = (2.0 / (9.0 * np.pi)) ** (1.0 / 3.0) * (
        fraction ** (1.0 / 3.0) * pressure_factor / denominator
    )
    splitting_term = diameter / (3.0 * lifetime * velocity)
    return coalescence_term - splitting_term
def orifice_volumetric_flow(u_0, area_m2, n_0) -> float:
    """单孔体积流量 V̇_0 = u_0 A / n_0。
    Source: Hamel Eq.3.37, S. 29
    """
    superficial = _positive("u_0", u_0)
    area = _positive("area_m2", area_m2)
    count = _positive("n_0", n_0)
    if count != int(count) or int(count) < 1:
        raise HydrodynamicsError("n_0 必须是正整数")
    return superficial * area / count
def initial_bubble_diameter(distributor, u_0, area_m2, n_0, u_mf) -> float:
    """技术布风器用式 (3.36)、(3.37)；多孔板用式 (3.38)。
    Source: Hamel Eq.3.36, S. 29
    Source: Hamel Eq.3.38, S. 29
    """
    if distributor == "technical":
        flow = orifice_volumetric_flow(u_0, area_m2, n_0)
        return 1.3 * (flow**2 / G) ** 0.2
    if distributor == "porous":
        superficial = _finite("u_0", u_0)
        minimum = _finite("u_mf", u_mf)
        if superficial <= minimum:
            raise HydrodynamicsError("u_0 必须大于 u_mf")
        return 2.78 / G * (superficial - minimum) ** 2
    raise HydrodynamicsError("布风器必须是 technical 或 porous")
def growth_start_height(distributor, d_b0, h_a) -> float:
    """技术布风器 h_0 = h_A + d_b,0/2；多孔板 h_0 = d_b,0/2。
    Source: Hamel Eq.3.39, S. 29
    Source: Hamel Eq.3.40, S. 29
    """
    diameter = _positive("d_b0", d_b0)
    if distributor == "technical":
        orifice_height = _finite("h_a", h_a)
        if orifice_height < 0.0:
            raise HydrodynamicsError("h_A 不能为负")
        return orifice_height + 0.5 * diameter
    if distributor == "porous":
        return 0.5 * diameter
    raise HydrodynamicsError("布风器必须是 technical 或 porous")
def exchange_coefficient(u_br, diameter_m, diffusivity, epsilon_mf, u_b) -> float:
    """本文采用的相间交换系数，S. 35 的式 (3.50)，系数 3/2。
    Source: Hamel Eq.3.50, S. 35
    """
    through = _finite("u_br", u_br)
    diameter = _positive("diameter_m", diameter_m)
    diffusion = _positive("diffusivity", diffusivity)
    voidage = _positive("epsilon_mf", epsilon_mf)
    velocity = _positive("u_b", u_b)
    if not 0.0 < voidage < 1.0:
        raise HydrodynamicsError("ε_mf 必须在 (0, 1) 内")
    convective = 1.5 * through / diameter
    diffusive = (144.0 * diffusion * voidage * velocity / (np.pi * diameter**3)) ** 0.5
    return convective + diffusive
@dataclass(frozen=True)
class BedConditions:
    """一次床层计算的物性与装置条件。ρ_g、μ_g、ν_g、D_g 都由调用方传入。"""
    d_p: float
    epsilon_mf: float
    rho_p: float
    rho_g: float
    mu_g: float
    nu_g: float
    pressure: float
    diffusivity: float
    vessel_diameter: float
    distributor: str
    geldart: str
    three_dimensional: bool
@dataclass(frozen=True)
class LocalState:
    """一个高度上的采用式结果。"""
    u_d: float
    u_b: float
    u_br: float
    epsilon_b: float
    epsilon_d: float
    alpha: float
    xi_b: float
    d_b: float
    growth_per_height: float
    k_bd: float
def local_bed_state(u_0, d_b, u_mf, conditions: BedConditions, search: HydroSearch) -> LocalState:
    """给定高度上的速度、份额、长大速率和 K_bd。
    压力系数是 ψ_b·(P/P_0)^0,2，不是印刷式里的减 1。解释性选择：ADR 0011 第 10 条。
    α_b 用式 (3.6) 的采用读法。解释性选择：ADR 0011 第 7 条。
    Source: Hamel Eq.3.43, S. 31
    """
    superficial = _finite("u_0", u_0)
    diameter = _positive("d_b", d_b)
    minimum = _positive("u_mf", u_mf)
    if superficial <= minimum:
        raise HydrodynamicsError("u_0 必须大于 u_mf")
    suspension = suspension_superficial_velocity(
        minimum,
        superficial,
        conditions.distributor,
        conditions.geldart,
        conditions.three_dimensional,
    )
    through = throughflow_velocity(suspension, conditions.pressure)
    terminal = terminal_velocity(conditions.d_p, conditions.nu_g, conditions.rho_p, conditions.rho_g)
    n_rz = richardson_zaki_exponent(
        settling_reynolds(terminal, conditions.d_p, conditions.rho_g, conditions.mu_g)
    )
    voidage = suspension_voidage(conditions.epsilon_mf, suspension, minimum, n_rz)
    kinematic = suspension_kinematic_viscosity(
        suspension_dynamic_viscosity(conditions.d_p),
        conditions.rho_g,
        voidage,
    )
    factor = (conditions.pressure / P_0) ** 0.2
    excess = psi_b(conditions.distributor, conditions.geldart) * (superficial - minimum) * factor
    velocity = bubble_rise_velocity(
        excess,
        upsilon_b(conditions.vessel_diameter, conditions.geldart),
        diameter,
        kinematic,
        search,
    )
    fraction = visible_bubble_fraction(superficial, suspension, velocity, through)
    if fraction <= 0.0 or fraction >= 1.0:
        raise HydrodynamicsError("ε_b 不在 (0, 1) 内")
    characteristic = alpha_b(velocity, suspension, voidage)
    coalescence = coalescence_parameter(characteristic)
    lifetime = bubble_lifetime(minimum, conditions.pressure)
    growth = bubble_growth_rate(
        diameter, fraction, coalescence, lifetime, velocity, conditions.pressure
    )
    coefficient = exchange_coefficient(
        through, diameter, conditions.diffusivity, conditions.epsilon_mf, velocity
    )
    return LocalState(
        suspension,
        velocity,
        through,
        fraction,
        voidage,
        characteristic,
        coalescence,
        diameter,
        growth,
        coefficient,
    )
def integrate_growth(d_start, h_start, h_end, n_steps, derivative) -> float:
    """对 dd_b/dh 做定步长 RK4。derivative(d_b) 返回式 (3.41) 的右端。
    Source: Hamel Eq.3.41, S. 31
    """
    start = _positive("d_start", d_start)
    begin = _finite("h_start", h_start)
    end = _finite("h_end", h_end)
    steps = int(_finite("n_steps", n_steps))
    if steps < 1 or steps != _finite("n_steps", n_steps):
        raise HydrodynamicsError("积分步数必须是正整数")
    if end <= begin:
        return start
    step = (end - begin) / steps
    diameter = start
    for _ in range(steps):
        diameter = _rk4(diameter, step, derivative)
    return diameter
def _rk4(diameter: float, step: float, derivative) -> float:
    k1 = derivative(diameter)
    k2 = derivative(diameter + 0.5 * step * k1)
    k3 = derivative(diameter + 0.5 * step * k2)
    k4 = derivative(diameter + step * k3)
    updated = diameter + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if not np.isfinite(updated) or updated <= 0.0:
        raise HydrodynamicsError("气泡直径在积分中变为非正或非有限")
    return updated
def state_at_height(
    height_m,
    u_0,
    u_mf,
    area_m2,
    n_0,
    h_a,
    n_steps,
    conditions: BedConditions,
    search: HydroSearch,
) -> LocalState:
    """从布风器积分到指定高度。h ≤ h_0 时直径为 d_b,0。
    Source: Hamel Eq.3.41, S. 31
    Source: Hamel Eq.3.39, S. 29
    """
    height = _finite("height_m", height_m)
    if height < 0.0:
        raise HydrodynamicsError("高度从布风器算起，不能为负")
    superficial = _positive("u_0", u_0)
    minimum = _positive("u_mf", u_mf)
    d_b0 = initial_bubble_diameter(
        conditions.distributor, superficial, area_m2, n_0, minimum
    )
    h_0 = growth_start_height(conditions.distributor, d_b0, h_a)
    if height <= h_0:
        return local_bed_state(superficial, d_b0, minimum, conditions, search)
    def derivative(diameter: float) -> float:
        return local_bed_state(superficial, diameter, minimum, conditions, search).growth_per_height
    diameter = integrate_growth(d_b0, h_0, height, n_steps, derivative)
    return local_bed_state(superficial, diameter, minimum, conditions, search)
@dataclass(frozen=True)
class BedProfile:
    """下标 0 是最上层 cell，与 hamel.balances 一致。
    最上层的上边界是床面，不在这个数组里。
    """
    u_mf: float
    epsilon_b: np.ndarray
    u_b: np.ndarray
    d_b: np.ndarray
    epsilon_d: np.ndarray
    k_bd: np.ndarray
    solids_mass: np.ndarray
    internal_wake_upflow: np.ndarray
def bed_profile(
    cell_edges_m,
    u_0,
    u_mf,
    area_m2,
    n_0,
    h_a,
    intervals_per_cell,
    conditions: BedConditions,
    search: HydroSearch,
) -> BedProfile:
    """各 cell 上边界的流动量，以及按 cell 积分的几何存量。
    cell_edges_m 从布风器向上，长度为 N+1。u_0[0] 是最上层。
    截面积在床区内取常数。存量解释性选择：ADR 0009；Wirsum (1998) 式 (3.26)。
    Source: Hamel Eq.3.10, S. 24
    Source: Hamel Eq.3.41, S. 31
    """
    edges = _edges(cell_edges_m)
    velocities = _velocities(u_0, edges.size - 1)
    count = edges.size - 1
    intervals = int(_finite("intervals_per_cell", intervals_per_cell))
    if intervals < 1 or intervals != _finite("intervals_per_cell", intervals_per_cell):
        raise HydrodynamicsError("每个 cell 的区间数必须是正整数")
    area = _positive("area_m2", area_m2)
    minimum = _positive("u_mf", u_mf)
    bottom_velocity = float(velocities[-1])
    d_b0 = initial_bubble_diameter(
        conditions.distributor, bottom_velocity, area, n_0, minimum
    )
    h_0 = growth_start_height(conditions.distributor, d_b0, h_a)
    if h_0 < float(edges[0]):
        raise HydrodynamicsError("h_0 低于给定床区，起始直径到床底之间的长大没有格子可积")
    epsilon_b = np.empty(count)
    u_b = np.empty(count)
    d_b = np.empty(count)
    epsilon_d = np.empty(count)
    k_bd = np.empty(count)
    solids = np.empty(count)
    entering = d_b0
    for index in range(count - 1, -1, -1):
        lo = float(edges[count - 1 - index])
        hi = float(edges[count - index])
        samples, states, entering = _march_cell(
            lo, hi, intervals, h_0, d_b0, entering, float(velocities[index]), minimum, conditions, search
        )
        top = states[-1]
        epsilon_b[index] = top.epsilon_b
        u_b[index] = top.u_b
        d_b[index] = top.d_b
        epsilon_d[index] = top.epsilon_d
        k_bd[index] = top.k_bd
        fraction = np.array([(1.0 - item.epsilon_b) * (1.0 - item.epsilon_d) for item in states])
        solids[index] = conditions.rho_p * area * float(np.trapz(fraction, samples))
        if solids[index] <= 0.0:
            raise HydrodynamicsError(f"cell {index} 的几何存量不是正数")
    internal = wake_solids_upflow(
        epsilon_b[1:], u_b[1:], conditions.rho_p, epsilon_d[1:], area
    )
    return BedProfile(minimum, epsilon_b, u_b, d_b, epsilon_d, k_bd, solids, internal)
def _march_cell(lo, hi, intervals, h_0, d_b0, entering, u_0, u_mf, conditions, search):
    if hi <= lo:
        raise HydrodynamicsError("cell 上边界必须高于下边界")
    nodes = np.linspace(lo, hi, intervals + 1)
    diameters = np.empty(intervals + 1)
    diameters[0] = d_b0 if lo <= h_0 else entering
    for step in range(intervals):
        begin = float(nodes[step])
        end = float(nodes[step + 1])
        if end <= h_0:
            diameters[step + 1] = d_b0
            continue
        start_height = begin
        start_diameter = float(diameters[step])
        if begin < h_0:
            start_height = h_0
            start_diameter = d_b0
        def derivative(diameter, velocity=u_0, minimum=u_mf, setup=conditions, used=search):
            return local_bed_state(velocity, diameter, minimum, setup, used).growth_per_height
        diameters[step + 1] = _rk4(start_diameter, end - start_height, derivative)
    states = [
        local_bed_state(u_0, float(diameters[point]), u_mf, conditions, search)
        for point in range(intervals + 1)
    ]
    return nodes, states, float(diameters[-1])
def _edges(values) -> np.ndarray:
    try:
        edges = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise HydrodynamicsError("cell_edges_m 不能转换成浮点数组") from exc
    if edges.ndim != 1 or edges.size < 2:
        raise HydrodynamicsError("cell_edges_m 至少要有两个高度")
    if not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0.0):
        raise HydrodynamicsError("cell_edges_m 必须严格递增且有限")
    return edges
def _velocities(values, count: int) -> np.ndarray:
    try:
        velocity = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise HydrodynamicsError("u_0 不能转换成浮点数组") from exc
    if velocity.shape != (count,):
        raise HydrodynamicsError("u_0 的长度必须等于 cell 数，下标 0 为最上层")
    if not np.isfinite(velocity).all():
        raise HydrodynamicsError("u_0 含 NaN 或无穷大")
    return velocity
def wake_solids_upflow(epsilon_b, u_b, rho_p, epsilon_d, area_m2) -> np.ndarray:
    """床内边界的尾涡向上固体流。f_w = 0,25（S. 39）。不要用于床面。
    解释性选择：ADR 0009。床面出流是 F_0·A，由调用方另行传入。
    Source: Hamel Eq.3.54, S. 37
    """
    fraction = np.asarray(epsilon_b, dtype=float)
    velocity = np.asarray(u_b, dtype=float)
    voidage = np.asarray(epsilon_d, dtype=float)
    density = _positive("rho_p", rho_p)
    area = _positive("area_m2", area_m2)
    if fraction.shape != velocity.shape or fraction.shape != voidage.shape:
        raise HydrodynamicsError("ε_b、u_b、ε_d 的形状必须相同")
    if fraction.ndim != 1:
        raise HydrodynamicsError("床内边界量必须是一维数组")
    if not np.isfinite(fraction).all() or not np.isfinite(velocity).all() or not np.isfinite(voidage).all():
        raise HydrodynamicsError("尾涡流的输入必须有限")
    if np.any(fraction <= 0.0):
        raise HydrodynamicsError("床内边界的 ε_b 必须大于 0")
    return F_W * fraction * velocity * density * (1.0 - voidage) * area
def bed_upflow(bed_surface_kg_s, internal_wake_kg_s) -> np.ndarray:
    """完整的向上固体流。下标 0 是调用方给出的床面值，后面是床内尾涡流。
    Source: Hamel Eq.3.54, S. 37
    """
    surface = _finite("bed_surface_kg_s", bed_surface_kg_s)
    internal = np.asarray(internal_wake_kg_s, dtype=float)
    if internal.ndim != 1 or not np.isfinite(internal).all():
        raise HydrodynamicsError("床内尾涡流必须是有限的一维数组")
    full = np.empty(internal.size + 1)
    full[0] = surface
    full[1:] = internal
    return full
def particle_density_from_bed_mass(
    bed_mass_kg,
    cell_edges_m,
    u_0,
    area_m2,
    n_0,
    h_a,
    intervals_per_cell,
    phi_s,
    conditions: BedConditions,
    search: HydroSearch,
) -> float:
    """使床区几何存量之和等于给定床料质量的 ρ_P。
    ε_b、ε_d 随 ρ_P 变化。扫描范围、轮数和闭合容差由调用方传入。求解手段：ADR 0020。
    解释性选择：ADR 0014 A4；ADR 0009。
    Source: Hamel Eq.3.4, S. 22
    Source: Hamel Eq.3.10, S. 24
    """
    target = _positive("bed_mass_kg", bed_mass_kg)
    sphericity = _positive("phi_s", phi_s)
    used = _use_search(search)
    def solids_mass(density: float):
        updated = _with_density(conditions, density)
        minimum = minimum_fluidization_velocity(
            updated.epsilon_mf,
            updated.nu_g,
            sphericity,
            updated.d_p,
            density,
            updated.rho_g,
        )
        profile = bed_profile(
            cell_edges_m,
            u_0,
            minimum,
            area_m2,
            n_0,
            h_a,
            intervals_per_cell,
            updated,
            used,
        )
        return float(np.sum(profile.solids_mass))
    found = []
    density = conditions.rho_g * used.density_start_factor
    for _ in range(used.density_scan_limit):
        if density > used.density_cap:
            break
        try:
            found.append((density, solids_mass(density)))
        except HydrodynamicsError:
            pass
        density *= used.density_growth_factor
    below = [item for item in found if item[1] <= target]
    above = [item for item in found if item[1] >= target]
    if not below or not above:
        raise HydrodynamicsError(
            "在登记的密度扫描范围内没有跨过给定床料质量。这是数值搜索失败，不能据此断言质量闭合没有物理解"
        )
    low, low_mass = max(below, key=lambda item: item[0])
    high, high_mass = min(above, key=lambda item: item[0])
    if abs(low_mass - target) <= used.density_mass_tolerance * target:
        return low
    for _ in range(used.density_bisection_steps):
        mid = 0.5 * (low + high)
        try:
            mid_mass = solids_mass(mid)
        except HydrodynamicsError as exc:
            raise HydrodynamicsError("密度区间内出现不能按鼓泡床计算的点") from exc
        if mid_mass < target:
            low, low_mass = mid, mid_mass
        else:
            high, high_mass = mid, mid_mass
    final = 0.5 * (low + high)
    final_mass = solids_mass(final)
    if abs(final_mass - target) > used.density_mass_tolerance * target:
        raise HydrodynamicsError(
            "密度二分后床料质量残差仍大于登记容差。这是数值搜索失败，不能据此断言质量闭合没有物理解"
        )
    return final
def _with_density(conditions: BedConditions, density: float) -> BedConditions:
    return BedConditions(
        conditions.d_p,
        conditions.epsilon_mf,
        density,
        conditions.rho_g,
        conditions.mu_g,
        conditions.nu_g,
        conditions.pressure,
        conditions.diffusivity,
        conditions.vessel_diameter,
        conditions.distributor,
        conditions.geldart,
        conditions.three_dimensional,
    )
def geldart_group(surface_volume_diameter_m: float, rho_s_kg_m3: float, rho_f_kg_m3: float) -> str:
    """常压 Geldart（1973）式 (6)、(8)。密度差用 g/cm³，粒径用 μm。

    Hamel 没有印出这两条分界，只在采用式里要求 Gruppe-B。式 (3.14) 不是分界式。
    ADR 0024 第 4 项。不是 B 类由调用方停止。
    Source: Hamel Eq.3.14, S. 26
    """
    diameter = _positive("d_sv", surface_volume_diameter_m)
    solid = _finite("rho_s", rho_s_kg_m3)
    fluid = _finite("rho_f", rho_f_kg_m3)
    difference = (solid - fluid) / 1000.0
    if difference <= 0.0:
        raise HydrodynamicsError("Geldart 判据的密度差不是正数")
    micron = diameter * 1.0e6
    if difference * micron <= 225.0:
        return "A"
    if difference * micron**2 >= 1.0e6:
        return "D"
    return "B"
def _require_technical_three_dimensional(distributor: str, three_dimensional: bool) -> None:
    if distributor != "technical":
        raise HydrodynamicsError("式 (3.11) 只用于技术布风器；多孔板原文没有采用式")
    if three_dimensional is not True:
        raise HydrodynamicsError("式 (3.11) 只用于三维床")
