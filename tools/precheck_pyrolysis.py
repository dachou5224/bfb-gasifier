#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单颗粒干燥与热解的平均释放时间预检。不进入 hamel/。

T_w = ∫ Y_w dt，T_v = ∫ Y_v dt（ADR 0015 第 10 项第 4 步）。
Y_w 是剩余水分份额，Y_v = 1 − V/V_∞ 是剩余挥发分份额。
有效域是 T_v ≥ T_w。这里只报告是否成立，不把 τ_pyr 改成 0。

式 (4.3) 取负号，式 (4.11) 的指数取负号，式 (4.12) 的前置因子取 3/R₀³。
这三处是 ADR 0015 第 5–7 条的读法。湿核内不热解：§4.1 写明热解只发生在已干壳层。

同步干燥与热解的温度来自已干壳层的非稳态导热（式 4.1，边界 4.2、4.3）。
S. 53 写明这一段不用准稳态剖面；S. 54 写明式 (4.12) 的温度由该非稳态模型给定。
原文说把式 (4.1) 变换后假设二次温度函数，系数没有印出，这里不补系数。

物性在 T21 未合入时用表 4.1 和临时气体数据，出处写在 INPUT_NOTES。
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

G_NOTE = "预检，不是热解闭合的裁定"


def release_time(times, values, tail: float = 0.0) -> float:
    """∫ Y dt。times 从过程起点计，tail 是最后一个采样之后到无穷的积分。

    指数衰减 Y = exp(−t/τ) 在采样末端的尾部是 τ·Y_end，整段积分等于 τ。
    """
    grid = np.asarray(times, dtype=float)
    series = np.asarray(values, dtype=float)
    if grid.ndim != 1 or series.shape != grid.shape or grid.size < 2:
        raise ValueError("释放时间需要至少两个采样")
    if np.any(np.diff(grid) <= 0.0):
        raise ValueError("时间采样必须严格递增")
    return float(np.trapz(series, grid) + tail)


def exponential_tail(remaining: float, time_constant: float) -> float:
    """Y(t1)·τ，即 exp(−(t−t1)/τ) 从 t1 积到无穷。"""
    if time_constant <= 0.0:
        raise ValueError("时间常数必须为正")
    return float(remaining * time_constant)


@dataclass(frozen=True)
class ParticleProperties:
    """单颗粒预检物性。除注明原文表的项外都是临时值。"""

    temperature_k: float
    pressure_pa: float
    moisture_dry_basis: float
    initial_temperature_k: float
    evaporation_temperature_k: float
    latent_heat_j_kg: float
    water_heat_j_kg_k: float
    solid_heat_j_kg_k: float
    density_kg_m3: float
    diffusivity_m2_s: float
    activation_mean_j_mol: float
    activation_sigma_j_mol: float
    frequency_s: float
    gas_conductivity_w_m_k: float
    gas_viscosity_pa_s: float
    gas_density_kg_m3: float
    prandtl: float
    superficial_m_s: float

    @property
    def conductivity(self) -> float:
        """λ_s = a ρ c。表 4.1 给了 a、ρ_P、c_P，没有单独的 λ。临时把三者当成同一套。"""
        return self.diffusivity_m2_s * self.density_kg_m3 * self.solid_heat_j_kg_k

    @property
    def corrected_latent_heat(self) -> float:
        """式 (4.4)。Source: Hamel Eq.4.4, S. 52"""
        water = self.moisture_dry_basis
        capacity = self.water_heat_j_kg_k + self.solid_heat_j_kg_k / water
        return self.latent_heat_j_kg + capacity * (self.evaporation_temperature_k - self.initial_temperature_k)


def nusselt(particle: ParticleProperties, diameter_m: float) -> float:
    """Nu = 2 + 1,2 Re^{1/2} Pr^{1/3}。特征速度临时取床层表观气速。

    Source: Hamel Eq.4.9, S. 53
    """
    reynolds = (
        particle.gas_density_kg_m3 * particle.superficial_m_s * diameter_m / particle.gas_viscosity_pa_s
    )
    return 2.0 + 1.2 * math.sqrt(max(reynolds, 0.0)) * particle.prandtl ** (1.0 / 3.0)


def heat_transfer_coefficient(particle: ParticleProperties, diameter_m: float) -> float:
    return nusselt(particle, diameter_m) * particle.gas_conductivity_w_m_k / diameter_m


def front_rate(phi: float, particle: ParticleProperties, radius_m: float) -> float:
    """dφ/dt，φ = r_e/R。式 (4.3) 按 ADR 0015 取负号，使锋面内移时热流为正。

    φ = 1 处的极限是 −α (T_a − T_e) / (R ρ h'_v w)。
    这是式 (4.5) 的准稳态降阶，只用来核对表面热流极限，不进入同步阶段的温度。
    Source: Hamel Eq.4.3, S. 52
    Source: Hamel Eq.4.5, S. 53
    """
    if phi <= 0.0:
        return 0.0
    alpha = heat_transfer_coefficient(particle, 2.0 * radius_m)
    biot = alpha * radius_m / particle.conductivity
    span = particle.temperature_k - particle.evaporation_temperature_k
    denominator = (
        radius_m**2
        * particle.corrected_latent_heat
        * particle.moisture_dry_basis
        * particle.density_kg_m3
        * phi
        * (phi + biot * (1.0 - phi))
    )
    return -particle.conductivity * biot * span / denominator


def _energy_nodes(particle: ParticleProperties, count: int = 12):
    nodes, weights = np.polynomial.hermite_e.hermegauss(count)
    scale = math.sqrt(2.0 * math.pi)
    return particle.activation_mean_j_mol + particle.activation_sigma_j_mol * nodes, weights / scale


def _remaining(exposure, frequency: float, weights) -> np.ndarray:
    decay = np.exp(np.clip(-frequency * exposure, -700.0, 0.0))
    return decay @ weights


def average_volatile(radius, remaining_profile, outer_m: float) -> float:
    """式 (4.12) 的体积平均。连续式的前置因子是 3/R³（ADR 0015）。

    径向用梯形公式。梯形对 r² 不是精确的，所以用同一套权重归一，均匀剩余份额得到 1。
    Source: Hamel Eq.4.12, S. 54
    """
    del outer_m
    weight = np.asarray(radius, dtype=float) ** 2
    profile = np.asarray(remaining_profile, dtype=float)
    return float(np.trapz(profile * weight, radius) / np.trapz(weight, radius))


# 初始干壳只为避开 ξ 变换在锋面贴着表面时的零厚度。临时离散，不是模型参数。
INITIAL_SHELL_FRACTION = 1.0e-3
FRONT_STOP_FRACTION = 1.0e-3
SHELL_INTERVALS = 12


def _shell_nodes(state, evaporation_k: float):
    """state[0] 是 r_e，其后是 ξ = 1/N … 1 上的温度。ξ = 0 固定为 T_e。"""
    values = np.empty(state.size)
    values[0] = evaporation_k
    values[1:] = state[1:]
    return values


def drying_rhs(time, state, particle: ParticleProperties, radius_m: float, alpha: float):
    """已干壳层上的式 (4.1)。ξ = (r − r_e)/(R − r_e) 把移动壳层固定到 [0, 1]。

    锋面速度来自式 (4.3) 的 ADR 0015 读法：导热通量为正时锋面内移。
    Source: Hamel Eq.4.1, S. 52
    Source: Hamel Eq.4.2, S. 52
    Source: Hamel Eq.4.3, S. 52
    """
    del time
    front = float(state[0])
    thickness = radius_m - front
    intervals = state.size - 1
    out = np.zeros_like(state)
    if thickness <= 0.0 or front <= FRONT_STOP_FRACTION * radius_m:
        return out
    values = _shell_nodes(state, particle.evaporation_temperature_k)
    xi_step = 1.0 / intervals
    surface_slope = (alpha * thickness / particle.conductivity) * (
        particle.temperature_k - values[-1]
    )
    ghost = values[-2] + 2.0 * xi_step * surface_slope
    front_slope = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (2.0 * xi_step)
    capacity = particle.corrected_latent_heat * particle.moisture_dry_basis * particle.density_kg_m3
    speed = -(particle.conductivity / capacity) * front_slope / thickness
    diffusivity = particle.diffusivity_m2_s
    for index in range(1, intervals + 1):
        left = values[index - 1]
        center = values[index]
        right = ghost if index == intervals else values[index + 1]
        dxi = (right - left) / (2.0 * xi_step)
        dxi2 = (right - 2.0 * center + left) / xi_step**2
        xi = index * xi_step
        radial = front + xi * thickness
        conduction = diffusivity * (dxi2 / thickness**2 + 2.0 * dxi / (radial * thickness))
        motion = speed * (1.0 - xi) / thickness * dxi
        out[index] = conduction + motion
    out[0] = speed
    return out


def integrate_drying(particle: ParticleProperties, radius_m: float, intervals: int = SHELL_INTERVALS):
    """把锋面从贴近表面积到接近中心。失败时返回 None，不改用准稳态剖面。"""
    if intervals < 4:
        raise ValueError("干壳离散至少要 4 个区间，锋面单侧差分才用得到两个内点")
    alpha = heat_transfer_coefficient(particle, 2.0 * radius_m)
    front0 = radius_m * (1.0 - INITIAL_SHELL_FRACTION)
    initial = np.full(intervals + 1, particle.evaporation_temperature_k)
    initial[0] = front0
    horizon = max(600.0, 800.0 * radius_m**2 / particle.diffusivity_m2_s)

    def center(_time, state):
        return float(state[0] - FRONT_STOP_FRACTION * radius_m)

    def rhs(time, state):
        return drying_rhs(time, state, particle, radius_m, alpha)

    center.terminal = True
    center.direction = -1
    tolerances = np.full(intervals + 1, 1.0e-3)
    tolerances[0] = 1.0e-10
    solved = solve_ivp(
        rhs,
        (0.0, horizon),
        initial,
        method="BDF",
        dense_output=True,
        events=center,
        rtol=1.0e-5,
        atol=tolerances,
    )
    if (not solved.success) or solved.sol is None or solved.t_events[0].size == 0:
        return None
    return solved


def temperature_on_radius(state, query, particle: ParticleProperties, radius_m: float) -> np.ndarray:
    """把壳层节点插到固定半径上。湿核（r < r_e）保持 T_e，不热解。"""
    front = float(state[0])
    values = _shell_nodes(state, particle.evaporation_temperature_k)
    nodes = np.linspace(0.0, 1.0, values.size)
    thickness = radius_m - front
    temperature = np.empty(np.shape(query), dtype=float)
    flat_query = np.asarray(query, dtype=float).reshape(-1)
    flat_temperature = temperature.reshape(-1)
    for index, radial in enumerate(flat_query):
        if radial <= front or thickness <= 0.0:
            flat_temperature[index] = particle.evaporation_temperature_k
        else:
            flat_temperature[index] = float(np.interp((radial - front) / thickness, nodes, values))
    return temperature


def drying_energy_balance(particle: ParticleProperties, radius_m: float, solved, samples: int = 40) -> dict:
    """同步阶段的积分能量：壳层显热增量 + 锋面潜热 = 表面传入的热。

    潜热用锋面扫过的干基水质量乘 h'_v。显热以 T_e 为参考，式 (4.4) 已经把
    湿核从初温升到 T_e 的热放进 h'_v。准稳态剖面会把显热增量漏掉，这项不平衡。
    """
    grid = np.linspace(float(solved.t[0]), float(solved.t[-1]), samples)
    alpha = heat_transfer_coefficient(particle, 2.0 * radius_m)
    area = 4.0 * math.pi * radius_m**2
    sensible = []
    power = []
    for instant in grid:
        state = solved.sol(instant)
        front = float(state[0])
        values = _shell_nodes(state, particle.evaporation_temperature_k)
        nodes = np.linspace(0.0, 1.0, values.size)
        radial = front + nodes * (radius_m - front)
        excess = (values - particle.evaporation_temperature_k) * radial**2
        sensible.append(
            particle.density_kg_m3
            * particle.solid_heat_j_kg_k
            * 4.0
            * math.pi
            * float(np.trapz(excess, radial))
        )
        power.append(area * alpha * (particle.temperature_k - float(values[-1])))
    entered = np.concatenate([[0.0], np.cumsum(0.5 * (np.asarray(power[:-1]) + np.asarray(power[1:])) * np.diff(grid))])
    latent_factor = (
        particle.corrected_latent_heat
        * particle.moisture_dry_basis
        * particle.density_kg_m3
        * 4.0
        * math.pi
        / 3.0
    )
    front0 = float(solved.sol(grid[0])[0])
    latent = [latent_factor * (front0**3 - float(solved.sol(instant)[0]) ** 3) for instant in grid]
    sensible = np.asarray(sensible)
    latent = np.asarray(latent)
    residual = sensible - sensible[0] + latent - entered
    scale = max(float(np.max(entered)), 1.0)
    return {
        "time_s": grid,
        "relative": float(np.max(np.abs(residual)) / scale),
        "entered_j": entered,
        "sensible_j": sensible,
        "latent_j": latent,
    }


def _heat_rhs(temperature, particle: ParticleProperties, radius_m: float, alpha: float):
    count = temperature.size
    step = radius_m / (count - 1)
    rate = np.empty(count)
    diffusivity = particle.diffusivity_m2_s
    rate[0] = 6.0 * diffusivity * (temperature[1] - temperature[0]) / step**2
    for index in range(1, count - 1):
        radial = index * step
        second = (temperature[index + 1] - 2.0 * temperature[index] + temperature[index - 1]) / step**2
        first = (temperature[index + 1] - temperature[index - 1]) / (2.0 * step)
        rate[index] = diffusivity * (second + 2.0 * first / radial)
    ghost = temperature[count - 2] + (2.0 * step * alpha / particle.conductivity) * (
        particle.temperature_k - temperature[count - 1]
    )
    radial = radius_m
    second = (ghost - 2.0 * temperature[count - 1] + temperature[count - 2]) / step**2
    first = (ghost - temperature[count - 2]) / (2.0 * step)
    rate[count - 1] = diffusivity * (second + 2.0 * first / radial)
    return rate


def single_particle(diameter_m: float, particle: ParticleProperties, radial_count: int = 21) -> dict:
    """在床温 T_a = T_WS 下积分单颗粒，返回 T_w、T_v 和 T_v ≥ T_w 是否成立。

    环境温度用床温：图 4.1 的床温是 T_WS，式 (4.2) 写 T_a。预检把二者当成同一个数，临时。
    """
    if diameter_m <= 0.0:
        raise ValueError("粒径必须为正")
    radius_m = 0.5 * diameter_m
    radius = np.linspace(0.0, radius_m, radial_count)
    energies, weights = _energy_nodes(particle)
    gas_constant = 8.314462618
    exposure = np.zeros((radial_count, energies.size))
    times = [0.0]
    # 初始干壳是数值起点。对应的水分份额从 φ³ 起算，不把这层壳写成尚未干燥。
    water = [(1.0 - INITIAL_SHELL_FRACTION) ** 3]
    volatile = [1.0]
    elapsed = 0.0
    # 湿核推进期间联立非稳态导热。准稳态式 (4.5) 不进入这段温度历程。
    dried = integrate_drying(particle, radius_m)
    if dried is None:
        return _unresolved(diameter_m, "非稳态干壳没有积分到锋面接近中心")
    substeps = 8
    for left, right in zip(dried.t[:-1], dried.t[1:]):
        width = float(right - left) / substeps
        if width <= 0.0:
            continue
        for step_index in range(substeps):
            instant = float(left) + (step_index + 1) * width
            state = dried.sol(instant)
            temperature = temperature_on_radius(state, radius, particle, radius_m)
            dry = radius >= float(state[0])
            exposure[dry] += np.exp(
                -energies / (gas_constant * np.maximum(temperature[dry, None], 250.0))
            ) * width
            profile = np.ones(radial_count)
            profile[dry] = _remaining(exposure[dry], particle.frequency_s, weights)
            times.append(instant)
            water.append((float(state[0]) / radius_m) ** 3)
            volatile.append(average_volatile(radius, profile, radius_m))
    elapsed = float(dried.t[-1])

    alpha = heat_transfer_coefficient(particle, diameter_m)
    initial = temperature_on_radius(dried.y[:, -1], radius, particle, radius_m)

    def rhs(_time, temperature):
        return _heat_rhs(temperature, particle, radius_m, alpha)

    def hot(_time, temperature):
        return float(np.max(np.abs(temperature - particle.temperature_k)) - 0.5)

    hot.terminal = True
    hot.direction = -1
    horizon = 40.0 * radius_m**2 / particle.diffusivity_m2_s
    heated = solve_ivp(
        rhs,
        (0.0, horizon),
        initial,
        method="BDF",
        dense_output=True,
        events=hot,
        rtol=1.0e-4,
        atol=1.0e-4,
        max_step=max(radius_m**2 / particle.diffusivity_m2_s / 40.0, 1.0e-4),
    )
    if not heated.success or heated.t[-1] <= 0.0:
        return _unresolved(diameter_m, "颗粒升温没有积分成功")
    if heated.t_events[0].size == 0:
        return _unresolved(diameter_m, "颗粒在导热时间尺度内没有升到床温")
    heat_times = np.linspace(0.0, float(heated.t[-1]), 80)
    for left, right in zip(heat_times[:-1], heat_times[1:]):
        midpoint = heated.sol(0.5 * (left + right))
        exposure += np.exp(-energies / (gas_constant * np.maximum(midpoint[:, None], 250.0))) * (right - left)
        elapsed_heat = elapsed + right
        profile = _remaining(exposure, particle.frequency_s, weights)
        times.append(float(elapsed_heat))
        water.append(0.0)
        volatile.append(average_volatile(radius, profile, radius_m))
    tail = _volatile_tail(radius, exposure, particle, energies, weights, gas_constant)
    water_time = release_time(times, water, 0.0)
    volatile_time = release_time(times, volatile, tail)
    return {
        "diameter_m": diameter_m,
        "T_w_s": water_time,
        "T_v_s": volatile_time,
        "holds": volatile_time >= water_time,
        "status": "成立" if volatile_time >= water_time else "不成立",
        "note": "",
    }


def _volatile_tail(radius, exposure, particle, energies, weights, gas_constant: float) -> float:
    """温度已到 T_WS 之后，各节点、各活化能的剩余挥发分按指数衰减积到无穷。"""
    coefficient = particle.frequency_s * np.exp(-energies / (gas_constant * particle.temperature_k))
    remaining = np.exp(np.clip(-particle.frequency_s * exposure, -700.0, 0.0))
    per_node = remaining @ (weights / coefficient)
    return average_volatile(radius, per_node, float(radius[-1]))


def _unresolved(diameter_m: float, note: str) -> dict:
    return {
        "diameter_m": diameter_m,
        "T_w_s": math.nan,
        "T_v_s": math.nan,
        "holds": None,
        "status": "未检验",
        "note": note,
    }


# 工况 1。床温与压力按任务卡；其余见 INPUT_NOTES。
BED_TEMPERATURE_K = 819.0 + 273.15
PRESSURE_PA = 500_000.0
MOISTURE_AS_RECEIVED = 0.118  # 表 7.2，11,8 %
PRIMARY_AIR_MOL_S = 0.6419598861822472
STEAM_KG_S = 0.001888888888888889
PURGE_N2_MOL_S = 0.07683689757393693
WATER_MOLAR_MASS = 0.01801528
AIR_MOLAR_MASS = 0.0289647
N2_MOLAR_MASS = 0.0280134
GAS_CONSTANT = 8.314462618
BED_AREA_M2 = math.pi * 0.075**2

INPUT_NOTES = (
    "床温 819 °C、压力 0,5 MPa：KURKELA_1992 表 5 的床温只作预检；压力来自表 7.2。",
    "T_a 取床温。式 (4.2) 写的是 T_a，图 4.1 的床温是 T_WS，原文没有写明二者相同。临时。",
    "水分 11,8 %（收到基，表 7.2）。w_0,tr = 0,118/(1−0,118)，临时换算成干基。",
    "燃料初温 298,15 K：ADR 0014 C2。",
    "T_e = 425 K、蒸发焓 2,108·10⁶ J/kg：0,5 MPa 下水的饱和温度和汽化焓，临时，不是 Hamel 给出的数。",
    "液态水比热 4180 J/(kg·K)：临时。",
    "E_0、σ_E、k_0,P、ρ_P、c_P、a：表 4.1 湿褐煤（Dersch 1994，S. 54）。"
    "临时用于莱茵褐煤。表中的 k_0,P、ρ_P、c_P 与式中的 k_0、ρ_s、c_s 原文没有写成同一组量。",
    "λ_s = a ρ_P c_P：临时，由表 4.1 的导温系数反推。",
    "式 (4.9) 的特征速度原文未说明，临时取床层表观气速（一次风 + 蒸汽 + 吹扫氮，按床温与 0,5 MPa）。",
    "气体粘度用空气 Sutherland 公式，导热用 λ = 2,334·10⁻³ T^{1.5}/(T+164,54) W/(m·K)，Pr = 0,70：临时。",
    "干空气、水、氮的摩尔质量与 R：临时。吹扫氮按 ADR 0006 暂定计入。",
    "同步阶段把式 (4.1) 变到 ξ = (r−r_e)/(R−r_e) 后用 12 个等距区间离散。"
    "二次温度函数的系数原文没有印出，不补。初始干壳厚度和锋面停止位置都是半径的 10⁻³，临时。"
    "锋面到达停止位置后，整颗颗粒继续按式 (4.1)、(4.2) 升温。活化能 12 个 Gauss–Hermite 节点，挥发分平均用 21 个固定半径。",
    "粒级用 kurkela1992_fig1.json 的筛孔上限与下限。<0,1 mm 没有下限，该端未检验。",
)


def case1_properties() -> ParticleProperties:
    steam = STEAM_KG_S / WATER_MOLAR_MASS
    amount = PRIMARY_AIR_MOL_S + steam + PURGE_N2_MOL_S
    mass = PRIMARY_AIR_MOL_S * AIR_MOLAR_MASS + STEAM_KG_S + PURGE_N2_MOL_S * N2_MOLAR_MASS
    molar_mass = mass / amount
    density = PRESSURE_PA * molar_mass / (GAS_CONSTANT * BED_TEMPERATURE_K)
    viscosity = 1.716e-5 * (BED_TEMPERATURE_K / 273.15) ** 1.5 * (273.15 + 110.4) / (BED_TEMPERATURE_K + 110.4)
    conductivity = 2.334e-3 * BED_TEMPERATURE_K**1.5 / (BED_TEMPERATURE_K + 164.54)
    volume = amount * GAS_CONSTANT * BED_TEMPERATURE_K / PRESSURE_PA
    return ParticleProperties(
        temperature_k=BED_TEMPERATURE_K,
        pressure_pa=PRESSURE_PA,
        moisture_dry_basis=MOISTURE_AS_RECEIVED / (1.0 - MOISTURE_AS_RECEIVED),
        initial_temperature_k=298.15,
        evaporation_temperature_k=425.0,
        latent_heat_j_kg=2.108e6,
        water_heat_j_kg_k=4180.0,
        solid_heat_j_kg_k=1256.0,
        density_kg_m3=1250.0,
        diffusivity_m2_s=1.0e-7,
        activation_mean_j_mol=192_000.0,
        activation_sigma_j_mol=40_000.0,
        frequency_s=1.67e13,
        gas_conductivity_w_m_k=conductivity,
        gas_viscosity_pa_s=viscosity,
        gas_density_kg_m3=density,
        prandtl=0.70,
        superficial_m_s=volume / BED_AREA_M2,
    )


def case1_diameters():
    """每个有限筛孔只算一次。最细一档没有下限，单独记为未检验。"""
    rows = [(None, "未检验：最细一档没有下限")]
    for bound in (0.1, 0.25, 0.5, 0.71, 1.0, 2.0, 3.15):
        rows.append((bound * 1.0e-3, f"筛孔 {bound} mm"))
    return rows


def case1_pyrolysis() -> list[dict]:
    particle = case1_properties()
    results = []
    for diameter, label in case1_diameters():
        if diameter is None:
            item = _unresolved(math.nan, label)
        else:
            item = single_particle(diameter, particle)
        item["label"] = label
        results.append(item)
    return results


def format_table(rows: list[dict]) -> str:
    header = "| 粒级 | T_w / s | T_v / s | T_v ≥ T_w |"
    rule = "|---|---|---|---|"
    body = []
    for row in rows:
        if row["status"] == "未检验" and not math.isfinite(row["T_w_s"]):
            body.append(f"| {row['label']} | — | — | 未检验 |")
            continue
        body.append(
            f"| {row['label']} | {row['T_w_s']:.4g} | {row['T_v_s']:.4g} | {row['status']} |"
        )
    return "\n".join([header, rule, *body])


def pyrolysis_sentence(rows: list[dict]) -> str:
    flags = [row["status"] for row in rows]
    if any(flag == "不成立" for flag in flags):
        head = "有粒级不成立"
    elif all(flag == "成立" for flag in flags):
        head = "成立"
    elif any(flag == "成立" for flag in flags) and any(flag == "未检验" for flag in flags):
        head = "已检验的粒级成立，最细档下限未检验"
    else:
        head = "未检验"
    return f"热解有效域预检（T_v ≥ T_w）的结果是「{head}」。{G_NOTE}。"


def main() -> None:
    rows = case1_pyrolysis()
    print("工况 1 热解有效域预检")
    print(format_table(rows))
    print()
    print(pyrolysis_sentence(rows))
    print()
    print("输入：")
    for note in INPUT_NOTES:
        print(f"- {note}")


if __name__ == "__main__":
    main()
