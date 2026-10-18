#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单颗粒干燥与热解的平均释放时间预检。不进入 hamel/。

T_w = ∫ Y_w dt，T_v = ∫ Y_v dt（ADR 0015 第 10 项第 4 步）。
Y_w 是剩余水分份额，Y_v = 1 − V/V_∞ 是剩余挥发分份额。
有效域是 T_v ≥ T_w。这里只报告是否成立，不把 τ_pyr 改成 0。

式 (4.3) 取负号，式 (4.11) 的指数取负号，式 (4.12) 的前置因子取 3/R₀³。
这三处是 ADR 0015 第 5–7 条的读法。湿核内不热解：§4.1 写明热解只发生在已干壳层。

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


def shell_temperature(radius_m: np.ndarray, front_m: float, outer_m: float, core_k: float, surface_k: float):
    """准稳态干壳温度。r ≤ r_e 保持 T_e。Source: Hamel Eq.4.5, S. 53"""
    temperature = np.full(radius_m.shape, core_k)
    dry = radius_m > front_m
    if not np.any(dry):
        return temperature
    radial = radius_m[dry]
    factor = (1.0 - front_m / radial) / (1.0 - front_m / outer_m)
    temperature[dry] = core_k + (surface_k - core_k) * factor
    return temperature


def surface_temperature(phi: float, particle: ParticleProperties, radius_m: float) -> float:
    alpha = heat_transfer_coefficient(particle, 2.0 * radius_m)
    biot = alpha * radius_m / particle.conductivity
    if phi >= 1.0:
        return particle.evaporation_temperature_k
    ratio = biot * (1.0 - phi) / phi
    core = particle.evaporation_temperature_k
    ambient = particle.temperature_k
    return (core + ratio * ambient) / (1.0 + ratio)


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
    water = [1.0]
    volatile = [1.0]
    phi = 1.0
    elapsed = 0.0
    # 锋面从表面推进到接近中心。φ³ 在 10⁻³ 时已是 10⁻⁹，后面的水分积分可以忽略。
    for _ in range(20_000):
        if phi <= 1.0e-3:
            break
        rate = front_rate(phi, particle, radius_m)
        if rate >= 0.0:
            raise RuntimeError("干燥锋面没有向内移动")
        step = min(0.05, max(1.0e-6, 0.01 / abs(rate)))
        next_phi = phi + step * rate
        if next_phi < 1.0e-3:
            step = (1.0e-3 - phi) / rate
            next_phi = 1.0e-3
        front = phi * radius_m
        surface = surface_temperature(phi, particle, radius_m)
        temperature = shell_temperature(
            radius, front, radius_m, particle.evaporation_temperature_k, surface
        )
        dry = radius >= front
        exposure[dry] += np.exp(
            -energies / (gas_constant * np.maximum(temperature[dry, None], 250.0))
        ) * step
        phi = next_phi
        elapsed += step
        profile = np.ones(radial_count)
        profile[dry] = _remaining(exposure[dry], particle.frequency_s, weights)
        times.append(elapsed)
        water.append(phi**3)
        volatile.append(average_volatile(radius, profile, radius_m))
    else:
        return _unresolved(diameter_m, "干燥锋面在步数上限内没有推进到中心")

    alpha = heat_transfer_coefficient(particle, diameter_m)
    initial = shell_temperature(
        radius,
        phi * radius_m,
        radius_m,
        particle.evaporation_temperature_k,
        surface_temperature(max(phi, 1.0e-3), particle, radius_m),
    )

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
    "活化能用 12 个 Gauss–Hermite 节点。径向 21 个节点。干燥推进到 φ = 10⁻³ 后改积导热方程。临时离散。",
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
