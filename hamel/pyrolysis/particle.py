"""单颗粒干燥与热解的平均释放时间。"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import simpson, solve_ivp

from hamel.cases.config import SolverConfig
from hamel.kinetics.arrhenius import arrhenius_rate_constant_array
from hamel.properties.data import gas_constant_j_mol_k, molar_mass_kg_mol, standard_pressure_pa
from hamel.properties.thermo import (
    chemical_potential_j_mol, density_kg_m3, enthalpy_j_mol, heat_capacity_j_mol_k,
    vaporization_enthalpy_j_mol,
)
from hamel.properties.transport import mixture_thermal_conductivity_w_m_k, mixture_viscosity_pa_s

class PyrolysisError(ValueError):
    """热解输入或积分不能按式 (4.1)–(4.12) 继续。"""

class PyrolysisDomainError(PyrolysisError):
    """T_v < T_w。解释性选择：ADR 0015 第 10 项。"""

@dataclass(frozen=True)
class ParticleProperties:
    """褐煤颗粒物性。表 4.1 的一组由 lignite_particle_properties 给出。"""

    activation_mean_j_mol: float
    activation_sigma_j_mol: float
    frequency_s: float
    density_kg_m3: float
    heat_capacity_j_kg_k: float
    diffusivity_m2_s: float

@dataclass(frozen=True)
class ReleaseTimes:
    """τ_dry = T_w，τ_pyr = T_v − T_w。tail_error_bound_s 是温度阈值造成的上界。"""

    t_w_s: float
    t_v_s: float
    tau_dry_s: float
    tau_pyr_s: float
    tail_error_bound_s: float

def lignite_particle_properties() -> ParticleProperties:
    """湿褐煤，Dersch 1994。k_0、ρ_s、c_s 用表内 k_0,P、ρ_P、c_P。Source: Hamel Tab.4.1, S. 54"""
    return ParticleProperties(192_000.0, 40_000.0, 1.67e13, 1250.0, 1256.0, 1.0e-7)

def pyrolysis_duration_s(water_time_s: float, volatile_time_s: float) -> float:
    """τ_pyr = T_v − T_w。T_v < T_w 超出两级串联的有效域。Source: Hamel Eq.4.10, S. 53"""
    if volatile_time_s < water_time_s:
        raise PyrolysisDomainError(
            f"T_v = {volatile_time_s} s < T_w = {water_time_s} s，超出两级串联的有效域"
        )
    return volatile_time_s - water_time_s

def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise PyrolysisError(f"{name} 必须是有限实数")
    number = float(value)
    if not math.isfinite(number):
        raise PyrolysisError(f"{name} 不是有限实数")
    return number

def _positive(name: str, value: object) -> float:
    number = _finite(name, value)
    if number <= 0.0:
        raise PyrolysisError(f"{name} 必须为正")
    return number


def conductivity_w_m_k(material: ParticleProperties) -> float:
    """λ_s = a ρ c。a 是式 (4.1) 的导温系数。物理定义。Source: Hamel Eq.4.1, S. 52"""
    return (
        _positive("diffusivity_m2_s", material.diffusivity_m2_s)
        * _positive("density_kg_m3", material.density_kg_m3)
        * _positive("heat_capacity_j_kg_k", material.heat_capacity_j_kg_k)
    )


def evaporation_temperature_k(pressure_pa: float) -> float:
    """床压下 μ(H2O,L) = μ(H2O,g)。凝相不作压力修正。Source: Hamel Eq.4.4, S. 52"""
    pressure, gas_constant, reference = _positive("pressure_pa", pressure_pa), gas_constant_j_mol_k(), standard_pressure_pa()

    def residual(temperature: float) -> float:
        return chemical_potential_j_mol("H2O_L", temperature) - (
            chemical_potential_j_mol("H2O", temperature) + gas_constant * temperature * math.log(pressure / reference)
        )

    low, high = 280.0, 590.0
    if residual(low) > 0.0 or residual(high) < 0.0:
        raise PyrolysisError("蒸发温度不在液态水多项式的可用区间内")
    for _ in range(60):
        mid = 0.5 * (low + high)
        low, high = (mid, high) if residual(mid) < 0.0 else (low, mid)
    return 0.5 * (low + high)


def corrected_vaporization_j_kg(moisture_mass_fraction: float, evaporation_k: float, initial_k: float, heat_capacity_j_kg_k: float) -> tuple[float, float]:
    """式 (4.4)。w_0,tr = x_w/(1−x_w)，x_w 为收到基。c_w(T_e−T_0) 取液态水焓差。Source: Hamel Eq.4.4, S. 52"""
    moisture = _finite("moisture_mass_fraction", moisture_mass_fraction)
    if moisture <= 0.0 or moisture >= 1.0:
        raise PyrolysisError("收到基水分必须在 (0, 1) 内才能计算修正蒸发焓")
    dry_basis, molar = moisture / (1.0 - moisture), molar_mass_kg_mol("H2O")
    latent = vaporization_enthalpy_j_mol(evaporation_k) / molar
    sensible_water = (enthalpy_j_mol("H2O_L", evaporation_k) - enthalpy_j_mol("H2O_L", initial_k)) / molar
    sensible_solid = _positive("heat_capacity_j_kg_k", heat_capacity_j_kg_k) * (evaporation_k - initial_k)
    return latent + sensible_water + sensible_solid / dry_basis, dry_basis


def _mixture_heat_capacity_j_kg_k(mole_fractions: dict[str, float], temperature_k: float) -> float:
    total = mass = 0.0
    for name, fraction in mole_fractions.items():
        share = float(fraction)
        if share < 0.0:
            raise PyrolysisError(f"{name} 的摩尔分数为负")
        total += share
        mass += share * molar_mass_kg_mol(name)
    if total <= 0.0:
        raise PyrolysisError("气体摩尔分数之和必须为正")
    capacity = sum(float(fraction) * heat_capacity_j_mol_k(name, temperature_k) for name, fraction in mole_fractions.items() if fraction)
    return capacity / mass


def heat_transfer_coefficient(diameter_m: float, temperature_k: float, surface_k: float, pressure_pa: float, mole_fractions: dict[str, float], relative_velocity_m_s: float, solid_conductivity: float) -> tuple[float, float]:
    """Nu 用式 (4.9)，Bi 用式 (4.8)。Re = u_r d/ν，Pr = c_p μ/λ，物性取膜温。

    u_r 与气体组成由调用方传入。解释性选择：ADR 0017 第 20 项。
    Source: Hamel Eq.4.8, S. 53
    Source: Hamel Eq.4.9, S. 53
    """
    diameter = _positive("diameter_m", diameter_m)
    velocity = _finite("relative_velocity_m_s", relative_velocity_m_s)
    if velocity < 0.0:
        raise PyrolysisError("相对速度为负时雷诺数的平方根没有实值")
    film = 0.5 * (_positive("temperature_k", temperature_k) + _positive("surface_k", surface_k))
    pressure = _positive("pressure_pa", pressure_pa)
    viscosity = mixture_viscosity_pa_s(mole_fractions, film, pressure)
    conductivity = mixture_thermal_conductivity_w_m_k(mole_fractions, film, pressure)
    prandtl = _mixture_heat_capacity_j_kg_k(mole_fractions, film) * viscosity / conductivity
    reynolds = density_kg_m3(film, pressure, mole_fractions) * velocity * diameter / viscosity
    nusselt = 2.0 + 1.2 * math.sqrt(reynolds) * prandtl ** (1.0 / 3.0)
    alpha = nusselt * conductivity / diameter
    return alpha, alpha * (0.5 * diameter) / _positive("solid_conductivity", solid_conductivity)


def remaining_time_integral(times, remaining, tail: float = 0.0) -> float:
    """T = ∫ Y dt。tail 是最后一个采样之后到无穷的积分。Source: Hamel Eq.4.10, S. 53"""
    grid, series = np.asarray(times, dtype=float), np.asarray(remaining, dtype=float)
    if grid.ndim != 1 or series.shape != grid.shape or grid.size < 2 or np.any(np.diff(grid) <= 0.0):
        raise PyrolysisError("释放时间至少要两个严格递增的采样")
    return float(np.trapz(series, grid) + tail)


def _activation_nodes(material: ParticleProperties, count: int):
    if count < 2:
        raise PyrolysisError("活化能节点至少为 2")
    nodes, weights = np.polynomial.hermite_e.hermegauss(count)
    energies = material.activation_mean_j_mol + material.activation_sigma_j_mol * nodes
    return energies, weights / math.sqrt(2.0 * math.pi)


def _volume_average(profile, radius: np.ndarray) -> float:
    """前置因子 3/R₀³。解释性选择：ADR 0015 第 7 项。Source: Hamel Eq.4.12, S. 54"""
    outer = float(radius[-1])
    return 3.0 * float(simpson(np.asarray(profile, dtype=float) * radius**2, x=radius)) / outer**3


def _use_config(config: SolverConfig) -> SolverConfig:
    intervals, nodes = config.pyrolysis_radial_intervals, config.pyrolysis_activation_nodes
    if isinstance(intervals, bool) or not isinstance(intervals, (int, np.integer)) or int(intervals) < 4 or int(intervals) % 2:
        raise PyrolysisError("径向区间数须为不小于 4 的偶数")
    _positive("pyrolysis_time_rtol", config.pyrolysis_time_rtol)
    if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)) or int(nodes) < 2:
        raise PyrolysisError("活化能节点数须为不小于 2 的整数")
    _positive("pyrolysis_tail_delta_t_k", config.pyrolysis_tail_delta_t_k)
    if _positive("pyrolysis_initial_shell_fraction", config.pyrolysis_initial_shell_fraction) >= 0.1:
        raise PyrolysisError("初始干壳份额必须小于 0.1")
    return config


def _shell_rhs(state, alpha, evaporation, bed, conductivity, diffusivity, capacity, outer):
    """ξ=(r−r_e)/(R−r_e) 上的式 (4.1)。式 (4.3) 取 −dr_e/dt。ADR 0015 第 5 项。"""
    front, thickness, intervals = float(state[0]), outer - float(state[0]), state.size - 1
    out = np.zeros_like(state)
    if thickness <= 0.0:
        return out
    values = np.empty(intervals + 1)
    values[0], values[1:] = evaporation, state[1:]
    step = 1.0 / intervals
    # Bi = α R / λ，∂T/∂ξ = (δ/R) Bi (T_a − T_s)。式 (4.2)、(4.8)。
    surface_slope = (thickness / outer) * (alpha * outer / conductivity) * (bed - values[-1])
    ghost = values[-2] + 2.0 * step * surface_slope
    front_slope = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (2.0 * step)
    speed = -(conductivity / capacity) * front_slope / thickness
    for index in range(1, intervals + 1):
        right = ghost if index == intervals else values[index + 1]
        dxi = (right - values[index - 1]) / (2.0 * step)
        dxi2 = (right - 2.0 * values[index] + values[index - 1]) / step**2
        xi = index * step
        conduction = diffusivity * (dxi2 / thickness**2 + 2.0 * dxi / ((front + xi * thickness) * thickness))
        out[index] = conduction + speed * (1.0 - xi) / thickness * dxi
    out[0] = speed
    return out


def _heat_rhs(temperature, alpha, bed, conductivity, diffusivity, outer):
    count, step = temperature.size, outer / (temperature.size - 1)
    rate = np.empty(count)
    rate[0] = 6.0 * diffusivity * (temperature[1] - temperature[0]) / step**2
    for index in range(1, count - 1):
        second = (temperature[index + 1] - 2.0 * temperature[index] + temperature[index - 1]) / step**2
        first = (temperature[index + 1] - temperature[index - 1]) / (2.0 * step)
        rate[index] = diffusivity * (second + 2.0 * first / (index * step))
    slope = (alpha / conductivity) * (bed - temperature[-1])
    ghost = temperature[-2] + 2.0 * step * slope
    second = (ghost - 2.0 * temperature[-1] + temperature[-2]) / step**2
    first = (ghost - temperature[-2]) / (2.0 * step)
    rate[-1] = diffusivity * (second + 2.0 * first / outer)
    return rate


def _map_shell(state, evaporation, radius, outer):
    front, intervals = float(state[0]), state.size - 1
    values = np.empty(intervals + 1)
    values[0], values[1:] = evaporation, state[1:]
    temperature = np.full(radius.shape, evaporation)
    dry = radius >= front
    thickness = outer - front
    if thickness > 0.0 and np.any(dry):
        temperature[dry] = np.interp((radius[dry] - front) / thickness, np.linspace(0.0, 1.0, values.size), values)
    # 湿核保持 T_e，不积累热解暴露。§4.1：热解只发生在已干壳层。
    return temperature, dry


def _exposure_rate(temperature, dry, energies, gas_constant):
    rate = np.zeros((temperature.size, energies.size))
    if np.any(dry):
        rate[dry] = arrhenius_rate_constant_array(1.0, temperature[dry, None], energies, gas_constant)
    return rate


def _remaining(theta, frequency, weights, radius) -> float:
    return _volume_average(np.exp(-frequency * theta) @ weights, radius)


def _march(theta, solved, temperature_of, radius, energies, gas_constant, frequency, weights, offset):
    """辛普森积分暴露量。Y 只记在求解器节点上，时间积分再对这些节点做梯形。"""
    times = [offset + float(solved.t[0])]
    values = [_remaining(theta, frequency, weights, radius)]
    for left, right in zip(solved.t[:-1], solved.t[1:]):
        rates = []
        for instant in (float(left), 0.5 * (float(left) + float(right)), float(right)):
            temperature, dry = temperature_of(solved.sol(instant))
            rates.append(_exposure_rate(temperature, dry, energies, gas_constant))
        theta = theta + (float(right) - float(left)) / 6.0 * (rates[0] + 4.0 * rates[1] + rates[2])
        times.append(offset + float(right))
        values.append(_remaining(theta, frequency, weights, radius))
    return theta, times, values


def _integrate(rhs, initial, horizon, event, rtol, atol):
    solved = solve_ivp(rhs, (0.0, horizon), initial, method="BDF", dense_output=True, events=event, rtol=rtol, atol=atol)
    if (not solved.success) or solved.sol is None:
        raise PyrolysisError(f"时间积分失败：{solved.message}")
    return solved


def drying_energy_relative_gap(diameter_m: float, temperature_k: float, pressure_pa: float, moisture_mass_fraction: float, initial_temperature_k: float, mole_fractions: dict[str, float], relative_velocity_m_s: float, material: ParticleProperties, config: SolverConfig) -> float:
    """同步阶段：壳层显热增量 + 锋面潜热，对照表面传入的热。显热以 T_e 为参考。

    Source: Hamel Eq.4.1, S. 52
    Source: Hamel Eq.4.3, S. 52
    Source: Hamel Eq.4.4, S. 52
    """
    solved, radius, evaporation, corrected, dry_basis, density, heat_capacity, alpha_of = _dry_particle(
        diameter_m, temperature_k, pressure_pa, moisture_mass_fraction, initial_temperature_k,
        mole_fractions, relative_velocity_m_s, material, config,
    )
    outer, entered, previous_power, previous_time = float(radius[-1]), 0.0, None, None
    latent_factor = corrected * dry_basis * density * 4.0 * math.pi / 3.0
    front_start, gap, heat, stored0 = float(solved.y[0, 0]), 0.0, 0.0, None
    for instant in solved.t:
        state = solved.sol(float(instant))
        temperature, _dry = _map_shell(state, evaporation, radius, outer)
        stored = density * heat_capacity * 4.0 * math.pi * float(simpson((temperature - evaporation) * radius**2, x=radius))
        power = 4.0 * math.pi * outer**2 * alpha_of(float(temperature[-1])) * (temperature_k - float(temperature[-1]))
        if previous_time is not None:
            entered += 0.5 * (previous_power + power) * (float(instant) - previous_time)
        previous_power, previous_time = power, float(instant)
        if stored0 is None:
            stored0 = stored
        latent = latent_factor * (front_start**3 - float(state[0]) ** 3)
        heat, gap = max(heat, entered), max(gap, abs(stored - stored0 + latent - entered))
    if heat <= 0.0:
        raise PyrolysisError("干燥阶段没有热量传入")
    return gap / heat


def _dry_particle(diameter_m, temperature_k, pressure_pa, moisture_mass_fraction, initial_temperature_k, mole_fractions, relative_velocity_m_s, material, config):
    config = _use_config(config)
    outer, bed = 0.5 * _positive("diameter_m", diameter_m), _positive("temperature_k", temperature_k)
    initial = _positive("initial_temperature_k", initial_temperature_k)
    if initial > bed:
        raise PyrolysisError("颗粒初温高于床温时，尾积分的误差界不再成立")
    evaporation = evaporation_temperature_k(pressure_pa)
    if bed <= evaporation:
        raise PyrolysisError("床温不高于蒸发温度时，干燥锋面不会向内移动")
    corrected, dry_basis = corrected_vaporization_j_kg(moisture_mass_fraction, evaporation, initial, material.heat_capacity_j_kg_k)
    solid_conductivity = conductivity_w_m_k(material)

    def alpha_of(surface):
        return heat_transfer_coefficient(2.0 * outer, bed, surface, pressure_pa, mole_fractions, relative_velocity_m_s, solid_conductivity)[0]

    intervals, fraction = int(config.pyrolysis_radial_intervals), float(config.pyrolysis_initial_shell_fraction)
    capacity = corrected * dry_basis * material.density_kg_m3
    state = np.full(intervals + 1, evaporation)
    state[0] = outer * (1.0 - fraction)
    horizon = 40.0 * outer * capacity / (alpha_of(evaporation) * (bed - evaporation))

    def rhs(_time, current):
        return _shell_rhs(current, alpha_of(float(current[-1])), evaporation, bed, solid_conductivity, material.diffusivity_m2_s, capacity, outer)

    def core(_time, current):
        return float(current[0] - fraction * outer)

    core.terminal, core.direction = True, -1
    tolerances = np.full(state.size, config.pyrolysis_time_rtol * bed)
    tolerances[0] = config.pyrolysis_time_rtol * outer
    solved = _integrate(rhs, state, horizon, core, config.pyrolysis_time_rtol, tolerances)
    if solved.t_events[0].size == 0:
        raise PyrolysisError("干燥锋面没有到达颗粒中心")
    return solved, np.linspace(0.0, outer, intervals + 1), evaporation, corrected, dry_basis, material.density_kg_m3, material.heat_capacity_j_kg_k, alpha_of


def _tail(theta, frequency, weights, energies, gas_constant, bed, delta):
    coefficient = arrhenius_rate_constant_array(frequency, bed, energies, gas_constant)
    slower = arrhenius_rate_constant_array(frequency, bed - delta, energies, gas_constant)
    remain = np.exp(-frequency * theta)
    return remain @ (weights / coefficient), remain @ (weights * (1.0 / slower - 1.0 / coefficient))


def _within_band(temperature, bed, delta) -> bool:
    return float(np.min(temperature)) >= bed - delta and float(np.max(temperature)) <= bed + delta


def particle_release_times(diameter_m: float, temperature_k: float, pressure_pa: float, moisture_mass_fraction: float, initial_temperature_k: float, mole_fractions: dict[str, float], relative_velocity_m_s: float, ultimate_volatile_mass_fraction: float, material: ParticleProperties, config: SolverConfig) -> ReleaseTimes:
    """在 T_WS 下积分单颗粒。τ_dry = T_w，τ_pyr = T_v − T_w。环境温度取床温。r_0 与 R_0 都取 d/2。

    零水分时 Y_w ≡ 0；V_∞ = 0 时 Y_v ≡ 0。这两种情形不作除法。解释性选择：ADR 0015 第 10 项、ADR 0027。
    Source: Hamel Eq.4.1, S. 52
    Source: Hamel Eq.4.10, S. 53
    Source: Hamel Eq.4.12, S. 54
    """
    config = _use_config(config)
    moisture = _finite("moisture_mass_fraction", moisture_mass_fraction)
    volatile = _finite("ultimate_volatile_mass_fraction", ultimate_volatile_mass_fraction)
    if moisture < 0.0 or moisture >= 1.0 or volatile < 0.0:
        raise PyrolysisError("水分必须在 [0, 1) 内，最终挥发分不能为负")
    if volatile == 0.0 and moisture == 0.0:
        return ReleaseTimes(0.0, 0.0, 0.0, 0.0, 0.0)
    bed = _positive("temperature_k", temperature_k)
    initial = _positive("initial_temperature_k", initial_temperature_k)
    gas_constant, frequency = gas_constant_j_mol_k(), _positive("frequency_s", material.frequency_s)
    energies, weights = _activation_nodes(material, int(config.pyrolysis_activation_nodes))
    volatile_times, volatile_values = [0.0], [1.0]
    if moisture > 0.0:
        solved, radius, evaporation, corrected, dry_basis, _density, _heat, alpha_of = _dry_particle(
            diameter_m, bed, pressure_pa, moisture, initial, mole_fractions, relative_velocity_m_s, material, config,
        )
        outer, fraction = float(radius[-1]), float(config.pyrolysis_initial_shell_fraction)
        prelude = fraction * outer * material.density_kg_m3 * dry_basis * corrected / (alpha_of(evaporation) * (bed - evaporation))
        # 壳厚低于登记份额时表面仍在 T_e。水分从 1 降到 (1−ε)^3，挥发分保持 1。
        water_times, water_values = [0.0, prelude], [1.0, (1.0 - fraction) ** 3]
        for instant in solved.t[1:]:
            water_times.append(prelude + float(instant))
            water_values.append((float(solved.sol(float(instant))[0]) / outer) ** 3)
        water_time = remaining_time_integral(water_times, water_values)
        theta = np.zeros((radius.size, energies.size))
        theta, volatile_times, volatile_values = _march(
            theta, solved, lambda state: _map_shell(state, evaporation, radius, outer),
            radius, energies, gas_constant, frequency, weights, prelude,
        )
        volatile_times, volatile_values = [0.0, *volatile_times], [1.0, *volatile_values]
        temperature = _map_shell(solved.y[:, -1], evaporation, radius, outer)[0]
    else:
        water_time, outer = 0.0, 0.5 * _positive("diameter_m", diameter_m)
        radius = np.linspace(0.0, outer, int(config.pyrolysis_radial_intervals) + 1)
        temperature, theta = np.full(radius.shape, initial), np.zeros((radius.size, energies.size))
    if volatile == 0.0:
        return ReleaseTimes(water_time, 0.0, water_time, 0.0, 0.0)
    delta = float(config.pyrolysis_tail_delta_t_k)
    if delta >= bed:
        raise PyrolysisError("尾积分温度阈值必须小于床温")
    solid_conductivity = conductivity_w_m_k(material)

    def alpha_of(surface: float) -> float:
        return heat_transfer_coefficient(2.0 * outer, bed, surface, pressure_pa, mole_fractions, relative_velocity_m_s, solid_conductivity)[0]

    if not _within_band(temperature, bed, delta):
        def rhs(_time, current):
            return _heat_rhs(current, alpha_of(float(current[-1])), bed, solid_conductivity, material.diffusivity_m2_s, outer)

        def hot(_time, current):
            return float(np.min(current) - (bed - delta))

        hot.terminal, hot.direction = True, 1
        fourier = outer**2 / material.diffusivity_m2_s
        lumped = material.density_kg_m3 * material.heat_capacity_j_kg_k * outer / (3.0 * alpha_of(float(temperature[-1])))
        heated = _integrate(rhs, temperature, 80.0 * max(fourier, lumped), hot, config.pyrolysis_time_rtol, config.pyrolysis_time_rtol * bed)
        if heated.t_events[0].size == 0:
            raise PyrolysisError("颗粒没有进入尾积分的温度带")
        theta, heat_times, heat_values = _march(
            theta, heated, lambda state: (np.asarray(state, dtype=float), np.ones(radius.shape, dtype=bool)),
            radius, energies, gas_constant, frequency, weights, volatile_times[-1],
        )
        volatile_times.extend(heat_times[1:])
        volatile_values.extend(heat_values[1:])
    mean, bound = _tail(theta, frequency, weights, energies, gas_constant, bed, delta)
    history = 0.0 if len(volatile_times) < 2 else remaining_time_integral(volatile_times, volatile_values)
    volatile_time, bound_time = history + _volume_average(mean, radius), _volume_average(bound, radius)
    return ReleaseTimes(water_time, volatile_time, water_time, pyrolysis_duration_s(water_time, volatile_time), bound_time)
