"""异相反应 R1–R4。速率是有符号进度；R1 的印刷负号已按进度约定去掉。"""
from __future__ import annotations

import math

from hamel.kinetics.arrhenius import arrhenius_rate_constant
from hamel.kinetics.common import (
    KineticsInputError,
    ReactionRate,
    _finite,
    _iteration_count,
    _particle_count,
    _positive,
    _real_power,
    _sphere_area_m2,
)

# 式 (5.20) 的参考压力，印刷为 101.300 Pa。不是 atm 的 101325 Pa。
_FILM_REFERENCE_PA = 101_300.0
# 牛顿法数值导数的相对步长 1e-6。求解器参数：ADR 0040。不是物性，也不改方程。
_NEWTON_RELATIVE_STEP = 1.0e-6
_LIGNITE_K0_M_S_K = 1.22
_LIGNITE_E_OVER_R_K = 10_300.0


def mechanism_factor(d_c_m: float, t_c_k: float) -> float:
    """式 (5.9)、(5.10) 的 Φ_C。

    d_c 与 T_c 按印刷符号传入。原文未说明它们是否等于 d_P、T_P。

    Source: Hamel Eq.5.9, S. 86
    Source: Hamel Eq.5.10, S. 86
    """
    diameter = _positive("d_c_m", d_c_m)
    temperature = _positive("t_c_k", t_c_k)
    p_c = arrhenius_rate_constant(2_500.0, temperature, activation_temperature_k=6_240.0)
    if diameter <= 0.05e-3:
        return (2.0 * p_c + 2.0) / (p_c + 2.0)
    if diameter <= 1.0e-3:
        removed = (p_c / 0.95) * (1.0e3 * diameter - 0.05)
        return (2.0 * p_c + 2.0 - removed) / (p_c + 2.0)
    return 1.0


def _film_diffusivity_m2_s(film_temperature_k: float, pressure_pa: float) -> float:
    """式 (5.20)。1,75 次只作用在温度比上。"""
    temperature_ratio = film_temperature_k / 1_500.0
    return 3.13e-4 * temperature_ratio**1.75 * (_FILM_REFERENCE_PA / pressure_pa)


def _mass_transfer_m_s(values, t_p_k: float) -> float:
    """式 (5.16)–(5.21)。式 (5.18) 的分母是运动粘度（ADR 0010 第 6 条）。"""
    film_temperature = 0.5 * (values["t_inf_k"] + t_p_k)
    diffusivity = _film_diffusivity_m2_s(film_temperature, values["pressure_pa"])
    kinematic = values["mu_g_pa_s"] / values["rho_g_kg_m3"]
    reynolds = values["u_r_m_s"] * values["d_p_m"] / kinematic
    schmidt = values["mu_g_pa_s"] / (values["rho_g_kg_m3"] * diffusivity)
    sherwood = 2.0 + 0.6 * _real_power("Re", reynolds, 0.5) * _real_power("Sc", schmidt, 0.33)
    return sherwood * diffusivity / values["d_p_m"]


def _r1_progress_per_particle(values, d_c_m, t_c_k, t_r1_k, t_p_k):
    """单颗粒进度 mol/s，以及外表面积。印刷式 (5.14) 的前导负号已去掉。"""
    phi = mechanism_factor(d_c_m, t_c_k)
    k_chemical = arrhenius_rate_constant(
        _LIGNITE_K0_M_S_K * t_r1_k,
        t_r1_k,
        activation_temperature_k=_LIGNITE_E_OVER_R_K,
    )
    k_diffusion = _mass_transfer_m_s(values, t_p_k)
    resistance = 1.0 / (phi * k_diffusion) + 1.0 / k_chemical
    if resistance == 0.0:
        raise KineticsInputError("式 (5.15) 的阻力为 0")
    area = _sphere_area_m2(values["d_p_m"])
    progress = area * (1.0 / resistance) * values["c_o2_mol_m3"]
    return progress, area


def _r1_numbers(
    d_k_m: float,
    t_inf_k: float,
    pressure_pa: float,
    u_r_m_s: float,
    mu_g_pa_s: float,
    rho_g_kg_m3: float,
    c_o2_mol_m3: float,
) -> dict[str, float]:
    return {
        "d_p_m": _positive("d_k_m", d_k_m),
        "t_inf_k": _positive("t_inf_k", t_inf_k),
        "pressure_pa": _positive("pressure_pa", pressure_pa),
        "u_r_m_s": _finite("u_r_m_s", u_r_m_s),
        "mu_g_pa_s": _positive("mu_g_pa_s", mu_g_pa_s),
        "rho_g_kg_m3": _positive("rho_g_kg_m3", rho_g_kg_m3),
        "c_o2_mol_m3": _finite("c_o2_mol_m3", c_o2_mol_m3),
    }


def r1_combustion_rate(
    m_char_kg, rho_p_kg_m3, d_k_m, d_c_m, t_c_k, t_r1_k, t_p_k,
    t_inf_k, c_o2_mol_m3, pressure_pa, u_r_m_s, mu_g_pa_s, rho_g_kg_m3,
) -> ReactionRate:
    """R1 单颗粒进度乘颗粒数。莱茵褐煤用表 5.2 的 Lignite A。

    N_k = m_C,k / (ρ_P · π d_k³ / 6)：反应颗粒数按焦炭质量取，ρ_P 作焦炭的表观密度。
    ``t_r1_k`` 是式 (5.23) 的 T，``d_c_m`` 与 ``t_c_k`` 是式 (5.9)、(5.10) 的符号。
    原文未说明它们是否等于颗粒温度和 d_P，这里不替换。

    解释性选择：ADR 0044

    Source: Hamel Eq.5.14, S. 86
    Source: Hamel Eq.5.15, S. 87
    Source: Hamel Eq.5.23, S. 88
    Source: Hamel Tab.5.2, S. 88
    """
    particle_mass = _finite("m_char_kg", m_char_kg)
    density = _positive("rho_p_kg_m3", rho_p_kg_m3)
    values = _r1_numbers(
        d_k_m, t_inf_k, pressure_pa, u_r_m_s, mu_g_pa_s, rho_g_kg_m3, c_o2_mol_m3
    )
    if values["u_r_m_s"] < 0.0:
        raise KineticsInputError("u_r_m_s 为负时雷诺数的平方根没有实值")
    _positive("t_r1_k", t_r1_k)
    _positive("t_p_k", t_p_k)
    progress, _area = _r1_progress_per_particle(values, d_c_m, t_c_k, t_r1_k, t_p_k)
    count = _particle_count(particle_mass, density, values["d_p_m"])
    return ReactionRate(specific=progress, extent_mol_s=progress * count)


def particle_temperature(
    d_k_m, d_c_m, t_c_k, t_r1_k, t_inf_k, c_o2_mol_m3, pressure_pa, u_r_m_s,
    mu_g_pa_s, rho_g_kg_m3, delta_h_j_mol, alpha_p_g_w_m2_k, epsilon_p,
    sigma_w_m2_k4, t_p_guess_k, tolerance_w, max_iterations,
) -> float:
    """由式 (5.22) 求 T_P。α_P,g、ε_P、σ 不给默认值，取值见 ADR 0017。

    容差是式 (5.22) 左右两边之差的绝对值，单位 W。牛顿步若离开 T>0，
    只把步长减半以留在定义域，不改方程，也不把温度截到某个下限。
    式 (5.23) 的 T 与式 (5.10) 的 T_c 在迭代中保持传入值；随 T_P 变化的是式 (5.21) 的膜温。

    Source: Hamel Eq.5.22, S. 87
    """
    values = _r1_numbers(
        d_k_m, t_inf_k, pressure_pa, u_r_m_s, mu_g_pa_s, rho_g_kg_m3, c_o2_mol_m3
    )
    if values["u_r_m_s"] < 0.0:
        raise KineticsInputError("u_r_m_s 为负时雷诺数的平方根没有实值")
    _positive("t_r1_k", t_r1_k)
    enthalpy = _finite("delta_h_j_mol", delta_h_j_mol)
    alpha = _finite("alpha_p_g_w_m2_k", alpha_p_g_w_m2_k)
    emissivity = _finite("epsilon_p", epsilon_p)
    sigma = _finite("sigma_w_m2_k4", sigma_w_m2_k4)
    tolerance = _finite("tolerance_w", tolerance_w)
    if tolerance < 0.0:
        raise KineticsInputError("tolerance_w 不能为负")
    rounds = _iteration_count("max_iterations", max_iterations)
    temperature = _positive("t_p_guess_k", t_p_guess_k)

    def residual(particle_k: float) -> float:
        progress, area = _r1_progress_per_particle(values, d_c_m, t_c_k, t_r1_k, particle_k)
        # 印刷式 R_C 为负的碳消耗速率，R_C = −进度。
        left = -progress * enthalpy
        convection = area * alpha * (particle_k - values["t_inf_k"])
        radiation = area * sigma * emissivity * (particle_k**4 - values["t_inf_k"] ** 4)
        return left - convection - radiation

    for _ in range(rounds):
        mismatch = residual(temperature)
        # 收敛判据是能量残差的大小，不作用在反应速率上。
        if abs(mismatch) <= tolerance:
            return temperature
        step_k = _NEWTON_RELATIVE_STEP * temperature
        above = residual(temperature + step_k)
        if temperature - step_k <= 0.0:
            derivative = (above - mismatch) / step_k
        else:
            below = residual(temperature - step_k)
            derivative = (above - below) / (2.0 * step_k)
        if derivative == 0.0 or not math.isfinite(derivative):
            raise KineticsInputError("式 (5.22) 的牛顿步没有有限导数")
        step = mismatch / derivative
        updated = temperature - step
        shrinks = 0
        while updated <= 0.0 and shrinks < 40:
            step *= 0.5
            updated = temperature - step
            shrinks += 1
        if updated <= 0.0 or not math.isfinite(updated):
            raise KineticsInputError("式 (5.22) 的迭代离开 T>0")
        temperature = updated
    if abs(residual(temperature)) <= tolerance:
        return temperature
    raise KineticsInputError(f"式 (5.22) 在 {rounds} 轮内残差未降到容差以内")


def _surface_progress(
    m_char_kg, rho_p_kg_m3, diameter_m, temperature_k,
    concentration_mol_m3, prefactor_m_s_k, activation_temperature_k,
) -> ReactionRate:
    """R2、R4：单颗粒 A_C·k·C，再乘按焦炭质量计的颗粒数。A_C 取颗粒外表面。"""
    mass = _finite("m_char_kg", m_char_kg)
    density = _positive("rho_p_kg_m3", rho_p_kg_m3)
    diameter = _positive("d_k_m", diameter_m)
    temperature = _positive("temperature_k", temperature_k)
    concentration = _finite("concentration_mol_m3", concentration_mol_m3)
    coefficient = arrhenius_rate_constant(
        prefactor_m_s_k * temperature,
        temperature,
        activation_temperature_k=activation_temperature_k,
    )
    per_particle = _sphere_area_m2(diameter) * coefficient * concentration
    count = _particle_count(mass, density, diameter)
    return ReactionRate(specific=per_particle, extent_mol_s=per_particle * count)


def r2_steam_gasification_rate(
    m_char_kg, rho_p_kg_m3, d_k_m, c_h2o_mol_m3, temperature_k,
) -> ReactionRate:
    """R2。式 (5.25) 的 T 原文未说明是哪一个温度，由调用方传入。

    颗粒数同 R1，按焦炭质量与 ρ_P 计（ADR 0044）。
    冻结稿记 A_C 的几何意义原文未写。式面是面积乘浓度，不是按粒级质量给出的比速率。

    解释性选择：ADR 0044

    Source: Hamel Eq.5.24, S. 88
    Source: Hamel Eq.5.25, S. 88
    """
    return _surface_progress(
        m_char_kg, rho_p_kg_m3, d_k_m, temperature_k, c_h2o_mol_m3, 3.42, 15_600.0
    )


def r3_boudouard_rate(
    m_char_kg: float,
    molar_mass_c_kg_mol: float,
    p_co2_pa: float,
    p_co_pa: float,
    temperature_k: float,
    gas_constant_j_mol_k: float,
) -> ReactionRate:
    """R3。比速率为 1/s，cell 速率只乘 m_j,k/M_C 一次。

    只乘焦炭质量，不乘颗粒内灰，也不用粒径和 ρ_P。T 取 cell 温度，R_g 由调用方传入。
    解释性选择：ADR 0010 第 8 条。负活化能保留。

    Source: Hamel Eq.5.26, S. 89
    Source: Hamel Eq.5.33, S. 89
    """
    mass = _finite("m_char_kg", m_char_kg)
    molar_mass = _positive("molar_mass_c_kg_mol", molar_mass_c_kg_mol)
    pressure_co2 = _finite("p_co2_pa", p_co2_pa)
    pressure_co = _finite("p_co_pa", p_co_pa)
    temperature = _positive("temperature_k", temperature_k)
    gas_constant = _positive("gas_constant_j_mol_k", gas_constant_j_mol_k)
    k_a = arrhenius_rate_constant(
        2.72, temperature, activation_energy_j_mol=188_000.0, gas_constant_j_mol_k=gas_constant
    )
    k_b = arrhenius_rate_constant(
        8.31e-7, temperature, activation_energy_j_mol=-38_900.0, gas_constant_j_mol_k=gas_constant
    )
    k_c = arrhenius_rate_constant(
        2.088e-14,
        temperature,
        activation_energy_j_mol=-24_500.0,
        gas_constant_j_mol_k=gas_constant,
    )
    denominator = 1.0 + k_b * pressure_co2 + k_c * pressure_co
    if denominator == 0.0:
        raise KineticsInputError("式 (5.26) 的分母为 0")
    specific = k_a * pressure_co2 / denominator
    return ReactionRate(specific=specific, extent_mol_s=specific * mass / molar_mass)


def r4_hydrogasification_rate(
    m_char_kg, rho_p_kg_m3, d_k_m, c_h2_mol_m3, temperature_k,
) -> ReactionRate:
    """R4。对 H₂ 是一次，不是计量系数 2 的二次。颗粒数同 R2。

    解释性选择：ADR 0044

    Source: Hamel Eq.5.34, S. 89
    Source: Hamel Eq.5.35, S. 89
    """
    return _surface_progress(
        m_char_kg, rho_p_kg_m3, d_k_m, temperature_k, c_h2_mol_m3, 3.42e-3, 15_600.0
    )
