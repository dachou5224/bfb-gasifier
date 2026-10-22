"""均相反应 R5–R11。进度速率的正值表示反应式正向。"""
from __future__ import annotations

import math

from hamel.kinetics.arrhenius import arrhenius_rate_constant
from hamel.kinetics.common import KineticsInputError, ReactionRate, _finite, _positive, _real_power

# 解释性选择：ADR 0010 第 2 条。式 (5.45) 的 p 以 atm 计。
_ATM_PA = 101_325.0
# 解释性选择（外推）：ADR 0013。Wen und Chaung 气流床所用，不是另测的流化床值。
_A_R8 = 0.02
# 式 (5.45) 第二个指数里与温度无关的项。
_R8_SHIFT = -8.91
# 解释性选择：ADR 0010 第 9 条。740 °C 归高温段。表 5.4 的区间端点含在内。
_CORELLA_LOW_K = 600.0 + 273.15
_CORELLA_SPLIT_K = 740.0 + 273.15
_CORELLA_HIGH_K = 900.0 + 273.15


def _extent(specific: float, factor: float) -> ReactionRate:
    return ReactionRate(specific=specific, extent_mol_s=specific * factor)


def r5_co_oxidation_rate(
    phase: str,
    c_co_mol_m3: float,
    c_o2_mol_m3: float,
    temperature_k: float,
    volume_m3: float,
    c_h2o_mol_m3: float | None = None,
) -> ReactionRate:
    """R5。悬浮相的短横线表示没有该因子：不乘 exp，也不乘水的浓度。

    解释性选择：ADR 0010 第 5 条。气泡相必须传入水的浓度。
    悬浮相即使传入 ``c_h2o_mol_m3`` 也不使用。

    Source: Hamel Eq.5.36, S. 90
    Source: Hamel Tab.5.3, S. 91
    """
    if phase not in ("bubble", "suspension"):
        raise KineticsInputError("phase 只能是 bubble 或 suspension")
    carbon_monoxide = _finite("c_co_mol_m3", c_co_mol_m3)
    oxygen = _real_power("c_o2_mol_m3", c_o2_mol_m3, 0.5 if phase == "suspension" else 0.3)
    temperature = _positive("temperature_k", temperature_k)
    volume = _finite("volume_m3", volume_m3)
    if phase == "bubble":
        if c_h2o_mol_m3 is None:
            raise KineticsInputError("气泡相 R5 必须传入 c_h2o_mol_m3")
        water = _real_power("c_h2o_mol_m3", c_h2o_mol_m3, 0.5)
        coefficient = arrhenius_rate_constant(
            1.91e6, temperature, activation_temperature_k=8_056.0
        )
        specific = coefficient * carbon_monoxide * oxygen * water
    else:
        specific = 1.8 * temperature**0.5 * carbon_monoxide * oxygen
    return _extent(specific, volume)


def r6_hydrogen_oxidation_rate(
    c_h2_mol_m3: float,
    c_o2_mol_m3: float,
    temperature_k: float,
    volume_m3: float,
) -> ReactionRate:
    """R6。印刷式没有前导负号，正的乘积即为进度速率。

    Source: Hamel Eq.5.38, S. 91
    Source: Hamel Eq.5.39, S. 91
    """
    hydrogen = _real_power("c_h2_mol_m3", c_h2_mol_m3, 1.5)
    oxygen = _finite("c_o2_mol_m3", c_o2_mol_m3)
    temperature = _positive("temperature_k", temperature_k)
    volume = _finite("volume_m3", volume_m3)
    coefficient = arrhenius_rate_constant(
        1.631e9 * temperature**1.5,
        temperature,
        activation_temperature_k=3_420.0,
    )
    return _extent(coefficient * hydrogen * oxygen, volume)


def r7_methane_oxidation_rate(
    c_ch4_mol_m3: float,
    c_o2_mol_m3: float,
    temperature_k: float,
    volume_m3: float,
) -> ReactionRate:
    """R7。两浓度都是一次。正的乘积即为进度速率。

    Source: Hamel Eq.5.40, S. 91
    Source: Hamel Eq.5.41, S. 91
    """
    methane = _finite("c_ch4_mol_m3", c_ch4_mol_m3)
    oxygen = _finite("c_o2_mol_m3", c_o2_mol_m3)
    temperature = _positive("temperature_k", temperature_k)
    volume = _finite("volume_m3", volume_m3)
    coefficient = arrhenius_rate_constant(
        3.552e11 / temperature,
        temperature,
        activation_temperature_k=15_700.0,
    )
    return _extent(coefficient * methane * oxygen, volume)


def r8_equilibrium_constant(
    mu_co_j_mol: float,
    mu_h2o_j_mol: float,
    mu_h2_j_mol: float,
    mu_co2_j_mol: float,
    temperature_k: float,
    gas_constant_j_mol_k: float,
) -> float:
    """式 (5.43)。化学势是合成或物性模块的输入，本函数不带数据。

    计量按 CO + H₂O → H₂ + CO₂：Σν_i μ_i⁰ = −μ_CO − μ_H₂O + μ_H₂ + μ_CO₂。
    解释性选择：ADR 0010 第 1 条。不使用式 (5.47)。

    Source: Hamel Eq.5.43, S. 92
    """
    chemical_potentials = (
        _finite("mu_co_j_mol", mu_co_j_mol),
        _finite("mu_h2o_j_mol", mu_h2o_j_mol),
        _finite("mu_h2_j_mol", mu_h2_j_mol),
        _finite("mu_co2_j_mol", mu_co2_j_mol),
    )
    temperature = _positive("temperature_k", temperature_k)
    gas_constant = _positive("gas_constant_j_mol_k", gas_constant_j_mol_k)
    summed = -chemical_potentials[0] - chemical_potentials[1] + chemical_potentials[2] + chemical_potentials[3]
    argument = -summed / (gas_constant * temperature)
    return math.exp(argument)


def r8_water_gas_shift_rate(
    y_co, y_h2o, y_h2, y_co2, mu_co_j_mol, mu_h2o_j_mol, mu_h2_j_mol,
    mu_co2_j_mol, temperature_k, pressure_pa, gas_constant_j_mol_k, m_ash_kg,
) -> ReactionRate:
    """R8 用式 (5.45)，不用式 (5.49)。比速率为 mol/(s·kg 灰)。

    y_CO > y_eq 时进度为正，CO 被正向消耗。y_CO < y_eq 时进度为负，保留逆反应。
    p 由 Pa 换算为 atm。a_R8 = 0,02。灰质量不是 Corella 的催化剂质量。
    解释性选择：ADR 0010 第 1、2、11 条；解释性选择（外推）：ADR 0013。

    Source: Hamel Eq.5.45, S. 94
    Source: Hamel Eq.5.46, S. 94
    """
    mole_co = _finite("y_co", y_co)
    mole_water = _finite("y_h2o", y_h2o)
    mole_hydrogen = _finite("y_h2", y_h2)
    mole_co2 = _finite("y_co2", y_co2)
    temperature = _positive("temperature_k", temperature_k)
    pressure = _positive("pressure_pa", pressure_pa)
    ash = _finite("m_ash_kg", m_ash_kg)
    constant = r8_equilibrium_constant(
        mu_co_j_mol,
        mu_h2o_j_mol,
        mu_h2_j_mol,
        mu_co2_j_mol,
        temperature,
        gas_constant_j_mol_k,
    )
    denominator = constant * mole_water
    if denominator == 0.0:
        raise KineticsInputError("式 (5.46) 的分母为 0")
    y_equilibrium = mole_co2 * mole_hydrogen / denominator
    pressure_atm = pressure / _ATM_PA
    first = arrhenius_rate_constant(2.77e8, temperature, activation_temperature_k=13_971.0)
    second = arrhenius_rate_constant(
        math.exp(_R8_SHIFT),
        temperature,
        activation_temperature_k=-5_553.0,
    )
    exponent = 0.5 - pressure_atm / 250.0
    specific = _A_R8 * first * (mole_co - y_equilibrium) * pressure_atm**exponent * second
    return _extent(specific, ash)


def r9_methane_reforming_rate(
    p_ch4_pa: float,
    temperature_k: float,
    gas_constant_j_mol_k: float,
    m_ash_kg: float,
) -> ReactionRate:
    """R9。印刷式是 CH₄ 的消耗速率，进度取其相反数。P_CH₄ 以 Pa 计。

    原文写明只在有固体的悬浮相计入。这里不按相把速率置零，灰质量由调用方传入。
    解释性选择：ADR 0010 第 4、11 条。

    Source: Hamel Eq.5.51, S. 95
    """
    pressure = _finite("p_ch4_pa", p_ch4_pa)
    temperature = _positive("temperature_k", temperature_k)
    gas_constant = _positive("gas_constant_j_mol_k", gas_constant_j_mol_k)
    ash = _finite("m_ash_kg", m_ash_kg)
    coefficient = arrhenius_rate_constant(
        6.113e-2,
        temperature,
        activation_energy_j_mol=137_327.0,
        gas_constant_j_mol_k=gas_constant,
    )
    return _extent(coefficient * pressure, ash)


def r10_tar_oxidation_rate(
    c_tar_mol_m3: float,
    c_o2_mol_m3: float,
    temperature_k: float,
    pressure_pa: float,
    volume_m3: float,
) -> ReactionRate:
    """R10，苯和萘都用芳烃参数。印刷式的前导负号已去掉。十六烷不实现。

    Source: Hamel Eq.5.59, S. 96
    Source: Hamel Eq.5.61, S. 96
    """
    tar = _real_power("c_tar_mol_m3", c_tar_mol_m3, 0.5)
    oxygen = _finite("c_o2_mol_m3", c_o2_mol_m3)
    temperature = _positive("temperature_k", temperature_k)
    pressure = _positive("pressure_pa", pressure_pa)
    volume = _finite("volume_m3", volume_m3)
    coefficient = arrhenius_rate_constant(20_700.0, temperature, activation_temperature_k=9_650.0)
    specific = coefficient * temperature * pressure**0.3 * tar * oxygen
    return _extent(specific, volume)


def r11_serio_rate(
    c_tar_mol_m3: float,
    temperature_k: float,
    gas_constant_j_mol_k: float,
    volume_m3: float,
    *,
    frequency_per_s: float = 5.42e4,
    activation_energy_j_mol: float = 100.5e3,
) -> ReactionRate:
    """气泡相焦油二次反应。默认表 5.4 的第一对，调用方可传入同一行的另一对。一级，乘气体体积。

    解释性选择：ADR 0010 第 9 条。对水没有级数。

    Source: Hamel Tab.5.4, S. 97
    """
    tar = _finite("c_tar_mol_m3", c_tar_mol_m3)
    temperature = _positive("temperature_k", temperature_k)
    gas_constant = _positive("gas_constant_j_mol_k", gas_constant_j_mol_k)
    volume = _finite("volume_m3", volume_m3)
    frequency = _positive("frequency_per_s", frequency_per_s)
    activation = _positive("activation_energy_j_mol", activation_energy_j_mol)
    coefficient = arrhenius_rate_constant(
        frequency,
        temperature,
        activation_energy_j_mol=activation,
        gas_constant_j_mol_k=gas_constant,
    )
    return _extent(coefficient * tar, volume)


def r11_corella_rate(
    c_tar_mol_m3: float,
    temperature_k: float,
    gas_constant_j_mol_k: float,
    m_catalyst_kg: float,
) -> ReactionRate:
    """悬浮相焦油二次反应。频率因子单位 m³/(s·kg 催化剂)，再乘催化剂质量。

    ``m_catalyst_kg`` 是焦炭加灰的质量（ADR 0014 A6），与 R8、R9 的灰质量不是同一个参数。
    740 °C 用高温段。低于 600 °C 或高于 900 °C 时拒绝计算：ADR 0010 第 9 条没有批准外推。
    解释性选择：ADR 0010 第 9 条。

    Source: Hamel Tab.5.4, S. 97
    """
    tar = _finite("c_tar_mol_m3", c_tar_mol_m3)
    temperature = _positive("temperature_k", temperature_k)
    gas_constant = _positive("gas_constant_j_mol_k", gas_constant_j_mol_k)
    catalyst = _finite("m_catalyst_kg", m_catalyst_kg)
    if temperature < _CORELLA_LOW_K or temperature > _CORELLA_HIGH_K:
        raise KineticsInputError("Corella 参数只覆盖 600–900 °C，ADR 0010 第 9 条没有批准外推")
    if temperature >= _CORELLA_SPLIT_K:
        prefactor = 0.7
        energy = 63.1e3
    else:
        prefactor = 7.5e-4
        energy = 5.9e3
    coefficient = arrhenius_rate_constant(
        prefactor,
        temperature,
        activation_energy_j_mol=energy,
        gas_constant_j_mol_k=gas_constant,
    )
    return _extent(coefficient * tar, catalyst)
