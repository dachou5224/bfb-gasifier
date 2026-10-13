"""均相气体反应动力学（R5–R9）。

R5: CO 氧化（气泡/悬浮相两机理）
R6: CH4 氧化
R7: CH4 水蒸气重整
R8: 水煤气变换 (WGSR)
R9: H2S 氧化

Source: ``docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md``;
Hamel (1999) Eq. 5.36, 5.45-5.47, 5.51
"""

from __future__ import annotations

import math
import numpy as np

from src.core.constants import BAR_PA, P0, Rg
from src.kinetics.arrhenius import k_jensen_h2, k_jensen_r7, k_standard
from src.kinetics.hamel_r6_h2_reference import (
    HAMEL_R6_H2_A,
    HAMEL_R6_H2_E_OVER_RG_K,
    HAMEL_R6_H2_H2_ORDER,
    HAMEL_R6_H2_O2_ORDER,
)
from src.kinetics.hamel_r7_ch4_reference import HAMEL_R7_CH4_A, HAMEL_R7_CH4_E_OVER_RG_K

# -----------------------------------------------------------------------
# 动力学常数（Source: Hamel 1999 Table 5.4）
# -----------------------------------------------------------------------

# R5 气泡相 (Hottel 1965; Hamel Table 5.3 / Eq. 5.36–5.37)
R5_bubble_k0: float = 1.91e6    # 与 Eq. 5.36 浓度幂次配套
R5_bubble_E_Rg: float = 8056.0  # [K]
R5_bubble_beta: float = 0.3     # [-] Table 5.3 O2 级数
R5_bubble_gamma: float = 0.5    # [-] Table 5.3 H2O 级数

# R5 悬浮相 (Hayhurst & Tucker 1990; Hamel Table 5.3)
# k = 1.8 * T^(-0.5); r = k * p_CO * p_O2^0.5 * p_H2O^0.5 [atm 基准]
R5_suspension_k_T05: float = 1.8

# R6 CH4 氧化 (de Souza-Santos 1989; Hamel Eq. 5.40–5.41, thesis R7)
R6_k0: float = HAMEL_R7_CH4_A
R6_E_Rg: float = HAMEL_R7_CH4_E_OVER_RG_K

# R7 CH4 重整（仓库 R7 = thesis R9；Hamel Eq. 5.51）
# 原式 mol/(s·kg_Asche)，p_CH4 按 bar→Pa：k0_Pa = 6.113e-2 / BAR_PA。
# 体积源 r_V = ρ_ash * R_m；无灰时为零。仅密相/载固悬浮相。
R7_A: float = 1.17e12
R7_E_Rg: float = 24000.0
R7_k0_mass_bar: float = 6.113e-2  # [mol/(s·kg_ash·bar)] Hamel Eq. 5.51
R7_k0_Pa: float = R7_k0_mass_bar / BAR_PA  # [mol/(s·kg_ash·Pa)]
R7_E_J: float = 137_327.0       # [J/mol] 德语 137.327 为千分符

# R8 WGSR (Hamel 1999 Eq. 5.45–5.49; Wen & Chaung 1979)
R8_k0: float = 2.77e8           # [1/s] per kg_ash，与 (y_CO - y_eq,CO) 相乘
R8_E_Rg: float = 13971.0        # [K] exp(-13971/T_g)
R8_p_corr_T: float = 5553.0     # [K] exp(-8.91 + 5553/T_g)
R8_a_R8: float = 0.02           # [-] Eq. 5.48 煤灰催化因子
R8_rho_ash_ref: float = 75.0    # [kg/m³] 无 cell 灰含率时的参考灰密度
# Eq. (5.49) SI 形式（page_57.png）
R8_k0_si: float = 8.82e-7       # [mol/(s·kg_ash·Pa²)] 与分压乘积项配套
R8_E_si_over_Rg: float = 90.853  # [K] exp(-90.853/(Rg·T))

# R9 H2S oxidation (占位参数；后续可由 Gibbs 替代)
R9_k0: float = 1.0e4
R9_E_Rg: float = 10_000.0

# R12 H2 氧化 (Jensen 1995; Hamel Eq. 5.38–5.39, thesis R6)
R12_A: float = HAMEL_R6_H2_A
R12_E_Rg: float = HAMEL_R6_H2_E_OVER_RG_K
R12_h2_order: float = HAMEL_R6_H2_H2_ORDER
R12_o2_order: float = HAMEL_R6_H2_O2_ORDER
R12_k0: float = R12_A  # 兼容旧名


def _clamp_exp_arg(x: float) -> float:
    val = float(x)
    if val < -100.0:
        return -100.0
    if val > 100.0:
        return 100.0
    return val


# -----------------------------------------------------------------------
# R5 CO 氧化
# -----------------------------------------------------------------------

def rate_R5_bubble(
    T: float,
    C_CO: float,
    C_O2: float,
    C_H2O: float,
    P: float,
    y: np.ndarray,
    species_index: dict[str, int] | None = None,
) -> float:
    """R5 气泡相 CO 氧化速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.36–5.37; Table 5.3

    r = k_ch,R5 * C_CO * C_O2^0.3 * C_H2O^0.5

    Hamel 原始论文对 CO 氧化写成纯正向动力学，不附加通用
    ``(1 - Q_p/K_eq)`` 热力学壳。保留 ``P/y/species_index`` 仅为兼容旧接口。
    """
    _ = (P, y, species_index)
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    C_H2O = max(C_H2O, 0.0)
    E_J = R5_bubble_E_Rg * Rg
    k = k_standard(R5_bubble_k0, E_J, T)
    return k * C_CO * (C_O2**R5_bubble_beta) * (C_H2O**R5_bubble_gamma)


def rate_R5_suspension(
    T: float,
    C_CO: float,
    C_O2: float,
    C_H2O: float,
    P: float,
    y: np.ndarray,
    species_index: dict[str, int] | None = None,
) -> float:
    """R5 悬浮相 CO 氧化速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Table 5.3; Hayhurst & Tucker (1990)

    Hamel 原始论文对 CO 氧化写成纯正向动力学，不附加通用
    ``(1 - Q_p/K_eq)`` 热力学壳。保留 ``P/y/species_index`` 仅为兼容旧接口。
    论文原式显式依赖 ``p_H2O^0.5``；此处不再保留项目历史上的蒸汽分压 floor。
    """
    _ = (P, y, species_index)
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    C_H2O = max(C_H2O, 0.0)

    p_CO = (C_CO * Rg * T) / P0
    p_O2 = (C_O2 * Rg * T) / P0
    p_H2O = (C_H2O * Rg * T) / P0

    k = R5_suspension_k_T05 * (max(T, 300.0) ** -0.5)
    return k * p_CO * (max(p_O2, 0.0) ** 0.5) * (max(p_H2O, 0.0) ** 0.5)


# -----------------------------------------------------------------------
# R6 CH4 氧化
# -----------------------------------------------------------------------

def rate_R6(T: float, C_CH4: float, C_O2: float) -> float:
    """R6 CH4 氧化速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.40–5.41; de Souza-Santos (1989)

    仓库编号 ``R6`` = thesis **R7**；r = k_ch,R7 * C_CH4 * C_O2。
    """
    C_CH4 = max(C_CH4, 0.0)
    C_O2 = max(C_O2, 0.0)
    k = k_jensen_r7(R6_k0, R6_E_Rg, max(T, 300.0))
    return k * C_CH4 * C_O2


# -----------------------------------------------------------------------
# R7 CH4 水蒸气重整
# -----------------------------------------------------------------------

def rate_R7(
    T: float,
    C_CH4: float,
    C_H2O: float,
    C_CO: float = 0.0,
    C_H2: float = 0.0,
    rho_ash_kg_m3: float | None = None,
) -> float:
    """R7 CH4 蒸汽重整体积速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.50–5.51；PDF 57 / 印刷 94–95

    仓库 ``R7`` = thesis R9。原式
    ``R_CH4 = -6.113e-2 exp(-137327/(Rg T)) p_CH4``，单位 mol/(s·kg_Asche)。
    体积源 ``r_V = ρ_ash R_m``，``R_m`` 取消耗为正。``ρ_ash<=0`` 或未传灰密度时为 0。
    原文：只在载固悬浮相计算；气泡与无灰自由空间为 0。
    原式依赖 CH4 分压，不含 ``(1 - Q_p/K_eq)``。
    ``p_CH4 = C_CH4 Rg T`` [Pa]；``k0`` 由 bar 基准除以 ``BAR_PA``。
    """
    _ = (C_H2O, C_CO, C_H2)
    rho_ash = 0.0 if rho_ash_kg_m3 is None else max(float(rho_ash_kg_m3), 0.0)
    if rho_ash <= 0.0:
        return 0.0
    T_safe = max(float(T), 300.0)
    p_ch4 = max(float(C_CH4), 0.0) * Rg * T_safe  # [Pa]
    k_mass = k_standard(R7_k0_Pa, R7_E_J, T_safe)  # [mol/(s·kg_ash·Pa)]
    return rho_ash * k_mass * p_ch4


# -----------------------------------------------------------------------
# R8 WGSR
# -----------------------------------------------------------------------

def rate_R8(
    T: float,
    P: float,
    y_CO: float,
    y_H2O: float,
    y_CO2: float,
    y_H2: float,
    rho_ash_kg_m3: float | None = None,
) -> float:
    """R8 WGSR 净速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.45–5.47

    论文原式单位 mol CO/(s·kg_灰)，密相体积速率 ``ρ_ash * R_mass``。
    逆 WGSR 由 ``(y_CO - y_eq,CO)`` 自动体现，无单独逆反应参数。
    压力项内部使用 atm（``P_atm = P / P0``）。
    """
    T_safe = max(float(T), 300.0)
    rho_ash = max(float(R8_rho_ash_ref if rho_ash_kg_m3 is None else rho_ash_kg_m3), 0.0)
    if rho_ash <= 0.0:
        return 0.0

    K_eq = wgsr_equilibrium_constant(T_safe)
    P_atm = max(float(P), 0.0) / P0
    k_arr = k_standard(R8_k0, R8_E_Rg * Rg, T_safe)
    exp_corr = math.exp(_clamp_exp_arg(-8.91 + R8_p_corr_T / T_safe))
    p_exp = 0.5 - P_atm / 250.0  # Hamel Eq. (5.45) 原式，不做论文外裁剪
    k_mass = R8_a_R8 * k_arr * (P_atm ** p_exp) * exp_corr  # [1/s] per kg_ash

    y_CO_s = max(float(y_CO), 0.0)
    y_H2O_s = max(float(y_H2O), 0.0)
    y_CO2_s = max(float(y_CO2), 0.0)
    y_H2_s = max(float(y_H2), 0.0)
    y_eq_co = wgsr_equilibrium_y_co(y_H2O_s, y_CO2_s, y_H2_s, K_eq)

    return rho_ash * k_mass * (y_CO_s - y_eq_co)


def rate_R8_si(
    T: float,
    P: float,
    p_CO: float,
    p_H2O: float,
    p_CO2: float,
    p_H2: float,
    rho_ash_kg_m3: float | None = None,
) -> float:
    """R8 WGSR 速率 SI 形式 [mol/(m³·s)]（Hamel Eq. 5.49）。

    R_CO = k0 · exp(-E/(Rg·T)) · [p_CO·p_H2O - p_H2·p_CO2 / k_eq]

    分压 ``p_*`` 与 ``P`` 均为 Pa；``k_eq`` 来自 Eq. (5.47) 修正转录（见 ``wgsr_equilibrium_constant``）。

    Source: ``docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md``; Hamel Eq. 5.49
    """
    T_safe = max(float(T), 300.0)
    rho_ash = max(float(R8_rho_ash_ref if rho_ash_kg_m3 is None else rho_ash_kg_m3), 0.0)
    if rho_ash <= 0.0:
        return 0.0

    K_eq = wgsr_equilibrium_constant(T_safe)
    k_arr = k_standard(R8_k0_si, R8_E_si_over_Rg * Rg, T_safe)
    drive_pa2 = max(float(p_CO), 0.0) * max(float(p_H2O), 0.0) - (
        max(float(p_H2), 0.0) * max(float(p_CO2), 0.0) / max(K_eq, 1e-30)
    )
    rate_mass = k_arr * drive_pa2  # [mol/(s·kg_ash)]
    return rho_ash * rate_mass


def wgsr_equilibrium_constant(T: float) -> float:
    """R8 WGSR 平衡常数 K_eq(T)（Hamel Eq. 5.47）。

    PDF page_57 字面 OCR 为 ``exp(-36.893 + 4.019/T_g)``，物理上不合理；
    德语数字 ``4.019`` 为千分符 → ``4019/T``，``-36,893`` 疑为 ``-3,6893``。
    实现采用 ``exp(-3.6893 + 4019/T)``，与 ``src/thermodynamics/equilibrium.py`` 一致。

    Ref: Hamel (1999) Eq. 5.47; ``data/hamel_pages_png/page_57.png``
    """
    exponent = _clamp_exp_arg(-3.6893 + 4019.0 / max(T, 300.0))
    return float(math.exp(exponent))


def wgsr_equilibrium_y_co(y_H2O: float, y_CO2: float, y_H2: float, K_eq: float) -> float:
    """WGSR 平衡 CO 摩尔分数 ``y_eq,CO``。

    Ref: Hamel (1999) Eq. 5.46
    """
    denom = max(float(K_eq) * max(float(y_H2O), 1e-30), 1e-30)
    return max(float(y_CO2), 0.0) * max(float(y_H2), 0.0) / denom


# -----------------------------------------------------------------------
# R9 H2S 氧化
# -----------------------------------------------------------------------

def rate_R9(T: float, C_H2S: float, C_O2: float) -> float:
    """R9 H2S 氧化速率 [mol/(m³·s)]。

    占位一级动力学；当 `use_gibbs_minor=True` 时由 Gibbs 微量组分平衡替代。
    """
    C_H2S = max(C_H2S, 0.0)
    C_O2 = max(C_O2, 0.0)
    k = k_standard(R9_k0, R9_E_Rg * Rg, max(T, 300.0))
    return k * C_H2S * C_O2


# -----------------------------------------------------------------------
# R12 H2 氧化
# -----------------------------------------------------------------------

def rate_R12(
    T: float,
    C_H2: float,
    C_O2: float,
    P: float,
    y: np.ndarray,
    species_index: dict[str, int] | None = None,
) -> float:
    """R12 H2 + 0.5 O2 -> H2O 速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.38–5.39; Jensen et al. (1995)

    仓库编号 ``R12`` = thesis **R6**；
    ``r = k_ch,M * C_H2^{3/2} * C_O2``（扫描原式，2026-09-14 二次诊断）。
    """
    _ = (P, y, species_index)
    C_H2 = max(C_H2, 0.0)
    C_O2 = max(C_O2, 0.0)
    k = k_jensen_h2(R12_A, R12_E_Rg, max(T, 300.0))
    return k * (C_H2**R12_h2_order) * (C_O2**R12_o2_order)
