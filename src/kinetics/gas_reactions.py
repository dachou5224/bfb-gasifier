"""均相气体反应动力学（R5–R9）。

R5: CO 氧化（气泡/悬浮相两机理）
R6: CH4 氧化
R7: CH4 水蒸气重整
R8: 水煤气变换 (WGSR)
R9: H2S 氧化

Source: Hamel (1999) Eq. 5.36, 5.45-5.47, 5.51; BFB_TechSpec_v11 §5.4
"""

from __future__ import annotations

import numpy as np

from src.core.constants import BAR_PA, P0, Rg
from src.kinetics.arrhenius import k_standard, k_jensen_r7

# -----------------------------------------------------------------------
# 动力学常数（Source: Hamel 1999 Table 5.4）
# -----------------------------------------------------------------------

# R5 气泡相 (Hottel 1965)
R5_bubble_k0: float = 1.91e6    # [ (m³/mol)^0.5 / s ]
R5_bubble_E_Rg: float = 8056.0  # [K]

# R5 悬浮相 (Hayhurst & Tucker 1990)
# r = k * p_CO * p_O2^0.5 * p_H2O^0.5 [atm 基准]
R5_suspension_k_T05: float = 1.8  # 恢复原始物理值

# R6 CH4 氧化 (de Souza-Santos 1989)
R6_k0: float = 3.552e11         # [ K/s ] -> k = k0/T * exp(-E/RT)
R6_E_Rg: float = 15700.0        # [K]

# R7 CH4 重整（保留旧常数名供审计脚本/历史对照使用；主实现改按 Hamel Eq.5.51）
R7_A: float = 1.17e12
R7_E_Rg: float = 24000.0
R7_k0_Pa: float = 6.113e-7      # [mol/(m³·s·Pa)] = 6.113e-2 / bar
R7_E_J: float = 137_327.0       # [J/mol]

# R8 WGSR (Hamel 1999)
R8_k0: float = 2.78e3           # [ (mol/(m³·s)) / bar^n ]
R8_E_Rg: float = 1510.0         # [K]
# Hamel 本地提取记录在煤灰/焦场景多处指向 a_R8 = 0.02。
# 当前 Phase-1 LU/HTW 目标工况为加压褐煤，因此先按该口径对齐。
R8_a_R8: float = 0.02

# R9 H2S oxidation (占位参数；后续可由 Gibbs 替代)
R9_k0: float = 1.0e4
R9_E_Rg: float = 10_000.0

# R12 H2 oxidation (代码扩展，不在严格 R1-R11 基线内)
R12_k0: float = 2.0e7
R12_E_Rg: float = 15_000.0


# -----------------------------------------------------------------------
# R5 CO 氧化
# -----------------------------------------------------------------------

def rate_R5_bubble(
    T: float,
    C_CO: float,
    C_O2: float,
    P: float,
    y: np.ndarray,
    species_index: dict[str, int] | None = None,
) -> float:
    """R5 气泡相 CO 氧化速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.36; Table 5.4

    Hamel 原始论文对 CO 氧化写成纯正向动力学，不附加通用
    ``(1 - Q_p/K_eq)`` 热力学壳。保留 ``P/y/species_index`` 仅为兼容旧接口。
    """
    _ = (P, y, species_index)
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    E_J = R5_bubble_E_Rg * Rg
    k = k_standard(R5_bubble_k0, E_J, T)
    return k * C_CO * C_O2**0.5


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

    Ref: Hamel (1999) Eq. 5.36; Table 5.4

    Hamel 原始论文对 CO 氧化写成纯正向动力学，不附加通用
    ``(1 - Q_p/K_eq)`` 热力学壳。保留 ``P/y/species_index`` 仅为兼容旧接口。
    """
    _ = (P, y, species_index)
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    C_H2O = max(C_H2O, 0.0)
    
    p_CO = (C_CO * Rg * T) / P0
    p_O2 = (C_O2 * Rg * T) / P0
    p_H2O = (C_H2O * Rg * T) / P0
    
    k = R5_suspension_k_T05 * np.sqrt(max(T, 300.0))
    # 引入 1e-4 的 H2O 兜底，防止 R2 耗尽水分导致 CO 氧化停止
    return k * p_CO * (max(p_O2, 0.0)**0.5) * (max(p_H2O, 1e-4)**0.5)


# -----------------------------------------------------------------------
# R6 CH4 氧化
# -----------------------------------------------------------------------

def rate_R6(T: float, C_CH4: float, C_O2: float) -> float:
    """R6 CH4 氧化速率 [mol/(m³·s)]。"""
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
) -> float:
    """R7 CH4 蒸汽重整正向速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.51

    注意：源码沿用本仓库历史编号 ``R7``，但这里实现的是 Hamel 论文中
    甲烷重整那条**单向**动力学。原式依赖 CH4 分压，不显式包含
    ``(1 - Q_p/K_eq)`` 或逆向甲烷化项。
    """
    _ = (C_H2O, C_CO, C_H2)
    T_safe = max(float(T), 300.0)
    p_ch4 = max(float(C_CH4), 0.0) * Rg * T_safe
    k = k_standard(R7_k0_Pa, R7_E_J, T_safe)
    return k * p_ch4


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
) -> float:
    """R8 WGSR 净速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.45-5.47

    按论文原貌使用 ``(y_CO - y_eq,CO)`` 作为平衡驱动力；在数值上与
    ``C_CO*C_H2O - C_CO2*C_H2/K_eq`` 等价，但更便于逐式追溯。
    """
    T_safe = max(float(T), 300.0)
    K_eq = wgsr_equilibrium_constant(T_safe)

    P_bar = P / BAR_PA
    k_std = k_standard(R8_k0, R8_E_Rg * Rg, T_safe)
    exp_corr = np.exp(np.clip(-8.91 + 5.553 / T_safe, -100.0, 100.0))
    P_factor = P_bar ** max(0.5 - P_bar / 250.0, 0.01)
    k_fwd = R8_a_R8 * k_std * P_factor * exp_corr  # [m³/(mol·s)] 校准值

    C_fac = P / (Rg * T_safe)
    y_CO_s = max(float(y_CO), 0.0)
    y_H2O_s = max(float(y_H2O), 0.0)
    y_CO2_s = max(float(y_CO2), 0.0)
    y_H2_s = max(float(y_H2), 0.0)
    y_eq_co = wgsr_equilibrium_y_co(y_H2O_s, y_CO2_s, y_H2_s, K_eq)

    return k_fwd * (C_fac**2) * y_H2O_s * (y_CO_s - y_eq_co)


def wgsr_equilibrium_constant(T: float) -> float:
    """R8 WGSR 平衡常数 K_eq(T)。

    Ref: Hamel (1999) Eq. 5.47
    """
    exponent = np.clip(-3.6893 + 4019.0 / max(T, 300.0), -100.0, 100.0)
    return float(np.exp(exponent))


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

    代码扩展反应，不属于严格的 R1-R11 基线；用于评估快速氢氧化对
    下部 O2 预算与温峰的影响。
    """
    C_H2 = max(C_H2, 0.0)
    C_O2 = max(C_O2, 0.0)
    k = k_standard(R12_k0, R12_E_Rg * Rg, max(T, 300.0))
    return k * C_H2 * C_O2**0.5
