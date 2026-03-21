"""均相气体反应动力学（R5–R9）。

R5: CO 氧化（气泡/悬浮相两机理）
R6: CH4 氧化
R7: CH4 水蒸气重整
R8: 水煤气变换 (WGSR)
R9: H2S 氧化

Source: Hamel (1999) §5.2; BFB_TechSpec_v11 §5.4
"""

from __future__ import annotations

import numpy as np

from src.core.constants import Rg, P0
from src.core.species import GAS_SPECIES_INDEX
from src.kinetics.arrhenius import k_standard, k_jensen_r7
from src.thermodynamics.equilibrium import calc_gibbs_driving_force

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

# R7 CH4 重整 (Jensen & Sørensen)
R7_A: float = 1.17e12           # [ K/(Pa·s) ] -> k = A/T * exp(-E/RT)
R7_E_Rg: float = 24000.0        # [K]

# R8 WGSR (Hamel 1999)
R8_k0: float = 2.78e3           # [ (mol/(m³·s)) / bar^n ]
R8_E_Rg: float = 1510.0         # [K]
R8_a_R8: float = 1.0            # 催化因子（褐煤默认 1.0）


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
    """R5 气泡相 CO 氧化速率 [mol/(m³·s)]。"""
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    E_J = R5_bubble_E_Rg * Rg
    k = k_standard(R5_bubble_k0, E_J, T)
    r_kinetic = k * C_CO * C_O2**0.5
    driving = calc_gibbs_driving_force(
        "R5", T, P, y, species_index or GAS_SPECIES_INDEX, clamp_irreversible=True
    )
    return r_kinetic * driving


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

    Source: Hayhurst & Tucker (1990); Hamel (1999) Table 5.4
    """
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    C_H2O = max(C_H2O, 0.0)
    
    p_CO = (C_CO * Rg * T) / 101325.0
    p_O2 = (C_O2 * Rg * T) / 101325.0
    p_H2O = (C_H2O * Rg * T) / 101325.0
    
    k = R5_suspension_k_T05 * np.sqrt(max(T, 300.0))
    # 引入 1e-4 的 H2O 兜底，防止 R2 耗尽水分导致 CO 氧化停止
    r_kinetic = k * p_CO * (max(p_O2, 0.0)**0.5) * (max(p_H2O, 1e-4)**0.5)
    
    driving = calc_gibbs_driving_force(
        "R5", T, P, y, species_index or GAS_SPECIES_INDEX, clamp_irreversible=True
    )
    return r_kinetic * driving


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
    P: float,
    y: np.ndarray,
    species_index: dict[str, int] | None = None,
) -> float:
    """R7 CH4 + H2O → CO + 3H2 速率 [mol/(m³·s)]。"""
    C_CH4 = max(C_CH4, 0.0)
    C_H2O = max(C_H2O, 0.0)
    k = k_jensen_r7(R7_A, R7_E_Rg, max(T, 300.0))
    driving = calc_gibbs_driving_force(
        "R7", T, P, y, species_index or GAS_SPECIES_INDEX, clamp_irreversible=True
    )
    return k * C_CH4 * C_H2O * driving


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
    """R8 WGSR 净速率 [mol/(m³·s)]。"""
    y_dict = {"CO": y_CO, "H2O": y_H2O, "CO2": y_CO2, "H2": y_H2}
    driving = calc_gibbs_driving_force(
        "R8", T, P, y_dict, clamp_irreversible=False
    )

    P_bar = P / 1e5
    k_std = k_standard(R8_k0, R8_E_Rg * Rg, T)
    exp_corr = np.exp(np.clip(-8.91 + 5.553 / max(T, 300.0), -100.0, 100.0))
    P_factor = P_bar ** max(0.5 - P_bar / 250.0, 0.01)

    k_eff = R8_a_R8 * k_std * P_factor * exp_corr
    C_total = P / (Rg * T) if T > 0 else 0.0
    return k_eff * C_total * driving


def wgsr_equilibrium_constant(T: float) -> float:
    """R8 WGSR 平衡常数 K_eq(T)。

    Source: Hamel (1999) Eq.5.46
    """
    return float(np.exp(-4.33 + 4577.8 / max(T, 300.0)))


# -----------------------------------------------------------------------
# R9 H2S 氧化
# -----------------------------------------------------------------------

def rate_R9(T: float, C_H2S: float, C_O2: float) -> float:
    """R9 H2S 氧化速率 [mol/(m³·s)]。"""
    C_H2S = max(C_H2S, 0.0)
    C_O2 = max(C_O2, 0.0)
    # 简化一级动力学
    k = 1.0e4 * np.exp(-10000.0 / max(T, 300.0))
    return k * C_H2S * C_O2
