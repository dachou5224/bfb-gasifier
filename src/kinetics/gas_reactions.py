"""R5–R9 均相气相反应。

R5: CO 氧化（气泡相 Hottel；悬浮相 Hayhurst & Tucker k=1.8·T^0.5）
R6: CH4 氧化（de Souza-Santos 1989）
R7: CH4 水蒸气重整（k_jensen_r7）
R8: WGSR（Hamel §5.2.4 催化速率，a_R8=0.02）
R9: H2S 氧化（简化处理）

Source: Hamel (1999) §5.2, Table 5.4; Hottel et al. (1965); Hayhurst & Tucker (1990);
        de Souza-Santos (1989); Jensen et al.
"""

from __future__ import annotations

import numpy as np

from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX
from src.kinetics.arrhenius import k_jensen_r7, k_standard
from src.thermodynamics.equilibrium import calc_gibbs_driving_force, get_K_eq


def wgsr_equilibrium_constant(T: float) -> float:
    """WGSR 平衡常数 K_eq（向后兼容，内部调用 equilibrium.get_K_eq）。"""
    return get_K_eq("R8", T)

# -----------------------------------------------------------------------
# 参数来源：Hamel (1999) §5.2, Table 5.4; docs/missing_parameters_summary.md §B
# -----------------------------------------------------------------------

# R5 CO 氧化
# 气泡相: Hottel et al. (1965), k = k0·exp(-E/Rg/T)
R5_bubble_k0: float = 1.91e6    # [1/s]
R5_bubble_E_Rg: float = 8056.0  # [K] (= E/Rg)
# 悬浮相: Hayhurst & Tucker (1990), k = 1.8·T^0.5（无 Arrhenius）
R5_suspension_k_T05: float = 1.8

# R6 CH4 氧化: de Souza-Santos (1989), k = 3.552e11/T · exp(-15700/T)
R6_k0: float = 3.552e11    # [1/s]
R6_E_Rg: float = 15_700.0  # [K]

# R7 CH4 + H2O 重整: Jensen & Sørensen（保留原参数）
R7_A: float = 3.0e5
R7_E_T: float = 15_000.0   # [K]

# R8 WGSR: Hamel §5.2.4 催化速率
# R_CO = a_R8 · 2.77e8 · (y_CO - y_eq) · exp(-13971/T) · P^(0.5-P/250) · exp(-8.91+5.553/T)
R8_a_R8: float = 0.02      # 煤灰催化因子
R8_k0: float = 2.77e8     # 前置系数
R8_E_Rg: float = 13_971.0 # [K]
# K_eq 由 thermodynamics.equilibrium.get_K_eq("R8", T) 提供（Benson 形式）

# R9 H2S 氧化（简化）
R9_k0: float = 1.0e8
R9_E: float = 80_000.0    # [J/mol]


# -----------------------------------------------------------------------
# R5 CO 氧化（TechSpec Eq. 5-1：驱动力 max(0, 1−Q_p/K_eq)）
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

    r = k · C_CO · C_O2^0.5 · max(0, 1 − Q_p/K_eq)
    k = 1.91e6 · exp(-8056/T)

    Source: Hottel et al. (1965); BFB_TechSpec_v11 Eq. 5-1
    """
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

    r = k · C_CO · C_O2^0.5 · C_H2O^0.5 · max(0, 1 − Q_p/K_eq)
    k = 1.8 · T^0.5（悬浮相表面催化，无 Arrhenius）

    Source: Hayhurst & Tucker (1990); BFB_TechSpec_v11 Eq. 5-1
    """
    C_CO = max(C_CO, 0.0)
    C_O2 = max(C_O2, 0.0)
    C_H2O = max(C_H2O, 0.0)
    k = R5_suspension_k_T05 * np.sqrt(max(T, 300.0))
    r_kinetic = k * C_CO * C_O2**0.5 * C_H2O**0.5
    driving = calc_gibbs_driving_force(
        "R5", T, P, y, species_index or GAS_SPECIES_INDEX, clamp_irreversible=True
    )
    return r_kinetic * driving


# -----------------------------------------------------------------------
# R6 CH4 氧化
# -----------------------------------------------------------------------

def rate_R6(T: float, C_CH4: float, C_O2: float) -> float:
    """R6 CH4 氧化速率 [mol/(m³·s)]。

    CH4 + 1.5 O2 -> CO + 2 H2O
    k = 3.552e11/T · exp(-15700/T)

    Source: de Souza-Santos (1989); Hamel (1999) §5.2, Table 5.4
    """
    C_CH4 = max(C_CH4, 0.0)
    C_O2 = max(C_O2, 0.0)
    k = (R6_k0 / max(T, 300.0)) * np.exp(np.clip(-R6_E_Rg / T, -100.0, 100.0))
    return k * C_CH4 * C_O2


# -----------------------------------------------------------------------
# R7 CH4 水蒸气重整（TechSpec Eq. 5-1：驱动力 (1−Q_p/K_eq)，可逆）
# -----------------------------------------------------------------------

def rate_R7(
    T: float,
    C_CH4: float,
    C_H2O: float,
    P: float,
    y: np.ndarray,
    species_index: dict[str, int] | None = None,
) -> float:
    """R7 CH4 + H2O → CO + 3H2 速率 [mol/(m³·s)]。

    r = k7 * C_CH4 * C_H2O * (1 − Q_p/K_eq)，可负表示逆反应

    Source: Jensen & Sørensen; BFB_TechSpec_v11 Eq. 5-1
    """
    C_CH4 = max(C_CH4, 0.0)
    C_H2O = max(C_H2O, 0.0)
    k = k_jensen_r7(R7_A, R7_E_T, T)
    r_kinetic = k * C_CH4 * C_H2O
    driving = calc_gibbs_driving_force(
        "R7", T, P, y, species_index or GAS_SPECIES_INDEX, clamp_irreversible=False
    )
    return r_kinetic * driving


# -----------------------------------------------------------------------
# R8 WGSR（TechSpec Eq. 5-1：驱动力 (1−Q_p/K_eq)，可逆）
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

    R = k_eff · C_total · (1 − Q_p/K_eq)，可负表示逆反应
    k_eff 含 Hamel 催化因子 a_R8、压力修正等

    Source: Hamel (1999) §5.2.4; BFB_TechSpec_v11 Eq. 5-1
    """
    y_dict = {"CO": y_CO, "H2O": y_H2O, "CO2": y_CO2, "H2": y_H2}
    driving = calc_gibbs_driving_force(
        "R8", T, P, y_dict, clamp_irreversible=False
    )

    P_bar = P / 1e5  # Pa -> bar
    exp_arrhenius = np.exp(np.clip(-R8_E_Rg / T, -100.0, 100.0))
    exp_corr = np.exp(-8.91 + 5.553 / max(T, 300.0))
    P_factor = P_bar ** max(0.5 - P_bar / 250.0, 0.01)

    k_eff = R8_a_R8 * R8_k0 * exp_arrhenius * P_factor * exp_corr
    C_total = P / (Rg * T) if T > 0 else 0.0
    return k_eff * C_total * driving


# -----------------------------------------------------------------------
# R9 H2S 氧化（简化一级）
# -----------------------------------------------------------------------

def rate_R9(T: float, C_H2S: float, C_O2: float) -> float:
    """R9 H2S 氧化速率 [mol/(m³·s)]。

    H2S + 1.5 O2 -> SO2 + H2O

    TODO: 简化一级动力学，待文献确认详细参数

    Source: 工程估值
    """
    C_H2S = max(C_H2S, 0.0)
    C_O2 = max(C_O2, 0.0)
    k = k_standard(R9_k0, R9_E, T)
    return k * C_H2S * C_O2**0.5
