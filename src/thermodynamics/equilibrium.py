"""平衡常数 K_eq 与反应商 Q_p 计算。

用于可逆反应驱动力 (1 − Q_p/K_eq)。
Source: BFB_TechSpec_v11 §5; Hamel (1999) §5.2
"""

from __future__ import annotations

from typing import Dict

import numpy as np

from src.core.constants import Rg
from src.core.species import gibbs_molar

# exp 安全区间，避免 inf/下溢为 0 导致 Q_p/K_eq 除零或 RuntimeWarning
_EXP_CLIP: float = 700.0
_K_EQ_FLOOR: float = 1e-300


def _exp_clipped(x: float) -> float:
    """``exp(x)``，将 ``x`` 裁剪到 ±700 附近，避免溢出/下溢为 0。"""
    return float(np.exp(np.clip(x, -_EXP_CLIP, _EXP_CLIP)))


# 反应 ID（与 gas_reactions 一致）
REACTION_R5 = "R5"   # 2 CO + O2 -> 2 CO2
REACTION_R7 = "R7"   # CH4 + H2O -> CO + 3 H2
REACTION_R8 = "R8"   # CO + H2O <-> CO2 + H2 (WGSR)


def _delta_g0(reaction_id: str, T: float) -> float:
    """反应标准 Gibbs 自由能变化 ΔG°(T) [J/mol]。

    以正向反应（消耗左侧、生成右侧）定义。
    """
    if reaction_id == REACTION_R5:
        # 2 CO + O2 -> 2 CO2
        dG = 2 * gibbs_molar("CO2", T) - 2 * gibbs_molar("CO", T) - gibbs_molar("O2", T)
    elif reaction_id == REACTION_R7:
        # CH4 + H2O -> CO + 3 H2
        dG = (
            gibbs_molar("CO", T) + 3 * gibbs_molar("H2", T)
            - gibbs_molar("CH4", T) - gibbs_molar("H2O", T)
        )
    elif reaction_id == REACTION_R8:
        # CO + H2O -> CO2 + H2
        dG = (
            gibbs_molar("CO2", T) + gibbs_molar("H2", T)
            - gibbs_molar("CO", T) - gibbs_molar("H2O", T)
        )
    else:
        raise ValueError(f"未知反应 ID: {reaction_id}")
    return dG


# R8 WGSR: Benson (1981) 拟合式，Hamel 给式量级与 sanity check 不符
R8_K_eq_A: float = -4.33
R8_K_eq_B: float = 4577.8  # K_eq = exp(4577.8/T - 4.33)


def get_K_eq(reaction_id: str, T: float, P: float | None = None) -> float:
    """平衡常数 K_eq（压力无关形式，K_p = K_eq * P^Δn）。

    R5/R7: ΔG° = −RT ln K_eq
    R8: Benson 拟合式 K_eq = exp(4577.8/T - 4.33)

    Parameters
    ----------
    reaction_id : str
        R5, R7, R8 之一
    T : float
        温度 [K]
    P : float | None
        压力 [Pa]，当前未用于 K_eq 计算

    Returns
    -------
    float
        K_eq（无量纲）
    """
    _ = P
    T_safe = max(T, 300.0)
    if reaction_id == REACTION_R8:
        return _exp_clipped(R8_K_eq_B / T_safe + R8_K_eq_A)
    dG = _delta_g0(reaction_id, T_safe)
    return _exp_clipped(-dG / (Rg * T_safe))


def calc_reaction_quotient(
    reaction_id: str,
    y: Dict[str, float] | np.ndarray,
    P: float,
    species_index: Dict[str, int] | None = None,
) -> float:
    """反应商 Q_p（基于分压）。

    Q_p = Π (P_j)^ν_j，产物为正、反应物为负。

    Parameters
    ----------
    reaction_id : str
        R5, R7, R8 之一
    y : dict or array
        摩尔分数。若为 dict，键为组分名；若为 array，需提供 species_index
    P : float
        总压 [Pa]
    species_index : dict | None
        {组分名: 索引}，仅当 y 为 array 时使用

    Returns
    -------
    float
        Q_p（无量纲，分压比）
    """
    if isinstance(y, np.ndarray) and species_index is not None:
        y_dict = {sp: float(y[i]) for sp, i in species_index.items()}
    elif isinstance(y, dict):
        y_dict = y
    else:
        raise ValueError("y 需为 dict 或 (array, species_index) 组合")

    def p(sp: str) -> float:
        return max(y_dict.get(sp, 0.0), 1e-30) * P

    if reaction_id == REACTION_R5:
        # 2 CO + O2 -> 2 CO2
        # Q_p = P_CO2^2 / (P_CO^2 * P_O2)
        return p("CO2") ** 2 / (p("CO") ** 2 * p("O2"))
    if reaction_id == REACTION_R7:
        # CH4 + H2O -> CO + 3 H2
        # Q_p = (P_CO * P_H2^3) / (P_CH4 * P_H2O)
        return (p("CO") * p("H2") ** 3) / (p("CH4") * p("H2O"))
    if reaction_id == REACTION_R8:
        # CO + H2O -> CO2 + H2
        # Q_p = (P_CO2 * P_H2) / (P_CO * P_H2O)
        return (p("CO2") * p("H2")) / (p("CO") * p("H2O"))

    raise ValueError(f"未知反应 ID: {reaction_id}")


def calc_gibbs_driving_force(
    reaction_id: str,
    T: float,
    P: float,
    y: Dict[str, float] | np.ndarray,
    species_index: Dict[str, int] | None = None,
    clamp_irreversible: bool = False,
) -> float:
    """计算驱动力 (1 − Q_p/K_eq)。

    Parameters
    ----------
    reaction_id : str
        R5, R7, R8 之一
    T, P : float
        温度 [K]、压力 [Pa]
    y : dict or array
        摩尔分数
    species_index : dict | None
        仅当 y 为 array 时使用
    clamp_irreversible : bool
        若 True，返回 max(0, 1 − Q_p/K_eq)（用于不可逆反应如 R5）

    Returns
    -------
    float
        驱动力，可负（表示逆反应）
    """
    K_eq = float(get_K_eq(reaction_id, T, P))
    if not np.isfinite(K_eq) or K_eq <= 0.0:
        K_eq = _K_EQ_FLOOR
    else:
        K_eq = max(K_eq, _K_EQ_FLOOR)
    Q_p = calc_reaction_quotient(reaction_id, y, P, species_index)
    driving = 1.0 - Q_p / K_eq
    if clamp_irreversible:
        return max(driving, 0.0)
    return driving
