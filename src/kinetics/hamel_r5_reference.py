"""Hamel (1999) Table 5.3 / Eq. 5.36–5.37 CO 氧化参考式。

Source: Hamel (1999) §5.2.1, Table 5.3; Hottel et al. (1965); Hayhurst & Tucker (1990).

说明
----
- 气泡相 (Blase / Hottel):
    k_ch,R5 = k0,R5 * exp(-(E/Rg)/T)                         … Eq. 5.37
    r [mol CO/(m³·s)] = k_ch,R5 * C_CO * C_O2^0.3 * C_H2O^0.5   … Eq. 5.36; Table 5.3 α=1, β=0.3, γ=0.5
    计量：CO + 1/2 O2 → CO2（dC_CO/dt）；耗氧 = 0.5 · r

- 悬浮相 (Suspension / Hayhurst & Tucker):
    k_ch,R5 = 1.8 * T^(-0.5)                                 … Table 5.3  ``1,8·T^(-0,5)``
    r [mol/(m³·s)] = k_ch,R5 * p_CO * p_O2^0.5 * p_H2O^0.5   … α=1, β=0.5；H2O 来自 Hayhurst

分压 p_i [atm] 由理想气体 C_i [mol/m³] 换算：p_i = C_i * Rg * T / P0。
"""

from __future__ import annotations

import math

from src.core.constants import P0, Rg

# Table 5.3 — Blase / Hottel et al. (1965)
HAMEL_R5_BUBBLE_K0: float = 1.91e6
HAMEL_R5_BUBBLE_E_OVER_RG_K: float = 8056.0
HAMEL_R5_BUBBLE_ALPHA: float = 1.0
HAMEL_R5_BUBBLE_BETA: float = 0.3
HAMEL_R5_BUBBLE_GAMMA: float = 0.5

# Table 5.3 — Suspension / Hayhurst & Tucker (1990)
HAMEL_R5_SUSPENSION_PREFACTOR: float = 1.8
HAMEL_R5_SUSPENSION_ALPHA: float = 1.0
HAMEL_R5_SUSPENSION_BETA: float = 0.5


def hamel_r5_bubble_k(T_K: float) -> float:
    """Eq. 5.37 气泡相速率常数 k_ch,R5（浓度形式，与 Eq. 5.36 配套）。"""
    t = float(max(T_K, 300.0))
    x = -HAMEL_R5_BUBBLE_E_OVER_RG_K / t
    x = float(min(max(x, -100.0), 100.0))
    return HAMEL_R5_BUBBLE_K0 * math.exp(x)


def hamel_r5_bubble_rate_mol_m3_s(
    T_K: float,
    C_CO_mol_m3: float,
    C_O2_mol_m3: float,
    C_H2O_mol_m3: float,
) -> float:
    """Hottel 气泡相 CO 氧化体积速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.36
    """
    k = hamel_r5_bubble_k(T_K)
    return (
        k
        * max(C_CO_mol_m3, 0.0)
        * (max(C_O2_mol_m3, 0.0) ** HAMEL_R5_BUBBLE_BETA)
        * (max(C_H2O_mol_m3, 0.0) ** HAMEL_R5_BUBBLE_GAMMA)
    )


def _partial_pressure_atm(C_mol_m3: float, T_K: float) -> float:
    return max(float(C_mol_m3), 0.0) * Rg * float(T_K) / P0


def hamel_r5_suspension_k(T_K: float, *, temperature_exponent: float = -0.5) -> float:
    """Table 5.3 悬浮相前置因子 1.8·T^exp。

    默认 exp=-0.5（OCR ``1,8·T^(-0,5)``）；``temperature_exponent=+0.5`` 用于对照当前实现。
    """
    t = float(max(T_K, 300.0))
    return HAMEL_R5_SUSPENSION_PREFACTOR * (t ** float(temperature_exponent))


def hamel_r5_suspension_rate_mol_m3_s(
    T_K: float,
    C_CO_mol_m3: float,
    C_O2_mol_m3: float,
    C_H2O_mol_m3: float,
    *,
    temperature_exponent: float = -0.5,
) -> float:
    """Hayhurst & Tucker 悬浮相 CO 氧化体积速率 [mol/(m³·s)]。"""
    k = hamel_r5_suspension_k(T_K, temperature_exponent=temperature_exponent)
    p_co = _partial_pressure_atm(C_CO_mol_m3, T_K)
    p_o2 = _partial_pressure_atm(C_O2_mol_m3, T_K)
    p_h2o = _partial_pressure_atm(C_H2O_mol_m3, T_K)
    return k * p_co * (p_o2**0.5) * (p_h2o**0.5)
