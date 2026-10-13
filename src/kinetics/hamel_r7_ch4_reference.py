"""Hamel (1999) Eq. 5.40–5.41 CH₄ 氧化参考式（thesis R7）。

Source: Hamel (1999) §5.2.3; de Souza-Santos (1989)

说明
----
仓库实现编号 ``R6`` 对应 thesis **R7**（CH₄ + 1.5O₂ → CO + 2H₂O）。

- 速率 [mol/(m³·s)]:
    r = k_ch,R7 * C_CH4 * C_O2                                       … Eq. 5.40

- 速率常数:
    k_ch,R7 = (3.552e11 / T) * exp(-15700 / T)                       … Eq. 5.41
"""

from __future__ import annotations

from src.kinetics.arrhenius import k_jensen_r7

HAMEL_R7_CH4_A: float = 3.552e11
HAMEL_R7_CH4_E_OVER_RG_K: float = 15700.0


def hamel_r7_ch4_k(T_K: float) -> float:
    """Eq. 5.41 速率常数 k_ch,R7 [1/s]。"""
    return k_jensen_r7(HAMEL_R7_CH4_A, HAMEL_R7_CH4_E_OVER_RG_K, T_K)


def hamel_r7_ch4_rate_mol_m3_s(
    T_K: float,
    C_CH4_mol_m3: float,
    C_O2_mol_m3: float,
) -> float:
    """de Souza-Santos CH₄ 氧化体积速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.40
    """
    k = hamel_r7_ch4_k(T_K)
    return k * max(C_CH4_mol_m3, 0.0) * max(C_O2_mol_m3, 0.0)
