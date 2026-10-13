"""Hamel (1999) Eq. 5.38–5.39 H₂ 氧化参考式（thesis R6）。

Source: Hamel (1999) §5.2.2 PDF 55 / 印刷 91；Jensen et al. (1995)

说明
----
仓库实现编号 ``R12`` 对应 thesis **R6**（H₂ + ½O₂ → H₂O）。

扫描原式（2026-09-14 二次诊断；同日 PDF 55 / 印刷 91 整页复核）：

- 速率 [mol/(m³·s)]:
    r = k_ch,R6 * C_H2^{3/2} * C_O2                                  … Eq. 5.38

- 速率常数:
    k_ch,R6 = 1.631e9 * T^{3/2} * exp(-(E/Rg)/T)                     … Eq. 5.39
    论文为 ``3.420``，按德文千分位记法读作 E/Rg = 3420 K（同 Eq. 5.41 的 15.700 → 15700）。

旧转录 ``r = k C_H2 C_O2^{1/2}`` 已撤回。温度幂保持正的 ``3/2``，不得改成负幂。
"""

from __future__ import annotations

from src.kinetics.arrhenius import k_jensen_h2

HAMEL_R6_H2_A: float = 1.631e9
HAMEL_R6_H2_E_OVER_RG_K: float = 3420.0
HAMEL_R6_H2_H2_ORDER: float = 1.5
HAMEL_R6_H2_O2_ORDER: float = 1.0


def hamel_r6_h2_k(T_K: float) -> float:
    """Eq. 5.39 速率常数 k_ch,M [ (m³/mol)^{3/2} / s ]。"""
    return k_jensen_h2(HAMEL_R6_H2_A, HAMEL_R6_H2_E_OVER_RG_K, T_K)


def hamel_r6_h2_rate_mol_m3_s(
    T_K: float,
    C_H2_mol_m3: float,
    C_O2_mol_m3: float,
) -> float:
    """Jensen H₂ 氧化体积速率 [mol/(m³·s)]。

    Ref: Hamel (1999) Eq. 5.38
    """
    k = hamel_r6_h2_k(T_K)
    c_h2 = max(C_H2_mol_m3, 0.0)
    c_o2 = max(C_O2_mol_m3, 0.0)
    return k * (c_h2 ** HAMEL_R6_H2_H2_ORDER) * (c_o2 ** HAMEL_R6_H2_O2_ORDER)
