"""R10–R11 焦油反应（两组分代理化学计量 + 速率函数）。

R10: 焦油氧化（均相气相，气泡/悬浮相同动力学，Eq.5.59）
R11: 焦油转化（裂解/重整）— 气泡相 Serio 均相热裂解；悬浮相 Corella 催化（× rho_cat）

Source:
- Hamel (1999) §5.2.6, Eq.5.59, Table 5.4 (pp.95–98)
- specs/species.md §2–3
- docs/CLAUDE.md Phase 3.3
"""

from __future__ import annotations

from typing import Dict, Literal

import numpy as np

from src.core.species import (
    TarComponent,
    TarFuelType,
    TarSurrogate,
    calc_tar_surrogate_fractions,
)
from src.kinetics.arrhenius import k_standard

TarReactionId = Literal["R10", "R11"]

# 约定：反应物为负，产物为正；键为气体组分名（与 species.py 主气相一致）
TAR_STOICH_R10: Dict[TarSurrogate, Dict[str, float]] = {
    "C6H6": {"O2": -3.0, "CO": 6.0, "H2": 3.0},
    "C10H8": {"O2": -5.0, "CO": 10.0, "H2": 4.0},
    "C16H34": {"O2": -8.0, "CO": 16.0, "H2": 17.0},
}

TAR_STOICH_R11: Dict[TarSurrogate, Dict[str, float]] = {
    "C6H6": {"H2O": -4.5, "CO": 4.5, "H2": 4.5, "CH4": 1.5},
    "C10H8": {"H2O": -(23.0 / 3.0), "CO": 23.0 / 3.0, "H2": 7.0, "CH4": 7.0 / 3.0},
    "C16H34": {"H2O": -10.5, "CO": 10.5, "H2": 16.5, "CH4": 5.5},
}

_TAR_STOICH_TABLE: Dict[TarReactionId, Dict[TarSurrogate, Dict[str, float]]] = {
    "R10": TAR_STOICH_R10,
    "R11": TAR_STOICH_R11,
}


def get_tar_stoichiometry(reaction_id: TarReactionId, surrogate: TarSurrogate) -> Dict[str, float]:
    """获取单一 tar 代理分子的反应化学计量系数。"""
    return dict(_TAR_STOICH_TABLE[reaction_id][surrogate])


def get_lumped_tar_stoichiometry(
    reaction_id: TarReactionId,
    fuel_type: TarFuelType,
    target_hc_ratio: float | None = None,
) -> Dict[str, float]:
    """将两组分 tar 代理按 H/C 配比折算为 lumped tar 的等效化学计量。"""
    fractions = calc_tar_surrogate_fractions(fuel_type, target_hc_ratio=target_hc_ratio)
    table = _TAR_STOICH_TABLE[reaction_id]

    lumped: Dict[str, float] = {"TAR1": 0.0, "TAR2": 0.0}
    non_zero_surrogates = [(s, x) for s, x in fractions.items() if x > 0.0]
    if len(non_zero_surrogates) != 2:
        raise ValueError("两组分 tar 模型要求恰好两个非零代理组分")

    for comp, (surrogate, x_s) in zip(("TAR1", "TAR2"), non_zero_surrogates):
        lumped[comp] -= x_s
        if x_s <= 0.0:
            continue
        for species, nu in table[surrogate].items():
            lumped[species] = lumped.get(species, 0.0) + x_s * nu
    return lumped


def get_tar_component_stoichiometry(
    reaction_id: TarReactionId,
    fuel_type: TarFuelType,
    component: TarComponent,
    target_hc_ratio: float | None = None,
) -> Dict[str, float]:
    """返回单一 tar 组分（TAR1/TAR2）对应的化学计量（含自身消耗 -1）。"""
    all_component_stoich = get_lumped_tar_stoichiometry(
        reaction_id=reaction_id,
        fuel_type=fuel_type,
        target_hc_ratio=target_hc_ratio,
    )
    frac_component = -all_component_stoich[component]
    if frac_component <= 0.0:
        raise ValueError(f"{component} 分数无效: {frac_component}")

    stoich: Dict[str, float] = {component: -1.0}
    for species, coeff in all_component_stoich.items():
        if species in {"TAR1", "TAR2"}:
            continue
        stoich[species] = coeff / frac_component
    return stoich


# ---------------------------------------------------------------------------
# R10 Eq.5.59, Table 5.4 — 气泡/悬浮相同热力学参数
# R = k10·exp(-E_Rg/T)·T·P^0.3·C_tar^0.5·C_O2（E/Rg 为表中活化温度 [K]）
# ---------------------------------------------------------------------------
R10_k0_AROM: float = 20_700.0   # 芳香（C6H6, C10H8）；Siminski (1972) 经 Hamel
R10_E_Rg_AROM: float = 9_650.0  # [K]
R10_k0_OLEF: float = 59.8      # 烯烃/烷烃（如 C16H34）
R10_E_Rg_OLEF: float = 12_200.0  # [K]


from src.kinetics.arrhenius import k_hobbs
from src.core.constants import Rg

def _r10_single_class(
    T: float,
    C_tar: float,
    C_O2: float,
    P: float,
    k10: float,
    E_Rg: float,
) -> float:
    """单类 tar（芳香或烯烃/烷烃）的 R10 速率 [mol/(m³·s)]。"""
    C_tar = max(C_tar, 0.0)
    C_O2 = max(C_O2, 0.0)
    # k = (
    #     k10
    #     * np.exp(np.clip(-E_Rg / max(T, 300.0), -100.0, 100.0))
    #     * max(T, 300.0)
    #     * (P ** 0.3)
    # )
    # k_hobbs: k0 * T * exp(-E/(Rg*T)) -> E/Rg is input E_Rg
    k = k_hobbs(k10, E_Rg * Rg, max(T, 300.0)) * (P ** 0.3)
    return k * (C_tar ** 0.5) * C_O2


def rate_R10(
    T: float,
    C_tar: float,
    C_O2: float,
    P: float,
    fuel_type: TarFuelType,
) -> float:
    """R10 焦油氧化 [mol/(m³·s)]，两组分代理按摩尔分率加权（芳香 vs 烯烃/烷烃）。

    Eq.5.59 (p.96)；气泡与悬浮相使用同一组热力学参数。

    Source: Hamel (1999) §5.2.6, Table 5.4
    """
    frac = calc_tar_surrogate_fractions(fuel_type)
    r = 0.0
    for surr, w in frac.items():
        if w <= 0.0:
            continue
        arom = surr in ("C6H6", "C10H8")
        k10, e_rg = (R10_k0_AROM, R10_E_Rg_AROM) if arom else (R10_k0_OLEF, R10_E_Rg_OLEF)
        r += w * _r10_single_class(T, C_tar, C_O2, P, k10, e_rg)
    return r


# ---------------------------------------------------------------------------
# R11 Table 5.4 — 气泡：Serio (1987)；悬浮：Corella (1991) × rho_cat
# ---------------------------------------------------------------------------
# 气泡相：均相热裂解，一级对 C_tar
R11_SERIO_K0: float = 5.42e4    # [1/s]
R11_SERIO_E: float = 100_500.0  # [J/mol]

# 悬浮相：k0 [m³/(s·kg_Kat)]，须乘局部催化剂质量密度 rho_cat [kg/m³]
R11_CORELLA_K0: float = 0.7     # [m³/(s·kg_Kat)]
R11_CORELLA_E: float = 63_100.0  # [J/mol]


def rate_R11_bubble(T: float, C_tar: float) -> float:
    """R11 气泡相：Serio et al. (1987) 均相热裂解 [mol/(m³·s)]。

    R = k0·exp(-E/(Rg·T))·C_tar

    Source: Hamel (1999) Table 5.4 p.98
    """
    C_tar = max(C_tar, 0.0)
    k = k_standard(R11_SERIO_K0, R11_SERIO_E, T)
    return k * C_tar


def rate_R11_suspension(T: float, C_tar: float, rho_cat: float) -> float:
    """R11 悬浮相：Corella et al. (1991) 催化转化 [mol/(m³·s)]。

    R = k0·rho_cat·exp(-E/(Rg·T))·C_tar，其中 rho_cat 为床料/灰等效质量浓度 [kg/m³]。

    Source: Hamel (1999) §5.2.6, Table 5.4
    """
    C_tar = max(C_tar, 0.0)
    rc = max(rho_cat, 0.0)
    k0_eff = R11_CORELLA_K0 * rc
    k = k_standard(k0_eff, R11_CORELLA_E, T)
    return k * C_tar


def rate_R11(T: float, C_tar: float, C_H2O: float = 0.0) -> float:
    """向后兼容：等价于 `rate_R11_bubble`（C_H2O 忽略）。"""
    _ = C_H2O
    return rate_R11_bubble(T, C_tar)
