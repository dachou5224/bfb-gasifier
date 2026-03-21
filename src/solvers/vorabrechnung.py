"""Vorabrechnung（预算）：在全局 NR 前生成与轴向温度一致的化学初值 x₀。

对应 Hamel (1999) Module 2 / §2.2 思想：先估计轴向温度与热解进度，再构造元素大致平衡的
气相摩尔流初值，降低 F(x₀) 与内层 NR 难度。

注：此为工程近似，非论文 Fortran 逐行复现。
"""

from __future__ import annotations

from typing import List

import numpy as np

from src.core.cell import Cell
from src.core.species import GAS_SPECIES_INDEX
from src.thermal.devolatilization import devolatilization_rate_for_cell


def estimate_axial_T_profile(
    n_cells: int,
    T_inlet: float,
    O2_feed: float,
    H2O_feed: float,
    N2_feed: float,
    fuel_feed_kg_s: float,
    C_dry: float,
    H_dry: float,
    moisture_wt: float,
    P: float,
) -> np.ndarray:
    """估计轴向温度剖面（底高顶低），用于 Vorabrechnung 初值。

    简化绝热温升 + 线性轴向衰减；顶温做下限使剖面落在常见床层/出口区间。
    """
    moisture_frac = moisture_wt / 100.0
    dry_fuel = fuel_feed_kg_s * (1.0 - moisture_frac)

    M_C = 12.011e-3
    n_C = dry_fuel * (C_dry / 100.0) / M_C

    dH_comb = 393_500.0  # J/mol O2（量级近似）
    total_mol_s = O2_feed + H2O_feed + N2_feed + max(n_C * 0.5, 1e-6)
    Cp_mix = 32.0  # J/(mol·K)

    dT_adiabatic = (O2_feed * dH_comb) / max(total_mol_s * Cp_mix, 1.0)
    T_bottom = float(np.clip(T_inlet + dT_adiabatic, 1050.0, 1495.0))
    T_top = float(np.clip(T_bottom * 0.82, 900.0, 1200.0))

    if n_cells <= 1:
        return np.array([T_bottom], dtype=np.float64)

    T_profile = np.array(
        [
            T_bottom + (T_top - T_bottom) * (i / max(n_cells - 1, 1))
            for i in range(n_cells)
        ],
        dtype=np.float64,
    )
    return T_profile


def generate_initial_x0(
    cells: List[Cell],
    O2_feed: float,
    H2O_feed: float,
    N2_feed: float,
    fuel_feed_kg_s: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    moisture_wt: float,
    ash_dry_wt: float,
    VM_daf: float,
    T_profile: np.ndarray,
    fuel_type: str = "coal",
) -> None:
    """按轴向高度与估计热解进度，写入各 cell 的 N_d/N_b/T/m_solid 初值（就地修改）。"""
    n = len(cells)
    idx = GAS_SPECIES_INDEX

    total_inlet = O2_feed + H2O_feed + N2_feed

    moisture_frac = moisture_wt / 100.0
    ash_frac = ash_dry_wt / 100.0
    vm_daf_frac = VM_daf / 100.0
    dry_feed = fuel_feed_kg_s * (1.0 - moisture_frac)
    daf_feed = dry_feed * (1.0 - ash_frac)

    M_C, M_H, M_O = 12.011e-3, 1.00794e-3, 15.999e-3
    c_dry = C_dry / 100.0
    h_dry = H_dry / 100.0
    o_dry = O_dry / 100.0
    to_daf = 1.0 / max(1.0 - ash_frac, 1e-9)
    c_daf = c_dry * to_daf
    h_daf = h_dry * to_daf
    o_daf = o_dry * to_daf

    daem_fuel = "brown_coal"
    if fuel_type in ("wood", "biomass"):
        daem_fuel = "wood"
    elif fuel_type in ("coal", "brown_coal"):
        daem_fuel = "brown_coal"

    for i, cell in enumerate(cells):
        T_i = float(T_profile[i])
        frac_height = (i + 0.5) / max(n, 1)

        o2_consumed_frac = min(1.0 - np.exp(-5.0 * frac_height), 1.0)
        o2_remaining = O2_feed * max(1.0 - o2_consumed_frac, 0.0)

        tau_est = cell.geo.dh / max(0.05, cell.u_mf if cell.u_mf > 0 else 0.05)
        X_vm, _ = devolatilization_rate_for_cell(
            T_bed=T_i,
            tau_cell=tau_est,
            T_init=max(float(T_profile[0]) * 0.3, 400.0),
            VM_daf=vm_daf_frac,
            fuel_type=daem_fuel,
        )

        X_vm_cum = min(X_vm * (frac_height + 0.1), 1.0)
        m_vm_released = daf_feed * vm_daf_frac * X_vm_cum

        nC_vm = m_vm_released * c_daf / M_C
        nH_vm = m_vm_released * h_daf / M_H / 2.0
        nO_vm = m_vm_released * o_daf / M_O / 2.0

        co2_frac = max(0.0, 1.0 - frac_height)
        co_frac = min(1.0, frac_height + 0.2)

        n_CO2 = nC_vm * co2_frac * 0.5
        n_CO = nC_vm * co_frac * 0.5
        n_CH4 = nC_vm * 0.05
        n_H2 = max(nH_vm - 2.0 * (nC_vm * co2_frac * 0.1), 0.01 * max(nH_vm, 1e-12))
        n_H2O = H2O_feed + (O2_feed * o2_consumed_frac * 0.5) - nO_vm * 0.5
        n_H2O = max(n_H2O, 0.01 * H2O_feed)
        n_N2 = N2_feed
        n_O2 = o2_remaining

        total = n_CO2 + n_CO + n_CH4 + n_H2 + n_H2O + n_N2 + n_O2
        total = max(total, total_inlet * 0.5)

        cell.N_d.fill(0.0)
        cell.N_b.fill(0.0)
        for sp, val in [
            ("CO2", n_CO2),
            ("CO", n_CO),
            ("CH4", n_CH4),
            ("H2", n_H2),
            ("H2O", n_H2O),
            ("N2", n_N2),
            ("O2", n_O2),
        ]:
            cell.N_d[idx[sp]] = max(val * 0.7, 1e-12)
            cell.N_b[idx[sp]] = max(val * 0.3, 1e-12)

        cell.T = T_i
        solid_frac = 0.3 + 0.7 * (1.0 - frac_height)
        nk = cell.solid.n_size_classes
        for k in range(nk):
            cell.m_solid[k] = max(fuel_feed_kg_s * solid_frac * 0.5 / max(nk, 1), 1e-12)
