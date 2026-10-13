"""Audit/benchmark 脚本公共辅助函数（Tier M 合并产物）。

仅包含逐 cell 分析前准备与炭反应速率计算等无状态纯函数，
避免 `audit_bed_co_co2_budget_lu.py` 与 `audit_bed_exit_to_reactor_exit_budget_lu.py`
重复相同样板。
"""

from __future__ import annotations

from typing import Any

from src.core.constants import Rg
from src.kinetics.char_reactions import (
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)


def prepare_cell_for_analysis(cell) -> None:
    """计算 cell 的流体力学、vorabrechnung 与相间交换（audit 逐 cell 分析前统一准备）。"""
    cell.calc_hydrodynamics()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau_val)
    cell.calc_exchange()


def compute_char_reaction_rates(cell, D_g: float) -> dict[str, Any]:
    """按 cell 当前状态计算炭异相反应 R1–R4 速率与相关辅助量。

    Returns
    -------
    dict with keys: r1, r2, r3, r4, alpha, total_area
    """
    areas = cell._calc_char_surface_area_per_class()
    total_area = float(__import__("numpy", fromlist=["sum"]).sum(areas))
    char_conversion = cell._compute_char_conversion()
    T = cell.T
    C_d = cell._concentrations("d")
    idx = __import__("src.core.species", fromlist=["GAS_SPECIES_INDEX"]).GAS_SPECIES_INDEX

    r1 = r2 = r3 = r4 = 0.0
    alpha = 1.0
    if total_area > 0.0:
        dc = d_core_from_spm_char_conversion(float(char_conversion), float(cell.solid.d_p))
        r1, alpha = rate_R1(T, C_d[idx["O2"]], float(cell.solid.d_p), D_g, dc, cell.fuel_type)
        r2 = rate_R2(T, C_d[idx["H2O"]], float(cell.solid.d_p), D_g, dc)
        r3 = rate_R3(T, C_d[idx["H2"]], float(cell.solid.d_p), D_g, dc)
        p_co2 = max(float(C_d[idx["CO2"]] * Rg * T), 0.0)
        p_co = max(float(C_d[idx["CO"]] * Rg * T), 0.0)
        r4 = cell.r4_scale * rate_R4_effective(T, p_co2, p_co, float(cell.solid.d_p), D_g, dc)

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "r4": r4,
        "alpha": alpha,
        "total_area": total_area,
    }
