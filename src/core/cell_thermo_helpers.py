"""Cell 热力学辅助：气相摩尔分数、浓度与焓流（与 cell_balances 一致）。

Source: ``docs/hamel_submodels/02_cell_balances_and_exchange.md``;
Hamel (1999) Eq. 2.5, 2.7
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from src.core.cell_balances import calc_gas_enthalpy_flow, calc_solid_enthalpy_flow
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX, N_GAS, cp_ash, cp_char, cp_sand

__all__ = [
    "calc_gas_enthalpy_flow",
    "calc_solid_enthalpy_flow",
    "gas_mole_fractions",
    "gas_concentrations",
    "solid_enthalpy_flow_for_fuel_props",
]


def gas_mole_fractions(
    N_b: npt.NDArray[np.float64],
    N_d: npt.NDArray[np.float64],
    phase: str,
) -> npt.NDArray[np.float64]:
    """气相摩尔分数 y_i。

    phase
        ``"b"`` / ``"d"``：单相；``"combined"`` / ``"total"`` / ``"mixed"``：
        气泡相与悬浮相摩尔流加权和。
    """
    if phase in ("combined", "total", "mixed"):
        work = np.maximum(N_b + N_d, 0.0)
    else:
        N = N_b if phase == "b" else N_d
        work = np.maximum(N, 0.0)
    total = np.sum(work)
    if total < 1e-12:
        out = np.zeros(N_GAS)
        out[GAS_SPECIES_INDEX["N2"]] = 1.0
        return out
    return work / total


def gas_concentrations(
    N_b: npt.NDArray[np.float64],
    N_d: npt.NDArray[np.float64],
    phase: str,
    P: float,  # [Pa]
    T: float,  # [K]
) -> npt.NDArray[np.float64]:
    """理想气体摩尔浓度 C_i = y_i * P / (R T)。单位：[mol/m^3]。"""
    return gas_mole_fractions(N_b, N_d, phase) * (float(P) / (Rg * float(T)))


def solid_enthalpy_flow_for_fuel_props(
    m: npt.ArrayLike,
    T: float,
    *,
    ash_dry_wt: float,
    VM_daf: float,
    h_f_dry: float,
) -> float:
    """固相总焓流（与 ``Cell._calc_solid_enthalpy_flow`` 等价，供薄委托）。"""
    return calc_solid_enthalpy_flow(
        m,
        T,
        ash_dry_wt=ash_dry_wt,
        VM_daf=VM_daf,
        h_f_dry=h_f_dry,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
    )
