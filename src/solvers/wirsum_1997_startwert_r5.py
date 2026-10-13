"""同相 R5 Startwert：x0 上消掉同一相里的 CO ∩ O₂。

Hamel §2.1 的 Startwertwahl，不是 F(x) 限幅，也不是工程 fastox。
计量 ``CO + 1/2 O2 → CO2``，按相独立。不做 R6 / R12 / R1，也不改 residual。

与 ``same_phase_r12_startwert``、``project_bed_holdup_fast_oxidation_closure`` 分名。

Source: Hamel (1999) §2.1 Startwertwahl; Table 5.3 / Eq. 5.36
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import GAS_SPECIES_INDEX

_IDX_CO = GAS_SPECIES_INDEX["CO"]
_IDX_O2 = GAS_SPECIES_INDEX["O2"]
_IDX_CO2 = GAS_SPECIES_INDEX["CO2"]


def same_phase_r5_overlap_mol_s(N: npt.ArrayLike) -> float:
    """同相可发生的 R5 进度 ξ = min(N_CO, 2 N_O2) [mol/s]。"""
    n = np.asarray(N, dtype=np.float64).reshape(-1)
    n_co = max(float(n[_IDX_CO]), 0.0)  # [mol/s]
    n_o2 = max(float(n[_IDX_O2]), 0.0)  # [mol/s]
    return float(min(n_co, 2.0 * n_o2))


def apply_same_phase_r5_startwert(N: npt.NDArray[np.float64]) -> dict[str, float]:
    """就地按 R5 计量消掉同相 CO ∩ O₂。只改 holdup，不改 F(x)。

    ξ = min(N_CO, 2 N_O2) [mol/s]：
    N_CO ← N_CO − ξ，N_O2 ← N_O2 − 0.5 ξ，N_CO2 ← N_CO2 + ξ。
    """
    n = np.asarray(N, dtype=np.float64)
    n_co = max(float(n[_IDX_CO]), 0.0)
    n_o2 = max(float(n[_IDX_O2]), 0.0)
    n_co2 = float(n[_IDX_CO2])
    xi = float(min(n_co, 2.0 * n_o2))  # [mol/s]
    n[_IDX_CO] = n_co - xi
    n[_IDX_O2] = n_o2 - 0.5 * xi
    n[_IDX_CO2] = n_co2 + xi
    N[:] = n
    return {
        "xi5_mol_s": xi,
        "N_CO_before": n_co,
        "N_O2_before": n_o2,
        "N_CO_after": float(n[_IDX_CO]),
        "N_O2_after": float(n[_IDX_O2]),
        "overlap_after_mol_s": same_phase_r5_overlap_mol_s(n),
    }


def apply_named_r5_startwert(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    """同相 R5 Startwert，并同步 FB/旋风气接口。不改 F(x) 定义。"""
    z = np.asarray(x, dtype=np.float64).copy()
    dense = apply_same_phase_r5_startwert(z[0:7])
    bubble = apply_same_phase_r5_startwert(z[7:14])
    n_tot = z[0:7] + z[7:14]
    z[16:23] = n_tot
    z[23:30] = n_tot
    z[30:37] = n_tot
    return {"x": z, "dense": dense, "bubble": bubble}
