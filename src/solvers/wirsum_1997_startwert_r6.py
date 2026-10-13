"""同相 R6 Startwert：x0 上消掉同一相里的 CH₄ ∩ O₂。

Hamel §2.1 的 Startwertwahl，不是 F(x) 限幅，也不是工程 fastox。
计量 ``CH4 + 1.5 O2 → CO + 2 H2O``，按相独立。不做 R5 / R12 / R7 / R1，
也不改 residual。R6 产的 CO 若与剩余 O2 同相，由调用方再跑 R5 Startwert。

与 ``same_phase_r12_startwert``、``same_phase_r5_startwert``、
``project_bed_holdup_fast_oxidation_closure`` 分名。

Source: Hamel (1999) §2.1 Startwertwahl; Eq. 5.40–5.41
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import GAS_SPECIES_INDEX

_IDX_CH4 = GAS_SPECIES_INDEX["CH4"]
_IDX_O2 = GAS_SPECIES_INDEX["O2"]
_IDX_CO = GAS_SPECIES_INDEX["CO"]
_IDX_H2O = GAS_SPECIES_INDEX["H2O"]


def same_phase_r6_overlap_mol_s(N: npt.ArrayLike) -> float:
    """同相可发生的 R6 进度 ξ = min(N_CH4, N_O2 / 1.5) [mol/s]。"""
    n = np.asarray(N, dtype=np.float64).reshape(-1)
    n_ch4 = max(float(n[_IDX_CH4]), 0.0)  # [mol/s]
    n_o2 = max(float(n[_IDX_O2]), 0.0)  # [mol/s]
    if n_o2 <= 0.0:
        return 0.0
    return float(min(n_ch4, n_o2 / 1.5))


def apply_same_phase_r6_startwert(N: npt.NDArray[np.float64]) -> dict[str, float]:
    """就地按 R6 计量消掉同相 CH₄ ∩ O₂。只改 holdup，不改 F(x)。

    ξ = min(N_CH4, N_O2/1.5) [mol/s]：
    N_CH4 ← N_CH4 − ξ，N_O2 ← N_O2 − 1.5 ξ，
    N_CO ← N_CO + ξ，N_H2O ← N_H2O + 2 ξ。
    """
    n = np.asarray(N, dtype=np.float64)
    n_ch4 = max(float(n[_IDX_CH4]), 0.0)
    n_o2 = max(float(n[_IDX_O2]), 0.0)
    n_co = float(n[_IDX_CO])
    n_h2o = float(n[_IDX_H2O])
    xi = float(min(n_ch4, n_o2 / 1.5)) if n_o2 > 0.0 else 0.0  # [mol/s]
    n[_IDX_CH4] = n_ch4 - xi
    n[_IDX_O2] = n_o2 - 1.5 * xi
    n[_IDX_CO] = n_co + xi
    n[_IDX_H2O] = n_h2o + 2.0 * xi
    N[:] = n
    return {
        "xi6_mol_s": xi,
        "N_CH4_before": n_ch4,
        "N_O2_before": n_o2,
        "N_CH4_after": float(n[_IDX_CH4]),
        "N_O2_after": float(n[_IDX_O2]),
        "N_CO_after": float(n[_IDX_CO]),
        "N_H2O_after": float(n[_IDX_H2O]),
        "overlap_after_mol_s": same_phase_r6_overlap_mol_s(n),
    }


def apply_named_r6_startwert(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    """同相 R6 Startwert，并同步 FB/旋风气接口。不改 F(x) 定义。"""
    z = np.asarray(x, dtype=np.float64).copy()
    dense = apply_same_phase_r6_startwert(z[0:7])
    bubble = apply_same_phase_r6_startwert(z[7:14])
    n_tot = z[0:7] + z[7:14]
    z[16:23] = n_tot
    z[23:30] = n_tot
    z[30:37] = n_tot
    return {"x": z, "dense": dense, "bubble": bubble}
