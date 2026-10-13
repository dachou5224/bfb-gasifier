"""R11 Startwert：计量消耗 holdup，不改 residual。

全量同相 ``H2O ∩ tar`` 会把 H2 写入 ``N_d``，撞上 R8/R12。延拓默认用
``apply_scaled_bubble_r11_predictor``：只在气泡相按 ``ξ=Δs·N_tar`` 预消耗，
多出的 H2 走反向 WGSR。``apply_named_r11_startwert`` 仍保留全量 holdup 形式。

Source: Hamel (1999) §2.1 Startwertwahl; Table 5.4; Eq. 5.45–5.46
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry
from src.solvers.wirsum_1997_named_f import R10_FUEL_TYPE
from src.solvers.wirsum_1997_stiffness_full_gas import IDX

_IDX_TAR1: int = 40
_IDX_TAR2: int = 41


def _stoich() -> dict[str, float]:
    return get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)


def apply_same_phase_r11_startwert(
    N: npt.NDArray[np.float64],
    n_tar1: float,
    n_tar2: float,
    sto: dict[str, float],
) -> dict[str, float]:
    """就地按 R11 计量消掉本相 H2O 与共享焦油。"""
    n = np.asarray(N, dtype=np.float64)
    nu_h2o = -float(sto.get("H2O", 0.0))
    n_h2o = max(float(n[IDX["H2O"]]), 0.0)
    n_t = max(float(n_tar1) + float(n_tar2), 0.0)
    xi = 0.0
    if n_h2o > 0.0 and n_t > 0.0 and nu_h2o > 0.0:
        xi = float(min(n_t, n_h2o / nu_h2o))
    for name, nu in sto.items():
        if name in IDX:
            n[IDX[name]] = float(n[IDX[name]]) + float(nu) * xi
    n[IDX["H2O"]] = max(float(n[IDX["H2O"]]), 0.0)
    n_t1 = max(float(n_tar1) + float(sto.get("TAR1", 0.0)) * xi, 0.0)
    n_t2 = max(float(n_tar2) + float(sto.get("TAR2", 0.0)) * xi, 0.0)
    N[:] = n
    return {
        "xi11_mol_s": xi,
        "N_H2O_before": n_h2o,
        "N_H2O_after": max(float(n[IDX["H2O"]]), 0.0),
        "N_tar1_after": max(n_t1, 0.0),
        "N_tar2_after": max(n_t2, 0.0),
        "N_CO_after": float(n[IDX["CO"]]),
        "N_H2_after": float(n[IDX["H2"]]),
        "N_CH4_after": float(n[IDX["CH4"]]),
    }


def _sync_named_gas_interfaces(z: npt.NDArray[np.float64]) -> None:
    n_tot = z[0:7] + z[7:14]
    z[16:23] = n_tot
    z[23:30] = n_tot
    z[30:37] = n_tot


def apply_scaled_bubble_r11_predictor(
    x: npt.NDArray[np.float64],
    *,
    ds: float,
) -> dict[str, Any]:
    """按 Δs 在气泡相预消耗 ξ=ds·N_tar，多出的 H2 走反向 WGSR。不改 F(x)。

    全量 holdup Startwert 会把 H2 写入 N_d，撞上 R8/R12。气泡 Serio 才是
    硬源；密相 Corella 留给 Newton。

    Source: Hamel (1999) §2.1 Startwertwahl; Table 5.4; Eq. 5.45–5.46
    """
    z = np.asarray(x, dtype=np.float64).copy()
    sto = _stoich()
    ds_use = float(np.clip(ds, 0.0, 1.0))
    n1 = max(float(z[_IDX_TAR1]), 0.0)
    n2 = max(float(z[_IDX_TAR2]), 0.0)
    tot = n1 + n2
    xi = ds_use * tot
    nu_h2o = -float(sto.get("H2O", 0.0))
    n_h2o_b = max(float(z[7 + IDX["H2O"]]), 0.0)
    if nu_h2o > 0.0:
        xi = min(xi, n_h2o_b / nu_h2o)
    Nb = z[7:14].copy()
    for name, nu in sto.items():
        if name in IDX:
            Nb[IDX[name]] = float(Nb[IDX[name]]) + float(nu) * xi
    Nb = np.maximum(Nb, 0.0)
    z[7:14] = Nb
    if tot > 0.0:
        frac = xi / tot
        z[_IDX_TAR1] = n1 * (1.0 - frac)
        z[_IDX_TAR2] = n2 * (1.0 - frac)
    xi8 = max(float(sto.get("H2", 0.0)) * xi, 0.0)
    n_co2 = float(z[7 + IDX["CO2"]])
    wgsr = False
    if xi8 > 0.0 and n_co2 >= xi8:
        z[7 + IDX["H2"]] = float(np.asarray(x, dtype=np.float64)[7 + IDX["H2"]])
        z[7 + IDX["CO"]] = float(z[7 + IDX["CO"]]) + xi8
        z[7 + IDX["CO2"]] = float(z[7 + IDX["CO2"]]) - xi8
        z[7 + IDX["H2O"]] = float(z[7 + IDX["H2O"]]) + xi8
        z[7:14] = np.maximum(z[7:14], 0.0)
        wgsr = True
    _sync_named_gas_interfaces(z)
    return {
        "x": z,
        "xi11_mol_s": float(xi),
        "ds": ds_use,
        "reverse_wgsr": wgsr,
        "N_tar_after": float(z[_IDX_TAR1] + z[_IDX_TAR2]),
        "N_b_H2O_after": float(z[7 + IDX["H2O"]]),
    }


def apply_named_r11_startwert(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    """共享焦油池：密相然后气泡。同步 FB/旋风气接口。不改 F(x)。"""
    z = np.asarray(x, dtype=np.float64).copy()
    sto = _stoich()
    n1 = max(float(z[_IDX_TAR1]), 0.0)
    n2 = max(float(z[_IDX_TAR2]), 0.0)
    dense = apply_same_phase_r11_startwert(z[0:7], n1, n2, sto)
    n1 = float(dense["N_tar1_after"])
    n2 = float(dense["N_tar2_after"])
    bubble = apply_same_phase_r11_startwert(z[7:14], n1, n2, sto)
    z[_IDX_TAR1] = float(bubble["N_tar1_after"])
    z[_IDX_TAR2] = float(bubble["N_tar2_after"])
    _sync_named_gas_interfaces(z)
    return {"x": z, "dense": dense, "bubble": bubble}
