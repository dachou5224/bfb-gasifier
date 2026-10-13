"""同相 R10 Startwert：共享焦油池上消掉各相 O₂ ∩ tar。

Hamel §2.1 的 Startwertwahl，不是 F(x) 限幅，也不是工程 fastox。
计量用生产 lumped ``get_lumped_tar_stoichiometry('R10')``。
先密相后气泡，共用 ``N_TAR1+N_TAR2``。不改 residual。

只消 holdup 重叠。气泡入口仍供氧时，O₂=0 且 H₂/CO>0 不是 s=1 盆地；
具名 R10 延拓默认从 E6 根走，不把本核当 s=1 初值。

Source: Hamel (1999) §2.1 Startwertwahl; Eq. 5.59
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
    return get_lumped_tar_stoichiometry("R10", R10_FUEL_TYPE)


def apply_same_phase_r10_startwert(
    N: npt.NDArray[np.float64],
    n_tar1: float,
    n_tar2: float,
    sto: dict[str, float],
) -> dict[str, float]:
    """就地按 R10 计量消掉本相 O₂ 与共享焦油。返回更新后的焦油 [mol/s]。"""
    n = np.asarray(N, dtype=np.float64)
    nu_o2 = -float(sto.get("O2", 0.0))
    n_o2 = max(float(n[IDX["O2"]]), 0.0)
    n_t = max(float(n_tar1) + float(n_tar2), 0.0)
    xi = 0.0
    if n_o2 > 0.0 and n_t > 0.0 and nu_o2 > 0.0:
        xi = float(min(n_t, n_o2 / nu_o2))
    for name, nu in sto.items():
        if name in IDX:
            n[IDX[name]] = float(n[IDX[name]]) + float(nu) * xi
    n[IDX["O2"]] = max(float(n[IDX["O2"]]), 0.0)
    n_t1 = max(float(n_tar1) + float(sto.get("TAR1", 0.0)) * xi, 0.0)
    n_t2 = max(float(n_tar2) + float(sto.get("TAR2", 0.0)) * xi, 0.0)
    N[:] = n
    return {
        "xi10_mol_s": xi,
        "N_O2_before": n_o2,
        "N_O2_after": max(float(n[IDX["O2"]]), 0.0),
        "N_tar1_after": max(n_t1, 0.0),
        "N_tar2_after": max(n_t2, 0.0),
        "N_CO_after": float(n[IDX["CO"]]),
        "N_H2_after": float(n[IDX["H2"]]),
    }


def apply_named_r10_startwert(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    """共享焦油池：密相然后气泡。同步 FB/旋风气接口。不改 F(x)。"""
    z = np.asarray(x, dtype=np.float64).copy()
    sto = _stoich()
    n1 = max(float(z[_IDX_TAR1]), 0.0)
    n2 = max(float(z[_IDX_TAR2]), 0.0)
    dense = apply_same_phase_r10_startwert(z[0:7], n1, n2, sto)
    n1 = float(dense["N_tar1_after"])
    n2 = float(dense["N_tar2_after"])
    bubble = apply_same_phase_r10_startwert(z[7:14], n1, n2, sto)
    z[_IDX_TAR1] = float(bubble["N_tar1_after"])
    z[_IDX_TAR2] = float(bubble["N_tar2_after"])
    n_tot = z[0:7] + z[7:14]
    z[16:23] = n_tot
    z[23:30] = n_tot
    z[30:37] = n_tot
    return {"x": z, "dense": dense, "bubble": bubble}
