"""自由板混合气 Startwert：只改 ``N_fb`` 上的燃料 ∩ O₂。

床层两相在 F_E6+R10+R11 根上可以各自无氧，混合成 ``N_fb=N_d+N_b``
后仍可能出现 H₂/CO/CH₄ ∩ O₂。这是 Hamel §2.1 的 Startwertwahl，
不是 F(x) 限幅，也不是生产 ``_freeboard_reaction_step``。

顺序：R12 → R5 → R6 → 再一刀 R5（R6 产 CO 若仍与 O₂ 重叠）。
同步旋风气拷贝 ``N_c1``/``N_c2``。不改床层 ``N_d``/``N_b``，不改 residual。

Source: Hamel (1999) §2.1 Startwertwahl; Eq. 5.36–5.41
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_startwert_r5 import apply_same_phase_r5_startwert
from src.solvers.wirsum_1997_startwert_r6 import apply_same_phase_r6_startwert
from src.solvers.wirsum_1997_startwert_r12 import apply_same_phase_r12_startwert


def apply_named_e7_fb_startwert(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    """只在 ``N_fb`` 上消氧化重叠，再拷到旋风气。不改 F(x)。"""
    z = np.asarray(x, dtype=np.float64).copy()
    n_fb = z[16:23].copy()
    r12 = apply_same_phase_r12_startwert(n_fb)
    r5_a = apply_same_phase_r5_startwert(n_fb)
    r6 = apply_same_phase_r6_startwert(n_fb)
    r5_b = apply_same_phase_r5_startwert(n_fb)
    z[16:23] = n_fb
    z[23:30] = n_fb
    z[30:37] = n_fb
    return {
        "x": z,
        "r12": r12,
        "r5_first": r5_a,
        "r6": r6,
        "r5_after_r6": r5_b,
        "xi12_mol_s": float(r12["xi12_mol_s"]),
        "xi5_mol_s": float(r5_a["xi5_mol_s"]) + float(r5_b["xi5_mol_s"]),
        "xi6_mol_s": float(r6["xi6_mol_s"]),
    }
