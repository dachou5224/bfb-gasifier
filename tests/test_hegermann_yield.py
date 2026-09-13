"""Hegermann Eq.4.21–4.22 与 tar 代理 atom-H/C 同构测试。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.species import (
    TAR_SURROGATE_FORMULA,
    calc_tar_surrogate_fractions,
    fuel_tar_target_hc_ratio,
)
from src.thermal.hegermann_yield import (
    tar_carbon_frac_of_volatile_c,
    xc_koks_hegermann,
    xc_teer_hegermann,
)


def test_hegermann_xc_trends_match_bild_4_8_4_9() -> None:
    """大气压高温：Koks~0.55、Teer~0.25；升压 Teer 降、Koks 升。"""
    k_lo = xc_koks_hegermann(1.0e5, 1173.0)
    t_lo = xc_teer_hegermann(1.0e5, 1173.0)
    k_hi = xc_koks_hegermann(2.5e6, 1173.0)
    t_hi = xc_teer_hegermann(2.5e6, 1173.0)
    assert 0.50 < k_lo < 0.60
    assert 0.20 < t_lo < 0.35
    assert k_hi > k_lo
    assert t_hi < t_lo
    assert 0.08 < t_hi < 0.20


def test_f_tar_of_volatile_equals_xc_teer_over_one_minus_koks() -> None:
    P, T = 2.5e6, 1173.0
    f = tar_carbon_frac_of_volatile_c(P, T)
    expected = xc_teer_hegermann(P, T) / max(1.0 - xc_koks_hegermann(P, T), 1e-12)
    assert f == pytest.approx(min(max(expected, 0.0), 0.9), rel=0.0, abs=1e-12)


def test_lu_fuel_hc_is_naphthalene_edge() -> None:
    """Wesseling LU C=61.5/H=4.1 在 C=12、H=1 下 H/C=0.8 → 纯萘。"""
    hc = fuel_tar_target_hc_ratio(61.5, 4.1)
    assert hc == pytest.approx(0.8, abs=1e-12)
    frac = calc_tar_surrogate_fractions("coal", target_hc_ratio=hc)
    assert frac["C6H6"] == pytest.approx(0.0, abs=1e-12)
    assert frac["C10H8"] == pytest.approx(1.0, abs=1e-12)


def test_table61_anchors_atom_hc() -> None:
    frac = calc_tar_surrogate_fractions("coal", target_hc_ratio=0.856)
    n_h = sum(frac[s] * TAR_SURROGATE_FORMULA[s][1] for s in frac)
    n_c = sum(frac[s] * TAR_SURROGATE_FORMULA[s][0] for s in frac)
    assert n_h / n_c == pytest.approx(0.856, abs=1e-12)
    assert frac["C6H6"] == pytest.approx(0.3932584269662921, rel=1e-6)
