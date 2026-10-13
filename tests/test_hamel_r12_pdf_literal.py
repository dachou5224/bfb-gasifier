"""Hamel Eq. 5.38–5.39：期望值按扫描原式手写，不调用生产参考函数。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.kinetics.gas_reactions import rate_R12

# Hamel (1999) §5.2.2 PDF 55 / 印刷 91 扫描式
_A = 1.631e9
_E_OVER_RG_K = 3420.0
_H2_ORDER = 1.5
_O2_ORDER = 1.0


def _pdf_literal_rate(T: float, C_H2: float, C_O2: float) -> float:
    t = max(float(T), 300.0)
    k = _A * (t**1.5) * math.exp(-_E_OVER_RG_K / t)
    return k * (max(C_H2, 0.0) ** _H2_ORDER) * (max(C_O2, 0.0) ** _O2_ORDER)


def test_rate_r12_matches_pdf_literal_not_old_half_order() -> None:
    T = 1325.0
    C_h2, C_o2 = 3.890713377968589, 18.54460093596481
    y = np.zeros(N_GAS)
    got = rate_R12(T, C_h2, C_o2, 2.5e6, y, GAS_SPECIES_INDEX)
    literal = _pdf_literal_rate(T, C_h2, C_o2)
    old_wrong = (
        _A
        * (T**1.5)
        * math.exp(-_E_OVER_RG_K / T)
        * C_h2
        * (C_o2**0.5)
    )
    assert got == pytest.approx(literal, rel=1e-12)
    assert got / old_wrong == pytest.approx(literal / old_wrong, rel=1e-12)
    assert got / old_wrong == pytest.approx(8.494, rel=2e-3)


def test_r12_oxygen_derivative_finite_at_trace_o2() -> None:
    T = 1325.0
    C_h2 = 3.890713377968589
    y = np.zeros(N_GAS)
    c1, c2 = 1.0e-8, 1.01e-8
    r1 = rate_R12(T, C_h2, c1, 2.5e6, y, GAS_SPECIES_INDEX)
    r2 = rate_R12(T, C_h2, c2, 2.5e6, y, GAS_SPECIES_INDEX)
    dr_dco2 = (r2 - r1) / (c2 - c1)
    assert math.isfinite(dr_dco2)
    expected_k = _A * (T**1.5) * math.exp(-_E_OVER_RG_K / T) * (C_h2**_H2_ORDER)
    assert dr_dco2 == pytest.approx(expected_k, rel=1e-3)
