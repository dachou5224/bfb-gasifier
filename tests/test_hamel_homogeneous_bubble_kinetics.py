"""Hamel Eq. 5.36–5.41 气泡/均相动力学与 ``gas_reactions`` 对照。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.kinetics import gas_reactions as gr
from src.kinetics.gas_reactions import rate_R12, rate_R6
from src.kinetics.hamel_r6_h2_reference import (
    HAMEL_R6_H2_A,
    HAMEL_R6_H2_E_OVER_RG_K,
    HAMEL_R6_H2_H2_ORDER,
    HAMEL_R6_H2_O2_ORDER,
    hamel_r6_h2_rate_mol_m3_s,
)
from src.kinetics.hamel_r7_ch4_reference import (
    HAMEL_R7_CH4_A,
    HAMEL_R7_CH4_E_OVER_RG_K,
    hamel_r7_ch4_rate_mol_m3_s,
)


class TestHamelR6H2Implementation:
    def test_r12_constants_match_eq_539(self):
        assert gr.R12_A == pytest.approx(HAMEL_R6_H2_A)
        assert gr.R12_E_Rg == pytest.approx(HAMEL_R6_H2_E_OVER_RG_K)
        assert gr.R12_h2_order == pytest.approx(HAMEL_R6_H2_H2_ORDER)
        assert gr.R12_o2_order == pytest.approx(HAMEL_R6_H2_O2_ORDER)
        assert gr.R12_h2_order == pytest.approx(1.5)
        assert gr.R12_o2_order == pytest.approx(1.0)

    def test_r12_rate_matches_eq_538(self):
        T = 1325.0
        C_h2, C_o2 = 0.3, 0.01
        y = np.zeros(N_GAS)
        r_impl = rate_R12(T, C_h2, C_o2, 2.5e6, y, GAS_SPECIES_INDEX)
        r_ref = hamel_r6_h2_rate_mol_m3_s(T, C_h2, C_o2)
        assert r_impl == pytest.approx(r_ref, rel=1e-12)

    def test_r12_not_legacy_k_standard_form(self):
        """旧实现 k_standard(2e7, 15000*Rg) 与 Jensen T^1.5 形式应显著不同。"""
        from src.kinetics.arrhenius import k_jensen_h2, k_standard
        from src.core.constants import Rg

        T = 1325.0
        k_legacy = k_standard(2.0e7, 15_000.0 * Rg, T)
        k_hamel = k_jensen_h2(HAMEL_R6_H2_A, HAMEL_R6_H2_E_OVER_RG_K, T)
        assert k_hamel / k_legacy > 1.0e4


class TestHamelR7Ch4Implementation:
    def test_r6_constants_match_eq_541(self):
        assert gr.R6_k0 == pytest.approx(HAMEL_R7_CH4_A)
        assert gr.R6_E_Rg == pytest.approx(HAMEL_R7_CH4_E_OVER_RG_K)

    def test_r6_rate_matches_eq_540(self):
        T = 1200.0
        C_ch4, C_o2 = 2.0, 0.02
        r_impl = rate_R6(T, C_ch4, C_o2)
        r_ref = hamel_r7_ch4_rate_mol_m3_s(T, C_ch4, C_o2)
        assert r_impl == pytest.approx(r_ref, rel=1e-12)


class TestBubbleO2Competition:
    """bed0 量级：有限 O₂ 在 R5 / R12 / R6 之间的分账（仅动力学，无 r5_scale）。"""

    def test_h2_can_compete_when_present(self):
        from src.kinetics.gas_reactions import rate_R5_bubble

        T = 1325.0
        P = 2.5e6
        C_o2 = 0.01
        C_h2o = 15.0
        y = np.zeros(N_GAS)

        # 富 H2、低 CO 热解气
        C_co, C_h2, C_ch4 = 0.5, 15.0, 1.0
        r5 = rate_R5_bubble(T, C_co, C_o2, C_h2o, P, y, GAS_SPECIES_INDEX)
        r12 = rate_R12(T, C_h2, C_o2, P, y, GAS_SPECIES_INDEX)
        r6 = rate_R6(T, C_ch4, C_o2)

        o2_r5 = 0.5 * r5
        o2_r12 = 0.5 * r12
        o2_r6 = 1.5 * r6
        total = o2_r5 + o2_r12 + o2_r6
        frac_h2 = o2_r12 / total
        assert frac_h2 > 0.5

    def test_h2_trace_rate_follows_pdf_orders(self):
        """微量 H₂ 下的速率由扫描原式决定，不再要求它压过 R5。"""
        T = 1325.0
        P = 2.5e6
        C_h2, C_o2 = 0.3, 0.01
        y = np.zeros(N_GAS)
        r12 = rate_R12(T, C_h2, C_o2, P, y, GAS_SPECIES_INDEX)
        r_ref = hamel_r6_h2_rate_mol_m3_s(T, C_h2, C_o2)
        assert r12 == pytest.approx(r_ref, rel=1e-12)

    def test_co_wins_without_h2(self):
        from src.kinetics.gas_reactions import rate_R5_bubble

        T = 1325.0
        P = 2.5e6
        C_o2 = 0.01
        C_h2o = 15.0
        y = np.zeros(N_GAS)
        r5_scale = 0.75

        C_co, C_ch4 = 12.0, 0.5
        r5 = r5_scale * rate_R5_bubble(T, C_co, C_o2, C_h2o, P, y, GAS_SPECIES_INDEX)
        r6 = rate_R6(T, C_ch4, C_o2)
        o2_r5 = 0.5 * r5
        o2_r6 = 1.5 * r6
        assert o2_r5 > o2_r6
