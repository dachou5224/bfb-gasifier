"""Hamel Table 5.3 / Eq. 5.36–5.37 与 ``gas_reactions.rate_R5_*`` 对照。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.constants import P0, Rg
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.kinetics import gas_reactions as gr
from src.kinetics.arrhenius import k_standard
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension
from src.kinetics.hamel_r5_reference import (
    HAMEL_R5_BUBBLE_BETA,
    HAMEL_R5_BUBBLE_E_OVER_RG_K,
    HAMEL_R5_BUBBLE_GAMMA,
    HAMEL_R5_BUBBLE_K0,
    hamel_r5_bubble_k,
    hamel_r5_bubble_rate_mol_m3_s,
    hamel_r5_suspension_rate_mol_m3_s,
)


class TestR5HamelTable53Constants:
    def test_bubble_k0_and_activation_match_table(self):
        assert gr.R5_bubble_k0 == pytest.approx(HAMEL_R5_BUBBLE_K0)
        assert gr.R5_bubble_E_Rg == pytest.approx(HAMEL_R5_BUBBLE_E_OVER_RG_K)
        assert gr.R5_bubble_beta == pytest.approx(HAMEL_R5_BUBBLE_BETA)
        assert gr.R5_bubble_gamma == pytest.approx(HAMEL_R5_BUBBLE_GAMMA)

    def test_eq_5_37_bubble_k_matches_k_standard(self):
        T_K = 1030.0
        e_j = HAMEL_R5_BUBBLE_E_OVER_RG_K * Rg
        k_ref = hamel_r5_bubble_k(T_K)
        k_impl = k_standard(gr.R5_bubble_k0, e_j, T_K)
        assert k_impl == pytest.approx(k_ref, rel=0.0, abs=0.0)


class TestR5BubbleImplementation:
    def test_bubble_rate_matches_eq_5_36(self):
        T = 1111.0
        P = 2.5e6
        C_co, C_o2, C_h2o = 12.0, 0.8, 15.0
        y = np.zeros(N_GAS)
        r_impl = rate_R5_bubble(T, C_co, C_o2, C_h2o, P, y, GAS_SPECIES_INDEX)
        r_ref = hamel_r5_bubble_rate_mol_m3_s(T, C_co, C_o2, C_h2o)
        assert r_impl == pytest.approx(r_ref, rel=1e-12)

    def test_bubble_uses_table53_beta_gamma_not_o2_half_order(self):
        """Eq. 5.36: C_O2^0.3 · C_H2O^0.5，非旧实现 C_O2^0.5。"""
        T = 1030.0
        C_co, C_o2, C_h2o = 10.0, 0.5, 15.0
        k = hamel_r5_bubble_k(T)
        r_eq536 = k * C_co * (C_o2**0.3) * (C_h2o**0.5)
        r_o2_05 = k * C_co * (C_o2**0.5) * (C_h2o**0.5)
        r_impl = rate_R5_bubble(T, C_co, C_o2, C_h2o, 2.5e6, np.zeros(N_GAS), GAS_SPECIES_INDEX)
        assert r_impl == pytest.approx(r_eq536, rel=1e-12)
        assert r_impl != pytest.approx(r_o2_05, rel=0.05)

    def test_bubble_vanishes_without_steam(self):
        y = np.zeros(N_GAS)
        r = rate_R5_bubble(1030.0, 1.0, 0.5, 0.0, 2.5e6, y, GAS_SPECIES_INDEX)
        assert r == pytest.approx(0.0)


class TestR5SuspensionImplementation:
    def test_suspension_matches_hamel_table53(self):
        T = 1030.0
        C_co, C_o2, C_h2o = 8.0, 0.2, 15.0
        y = np.zeros(N_GAS)
        r_impl = rate_R5_suspension(T, C_co, C_o2, C_h2o, 2.5e6, y, GAS_SPECIES_INDEX)
        r_hamel = hamel_r5_suspension_rate_mol_m3_s(T, C_co, C_o2, C_h2o, temperature_exponent=-0.5)
        assert r_impl == pytest.approx(r_hamel, rel=1e-12)

        k_code = gr.R5_suspension_k_T05 * (T ** -0.5)
        p_co = C_co * Rg * T / P0
        p_o2 = C_o2 * Rg * T / P0
        p_h2o = C_h2o * Rg * T / P0
        r_code_form = k_code * p_co * (p_o2**0.5) * (p_h2o**0.5)
        assert r_impl == pytest.approx(r_code_form, rel=1e-12)

    def test_suspension_vanishes_without_steam(self):
        y = np.zeros(N_GAS)
        r = rate_R5_suspension(1030.0, 1.0, 0.5, 0.0, 2.5e6, y, GAS_SPECIES_INDEX)
        assert r == pytest.approx(0.0)
