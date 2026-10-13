"""thermodynamics 模块单元测试。

测试 equilibrium、gibbs_minimizer、minor_species 及 gas_reactions 驱动力集成。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ===== equilibrium.py =====

class TestEquilibrium:
    """K_eq、Q_p、驱动力计算。"""

    def test_get_K_eq_r5_positive(self):
        from src.thermodynamics.equilibrium import get_K_eq
        K = get_K_eq("R5", 1200.0)
        # CO + 1/2 O2 口径的 K 约为旧 2CO+O2 口径的平方根，仍应 ≫1
        assert K > 1e6, "CO 氧化在 1200K 应强烈正向"

    def test_get_K_eq_r8_positive(self):
        from src.thermodynamics.equilibrium import get_K_eq
        K = get_K_eq("R8", 1200.0)
        assert K > 0 and np.isfinite(K), "WGSR K_eq 应为正有限值"

    def test_calc_reaction_quotient_r5(self):
        from src.thermodynamics.equilibrium import calc_reaction_quotient
        y = {"CO": 0.1, "CO2": 0.15, "O2": 0.05, "H2": 0.2, "H2O": 0.2, "CH4": 0.05, "N2": 0.25}
        Q = calc_reaction_quotient("R5", y, 2.5e6)
        assert Q > 0 and np.isfinite(Q)

    def test_calc_gibbs_driving_force_r5_clamp(self):
        from src.thermodynamics.equilibrium import calc_gibbs_driving_force
        y = {"CO": 0.1, "CO2": 0.15, "O2": 0.05, "H2": 0.2, "H2O": 0.2, "CH4": 0.05, "N2": 0.25}
        d = calc_gibbs_driving_force("R5", 1200, 2.5e6, y, clamp_irreversible=True)
        assert d >= 0.0 and d <= 1.0

    def test_calc_gibbs_driving_force_r8_reversible(self):
        from src.thermodynamics.equilibrium import calc_gibbs_driving_force
        y = {"CO": 0.2, "H2O": 0.2, "CO2": 0.1, "H2": 0.1, "N2": 0.4}
        d = calc_gibbs_driving_force("R8", 1200, 2.5e6, y, clamp_irreversible=False)
        assert np.isfinite(d)


# ===== gibbs_minimizer.py =====

class TestGibbsMinimizer:
    """Gibbs 最小化求解器。"""

    def test_solve_sulfur_simple(self):
        from src.thermodynamics.minor_species import solve_sulfur_distribution
        elements = {"S": 0.01, "H": 0.5, "O": 0.3, "C": 0.2}
        r = solve_sulfur_distribution(1200.0, 2.5e6, elements)
        assert "H2S" in r and "SO2" in r and "COS" in r
        assert all(np.isfinite(v) and v >= 0 for v in r.values())
        total_s = r["H2S"] + r["SO2"] + r["COS"]
        assert total_s >= 0, "总硫量应为非负"

    def test_solve_minor_species_supports_warm_start_and_diag(self):
        from src.thermodynamics.minor_species import solve_minor_species

        elements = {"C": 0.8, "H": 1.2, "O": 1.4, "N": 3.0}
        candidates = ["CO2", "CO", "CH4", "H2", "H2O", "O2", "N2"]
        r1, d1 = solve_minor_species(
            1200.0,
            2.5e6,
            elements,
            candidates,
            return_diag=True,
        )
        r2, d2 = solve_minor_species(
            1200.0,
            2.5e6,
            elements,
            candidates,
            lambda0=np.asarray(d1["lambda"]).tolist(),
            ln_N0=float(d1["ln_N"]),
            return_diag=True,
        )
        assert set(r1.keys()) == set(candidates)
        assert set(r2.keys()) == set(candidates)
        assert d1["n_iter"] >= 1 and d2["n_iter"] >= 1
        assert np.isfinite(d1["final_residual"])
        assert np.isfinite(d2["final_residual"])


# ===== gas_reactions 驱动力集成 =====

class TestGasReactionsDrivingForce:
    """Hamel 对齐后的 R5/R7/R8 动力学语义。"""

    def test_rate_r5_bubble_stays_forward_without_gibbs_shell(self):
        from src.kinetics.gas_reactions import rate_R5_bubble
        from src.core.species import GAS_SPECIES_INDEX, N_GAS
        y = np.zeros(N_GAS)
        y[GAS_SPECIES_INDEX["CO"]] = 0.1
        y[GAS_SPECIES_INDEX["O2"]] = 0.05
        y[GAS_SPECIES_INDEX["CO2"]] = 0.1
        r = rate_R5_bubble(1200.0, 1.0, 0.5, 0.2, 2.5e6, y, GAS_SPECIES_INDEX)
        assert np.isfinite(r) and r >= 0

    def test_rate_r5_suspension_vanishes_without_steam(self):
        from src.kinetics.gas_reactions import rate_R5_suspension
        from src.core.species import GAS_SPECIES_INDEX, N_GAS

        y = np.zeros(N_GAS)
        r = rate_R5_suspension(1200.0, 1.0, 0.5, 0.0, 2.5e6, y, GAS_SPECIES_INDEX)

        assert r == pytest.approx(0.0)

    def test_rate_r7_is_forward_only_under_product_rich_state(self):
        from src.kinetics.gas_reactions import rate_R7

        r = rate_R7(1200.0, 0.05, 0.20, 8.0, 12.0, rho_ash_kg_m3=75.0)

        assert np.isfinite(r)
        assert r > 0.0
        assert rate_R7(1200.0, 0.05, 0.20, 8.0, 12.0, rho_ash_kg_m3=0.0) == 0.0
        assert rate_R7(1200.0, 0.05, 0.20, 8.0, 12.0) == 0.0

    def test_rate_r8_zero_when_y_co_matches_hamel_equilibrium(self):
        from src.kinetics.gas_reactions import rate_R8, wgsr_equilibrium_constant, wgsr_equilibrium_y_co

        T = 1200.0
        y_h2o = 0.20
        y_co2 = 0.10
        y_h2 = 0.10
        y_co = wgsr_equilibrium_y_co(y_h2o, y_co2, y_h2, wgsr_equilibrium_constant(T))

        r = rate_R8(T, 2.5e6, y_co, y_h2o, y_co2, y_h2)

        assert abs(r) < 1e-12

    def test_rate_r8_allows_negative(self):
        from src.kinetics.gas_reactions import rate_R8
        # 高 CO2、H2，低 CO、H2O -> Q_p 大，驱动力可负
        r = rate_R8(1200.0, 2.5e6, 0.01, 0.01, 0.3, 0.3)
        assert np.isfinite(r)

    def test_rate_r8_si_zero_at_equilibrium_pressures(self):
        from src.kinetics.gas_reactions import rate_R8_si, wgsr_equilibrium_constant

        T = 1100.0
        P = 2.0e6
        keq = wgsr_equilibrium_constant(T)
        p_h2o = 0.15 * P
        p_co2 = 0.10 * P
        p_h2 = 0.12 * P
        p_co = p_co2 * p_h2 / (keq * p_h2o)
        r = rate_R8_si(T, P, p_co, p_h2o, p_co2, p_h2)
        assert abs(r) < 1e-9


class TestCharReactionsHamel51:
    """Hamel 5.1 异相炭反应的速率常数口径。"""

    def test_r1_constants_follow_hamel_table_52(self):
        from src.kinetics.char_reactions import r1_hamel_kinetic_constants
        from src.core.constants import Rg

        k0_coal, e_coal = r1_hamel_kinetic_constants("coal")
        k0_fallback, e_fallback = r1_hamel_kinetic_constants("biomass")

        assert k0_coal == pytest.approx(1.2)
        assert e_coal == pytest.approx(10_300.0 * Rg)
        assert k0_fallback == pytest.approx(2.3)
        assert e_fallback == pytest.approx(11_100.0 * Rg)

    def test_r1_n1_series_matches_hamel_eq_512_514(self):
        """§5.1.2：n=1，r = C / (1/k_ch + 1/(φ k_dg))，SPM 无灰层。"""
        from src.kinetics.arrhenius import k_hobbs
        from src.kinetics.char_reactions import (
            _k_diff_sphere,
            phi_c,
            r1_hamel_kinetic_constants,
            rate_R1,
            sherwood_number_eq517,
        )

        T = 1050.0
        C_O2 = 28.0  # [mol/m³]
        d_p = 3.0e-3  # [m]
        D_g = 6.8e-6  # [m²/s]
        u_slip = 1.0  # [m/s]
        nu_g = 1.5e-4  # [m²/s]
        Sh = sherwood_number_eq517(d_p, D_g, u_slip=u_slip, nu_g=nu_g)
        assert Sh == pytest.approx(2.0 + 0.6 * (u_slip * d_p / nu_g) ** 0.5 * (nu_g / D_g) ** (1.0 / 3.0))
        k0, e = r1_hamel_kinetic_constants("coal")
        k_ch = k_hobbs(k0, e, T)
        phi = phi_c(d_p, T)
        k_dg = _k_diff_sphere(D_g, d_p, Sh=Sh)
        expected = C_O2 / (1.0 / k_ch + 1.0 / (phi * k_dg))
        r, alpha = rate_R1(T, C_O2, d_p, D_g, None, "coal", u_slip=u_slip, nu_g=nu_g)
        assert alpha == pytest.approx(1.0 / phi)
        assert r == pytest.approx(expected, rel=1e-12)
        # 缺滑移 → Sh=2（Eq.5.17 Re→0 极限）
        r0, _ = rate_R1(T, C_O2, d_p, D_g, None, "coal")
        k_dg0 = _k_diff_sphere(D_g, d_p, Sh=2.0)
        assert r0 == pytest.approx(C_O2 / (1.0 / k_ch + 1.0 / (phi * k_dg0)), rel=1e-12)

    def test_r1_spm_ignores_d_core_ash_shrink(self):
        """主链 SPM 不因 X→d_core 引入灰层阻力。"""
        from src.kinetics.char_reactions import rate_R1

        T, C, d_p, D_g = 1100.0, 20.0, 2.5e-3, 1.0e-5
        r_spm, _ = rate_R1(T, C, d_p, D_g, d_p * 0.5, "coal", geometry="spm")
        r_none, _ = rate_R1(T, C, d_p, D_g, None, "coal", geometry="spm")
        assert r_spm == pytest.approx(r_none)
        r_ash, _ = rate_R1(T, C, d_p, D_g, d_p * 0.5, "coal", geometry="scm_ash")
        assert r_ash < r_spm

    def test_r4_boudouard_constants_follow_weeda_1995(self):
        from src.kinetics.char_reactions import rate_R4
        from src.kinetics.arrhenius import k_standard

        T = 1073.0
        p_co2 = 5.0e4
        p_co = 1.0e4
        expected = (
            k_standard(2.72, 188_000.0, T) * p_co2
            / (1.0 + k_standard(8.31e-7, -38_900.0, T) * p_co2 + k_standard(2.088e-14, -24_500.0, T) * p_co)
        )

        assert rate_R4(T, p_co2, p_co) == pytest.approx(expected)

    def test_r2_surface_rate_follows_hamel_eq_524_525(self):
        import math

        from src.core.constants import Rg
        from src.kinetics.arrhenius import k_hobbs
        from src.kinetics.char_reactions import R2_E_Rg, R2_k0, rate_R2

        T = 1100.0
        C_H2O = 62.87  # [mol/m³]
        k_ch = k_hobbs(R2_k0, R2_E_Rg * Rg, T)
        expected = k_ch * C_H2O
        assert rate_R2(T, C_H2O, 2.25e-3, 1.0e-5, 1.5e-3) == pytest.approx(expected)
        assert k_ch == pytest.approx(3.42 * T * math.exp(-15_600.0 / T))

    def test_r3_surface_rate_follows_hamel_eq_534_535(self):
        import math

        from src.core.constants import Rg
        from src.kinetics.arrhenius import k_hobbs
        from src.kinetics.char_reactions import R3_E_Rg, R3_k0, rate_R3

        T = 1100.0
        C_H2 = 12.0  # [mol/m³]
        k_ch = k_hobbs(R3_k0, R3_E_Rg * Rg, T)
        expected = k_ch * C_H2
        assert rate_R3(T, C_H2, 2.25e-3, 1.0e-5, 1.5e-3) == pytest.approx(expected)
        assert k_ch == pytest.approx(3.42e-3 * T * math.exp(-15_600.0 / T))

    def test_r4_effective_matches_hamel_eq_533_kinetic_rate(self):
        from src.kinetics.char_reactions import rate_R4, rate_R4_effective

        T = 1100.0
        p_co2 = 3.5e4
        p_co = 8.0e3
        assert rate_R4_effective(T, p_co2, p_co, 2.25e-3, 1.0e-5, 1.5e-3) == pytest.approx(
            rate_R4(T, p_co2, p_co)
        )

    def test_r4_extent_follows_hamel_eq_533_mass_basis(self):
        """Eq.5.33：R [mol C/s] = R* [1/s] · m / M_C，不得乘外表面积。"""
        from src.kinetics.char_reactions import R4_M_C, extent_R4, rate_R4

        T = 1100.0
        p_co2 = 3.5e4  # [Pa]
        p_co = 8.0e3  # [Pa]
        m_char = 2.5  # [kg]
        r_star = rate_R4(T, p_co2, p_co)
        assert extent_R4(T, p_co2, p_co, m_char) == pytest.approx(r_star * m_char / R4_M_C)


# ===== Cell 集成 =====

class TestCellGibbsIntegration:
    """Cell 中 Gibbs 相关方法。"""

    def test_calc_gibbs_correction(self):
        from src.core.cell import Cell
        from src.core.species import GAS_SPECIES_INDEX
        c = Cell()
        c.N_d[GAS_SPECIES_INDEX["CO"]] = 1.0
        c.N_d[GAS_SPECIES_INDEX["H2O"]] = 1.0
        c.N_d[GAS_SPECIES_INDEX["CO2"]] = 0.2
        c.N_d[GAS_SPECIES_INDEX["H2"]] = 0.1
        c.N_d[GAS_SPECIES_INDEX["N2"]] = 5.0
        d = c.calc_gibbs_correction("R8", "d", clamp_irreversible=False)
        assert np.isfinite(d)

    def test_get_element_release(self):
        from src.core.cell import Cell
        from src.core.species import GAS_SPECIES_INDEX
        c = Cell()
        c.N_d[GAS_SPECIES_INDEX["N2"]] = 5.0
        c.N_d[GAS_SPECIES_INDEX["NH3"]] = 0.2
        s = c.get_element_release("S")
        n = c.get_element_release("N")
        assert s >= 0 and n >= 0
        assert n == pytest.approx(0.2)

    def test_calc_minor_species_gibbs_empty(self):
        from src.core.cell import Cell
        c = Cell()
        r = c.calc_minor_species_gibbs(use_sulfur=True, use_nitrogen=True)
        assert isinstance(r, dict)

    def test_calc_minor_species_gibbs_clamps_nh3_to_minor_n_budget(self, monkeypatch):
        from src.core.cell import Cell
        from src.thermodynamics import minor_species as minor_species_mod

        c = Cell()
        monkeypatch.setattr(c, "get_element_release", lambda el: 0.2 if el.upper() == "N" else 0.0)
        monkeypatch.setattr(
            minor_species_mod,
            "solve_nitrogen_distribution",
            lambda T, P, elements: {"NH3": 1.0e12},
        )

        r = c.calc_minor_species_gibbs(use_sulfur=False, use_nitrogen=True)

        assert r["NH3"] == pytest.approx(0.2)

    def test_calc_minor_species_gibbs_reuses_exact_input_cache(self, monkeypatch):
        from src.core.cell import Cell
        from src.core.species import GAS_SPECIES_INDEX
        from src.thermodynamics import minor_species as minor_species_mod

        c = Cell()
        idx = GAS_SPECIES_INDEX
        c.N_d[idx["CO"]] = 1.0
        c.N_d[idx["H2O"]] = 1.0
        c.N_d[idx["N2"]] = 5.0

        calls = {"n": 0}

        def _fake_n_solver(T, P, elements, **kwargs):
            calls["n"] += 1
            if kwargs.get("return_diag"):
                return {"NH3": 0.05}, {"lambda": np.zeros(4), "ln_N": 0.0}
            return {"NH3": 0.05}

        monkeypatch.setattr(c, "get_element_release", lambda el: 0.2 if el.upper() == "N" else 0.0)
        monkeypatch.setattr(minor_species_mod, "solve_nitrogen_distribution", _fake_n_solver)

        r1 = c.calc_minor_species_gibbs(use_sulfur=False, use_nitrogen=True)
        r2 = c.calc_minor_species_gibbs(use_sulfur=False, use_nitrogen=True)

        assert calls["n"] == 1
        assert r1 == r2

    def test_calc_minor_species_gibbs_cache_invalidates_on_state_change(self, monkeypatch):
        from src.core.cell import Cell
        from src.core.species import GAS_SPECIES_INDEX
        from src.thermodynamics import minor_species as minor_species_mod

        c = Cell()
        idx = GAS_SPECIES_INDEX
        c.N_d[idx["CO"]] = 1.0
        c.N_d[idx["H2O"]] = 1.0
        c.N_d[idx["N2"]] = 5.0

        calls = {"n": 0}

        def _fake_n_solver(T, P, elements, **kwargs):
            calls["n"] += 1
            if kwargs.get("return_diag"):
                return {"NH3": 0.05}, {"lambda": np.zeros(4), "ln_N": 0.0}
            return {"NH3": 0.05}

        monkeypatch.setattr(c, "get_element_release", lambda el: 0.2 if el.upper() == "N" else 0.0)
        monkeypatch.setattr(minor_species_mod, "solve_nitrogen_distribution", _fake_n_solver)

        c.calc_minor_species_gibbs(use_sulfur=False, use_nitrogen=True)
        c.T += 1.0
        c.calc_minor_species_gibbs(use_sulfur=False, use_nitrogen=True)

        assert calls["n"] == 2
