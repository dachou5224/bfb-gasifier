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
        assert K > 1e10, "CO 氧化在 1200K 应强烈正向"

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
        from src.core.species import GAS_SPECIES_INDEX
        y = np.zeros(11)
        y[GAS_SPECIES_INDEX["CO"]] = 0.1
        y[GAS_SPECIES_INDEX["O2"]] = 0.05
        y[GAS_SPECIES_INDEX["CO2"]] = 0.1
        r = rate_R5_bubble(1200.0, 1.0, 0.5, 2.5e6, y, GAS_SPECIES_INDEX)
        assert np.isfinite(r) and r >= 0

    def test_rate_r7_is_forward_only_under_product_rich_state(self):
        from src.kinetics.gas_reactions import rate_R7

        r = rate_R7(1200.0, 0.05, 0.20, 8.0, 12.0)

        assert np.isfinite(r)
        assert r > 0.0

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
