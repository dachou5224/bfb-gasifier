from __future__ import annotations

import numpy as np
import pytest

from src.core.constants import P0, Rg
from src.thermodynamics.gibbs_hamel_reduced import ReducedHamelGibbsSolver


class _ToySpeciesDB:
    def __init__(self) -> None:
        self._mu0 = {
            "CO2": -3.0e5,
            "CO": -2.2e5,
            "CH4": -1.0e5,
            "H2": -1.0e4,
            "H2O": -2.4e5,
            "O2": 0.0,
            "N2": 0.0,
        }
        self._atoms = {
            ("CO2", "C"): 1,
            ("CO2", "O"): 2,
            ("CO", "C"): 1,
            ("CO", "O"): 1,
            ("CH4", "C"): 1,
            ("CH4", "H"): 4,
            ("H2", "H"): 2,
            ("H2O", "H"): 2,
            ("H2O", "O"): 1,
            ("O2", "O"): 2,
            ("N2", "N"): 2,
        }

    def get_mu0(self, species: str, T: float) -> float:
        return float(self._mu0[species])

    def get_atom_count(self, species: str, element: str) -> int:
        return int(self._atoms.get((species, element), 0))


def _solver() -> ReducedHamelGibbsSolver:
    return ReducedHamelGibbsSolver(_ToySpeciesDB())


def test_reduced_solver_rejects_single_element_problem() -> None:
    solver = _solver()
    with pytest.raises(ValueError, match="at least 2 tracked elements"):
        solver.solve(
            T=1200.0,
            P=101325.0,
            elements={"C": 1.0},
            candidates=["CO2", "CO"],
        )


def test_reduced_solver_rejects_empty_candidates() -> None:
    solver = _solver()
    with pytest.raises(ValueError, match="non-empty candidate species"):
        solver.solve(
            T=1200.0,
            P=101325.0,
            elements={"C": 1.0, "O": 1.0},
            candidates=[],
        )


def test_reduced_solver_toy_co_system_closes_elements() -> None:
    solver = _solver()
    out, diag = solver.solve(
        T=1200.0,
        P=101325.0,
        elements={"C": 1.0, "O": 1.0},
        candidates=["CO2", "CO", "O2"],
        return_diag=True,
    )
    assert diag["converged"] is True
    assert diag["max_abs_element_closure"] < 1e-10
    assert out["CO"] == pytest.approx(1.0, rel=0.0, abs=1e-10)
    assert out["CO2"] < 1e-10
    assert out["O2"] < 1e-10


def test_reference_element_selection_skips_largest_singleton_pool() -> None:
    solver = _solver()
    elem_names = ["C", "H", "N"]
    candidates = ["CH4", "CO", "N2"]
    a = solver._build_atom_matrix(elem_names, candidates)
    b = np.array([1.0, 2.0, 100.0], dtype=np.float64)
    ref_idx = solver._select_reference_element(a, b)
    assert elem_names[ref_idx] == "C"


def test_extract_fixed_singleton_species_extracts_n2_pool() -> None:
    solver = _solver()
    rem_elements, rem_candidates, fixed = solver._extract_fixed_singleton_species(
        elements={"C": 1.0, "O": 2.0, "N": 10.0},
        candidates=["CO2", "CO", "O2", "N2"],
    )
    assert fixed["N2"] == pytest.approx(5.0, abs=1e-12)
    assert "N2" not in rem_candidates
    assert "N" not in rem_elements


def test_reduced_jacobian_fd_matches_analytic_form() -> None:
    solver = _solver()
    elements = {"C": 1.2, "H": 2.4, "O": 1.7}
    candidates = ["CO2", "CO", "CH4", "H2", "H2O", "O2"]
    elem_names = list(elements.keys())
    a = solver._build_atom_matrix(elem_names, candidates)
    b = np.array([elements[e] for e in elem_names], dtype=np.float64)
    mu0 = np.array([solver.species_db.get_mu0(sp, 1200.0) for sp in candidates], dtype=np.float64)
    c = mu0 / (Rg * 1200.0) + np.log(101325.0 / P0)

    ref_idx = solver._select_reference_element(a, b)
    active_idx = np.array([i for i in range(len(elem_names)) if i != ref_idx], dtype=np.int64)
    lam = np.array([-1.5, 0.0, -2.0], dtype=np.float64)
    lam[ref_idx] = 0.0
    q = solver._calc_q(lam, a, c)
    g = solver._calc_g(a, q)

    jac_analytic = solver._reduced_jacobian(g=g, b=b, ref_idx=ref_idx, active_idx=active_idx)
    jac_fd = solver._reduced_jacobian_fd(
        lam=lam,
        a=a,
        c=c,
        b=b,
        ref_idx=ref_idx,
        active_idx=active_idx,
    )
    np.testing.assert_allclose(jac_fd, jac_analytic, rtol=1e-5, atol=1e-5)


def test_reduced_solver_closure_relaxed_flag_when_ls_residual_above_tol(monkeypatch: pytest.MonkeyPatch) -> None:
    solver = _solver()

    def _fake_ls_refine(*, lam, a, c, b, ref_idx, active_idx):
        return np.array(lam, dtype=np.float64), 5e-3

    monkeypatch.setattr(solver, "_least_squares_refine", _fake_ls_refine)

    out, diag = solver.solve(
        T=1200.0,
        P=101325.0,
        elements={"C": 1.0, "O": 1.0},
        candidates=["CO2", "CO", "O2"],
        lambda0=np.array([0.0, -40.0]),
        max_iter=0,
        tol=1e-10,
        return_diag=True,
    )
    assert out["CO"] > 0.0
    assert diag["used_least_squares"] is True
    assert diag["converged"] is True
    assert diag.get("converged_by") == "closure_relaxed"
    assert diag["final_residual"] == pytest.approx(5e-3, abs=1e-12)
