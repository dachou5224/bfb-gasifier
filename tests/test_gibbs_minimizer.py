from __future__ import annotations

import numpy as np

from src.thermodynamics.gibbs_minimizer import GibbsMinimizer, hamel_a40_initial_guess


class _DummySpeciesDB:
    def __init__(self) -> None:
        self._mu0 = {
            "H2": 0.0,
            "H2O": -220_000.0,
            "O2": 0.0,
        }
        self._atoms = {
            "H2": {"H": 2, "O": 0},
            "H2O": {"H": 2, "O": 1},
            "O2": {"H": 0, "O": 2},
        }

    def get_mu0(self, species: str, T: float) -> float:
        _ = T
        return float(self._mu0[species])

    def get_atom_count(self, species: str, element: str) -> int:
        return int(self._atoms.get(species, {}).get(element, 0))


def test_gibbs_minimizer_defaults_use_hamel_a40_and_holub_tol() -> None:
    from src.core.hamel_parity_defaults import (
        HAMEL_GIBBS_INITIAL_GUESS_MODE,
        HAMEL_GIBBS_TOL_MODE,
    )

    g = GibbsMinimizer(_DummySpeciesDB())
    _, diag = g.solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        return_diag=True,
    )
    assert diag["initial_guess_mode"] == HAMEL_GIBBS_INITIAL_GUESS_MODE
    assert diag["tol_mode"] == HAMEL_GIBBS_TOL_MODE
    assert diag["converged"] is True


def test_gibbs_minimizer_returns_nonnegative_flows():
    g = GibbsMinimizer(_DummySpeciesDB())
    out = g.solve(
        T=1200.0,
        P=2.5e6,
        elements={"H": 2.0, "O": 1.0},
        candidates=["H2", "H2O", "O2"],
    )
    vals = np.array([out["H2"], out["H2O"], out["O2"]], dtype=float)
    assert np.all(np.isfinite(vals))
    assert np.all(vals >= 0.0)
    assert float(np.sum(vals)) > 0.0


def test_gibbs_minimizer_holub_vonka_a34_mode_converges() -> None:
    g = GibbsMinimizer(_DummySpeciesDB())
    out, diag = g.solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        tol=1e-8,
        tol_mode="holub_vonka_a34",
        return_diag=True,
    )
    assert isinstance(out, dict)
    assert diag["tol_mode"] == "holub_vonka_a34"
    assert diag["converged"] is True
    assert diag.get("converged_by") == "holub_vonka_a34"
    assert float(diag["holub_vonka_a34"]) < 1e-8


def test_gibbs_minimizer_rejects_unknown_tol_mode() -> None:
    g = GibbsMinimizer(_DummySpeciesDB())
    try:
        g.solve(
            T=1100.0,
            P=2.0e6,
            elements={"H": 2.0, "O": 1.0},
            candidates=["H2", "H2O", "O2"],
            tol_mode="unknown_mode",
        )
    except ValueError as exc:
        assert "tol_mode" in str(exc)
    else:
        raise AssertionError("expected ValueError for unknown tol_mode")


def test_gibbs_minimizer_hamel_a40_initial_guess_mode_converges() -> None:
    g = GibbsMinimizer(_DummySpeciesDB())
    out, diag = g.solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        initial_guess_mode="hamel_a40",
        return_diag=True,
    )
    assert diag["initial_guess_mode"] == "hamel_a40"
    assert diag["converged"] is True
    vals = np.array([out["H2"], out["H2O"], out["O2"]], dtype=float)
    assert np.all(vals >= 0.0)
    assert float(np.sum(vals)) > 0.0


def test_hamel_a40_initial_guess_produces_finite_lambda() -> None:
    a = np.array([[2.0, 2.0, 0.0], [0.0, 1.0, 2.0]], dtype=np.float64)
    b = np.array([3.0, 1.5], dtype=np.float64)
    c = np.array([0.0, -1.0, 0.5], dtype=np.float64)
    lam, ln_n = hamel_a40_initial_guess(a, b, c)
    assert lam.shape == (2,)
    assert np.all(np.isfinite(lam))
    assert np.isfinite(ln_n)


def test_gibbs_minimizer_diag_mode_exposes_convergence_metadata():
    g = GibbsMinimizer(_DummySpeciesDB())
    out, diag = g.solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        return_diag=True,
    )
    assert isinstance(out, dict)
    assert "converged" in diag
    assert "n_iter" in diag
    assert "final_residual" in diag
    assert np.isfinite(float(diag["final_residual"]))
