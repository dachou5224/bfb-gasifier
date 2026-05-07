from __future__ import annotations

import numpy as np

from src.thermodynamics.gibbs_minimizer import GibbsMinimizer


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
