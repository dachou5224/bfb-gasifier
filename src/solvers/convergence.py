"""Convergence objects for Hamel-style dual-loop solve.

Check1: inner NR residual convergence (equation residual space).
Check2: outer Abgleich physical alignment (state-space consistency).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class InnerConvergence:
    """Inner NR convergence criterion (Check1)."""

    atol: float = 1e-8
    rtol: float = 1e-6

    def is_converged(self, residual_norm: float, x_norm: float) -> bool:
        threshold = float(self.atol) + float(self.rtol) * max(float(x_norm), 0.0)
        return float(residual_norm) <= threshold


@dataclass(frozen=True)
class OuterConvergence:
    """Outer physical alignment criterion (Check2)."""

    rel_tol: float = 1e-4
    abs_tol_T: float = 0.5  # [K]

    def is_aligned(self, dT_outer_max: float, T_ref_scale: float) -> bool:
        """Check physical alignment between Vorabrechnung and cell solve.

        Parameters
        ----------
        dT_outer_max:
            Max absolute axial temperature mismatch [K] after one outer update.
        T_ref_scale:
            Reference temperature scale [K] for relative tolerance.
        """
        t_scale = max(float(T_ref_scale), 1.0)
        tol_abs = max(float(self.abs_tol_T), float(self.rel_tol) * t_scale)
        return float(dT_outer_max) <= tol_abs

