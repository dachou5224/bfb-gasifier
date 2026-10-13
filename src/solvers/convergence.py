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
    abs_tol_T: float = 1.0  # [K]
    enable_gas_phase_split_gate: bool = False
    max_gas_phase_split_rms: float = 0.05
    enable_bed0_primary_gate: bool = False
    max_bed0_primary_rel_drift: float = 0.10

    def is_aligned(
        self,
        dT_outer_max: float,
        T_ref_scale: float,
        *,
        gas_phase_split_rms: float | None = None,
        bed0_primary_rel_drift: float | None = None,
    ) -> bool:
        """Check physical alignment between Vorabrechnung and cell solve.

        Temperature: max |ΔT| after one inner pass vs pre-refresh state.
        Gas (optional): bed major phase-split RMS and bed0 primary vs zu drift.
        """
        t_scale = max(float(T_ref_scale), 1.0)
        tol_abs = max(float(self.abs_tol_T), float(self.rel_tol) * t_scale)
        if float(dT_outer_max) > tol_abs:
            return False
        if self.enable_gas_phase_split_gate and gas_phase_split_rms is not None:
            if float(gas_phase_split_rms) > float(self.max_gas_phase_split_rms):
                return False
        if self.enable_bed0_primary_gate and bed0_primary_rel_drift is not None:
            if float(bed0_primary_rel_drift) > float(self.max_bed0_primary_rel_drift):
                return False
        return True


@dataclass(frozen=True)
class OuterMeritEarlyStop:
    """Composite merit early-stop for outer Abgleich (thesis Check1 merit).

    Merit = max(rms_scaled, gas_phase_split_rms, energy_rms) when the
    corresponding flags are enabled.  Single-factor RMS stall is avoided
    because energy and phase-split can move independently of bulk RMS.
    """

    enabled: bool = True
    min_outer_iters: int = 8
    regression_rel_tol: float = 0.25
    regression_patience: int = 1
    plateau_rel_tol: float = 0.002
    use_gas_phase_split_merit: bool = True
    use_energy_merit: bool = True

