"""预投影 trial / LSQ / restore 共享引擎。

Ref: docs/opt_implementation_plan.md OPT-002 Phase 1;
     Hamel (1999) NR Startwertwahl bridge patterns in connectivity.py.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np

from src.core.cell import Cell
from src.core.cell_residual_pipeline import execute_cell_residual_pipeline


def preproj_cell_residual_vector(
    cell: Cell,
    *groups: str,
    reaction_rate_multiplier: float | None = None,
) -> np.ndarray:
    """Partial/full cell residual view for preprojection LSQ (no copy)."""
    rm = (
        float(reaction_rate_multiplier)
        if reaction_rate_multiplier is not None
        else float(cell.nr_reaction_rate_multiplier)
    )
    return execute_cell_residual_pipeline(
        cell,
        rate_multiplier=rm,
        groups=frozenset(groups),
        copy=False,
    )


def preproj_cell_energy_residual_abs(cell: Cell) -> float:
    out = execute_cell_residual_pipeline(
        cell,
        rate_multiplier=cell.nr_reaction_rate_multiplier,
        groups=frozenset({"energy"}),
        copy=False,
    )
    return abs(float(out[0]))


def preproj_refresh_trial_bc(
    cell_index: int,
    *,
    apply_bc_fn: Callable[[], None] | None,
    apply_local_bc_fn: Callable[[int], None] | None = None,
    neighbor_indices: tuple[int, ...] = (),
) -> None:
    if apply_local_bc_fn is not None:
        apply_local_bc_fn(int(cell_index))
        for idx in neighbor_indices:
            apply_local_bc_fn(int(idx))
    elif apply_bc_fn is not None:
        apply_bc_fn()


def preproj_finalize_bc(apply_bc_fn: Callable[[], None] | None) -> None:
    if apply_bc_fn is not None:
        apply_bc_fn()


@dataclass(frozen=True)
class PreprojectionBcContext:
    """Trial 内 local BC refresh；accept/reject 后 full BC finalize。"""

    cell_index: int
    apply_bc_fn: Callable[[], None] | None = None
    apply_local_bc_fn: Callable[[int], None] | None = None
    neighbor_indices: tuple[int, ...] = ()

    def needs_bc(self) -> bool:
        return self.apply_bc_fn is not None or self.apply_local_bc_fn is not None

    def finalize(self) -> None:
        preproj_finalize_bc(self.apply_bc_fn)

    def refresh_trial(self) -> None:
        preproj_refresh_trial_bc(
            self.cell_index,
            apply_bc_fn=self.apply_bc_fn,
            apply_local_bc_fn=self.apply_local_bc_fn,
            neighbor_indices=self.neighbor_indices,
        )


@dataclass(frozen=True)
class LsqPreprojectionSpec:
    x0: np.ndarray  # [-]
    lower: np.ndarray  # [-]
    upper: np.ndarray  # [-]
    max_nfev: int = 32  # [-]
    xtol: float = 1.0e-10  # [-]
    ftol: float = 1.0e-10  # [-]
    gtol: float = 1.0e-10  # [-]


def run_improvement_gated_lsq(
    *,
    spec: LsqPreprojectionSpec,
    set_values: Callable[[np.ndarray], None],
    objective: Callable[[np.ndarray], np.ndarray],
    initial_metric: float,
    final_metric: Callable[[], float],
    improvement_factor: float,
    restore: Callable[[], None],
    sync_full_residuals: Callable[[], None],
    bc: PreprojectionBcContext | None = None,
    metric_accept: Callable[[float, float], bool] | None = None,
    finalize_on_accept: bool = True,
    finalize_on_reject: bool = True,
    sync_on_accept: bool = True,
    on_accept: Callable[[], None] | None = None,
) -> bool:
    """Bounded LSQ preprojection with improvement-factor accept gate.

    On failure (LSQ exception or gate reject): restore snapshot, optional full BC,
    sync full cell residuals, return False.
    """
    if bc is not None and bc.needs_bc():
        bc.finalize()

    if not np.isfinite(initial_metric):
        return False

    def _default_accept(final: float, initial: float) -> bool:
        return np.isfinite(final) and final <= improvement_factor * initial

    accept_fn = metric_accept or _default_accept

    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            objective,
            x0=np.clip(spec.x0, spec.lower, spec.upper),
            bounds=(spec.lower, spec.upper),
            max_nfev=int(spec.max_nfev),
            xtol=float(spec.xtol),
            ftol=float(spec.ftol),
            gtol=float(spec.gtol),
        )
        set_values(np.asarray(sol.x, dtype=np.float64))
        if bc is not None:
            bc.refresh_trial()
        final = float(final_metric())
    except Exception:
        restore()
        if bc is not None and finalize_on_reject:
            bc.finalize()
        sync_full_residuals()
        return False

    if accept_fn(final, initial_metric):
        if bc is not None and finalize_on_accept:
            bc.finalize()
        if on_accept is not None:
            on_accept()
        if sync_on_accept:
            sync_full_residuals()
        return True

    restore()
    if bc is not None and finalize_on_reject:
        bc.finalize()
    sync_full_residuals()
    return False


def scaled_rms(vector: np.ndarray, ref_scale: float) -> float:
    """Scaled RMS for preprojection merit gates."""
    ref = max(float(ref_scale), 1.0e-12)
    v = np.asarray(vector, dtype=np.float64)
    if v.size == 0:
        return 0.0
    return float(np.sqrt(np.mean((v / ref) ** 2)))


def run_lsq_preprojection(
    *,
    spec: LsqPreprojectionSpec,
    objective: Callable[[np.ndarray], np.ndarray],
    apply_solution: Callable[[np.ndarray], None],
    restore: Callable[[], None],
    accept_predicate: Callable[[], bool],
    sync_full_residuals: Callable[[], None],
    bc: PreprojectionBcContext | None = None,
    finalize_before: bool = False,
    finalize_on_reject: bool = True,
    sync_on_reject: bool = True,
    on_accept: Callable[[], None] | None = None,
) -> bool:
    """Bounded LSQ with arbitrary post-process and accept predicate.

    Ref: ``preproject_bed_gas_temperature_to_energy_closure`` gas/T bridge.
    """
    if finalize_before and bc is not None and bc.needs_bc():
        bc.finalize()

    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            objective,
            x0=np.clip(spec.x0, spec.lower, spec.upper),
            bounds=(spec.lower, spec.upper),
            max_nfev=int(spec.max_nfev),
            xtol=float(spec.xtol),
            ftol=float(spec.ftol),
            gtol=float(spec.gtol),
        )
        apply_solution(np.asarray(sol.x, dtype=np.float64))
    except Exception:
        restore()
        if bc is not None and finalize_on_reject:
            bc.finalize()
        if sync_on_reject:
            sync_full_residuals()
        return False

    if accept_predicate():
        if on_accept is not None:
            on_accept()
        return True

    restore()
    if bc is not None and finalize_on_reject:
        bc.finalize()
    if sync_on_reject:
        sync_full_residuals()
    return False


def run_relax_schedule(
    candidates: Iterable[float],
    *,
    restore: Callable[[], None],
    try_trial: Callable[[float], bool],
    prepare_trial: Callable[[], None] | None = None,
    finalize_failure: Callable[[], None] | None = None,
) -> tuple[bool, float]:
    """Restore → prepare → trial loop over relax candidates (halving schedule).

    Ref: ``run_bottom_joint_start_value_preprojection``; Hamel bed0 joint Startwertwahl.
    """
    for value in candidates:
        restore()
        if prepare_trial is not None:
            prepare_trial()
        if try_trial(float(value)):
            return True, float(value)
    restore()
    if finalize_failure is not None:
        finalize_failure()
    return False, 0.0
