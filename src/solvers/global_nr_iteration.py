"""Global NR main loop extracted from global_nr_solver (OPT-005).

Source: Hamel (1999) §2.3 Eq. 2.9 damped Newton iteration body.
"""

from __future__ import annotations

from collections.abc import Callable
from time import perf_counter
from typing import Any, List, Tuple

import numpy as np
import scipy.sparse as sp

from src.core.cell import Cell
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers import global_nr_solver as _gnr
from src.solvers.convergence import InnerConvergence
from src.solvers.nr_indexing import (
    cell_offsets,
    gas_mode_for_nr,
    has_temperature_var,
    main_nr_gas_idx_for_cell,
    n_var,
    pack_reactor,
    temperature_offset,
    unpack_reactor,
)

# Excluded-gas snapshot: (N_d[excl], N_b[excl], excl_indices)
_ExcludedGasSnap = Tuple[np.ndarray, np.ndarray, np.ndarray]


def _excluded_gas_indices_for_cell(cell: Cell) -> np.ndarray:
    """Gas species indices not packed into the NR state vector (e.g. H2S)."""
    packed = {int(i) for i in main_nr_gas_idx_for_cell(cell)}
    return np.asarray([i for i in range(N_GAS) if i not in packed], dtype=np.int64)


def _snapshot_excluded_gas_holdups(cells: List[Cell]) -> list[_ExcludedGasSnap]:
    """Snapshot NR-excluded gas holdups so rejected fastox trials can roll back."""
    snaps: list[_ExcludedGasSnap] = []
    for cell in cells:
        excl = _excluded_gas_indices_for_cell(cell)
        snaps.append(
            (
                np.asarray(cell.N_d[excl], dtype=np.float64).copy(),
                np.asarray(cell.N_b[excl], dtype=np.float64).copy(),
                excl,
            )
        )
    return snaps


def _restore_excluded_gas_holdups(
    cells: List[Cell],
    snaps: list[_ExcludedGasSnap],
) -> None:
    for cell, (nd, nb, excl) in zip(cells, snaps):
        if excl.size == 0:
            continue
        cell.N_d[excl] = nd
        cell.N_b[excl] = nb


def _bed_species_holdup_mol(cells: List[Cell], species: str) -> float:
    """Bed total holdup [mol] for one gas species."""
    idx = GAS_SPECIES_INDEX
    j = idx[species]
    total = 0.0
    for cell in cells:
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        total += float(max(cell.N_d[j], 0.0) + max(cell.N_b[j], 0.0))
    return float(total)


def _bed_syngas_holdup_mol(cells: List[Cell]) -> float:
    """Bed total CO+H₂+CH₄ holdup [mol] for LS inventory-collapse guard."""
    return float(
        _bed_species_holdup_mol(cells, "CO")
        + _bed_species_holdup_mol(cells, "H2")
        + _bed_species_holdup_mol(cells, "CH4")
    )


_SYNGAS_SPECIES = ("CO", "H2", "CH4")


def _syngas_packed_indices(cells: List[Cell]) -> list[int]:
    """Global packed-state indices of CO/H₂/CH₄ holdup unknowns."""
    syngas = {GAS_SPECIES_INDEX[sp] for sp in _SYNGAS_SPECIES}
    indices: list[int] = []
    offset = 0
    for cell in cells:
        nv = n_var(cell)
        gas_idx = main_nr_gas_idx_for_cell(cell)
        mode = gas_mode_for_nr(cell)
        if mode == "two_phase":
            packed_sp = list(gas_idx) + list(gas_idx)
        elif mode == "single":
            packed_sp = list(gas_idx)
        else:
            packed_sp = []
        for local, sp in enumerate(packed_sp):
            if int(sp) in syngas:
                indices.append(offset + local)
        offset += nv
    return indices


def _floor_syngas_holdups_in_state(
    x_trial: np.ndarray,
    x_base: np.ndarray,
    cells: List[Cell],
    ratio: float,
) -> int:
    """Floor packed CO/H₂/CH₄ so a trial cannot erase syngas before fastox/merit.

    Returns the number of components that were raised to the floor.
    """
    floor_ratio = float(np.clip(ratio, 0.0, 1.0))
    n_floored = 0
    for idx in _syngas_packed_indices(cells):
        if idx >= x_trial.size or idx >= x_base.size:
            continue
        floor = floor_ratio * max(float(x_base[idx]), 0.0)
        if float(x_trial[idx]) < floor:
            x_trial[idx] = floor
            n_floored += 1
    return int(n_floored)


def _recompute_residual_after_hydro_refresh(
    *,
    x: np.ndarray,
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    scale: np.ndarray,
    line_search_merit_enabled: bool,
    line_search_energy_merit_enabled: bool,
    timing: dict[str, float],
    counts: dict[str, int],
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Recompute F after accepted-step hydro refresh so next RHS/J share one mapping.

    Accept-step hydro refresh changes K_bd/ε etc. Keeping the pre-refresh trial
    residual as the next Newton RHS would pair a stale F with a refreshed J.
    """
    t0 = perf_counter()
    F = _gnr.global_residual(x, cells, apply_bc_fn)
    timing["line_search_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] = int(counts.get("global_residual_calls", 0)) + 1
    counts["hydrodynamics_refresh_residual_recompute"] = int(
        counts.get("hydrodynamics_refresh_residual_recompute", 0)
    ) + 1
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    norm_F = _gnr._rms_norm(F_hat)
    merit_F = _gnr._line_search_merit(
        F_hat,
        cells,
        use_gas_phase_split=line_search_merit_enabled,
        use_energy=line_search_energy_merit_enabled,
    )
    return F, F_hat, float(norm_F), float(merit_F)


def execute_global_nr_solve(
    cells: List[Cell],
    apply_bc_fn: Callable[[], None],
    apply_local_bc_fn: Callable[[int], None] | None = None,
    affected_residual_cells_fn: Callable[[int], tuple[int, ...]] | None = None,
    ref_gas_mol_s: float = 80.0,
    ref_solid_kg_s: float = 1.0,
    ref_energy_W: float = 2e7,
    max_iter: int = 25,
    tol_rms: float = 0.01,
    inner_convergence: InnerConvergence | None = None,
    lambda_init: float = 0.5,
    n_damp_halvings: int = 12,
    lambda_min: float = 1.0 / 1024.0,
    jacobian_strategy: str = "block_tridiag_structured",
    linear_solver_backend: str | None = None,
    jacobian_lag_steps: int = 1,
    allow_gd_fallback: bool = False,
    prefer_full_step: bool = False,
    step_model: str = "newton",
    lm_mu0: float = 1.0e-4,
    lm_mu_growth: float = 10.0,
    ptc_alpha0: float = 1.0,
    ptc_alpha_growth: float = 2.0,
    equil_iters: int = 0,
    equil_scale_clip: float = 1.0e3,
    line_search_max_trials: int = 0,
    nonmonotone_enabled: bool = False,
    nonmonotone_window: int = 5,
    nonmonotone_relax: float = 1.0,
    enforce_gas_phase_split_convergence: bool = False,
    tol_gas_phase_split_max: float | None = None,
    line_search_use_gas_phase_split_merit: bool = False,
    line_search_use_energy_merit: bool = False,
    refresh_hydrodynamics_on_accepted_step: bool = False,
    line_search_fast_oxidation_projection: bool = False,
    line_search_per_phase_fastox: bool = False,
    line_search_per_phase_fastox_bed_limit: int = 0,
    clip_same_phase_oxidizer_fuel_step: bool = False,
    allow_oxidizer_into_fuel_step: bool = False,
    line_search_syngas_collapse_guard: bool = True,
    line_search_syngas_collapse_ratio: float = 0.05,
    line_search_syngas_floor_ratio: float = 0.0,
    line_search_refresh_aware_accept: bool = False,
    line_search_split_constrained_step: bool = False,
    clip_solid_ftb_relative_floor: bool = False,
    inner_t_step_default_cap_K: float = 600.0,
    inner_split_dominant_t_step_cap_K: float | None = 150.0,
    inner_split_dominant_merit_ratio: float = 0.85,
    inner_t_cap_hold_split_dominant: bool = False,
    inner_t_cap_split_release_ratio: float = 0.70,
    inner_abs_dT_budget_K: float | None = None,
    connectivity_graph: Any | None = None,
    record_inner_group_history: bool = False,
    verbose: bool = False,
) -> dict:
    """Damped Newton–Raphson main loop (Hamel Eq. 2.9)."""
    solve_started = perf_counter()
    inner_group_history: list[dict[str, float | int | bool]] = []
    t_cap_hold_split = bool(inner_t_cap_hold_split_dominant)
    T_ref_abs = np.array(
        [float(c.T) if has_temperature_var(c) else float("nan") for c in cells],
        dtype=np.float64,
    )
    x = pack_reactor(cells)
    timing = {
        "initial_residual_s": 0.0,
        "jacobian_build_s": 0.0,
        "linear_solve_s": 0.0,
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "line_search_residual_s": 0.0,
        "final_residual_s": 0.0,
        "total_s": 0.0,
    }
    counts = {
        "global_residual_calls": 0,
        "jacobian_cell_residual_calls": 0,
        "jacobian_nnz_last": 0,
        "jacobian_zero_cols_last": 0,
        "jacobian_zero_rows_last": 0,
        "jacobian_rebuilds": 0,
        "jacobian_reuses": 0,
        "line_search_evaluations": 0,
        "line_search_backtracks": 0,
        "line_search_failures": 0,
        "line_search_retries": 0,
        "linear_regularized_solves": 0,
        "linear_lstsq_fallbacks": 0,
        "structured_fallbacks": 0,
        "nr_schur_size_last": 0,
    }
    accepted_lambda_history: list[float | None] = []
    line_search_trial_counts: list[int] = []
    clip_history: list[dict] = []
    group_norms = {"gas": float("inf"), "solid": float("inf"), "energy": float("inf")}
    group_max_abs = {"gas": float("inf"), "solid": float("inf"), "energy": float("inf")}

    t0 = perf_counter()
    F = _gnr.global_residual(x, cells, apply_bc_fn)
    timing["initial_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] += 1
    x = pack_reactor(cells)
    scale = _gnr.build_equation_scales(cells, ref_gas_mol_s, ref_solid_kg_s, ref_energy_W)

    n_total = len(F)

    F_hat = F / scale
    norm_F = _gnr._rms_norm(F_hat)
    line_search_merit_enabled = bool(line_search_use_gas_phase_split_merit)
    line_search_energy_merit_enabled = bool(line_search_use_energy_merit)
    refresh_hydro_on_accept = bool(refresh_hydrodynamics_on_accepted_step)
    ls_per_phase_fastox = bool(line_search_per_phase_fastox)
    ls_per_phase_bed_limit = int(line_search_per_phase_fastox_bed_limit)
    ls_fast_ox_proj = bool(line_search_fast_oxidation_projection) or ls_per_phase_fastox
    ls_syngas_guard = bool(line_search_syngas_collapse_guard)
    ls_syngas_ratio = float(np.clip(line_search_syngas_collapse_ratio, 0.0, 1.0))
    ls_syngas_floor_ratio = float(np.clip(line_search_syngas_floor_ratio, 0.0, 1.0))
    ls_syngas_floor = bool(ls_syngas_guard and ls_syngas_floor_ratio > 0.0)
    ls_refresh_aware = bool(line_search_refresh_aware_accept)
    ls_split_constrained = bool(line_search_split_constrained_step)
    clip_same_phase_ox_fuel = bool(clip_same_phase_oxidizer_fuel_step)
    allow_ox_into_fuel = bool(allow_oxidizer_into_fuel_step)
    merit_F = _gnr._line_search_merit(
        F_hat,
        cells,
        use_gas_phase_split=line_search_merit_enabled,
        use_energy=line_search_energy_merit_enabled,
    )
    history: List[float] = [merit_F]
    converged = False
    lambda_seed = 1.0 if bool(prefer_full_step) else float(lambda_init)
    consecutive_line_search_failures = 0
    iter_attempts = 0
    lag_steps = max(int(jacobian_lag_steps), 1)
    cached_J: sp.csr_matrix | None = None
    cached_jac_meta: dict | None = None
    cached_jac_iter = -1
    linear_solver_backend_last = ""
    if linear_solver_backend is None:
        from src.solvers.global_nr_solver import structured_jacobian_uses_direct_linear_solver

        linear_solver_backend = (
            "structured_direct"
            if structured_jacobian_uses_direct_linear_solver(jacobian_strategy)
            else "sparse_direct_fallback"
        )
    step_model_resolved = str(step_model).strip().lower()
    if step_model_resolved not in {"newton", "lm", "ptc", "equil_newton"}:
        raise ValueError(f"Unsupported step_model={step_model!r}; expected 'newton', 'lm', 'ptc', or 'equil_newton'")
    lm_mu0_resolved = float(max(lm_mu0, 1.0e-12))
    lm_mu_growth_resolved = float(max(lm_mu_growth, 1.0))
    ptc_alpha0_resolved = float(max(ptc_alpha0, 1.0e-12))
    ptc_alpha_growth_resolved = float(max(ptc_alpha_growth, 1.0))
    equil_iters_resolved = int(max(equil_iters, 0))
    equil_scale_clip_resolved = float(max(equil_scale_clip, 1.0))
    line_search_max_trials_resolved = int(max(line_search_max_trials, 0))
    nonmonotone_enabled_resolved = bool(nonmonotone_enabled)
    nonmonotone_window_resolved = int(max(nonmonotone_window, 1))
    nonmonotone_relax_resolved = float(max(nonmonotone_relax, 1.0))

    # Pre-compute T DOF indices (last DOF of each cell block).
    # Used by _solve_linear_step to pass Tikhonov regularisation info.
    _t_dof_idx: np.ndarray | None = None
    _t_dof_list: list[int] = []
    _off = 0
    for _c in cells:
        _nv = n_var(_c)
        if has_temperature_var(_c):
            _t_dof_list.append(_off + temperature_offset(_c))
        _off += _nv
    if _t_dof_list:
        _t_dof_idx = np.array(_t_dof_list, dtype=np.intp)

    for it in range(max_iter):
        # 总是按真实 cell-by-cell 布局计算三类残差 RMS，用于自适应 T 步长判断。
        group_metrics = _gnr._residual_group_metrics(F_hat, cells)
        group_norms = {
            "gas": group_metrics["gas_rms"],
            "solid": group_metrics["solid_rms"],
            "energy": group_metrics["energy_rms"],
        }
        group_max_abs = {
            "gas": group_metrics["gas_max_abs"],
            "solid": group_metrics["solid_max_abs"],
            "energy": group_metrics["energy_max_abs"],
        }
        gas_norm = group_norms["gas"]
        sol_norm = group_norms["solid"]
        ene_norm = group_norms["energy"]
        gas_phase_iter = _gnr._residual_gas_phase_metrics(F_hat, cells)
        split_rms_iter = float(gas_phase_iter.get("gas_phase_split_rms", 0.0))
        t_cap_hold_split = _gnr.update_t_cap_split_hysteresis(
            t_cap_hold_split,
            split_rms=split_rms_iter,
            merit_F=float(merit_F),
            merit_without_energy=float(max(norm_F, split_rms_iter)),
            engage_ratio=float(inner_split_dominant_merit_ratio),
            release_ratio=float(inner_t_cap_split_release_ratio),
        )
        _t_max_K_adaptive = _gnr.resolve_inner_t_max_step_K(
            ene_norm=float(ene_norm),
            tol_rms=float(tol_rms),
            split_rms=split_rms_iter,
            merit_F=float(merit_F),
            dominant_non_energy=float(max(gas_norm, sol_norm)),
            default_cap_K=float(inner_t_step_default_cap_K),
            split_dominant_cap_K=inner_split_dominant_t_step_cap_K,
            split_dominant_merit_ratio=float(inner_split_dominant_merit_ratio),
            line_search_split_enabled=line_search_merit_enabled,
            hold_split_cap=bool(t_cap_hold_split),
        )
        dT_so_far = 0.0
        if T_ref_abs.size:
            T_now = np.array(
                [float(c.T) if has_temperature_var(c) else float("nan") for c in cells],
                dtype=np.float64,
            )
            mask = np.isfinite(T_ref_abs) & np.isfinite(T_now)
            if np.any(mask):
                dT_so_far = float(np.max(np.abs(T_now[mask] - T_ref_abs[mask])))
        _t_max_K_adaptive = _gnr.resolve_t_step_with_abs_dT_budget(
            _t_max_K_adaptive,
            dT_so_far_K=dT_so_far,
            abs_dT_budget_K=inner_abs_dT_budget_K,
        )
        if record_inner_group_history:
            inner_group_history.append(
                {
                    "iter": int(it),
                    "merit": float(merit_F),
                    "rms_scaled": float(norm_F),
                    "gas_rms": float(gas_norm),
                    "solid_rms": float(sol_norm),
                    "energy_rms": float(ene_norm),
                    "split_rms": split_rms_iter,
                    "t_max_step_K": float(_t_max_K_adaptive),
                    "t_cap_hold_split": bool(t_cap_hold_split),
                    "dT_so_far_K": float(dT_so_far),
                    "converged_check": bool(
                        _gnr._inner_nr_check1_satisfied(
                            F_hat=F_hat,
                            cells=cells,
                            norm_F=norm_F,
                            x_rms=_gnr._state_rms_norm(x),
                            tol_rms=tol_rms,
                            inner_convergence=inner_convergence,
                            enforce_gas_phase_split_convergence=enforce_gas_phase_split_convergence,
                            line_search_merit_enabled=line_search_merit_enabled,
                            line_search_energy_merit_enabled=line_search_energy_merit_enabled,
                            tol_gas_phase_split_max=tol_gas_phase_split_max,
                        )
                    ),
                }
            )
        if verbose:
            print(
                f"  NR iter {it}: RMS={norm_F:.3e}  gas={gas_norm:.3e}  "
                f"solid={sol_norm:.3e}  energy={ene_norm:.3e}  "
                f"split={split_rms_iter:.3e}  t_cap={_t_max_K_adaptive:.1f}K"
            )

        x_rms = _gnr._state_rms_norm(x)
        if _gnr._inner_nr_check1_satisfied(
            F_hat=F_hat,
            cells=cells,
            norm_F=norm_F,
            x_rms=x_rms,
            tol_rms=tol_rms,
            inner_convergence=inner_convergence,
            enforce_gas_phase_split_convergence=enforce_gas_phase_split_convergence,
            line_search_merit_enabled=line_search_merit_enabled,
            line_search_energy_merit_enabled=line_search_energy_merit_enabled,
            tol_gas_phase_split_max=tol_gas_phase_split_max,
        ):
            converged = True
            break
        iter_attempts += 1

        # 构造 Jacobian。注意：这里 F 必须是与当前 x 对应的最新残差
        # 平稳阶段自适应 Jacobian lag：
        # 当上一轮 line-search 一次命中且阻尼不小，允许多复用一次 Jacobian
        # 以减少 FD 构造开销（不改变物理方程，仅改变线性化频率）。
        smooth_step = (
            len(line_search_trial_counts) > 0
            and line_search_trial_counts[-1] == 1
            and len(accepted_lambda_history) > 0
            and accepted_lambda_history[-1] is not None
            and float(accepted_lambda_history[-1]) >= 0.5
            and norm_F < max(10.0 * tol_rms, 0.2)
        )
        lag_steps_eff = max(lag_steps, 2) if smooth_step else lag_steps

        rebuild_jacobian = (
            cached_J is None
            or (it - cached_jac_iter) >= lag_steps_eff
        )
        if rebuild_jacobian:
            t0 = perf_counter()
            J, jac_meta = _gnr.build_jacobian_fd(
                x,
                F,
                cells,
                apply_bc_fn,
                scale,
                strategy=jacobian_strategy,
                apply_local_bc_fn=apply_local_bc_fn,
                affected_residual_cells_fn=affected_residual_cells_fn,
                connectivity_graph=connectivity_graph,
                verbose=verbose,
            )
            timing["jacobian_build_s"] += perf_counter() - t0
            cached_J = J
            cached_jac_meta = jac_meta
            cached_jac_iter = it
            counts["jacobian_rebuilds"] += 1
            counts["jacobian_cell_residual_calls"] += int(jac_meta.get("residual_cell_calls", 0))
            counts["jacobian_nnz_last"] = int(jac_meta.get("nnz", J.nnz))
            counts["jacobian_zero_cols_last"] = int(jac_meta.get("zero_cols", 0))
            counts["jacobian_zero_rows_last"] = int(jac_meta.get("zero_rows", 0))
            # 新 Jacobian 意味着全新牛顿方向；重置 lambda_seed 到初始值，
            # 以便线搜索从大步长开始探索，不沿用上次成功的小步长。
            lambda_seed = 1.0 if bool(prefer_full_step) else float(lambda_init)
        else:
            J = cached_J
            jac_meta = cached_jac_meta or {}
            counts["jacobian_reuses"] += 1

        # 使用稀疏解法求解 Newton 方向 Δx；对秩亏 Jacobian 做轻度对角正则化。
        t0 = perf_counter()
        if step_model_resolved == "lm":
            lm_mu_iter = lm_mu0_resolved * (lm_mu_growth_resolved ** max(consecutive_line_search_failures, 0))
            dx, linear_solver_used, linear_regularized, linear_diag = _gnr._solve_lm_step(
                J,
                F_hat,
                lm_mu=lm_mu_iter,
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
            counts["lm_steps"] = int(counts.get("lm_steps", 0)) + 1
        elif step_model_resolved == "ptc":
            ptc_alpha_iter = ptc_alpha0_resolved * (ptc_alpha_growth_resolved ** max(consecutive_line_search_failures, 0))
            dx, linear_solver_used, linear_regularized, linear_diag = _gnr._solve_ptc_step(
                J,
                F_hat,
                ptc_alpha=ptc_alpha_iter,
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
            counts["ptc_steps"] = int(counts.get("ptc_steps", 0)) + 1
        elif step_model_resolved == "equil_newton":
            dx, linear_solver_used, linear_regularized, linear_diag = _gnr._solve_equilibrated_newton_step(
                J,
                -F_hat,
                zero_rows=int(jac_meta.get("zero_rows", 0)),
                zero_cols=int(jac_meta.get("zero_cols", 0)),
                equil_iters=equil_iters_resolved,
                scale_clip=equil_scale_clip_resolved,
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
            counts["equil_steps"] = int(counts.get("equil_steps", 0)) + 1
        else:
            dx, linear_solver_used, linear_regularized, linear_diag = _gnr._solve_linear_step(
                J,
                -F_hat,
                zero_rows=int(jac_meta.get("zero_rows", 0)),
                zero_cols=int(jac_meta.get("zero_cols", 0)),
                structured_jacobian=jac_meta.get("structured_jacobian"),
                linear_solver_backend=str(linear_solver_backend),
                t_dof_indices=_t_dof_idx,
                t_max_step_K=_t_max_K_adaptive,
            )
        timing["linear_solve_s"] += perf_counter() - t0
        timing["band_lu_s"] = timing.get("band_lu_s", 0.0) + float(linear_diag.get("band_lu_s", 0.0))
        timing["side_update_s"] = timing.get("side_update_s", 0.0) + float(linear_diag.get("side_update_s", 0.0))
        if linear_regularized:
            counts["linear_regularized_solves"] += 1
        if linear_solver_used == "lstsq":
            counts["linear_lstsq_fallbacks"] += 1
        linear_solver_backend_last = str(linear_solver_used)
        counts["structured_fallbacks"] = int(counts.get("structured_fallbacks", 0)) + int(
            bool(linear_diag.get("fallback_used", False) and str(linear_solver_backend) == "structured_direct")
        )
        counts["nr_schur_size_last"] = int(linear_diag.get("schur_size", 0))

        # Energy-dominated plateau: optionally freeze gas holdup Newton components
        # so LS cannot trade energy↓ for dense–bubble split↑.
        if ls_split_constrained and _gnr._energy_dominates_merit(
            float(_gnr._residual_group_metrics(F_hat, cells).get("energy_rms", 0.0)),
            float(merit_F),
            energy_dominate_ratio=float(inner_split_dominant_merit_ratio),
        ):
            dx = _gnr._zero_gas_holdup_newton_step(dx, cells)
            counts["line_search_split_constrained_steps"] = int(
                counts.get("line_search_split_constrained_steps", 0)
            ) + 1

        # 分相切向：禁止同相增加 O₂∩合成气，使按相 fastox 近似恒等（handoff §5.51）。
        if clip_same_phase_ox_fuel:
            dx_before_ox = np.asarray(dx, dtype=np.float64)
            dx = _gnr._zero_same_phase_oxidizer_fuel_newton_step(
                dx, cells, allow_oxidizer_into_fuel=allow_ox_into_fuel
            )
            n_ox_zeroed = int(np.count_nonzero(np.abs(dx - dx_before_ox) > 0.0))
            if n_ox_zeroed > 0:
                counts["same_phase_oxidizer_fuel_clip_components"] = int(
                    counts.get("same_phase_oxidizer_fuel_clip_components", 0)
                ) + n_ox_zeroed

        # 对 dx 进行物理约束下的步长剪切
        dx_clipped = _gnr._clip_dx(
            dx,
            cells,
            ref_gas_mol_s,
            ref_solid_kg_s,
            t_step_limit_K=_t_max_K_adaptive,
            solid_ftb_relative_floor=bool(clip_solid_ftb_relative_floor),
        )
        clip_diag = _gnr._clip_dx_diagnostics(dx, dx_clipped, cells)
        newton_dir_diag = _gnr.diagnose_newton_step_consistency(
            J,
            F_hat,
            dx,
            dx_clipped,
        )
        clip_diag.update(
            {
                "newton_linear_defect_rel": float(newton_dir_diag.get("linear_defect_rel", np.nan)),
                "newton_clip_linear_defect_rel": float(
                    newton_dir_diag.get("clip_linear_defect_rel", np.nan)
                ),
                "newton_clip_displacement_rel": float(
                    newton_dir_diag.get("clip_displacement_rel", np.nan)
                ),
            }
        )

        # --- 阻尼线搜索 (Damped Newton) ---
        lambda_cap = 1.0 if bool(prefer_full_step) else max(float(lambda_init), lambda_min)
        # Hamel-style damped Newton starts fresh Newton directions from full step.
        # A line-search retry, however, is the same accepted state and the same
        # linearisation; restart it from the best smaller lambda instead of
        # rebuilding an identical FD Jacobian and trying lambda=1 again.
        if bool(prefer_full_step) and consecutive_line_search_failures == 0:
            lam = 1.0
        else:
            lam = float(np.clip(lambda_seed, lambda_min, lambda_cap))
        found_step = False
        trial_evals = 0
        best_trial_nt = float(np.inf)
        best_trial_merit = float(np.inf)
        best_trial_lambda: float | None = None
        nm_recent_max = float(max(history[-nonmonotone_window_resolved:])) if history else float(merit_F)
        nm_accept_upper = float(nm_recent_max * nonmonotone_relax_resolved)
        ls_trial_cap = (
            min(int(n_damp_halvings), line_search_max_trials_resolved)
            if line_search_max_trials_resolved > 0
            else int(n_damp_halvings)
        )
        curr_groups = _gnr._residual_group_metrics(F_hat, cells)
        curr_gas_phase = _gnr._residual_gas_phase_metrics(F_hat, cells)
        split_F_ls = float(curr_gas_phase.get("gas_phase_split_rms", 0.0))
        energy_F_ls = float(curr_groups.get("energy_rms", 0.0))
        # fastox 会按水力份额重分整条气体向量，含 NR 未打包组分（如 H2S）。
        # 拒绝试探时 pack/unpack 无法回滚这些隐藏状态，故在 LS 入口快照。
        excluded_gas_snap = (
            _snapshot_excluded_gas_holdups(cells) if ls_fast_ox_proj else None
        )
        syngas_base_mol = _bed_syngas_holdup_mol(cells) if ls_syngas_guard else 0.0
        co_base_mol = _bed_species_holdup_mol(cells, "CO") if ls_syngas_guard else 0.0

        for damp_step in range(ls_trial_cap):
            x_trial = x + lam * dx_clipped
            temp_bound_hits = 0
            # 对尝试步进行物理边界保护（必须正值，且 T 在范围内）
            offset = 0
            for cell in cells:
                nv = n_var(cell)
                if has_temperature_var(cell):
                    temp_local = temperature_offset(cell)
                    # 气/固质量流率 >= 0
                    x_trial[offset : offset + temp_local] = np.maximum(
                        x_trial[offset : offset + temp_local],
                        0.0,
                    )
                    # 温度保护
                    raw_T_trial = float(x_trial[offset + temp_local])
                    t_min, t_max = _gnr._temperature_bounds_for_cell(cell)
                    clipped_T_trial = float(np.clip(raw_T_trial, t_min, t_max))
                    temp_bound_hits += int(abs(clipped_T_trial - raw_T_trial) > 1e-12)
                    x_trial[offset + temp_local] = clipped_T_trial
                else:
                    x_trial[offset : offset + nv] = np.maximum(x_trial[offset : offset + nv], 0.0)
                offset += nv

            # Evaluate projected trial states even when temperature bounds are hit.
            # Rejecting all bound-hit trials can starve line-search (0 evaluations)
            # when one cell sits near a hard bound (e.g. T≈300 K).

            # Floor CO/H₂/CH₄ before fastox so Newton cannot wipe syngas, get a
            # fastox-repaired false merit drop, and lock the trajectory.
            if ls_syngas_floor:
                n_floored = _floor_syngas_holdups_in_state(
                    x_trial, x, cells, ls_syngas_floor_ratio
                )
                if n_floored > 0:
                    counts["line_search_syngas_floor_hits"] = int(
                        counts.get("line_search_syngas_floor_hits", 0)
                    ) + int(n_floored)

            # 快氧化流形投影：牛顿步易把 H₂ 与 O₂ 叠在同一相，裸 R12 打到 10⁸⁺。
            # per_phase：只消同相 overlap，保留气泡氧∩悬浮氢。合计 local_participating 为旧工程桥。
            if ls_fast_ox_proj:
                from types import SimpleNamespace

                from src.solvers.vorabrechnung.vorab_x0 import (
                    project_bed_holdup_fast_oxidation_closure,
                )

                if excluded_gas_snap is not None:
                    _restore_excluded_gas_holdups(cells, excluded_gas_snap)
                unpack_reactor(x_trial, cells)
                fuel_type = str(getattr(cells[0], "fuel_type", "coal") if cells else "coal")
                per_phase_limit = None
                if ls_per_phase_fastox and ls_per_phase_bed_limit > 0:
                    per_phase_limit = int(ls_per_phase_bed_limit)
                project_bed_holdup_fast_oxidation_closure(
                    cells,
                    SimpleNamespace(
                        vorab_transport_x0_fast_oxidation_closure_thesis=True,
                        vorab_transport_x0_char_oxidation_o2_closure_thesis=False,
                        fuel_type=fuel_type,
                    ),
                    bed_cell_limit=per_phase_limit,
                    phase_update="per_phase" if ls_per_phase_fastox else "local_participating",
                )
                x_trial = pack_reactor(cells)
                counts["line_search_fast_ox_projections"] = int(
                    counts.get("line_search_fast_ox_projections", 0)
                ) + 1
                # fastox can still burn CO/H₂; re-apply floor so merit is not
                # evaluated on a chemically wiped syngas state.
                if ls_syngas_floor:
                    n_floored = _floor_syngas_holdups_in_state(
                        x_trial, x, cells, ls_syngas_floor_ratio
                    )
                    if n_floored > 0:
                        unpack_reactor(x_trial, cells)
                        counts["line_search_syngas_floor_hits"] = int(
                            counts.get("line_search_syngas_floor_hits", 0)
                        ) + int(n_floored)

            # 计算尝试步的残差范数
            t0 = perf_counter()
            F_trial = _gnr.global_residual(x_trial, cells, apply_bc_fn)
            timing["line_search_residual_s"] += perf_counter() - t0
            counts["global_residual_calls"] += 1
            trial_evals += 1
            F_hat_trial = F_trial / scale
            nt = _gnr._rms_norm(F_hat_trial)
            merit_trial = _gnr._line_search_merit(
                F_hat_trial,
                cells,
                use_gas_phase_split=line_search_merit_enabled,
                use_energy=line_search_energy_merit_enabled,
            )
            split_trial = float(
                _gnr._residual_gas_phase_metrics(F_hat_trial, cells).get("gas_phase_split_rms", 0.0)
            )
            if np.isfinite(merit_trial) and merit_trial < best_trial_merit:
                best_trial_merit = float(merit_trial)
                best_trial_nt = float(nt)
                best_trial_lambda = float(lam)
            if verbose and trial_evals <= 3:
                print(
                    f"    LS trial {damp_step}: lam={lam:.4e}  rms={nt:.4e}  "
                    f"merit={merit_trial:.4e}  (merit_F={merit_F:.4e})"
                )

            # Guard: clipped Newton can drive CO → 0 while H₂ remains; fastox then
            # repairs R12 residuals so merit falls and the ruined trial is accepted.
            # Check CO alone (LU parity species) and total syngas.
            # Allow a tiny absolute slack so the packed-state floor (exact ratio×base)
            # is not rejected by floating-point roundoff.
            syngas_collapsed = False
            if ls_syngas_guard:
                co_floor = ls_syngas_ratio * co_base_mol
                syn_floor = ls_syngas_ratio * syngas_base_mol
                atol_co = 1.0e-9 * max(co_base_mol, 1.0)
                atol_syn = 1.0e-9 * max(syngas_base_mol, 1.0)
                if co_base_mol > 1.0e-6:
                    co_trial_mol = _bed_species_holdup_mol(cells, "CO")
                    if co_trial_mol + atol_co < co_floor:
                        syngas_collapsed = True
                if (not syngas_collapsed) and syngas_base_mol > 1.0e-6:
                    syngas_trial_mol = _bed_syngas_holdup_mol(cells)
                    if syngas_trial_mol + atol_syn < syn_floor:
                        syngas_collapsed = True
                if syngas_collapsed:
                    counts["line_search_syngas_collapse_rejects"] = int(
                        counts.get("line_search_syngas_collapse_rejects", 0)
                    ) + 1

            accepted_nonmonotone = bool(
                nonmonotone_enabled_resolved
                and consecutive_line_search_failures > 0
                and merit_trial >= merit_F
                and np.isfinite(merit_trial)
                and merit_trial <= nm_accept_upper
            )
            trial_ok = (not syngas_collapsed) and _gnr._line_search_trial_acceptable(
                merit_trial=float(merit_trial),
                merit_F=float(merit_F),
                split_trial=float(split_trial),
                split_F=float(split_F_ls),
                energy_F=float(energy_F_ls),
                use_gas_phase_split=line_search_merit_enabled,
                use_energy=line_search_energy_merit_enabled,
                accepted_nonmonotone=accepted_nonmonotone,
                energy_dominate_ratio=float(inner_split_dominant_merit_ratio),
            )
            refresh_aware_accepted = False
            merit_pre_refresh = float(merit_trial)
            merit_post_refresh = float(merit_trial)
            # Refresh-aware joint accept: energy↓ but max-merit↑ (split climb) may
            # still be OK after hydro refresh if post merit/split are not worse.
            if (
                (not trial_ok)
                and (not syngas_collapsed)
                and ls_refresh_aware
                and line_search_merit_enabled
                and line_search_energy_merit_enabled
            ):
                energy_trial = float(
                    _gnr._residual_group_metrics(F_hat_trial, cells).get("energy_rms", 0.0)
                )
                if _gnr._energy_dominates_merit(
                    float(energy_F_ls),
                    float(merit_F),
                    energy_dominate_ratio=float(inner_split_dominant_merit_ratio),
                ) and float(energy_trial) < float(energy_F_ls):
                    from src.solvers.vorabrechnung import (
                        backup_vorabrechnung_hydro_freeze,
                        refresh_cells_hydrodynamics_for_frozen_inner,
                        restore_vorabrechnung_hydro_freeze,
                    )

                    hydro_backup = backup_vorabrechnung_hydro_freeze(cells)
                    refresh_cells_hydrodynamics_for_frozen_inner(cells)
                    counts["line_search_refresh_aware_probes"] = int(
                        counts.get("line_search_refresh_aware_probes", 0)
                    ) + 1
                    t0 = perf_counter()
                    F_post = _gnr.global_residual(x_trial, cells, apply_bc_fn)
                    timing["line_search_residual_s"] += perf_counter() - t0
                    counts["global_residual_calls"] += 1
                    F_hat_post = F_post / scale
                    merit_post = _gnr._line_search_merit(
                        F_hat_post,
                        cells,
                        use_gas_phase_split=line_search_merit_enabled,
                        use_energy=line_search_energy_merit_enabled,
                    )
                    split_post = float(
                        _gnr._residual_gas_phase_metrics(F_hat_post, cells).get(
                            "gas_phase_split_rms", 0.0
                        )
                    )
                    if _gnr._refresh_aware_joint_acceptable(
                        energy_F=float(energy_F_ls),
                        energy_trial=float(energy_trial),
                        merit_F=float(merit_F),
                        merit_post=float(merit_post),
                        split_F=float(split_F_ls),
                        split_post=float(split_post),
                        energy_dominate_ratio=float(inner_split_dominant_merit_ratio),
                    ):
                        trial_ok = True
                        refresh_aware_accepted = True
                        F_trial = F_post
                        F_hat_trial = F_hat_post
                        nt = _gnr._rms_norm(F_hat_trial)
                        merit_trial = float(merit_post)
                        split_trial = float(split_post)
                        merit_post_refresh = float(merit_post)
                        counts["line_search_refresh_aware_accepts"] = int(
                            counts.get("line_search_refresh_aware_accepts", 0)
                        ) + 1
                    else:
                        restore_vorabrechnung_hydro_freeze(cells, hydro_backup)
                        counts["line_search_refresh_aware_rejects"] = int(
                            counts.get("line_search_refresh_aware_rejects", 0)
                        ) + 1

            if trial_ok:
                if verbose and damp_step > 0:
                    print(
                        f"    Line search OK at step {damp_step}: lam={lam:.4f}, "
                        f"rms={nt:.3e}, merit={merit_trial:.3e}"
                    )
                found_step = True
                consecutive_line_search_failures = 0
                prev_norm_F = float(norm_F)
                x = pack_reactor(cells)
                F = F_trial
                F_hat = F_hat_trial
                norm_F = nt
                accepted_merit = float(merit_trial)
                merit_F = accepted_merit
                merit_pre_refresh = float(merit_pre_refresh)
                merit_post_refresh = float(merit_post_refresh)
                # P1：默认 outer_fixed 不走此支路。accept_refresh opt-in 时必须重算 F。
                # Refresh-aware 接受已刷过 h 并持有 post residual，勿再刷一次。
                if refresh_hydro_on_accept and (not refresh_aware_accepted):
                    from src.solvers.vorabrechnung import (
                        refresh_cells_hydrodynamics_for_frozen_inner,
                    )

                    refresh_cells_hydrodynamics_for_frozen_inner(cells)
                    counts["hydrodynamics_refresh_on_accept"] = int(
                        counts.get("hydrodynamics_refresh_on_accept", 0)
                    ) + 1
                    # 水力刷新改变 K_bd/ε；必须用刷新后残差作为下一步 RHS，并丢弃 J。
                    F, F_hat, norm_F, merit_F = _recompute_residual_after_hydro_refresh(
                        x=x,
                        cells=cells,
                        apply_bc_fn=apply_bc_fn,
                        scale=scale,
                        line_search_merit_enabled=line_search_merit_enabled,
                        line_search_energy_merit_enabled=line_search_energy_merit_enabled,
                        timing=timing,
                        counts=counts,
                    )
                    merit_post_refresh = float(merit_F)
                    x = pack_reactor(cells)
                    cached_J = None
                    cached_jac_meta = None
                    cached_jac_iter = -1
                elif refresh_aware_accepted:
                    cached_J = None
                    cached_jac_meta = None
                    cached_jac_iter = -1
                counts["line_search_evaluations"] += int(trial_evals)
                counts["line_search_backtracks"] += max(int(trial_evals) - 1, 0)
                if accepted_nonmonotone:
                    counts["nonmonotone_accepts"] = int(counts.get("nonmonotone_accepts", 0)) + 1
                accepted_lambda_history.append(float(lam))
                line_search_trial_counts.append(int(trial_evals))
                clip_history.append(
                    {
                        "iter": int(it + 1),
                        "accepted_lambda": float(lam),
                        "line_search_trials": int(trial_evals),
                        "line_search_failed": False,
                        "nonmonotone_accepted": bool(accepted_nonmonotone),
                        "refresh_aware_accepted": bool(refresh_aware_accepted),
                        "best_trial_rms_scaled": float(best_trial_nt),
                        "best_trial_merit_scaled": float(best_trial_merit),
                        "best_trial_lambda": best_trial_lambda,
                        "merit_pre_refresh": float(merit_pre_refresh),
                        "merit_post_refresh": float(merit_post_refresh),
                        **clip_diag,
                    }
                )
                lambda_seed = _gnr._next_lambda_seed(
                    float(lambda_cap),
                    float(lam),
                    int(trial_evals),
                    float(accepted_merit / max(prev_norm_F, 1.0e-30)),
                )
                break
            
            lam *= 0.5
            if lam < lambda_min:
                break

        if not found_step:
            if excluded_gas_snap is not None:
                _restore_excluded_gas_holdups(cells, excluded_gas_snap)
            if line_search_max_trials_resolved > 0 and int(trial_evals) >= int(ls_trial_cap):
                counts["line_search_cap_hits"] = int(counts.get("line_search_cap_hits", 0)) + 1
            counts["line_search_evaluations"] += int(trial_evals)
            counts["line_search_backtracks"] += max(int(trial_evals) - 1, 0)
            counts["line_search_failures"] += 1
            accepted_lambda_history.append(None)
            line_search_trial_counts.append(int(trial_evals))
            clip_history.append(
                {
                    "iter": int(it + 1),
                    "accepted_lambda": None,
                    "line_search_trials": int(trial_evals),
                    "line_search_failed": True,
                    "best_trial_rms_scaled": float(best_trial_nt),
                    "best_trial_merit_scaled": float(best_trial_merit),
                    "best_trial_lambda": best_trial_lambda,
                    **clip_diag,
                }
            )
            if verbose:
                print(f"  NR iter {it}: Line search failed to reduce merit. Stopping.")
            consecutive_line_search_failures += 1
            retry_allowed = consecutive_line_search_failures <= 1 and (it + 1) < max_iter
            if retry_allowed:
                counts["line_search_retries"] += 1
                lam_retry_base = float(best_trial_lambda) if best_trial_lambda is not None else float(lambda_seed)
                lambda_seed = float(np.clip(0.5 * lam_retry_base, lambda_min, lambda_cap))
                # 恢复 cells 到接受状态（line search trial 会把 cells 留在最后一次试探）
                _ = _gnr.global_residual(x, cells, apply_bc_fn)
                continue

            if not bool(allow_gd_fallback):
                # Hamel-aligned damped Newton path: stop after line-search failure
                # (after one smaller-lambda retry) instead of switching to a
                # gradient-descent surrogate direction.
                _ = _gnr.global_residual(x, cells, apply_bc_fn)
                counts["global_residual_calls"] += 1
                break

            # --- 梯度下降兜底 (Gradient-descent fallback) ---
            # 当 NR 线搜索彻底失败（lstsq 方向在当前状态为反下降方向）时，
            # 使用 dx_gd = -J^T F_hat 作为替代方向。
            # 该方向保证是 ||F_hat||^2 的下降方向：
            #   d/dλ ||F_hat(x + λ dx_gd)||^2|_{λ=0} = -2||J^T F_hat||^2 ≤ 0
            # 即使 Jacobian 近奇异（holdup_transport + frozen Vorabrechnung），
            # 只要 J^T F_hat ≠ 0，梯度方向就是下降方向。
            # 恢复 cells 到接受状态（line search trial 会把 cells 留在最后一次试探状态）
            _ = _gnr.global_residual(x, cells, apply_bc_fn)
            J_arr = J.toarray() if sp.issparse(J) else np.asarray(J)
            dx_gd = -(J_arr.T @ F_hat)
            # 当能量已满足时，清零 T 分量：任何 T 扰动都会破坏已收敛的能量方程。
            _energy_satisfied_for_gd = ene_norm < tol_rms * 1e-3 and _energy_frac < 1e-3
            if _energy_satisfied_for_gd and _t_dof_idx is not None:
                dx_gd[_t_dof_idx] = 0.0
            gd_norm = float(np.linalg.norm(dx_gd))
            if verbose:
                print(f"  GD fallback: gd_norm={gd_norm:.3e}, norm_F={norm_F:.3e}, freeze_T={_energy_satisfied_for_gd}")
            if gd_norm > 1e-10:
                # 将方向归一化后按气相参考量缩放，确保步长物理合理
                gd_scale = 0.1 * float(ref_gas_mol_s) / gd_norm
                dx_gd_step = gd_scale * dx_gd
                if clip_same_phase_ox_fuel:
                    dx_gd_step = _gnr._zero_same_phase_oxidizer_fuel_newton_step(
                        dx_gd_step, cells, allow_oxidizer_into_fuel=allow_ox_into_fuel
                    )
                dx_gd_clipped = _gnr._clip_dx(
                    dx_gd_step,
                    cells,
                    ref_gas_mol_s,
                    ref_solid_kg_s,
                    zero_empty_solid=True,
                    t_step_limit_K=_t_max_K_adaptive,
                    solid_ftb_relative_floor=bool(clip_solid_ftb_relative_floor),
                )
                if verbose:
                    print(f"    gd_scale={gd_scale:.3e}, max|dx_gd_clip|={np.max(np.abs(dx_gd_clipped)):.3e}")
                clip_diag_gd = _gnr._clip_dx_diagnostics(dx_gd_clipped, dx_gd_clipped, cells)
                lam_gd = float(lambda_init)
                found_gd = False
                gd_trial_evals = 0
                # 追踪线搜索中最优试探（用于非单调兜底）
                best_gd_nt = float(np.inf)
                best_gd_x: np.ndarray | None = None
                best_gd_F: np.ndarray | None = None
                best_gd_F_hat: np.ndarray | None = None
                best_gd_lam: float = float(lambda_init)

                for _ in range(n_damp_halvings * 2):
                    x_trial_gd = x + lam_gd * dx_gd_clipped
                    offset = 0
                    gd_temp_hits = 0
                    for cell in cells:
                        nv = n_var(cell)
                        if has_temperature_var(cell):
                            temp_local = temperature_offset(cell)
                            x_trial_gd[offset : offset + temp_local] = np.maximum(
                                x_trial_gd[offset : offset + temp_local], 0.0
                            )
                            raw_T = float(x_trial_gd[offset + temp_local])
                            clipped_T = float(np.clip(raw_T, 300.0, 2500.0))
                            gd_temp_hits += int(abs(clipped_T - raw_T) > 1e-12)
                            x_trial_gd[offset + temp_local] = clipped_T
                        else:
                            x_trial_gd[offset : offset + nv] = np.maximum(
                                x_trial_gd[offset : offset + nv], 0.0
                            )
                        offset += nv
                    if gd_temp_hits > 0:
                        if lam_gd <= lambda_min:
                            break
                        lam_gd *= 0.5
                        continue
                    F_gd = _gnr.global_residual(x_trial_gd, cells, apply_bc_fn)
                    counts["global_residual_calls"] += 1
                    gd_trial_evals += 1
                    F_hat_gd = F_gd / scale
                    nt_gd = _gnr._rms_norm(F_hat_gd)
                    if np.isfinite(nt_gd) and nt_gd < best_gd_nt:
                        best_gd_nt = float(nt_gd)
                        best_gd_x = x_trial_gd.copy()
                        best_gd_F = F_gd.copy()
                        best_gd_F_hat = F_hat_gd.copy()
                        best_gd_lam = float(lam_gd)
                    if verbose:
                        print(f"    GD lam={lam_gd:.4f}  rms={nt_gd:.4e}  {'<<ACCEPT' if np.isfinite(nt_gd) and nt_gd < norm_F else ''}")
                        if gd_trial_evals == 1:
                            # First GD trial: diagnose which residual component is largest
                            big_idx = int(np.argmax(np.abs(F_hat_gd)))
                            print(
                                f"    [GD diag] biggest F_hat: idx={big_idx} val={F_hat_gd[big_idx]:.3e} "
                                f"x={x_trial_gd[big_idx]:.3e} prev_F_hat={F_hat[big_idx]:.3e}"
                            )
                            offsets_diag = cell_offsets(cells)
                            cell_idx = max(0, int(np.searchsorted(offsets_diag, big_idx, side="right") - 1)) if cells else 0
                            print(f"    [GD diag] cell≈{cell_idx}, change_at_idx={x_trial_gd[big_idx]-x[big_idx]:.3e}")
                    if np.isfinite(nt_gd) and nt_gd < norm_F:
                        # 严格下降
                        x = x_trial_gd
                        F = F_gd
                        F_hat = F_hat_gd
                        norm_F = nt_gd
                        merit_F = _gnr._line_search_merit(
                            F_hat_gd,
                            cells,
                            use_gas_phase_split=line_search_merit_enabled,
                            use_energy=line_search_energy_merit_enabled,
                        )
                        found_gd = True
                        consecutive_line_search_failures = 0
                        counts.setdefault("gd_fallback_accepts", 0)
                        counts["gd_fallback_accepts"] += 1
                        accepted_lambda_history.append(float(lam_gd))
                        line_search_trial_counts.append(int(gd_trial_evals))
                        clip_history.append(
                            {
                                "iter": int(it + 1),
                                "accepted_lambda": float(lam_gd),
                                "line_search_trials": int(gd_trial_evals),
                                "line_search_failed": False,
                                "gd_fallback": True,
                                "best_trial_rms_scaled": float(nt_gd),
                                "best_trial_lambda": float(lam_gd),
                                **clip_diag_gd,
                            }
                        )
                        lambda_seed = float(lambda_min)
                        break
                    lam_gd *= 0.5
                    if lam_gd < lambda_min:
                        break

                if not found_gd and best_gd_x is not None:
                    # 非单调有界兜底：若 GD 最优试探比当前 rms 仅高 ≤5%，
                    # 则接受之（有界 non-monotone，防止发散），以便外层继续推进。
                    nm_rel_tol = 1.05
                    if best_gd_nt < nm_rel_tol * norm_F:
                        _old_norm_F = norm_F
                        x = best_gd_x
                        F = best_gd_F
                        F_hat = best_gd_F_hat  # type: ignore[assignment]
                        norm_F = best_gd_nt
                        merit_F = _gnr._line_search_merit(
                            F_hat,
                            cells,
                            use_gas_phase_split=line_search_merit_enabled,
                            use_energy=line_search_energy_merit_enabled,
                        )
                        found_gd = True
                        consecutive_line_search_failures = 0
                        counts.setdefault("gd_nm_accepts", 0)
                        counts["gd_nm_accepts"] += 1
                        accepted_lambda_history.append(float(best_gd_lam))
                        line_search_trial_counts.append(int(gd_trial_evals))
                        clip_history.append(
                            {
                                "iter": int(it + 1),
                                "accepted_lambda": float(best_gd_lam),
                                "line_search_trials": int(gd_trial_evals),
                                "line_search_failed": False,
                                "gd_fallback": True,
                                "gd_nonmonotone": True,
                                "best_trial_rms_scaled": float(best_gd_nt),
                                "best_trial_lambda": float(best_gd_lam),
                                **clip_diag_gd,
                            }
                        )
                        if verbose:
                            print(
                                f"  GD non-monotone accept: rms {norm_F:.4e} "
                                f"(ratio={norm_F / _old_norm_F:.3f})"
                            )
                        lambda_seed = float(lambda_min)
                    if found_gd:
                        if refresh_hydro_on_accept:
                            from src.solvers.vorabrechnung import (
                                refresh_cells_hydrodynamics_for_frozen_inner,
                            )

                            refresh_cells_hydrodynamics_for_frozen_inner(cells)
                            counts["hydrodynamics_refresh_on_accept"] = int(
                                counts.get("hydrodynamics_refresh_on_accept", 0)
                            ) + 1
                            F, F_hat, norm_F, merit_F = _recompute_residual_after_hydro_refresh(
                                x=x,
                                cells=cells,
                                apply_bc_fn=apply_bc_fn,
                                scale=scale,
                                line_search_merit_enabled=line_search_merit_enabled,
                                line_search_energy_merit_enabled=line_search_energy_merit_enabled,
                                timing=timing,
                                counts=counts,
                            )
                            x = pack_reactor(cells)
                            cached_J = None
                            cached_jac_meta = None
                            cached_jac_iter = -1
                        history.append(merit_F)
                        continue  # 继续下一个 NR 迭代
            break

        history.append(merit_F)

    # 结果解包回反应器对象
    t0 = perf_counter()
    F_final = _gnr.global_residual(x, cells, apply_bc_fn)
    timing["final_residual_s"] += perf_counter() - t0
    counts["global_residual_calls"] += 1
    x = pack_reactor(cells)
    final_norm = float(np.linalg.norm(F_final))
    F_hat_final = np.asarray(F_final, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    final_group_metrics = _gnr._residual_group_metrics(F_hat_final, cells)
    final_gas_phase_metrics = _gnr._residual_gas_phase_metrics(F_hat_final, cells)
    final_merit = _gnr._line_search_merit(
        F_hat_final,
        cells,
        use_gas_phase_split=bool(
            line_search_use_gas_phase_split_merit and enforce_gas_phase_split_convergence
        ),
        use_energy=bool(line_search_use_energy_merit),
    )
    timing["total_s"] = perf_counter() - solve_started

    return {
        "converged": converged,
        "n_iter": len(history),
        "residual": final_norm,
        "rms_scaled_final": float(final_merit),
        "rms_scaled_final_last_iter": float(norm_F),
        "rms_scaled_gas_final": final_group_metrics["gas_rms"],
        "rms_scaled_solid_final": final_group_metrics["solid_rms"],
        "rms_scaled_energy_final": final_group_metrics["energy_rms"],
        "rms_scaled_gas_combined_final": final_gas_phase_metrics["gas_combined_rms"],
        "rms_scaled_gas_phase_split_final": final_gas_phase_metrics["gas_phase_split_rms"],
        "rms_scaled_component_max_final": max(
            final_group_metrics["gas_rms"],
            final_group_metrics["solid_rms"],
            final_group_metrics["energy_rms"],
        ),
        "max_abs_scaled_gas_final": final_group_metrics["gas_max_abs"],
        "max_abs_scaled_solid_final": final_group_metrics["solid_max_abs"],
        "max_abs_scaled_energy_final": final_group_metrics["energy_max_abs"],
        "max_abs_scaled_gas_combined_final": final_gas_phase_metrics["gas_combined_max_abs"],
        "max_abs_scaled_gas_phase_split_final": final_gas_phase_metrics["gas_phase_split_max_abs"],
        "max_abs_scaled_final": max(
            final_group_metrics["gas_max_abs"],
            final_group_metrics["solid_max_abs"],
            final_group_metrics["energy_max_abs"],
        ),
        "norm_history": history,
        "n_newton_iters_attempted": int(iter_attempts),
        "jacobian_strategy": jacobian_strategy,
        "nr_layout_audit": _gnr._nr_layout_audit(cells, F_final),
        "jacobian_structure": (
            cached_jac_meta.get("jacobian_structure")
            if cached_jac_meta is not None
            else {
                "main_chain_blocks": int(len(cells)),
                "band_block_count": 0,
                "side_element_count": 0,
                "side_tail_block_count": 0,
                "side_head_block_count": 0,
                "unexpected_block_count": 0,
                "unexpected_blocks": [],
                "structure_validation_ok": True,
            }
        ),
        "linear_solver_backend": str(linear_solver_backend),
        "linear_solver_backend_last": linear_solver_backend_last,
        "step_model": step_model_resolved,
        "lm_mu0": lm_mu0_resolved,
        "lm_mu_growth": lm_mu_growth_resolved,
        "ptc_alpha0": ptc_alpha0_resolved,
        "ptc_alpha_growth": ptc_alpha_growth_resolved,
        "equil_iters": equil_iters_resolved,
        "equil_scale_clip": equil_scale_clip_resolved,
        "line_search_max_trials": line_search_max_trials_resolved,
        "nonmonotone_enabled": nonmonotone_enabled_resolved,
        "nonmonotone_window": nonmonotone_window_resolved,
        "nonmonotone_relax": nonmonotone_relax_resolved,
        "line_search_use_gas_phase_split_merit": line_search_merit_enabled,
        "line_search_refresh_aware_accept": ls_refresh_aware,
        "line_search_split_constrained_step": ls_split_constrained,
        "refresh_hydrodynamics_on_accepted_step": refresh_hydro_on_accept,
        "inner_t_cap_hold_split_dominant": bool(t_cap_hold_split),
        "jacobian_lag_steps": int(lag_steps),
        "timing": timing,
        "counts": counts,
        "accepted_lambda_history": accepted_lambda_history,
        "line_search_trial_counts": line_search_trial_counts,
        "clip_history": clip_history,
        "inner_group_history": inner_group_history,
    }
