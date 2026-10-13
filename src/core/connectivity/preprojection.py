"""Connectivity submodule — see package ``__init__`` for public API."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from src.core.cell import Cell, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_gas_enthalpy_flow, calc_solid_enthalpy_flow
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX, cp_ash, cp_char, cp_sand
from src.core.connectivity_graph import (
    ConnectivityBlock,
    ConnectivityEdge,
    ConnectivityGraph,
    bed_top_freeboard_ab_source_block_id,
    block_id_bed_index,
    build_connectivity_graph_from_reactor_params,
    cyclone_inlet_source_block_id,
    axial_upstream_block_id,
    freeboard_chain_inlet_source_block_id,
    infer_stream_class,
    ordered_freeboard_block_ids,
    resolve_bed_neighbor_block_ids,
    resolve_registry_cell,
    solid_ab_freeboard_downstream_block_id,
    solid_ab_source_block_id,
    solid_auf_source_block_id,
)
from src.solvers.nr_indexing import (
    MAJOR_NR_GAS_SPECIES,
    nr_dense_only_pyrolysis_support_species,
    nr_projection_gas_species,
)

_F_W = 0.25
_ZETA_W = 0.40


def _preproj_gas_reaction_rate(cfg: Any) -> float:
    """Gas/energy preprojection: freeze main kinetics (Hamel Vorab Startwert semantics)."""
    return float(getattr(cfg, "preproj_gas_reaction_rate_multiplier_thesis", 0.0))


def _preproj_phase_split_reaction_rate(cfg: Any) -> float:
    """Phase-split LSQ uses full kinetics so closure matches inner NR merit."""
    return float(getattr(cfg, "preproj_phase_split_reaction_rate_multiplier_thesis", 1.0))

from src.core.preprojection_engine import (
    LsqPreprojectionSpec,
    PreprojectionBcContext,
    preproj_cell_energy_residual_abs as _preproj_cell_energy_residual_abs,
    preproj_cell_residual_vector as _preproj_cell_residual_vector,
    preproj_finalize_bc as _preproj_finalize_bc,
    preproj_refresh_trial_bc as _preproj_refresh_trial_bc,
    run_improvement_gated_lsq,
    run_lsq_preprojection,
    run_relax_schedule,
    scaled_rms,
)
from src.core.connectivity.joint_gate import (
    _blend_bed0_state_from_snapshot,
    _bottom_joint_global_gate_pass,
    _bottom_joint_relax_candidates,
    _bottom_neighbor_energy_residual_abs,
    _measure_nr_init_scaled_metrics,
    resolve_preinner_tier_plan,
)
from src.core.connectivity.solid_transport import (
    max_upper_bed_solid_auf_flux_ratio_error,
    max_upper_bed_solid_stream_imbalance,
    reconcile_upper_bed_solid_holdup_for_init,
)

def _preinner_split_only_fallback_gate_pass(
    metrics_before: dict[str, float],
    metrics_after: dict[str, float],
    cfg: Any,
) -> bool:
    """Strict merit gate for T2-only retry after combined pre-inner rollback."""
    return _bottom_joint_global_gate_pass(metrics_before, metrics_after, cfg)


def preproject_bottom_primary_gas_state_to_exchange_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    improvement_factor: float = 0.75,
) -> bool:
    """Preproject bed0 primary gas split against the local K_bd exchange closure.

    Hamel gives the Vorabrechnung module responsibility for hydrodynamics and
    the outer Abgleich with the cell model, but does not prescribe an explicit
    start-value projection.  This bridge is therefore deliberately narrow:
    only O2/H2O/N2 in the bottom bed cell are repartitioned, each species'
    total molar flow is conserved, and reaction products are left untouched.
    """
    if not cells:
        return False
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    gas_rm = _preproj_gas_reaction_rate(cfg)

    idxs = [GAS_SPECIES_INDEX[sp] for sp in ("O2", "H2O", "N2")]
    totals = np.array([float(max(cell.N_d[j] + cell.N_b[j], 0.0)) for j in idxs], dtype=np.float64)
    if np.any(~np.isfinite(totals)) or np.any(totals <= 1e-10):
        return False

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()

    def _invalidate_gas_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _set_dense_values(values: np.ndarray) -> None:
        clipped = np.clip(np.asarray(values, dtype=np.float64), 1e-12, totals - 1e-12)
        for value, total, j in zip(clipped, totals, idxs):
            cell.N_d[j] = float(value)
            cell.N_b[j] = float(total - value)
        _invalidate_gas_state_cache()

    def _primary_dense_residual(values: np.ndarray) -> np.ndarray:
        _set_dense_values(values)
        res = _preproj_cell_residual_vector(cell, "gas", reaction_rate_multiplier=gas_rm)
        return res[idxs].copy()

    x0 = np.array([float(cell.N_d[j]) for j in idxs], dtype=np.float64)
    initial_res = _primary_dense_residual(x0)
    initial_norm = float(np.linalg.norm(initial_res, ord=2))
    if not np.isfinite(initial_norm) or initial_norm <= 1e-12:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()
        return False

    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            _primary_dense_residual,
            x0=np.clip(x0, 1e-12, totals - 1e-12),
            bounds=(np.full_like(totals, 1e-12), totals - 1e-12),
            max_nfev=12,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )
        candidate = np.asarray(sol.x, dtype=np.float64)
        final_res = _primary_dense_residual(candidate)
    except Exception:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()
        cell.residuals()
        return False

    final_norm = float(np.linalg.norm(final_res, ord=2))
    if np.isfinite(final_norm) and final_norm <= improvement_factor * initial_norm:
        _set_dense_values(candidate)
        cell.residuals()
        return True

    cell.N_d[:] = original_d
    cell.N_b[:] = original_b
    _invalidate_gas_state_cache()
    cell.residuals()
    return False

def preproject_bottom_major_gas_state_to_local_balance(
    cells: list[Cell],
    cfg: Any,
    *,
    improvement_factor: float = 0.75,
) -> bool:
    """Preproject bed0 major gas x0 against the local two-phase gas balances.

    This is a start-value bridge for the bottom oxidation/devolatilization zone:
    primary species totals (O2/H2O/N2) remain fixed, while major kinetic products
    (CO/H2/CH4/CO2) may receive dense and bubble outlet support.  Tar species are
    intentionally excluded because the current tar source is a separate lumped
    pyrolysis carrier and including it worsened the O2/H2O split in static probes.
    """
    if not cells:
        return False
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    gas_rm = _preproj_gas_reaction_rate(cfg)

    primary_names = ("O2", "H2O", "N2")
    product_names = ("CO", "H2", "CH4", "CO2")
    primary_idxs = [GAS_SPECIES_INDEX[sp] for sp in primary_names]
    product_idxs = [GAS_SPECIES_INDEX[sp] for sp in product_names]
    objective_idxs = primary_idxs + product_idxs
    objective_rows = objective_idxs + [j + cell.N_d.shape[0] for j in objective_idxs]
    primary_totals = np.array(
        [float(max(cell.N_d[j] + cell.N_b[j], 0.0)) for j in primary_idxs],
        dtype=np.float64,
    )
    if np.any(~np.isfinite(primary_totals)) or np.any(primary_totals <= 1e-10):
        return False

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()

    def _invalidate_gas_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _set_values(values: np.ndarray) -> None:
        x = np.asarray(values, dtype=np.float64)
        k = 0
        for value, total, j in zip(x[k : k + len(primary_idxs)], primary_totals, primary_idxs):
            dense = float(np.clip(value, 1e-12, total - 1e-12))
            cell.N_d[j] = dense
            cell.N_b[j] = float(total - dense)
        k += len(primary_idxs)
        for value, j in zip(x[k : k + len(product_idxs)], product_idxs):
            cell.N_d[j] = max(float(value), 1e-12)
        k += len(product_idxs)
        for value, j in zip(x[k : k + len(product_idxs)], product_idxs):
            cell.N_b[j] = max(float(value), 1e-12)
        _invalidate_gas_state_cache()

    def _residual_objective(values: np.ndarray) -> np.ndarray:
        _set_values(values)
        res = _preproj_cell_residual_vector(cell, "gas", reaction_rate_multiplier=gas_rm)
        return res[objective_rows].copy()

    initial_res = _preproj_cell_residual_vector(cell, "gas", reaction_rate_multiplier=gas_rm)
    initial_gas_max = float(
        max(
            np.max(np.abs(initial_res[: cell.N_d.shape[0]])),
            np.max(np.abs(initial_res[cell.N_d.shape[0] : 2 * cell.N_d.shape[0]])),
        )
    )
    if not np.isfinite(initial_gas_max) or initial_gas_max <= 1e-12:
        return False

    x0 = np.concatenate(
        (
            np.array([float(cell.N_d[j]) for j in primary_idxs], dtype=np.float64),
            np.array([max(float(cell.N_d[j]), 1e-9) for j in product_idxs], dtype=np.float64),
            np.array([max(float(cell.N_b[j]), 1e-9) for j in product_idxs], dtype=np.float64),
        )
    )
    lower = np.full_like(x0, 1e-12)
    upper = np.concatenate(
        (
            primary_totals - 1e-12,
            np.full(len(product_idxs), 80.0, dtype=np.float64),
            np.full(len(product_idxs), 80.0, dtype=np.float64),
        )
    )

    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            _residual_objective,
            x0=np.clip(x0, lower, upper),
            bounds=(lower, upper),
            max_nfev=20,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )
        candidate = np.asarray(sol.x, dtype=np.float64)
        _set_values(candidate)
        final_res = _preproj_cell_residual_vector(cell, "gas", reaction_rate_multiplier=gas_rm)
    except Exception:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()
        cell.residuals()
        return False

    final_gas_max = float(
        max(
            np.max(np.abs(final_res[: cell.N_d.shape[0]])),
            np.max(np.abs(final_res[cell.N_d.shape[0] : 2 * cell.N_d.shape[0]])),
        )
    )
    if np.isfinite(final_gas_max) and final_gas_max <= improvement_factor * initial_gas_max:
        return True

    cell.N_d[:] = original_d
    cell.N_b[:] = original_b
    _invalidate_gas_state_cache()
    cell.residuals()
    return False

def _bed_major_gas_phase_split_vector(cell: Cell, *, reaction_rate_multiplier: float | None = None) -> np.ndarray:
    """Dense minus bubble residuals for the seven major NR gas species."""
    res = _preproj_cell_residual_vector(
        cell,
        "gas",
        reaction_rate_multiplier=reaction_rate_multiplier,
    )
    n_gas = cell.N_d.shape[0]
    return np.array(
        [float(res[GAS_SPECIES_INDEX[sp]] - res[n_gas + GAS_SPECIES_INDEX[sp]]) for sp in MAJOR_NR_GAS_SPECIES],
        dtype=np.float64,
    )

def _bed_major_gas_phase_split_rms(cell: Cell, *, reaction_rate_multiplier: float | None = None) -> float:
    split = _bed_major_gas_phase_split_vector(cell, reaction_rate_multiplier=reaction_rate_multiplier)
    if split.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(split**2)))

def preproject_bed_major_gas_phase_split_to_exchange_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    cell_index: int = 0,
    improvement_factor: float = 0.75,
    preserve_primary_inlet_totals: bool | None = None,
    apply_bc_fn: Callable[[], None] | None = None,
    apply_local_bc_fn: Callable[[int], None] | None = None,
) -> bool:
    """Preproject one bed cell's major gas holdup against dense/bubble phase-split closure.

    Ref: Hamel (1999) two-phase gas balances; thesis NR split merit gate.
    """
    if not cells or cell_index < 0 or cell_index >= len(cells):
        return False
    cell = cells[cell_index]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    gas_rm = _preproj_phase_split_reaction_rate(cfg)

    if preserve_primary_inlet_totals is None:
        preserve_primary_inlet_totals = int(cell_index) == 0

    ref_gas = max(
        float(
            np.sum(
                np.maximum(
                    cell.N_zu_d
                    + cell.N_zu_b
                    + cell.N_d_in
                    + cell.N_b_in
                    + cell.N_rez_d
                    + cell.N_rez_b,
                    0.0,
                )
            )
        ),
        float(np.sum(np.maximum(cell.N_d + cell.N_b, 0.0))),
        1.0,
    )

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()

    def _invalidate_gas_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _scaled_split_vector() -> np.ndarray:
        return _bed_major_gas_phase_split_vector(cell, reaction_rate_multiplier=gas_rm) / ref_gas

    def _scaled_split_rms() -> float:
        split = _scaled_split_vector()
        if split.size == 0:
            return 0.0
        return float(np.sqrt(np.mean(split**2)))

    if preserve_primary_inlet_totals:
        primary_names = ("O2", "H2O", "N2")
        product_names = tuple(sp for sp in MAJOR_NR_GAS_SPECIES if sp not in primary_names)
        primary_idxs = [GAS_SPECIES_INDEX[sp] for sp in primary_names]
        product_idxs = [GAS_SPECIES_INDEX[sp] for sp in product_names]
        primary_totals = np.array(
            [float(max(cell.N_d[j] + cell.N_b[j], 0.0)) for j in primary_idxs],
            dtype=np.float64,
        )
        if np.any(~np.isfinite(primary_totals)) or np.any(primary_totals <= 1e-10):
            return False

        def _set_values(values: np.ndarray) -> None:
            x = np.asarray(values, dtype=np.float64)
            k = 0
            for value, total, j in zip(x[k : k + len(primary_idxs)], primary_totals, primary_idxs):
                dense = float(np.clip(value, 1e-12, total - 1e-12))
                cell.N_d[j] = dense
                cell.N_b[j] = float(total - dense)
            k += len(primary_idxs)
            for value, j in zip(x[k : k + len(product_idxs)], product_idxs):
                cell.N_d[j] = max(float(value), 1e-12)
            k += len(product_idxs)
            for value, j in zip(x[k : k + len(product_idxs)], product_idxs):
                cell.N_b[j] = max(float(value), 1e-12)
            _invalidate_gas_state_cache()

        x0 = np.concatenate(
            (
                np.array([float(cell.N_d[j]) for j in primary_idxs], dtype=np.float64),
                np.array([max(float(cell.N_d[j]), 1e-9) for j in product_idxs], dtype=np.float64),
                np.array([max(float(cell.N_b[j]), 1e-9) for j in product_idxs], dtype=np.float64),
            )
        )
        lower = np.full_like(x0, 1e-12)
        upper = np.concatenate(
            (
                primary_totals - 1e-12,
                np.full(len(product_idxs), 80.0, dtype=np.float64),
                np.full(len(product_idxs), 80.0, dtype=np.float64),
            )
        )
    else:
        major_idxs = [GAS_SPECIES_INDEX[sp] for sp in MAJOR_NR_GAS_SPECIES]
        species_totals = np.array(
            [float(max(cell.N_d[j] + cell.N_b[j], 0.0)) for j in major_idxs],
            dtype=np.float64,
        )
        if np.any(~np.isfinite(species_totals)):
            return False
        active = species_totals > 1.0e-10
        if not np.any(active):
            return False
        active_idxs = [major_idxs[i] for i in range(len(major_idxs)) if active[i]]
        active_totals = species_totals[active]

        def _set_values(values: np.ndarray) -> None:
            x = np.asarray(values, dtype=np.float64).reshape(-1)
            for value, total, j in zip(x, active_totals, active_idxs):
                dense = float(np.clip(value, 1e-12, total - 1e-12))
                cell.N_d[j] = dense
                cell.N_b[j] = float(total - dense)
            _invalidate_gas_state_cache()

        x0 = np.array([float(cell.N_d[j]) for j in active_idxs], dtype=np.float64)
        lower = np.full_like(x0, 1e-12)
        upper = active_totals - 1e-12

    bc_ctx = PreprojectionBcContext(
        cell_index=int(cell_index),
        apply_bc_fn=apply_bc_fn,
        apply_local_bc_fn=apply_local_bc_fn,
    )
    if bc_ctx.needs_bc():
        bc_ctx.finalize()
    initial_split_rms = _scaled_split_rms()
    if not np.isfinite(initial_split_rms) or initial_split_rms <= 1.0e-12:
        return False

    def _restore() -> None:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()

    def _phase_split_objective(values: np.ndarray) -> np.ndarray:
        _set_values(values)
        bc_ctx.refresh_trial()
        return _scaled_split_vector()

    def _record_accept_diag() -> None:
        cell._bottom_major_phase_split_preprojection_diag = {
            "cell_index": int(cell_index),
            "initial_split_rms_scaled": float(initial_split_rms),
            "final_split_rms_scaled": float(_scaled_split_rms()),
        }

    return run_improvement_gated_lsq(
        spec=LsqPreprojectionSpec(
            x0=x0,
            lower=lower,
            upper=upper,
            max_nfev=32,
            xtol=1.0e-10,
            ftol=1.0e-10,
            gtol=1.0e-10,
        ),
        set_values=_set_values,
        objective=_phase_split_objective,
        initial_metric=initial_split_rms,
        final_metric=_scaled_split_rms,
        improvement_factor=float(improvement_factor),
        restore=_restore,
        sync_full_residuals=cell.residuals,
        bc=bc_ctx,
        finalize_on_accept=False,
        on_accept=_record_accept_diag,
    )

def preproject_bottom_major_gas_phase_split_to_exchange_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    improvement_factor: float = 0.75,
    apply_bc_fn: Callable[[], None] | None = None,
) -> bool:
    """Preproject bed0 major gas holdup against dense/bubble phase-split closure."""
    return preproject_bed_major_gas_phase_split_to_exchange_closure(
        cells,
        cfg,
        cell_index=0,
        improvement_factor=improvement_factor,
        preserve_primary_inlet_totals=True,
        apply_bc_fn=apply_bc_fn,
    )

def _per_bed_nr_phase_split_rms_map(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
) -> dict[int, float]:
    """Per-bed phase-split RMS aligned with thesis global NR merit."""
    from src.solvers.global_nr_solver import build_equation_scales, global_residual, pack_reactor
    from src.solvers.nr_indexing import gas_mode_for_nr, n_gas_var, n_var
    from src.solvers.nr_indexing import main_nr_gas_phase_split_positions_for_cell

    apply_bc_fn()
    ref_gas = max(float(getattr(cfg, "O2_feed", 0.0) + getattr(cfg, "H2O_feed", 0.0) + getattr(cfg, "N2_feed", 0.0)), 1.0)
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    F = global_residual(pack_reactor(cells), cells, apply_bc_fn)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)

    out: dict[int, float] = {}
    offset = 0
    for i, cell in enumerate(cells):
        nv = n_var(cell)
        gas_mode = gas_mode_for_nr(cell)
        ng = n_gas_var(cell)
        if gas_mode == "two_phase" and ng > 0:
            half = ng // 2
            dense = F_hat[offset : offset + half]
            bubble = F_hat[offset + half : offset + ng]
            split_pos = main_nr_gas_phase_split_positions_for_cell(cell)
            if split_pos.size:
                terms = dense[split_pos] - bubble[split_pos]
                out[int(i)] = float(np.sqrt(np.mean(terms**2)))
        offset += nv
    return out


def _select_outer_preinner_phase_split_bed_indices(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
) -> list[int]:
    max_bed_index = int(max(getattr(cfg, "nr_outer_preinner_phase_split_max_bed_index_thesis", 0), 0))
    include_bed1 = bool(getattr(cfg, "nr_outer_preinner_bed1_phase_split_thesis", False))
    tau_upper = float(getattr(cfg, "nr_outer_preinner_phase_split_upper_bed_min_rms_thesis", 0.06))
    candidates = _outer_preinner_phase_split_bed_indices(
        cells,
        max_bed_index=max_bed_index,
        include_bed1=include_bed1,
    )
    split_by_bed = _per_bed_nr_phase_split_rms_map(cells, cfg, apply_bc_fn=apply_bc_fn)
    bed0_split = float(split_by_bed.get(0, 0.0))
    global_split = float(max(split_by_bed.values(), default=0.0)) if split_by_bed else 0.0
    # 平台期：放宽上层阈值，让 bed3+ 进入 T2（oneshot 0.037→0.022 需多格）
    plateau = bool(global_split >= 0.03 or bed0_split >= 0.05)
    tau_eff = min(float(tau_upper), 0.015) if plateau else float(tau_upper)
    selected: list[int] = []
    for bed_index in candidates:
        local_split = float(split_by_bed.get(bed_index, 0.0))
        if bed_index == 0:
            selected.append(0)
        elif bed_index == 1 and include_bed1 and local_split >= min(tau_eff, 0.02):
            # 显式 bed1 旗标：只要局部 split 非冷即试（不再要求严热于 bed0）
            selected.append(1)
        elif plateau and local_split >= tau_eff:
            selected.append(bed_index)
        elif local_split >= tau_eff and local_split >= 0.95 * max(bed0_split, tau_eff):
            # 非平台：达到阈值且与 bed0 同量级
            selected.append(bed_index)
    return selected


def _outer_preinner_phase_split_bed_indices(
    cells: list[Cell],
    *,
    max_bed_index: int,
    include_bed1: bool,
) -> list[int]:
    """Bed indices for outer T2 split LSQ (bed1 optional, bed2+ via max_bed_index).

    ``include_bed1=True`` must include bed1 even when ``max_bed_index=0``
    (Phase2 / VTT contract: flag alone enables bed1 preprojection).
    """
    indices: list[int] = []
    for bed_index in range(int(max_bed_index) + 1):
        if bed_index >= len(cells):
            break
        if bed_index == 1 and not include_bed1:
            continue
        if str(getattr(cells[bed_index], "cell_type", "bed")) != "bed":
            continue
        indices.append(bed_index)
    if include_bed1 and 1 not in indices and len(cells) > 1:
        if str(getattr(cells[1], "cell_type", "bed")) == "bed":
            indices.append(1)
            indices.sort()
    return indices


def _bed1_hotspot_followup_should_run(
    split_by_bed: dict[int, float],
    cfg: Any,
    *,
    already_selected: bool,
    global_split_rms: float | None = None,
) -> bool:
    """True when bed1 is the post-bed0 late-floor hotspot (not mid-solve).

    Ref: thesis NR split merit; late-floor bed1 H₂ after bed0 CO₂ T2.
    Mid-solve bed1 T2 raises energy (8-outer A/B); require bed0 cold + low global S.
    """
    if already_selected:
        return False
    if not bool(getattr(cfg, "nr_outer_preinner_bed1_hotspot_followup_thesis", True)):
        return False
    tau = float(getattr(cfg, "nr_outer_preinner_bed1_hotspot_min_rms_thesis", 0.05))
    tau_bed0 = float(getattr(cfg, "nr_outer_preinner_bed1_hotspot_max_bed0_rms_thesis", 0.01))
    tau_global = float(getattr(cfg, "nr_outer_preinner_bed1_hotspot_max_global_rms_thesis", 0.04))
    bed1_split = float(split_by_bed.get(1, 0.0))
    bed0_split = float(split_by_bed.get(0, 0.0))
    if not (np.isfinite(bed1_split) and np.isfinite(bed0_split)):
        return False
    if bed0_split > tau_bed0:
        return False
    if global_split_rms is not None and float(global_split_rms) > tau_global:
        return False
    return bool(bed1_split >= tau and bed1_split > bed0_split)


def oneshot_bed_dense_species_sum_absorb(
    cell: Cell,
    *,
    species: tuple[str, ...] = ("CO", "H2"),
    fraction: float = 1.0,
    mode: str = "dense",
    bidirectional: bool = False,
) -> dict[str, float]:
    """把物种总量残差吸入/放出 holdup（诊断/Pre-Inner 桥）。

    ``res_sum = in_d+in_b+R_d+R_b-(N_d+N_b)``。
    默认仅当 ``res_sum>0`` 时吸入 ``fraction * res_sum``（历史 CO 桥）。
    ``bidirectional=True``：``res_sum<0`` 时按相份额放出 holdup（bed3 N2 过填；§5.25）。

    - ``mode="dense"``：全部写入/放出密相（会抬 split；handoff §5.5）
    - ``mode="phase_share"``：按当前 ``N_d/(N_d+N_b)`` 分给两相（末态证据可到
      y_CO≈0.145 且不抬 split；handoff §5.11）

    Ref: Hamel (1999) Eq.2-4/2-5 总量形式；handoff §5.4 / §5.11 / §5.25.
    """
    frac = float(np.clip(fraction, 0.0, 1.0))
    absorb_mode = str(mode or "dense").strip().lower()
    added: dict[str, float] = {}
    for sp in species:
        if sp not in GAS_SPECIES_INDEX:
            continue
        j = int(GAS_SPECIES_INDEX[sp])
        in_d = float(cell.N_zu_d[j] + cell.N_rez_d[j] + cell.N_d_in[j])
        in_b = float(cell.N_zu_b[j] + cell.N_rez_b[j] + cell.N_b_in[j])
        R_d = float(cell.R_gas_d[j])
        R_b = float(cell.R_gas_b[j])
        N_d = float(cell.N_d[j])
        N_b = float(cell.N_b[j])
        res_sum = in_d + in_b + R_d + R_b - N_d - N_b
        if not np.isfinite(res_sum) or frac <= 0.0:
            added[sp] = 0.0
            continue
        if res_sum > 1.0e-12:
            delta = float(frac * res_sum)
        elif bidirectional and res_sum < -1.0e-12:
            delta = float(frac * res_sum)  # 负：放出
        else:
            added[sp] = 0.0
            continue
        if absorb_mode == "phase_share":
            tot = max(N_d + N_b, 1e-12)
            f_d = float(N_d / tot)
            cell.N_d[j] = max(N_d + f_d * delta, 1e-12)
            cell.N_b[j] = max(N_b + (1.0 - f_d) * delta, 1e-12)
        else:
            # dense：放出时优先从密相扣，不足再动气泡相
            if delta >= 0.0:
                cell.N_d[j] = max(N_d + delta, 1e-12)
            else:
                take = min(N_d - 1e-12, -delta)
                cell.N_d[j] = max(N_d - take, 1e-12)
                rem = -delta - take
                if rem > 0.0:
                    cell.N_b[j] = max(N_b - rem, 1e-12)
        added[sp] = delta
    if any(abs(float(v)) > 0.0 for v in added.values()):
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False
    return added


def run_outer_preinner_upper_dense_syngas_absorb(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
) -> dict[str, bool | float | list[int]]:
    """上段床密相合成气总量吸收（opt-in）。

    非 bed0 的 phase-split LSQ 冻住 ``N_d+N_b``，无法关 ``res_sum``；本桥直接把
    正总量残差写入 ``N_d``。默认关。

    Ref: handoff §5.4; Hamel Eq.2-4/2-5.
    """
    if not bool(getattr(cfg, "nr_outer_preinner_upper_dense_syngas_absorb_thesis", False)):
        return {"enabled": False}
    min_bed = int(max(getattr(cfg, "nr_outer_preinner_upper_dense_absorb_min_bed_index_thesis", 6), 0))
    raw_sp = getattr(cfg, "nr_outer_preinner_upper_dense_absorb_species_thesis", ("CO", "H2"))
    species = tuple(str(s) for s in raw_sp)
    fraction = float(getattr(cfg, "nr_outer_preinner_upper_dense_absorb_fraction_thesis", 1.0))
    absorb_mode = str(
        getattr(cfg, "nr_outer_preinner_upper_dense_absorb_mode_thesis", "dense") or "dense"
    ).strip().lower()
    apply_bc_fn()
    metrics_before = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    touched: list[int] = []
    added_CO = 0.0
    added_H2 = 0.0
    snaps: list[tuple[int, np.ndarray, np.ndarray]] = []
    for i, cell in enumerate(cells):
        if i < min_bed:
            continue
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        cell.residuals()
        snap_d = np.asarray(cell.N_d, dtype=np.float64).copy()
        snap_b = np.asarray(cell.N_b, dtype=np.float64).copy()
        added = oneshot_bed_dense_species_sum_absorb(
            cell, species=species, fraction=fraction, mode=absorb_mode
        )
        if any(float(v) > 0.0 for v in added.values()):
            snaps.append((i, snap_d, snap_b))
            touched.append(i)
            added_CO += float(added.get("CO", 0.0))
            added_H2 += float(added.get("H2", 0.0))
    if not touched:
        return {
            "enabled": True,
            "beds_touched": [],
            "added_CO_mol_s": 0.0,
            "added_H2_mol_s": 0.0,
            "accepted": False,
            "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
            "split_rms_after": float(metrics_before["gas_phase_split_rms"]),
            "merit_before": float(metrics_before["merit_scaled"]),
            "merit_after": float(metrics_before["merit_scaled"]),
        }
    apply_bc_fn()
    metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    skip_local = bool(
        getattr(cfg, "nr_outer_preinner_upper_dense_absorb_skip_local_gate_thesis", False)
    )
    # 局部吸收不得抬高全局 split / merit（与 bed2+ LSQ 全局守卫同口径）
    split_ok = float(metrics_after["gas_phase_split_rms"]) <= float(
        metrics_before["gas_phase_split_rms"]
    ) * 1.001
    merit_ok = float(metrics_after["merit_scaled"]) <= float(metrics_before["merit_scaled"]) * 1.001
    rms_ok = float(metrics_after["rms_scaled"]) <= float(metrics_before["rms_scaled"]) * 1.05
    # phase_share 末态证据：split/rms 稳，但 energy 主导的 merit 可能微升（§5.11）
    if absorb_mode == "phase_share" and not skip_local:
        accepted = bool(split_ok and rms_ok)
    else:
        accepted = bool(skip_local or (split_ok and merit_ok))
    if not accepted:
        for i, snap_d, snap_b in snaps:
            cells[i].N_d[:] = snap_d
            cells[i].N_b[:] = snap_b
            cells[i]._thermo_cache_valid = False
            cells[i]._hydro_cache_valid = False
        apply_bc_fn()
        metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    return {
        "enabled": True,
        "beds_touched": [int(i) for i in touched],
        "added_CO_mol_s": float(added_CO),
        "added_H2_mol_s": float(added_H2),
        "accepted": accepted,
        "absorb_mode": absorb_mode,
        "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
        "split_rms_after": float(metrics_after["gas_phase_split_rms"]),
        "merit_before": float(metrics_before["merit_scaled"]),
        "merit_after": float(metrics_after["merit_scaled"]),
    }


def run_outer_preinner_upper_co_local2x2(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
) -> dict[str, bool | float | list[int]]:
    """上段床 CO 局部 2×2 Newton（opt-in）：同时动 ``N_d,CO`` 与 ``N_b,CO``。

    全系统 Newton 对 dense-CO 只给出 ~4% 理想吸收步；同点局部 2×2 约 50% 且
    ``frac=1`` oneshot 可降 split（handoff §5.7）。默认关。

    Ref: Hamel (1999) Eq.2-4/2-5；局部缩并 NR。
    """
    if not bool(getattr(cfg, "nr_outer_preinner_upper_co_local2x2_thesis", False)):
        return {"enabled": False}

    from src.solvers.global_nr_solver import (
        build_equation_scales,
        global_residual,
        pack_reactor,
        unpack_reactor,
    )
    from src.solvers.nr_indexing import (
        gas_mode_for_nr,
        install_main_nr_gas_indexing_on_cells,
        main_nr_gas_idx_for_cell,
        n_gas_var,
        n_var,
    )
    from src.core.species import GAS_SPECIES

    min_bed = int(max(getattr(cfg, "nr_outer_preinner_upper_co_local2x2_min_bed_index_thesis", 6), 0))
    frac = float(np.clip(getattr(cfg, "nr_outer_preinner_upper_co_local2x2_fraction_thesis", 1.0), 0.0, 1.0))
    if frac <= 0.0:
        return {"enabled": True, "accepted": False, "beds_touched": [], "fraction": frac}

    install_main_nr_gas_indexing_on_cells(cells, cfg)
    apply_bc_fn()
    metrics_before = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    ref_gas = max(
        float(getattr(cfg, "O2_feed", 0.0) + getattr(cfg, "H2O_feed", 0.0) + getattr(cfg, "N2_feed", 0.0)),
        1.0,
    )
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    x0 = pack_reactor(cells)
    F0 = global_residual(x0, cells, apply_bc_fn)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F0, dtype=np.float64) / np.asarray(scale, dtype=np.float64)

    # DOF map: bed -> (dense_dof, bubble_dof, dense_res, bubble_res)
    dof_map: dict[int, tuple[int, int, int, int]] = {}
    offset = 0
    for i, cell in enumerate(cells):
        nv = int(n_var(cell))
        ng = int(n_gas_var(cell))
        if gas_mode_for_nr(cell) == "two_phase" and ng > 0 and i >= min_bed:
            if str(getattr(cell, "cell_type", "bed")) == "bed":
                gas_idx = main_nr_gas_idx_for_cell(cell)
                half = ng // 2
                for p, gi in enumerate(gas_idx):
                    if str(GAS_SPECIES[int(gi)]) == "CO":
                        dof_map[i] = (offset + p, offset + half + p, offset + p, offset + half + p)
                        break
        offset += nv

    if not dof_map:
        return {
            "enabled": True,
            "accepted": False,
            "beds_touched": [],
            "fraction": frac,
            "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
            "split_rms_after": float(metrics_before["gas_phase_split_rms"]),
            "merit_before": float(metrics_before["merit_scaled"]),
            "merit_after": float(metrics_before["merit_scaled"]),
        }

    snaps: list[tuple[int, np.ndarray, np.ndarray]] = []
    touched: list[int] = []
    steps: dict[str, list[float]] = {}
    eps = 1.0e-4
    j_co = int(GAS_SPECIES_INDEX["CO"])
    # FD columns then apply — restore x0 between column probes; hold snaps of pre-edit state.
    for bed_index, (md, mb, rd, rb) in sorted(dof_map.items()):
        A = np.zeros((2, 2), dtype=np.float64)
        for col, dof in enumerate((md, mb)):
            xp = np.asarray(x0, dtype=np.float64).copy()
            xp[int(dof)] += eps
            unpack_reactor(xp, cells)
            apply_bc_fn()
            Fp = global_residual(pack_reactor(cells), cells, apply_bc_fn)
            Fph = np.asarray(Fp, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
            A[0, col] = (float(Fph[rd]) - float(F_hat[rd])) / eps
            A[1, col] = (float(Fph[rb]) - float(F_hat[rb])) / eps
        unpack_reactor(x0, cells)
        apply_bc_fn()
        f2 = np.array([float(F_hat[rd]), float(F_hat[rb])], dtype=np.float64)
        try:
            dx2 = np.linalg.solve(A, -f2)
        except np.linalg.LinAlgError:
            dx2, *_ = np.linalg.lstsq(A, -f2, rcond=None)
        if not np.all(np.isfinite(dx2)):
            continue
        steps[str(bed_index)] = [float(dx2[0]), float(dx2[1])]

    unpack_reactor(x0, cells)
    apply_bc_fn()
    for bed_index, (md, mb, _rd, _rb) in sorted(dof_map.items()):
        if str(bed_index) not in steps:
            continue
        cell = cells[bed_index]
        snap_d = np.asarray(cell.N_d, dtype=np.float64).copy()
        snap_b = np.asarray(cell.N_b, dtype=np.float64).copy()
        dnd, dnb = steps[str(bed_index)]
        cell.N_d[j_co] = max(float(x0[md]) + frac * dnd, 1e-12)
        cell.N_b[j_co] = max(float(x0[mb]) + frac * dnb, 1e-12)
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False
        snaps.append((bed_index, snap_d, snap_b))
        touched.append(bed_index)

    if not touched:
        return {
            "enabled": True,
            "accepted": False,
            "beds_touched": [],
            "fraction": frac,
            "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
            "split_rms_after": float(metrics_before["gas_phase_split_rms"]),
            "merit_before": float(metrics_before["merit_scaled"]),
            "merit_after": float(metrics_before["merit_scaled"]),
        }

    apply_bc_fn()
    metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    split_ok = float(metrics_after["gas_phase_split_rms"]) <= float(
        metrics_before["gas_phase_split_rms"]
    ) * 1.001
    merit_ok = float(metrics_after["merit_scaled"]) <= float(metrics_before["merit_scaled"]) * 1.001
    # rms 单独守卫：merit 可能被 energy 主导而放过 rms 爆炸
    rms_ok = float(metrics_after["rms_scaled"]) <= float(metrics_before["rms_scaled"]) * 1.001
    accepted = bool(split_ok and merit_ok and rms_ok)
    if not accepted:
        for bed_index, snap_d, snap_b in snaps:
            cells[bed_index].N_d[:] = snap_d
            cells[bed_index].N_b[:] = snap_b
            cells[bed_index]._thermo_cache_valid = False
            cells[bed_index]._hydro_cache_valid = False
        apply_bc_fn()
        metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)

    return {
        "enabled": True,
        "accepted": accepted,
        "beds_touched": [int(i) for i in touched],
        "fraction": frac,
        "dx2_by_bed": steps,
        "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
        "split_rms_after": float(metrics_after["gas_phase_split_rms"]),
        "merit_before": float(metrics_before["merit_scaled"]),
        "merit_after": float(metrics_after["merit_scaled"]),
        "rms_before": float(metrics_before["rms_scaled"]),
        "rms_after": float(metrics_after["rms_scaled"]),
        "rms_ok": bool(rms_ok),
        "split_ok": bool(split_ok),
        "merit_ok": bool(merit_ok),
    }


def run_outer_preinner_gas_phase_split_preprojection(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    apply_local_bc_fn: Callable[[int], None] | None = None,
    split_shock: bool = False,
) -> dict[str, bool | float]:
    """Major phase-split bridge after Vorabrechnung refresh, before inner global NR.

    Ref: Hamel outer Abgleich Startwertwahl; thesis NR split merit gate.
    """
    if not bool(getattr(cfg, "nr_outer_preinner_phase_split_preprojection_thesis", True)):
        return {"enabled": False}
    max_bed_index = int(max(getattr(cfg, "nr_outer_preinner_phase_split_max_bed_index_thesis", 0), 0))
    include_bed1 = bool(getattr(cfg, "nr_outer_preinner_bed1_phase_split_thesis", False))
    bed_indices = _select_outer_preinner_phase_split_bed_indices(
        cells,
        cfg,
        apply_bc_fn=apply_bc_fn,
    )

    metrics_before = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    bed0_factor = float(
        getattr(cfg, "nr_outer_preinner_split_shock_bed0_improvement_thesis", 0.65)
        if split_shock
        else getattr(cfg, "nr_outer_preinner_split_plateau_bed0_improvement_thesis", 0.85)
    )
    accepted: dict[str, bool] = {}
    global_split_rollback = 0
    global_residual_rollback = 0
    running_split = float(metrics_before["gas_phase_split_rms"])
    running_rms = float(metrics_before["rms_scaled"])
    running_max_abs = float(metrics_before["max_abs_scaled"])
    # 局部 split 向量可降，但全系统 NR 残差爆炸（bed0 bubble/H2O，handoff §5.9）
    residual_guard = bool(
        getattr(cfg, "nr_outer_preinner_split_global_residual_guard_thesis", True)
    )
    rms_factor = float(
        getattr(cfg, "nr_outer_preinner_split_global_rms_factor_thesis", 1.25)
    )
    max_abs_factor = float(
        getattr(cfg, "nr_outer_preinner_split_global_max_abs_factor_thesis", 1.25)
    )
    max_abs_delta = float(
        getattr(cfg, "nr_outer_preinner_split_global_max_abs_delta_thesis", 0.5)
    )
    for bed_index in bed_indices:
        if bed_index >= len(cells):
            continue
        cell = cells[bed_index]
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        # 上层床 LSQ：0.95 过严导致 stall 上 bed2+ 拒收；0.99 与 oneshot 证据对齐
        bed_factor = bed0_factor if bed_index == 0 else 0.99
        snap_d = np.asarray(cell.N_d, dtype=np.float64).copy()
        snap_b = np.asarray(cell.N_b, dtype=np.float64).copy()
        # bed0 旧模式：冻 O2/H2O/N2 总量、产物 Nd/Nb 自由≤80 → 可凭空造 CO 并引爆 H2O 残差（§5.10）
        bed0_phase_only = bool(
            getattr(cfg, "nr_outer_preinner_bed0_phase_only_majors_thesis", True)
        )
        preserve_primary = bool(bed_index == 0) and not bed0_phase_only
        ok = bool(
            preproject_bed_major_gas_phase_split_to_exchange_closure(
                cells,
                cfg,
                cell_index=bed_index,
                improvement_factor=bed_factor,
                preserve_primary_inlet_totals=preserve_primary,
                apply_bc_fn=apply_bc_fn,
                apply_local_bc_fn=apply_local_bc_fn,
            )
        )
        if ok:
            apply_bc_fn()
            m_now = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
            split_now = float(m_now["gas_phase_split_rms"])
            rms_now = float(m_now["rms_scaled"])
            max_now = float(m_now["max_abs_scaled"])
            split_worse = split_now > running_split * 1.001
            residual_blow = bool(residual_guard) and (
                rms_now > running_rms * max(rms_factor, 1.0)
                or max_now > max(running_max_abs * max(max_abs_factor, 1.0), running_max_abs + max_abs_delta)
            )
            # 局部 LSQ 可接受但抬高全局 split / 引爆全系统残差 → 回滚该格
            if split_worse or residual_blow:
                cell.N_d[:] = snap_d
                cell.N_b[:] = snap_b
                cell._thermo_cache_valid = False
                cell._hydro_cache_valid = False
                apply_bc_fn()
                ok = False
                global_split_rollback += 1
                if residual_blow and not split_worse:
                    global_residual_rollback += 1
            else:
                running_split = split_now
                running_rms = rms_now
                running_max_abs = max_now
        accepted[f"bed{bed_index}_accepted"] = ok

    # After bed0 T2, bed1 may become the sole late-floor hotspot.
    bed1_followup = False
    if len(cells) > 1 and str(getattr(cells[1], "cell_type", "bed")) == "bed":
        apply_bc_fn()
        metrics_mid = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
        split_after = _per_bed_nr_phase_split_rms_map(cells, cfg, apply_bc_fn=apply_bc_fn)
        if _bed1_hotspot_followup_should_run(
            split_after,
            cfg,
            already_selected=(1 in bed_indices),
            global_split_rms=float(metrics_mid["gas_phase_split_rms"]),
        ):
            bed1_followup = True
            cell1 = cells[1]
            snap1_d = np.asarray(cell1.N_d, dtype=np.float64).copy()
            snap1_b = np.asarray(cell1.N_b, dtype=np.float64).copy()
            ok1 = bool(
                preproject_bed_major_gas_phase_split_to_exchange_closure(
                    cells,
                    cfg,
                    cell_index=1,
                    improvement_factor=0.95,
                    preserve_primary_inlet_totals=False,
                    apply_bc_fn=apply_bc_fn,
                    apply_local_bc_fn=apply_local_bc_fn,
                )
            )
            if ok1 and residual_guard:
                apply_bc_fn()
                m1 = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
                if (
                    float(m1["rms_scaled"]) > running_rms * max(rms_factor, 1.0)
                    or float(m1["max_abs_scaled"])
                    > max(
                        running_max_abs * max(max_abs_factor, 1.0),
                        running_max_abs + max_abs_delta,
                    )
                ):
                    cell1.N_d[:] = snap1_d
                    cell1.N_b[:] = snap1_b
                    cell1._thermo_cache_valid = False
                    cell1._hydro_cache_valid = False
                    apply_bc_fn()
                    ok1 = False
                    global_split_rollback += 1
                    global_residual_rollback += 1
            accepted["bed1_accepted"] = ok1

    apply_bc_fn()
    metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    return {
        "enabled": True,
        "max_bed_index": float(max_bed_index),
        "include_bed1": include_bed1,
        "bed_indices": [int(i) for i in bed_indices],
        "bed1_hotspot_followup": bool(bed1_followup),
        "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
        "split_rms_after": float(metrics_after["gas_phase_split_rms"]),
        "merit_before": float(metrics_before["merit_scaled"]),
        "merit_after": float(metrics_after["merit_scaled"]),
        "split_shock": bool(split_shock),
        "bed0_improvement_factor": float(bed0_factor),
        "bed0_accepted": bool(accepted.get("bed0_accepted", False)),
        "bed1_accepted": bool(accepted.get("bed1_accepted", False)),
        "global_split_rollbacks": float(global_split_rollback),
        "global_residual_rollbacks": float(global_residual_rollback),
        "rms_scaled_before": float(metrics_before["rms_scaled"]),
        "rms_scaled_after": float(metrics_after["rms_scaled"]),
        "max_abs_scaled_before": float(metrics_before["max_abs_scaled"]),
        "max_abs_scaled_after": float(metrics_after["max_abs_scaled"]),
    }

def preproject_bed_gas_temperature_to_energy_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    cell_index: int = 0,
    gas_relax: float = 1.0,
    temperature_bounds_K: tuple[float, float] = (650.0, 1600.0),
    improvement_factor: float = 0.75,
    couple_neighbor_cell_index: int | None = None,
    apply_bc_fn: Callable[[], None] | None = None,
    apply_local_bc_fn: Callable[[int], None] | None = None,
    neighbor_energy_weight: float = 2.0,
) -> bool:
    """Bed-cell x0 projection for total gas closure and energy accessibility.

    This is an initialization-only bridge for stiff bottom / neighbor beds.
    It does not modify reaction sources and does not put drying/pyrolysis gas
    into the bubble source.  Instead:

    - major gas dense/bubble outlet values are solved with bounded least squares
      against total gas closure, phase residuals, and energy closure;
    - TAR/NH3 source support is added only to dense outlet variables because the
      direct pyrolysis source remains a suspension-phase source;
    - temperature is bounded to keep the projection in a plausible x0 range.
    - optional ``couple_neighbor_cell_index`` adds the downstream bed temperature
      as an unknown and refreshes BC before evaluating the neighbor energy row.

    Ref: status_report.md bottom-zone Startwertwahl bridge; Hamel Vorabrechnung.
    """
    if not cells:
        return False
    if cell_index < 0 or cell_index >= len(cells):
        return False
    cell = cells[cell_index]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    gas_rm = _preproj_gas_reaction_rate(cfg)

    neighbor: Cell | None = None
    couple_neighbor = couple_neighbor_cell_index is not None
    if couple_neighbor:
        if apply_bc_fn is None:
            return False
        n_idx = int(couple_neighbor_cell_index)
        if n_idx < 0 or n_idx >= len(cells):
            return False
        neighbor = cells[n_idx]
        if str(getattr(neighbor, "cell_type", "bed")) != "bed":
            return False

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()
    original_T = float(cell.T)
    original_neighbor_T: float | None = None
    if neighbor is not None:
        original_neighbor_T = float(neighbor.T)

    def _invalidate_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _gas_max(res: np.ndarray) -> float:
        n = cell.N_d.shape[0]
        return float(max(np.max(np.abs(res[:n])), np.max(np.abs(res[n : 2 * n]))))

    def _combined_gas_max(res: np.ndarray) -> float:
        n = cell.N_d.shape[0]
        return float(np.max(np.abs(res[:n] + res[n : 2 * n])))

    n_gas = cell.N_d.shape[0]
    projection_omega = float(np.clip(gas_relax, 0.0, 1.0))
    projection_names = nr_projection_gas_species(cfg)
    dense_idxs = [GAS_SPECIES_INDEX[sp] for sp in projection_names]
    bubble_idxs = [GAS_SPECIES_INDEX[sp] for sp in MAJOR_NR_GAS_SPECIES if sp in projection_names]
    dense_only_source_names = nr_dense_only_pyrolysis_support_species(cfg)
    n_dense = len(dense_idxs)
    n_bubble = len(bubble_idxs)
    # Metrics/gates use the same species subset as the LSQ unknown vector.
    gate_idxs = dense_idxs

    def _projection_gas_max(res: np.ndarray) -> float:
        dense_part = max(abs(float(res[j])) for j in dense_idxs)
        bubble_part = max(abs(float(res[n_gas + j])) for j in bubble_idxs) if bubble_idxs else 0.0
        return float(max(dense_part, bubble_part))

    def _combined_projection_gas_max(res: np.ndarray) -> float:
        return float(max(abs(float(res[j] + res[n_gas + j])) for j in gate_idxs))

    initial_res = _preproj_cell_residual_vector(
        cell,
        "gas",
        "energy",
        reaction_rate_multiplier=gas_rm,
    )
    initial_energy_abs = abs(float(initial_res[-1]))
    initial_gas_max = _projection_gas_max(initial_res)
    initial_combined_gas_max = _combined_projection_gas_max(initial_res)
    if not np.isfinite(initial_energy_abs) or initial_energy_abs <= 1e-9:
        return False

    fixed_dense_support = original_d.copy()
    fixed_bubble_support = original_b.copy()
    for sp in dense_only_source_names:
        j = GAS_SPECIES_INDEX[sp]
        dense_support = max(float(initial_res[j]), 0.0)
        # Seed pyrolysis support at the full residual; ``projection_omega`` only
        # scales the major-gas/T move so partial relax does not desync LSQ closure.
        fixed_dense_support[j] = max(float(fixed_dense_support[j] + dense_support), 0.0)
        # Keep bubble untouched: direct drying/pyrolysis source remains dense-only.
        fixed_bubble_support[j] = max(float(fixed_bubble_support[j]), 0.0)

    def _h_gas(flow: np.ndarray, temperature: float) -> float:
        return calc_gas_enthalpy_flow(np.asarray(flow, dtype=np.float64), float(temperature), cell._h_cache)

    def _h_solid(flow: np.ndarray, temperature: float) -> float:
        return calc_solid_enthalpy_flow(
            np.asarray(flow, dtype=np.float64),
            float(temperature),
            ash_dry_wt=cell.solid.ash_dry_wt,
            VM_daf=cell.solid.VM_daf,
            h_f_dry=cell.solid.h_f_dry,
            cp_char_fn=cp_char,
            cp_ash_fn=cp_ash,
            cp_sand_fn=cp_sand,
        )

    attr_lo = float(getattr(cell, "_nr_temperature_min_K", temperature_bounds_K[0]))
    attr_hi = float(getattr(cell, "_nr_temperature_max_K", temperature_bounds_K[1]))
    if not np.isfinite(attr_lo):
        attr_lo = float(temperature_bounds_K[0])
    if not np.isfinite(attr_hi):
        attr_hi = float(temperature_bounds_K[1])
    t_lo = max(float(temperature_bounds_K[0]), attr_lo, 300.0)
    t_hi = min(float(temperature_bounds_K[1]), attr_hi, 2500.0)
    t_hi = max(t_hi, t_lo + 1.0)

    neighbor_t_lo = t_lo
    neighbor_t_hi = t_hi
    if neighbor is not None:
        n_attr_lo = float(getattr(neighbor, "_nr_temperature_min_K", temperature_bounds_K[0]))
        n_attr_hi = float(getattr(neighbor, "_nr_temperature_max_K", temperature_bounds_K[1]))
        if np.isfinite(n_attr_lo):
            neighbor_t_lo = max(float(temperature_bounds_K[0]), n_attr_lo, 300.0)
        if np.isfinite(n_attr_hi):
            neighbor_t_hi = min(float(temperature_bounds_K[1]), n_attr_hi, 2500.0)
        neighbor_t_hi = max(neighbor_t_hi, neighbor_t_lo + 1.0)

    initial_neighbor_energy_abs = 0.0
    if neighbor is not None:
        _preproj_finalize_bc(apply_bc_fn)
        initial_neighbor_energy_abs = _preproj_cell_energy_residual_abs(neighbor)

    neighbor_indices = ()
    if neighbor is not None and couple_neighbor_cell_index is not None:
        neighbor_indices = (int(cell_index), int(couple_neighbor_cell_index))

    def _refresh_trial_bc() -> None:
        _preproj_refresh_trial_bc(
            cell_index,
            apply_bc_fn=apply_bc_fn,
            apply_local_bc_fn=apply_local_bc_fn,
            neighbor_indices=neighbor_indices,
        )

    def _set_candidate(x: np.ndarray) -> None:
        arr = np.asarray(x, dtype=np.float64)
        cell.N_d[:] = fixed_dense_support
        cell.N_b[:] = fixed_bubble_support
        for value, j in zip(arr[:n_dense], dense_idxs):
            cell.N_d[j] = max(float(value), 1e-12)
        for value, j in zip(arr[n_dense : n_dense + n_bubble], bubble_idxs):
            cell.N_b[j] = max(float(value), 1e-12)
        if couple_neighbor:
            cell.T = float(np.clip(arr[-2], t_lo, t_hi))
            if neighbor is not None:
                neighbor.T = float(np.clip(arr[-1], neighbor_t_lo, neighbor_t_hi))
        else:
            cell.T = float(np.clip(arr[-1], t_lo, t_hi))
        _invalidate_state_cache()
        if neighbor is not None:
            neighbor._thermo_cache_valid = False
            neighbor._hydro_cache_valid = False
        _refresh_trial_bc()

    ref_gas = max(
        float(
            np.sum(
                np.maximum(
                    cell.N_zu_d
                    + cell.N_zu_b
                    + cell.N_d_in
                    + cell.N_b_in
                    + cell.N_rez_d
                    + cell.N_rez_b,
                    0.0,
                )
            )
        ),
        1.0,
    )
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)

    def _objective(x: np.ndarray) -> np.ndarray:
        _set_candidate(x)
        res = _preproj_cell_residual_vector(
            cell,
            "gas",
            "energy",
            reaction_rate_multiplier=gas_rm,
        )
        combined = (res[dense_idxs] + res[[n_gas + j for j in dense_idxs]]) / ref_gas
        phase = np.concatenate(
            (
                res[dense_idxs],
                res[[n_gas + j for j in bubble_idxs]],
            )
        ) / ref_gas
        energy = np.array([float(res[-1]) / ref_energy], dtype=np.float64)
        terms = [combined, 0.4 * phase, 3.0 * energy]
        if neighbor is not None:
            neighbor_energy = np.array(
                [float(neighbor_energy_weight) * _preproj_cell_energy_residual_abs(neighbor) / ref_energy],
                dtype=np.float64,
            )
            terms.append(neighbor_energy)
        # Total closure is the priority; phase rows remain in the objective so
        # the projection cannot "solve" total gas by destroying two-phase balance.
        return np.concatenate(terms)

    x0 = np.concatenate(
        (
            np.array([max(float(original_d[j]), 1e-12) for j in dense_idxs], dtype=np.float64),
            np.array([max(float(original_b[j]), 1e-12) for j in bubble_idxs], dtype=np.float64),
            np.array([float(np.clip(original_T, t_lo, t_hi))], dtype=np.float64),
        )
    )
    if neighbor is not None:
        assert original_neighbor_T is not None
        x0 = np.concatenate(
            (
                x0,
                np.array([float(np.clip(original_neighbor_T, neighbor_t_lo, neighbor_t_hi))], dtype=np.float64),
            )
        )
    lower = np.full_like(x0, 1e-12)
    upper = np.full_like(x0, 120.0)
    if couple_neighbor:
        lower[-2] = t_lo
        upper[-2] = t_hi
        lower[-1] = neighbor_t_lo
        upper[-1] = neighbor_t_hi
    else:
        lower[-1] = t_lo
        upper[-1] = t_hi

    def _restore_gas_t_state() -> None:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        cell.T = original_T
        if neighbor is not None and original_neighbor_T is not None:
            neighbor.T = original_neighbor_T
            neighbor._thermo_cache_valid = False
            neighbor._hydro_cache_valid = False
        _invalidate_state_cache()

    final_res_holder: dict[str, np.ndarray | None] = {"res": None}

    def _apply_lsq_solution(raw_x: np.ndarray) -> None:
        relaxed_candidate = np.asarray(x0 + projection_omega * (raw_x - x0), dtype=np.float64)
        _set_candidate(relaxed_candidate)
        final_res_holder["res"] = _preproj_cell_residual_vector(
            cell,
            "gas",
            "energy",
            reaction_rate_multiplier=gas_rm,
        )

    def _accept_gas_t_projection() -> bool:
        final_res = final_res_holder["res"]
        if final_res is None:
            return False
        final_energy_abs = abs(float(final_res[-1]))
        final_combined_gas_max = _combined_projection_gas_max(final_res)
        final_gas_max = _projection_gas_max(final_res)
        required_energy = initial_energy_abs - projection_omega * (
            initial_energy_abs - improvement_factor * initial_energy_abs
        )
        required_combined = initial_combined_gas_max - projection_omega * (
            initial_combined_gas_max - improvement_factor * max(initial_combined_gas_max, 1e-9)
        )
        energy_improved = np.isfinite(final_energy_abs) and final_energy_abs <= required_energy
        if projection_omega < 1.0 - 1e-12:
            gas_improved = (
                np.isfinite(final_combined_gas_max)
                and final_combined_gas_max <= required_combined
                and np.isfinite(final_gas_max)
            )
        else:
            gas_improved = (
                np.isfinite(final_combined_gas_max)
                and final_combined_gas_max <= required_combined
                and np.isfinite(final_gas_max)
                and final_gas_max <= max(initial_gas_max, 1e-9)
            )
        if not (energy_improved and gas_improved):
            return False
        if neighbor is not None:
            neighbor_energy_abs = _preproj_cell_energy_residual_abs(neighbor)
            guard_ratio = float(max(getattr(cfg, "nr_bottom_joint_bed1_energy_guard_ratio_thesis", 1.25), 1.0))
            max_increase = max((guard_ratio - 1.0) * max(initial_neighbor_energy_abs, 1.0e3), 1.0e6)
            if neighbor_energy_abs > initial_neighbor_energy_abs + max_increase:
                return False
        return True

    def _record_gas_t_accept_diag() -> None:
        final_res = final_res_holder["res"]
        assert final_res is not None
        final_energy_abs = abs(float(final_res[-1]))
        final_combined_gas_max = _combined_projection_gas_max(final_res)
        final_gas_max = _projection_gas_max(final_res)
        neighbor_energy_abs = _preproj_cell_energy_residual_abs(neighbor) if neighbor is not None else None
        diag = {
            "cell_index": int(cell_index),
            "initial_energy_abs_W": float(initial_energy_abs),
            "final_energy_abs_W": float(final_energy_abs),
            "initial_gas_max_mol_s": float(initial_gas_max),
            "final_gas_max_mol_s": float(final_gas_max),
            "initial_combined_gas_max_mol_s": float(initial_combined_gas_max),
            "final_combined_gas_max_mol_s": float(final_combined_gas_max),
            "initial_full_gas_max_mol_s": float(_gas_max(initial_res)),
            "final_full_gas_max_mol_s": float(_gas_max(final_res)),
            "candidate_temperature_K": float(cell.T),
            "bc_aware_neighbor_index": int(couple_neighbor_cell_index) if couple_neighbor else None,
            "initial_neighbor_energy_abs_W": float(initial_neighbor_energy_abs) if neighbor is not None else None,
            "final_neighbor_energy_abs_W": float(neighbor_energy_abs) if neighbor_energy_abs is not None else None,
        }
        if cell_index == 0:
            cell._bottom_joint_preprojection_diag = diag
        else:
            cell._bottom_neighbor_preprojection_diag = diag

    return run_lsq_preprojection(
        spec=LsqPreprojectionSpec(
            x0=x0,
            lower=lower,
            upper=upper,
            max_nfev=40,
            xtol=1.0e-8,
            ftol=1.0e-8,
            gtol=1.0e-8,
        ),
        objective=_objective,
        apply_solution=_apply_lsq_solution,
        restore=_restore_gas_t_state,
        accept_predicate=_accept_gas_t_projection,
        sync_full_residuals=cell.residuals,
        bc=PreprojectionBcContext(
            cell_index=int(cell_index),
            apply_bc_fn=apply_bc_fn,
            apply_local_bc_fn=apply_local_bc_fn,
            neighbor_indices=neighbor_indices,
        ),
        on_accept=_record_gas_t_accept_diag,
    )

def preproject_bed_solid_holdup_to_transport_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    cell_index: int,
    apply_bc_fn: Callable[[], None] | None = None,
    apply_local_bc_fn: Callable[[int], None] | None = None,
    improvement_factor: float = 0.75,
    holdup_upper_scale: float = 3.0,
    skip_rms_scaled: float | None = None,
    max_nfev: int | None = None,
) -> bool:
    """Preproject one bed cell's solid holdup against local Eq. 2.6 transport closure.

    Adjusts ``m_solid`` only; reactions/BC remain live through ``cell.residuals()``.
    Ref: Hamel (1999) Eq. 2.6; thesis holdup_transport bottom-zone Startwertwahl bridge.
    """
    if cell_index < 0 or cell_index >= len(cells):
        return False
    cell = cells[cell_index]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False
    if str(getattr(cell, "solid_state_model", "legacy_stream")) != "holdup_transport":
        return False

    bc_ctx = PreprojectionBcContext(
        cell_index=int(cell_index),
        apply_bc_fn=apply_bc_fn,
        apply_local_bc_fn=apply_local_bc_fn,
    )
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    original = np.asarray(cell.m_solid, dtype=np.float64).copy()
    n_solid = int(original.size)
    if n_solid <= 0:
        return False

    if bc_ctx.needs_bc():
        bc_ctx.finalize()
    holdup_rm = float(
        getattr(cfg, "preproj_solid_holdup_reaction_rate_multiplier_thesis", 0.0)
    )
    initial_solid = _preproj_cell_residual_vector(
        cell,
        "solid",
        reaction_rate_multiplier=holdup_rm,
    )
    initial_rms = scaled_rms(initial_solid, ref_solid)
    if not np.isfinite(initial_rms) or initial_rms <= 1.0e-12:
        return False
    skip_threshold = (
        float(skip_rms_scaled)
        if skip_rms_scaled is not None
        else float(getattr(cfg, "nr_bottom_zone_solid_holdup_skip_rms_scaled_thesis", 0.0))
    )
    if skip_threshold > 0.0 and initial_rms <= skip_threshold:
        cell._bottom_zone_solid_preprojection_diag = {
            "cell_index": int(cell_index),
            "skipped_well_converged": True,
            "initial_solid_rms_kg_s": float(initial_rms * ref_solid),
        }
        return False
    lsq_max_nfev = int(
        max_nfev
        if max_nfev is not None
        else int(getattr(cfg, "nr_bottom_zone_solid_holdup_lsq_max_nfev_thesis", 12))
    )
    lsq_max_nfev = max(lsq_max_nfev, 4)

    def _invalidate() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _restore() -> None:
        cell.m_solid[:] = original
        _invalidate()

    def _set_holdup(values: np.ndarray) -> None:
        cell.m_solid[:] = np.maximum(np.asarray(values, dtype=np.float64).reshape(original.shape), 1.0e-12)
        _invalidate()

    def _solid_objective(values: np.ndarray) -> np.ndarray:
        _set_holdup(values)
        bc_ctx.refresh_trial()
        solid = _preproj_cell_residual_vector(
            cell,
            "solid",
            reaction_rate_multiplier=holdup_rm,
        )
        return (solid / ref_solid).copy()

    def _final_rms() -> float:
        solid = _preproj_cell_residual_vector(
            cell,
            "solid",
            reaction_rate_multiplier=holdup_rm,
        )
        return scaled_rms(solid, ref_solid)

    x0 = np.maximum(original.reshape(-1), 1.0e-12)
    upper_scale = float(max(holdup_upper_scale, 1.05))
    upper = np.maximum(x0 * upper_scale, 1.0e-6)
    lower = np.full_like(x0, 1.0e-12)
    initial_rms_kg_s = float(initial_rms * ref_solid)

    def _record_accept_diag() -> None:
        cell._bottom_zone_solid_preprojection_diag = {
            "cell_index": int(cell_index),
            "initial_solid_rms_kg_s": initial_rms_kg_s,
            "final_solid_rms_kg_s": float(_final_rms() * ref_solid),
        }

    return run_improvement_gated_lsq(
        spec=LsqPreprojectionSpec(
            x0=x0,
            lower=lower,
            upper=upper,
            max_nfev=lsq_max_nfev,
            xtol=1.0e-8,
            ftol=1.0e-8,
            gtol=1.0e-8,
        ),
        set_values=_set_holdup,
        objective=_solid_objective,
        initial_metric=initial_rms,
        final_metric=_final_rms,
        improvement_factor=float(improvement_factor),
        restore=_restore,
        sync_full_residuals=cell.residuals,
        bc=bc_ctx,
        on_accept=_record_accept_diag,
    )

def align_downstream_bed_temperatures_to_upstream(
    cells: list[Cell],
    cfg: Any,
    *,
    upstream_index: int = 0,
    downstream_indices: tuple[int, ...] = (1, 2),
    blend_upstream: float | None = None,
) -> bool:
    """Blend downstream bed temperatures toward an upstream neighbor after bottom gas/T moves."""
    if upstream_index < 0 or upstream_index >= len(cells):
        return False
    upstream = cells[upstream_index]
    if str(getattr(upstream, "cell_type", "bed")) != "bed":
        return False
    if blend_upstream is None:
        blend_upstream = float(getattr(cfg, "nr_bottom_zone_downstream_temperature_blend_thesis", 0.75))
    blend = float(np.clip(blend_upstream, 0.0, 1.0))
    changed = False
    t_up = float(upstream.T)
    for j in downstream_indices:
        if j < 0 or j >= len(cells):
            continue
        cell = cells[j]
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        t_lo = float(getattr(cell, "_nr_temperature_min_K", 300.0))
        t_hi = float(getattr(cell, "_nr_temperature_max_K", 2500.0))
        if not np.isfinite(t_lo):
            t_lo = 300.0
        if not np.isfinite(t_hi):
            t_hi = 2500.0
        t_hi = max(t_hi, t_lo + 1.0)
        new_t = float(np.clip((1.0 - blend) * float(cell.T) + blend * t_up, t_lo, t_hi))
        if abs(new_t - float(cell.T)) > 1.0e-9:
            cell.T = new_t
            cell._thermo_cache_valid = False
            cell._hydro_cache_valid = False
            changed = True
    return changed

def run_bottom_zone_solid_energy_preprojection(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None] | None = None,
    apply_local_bc_fn: Callable[[int], None] | None = None,
    solid_repass: bool | None = None,
    solid_holdup_skip_rms_scaled: float | None = None,
) -> dict[str, bool | float]:
    """Bottom-zone solid holdup + bed1 energy response after guarded joint x0.

    Ref: Hamel (1999) Eq. 2.6/2.7; thesis bottom-zone Startwertwahl bridge.
    Order: optional downstream-T blend → solid LSQ → bed1 fence T projection → solid re-pass.
    """
    max_bed_index = int(max(getattr(cfg, "nr_bottom_zone_preprojection_max_bed_index_thesis", 2), 0))
    align_downstream = bool(getattr(cfg, "nr_bottom_zone_align_downstream_temperature_thesis", False))
    bed1_energy_t = bool(getattr(cfg, "nr_bottom_zone_bed1_energy_temperature_preprojection_thesis", True))
    repass_enabled = (
        bool(solid_repass)
        if solid_repass is not None
        else bool(getattr(cfg, "nr_bottom_zone_solid_repass_thesis", True))
    )
    diag: dict[str, bool | float] = {
        "enabled": True,
        "max_bed_index": float(max_bed_index),
        "solid_preprojections_accepted": 0.0,
        "downstream_temperature_aligned": False,
        "bed1_energy_temperature_preprojection_accepted": False,
        "solid_repass_enabled": bool(repass_enabled),
    }
    if not cells:
        return diag

    if align_downstream:
        diag["downstream_temperature_aligned"] = bool(
            align_downstream_bed_temperatures_to_upstream(
                cells,
                cfg,
                upstream_index=0,
                downstream_indices=tuple(
                    j
                    for j in range(1, min(max_bed_index + 1, len(cells)))
                    if str(getattr(cells[j], "cell_type", "bed")) == "bed"
                ),
            )
        )
        if apply_bc_fn is not None:
            apply_bc_fn()

    solid_ok: list[bool] = []
    for i in range(min(max_bed_index + 1, len(cells))):
        solid_ok.append(
            preproject_bed_solid_holdup_to_transport_closure(
                cells,
                cfg,
                cell_index=i,
                apply_bc_fn=apply_bc_fn,
                apply_local_bc_fn=apply_local_bc_fn,
                skip_rms_scaled=solid_holdup_skip_rms_scaled,
            )
        )
    if apply_bc_fn is not None:
        apply_bc_fn()

    if (
        bed1_energy_t
        and len(cells) >= 2
        and str(getattr(cells[1], "cell_type", "bed")) == "bed"
        and apply_bc_fn is not None
    ):
        diag["bed1_energy_temperature_preprojection_accepted"] = bool(
            preproject_bed_temperature_to_local_energy_closure(
                cells,
                cfg,
                cell_index=1,
                apply_bc_fn=apply_bc_fn,
                apply_local_bc_fn=apply_local_bc_fn,
            )
        )
        if apply_bc_fn is not None:
            apply_bc_fn()

    if repass_enabled:
        for i in range(1, min(max_bed_index + 1, len(cells))):
            if str(getattr(cells[i], "cell_type", "bed")) != "bed":
                continue
            solid_ok.append(
                preproject_bed_solid_holdup_to_transport_closure(
                    cells,
                    cfg,
                    cell_index=i,
                    apply_bc_fn=apply_bc_fn,
                    apply_local_bc_fn=apply_local_bc_fn,
                    skip_rms_scaled=solid_holdup_skip_rms_scaled,
                )
            )
    diag["solid_preprojections_accepted"] = float(sum(1 for ok in solid_ok if ok))
    return diag

def _snapshot_nr_bottom_bed_states(
    cells: list[Cell],
    max_bed_index: int,
) -> list[dict[str, np.ndarray | float]]:
    """Capture bottom-bed NR state for outer pre-inner rollback."""
    snapshots: list[dict[str, np.ndarray | float]] = []
    for i in range(min(max_bed_index + 1, len(cells))):
        cell = cells[i]
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        snapshots.append(
            {
                "cell_index": float(i),
                "N_d": np.asarray(cell.N_d, dtype=np.float64).copy(),
                "N_b": np.asarray(cell.N_b, dtype=np.float64).copy(),
                "T": float(cell.T),
                "m_solid": np.asarray(cell.m_solid, dtype=np.float64).copy(),
            }
        )
    return snapshots

def _restore_nr_bottom_bed_states(
    cells: list[Cell],
    snapshots: list[dict[str, np.ndarray | float]],
) -> None:
    """Restore bottom-bed NR state captured by ``_snapshot_nr_bottom_bed_states``."""
    for item in snapshots:
        i = int(item["cell_index"])
        if i < 0 or i >= len(cells):
            continue
        cell = cells[i]
        cell.N_d[:] = np.asarray(item["N_d"], dtype=np.float64)
        cell.N_b[:] = np.asarray(item["N_b"], dtype=np.float64)
        cell.T = float(item["T"])
        cell.m_solid[:] = np.asarray(item["m_solid"], dtype=np.float64)
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

def run_outer_preinner_start_value_bridge(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    apply_local_bc_fn: Callable[[int], None] | None = None,
    top_above_cell: Cell | None = None,
    check1_continuation: bool = False,
) -> dict[str, bool | float]:
    """Outer Abgleich Startwert bridge after refresh and before inner NR.

    Merit-driven sequence (when ``nr_outer_preinner_merit_driven_thesis``):
    optional holdup reconcile (T1) → phase-split LSQ (T2) → zone solid/energy (T3)
    → optional bed0 resplit, with global merit gate rollback.

    On Check1 continuation (``check1_continuation=True``), holdup reconcile is
    skipped by default (Bild 2.2: Check1 fail → stay in NR, not Node-A remap).
    Gate rollback restores the pre-mutation snapshot (including upper-bed holdup)
    and remeasures metrics (no sticky inventory / no fake Δmerit=0).

    Ref: Hamel (1999) Bild 2.2 Check1/Check2; engineering Pre-Inner bridge.
    """
    split_cfg = bool(getattr(cfg, "nr_outer_preinner_phase_split_preprojection_thesis", True))
    zone_cfg = bool(getattr(cfg, "nr_outer_preinner_solid_energy_preprojection_thesis", True))
    absorb_cfg = bool(getattr(cfg, "nr_outer_preinner_upper_dense_syngas_absorb_thesis", False))
    local2_cfg = bool(getattr(cfg, "nr_outer_preinner_upper_co_local2x2_thesis", False))
    bed_top_solid_cfg = bool(
        getattr(cfg, "nr_outer_preinner_bed_top_solid_holdup_thesis", False)
    )
    targeted_split_beds = tuple(
        int(i) for i in (getattr(cfg, "nr_outer_preinner_targeted_phase_split_beds_thesis", ()) or ())
    )
    targeted_split_cfg = bool(targeted_split_beds)
    targeted_sum_beds = tuple(
        int(i) for i in (getattr(cfg, "nr_outer_preinner_targeted_species_sum_beds_thesis", ()) or ())
    )
    targeted_sum_cfg = bool(targeted_sum_beds)
    # CO 带宽 / local2 接受次数后冻结 upper 库存桥（§5.13）
    upper_bridge_frozen = bool(getattr(cfg, "_upper_co_bridge_frozen_runtime", False))
    y_co_top = float("nan")
    local2_accept_count = int(getattr(cfg, "_upper_co_local2_accept_count_runtime", 0) or 0)
    freeze_after_n = int(
        max(getattr(cfg, "nr_outer_freeze_upper_co_bridge_after_local2_accepts_thesis", 0) or 0, 0)
    )
    if freeze_after_n > 0 and local2_accept_count >= freeze_after_n:
        setattr(cfg, "_upper_co_bridge_frozen_runtime", True)
        upper_bridge_frozen = True
    if bool(getattr(cfg, "nr_outer_freeze_upper_co_bridge_after_y_co_thesis", False)):
        bed_cells_for_y = [c for c in cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if bed_cells_for_y:
            top = bed_cells_for_y[-1]
            j_co = int(GAS_SPECIES_INDEX["CO"])
            n_co = float(top.N_d[j_co] + top.N_b[j_co])
            n_gas = float(sum(top.N_d) + sum(top.N_b))
            y_co_top = n_co / max(n_gas, 1e-30)
            thr = float(getattr(cfg, "nr_outer_freeze_upper_co_bridge_y_co_threshold_thesis", 0.14))
            if upper_bridge_frozen or y_co_top >= thr:
                setattr(cfg, "_upper_co_bridge_frozen_runtime", True)
                upper_bridge_frozen = True
    if upper_bridge_frozen:
        absorb_cfg = False
        local2_cfg = False
    if (
        not split_cfg
        and not zone_cfg
        and not absorb_cfg
        and not local2_cfg
        and not bed_top_solid_cfg
        and not targeted_split_cfg
        and not targeted_sum_cfg
    ):
        return {
            "enabled": False,
            "upper_co_bridge_frozen": bool(upper_bridge_frozen),
            "y_co_top": float(y_co_top) if np.isfinite(y_co_top) else float("nan"),
            "local2_accept_count": int(local2_accept_count),
        }

    split_max = int(max(getattr(cfg, "nr_outer_preinner_phase_split_max_bed_index_thesis", 0), 0))
    solid_max = int(
        max(getattr(cfg, "nr_outer_preinner_solid_energy_max_bed_index_thesis", 2), 0)
        if zone_cfg
        else 0
    )
    max_bed_index = max(split_max, solid_max)
    use_global_gate = bool(getattr(cfg, "nr_outer_preinner_global_merit_gate_thesis", True))

    apply_bc_fn()
    metrics_before = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    bed_cells = [c for c in cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    flux_mismatch = max_upper_bed_solid_auf_flux_ratio_error(bed_cells, min_bed_index=1)
    stream_imbalance = max_upper_bed_solid_stream_imbalance(bed_cells, min_bed_index=1)
    tier_plan = resolve_preinner_tier_plan(
        metrics_before,
        cfg,
        holdup_flux_mismatch=flux_mismatch,
        holdup_stream_imbalance=stream_imbalance,
    )
    # F4: Bild 2.2 Check1 continuation — do not algebraically remap holdup.
    holdup_skipped_continuation = False
    allow_holdup_on_continuation = bool(
        getattr(cfg, "nr_outer_preinner_holdup_on_check1_continuation_thesis", False)
    )
    if bool(check1_continuation) and not allow_holdup_on_continuation and bool(
        tier_plan.get("run_holdup_reconcile", False)
    ):
        tier_plan = dict(tier_plan)
        tier_plan["run_holdup_reconcile"] = False
        tier_plan["holdup_skipped_check1_continuation"] = True
        holdup_skipped_continuation = True

    split_on = bool(tier_plan["run_split_tier"])
    solid_energy_on = bool(tier_plan["run_solid_energy_tier"])
    run_holdup = bool(tier_plan.get("run_holdup_reconcile", False))
    run_absorb = bool(absorb_cfg)
    run_local2 = bool(local2_cfg)
    if (
        not split_on
        and not solid_energy_on
        and not run_holdup
        and not run_absorb
        and not run_local2
        and not bed_top_solid_cfg
        and not targeted_split_cfg
        and not targeted_sum_cfg
    ):
        return {
            "enabled": False,
            "skipped_all_tiers": True,
            "tier_plan": tier_plan,
            "holdup_skipped_check1_continuation": bool(holdup_skipped_continuation),
            "merit_before": float(metrics_before["merit_scaled"]),
            "merit_after": float(metrics_before["merit_scaled"]),
            "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
            "split_rms_after": float(metrics_before["gas_phase_split_rms"]),
            "solid_rms_before": float(metrics_before["solid_rms"]),
            "solid_rms_after": float(metrics_before["solid_rms"]),
            "energy_rms_before": float(metrics_before["energy_rms"]),
            "energy_rms_after": float(metrics_before["energy_rms"]),
        }

    # F1: snapshot BEFORE mutations. Holdup edits bed1+; include all beds when T1 runs.
    snapshot_max = max_bed_index
    if (run_holdup or run_absorb or run_local2) and bed_cells:
        snapshot_max = max(snapshot_max, len(bed_cells) - 1)
    snapshot = _snapshot_nr_bottom_bed_states(cells, snapshot_max)

    split_shock = bool(tier_plan.get("split_shock", False))
    holdup_before_split = run_holdup and not split_shock
    holdup_after_split = run_holdup and split_shock

    holdup_reconciled = False

    def _run_holdup_reconcile() -> bool:
        return bool(
            reconcile_upper_bed_solid_holdup_for_init(
                bed_cells,
                cfg,
                top_above_cell=top_above_cell,
                allow_reduce=True,
            )
        )

    if holdup_before_split:
        holdup_reconciled = _run_holdup_reconcile()
        if holdup_reconciled:
            apply_bc_fn()
            # Do not refresh snapshot: gate rollback must undo holdup (F1).

    split_diag = run_outer_preinner_gas_phase_split_preprojection(
        cells,
        cfg,
        apply_bc_fn=apply_bc_fn,
        apply_local_bc_fn=apply_local_bc_fn,
        split_shock=split_shock,
    ) if split_on else {"enabled": False}

    absorb_post_gate = bool(
        getattr(cfg, "nr_outer_preinner_upper_dense_absorb_post_gate_thesis", False)
    )
    absorb_diag: dict[str, bool | float | list[int]] = {"enabled": False}
    if run_absorb and not absorb_post_gate:
        absorb_diag = run_outer_preinner_upper_dense_syngas_absorb(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
        )

    # local2 默认 post-gate：与 split/zone 同捆时，后期 split LSQ 会把 rms 打到 O(10)
    # 并连带回滚已接受的 local2（handoff §5.8）。bundle 模式仅 opt-in。
    local2_post_gate = bool(
        getattr(cfg, "nr_outer_preinner_upper_co_local2x2_post_gate_thesis", True)
    )
    local2_diag: dict[str, bool | float | list[int]] = {"enabled": False}
    if run_local2 and not local2_post_gate:
        local2_diag = run_outer_preinner_upper_co_local2x2(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
        )

    if holdup_after_split:
        holdup_after = _run_holdup_reconcile()
        holdup_reconciled = bool(holdup_reconciled or holdup_after)
        if holdup_after:
            apply_bc_fn()

    zone_diag: dict[str, bool | float] = {"enabled": False}
    if solid_energy_on:
        zone_diag = run_bottom_zone_solid_energy_preprojection(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
            apply_local_bc_fn=apply_local_bc_fn,
            solid_repass=bool(getattr(cfg, "nr_outer_preinner_solid_repass_thesis", False)),
            solid_holdup_skip_rms_scaled=float(
                getattr(cfg, "nr_outer_preinner_solid_holdup_skip_rms_scaled_thesis", 0.01)
            ),
        )

    bed0_resplit = False
    if bool(tier_plan.get("run_bed0_resplit", False)):
        bed0_resplit = bool(
            preproject_bottom_major_gas_phase_split_to_exchange_closure(
                cells,
                cfg,
                improvement_factor=0.85,
                apply_bc_fn=apply_bc_fn,
            )
        )
        apply_bc_fn()

    apply_bc_fn()
    metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    global_gate_pass = True
    split_only_fallback = False
    if use_global_gate and not _bottom_joint_global_gate_pass(metrics_before, metrics_after, cfg):
        _restore_nr_bottom_bed_states(cells, snapshot)
        apply_bc_fn()
        # F2: remeasure after rollback (do not copy metrics_before).
        metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
        global_gate_pass = False
        holdup_reconciled = False
        if split_on:
            split_retry = run_outer_preinner_gas_phase_split_preprojection(
                cells,
                cfg,
                apply_bc_fn=apply_bc_fn,
                apply_local_bc_fn=apply_local_bc_fn,
                split_shock=split_shock,
            )
            apply_bc_fn()
            metrics_retry = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
            if _preinner_split_only_fallback_gate_pass(metrics_before, metrics_retry, cfg):
                metrics_after = metrics_retry
                split_diag = split_retry
                global_gate_pass = True
                split_only_fallback = True
            else:
                _restore_nr_bottom_bed_states(cells, snapshot)
                apply_bc_fn()
                metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
                split_diag = {"enabled": bool(split_on)}

    accepted = bool(global_gate_pass)
    # bundle 路径：split_only_fallback 会清掉 snapshot 内的 local2，不得虚报 accepted
    local2_bundle_kept = bool(
        local2_diag.get("accepted", False) and accepted and not split_only_fallback
    )

    if run_local2 and local2_post_gate:
        local2_diag = run_outer_preinner_upper_co_local2x2(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
        )
        if bool(local2_diag.get("accepted", False)):
            apply_bc_fn()
            metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
        local2_kept = bool(local2_diag.get("accepted", False))
    else:
        local2_kept = bool(local2_bundle_kept)

    if local2_kept:
        local2_accept_count = int(local2_accept_count) + 1
        setattr(cfg, "_upper_co_local2_accept_count_runtime", int(local2_accept_count))
        if freeze_after_n > 0 and local2_accept_count >= freeze_after_n:
            setattr(cfg, "_upper_co_bridge_frozen_runtime", True)
            upper_bridge_frozen = True

    absorb_kept = bool(absorb_diag.get("accepted", False)) and accepted and not split_only_fallback
    if run_absorb and absorb_post_gate:
        absorb_diag = run_outer_preinner_upper_dense_syngas_absorb(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
        )
        if bool(absorb_diag.get("accepted", False)):
            apply_bc_fn()
            metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
        absorb_kept = bool(absorb_diag.get("accepted", False))

    # §5.21/§5.23：床顶 solid holdup LSQ（post-gate；可自顶向下多床）
    bed_top_solid_diag: dict[str, bool | float | int | list[int]] = {
        "enabled": False,
        "accepted": False,
        "armed": False,
        "cell_index": -1,
        "n_beds": 1,
        "accepted_count": 0,
        "beds": [],
    }
    if bed_top_solid_cfg and bed_cells:
        if not np.isfinite(y_co_top):
            top = bed_cells[-1]
            j_co = int(GAS_SPECIES_INDEX["CO"])
            n_co = float(top.N_d[j_co] + top.N_b[j_co])
            n_gas = float(sum(top.N_d) + sum(top.N_b))
            y_co_top = n_co / max(n_gas, 1e-30)
        defer_bt = bool(
            getattr(cfg, "nr_outer_preinner_bed_top_solid_holdup_after_y_co_thesis", False)
        )
        thr_bt = float(
            getattr(cfg, "nr_outer_preinner_bed_top_solid_holdup_y_co_threshold_thesis", 0.14)
        )
        n_top = int(
            max(getattr(cfg, "nr_outer_preinner_bed_top_solid_holdup_n_beds_thesis", 1) or 1, 1)
        )
        armed_bt = (not defer_bt) or (np.isfinite(y_co_top) and y_co_top >= thr_bt)
        bed_top_solid_diag["enabled"] = True
        bed_top_solid_diag["armed"] = bool(armed_bt)
        bed_top_solid_diag["n_beds"] = int(n_top)
        bed_top_solid_diag["y_co_top"] = float(y_co_top) if np.isfinite(y_co_top) else float("nan")
        if armed_bt:
            bed_idxs = [
                i
                for i, c in enumerate(cells)
                if str(getattr(c, "cell_type", "bed")) == "bed"
            ]
            target_idxs = bed_idxs[-n_top:] if bed_idxs else []
            bed_top_solid_diag["beds"] = [int(i) for i in target_idxs]
            bed_top_solid_diag["cell_index"] = int(target_idxs[-1]) if target_idxs else -1
            accepted_count = 0
            for idx in target_idxs:
                ok_bt = preproject_bed_solid_holdup_to_transport_closure(
                    cells,
                    cfg,
                    cell_index=int(idx),
                    apply_bc_fn=apply_bc_fn,
                    apply_local_bc_fn=apply_local_bc_fn,
                    skip_rms_scaled=0.0,
                    improvement_factor=0.5,
                )
                if ok_bt:
                    accepted_count += 1
                    apply_bc_fn()
            bed_top_solid_diag["accepted_count"] = int(accepted_count)
            bed_top_solid_diag["accepted"] = bool(accepted_count > 0)
            if accepted_count > 0:
                metrics_after = _measure_nr_init_scaled_metrics(
                    cells, cfg, apply_bc_fn=apply_bc_fn
                )

    # §5.24：指定 bed 的 phase-split LSQ（post-gate；避免扩 max_bed 连带 bed0 爆炸）
    targeted_split_diag: dict[str, bool | float | list[int]] = {
        "enabled": False,
        "armed": False,
        "accepted": False,
        "accepted_count": 0.0,
        "beds": [],
        "beds_accepted": [],
    }
    if targeted_split_cfg:
        if not np.isfinite(y_co_top) and bed_cells:
            top = bed_cells[-1]
            j_co = int(GAS_SPECIES_INDEX["CO"])
            n_co = float(top.N_d[j_co] + top.N_b[j_co])
            n_gas = float(sum(top.N_d) + sum(top.N_b))
            y_co_top = n_co / max(n_gas, 1e-30)
        defer_ts = bool(
            getattr(cfg, "nr_outer_preinner_targeted_phase_split_after_y_co_thesis", False)
        )
        thr_ts = float(
            getattr(cfg, "nr_outer_preinner_targeted_phase_split_y_co_threshold_thesis", 0.14)
        )
        armed_ts = (not defer_ts) or (np.isfinite(y_co_top) and y_co_top >= thr_ts)
        targeted_split_diag["enabled"] = True
        targeted_split_diag["armed"] = bool(armed_ts)
        targeted_split_diag["beds"] = [int(i) for i in targeted_split_beds]
        targeted_split_diag["y_co_top"] = (
            float(y_co_top) if np.isfinite(y_co_top) else float("nan")
        )
        if armed_ts:
            beds_ok: list[int] = []
            for bed_i in targeted_split_beds:
                if bed_i < 0 or bed_i >= len(cells):
                    continue
                if str(getattr(cells[bed_i], "cell_type", "bed")) != "bed":
                    continue
                ok_ts = preproject_bed_major_gas_phase_split_to_exchange_closure(
                    cells,
                    cfg,
                    cell_index=int(bed_i),
                    apply_bc_fn=apply_bc_fn,
                    apply_local_bc_fn=apply_local_bc_fn,
                    improvement_factor=0.5,
                )
                if ok_ts:
                    beds_ok.append(int(bed_i))
                    apply_bc_fn()
            targeted_split_diag["beds_accepted"] = beds_ok
            targeted_split_diag["accepted_count"] = float(len(beds_ok))
            targeted_split_diag["accepted"] = bool(beds_ok)
            if beds_ok:
                metrics_after = _measure_nr_init_scaled_metrics(
                    cells, cfg, apply_bc_fn=apply_bc_fn
                )

    # §5.25：指定 bed 物种总量闭合（含负 res_sum 放出；post-gate）
    targeted_sum_diag: dict[str, bool | float | list[int] | list[str]] = {
        "enabled": False,
        "armed": False,
        "accepted": False,
        "accepted_count": 0.0,
        "beds": [],
        "beds_accepted": [],
        "added_N2_mol_s": 0.0,
    }
    if targeted_sum_cfg:
        if not np.isfinite(y_co_top) and bed_cells:
            top = bed_cells[-1]
            j_co = int(GAS_SPECIES_INDEX["CO"])
            n_co = float(top.N_d[j_co] + top.N_b[j_co])
            n_gas = float(sum(top.N_d) + sum(top.N_b))
            y_co_top = n_co / max(n_gas, 1e-30)
        defer_su = bool(
            getattr(cfg, "nr_outer_preinner_targeted_species_sum_after_y_co_thesis", False)
        )
        thr_su = float(
            getattr(cfg, "nr_outer_preinner_targeted_species_sum_y_co_threshold_thesis", 0.14)
        )
        armed_su = (not defer_su) or (np.isfinite(y_co_top) and y_co_top >= thr_su)
        raw_sp_su = getattr(
            cfg, "nr_outer_preinner_targeted_species_sum_species_thesis", ("N2",)
        )
        species_su = tuple(str(s) for s in (raw_sp_su or ("N2",)))
        frac_su = float(
            getattr(cfg, "nr_outer_preinner_targeted_species_sum_fraction_thesis", 1.0)
        )
        mode_su = str(
            getattr(cfg, "nr_outer_preinner_targeted_species_sum_mode_thesis", "phase_share")
            or "phase_share"
        )
        bidir_su = bool(
            getattr(cfg, "nr_outer_preinner_targeted_species_sum_bidirectional_thesis", True)
        )
        targeted_sum_diag["enabled"] = True
        targeted_sum_diag["armed"] = bool(armed_su)
        targeted_sum_diag["beds"] = [int(i) for i in targeted_sum_beds]
        targeted_sum_diag["species"] = list(species_su)
        if armed_su:
            beds_ok_su: list[int] = []
            added_n2 = 0.0
            for bed_i in targeted_sum_beds:
                if bed_i < 0 or bed_i >= len(cells):
                    continue
                cell_su = cells[bed_i]
                if str(getattr(cell_su, "cell_type", "bed")) != "bed":
                    continue
                cell_su.residuals()
                added_su = oneshot_bed_dense_species_sum_absorb(
                    cell_su,
                    species=species_su,
                    fraction=frac_su,
                    mode=mode_su,
                    bidirectional=bidir_su,
                )
                if any(abs(float(v)) > 0.0 for v in added_su.values()):
                    beds_ok_su.append(int(bed_i))
                    added_n2 += float(added_su.get("N2", 0.0))
                    apply_bc_fn()
            targeted_sum_diag["beds_accepted"] = beds_ok_su
            targeted_sum_diag["accepted_count"] = float(len(beds_ok_su))
            targeted_sum_diag["accepted"] = bool(beds_ok_su)
            targeted_sum_diag["added_N2_mol_s"] = float(added_n2)
            if beds_ok_su:
                metrics_after = _measure_nr_init_scaled_metrics(
                    cells, cfg, apply_bc_fn=apply_bc_fn
                )

    return {
        "enabled": True,
        "max_bed_index": float(max_bed_index),
        "global_gate_pass": bool(global_gate_pass),
        "split_only_fallback": bool(split_only_fallback),
        "merit_driven": bool(tier_plan.get("merit_driven", False)),
        "tier_plan": tier_plan,
        "holdup_reconciled": bool(holdup_reconciled and accepted),
        "holdup_skipped_check1_continuation": bool(holdup_skipped_continuation),
        "holdup_flux_mismatch_before": float(flux_mismatch),
        "merit_before": float(metrics_before["merit_scaled"]),
        "merit_after": float(metrics_after["merit_scaled"]),
        "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
        "split_rms_after": float(metrics_after["gas_phase_split_rms"]),
        "solid_rms_before": float(metrics_before["solid_rms"]),
        "solid_rms_after": float(metrics_after["solid_rms"]),
        "energy_rms_before": float(metrics_before["energy_rms"]),
        "energy_rms_after": float(metrics_after["energy_rms"]),
        "bed0_accepted": bool(split_diag.get("bed0_accepted", False)) and accepted,
        "bed1_split_accepted": bool(split_diag.get("bed1_accepted", False)) and accepted,
        "bed1_energy_t_accepted": bool(
            zone_diag.get("bed1_energy_temperature_preprojection_accepted", False)
        )
        and accepted,
        "solid_preprojections_accepted": float(zone_diag.get("solid_preprojections_accepted", 0.0))
        if accepted
        else 0.0,
        "bed0_resplit_after_solid_energy": bool(bed0_resplit and accepted),
        "upper_dense_absorb_enabled": bool(absorb_diag.get("enabled", False)) or bool(run_absorb),
        "upper_dense_absorb_accepted": bool(absorb_kept),
        "upper_dense_absorb_post_gate": bool(absorb_post_gate),
        "upper_dense_absorb_mode": str(absorb_diag.get("absorb_mode", "") or ""),
        "upper_dense_absorb_beds": list(absorb_diag.get("beds_touched", []) or []),
        "upper_dense_absorb_added_CO_mol_s": float(absorb_diag.get("added_CO_mol_s", 0.0) or 0.0),
        "upper_dense_absorb_added_H2_mol_s": float(absorb_diag.get("added_H2_mol_s", 0.0) or 0.0),
        "upper_co_local2x2_enabled": bool(local2_diag.get("enabled", False)) or bool(run_local2),
        "upper_co_local2x2_accepted": bool(local2_kept),
        "upper_co_local2x2_post_gate": bool(local2_post_gate),
        "upper_co_local2x2_beds": list(local2_diag.get("beds_touched", []) or []),
        "upper_co_local2x2_fraction": float(local2_diag.get("fraction", 0.0) or 0.0),
        "bed_top_solid_holdup_enabled": bool(bed_top_solid_diag.get("enabled", False)),
        "bed_top_solid_holdup_armed": bool(bed_top_solid_diag.get("armed", False)),
        "bed_top_solid_holdup_accepted": bool(bed_top_solid_diag.get("accepted", False)),
        "bed_top_solid_holdup_cell": float(bed_top_solid_diag.get("cell_index", -1)),
        "bed_top_solid_holdup_n_beds": float(bed_top_solid_diag.get("n_beds", 1)),
        "bed_top_solid_holdup_accepted_count": float(
            bed_top_solid_diag.get("accepted_count", 0)
        ),
        "bed_top_solid_holdup_beds": list(bed_top_solid_diag.get("beds") or []),
        "targeted_phase_split_enabled": bool(targeted_split_diag.get("enabled", False)),
        "targeted_phase_split_armed": bool(targeted_split_diag.get("armed", False)),
        "targeted_phase_split_accepted": bool(targeted_split_diag.get("accepted", False)),
        "targeted_phase_split_accepted_count": float(
            targeted_split_diag.get("accepted_count", 0.0)
        ),
        "targeted_phase_split_beds": list(targeted_split_diag.get("beds") or []),
        "targeted_phase_split_beds_accepted": list(
            targeted_split_diag.get("beds_accepted") or []
        ),
        "targeted_species_sum_enabled": bool(targeted_sum_diag.get("enabled", False)),
        "targeted_species_sum_armed": bool(targeted_sum_diag.get("armed", False)),
        "targeted_species_sum_accepted": bool(targeted_sum_diag.get("accepted", False)),
        "targeted_species_sum_accepted_count": float(
            targeted_sum_diag.get("accepted_count", 0.0)
        ),
        "targeted_species_sum_beds": list(targeted_sum_diag.get("beds") or []),
        "targeted_species_sum_beds_accepted": list(
            targeted_sum_diag.get("beds_accepted") or []
        ),
        "targeted_species_sum_added_N2_mol_s": float(
            targeted_sum_diag.get("added_N2_mol_s", 0.0)
        ),
        "upper_co_bridge_frozen": bool(upper_bridge_frozen),
        "y_co_top": float(y_co_top) if np.isfinite(y_co_top) else float("nan"),
        "local2_accept_count": int(local2_accept_count),
    }

def preproject_bottom_bc_aware_bed0_bed1_joint_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    gas_relax: float = 1.0,
    apply_bc_fn: Callable[[], None] | None = None,
    temperature_bounds_K: tuple[float, float] = (650.0, 1600.0),
    improvement_factor: float = 0.75,
) -> bool:
    """BC-aware bottom-zone joint x0 via alpha-blended bed0 projection + bed1 T response."""
    if apply_bc_fn is None or len(cells) < 2:
        return preproject_bottom_total_gas_and_temperature_to_energy_closure(
            cells,
            cfg,
            gas_relax=gas_relax,
            temperature_bounds_K=temperature_bounds_K,
            improvement_factor=improvement_factor,
        )
    bed0 = cells[0]
    snap = (bed0.N_d.copy(), bed0.N_b.copy(), float(bed0.T))
    snap_bed1_t = float(cells[1].T)

    def _restore() -> None:
        bed0.N_d[:] = snap[0]
        bed0.N_b[:] = snap[1]
        bed0.T = snap[2]
        bed0._thermo_cache_valid = False
        bed0._hydro_cache_valid = False
        cells[1].T = snap_bed1_t
        cells[1]._thermo_cache_valid = False
        cells[1]._hydro_cache_valid = False

    metrics_before = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
    ok, _alpha = _attempt_bc_aware_bottom_joint_projection(
        cells,
        cfg,
        apply_bc_fn=apply_bc_fn,
        metrics_before=metrics_before,
        restore_fn=_restore,
        snap_bed0=snap,
        bed1_energy_before=_bottom_neighbor_energy_residual_abs(cells),
        guard_ratio=float(max(getattr(cfg, "nr_bottom_joint_bed1_energy_guard_ratio_thesis", 1.25), 1.0)),
    )
    if not ok:
        _restore()
        if apply_bc_fn is not None:
            apply_bc_fn()
    return bool(ok)

def preproject_bed_temperature_to_local_energy_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    cell_index: int,
    apply_bc_fn: Callable[[], None],
    apply_local_bc_fn: Callable[[int], None] | None = None,
    improvement_factor: float = 0.75,
    temperature_bounds_K: tuple[float, float] = (650.0, 1600.0),
) -> bool:
    """Single-cell temperature projection after BC refresh (downstream energy response)."""
    if cell_index < 0 or cell_index >= len(cells):
        return False
    cell = cells[cell_index]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    bc_ctx = PreprojectionBcContext(
        cell_index=int(cell_index),
        apply_bc_fn=apply_bc_fn,
        apply_local_bc_fn=apply_local_bc_fn,
    )
    bc_ctx.finalize()
    original_T = float(cell.T)
    initial_energy_abs = _preproj_cell_energy_residual_abs(cell)
    if not np.isfinite(initial_energy_abs) or initial_energy_abs <= 1e-9:
        return False

    attr_lo = float(getattr(cell, "_nr_temperature_min_K", temperature_bounds_K[0]))
    attr_hi = float(getattr(cell, "_nr_temperature_max_K", temperature_bounds_K[1]))
    t_lo = max(float(temperature_bounds_K[0]), attr_lo, 300.0)
    t_hi = min(float(temperature_bounds_K[1]), attr_hi, 2500.0)
    t_hi = max(t_hi, t_lo + 1.0)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    bed1_energy_lsq_max_nfev = int(
        max(getattr(cfg, "nr_bottom_zone_bed1_energy_lsq_max_nfev_thesis", 12), 4)
    )

    def _restore() -> None:
        cell.T = original_T
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _apply_temperature(temperature_K: float) -> None:
        cell.T = float(np.clip(temperature_K, t_lo, t_hi))
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _set_values(values: np.ndarray) -> None:
        _apply_temperature(float(values[0]))

    def _objective(x: np.ndarray) -> np.ndarray:
        _apply_temperature(float(x[0]))
        bc_ctx.refresh_trial()
        return np.array([_preproj_cell_energy_residual_abs(cell) / ref_energy], dtype=np.float64)

    x0 = np.array([float(np.clip(original_T, t_lo, t_hi))], dtype=np.float64)

    def _record_accept_diag() -> None:
        cell._bottom_neighbor_preprojection_diag = {
            "cell_index": int(cell_index),
            "initial_energy_abs_W": float(initial_energy_abs),
            "final_energy_abs_W": float(_preproj_cell_energy_residual_abs(cell)),
            "candidate_temperature_K": float(cell.T),
        }

    return run_improvement_gated_lsq(
        spec=LsqPreprojectionSpec(
            x0=x0,
            lower=np.array([t_lo], dtype=np.float64),
            upper=np.array([t_hi], dtype=np.float64),
            max_nfev=bed1_energy_lsq_max_nfev,
            xtol=1.0e-8,
            ftol=1.0e-8,
            gtol=1.0e-8,
        ),
        set_values=_set_values,
        objective=_objective,
        initial_metric=initial_energy_abs,
        final_metric=lambda: _preproj_cell_energy_residual_abs(cell),
        improvement_factor=float(improvement_factor),
        restore=_restore,
        sync_full_residuals=cell.residuals,
        bc=bc_ctx,
        finalize_on_accept=False,
        sync_on_accept=False,
        on_accept=_record_accept_diag,
    )

def preproject_bottom_total_gas_and_temperature_to_energy_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    gas_relax: float = 1.0,
    temperature_bounds_K: tuple[float, float] = (650.0, 1600.0),
    improvement_factor: float = 0.75,
) -> bool:
    """Joint bed0 x0 projection wrapper around ``preproject_bed_gas_temperature_to_energy_closure``."""
    return preproject_bed_gas_temperature_to_energy_closure(
        cells,
        cfg,
        cell_index=0,
        gas_relax=gas_relax,
        temperature_bounds_K=temperature_bounds_K,
        improvement_factor=improvement_factor,
    )

def _attempt_bc_aware_bottom_joint_projection(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    metrics_before: dict[str, float] | None,
    restore_fn: Callable[[], None],
    snap_bed0: tuple[np.ndarray, np.ndarray, float],
    bed1_energy_before: float | None,
    guard_ratio: float,
) -> tuple[bool, float]:
    """Alpha-blend bed0 projection direction under BC-refreshed global/bed1 guards."""
    bed0 = cells[0]
    base_d, base_b, base_T = snap_bed0
    restore_fn()
    apply_bc_fn()
    if not preproject_bottom_total_gas_and_temperature_to_energy_closure(cells, cfg, gas_relax=1.0):
        restore_fn()
        apply_bc_fn()
        return False, 0.0
    target_d, target_b, target_T = bed0.N_d.copy(), bed0.N_b.copy(), float(bed0.T)
    restore_fn()
    apply_bc_fn()

    max_increase = max((guard_ratio - 1.0) * max(bed1_energy_before or 0.0, 1.0e3), 1.0e6)

    def _try_bc_aware_omega(omega: float) -> bool:
        _blend_bed0_state_from_snapshot(
            bed0,
            base_d=base_d,
            base_b=base_b,
            base_T=base_T,
            target_d=target_d,
            target_b=target_b,
            target_T=target_T,
            omega=omega,
        )
        apply_bc_fn()
        preproject_bed_temperature_to_local_energy_closure(
            cells,
            cfg,
            cell_index=1,
            apply_bc_fn=apply_bc_fn,
        )
        apply_bc_fn()
        bed1_energy_after = _bottom_neighbor_energy_residual_abs(cells)
        guard_ok = (
            bed1_energy_before is None
            or bed1_energy_after is None
            or bed1_energy_after <= bed1_energy_before + max_increase
        )
        if metrics_before is not None:
            metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=None)
            gate_ok = _bottom_joint_global_gate_pass(metrics_before, metrics_after, cfg)
            if gate_ok and guard_ok:
                bed0._bottom_joint_preprojection_diag = {
                    **(getattr(bed0, "_bottom_joint_preprojection_diag", None) or {}),
                    "bc_aware_neighbor_index": 1,
                    "bc_aware_blend_alpha": float(omega),
                    "initial_neighbor_energy_abs_W": float(bed1_energy_before or 0.0),
                    "final_neighbor_energy_abs_W": float(bed1_energy_after or 0.0),
                    "initial_merit_scaled": float(metrics_before.get("merit_scaled", metrics_before["rms_scaled"])),
                    "final_merit_scaled": float(metrics_after.get("merit_scaled", metrics_after["rms_scaled"])),
                }
                return True
            return False
        return guard_ok

    ok, alpha = run_relax_schedule(
        _bottom_joint_relax_candidates(cfg),
        restore=restore_fn,
        prepare_trial=apply_bc_fn,
        try_trial=_try_bc_aware_omega,
        finalize_failure=lambda: (restore_fn(), apply_bc_fn()),
    )
    if ok:
        return True, alpha
    return False, 0.0

def run_bottom_joint_start_value_preprojection(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Any | None = None,
) -> dict[str, bool | float]:
    """Joint bed0 gas/T preprojection with relax schedule and optional global RMS gate."""
    empty_diag = {
        "joint": False,
        "joint_gas_relax": 0.0,
        "bed1_guard_tripped": False,
        "joint_global_gate_pass": False,
        "joint_rms_scaled_before": 0.0,
        "joint_rms_scaled_after": 0.0,
        "joint_merit_scaled_before": 0.0,
        "joint_merit_scaled_after": 0.0,
        "joint_split_rms_before": 0.0,
        "joint_split_rms_after": 0.0,
        "joint_split_aware_gate": False,
        "bc_aware_bed0_bed1_joint": False,
        "bed1_energy_abs_before_W": 0.0,
        "bed1_energy_abs_after_W": 0.0,
    }
    if not cells:
        return empty_diag
    cell = cells[0]
    snapshot_d = cell.N_d.copy()
    snapshot_b = cell.N_b.copy()
    snapshot_T = float(cell.T)
    snapshot_bed1_T: float | None = None
    if len(cells) >= 2 and str(getattr(cells[1], "cell_type", "bed")) == "bed":
        snapshot_bed1_T = float(cells[1].T)

    use_bc_aware = bool(getattr(cfg, "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis", True))

    def _restore_bed0_snapshot() -> None:
        cell.N_d[:] = snapshot_d
        cell.N_b[:] = snapshot_b
        cell.T = snapshot_T
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False
        if snapshot_bed1_T is not None:
            cells[1].T = snapshot_bed1_T
            cells[1]._thermo_cache_valid = False
            cells[1]._hydro_cache_valid = False

    snap_bed0 = (snapshot_d, snapshot_b, snapshot_T)

    use_global_gate = bool(getattr(cfg, "nr_bottom_joint_global_rms_gate_thesis", True))
    use_split_aware_gate = bool(getattr(cfg, "nr_bottom_joint_split_aware_gate_thesis", True))
    metrics_before = (
        _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=apply_bc_fn)
        if use_global_gate
        else None
    )
    bed1_energy_before = _bottom_neighbor_energy_residual_abs(cells)
    guard_ratio = float(max(getattr(cfg, "nr_bottom_joint_bed1_energy_guard_ratio_thesis", 1.25), 1.0))
    joint_ok = False
    accepted_relax = 0.0
    bed1_guard_tripped = False
    global_gate_pass = False
    metrics_after = metrics_before

    if use_bc_aware and apply_bc_fn is not None and len(cells) >= 2:
        _restore_bed0_snapshot()
        if apply_bc_fn is not None:
            apply_bc_fn()
        bc_ok, bc_alpha = _attempt_bc_aware_bottom_joint_projection(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
            metrics_before=metrics_before if use_global_gate else None,
            restore_fn=_restore_bed0_snapshot,
            snap_bed0=snap_bed0,
            bed1_energy_before=bed1_energy_before,
            guard_ratio=guard_ratio,
        )
        if bc_ok:
            joint_ok = True
            accepted_relax = bc_alpha
            global_gate_pass = bool(use_global_gate)
            if use_global_gate:
                metrics_after = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=None)
        else:
            _restore_bed0_snapshot()
            if apply_bc_fn is not None:
                apply_bc_fn()

    if not joint_ok:
        trial_state = {
            "global_gate_pass": False,
            "bed1_guard_tripped": False,
            "metrics_after": metrics_before,
        }

        def _prepare_joint_trial() -> None:
            if apply_bc_fn is not None:
                apply_bc_fn()

        def _try_joint_relax(trial_relax: float) -> bool:
            if not preproject_bottom_total_gas_and_temperature_to_energy_closure(
                cells,
                cfg,
                gas_relax=trial_relax,
            ):
                return False
            if apply_bc_fn is not None:
                apply_bc_fn()
            if use_global_gate:
                metrics_after_trial = _measure_nr_init_scaled_metrics(cells, cfg, apply_bc_fn=None)
                trial_state["metrics_after"] = metrics_after_trial
                if _bottom_joint_global_gate_pass(metrics_before, metrics_after_trial, cfg):
                    trial_state["global_gate_pass"] = True
                    return True
                return False
            bed1_energy_after = _bottom_neighbor_energy_residual_abs(cells)
            max_increase = max((guard_ratio - 1.0) * max(bed1_energy_before or 0.0, 1.0e3), 1.0e6)
            guard_ok = (
                bed1_energy_before is None
                or bed1_energy_after is None
                or bed1_energy_after <= bed1_energy_before + max_increase
            )
            if guard_ok:
                return True
            if trial_relax >= 1.0 - 1.0e-12:
                trial_state["bed1_guard_tripped"] = True
                return True
            return False

        def _finalize_joint_failure() -> None:
            if apply_bc_fn is not None:
                apply_bc_fn()

        relax_ok, trial_relax = run_relax_schedule(
            _bottom_joint_relax_candidates(cfg),
            restore=_restore_bed0_snapshot,
            prepare_trial=_prepare_joint_trial,
            try_trial=_try_joint_relax,
            finalize_failure=_finalize_joint_failure,
        )
        if relax_ok:
            joint_ok = True
            accepted_relax = trial_relax
            global_gate_pass = bool(trial_state["global_gate_pass"])
            bed1_guard_tripped = bool(trial_state["bed1_guard_tripped"])
            if use_global_gate and trial_state["metrics_after"] is not None:
                metrics_after = trial_state["metrics_after"]
    if not joint_ok:
        _restore_bed0_snapshot()
        if apply_bc_fn is not None:
            apply_bc_fn()
        metrics_after = metrics_before if metrics_before is not None else metrics_after
    return {
        "joint": bool(joint_ok),
        "joint_gas_relax": float(accepted_relax if joint_ok else 0.0),
        "bed1_guard_tripped": bool(bed1_guard_tripped),
        "joint_global_gate_pass": bool(global_gate_pass),
        "joint_rms_scaled_before": float((metrics_before or {}).get("rms_scaled", 0.0)),
        "joint_rms_scaled_after": float((metrics_after or {}).get("rms_scaled", 0.0)),
        "joint_merit_scaled_before": float((metrics_before or {}).get("merit_scaled", 0.0)),
        "joint_merit_scaled_after": float((metrics_after or {}).get("merit_scaled", 0.0)),
        "joint_split_rms_before": float((metrics_before or {}).get("gas_phase_split_rms", 0.0)),
        "joint_split_rms_after": float((metrics_after or {}).get("gas_phase_split_rms", 0.0)),
        "joint_split_aware_gate": bool(use_split_aware_gate and use_global_gate),
        "bc_aware_bed0_bed1_joint": bool(use_bc_aware and apply_bc_fn is not None and len(cells) >= 2),
        "bed1_energy_abs_before_W": float(bed1_energy_before or 0.0),
        "bed1_energy_abs_after_W": float(_bottom_neighbor_energy_residual_abs(cells) or 0.0),
    }

def run_bottom_start_value_preprojections(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Any | None = None,
) -> dict[str, bool | float]:
    """Run all bottom-bed x0 bridges (primary, major, guarded joint)."""
    primary_ok = preproject_bottom_primary_gas_state_to_exchange_closure(cells, cfg)
    major_ok = preproject_bottom_major_gas_state_to_local_balance(cells, cfg)
    joint = run_bottom_joint_start_value_preprojection(cells, cfg, apply_bc_fn=apply_bc_fn)
    return {
        "primary": bool(primary_ok),
        "major": bool(major_ok),
        **joint,
    }

