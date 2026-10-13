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

from src.core.preprojection_engine import (
    preproj_cell_energy_residual_abs as _preproj_cell_energy_residual_abs,
)

def _bottom_joint_gate_merit(metrics: dict[str, float], cfg: Any) -> float:
    """Scalar gate merit aligned with thesis NR line-search (RMS, split, energy max)."""
    merit = float(metrics["rms_scaled"])
    if bool(getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)):
        merit = max(merit, float(metrics.get("gas_phase_split_rms", 0.0)))
    if bool(getattr(cfg, "nr_bottom_joint_energy_aware_gate_thesis", True)):
        merit = max(merit, float(metrics.get("energy_rms", 0.0)))
    return merit

def resolve_preinner_tier_plan(
    metrics: dict[str, float],
    cfg: Any,
    *,
    holdup_flux_mismatch: float = 0.0,
    holdup_stream_imbalance: float = 0.0,
) -> dict[str, bool | float]:
    """Decide which pre-inner projection tiers to run from current scaled metrics.

    Tier T2: major phase-split LSQ when ``gas_phase_split_rms`` or
    ``gas_phase_split_max_abs`` exceeds threshold.
    Tier T1: algebraic upper-bed holdup reconcile when energy merit is elevated,
    a single-bed energy hotspot exists, or auf-flux closure is violated on bed1+.
    Tier T3: bottom-zone solid/energy LSQ when energy or solid RMS exceeds threshold
    (skipped during refresh split shock when split alone dominates merit).

    Ref: docs/vorab_init_abgleich_vs_outer_preinner.md; merit-driven pre-inner policy.
    """
    split_cfg_on = bool(getattr(cfg, "nr_outer_preinner_phase_split_preprojection_thesis", True))
    zone_cfg_on = bool(getattr(cfg, "nr_outer_preinner_solid_energy_preprojection_thesis", True))
    merit_driven = bool(getattr(cfg, "nr_outer_preinner_merit_driven_thesis", True))

    split_rms = float(metrics.get("gas_phase_split_rms", 0.0))
    split_max_abs = float(metrics.get("gas_phase_split_max_abs", 0.0))
    energy_rms = float(metrics.get("energy_rms", 0.0))
    solid_rms = float(metrics.get("solid_rms", 0.0))
    merit = float(metrics.get("merit_scaled", 0.0))
    tau_split = float(getattr(cfg, "nr_outer_preinner_split_tier_min_rms_thesis", 0.03))
    tau_split_max = float(
        getattr(cfg, "nr_outer_preinner_split_tier_min_max_abs_thesis", 0.05)
    )
    tau_energy = float(getattr(cfg, "nr_outer_preinner_energy_tier_min_rms_thesis", 0.05))
    tau_solid = float(getattr(cfg, "nr_outer_preinner_solid_tier_min_rms_thesis", 0.01))
    tau_split_shock = float(getattr(cfg, "nr_outer_preinner_split_shock_rms_thesis", 0.15))
    tau_max_bed_energy = float(
        getattr(cfg, "nr_outer_preinner_max_bed_energy_scaled_thesis", 0.08)
    )
    tau_flux_mismatch = float(
        getattr(cfg, "nr_outer_preinner_holdup_flux_mismatch_min_thesis", 0.15)
    )
    tau_stream_imbalance = float(
        getattr(cfg, "nr_outer_preinner_holdup_stream_imbalance_min_thesis", 0.15)
    )
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    max_bed_energy_scaled = float(metrics.get("max_bed_energy_abs_W", 0.0)) / ref_energy
    flux_mismatch = float(max(0.0, holdup_flux_mismatch))
    stream_imbalance = float(max(0.0, holdup_stream_imbalance))

    if not merit_driven:
        return {
            "merit_driven": False,
            "run_split_tier": bool(split_cfg_on),
            "run_solid_energy_tier": bool(zone_cfg_on),
            "run_holdup_reconcile": False,
            "run_bed0_resplit": bool(split_cfg_on and zone_cfg_on),
            "split_rms": split_rms,
            "gas_phase_split_max_abs": split_max_abs,
            "energy_rms": energy_rms,
            "solid_rms": solid_rms,
            "merit_scaled": merit,
        }

    energy_elevated = energy_rms >= tau_energy or (
        energy_rms >= max(0.5 * merit, tau_solid) and energy_rms >= tau_solid
    )
    bed_energy_hotspot = max_bed_energy_scaled >= tau_max_bed_energy
    holdup_flux_elevated = flux_mismatch >= tau_flux_mismatch
    stream_imbalance_elevated = stream_imbalance >= tau_stream_imbalance
    holdup_transport_mismatch = holdup_flux_elevated and stream_imbalance_elevated
    split_shock = split_rms >= tau_split_shock
    split_dominates_merit = merit > 1.0e-12 and split_rms >= 0.85 * merit
    if split_shock and split_dominates_merit:
        holdup_trigger = holdup_transport_mismatch
    else:
        holdup_trigger = energy_elevated or bed_energy_hotspot or holdup_transport_mismatch
    run_split = bool(
        split_cfg_on and (split_rms >= tau_split or split_max_abs >= tau_split_max)
    )
    run_reconcile = bool(
        getattr(cfg, "nr_outer_preinner_reconcile_holdup_before_zone_thesis", True)
        and holdup_trigger
    )
    run_zone = bool(
        zone_cfg_on
        and (energy_elevated or solid_rms >= tau_solid or bed_energy_hotspot)
        and not split_shock
        and not (split_dominates_merit and split_rms > energy_rms)
    )
    return {
        "merit_driven": True,
        "run_split_tier": run_split,
        "run_solid_energy_tier": run_zone,
        "run_holdup_reconcile": run_reconcile,
        "run_bed0_resplit": bool(run_split and run_zone),
        "split_rms": split_rms,
        "gas_phase_split_max_abs": split_max_abs,
        "energy_rms": energy_rms,
        "solid_rms": solid_rms,
        "merit_scaled": merit,
        "max_bed_energy_scaled": max_bed_energy_scaled,
        "holdup_flux_mismatch": flux_mismatch,
        "holdup_stream_imbalance": stream_imbalance,
        "holdup_transport_mismatch": bool(holdup_transport_mismatch),
        "split_shock": bool(split_shock),
        "split_dominates_merit": bool(split_dominates_merit),
        "tau_split": tau_split,
        "tau_split_max_abs": tau_split_max,
        "tau_energy": tau_energy,
        "tau_solid": tau_solid,
    }


def _bottom_joint_global_gate_pass(
    metrics_before: dict[str, float],
    metrics_after: dict[str, float],
    cfg: Any,
) -> bool:
    """Accept bottom joint x0 when global merit does not increase (split-aware in thesis)."""
    eps = 1.0e-12
    use_split_aware = bool(getattr(cfg, "nr_bottom_joint_split_aware_gate_thesis", True))
    if use_split_aware and bool(getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)):
        merit_before = _bottom_joint_gate_merit(metrics_before, cfg)
        merit_after = _bottom_joint_gate_merit(metrics_after, cfg)
        return np.isfinite(merit_after) and merit_after <= merit_before + eps
    return (
        metrics_after["rms_scaled"] <= metrics_before["rms_scaled"] + eps
        and metrics_after["max_abs_scaled"] <= metrics_before["max_abs_scaled"] + eps
    )

def _measure_nr_init_scaled_metrics(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Any | None,
) -> dict[str, float]:
    """Measure scaled NR residual metrics after BC refresh (matches inner NR entry)."""
    from src.solvers.global_nr_solver import (
        _residual_gas_phase_metrics,
        _residual_group_metrics,
        _rms_norm,
        build_equation_scales,
        global_residual,
        pack_reactor,
    )

    if apply_bc_fn is not None:
        apply_bc_fn()
    ref_gas = max(float(getattr(cfg, "O2_feed", 0.0) + getattr(cfg, "H2O_feed", 0.0) + getattr(cfg, "N2_feed", 0.0)), 1.0)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    F = global_residual(pack_reactor(cells), cells, apply_bc_fn or (lambda: None))
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    gas_phase = _residual_gas_phase_metrics(F_hat, cells)
    group_metrics = _residual_group_metrics(F_hat, cells)
    max_bed_energy_W = 0.0
    for cell in cells:
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        val = _preproj_cell_energy_residual_abs(cell)
        if np.isfinite(val):
            max_bed_energy_W = max(max_bed_energy_W, val)
    metrics = {
        "rms_scaled": float(_rms_norm(F_hat)),
        "max_abs_scaled": float(np.max(np.abs(F_hat))) if F_hat.size else 0.0,
        "max_bed_energy_abs_W": float(max_bed_energy_W),
        "gas_phase_split_rms": float(gas_phase["gas_phase_split_rms"]),
        "gas_phase_split_max_abs": float(gas_phase["gas_phase_split_max_abs"]),
        "solid_rms": float(group_metrics["solid_rms"]),
        "energy_rms": float(group_metrics["energy_rms"]),
    }
    metrics["merit_scaled"] = float(_bottom_joint_gate_merit(metrics, cfg))
    metrics["max_bed_energy_scaled"] = float(max_bed_energy_W / ref_energy)
    return metrics

def _bottom_neighbor_energy_residual_abs(cells: list[Cell]) -> float | None:
    if len(cells) < 2:
        return None
    neighbor = cells[1]
    if str(getattr(neighbor, "cell_type", "bed")) != "bed":
        return None
    val = _preproj_cell_energy_residual_abs(neighbor)
    return val if np.isfinite(val) else None

def _blend_bed0_state_from_snapshot(
    cell: Cell,
    *,
    base_d: np.ndarray,
    base_b: np.ndarray,
    base_T: float,
    target_d: np.ndarray,
    target_b: np.ndarray,
    target_T: float,
    omega: float,
) -> None:
    """Convex blend between baseline and projected bed0 state."""
    blend = float(np.clip(omega, 0.0, 1.0))
    cell.N_d[:] = np.maximum(base_d + blend * (target_d - base_d), 0.0)
    cell.N_b[:] = np.maximum(base_b + blend * (target_b - base_b), 0.0)
    cell.T = float(base_T + blend * (target_T - base_T))
    cell._thermo_cache_valid = False
    cell._hydro_cache_valid = False

def _resolve_bottom_joint_gas_relax_min(cfg: Any) -> float:
    if not bool(getattr(cfg, "thesis_mode", False)):
        return 1.0
    return float(
        np.clip(
            getattr(cfg, "nr_bottom_joint_preprojection_gas_relax_min_thesis", 0.05),
            0.01,
            1.0,
        )
    )

def _bottom_joint_relax_candidates(cfg: Any) -> list[float]:
    """Relax schedule: full projection first, then halving down to thesis minimum."""
    min_omega = _resolve_bottom_joint_gas_relax_min(cfg)
    candidates: list[float] = []
    trial = 1.0
    while trial >= min_omega - 1.0e-12:
        candidates.append(float(trial))
        trial *= 0.5
    if min_omega not in candidates:
        candidates.append(float(min_omega))
    return sorted(set(candidates), reverse=True)

