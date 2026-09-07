"""Vorabrechnung — vorab refresh segment (OPT-004)."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Dict, List

import numpy as np

from src.core.cell import Cell, S_CHAR, S_VM, S_MOISTURE, S_ASH
from src.core.connectivity import cell_total_solid_holdup
from src.core.species import GAS_SPECIES_INDEX, get_atom_count, gibbs_molar
from src.thermal.devolatilization import (
    devolatilization_rate_for_cell,
    vm_devolatilization_zone_cumulative_fraction,
    vm_devolatilization_zone_increment_fraction,
)
from src.thermal.drying import solve_drying_CN
from src.thermodynamics.gibbs_hamel_reduced import ReducedHamelGibbsSolver
from src.thermodynamics.gibbs_minimizer import GibbsMinimizer
from src.solvers.vorabrechnung.vorab_budgets import vorabrechnung_tau_for_cell

def refresh_cell_vorabrechnung(cell: Cell, *, force: bool = False) -> None:
    """按 Hamel 外层 Abgleich 语义刷新单格 Vorabrechnung。

    顺序约束：
    1. 先重算 hydrodynamics（u_mf, u_b, eps_b, K_bd...）
    2. 再按最新流体力学停留时间计算 drying/pyrolysis 缓存
    3. inner NR 期间只读这些缓存，不在残差扰动中重复重算

    thesis single-shot 冻结后仅刷新水力学快照，不重算干燥/DAEM。
    """
    if bool(getattr(cell, "_vorab_drying_pyro_sources_frozen", False)):
        refresh_cell_hydrodynamics_for_frozen_inner(cell)
        return
    cell.calc_hydrodynamics()
    if force:
        cell.invalidate_vorabrechnung_cache()
    cell.compute_vorabrechnung(vorabrechnung_tau_for_cell(cell))


def refresh_cell_hydrodynamics_for_frozen_inner(cell: Cell) -> None:
    """Refresh only hydrodynamics scalars and snapshot for inner-NR freeze.

    Does **not** recompute bed ``K_solid_*`` (needs the full bed chain). Prefer
    ``refresh_cells_hydrodynamics_for_frozen_inner`` on accept-step / outer hydro
    refresh so transport coefficients stay consistent with the new hydro state.

    This path intentionally keeps already prepared drying/pyrolysis source caches intact.
    """
    cell.calc_hydrodynamics()
    cell._snapshot_vorabrechnung_hydrodynamics()


def _update_bed_solid_transport_after_hydro_refresh(cells: List[Cell], cfg: Any | None = None) -> None:
    """Recompute bed ``K_auf`` / ``K_ab`` from the just-refreshed hydro fields.

    Under inner-NR freeze, ``apply_bc`` restores snapshot ``K_solid`` and skips
    ``update_bed_solid_transport_coefficients``. Hydro-only refresh therefore left
    stale transport coefficients in the freeze cache (seed850 O1: freeze E≈0.23 vs
    live E≈0.44 after large ΔT).

    Ref: Bild 2.2 outer Abgleich — hydro/transport must track accepted state updates.
    """
    from src.core.connectivity.solid_transport import (
        update_bed_solid_transport_coefficients,
        update_side_block_solid_transport_coefficients,
    )

    bed_cells = [c for c in cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    if not bed_cells:
        return
    freeboard_cells = [c for c in cells if str(getattr(c, "cell_type", "")) == "freeboard"]
    top_downflow_total = 0.0
    top_downflow_components: np.ndarray | None = None
    if freeboard_cells:
        top_downflow_matrix = np.maximum(freeboard_cells[0]._solid_downflow_rates(), 0.0)
        top_downflow_total = float(np.sum(top_downflow_matrix))
        top_downflow_components = np.array(
            [
                float(np.sum(top_downflow_matrix[:, S_CHAR])),
                float(np.sum(top_downflow_matrix[:, S_ASH])),
            ],
            dtype=np.float64,
        )
    update_bed_solid_transport_coefficients(
        bed_cells,
        cfg if cfg is not None else object(),
        top_downflow_total=top_downflow_total,
        top_downflow_components=top_downflow_components,
    )
    cyclone_cells = [c for c in cells if str(getattr(c, "cell_type", "")) == "cyclone"]
    return_cells = [c for c in cells if str(getattr(c, "cell_type", "")) == "return_leg"]
    if cyclone_cells and return_cells and cfg is not None:
        update_side_block_solid_transport_coefficients(
            cyclone_cells[0],
            return_cells[0],
            cfg,
        )


def refresh_cells_hydrodynamics_for_frozen_inner(
    cells: List[Cell],
    *,
    cfg: Any | None = None,
) -> None:
    """Refresh hydrodynamics + bed solid transport K, then snapshot for freeze.

    Call this on accepted Newton steps and outer hydro-only refresh. Per-cell
    ``refresh_cell_hydrodynamics_for_frozen_inner`` alone is insufficient because
    ``K_solid_*`` are chain closures over the bed, not single-cell hydro outputs.
    """
    for cell in cells:
        cell.calc_hydrodynamics()
    _update_bed_solid_transport_after_hydro_refresh(cells, cfg)
    for cell in cells:
        cell._snapshot_vorabrechnung_hydrodynamics()


def backup_vorabrechnung_hydro_freeze(
    cells: List[Cell],
) -> list[tuple[dict[str, float], np.ndarray, np.ndarray, bool, bool]]:
    """Deep-copy freeze-cache hydro/transport for refresh-aware LS probe rollback."""
    backups: list[tuple[dict[str, float], np.ndarray, np.ndarray, bool, bool]] = []
    for cell in cells:
        cache = {str(k): float(v) for k, v in dict(cell._vorab_hydro_cache).items()}
        backups.append(
            (
                cache,
                np.array(cell._vorab_K_solid_auf, dtype=np.float64, copy=True),
                np.array(cell._vorab_K_solid_ab, dtype=np.float64, copy=True),
                bool(cell._vorab_hydro_cache_valid),
                bool(cell._vorab_solid_transport_cache_valid),
            )
        )
    return backups


def restore_vorabrechnung_hydro_freeze(
    cells: List[Cell],
    backups: list[tuple[dict[str, float], np.ndarray, np.ndarray, bool, bool]],
) -> None:
    """Restore freeze cache and live hydro fields after a rejected refresh probe."""
    if len(backups) != len(cells):
        raise ValueError("hydro freeze backup length must match cells")
    for cell, (cache, k_auf, k_ab, hydro_ok, solid_ok) in zip(cells, backups):
        cell._vorab_hydro_cache.clear()
        cell._vorab_hydro_cache.update(cache)
        cell._vorab_K_solid_auf[:, :] = k_auf
        cell._vorab_K_solid_ab[:, :] = k_ab
        cell._vorab_hydro_cache_valid = bool(hydro_ok)
        cell._vorab_solid_transport_cache_valid = bool(solid_ok)
        if bool(hydro_ok):
            cell._apply_frozen_vorabrechnung_hydrodynamics()
        cell._hydro_cache_valid = False


def _hydro_field_vector_from_mapping(
    *,
    u_mf: float,
    u_b: float,
    d_b: float,
    eps_b: float,
    eps_d: float,
    K_bd: float,
    V_b: float,
    V_d: float,
    u0: float,
    u_d: float,
    K_solid_auf_abs: float,
    K_solid_ab_abs: float,
) -> list[float]:
    """Shared layout for freeze/live hydro fingerprints."""
    return [
        float(u_mf),
        float(u_b),
        float(d_b),
        float(eps_b),
        float(eps_d),
        float(K_bd),
        float(V_b),
        float(V_d),
        float(u0),
        float(u_d),
        float(K_solid_auf_abs),
        float(K_solid_ab_abs),
    ]


def hydrodynamics_freeze_fingerprint(cells: List[Cell]) -> dict[str, float]:
    """Fingerprint of freeze-cache hydro fields (Inner NR mapping ``h``)."""
    vals: list[float] = []
    for cell in cells:
        cache = getattr(cell, "_vorab_hydro_cache", None) or {}
        vals.extend(
            _hydro_field_vector_from_mapping(
                u_mf=float(cache.get("u_mf", 0.0)),
                u_b=float(cache.get("u_b", 0.0)),
                d_b=float(cache.get("d_b", 0.0)),
                eps_b=float(cache.get("eps_b", 0.0)),
                eps_d=float(cache.get("eps_d", 0.0)),
                K_bd=float(cache.get("K_bd", 0.0)),
                V_b=float(cache.get("V_b", 0.0)),
                V_d=float(cache.get("V_d", 0.0)),
                u0=float(cache.get("u0", 0.0)),
                u_d=float(cache.get("u_d", 0.0)),
                K_solid_auf_abs=float(np.sum(np.abs(getattr(cell, "_vorab_K_solid_auf", 0.0)))),
                K_solid_ab_abs=float(np.sum(np.abs(getattr(cell, "_vorab_K_solid_ab", 0.0)))),
            )
        )
    arr = np.asarray(vals, dtype=np.float64)
    return {
        "n": float(arr.size),
        "sum": float(np.sum(arr)) if arr.size else 0.0,
        "sumsq": float(np.dot(arr, arr)) if arr.size else 0.0,
        "valid_count": float(
            sum(1 for c in cells if bool(getattr(c, "_vorab_hydro_cache_valid", False)))
        ),
    }


def live_hydrodynamics_fingerprint(cells: List[Cell]) -> dict[str, float]:
    """Fingerprint of live hydro fields currently on the cells."""
    vals: list[float] = []
    for cell in cells:
        vals.extend(
            _hydro_field_vector_from_mapping(
                u_mf=float(cell.u_mf),
                u_b=float(cell.u_b),
                d_b=float(cell.d_b),
                eps_b=float(cell.eps_b),
                eps_d=float(cell.eps_d),
                K_bd=float(cell.K_bd),
                V_b=float(cell.V_b),
                V_d=float(cell.V_d),
                u0=float(cell.u0),
                u_d=float(cell.u_d),
                K_solid_auf_abs=float(np.sum(np.abs(cell.K_solid_auf))),
                K_solid_ab_abs=float(np.sum(np.abs(cell.K_solid_ab))),
            )
        )
    arr = np.asarray(vals, dtype=np.float64)
    return {
        "n": float(arr.size),
        "sum": float(np.sum(arr)) if arr.size else 0.0,
        "sumsq": float(np.dot(arr, arr)) if arr.size else 0.0,
    }


def prepare_inner_nr_hydrodynamics_freeze(
    cells: List[Cell],
    *,
    cfg: Any | None = None,
) -> dict[str, float]:
    """Recompute+snapshot hydro at Inner entry so freeze matches Pre-Inner state.

    P1 Hamel ``outer_fixed`` semantics: hydro mapping ``h`` is fixed for the whole
    Inner NR solve. Updates happen at outer refresh and here at Inner entry (after
    Pre-Inner). Mid-Inner accept-step refresh is opt-in only
    (``nr_refresh_hydrodynamics_on_accepted_step_thesis``).

    Outer refresh snapshots ``h`` before Pre-Inner. Pre-Inner may mutate N/T/holdup
    and gate on live hydro residuals. Without this prepare step, Inner would solve
    ``F(x; h_stale)`` under the pre-PreInner freeze cache.

    Returns pre/post freeze fingerprints for diagnostics.
    """
    fp_before = hydrodynamics_freeze_fingerprint(cells)
    refresh_cells_hydrodynamics_for_frozen_inner(cells, cfg=cfg)
    fp_after = hydrodynamics_freeze_fingerprint(cells)
    live_fp = live_hydrodynamics_fingerprint(cells)
    return {
        "freeze_sum_before": float(fp_before["sum"]),
        "freeze_sum_after": float(fp_after["sum"]),
        "live_sum_after": float(live_fp["sum"]),
        "freeze_shifted": float(
            abs(float(fp_after["sum"]) - float(fp_before["sum"])) > 1e-12
            or abs(float(fp_after["sumsq"]) - float(fp_before["sumsq"])) > 1e-12
        ),
        "freeze_matches_live": float(
            abs(float(fp_after["sum"]) - float(live_fp["sum"])) <= 1e-9
            and abs(float(fp_after["sumsq"]) - float(live_fp["sumsq"])) <= 1e-9
        ),
    }


def initialize_fixed_vorabrechnung_sources_for_cells(
    cells: List[Cell],
    *,
    T_reference: float,
    per_cell_temperature: bool = False,
    vm_zone_cell_indices: list[int] | None = None,
) -> None:
    """Initialize fixed drying/pyrolysis Vorabrechnung sources at a reference temperature.

    Hamel-style strict mode: drying/DAEM are budgeted once (Vorabrechnung) and then
    treated as fixed sources for subsequent outer/inner iterations.

    ``per_cell_temperature=True`` 时各 cell 用自身 ``T`` 计算源项（VM 脱挥发分区）。
    """
    T_ref = float(max(T_reference, 1.0))
    zone_set = set(vm_zone_cell_indices or [])
    for i, cell in enumerate(cells):
        T_saved = float(cell.T)
        use_local_t = bool(per_cell_temperature) or (zone_set and i in zone_set)
        cell.T = float(T_saved if use_local_t else T_ref)
        cell.invalidate_vorabrechnung_cache()
        refresh_cell_vorabrechnung(cell, force=False)
        cell.T = T_saved
        refresh_cell_hydrodynamics_for_frozen_inner(cell)


def refresh_vorabrechnung_for_cells(
    cells: List[Cell],
    *,
    force: bool = False,
    refresh_sources: bool = True,
) -> None:
    """批量刷新 Vorabrechnung（outer Abgleich 调用入口）。

    Parameters
    ----------
    refresh_sources:
        True: refresh hydrodynamics + drying/pyrolysis sources.
        False: refresh hydrodynamics/snapshots only, keep drying/pyrolysis fixed.
    """
    if not refresh_sources:
        refresh_cells_hydrodynamics_for_frozen_inner(cells)
        return
    for cell in cells:
        refresh_cell_vorabrechnung(cell, force=force)
    # After source refresh, recompute bed K from the new hydro and resnapshot so
    # the next inner freeze does not carry stale transport coefficients.
    _update_bed_solid_transport_after_hydro_refresh(cells)
    for cell in cells:
        cell._snapshot_vorabrechnung_hydrodynamics()


