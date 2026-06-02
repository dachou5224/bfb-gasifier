"""Legacy Gauss–Seidel 诊断与软回滚辅助（与 NR 主路径隔离）。

仅被 ``Reactor._solve_gauss_seidel`` 等遗留路径使用；全局 NR 不应依赖本模块。
"""

from __future__ import annotations

import numpy as np

from src.core.cell import Cell, S_CHAR, S_MOISTURE, S_VM
from src.core.species import GAS_SPECIES_INDEX


def _positive_rebound_penalty(values: np.ndarray, scale: float) -> float:
    """Return normalized positive-rebound amount along an axial profile."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.size <= 1:
        return 0.0
    diffs = np.diff(arr)
    rebound = float(np.sum(np.maximum(diffs, 0.0)))
    return rebound / max(float(scale), 1e-12)


def _gs_profile_metrics(
    cells: list[Cell],
    o2_feed: float,
    T_ref: np.ndarray | None = None,
    m_char_in: float | None = None,
) -> dict[str, float]:
    """Build GS profile-quality metrics for soft rollback decisions."""
    idx_o2 = GAS_SPECIES_INDEX["O2"]
    o2_prof = np.array([max(float(c.N_d[idx_o2] + c.N_b[idx_o2]), 0.0) for c in cells], dtype=np.float64)
    vm_prof = np.array([max(float(np.sum(c.m_solid[:, S_VM])), 0.0) for c in cells], dtype=np.float64)
    moist_prof = np.array([max(float(np.sum(c.m_solid[:, S_MOISTURE])), 0.0) for c in cells], dtype=np.float64)

    o2_rebound = _positive_rebound_penalty(o2_prof, scale=max(float(o2_feed), 1e-12))
    vm_rebound = _positive_rebound_penalty(vm_prof, scale=max(float(vm_prof[0]) if vm_prof.size > 0 else 0.0, 1e-12))
    moist_rebound = _positive_rebound_penalty(
        moist_prof, scale=max(float(moist_prof[0]) if moist_prof.size > 0 else 0.0, 1e-12)
    )
    o2_slip = float(o2_prof[-1] / max(float(o2_feed), 1e-12)) if o2_prof.size > 0 else 0.0

    runaway_penalty = 0.0
    if T_ref is not None and m_char_in is not None and len(cells) > 0:
        T_ref_arr = np.asarray(T_ref, dtype=np.float64)
        T_cur = np.array([float(c.T) for c in cells], dtype=np.float64)
        hot_excess = np.maximum(T_cur - T_ref_arr[: len(cells)], 0.0)
        hot_rel = float(np.mean(hot_excess / np.maximum(T_ref_arr[: len(cells)], 1.0)))
        m_char_out = float(np.sum(np.maximum(cells[-1].m_solid[:, S_CHAR], 0.0)))
        x_char = float(np.clip(1.0 - m_char_out / max(float(m_char_in), 1e-12), 0.0, 1.0))
        runaway_penalty = hot_rel * max(x_char - 0.90, 0.0)

    penalty = 2.0 * o2_slip + o2_rebound + vm_rebound + moist_rebound + runaway_penalty
    return {
        "o2_rebound": float(o2_rebound),
        "vm_rebound": float(vm_rebound),
        "moist_rebound": float(moist_rebound),
        "o2_slip": float(o2_slip),
        "runaway_penalty": float(runaway_penalty),
        "penalty": float(penalty),
    }


def _clamp_soft_rollback_monotonic(cell: Cell) -> None:
    """Clamp O2/VM/moisture to inlet budgets to avoid upward rebound."""
    idx_o2 = GAS_SPECIES_INDEX["O2"]
    o2_cap = max(float(cell.N_d_in[idx_o2] + cell.N_b_in[idx_o2]), 0.0)
    o2_cur = max(float(cell.N_d[idx_o2] + cell.N_b[idx_o2]), 0.0)
    if o2_cur > o2_cap + 1e-12:
        fac = o2_cap / max(o2_cur, 1e-12)
        cell.N_d[idx_o2] *= fac
        cell.N_b[idx_o2] *= fac

    for comp_idx in (S_VM, S_MOISTURE):
        cap = max(float(np.sum(cell.m_solid_in[:, comp_idx] + cell.m_solid_zu[:, comp_idx])), 0.0)
        cur = max(float(np.sum(cell.m_solid[:, comp_idx])), 0.0)
        if cur > cap + 1e-12:
            cell.m_solid[:, comp_idx] *= cap / max(cur, 1e-12)


def _capture_cells_state(cells: list[Cell]) -> list[dict[str, np.ndarray | float]]:
    """Capture a lightweight mutable state snapshot for rollback."""
    snapshots: list[dict[str, np.ndarray | float]] = []
    for c in cells:
        snapshots.append(
            {
                "N_d": np.array(c.N_d, copy=True),
                "N_b": np.array(c.N_b, copy=True),
                "m_solid": np.array(c.m_solid, copy=True),
                "T": float(c.T),
            }
        )
    return snapshots


def _restore_cells_state(cells: list[Cell], snapshots: list[dict[str, np.ndarray | float]]) -> None:
    """Restore mutable fields from snapshots."""
    for c, snap in zip(cells, snapshots):
        c.N_d[:] = np.asarray(snap["N_d"], dtype=np.float64)
        c.N_b[:] = np.asarray(snap["N_b"], dtype=np.float64)
        c.m_solid[:, :] = np.asarray(snap["m_solid"], dtype=np.float64)
        c.T = float(snap["T"])
