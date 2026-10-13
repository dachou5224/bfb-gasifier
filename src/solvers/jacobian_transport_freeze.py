"""Freeze ``K_solid_auf`` / ``K_solid_ab`` during global-NR Jacobian column FD.

Ref: Hamel (1999) Eq. (2-6) holdup transport. Transport coefficients are not packed
NR unknowns but must stay fixed at the linearization point while ``m`` columns are
perturbed; ``m_solid_auf_in = K_auf * m`` may still refresh from perturbed holdup.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from src.core.cell import Cell


def enter_jacobian_transport_k_freeze(cells: Iterable[Cell]) -> None:
    """Snapshot current transport coefficients on all solver cells."""
    for cell in cells:
        cell._jacobian_frozen_transport_k = True
        cell._jacobian_K_solid_auf_snapshot = np.array(cell.K_solid_auf, dtype=np.float64, copy=True)
        cell._jacobian_K_solid_ab_snapshot = np.array(cell.K_solid_ab, dtype=np.float64, copy=True)


def exit_jacobian_transport_k_freeze(cells: Iterable[Cell]) -> None:
    """Clear Jacobian transport freeze flags and snapshots."""
    for cell in cells:
        cell._jacobian_frozen_transport_k = False
        cell._jacobian_K_solid_auf_snapshot = None
        cell._jacobian_K_solid_ab_snapshot = None


def restore_jacobian_frozen_transport_k(cells: Iterable[Cell]) -> bool:
    """Replay snapshotted ``K_solid`` if Jacobian freeze is active."""
    restored = False
    for cell in cells:
        if not bool(getattr(cell, "_jacobian_frozen_transport_k", False)):
            continue
        snap_auf = getattr(cell, "_jacobian_K_solid_auf_snapshot", None)
        snap_ab = getattr(cell, "_jacobian_K_solid_ab_snapshot", None)
        if snap_auf is not None:
            cell.K_solid_auf[:, :] = snap_auf
        if snap_ab is not None:
            cell.K_solid_ab[:, :] = snap_ab
        restored = True
    return restored


def jacobian_transport_k_is_frozen(cells: Iterable[Cell]) -> bool:
    return any(bool(getattr(cell, "_jacobian_frozen_transport_k", False)) for cell in cells)
