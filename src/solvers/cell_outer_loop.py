"""Cell-level preparation stage for GS outer loop.

This isolates Zellenmodell pre-step logic from Reactor orchestration:
- inlet-based prefill
- relaxed warm-start blend
- Vorabrechnung refresh
- light pre-reaction nudging
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.core.cell import Cell
from src.solvers.vorabrechnung import refresh_cell_vorabrechnung


@dataclass
class CellOuterPrepState:
    seed_relax: float
    pre_relax: float


def prepare_cell_for_gs_outer_iteration(
    *,
    cell: Cell,
    g_iter: int,
    T_old: float,
    N_d_old: np.ndarray,
    N_b_old: np.ndarray,
    m_old: np.ndarray,
    T_est: float,
) -> CellOuterPrepState:
    """Prepare one cell state before `solve_cell(...)` in GS path."""
    N_d_prefill = np.maximum(cell.N_d_in + cell.N_zu_d + cell.N_rez_d, 0.0)
    N_b_prefill = np.maximum(cell.N_b_in + cell.N_zu_b + cell.N_rez_b, 0.0)
    m_prefill = np.maximum(cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez, 0.0)

    if g_iter >= 1:
        seed_relax = 0.55
        cell.N_d[:] = np.maximum(seed_relax * N_d_old + (1.0 - seed_relax) * N_d_prefill, 0.0)
        cell.N_b[:] = np.maximum(seed_relax * N_b_old + (1.0 - seed_relax) * N_b_prefill, 0.0)
        cell.m_solid[:] = np.maximum(seed_relax * m_old + (1.0 - seed_relax) * m_prefill, 0.0)
        cell.T = seed_relax * T_old + (1.0 - seed_relax) * float(T_est)
    else:
        seed_relax = 0.0
        cell.N_d[:] = N_d_prefill
        cell.N_b[:] = N_b_prefill
        cell.m_solid[:] = m_prefill

    refresh_cell_vorabrechnung(cell, force=False)
    cell.calc_reactions()

    pre_relax = 0.15 if g_iter == 0 else 0.05
    cell.N_d[:] = np.maximum(
        cell.N_d_in + cell.N_zu_d + cell.N_rez_d + pre_relax * cell.R_gas_d,
        0.0,
    )
    cell.N_b[:] = np.maximum(
        cell.N_b_in + cell.N_zu_b + cell.N_rez_b + pre_relax * cell.R_gas_b,
        0.0,
    )
    cell.m_solid[:] = np.maximum(
        cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez + pre_relax * cell.R_solid,
        0.0,
    )

    return CellOuterPrepState(seed_relax=float(seed_relax), pre_relax=float(pre_relax))

