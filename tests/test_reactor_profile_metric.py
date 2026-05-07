from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import Cell, CellGeometry, SolidProps, S_MOISTURE, S_VM
from src.core.reactor import _clamp_soft_rollback_monotonic, _gs_profile_metrics, _positive_rebound_penalty
from src.core.species import GAS_SPECIES_INDEX


def _make_cell() -> Cell:
    geo = CellGeometry(D_bed=0.6, dh=0.5)
    solid = SolidProps(
        rho_s=1400.0,
        d_p=2.25e-3,
        n_size_classes=1,
        d_p_classes=np.array([2.25e-3]),
        mass_fractions=np.array([1.0]),
        phi_s=0.86,
        eps_mf=0.45,
    )
    cell = Cell(geo=geo, solid=solid)
    cell.T = 1100.0
    cell.P = 2.5e6
    return cell


def test_positive_rebound_penalty_zero_for_monotone_decrease():
    values = np.array([1.0, 0.8, 0.3, 0.1], dtype=float)
    assert _positive_rebound_penalty(values, scale=1.0) == 0.0


def test_positive_rebound_penalty_counts_only_upward_steps():
    values = np.array([1.0, 0.7, 0.9, 0.4, 0.5], dtype=float)
    assert _positive_rebound_penalty(values, scale=1.0) == pytest.approx(0.3)


def test_gs_profile_metrics_penalize_o2_slip_and_solid_rebound():
    cells = [_make_cell() for _ in range(3)]
    idx_o2 = GAS_SPECIES_INDEX["O2"]

    cells[0].N_d[idx_o2] = 1.0
    cells[1].N_d[idx_o2] = 0.2
    cells[2].N_d[idx_o2] = 0.4

    cells[0].m_solid[0, S_VM] = 0.10
    cells[1].m_solid[0, S_VM] = 0.02
    cells[2].m_solid[0, S_VM] = 0.03

    cells[0].m_solid[0, S_MOISTURE] = 0.05
    cells[1].m_solid[0, S_MOISTURE] = 0.01
    cells[2].m_solid[0, S_MOISTURE] = 0.015

    metrics = _gs_profile_metrics(cells, o2_feed=1.0)

    assert metrics["o2_rebound"] == pytest.approx(0.2)
    assert metrics["vm_rebound"] == pytest.approx(0.1)
    assert metrics["moist_rebound"] == pytest.approx(0.1)
    assert metrics["o2_slip"] == pytest.approx(0.4)
    assert metrics["runaway_penalty"] == pytest.approx(0.0)
    assert metrics["penalty"] == pytest.approx(1.2)


def test_gs_profile_metrics_penalize_hot_high_conversion_runaway():
    cells = [_make_cell() for _ in range(3)]
    for cell in cells:
        cell.m_solid[0, S_VM] = 0.0
        cell.m_solid[0, S_MOISTURE] = 0.0
    cells[0].m_solid[0, 0] = 0.04
    cells[1].m_solid[0, 0] = 0.01
    cells[2].m_solid[0, 0] = 0.001
    cells[0].T = 1380.0
    cells[1].T = 1340.0
    cells[2].T = 1310.0

    metrics = _gs_profile_metrics(
        cells,
        o2_feed=1.0,
        T_ref=np.array([1260.0, 1210.0, 1130.0]),
        m_char_in=0.04,
    )

    assert metrics["runaway_penalty"] > 0.0
    assert metrics["penalty"] >= metrics["runaway_penalty"]


def test_clamp_soft_rollback_monotonic_limits_o2_vm_and_moisture_to_inlet():
    cell = _make_cell()
    idx_o2 = GAS_SPECIES_INDEX["O2"]

    cell.N_d_in[idx_o2] = 0.2
    cell.N_b_in[idx_o2] = 0.1
    cell.N_d[idx_o2] = 0.5
    cell.N_b[idx_o2] = 0.4

    cell.m_solid_in[0, S_VM] = 0.05
    cell.m_solid_zu[0, S_VM] = 0.01
    cell.m_solid[0, S_VM] = 0.20

    cell.m_solid_in[0, S_MOISTURE] = 0.03
    cell.m_solid_zu[0, S_MOISTURE] = 0.01
    cell.m_solid[0, S_MOISTURE] = 0.12

    _clamp_soft_rollback_monotonic(cell)

    assert cell.N_d[idx_o2] + cell.N_b[idx_o2] == pytest.approx(0.3)
    assert np.sum(cell.m_solid[:, S_VM]) == pytest.approx(0.06)
    assert np.sum(cell.m_solid[:, S_MOISTURE]) == pytest.approx(0.04)
