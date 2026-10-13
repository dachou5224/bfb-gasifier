"""床顶 Eq.2.7 固体焓出口透传。"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import GAS_SPECIES, N_GAS
from src.solvers.global_nr_solver import (
    arm_bed_top_solid_energy_passthrough_if_y_co,
    install_bed_top_solid_energy_passthrough_on_cells,
)


def _bare_bed() -> Cell:
    cell = Cell.__new__(Cell)
    nk = 2
    cell.solid_state_model = "holdup_transport"
    cell.m_solid = np.ones((nk, N_SOLID_COMP))
    cell.K_solid_auf = np.full((nk, N_SOLID_COMP), 0.5)
    cell.K_solid_ab = np.full((nk, N_SOLID_COMP), 0.1)
    cell.m_solid_in = np.zeros((nk, N_SOLID_COMP))
    cell.m_solid_auf_in = np.full((nk, N_SOLID_COMP), 0.2)
    cell.m_solid_ab_in = np.full((nk, N_SOLID_COMP), 0.05)
    cell.m_solid_zu = np.zeros((nk, N_SOLID_COMP))
    cell.m_solid_rez = np.zeros((nk, N_SOLID_COMP))
    cell._solid_energy_passthrough = False
    return cell


def test_passthrough_uses_inlet_mass_not_Km():
    cell = _bare_bed()
    k_out = cell._solid_outflow_rates()
    assert np.allclose(cell._solid_energy_outlet_rates(), k_out)
    cell._solid_energy_passthrough = True
    e_out = cell._solid_energy_outlet_rates()
    expect = cell.m_solid_auf_in + cell.m_solid_ab_in
    assert np.allclose(e_out, expect)
    assert not np.allclose(e_out, k_out)


def test_freeboard_closure_still_passthrough():
    cell = _bare_bed()
    cell.solid_state_model = "freeboard_closure"
    cell._solid_energy_passthrough = False
    e_out = cell._solid_energy_outlet_rates()
    expect = cell.m_solid_auf_in + cell.m_solid_ab_in
    assert np.allclose(e_out, expect)


def _fake_beds(n: int = 3) -> list[Cell]:
    cells: list[Cell] = []
    for _ in range(n):
        c = Cell.__new__(Cell)
        c.cell_type = "bed"
        c.N_d = np.zeros(N_GAS)
        c.N_b = np.zeros(N_GAS)
        c._solid_energy_passthrough = False
        cells.append(c)
    return cells


def test_install_immediate_arms_top_only():
    cells = _fake_beds(3)
    cfg = SimpleNamespace(
        nr_bed_top_solid_energy_passthrough_thesis=True,
        nr_bed_top_solid_energy_passthrough_after_y_co_thesis=False,
    )
    install_bed_top_solid_energy_passthrough_on_cells(cells, cfg)
    assert cells[0]._solid_energy_passthrough is False
    assert cells[1]._solid_energy_passthrough is False
    assert cells[2]._solid_energy_passthrough is True


def test_deferred_arms_after_y_co():
    cells = _fake_beds(3)
    cfg = SimpleNamespace(
        nr_bed_top_solid_energy_passthrough_thesis=True,
        nr_bed_top_solid_energy_passthrough_after_y_co_thesis=True,
        nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis=0.14,
    )
    install_bed_top_solid_energy_passthrough_on_cells(cells, cfg)
    assert cells[2]._solid_energy_passthrough is False
    assert cells[2]._solid_energy_passthrough_deferred is True
    j_co = GAS_SPECIES.index("CO")
    j_n2 = GAS_SPECIES.index("N2")
    cells[2].N_d[j_co] = 0.10
    cells[2].N_d[j_n2] = 0.90  # y_CO=0.10
    r1 = arm_bed_top_solid_energy_passthrough_if_y_co(cells, cfg)
    assert r1["armed"] is False
    cells[2].N_d[j_co] = 0.20
    cells[2].N_d[j_n2] = 0.80  # y_CO=0.20
    r2 = arm_bed_top_solid_energy_passthrough_if_y_co(cells, cfg)
    assert r2["newly_armed"] is True
    assert cells[2]._solid_energy_passthrough is True
    assert cells[0]._solid_energy_passthrough is False
