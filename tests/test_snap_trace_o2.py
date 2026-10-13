"""微量 O₂ holdup 粘滞清零。"""

from __future__ import annotations

import numpy as np

from src.core.cell import Cell
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.global_nr_solver import install_trace_o2_snap_on_cells, snap_trace_o2_holdups


def _bare_cell() -> Cell:
    cell = Cell.__new__(Cell)
    cell.N_d = np.zeros(N_GAS, dtype=np.float64)
    cell.N_b = np.zeros(N_GAS, dtype=np.float64)
    cell.N_d_in = np.zeros(N_GAS, dtype=np.float64)
    cell.N_b_in = np.zeros(N_GAS, dtype=np.float64)
    cell.N_zu_d = np.zeros(N_GAS, dtype=np.float64)
    cell.N_zu_b = np.zeros(N_GAS, dtype=np.float64)
    return cell


def test_snap_trace_o2_clears_dust_without_inlet():
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-16
    cell.N_b[j] = 2.0e-16
    install_trace_o2_snap_on_cells(
        [cell],
        type("Cfg", (), {"nr_snap_trace_o2_holdup_thesis": True, "nr_snap_trace_o2_atol_mol_s_thesis": 1e-12})(),
    )
    n = snap_trace_o2_holdups([cell])
    assert n == 2
    assert cell.N_d[j] == 0.0
    assert cell.N_b[j] == 0.0


def test_snap_trace_o2_keeps_o2_when_inlet_present():
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-16
    cell.N_d_in[j] = 1.0e-6
    install_trace_o2_snap_on_cells(
        [cell],
        type("Cfg", (), {"nr_snap_trace_o2_holdup_thesis": True, "nr_snap_trace_o2_atol_mol_s_thesis": 1e-12})(),
    )
    n = snap_trace_o2_holdups([cell])
    assert n == 0
    assert cell.N_d[j] == 1.0e-16


def test_snap_trace_o2_disabled_by_default():
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-16
    install_trace_o2_snap_on_cells([cell], type("Cfg", (), {})())
    assert snap_trace_o2_holdups([cell]) == 0
    assert cell.N_d[j] == 1.0e-16


def test_snap_trace_o2_min_bed_index_skips_bed0():
    cells = [_bare_cell(), _bare_cell()]
    for c in cells:
        c.cell_type = "bed"
    j = int(GAS_SPECIES_INDEX["O2"])
    cells[0].N_d[j] = 1.0e-16
    cells[1].N_d[j] = 1.0e-16
    install_trace_o2_snap_on_cells(
        cells,
        type(
            "Cfg",
            (),
            {
                "nr_snap_trace_o2_holdup_thesis": True,
                "nr_snap_trace_o2_atol_mol_s_thesis": 1e-12,
                "nr_snap_trace_o2_min_bed_index_thesis": 1,
            },
        )(),
    )
    assert cells[0]._nr_snap_trace_o2 is False
    assert cells[1]._nr_snap_trace_o2 is True
    n = snap_trace_o2_holdups(cells)
    assert n == 1
    assert cells[0].N_d[j] == 1.0e-16
    assert cells[1].N_d[j] == 0.0
