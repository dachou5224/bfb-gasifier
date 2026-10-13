"""微量 O₂ 速率门控：不动 holdup，仅清零速率评价用的 C/y。"""

from __future__ import annotations

import numpy as np

from src.core.cell import Cell
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.global_nr_solver import (
    gate_trace_o2_for_rate_eval,
    install_trace_o2_rate_gate_on_cells,
)


def _bare_cell() -> Cell:
    cell = Cell.__new__(Cell)
    cell.N_d = np.zeros(N_GAS, dtype=np.float64)
    cell.N_b = np.zeros(N_GAS, dtype=np.float64)
    cell.N_d_in = np.zeros(N_GAS, dtype=np.float64)
    cell.N_b_in = np.zeros(N_GAS, dtype=np.float64)
    cell.N_zu_d = np.zeros(N_GAS, dtype=np.float64)
    cell.N_zu_b = np.zeros(N_GAS, dtype=np.float64)
    cell.cell_type = "bed"
    return cell


def test_rate_gate_zeros_c_y_keeps_holdup():
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-16
    cell.N_b[j] = 2.0e-16
    install_trace_o2_rate_gate_on_cells(
        [cell],
        type(
            "Cfg",
            (),
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_atol_mol_s_thesis": 1e-12,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 0,
            },
        )(),
    )
    C_b = np.ones(N_GAS)
    C_d = np.ones(N_GAS)
    y_b = np.ones(N_GAS) * 0.1
    y_d = np.ones(N_GAS) * 0.1
    C_b[j] = 1e-9
    C_d[j] = 1e-9
    y_b[j] = 1e-12
    y_d[j] = 1e-12
    Cb2, Cd2, yb2, yd2 = gate_trace_o2_for_rate_eval(cell, C_b, C_d, y_b, y_d)
    assert Cb2[j] == 0.0 and Cd2[j] == 0.0
    assert yb2[j] == 0.0 and yd2[j] == 0.0
    assert cell.N_d[j] == 1.0e-16
    assert cell.N_b[j] == 2.0e-16
    # 原数组不被原地改写
    assert C_b[j] == 1e-9


def test_rate_gate_ignores_inlet_uses_local_holdup():
    """上游 N_in 不放行；局地灰尘仍门控（bed2 典型）。"""
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-16
    cell.N_d_in[j] = 1.0e-6
    install_trace_o2_rate_gate_on_cells(
        [cell],
        type(
            "Cfg",
            (),
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 0,
            },
        )(),
    )
    C = np.ones(N_GAS)
    C[j] = 1e-9
    y = np.ones(N_GAS) * 0.1
    y[j] = 1e-12
    Cb2, _, _, _ = gate_trace_o2_for_rate_eval(cell, C, C, y, y)
    assert Cb2[j] == 0.0


def test_rate_gate_keeps_real_local_o2():
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-3  # [mol/s] 真实局地库存
    install_trace_o2_rate_gate_on_cells(
        [cell],
        type(
            "Cfg",
            (),
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 0,
            },
        )(),
    )
    C = np.ones(N_GAS)
    C[j] = 1.0  # [mol/m³]
    y = np.ones(N_GAS) * 0.1
    y[j] = 0.05
    Cb2, _, _, _ = gate_trace_o2_for_rate_eval(cell, C, C, y, y)
    assert Cb2[j] == 1.0


def test_rate_gate_default_off():
    cell = _bare_cell()
    j = int(GAS_SPECIES_INDEX["O2"])
    cell.N_d[j] = 1.0e-16
    install_trace_o2_rate_gate_on_cells([cell], type("Cfg", (), {})())
    C = np.ones(N_GAS)
    C[j] = 1e-9
    y = np.ones(N_GAS)
    Cb2, _, _, _ = gate_trace_o2_for_rate_eval(cell, C, C, y, y)
    assert Cb2[j] == 1e-9


def test_rate_gate_min_bed_index():
    cells = [_bare_cell(), _bare_cell(), _bare_cell()]
    install_trace_o2_rate_gate_on_cells(
        cells,
        type(
            "Cfg",
            (),
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 2,
            },
        )(),
    )
    assert cells[0]._nr_trace_o2_rate_gate is False
    assert cells[1]._nr_trace_o2_rate_gate is False
    assert cells[2]._nr_trace_o2_rate_gate is True


def test_rate_gate_deferred_arms_after_y_co():
    from src.solvers.global_nr_solver import arm_trace_o2_rate_gate_if_y_co

    cells = [_bare_cell(), _bare_cell(), _bare_cell()]
    j_co = int(GAS_SPECIES_INDEX["CO"])
    j_n2 = int(GAS_SPECIES_INDEX["N2"])
    # 床顶 CO 先低后高
    for c in cells:
        c.N_d[j_n2] = 1.0
        c.N_b[j_n2] = 0.0
    cells[-1].N_d[j_co] = 0.05  # y_CO≈0.05
    cfg = type(
        "Cfg",
        (),
        {
            "nr_trace_o2_rate_gate_thesis": True,
            "nr_trace_o2_rate_gate_min_bed_index_thesis": 2,
            "nr_trace_o2_rate_gate_after_y_co_thesis": True,
            "nr_trace_o2_rate_gate_y_co_threshold_thesis": 0.14,
        },
    )()
    install_trace_o2_rate_gate_on_cells(cells, cfg)
    assert cells[2]._nr_trace_o2_rate_gate is False
    r1 = arm_trace_o2_rate_gate_if_y_co(cells, cfg)
    assert r1["armed"] is False
    cells[-1].N_d[j_co] = 0.20  # y_CO≈0.167
    r2 = arm_trace_o2_rate_gate_if_y_co(cells, cfg)
    assert r2["newly_armed"] is True
    assert cells[2]._nr_trace_o2_rate_gate is True
    assert cells[0]._nr_trace_o2_rate_gate is False
