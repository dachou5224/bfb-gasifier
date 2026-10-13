"""物种总量双向闭合（负 res_sum 放出；§5.25）。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import Cell
from src.core.connectivity.preprojection import oneshot_bed_dense_species_sum_absorb
from src.core.species import GAS_SPECIES, N_GAS


def _cell_n2(*, n_hold: float, n_in: float) -> Cell:
    c = Cell.__new__(Cell)
    c.N_d = np.zeros(N_GAS)
    c.N_b = np.zeros(N_GAS)
    c.N_zu_d = np.zeros(N_GAS)
    c.N_zu_b = np.zeros(N_GAS)
    c.N_rez_d = np.zeros(N_GAS)
    c.N_rez_b = np.zeros(N_GAS)
    c.N_d_in = np.zeros(N_GAS)
    c.N_b_in = np.zeros(N_GAS)
    c.R_gas_d = np.zeros(N_GAS)
    c.R_gas_b = np.zeros(N_GAS)
    j = GAS_SPECIES.index("N2")
    c.N_d[j] = 0.6 * n_hold
    c.N_b[j] = 0.4 * n_hold
    c.N_d_in[j] = 0.6 * n_in
    c.N_b_in[j] = 0.4 * n_in
    return c


def test_default_ignores_negative_res_sum():
    c = _cell_n2(n_hold=20.0, n_in=10.0)
    j = GAS_SPECIES.index("N2")
    before = float(c.N_d[j] + c.N_b[j])
    added = oneshot_bed_dense_species_sum_absorb(c, species=("N2",), mode="phase_share")
    assert added["N2"] == 0.0
    assert float(c.N_d[j] + c.N_b[j]) == pytest.approx(before)


def test_bidirectional_reduces_overfilled_holdup():
    c = _cell_n2(n_hold=20.0, n_in=10.0)
    j = GAS_SPECIES.index("N2")
    added = oneshot_bed_dense_species_sum_absorb(
        c, species=("N2",), mode="phase_share", bidirectional=True, fraction=1.0
    )
    assert added["N2"] == pytest.approx(-10.0)
    assert float(c.N_d[j] + c.N_b[j]) == pytest.approx(10.0)
    # phase share preserved
    assert float(c.N_d[j] / (c.N_d[j] + c.N_b[j])) == pytest.approx(0.6)
