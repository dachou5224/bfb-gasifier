"""phase_share 总量吸收：按相份额分给两相。"""

from __future__ import annotations

import numpy as np

from src.core.cell import Cell
from src.core.connectivity.preprojection import oneshot_bed_dense_species_sum_absorb
from src.core.species import GAS_SPECIES_INDEX


def test_phase_share_absorb_preserves_phase_fraction():
    cell = Cell.__new__(Cell)
    n = len(GAS_SPECIES_INDEX)
    cell.N_d = np.zeros(n, dtype=np.float64)
    cell.N_b = np.zeros(n, dtype=np.float64)
    cell.N_zu_d = np.zeros(n, dtype=np.float64)
    cell.N_zu_b = np.zeros(n, dtype=np.float64)
    cell.N_rez_d = np.zeros(n, dtype=np.float64)
    cell.N_rez_b = np.zeros(n, dtype=np.float64)
    cell.N_d_in = np.zeros(n, dtype=np.float64)
    cell.N_b_in = np.zeros(n, dtype=np.float64)
    cell.R_gas_d = np.zeros(n, dtype=np.float64)
    cell.R_gas_b = np.zeros(n, dtype=np.float64)
    cell._thermo_cache_valid = True
    cell._hydro_cache_valid = True

    j = int(GAS_SPECIES_INDEX["CO"])
    cell.N_d[j] = 3.0
    cell.N_b[j] = 1.0
    # res_sum = 0+0+4+0 - 4 = +4? set R so res_sum=+4
    cell.R_gas_d[j] = 4.0
    # in+R - N = 4 - 4 = 0; bump R
    cell.R_gas_d[j] = 8.0  # res_sum = 8 - 4 = 4

    f0 = 3.0 / 4.0
    added = oneshot_bed_dense_species_sum_absorb(
        cell, species=("CO",), fraction=1.0, mode="phase_share"
    )
    assert added["CO"] == 4.0
    tot = float(cell.N_d[j] + cell.N_b[j])
    assert abs(tot - 8.0) < 1e-9
    assert abs(float(cell.N_d[j]) / tot - f0) < 1e-9
