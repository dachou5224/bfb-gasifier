"""上段密相合成气总量吸收 oneshot 单测。"""

from __future__ import annotations

import numpy as np

from src.core.connectivity.preprojection import oneshot_bed_dense_species_sum_absorb
from src.core.species import GAS_SPECIES_INDEX, N_GAS


class _FakeCell:
    def __init__(self) -> None:
        self.N_zu_d = np.zeros(N_GAS)
        self.N_zu_b = np.zeros(N_GAS)
        self.N_rez_d = np.zeros(N_GAS)
        self.N_rez_b = np.zeros(N_GAS)
        self.N_d_in = np.zeros(N_GAS)
        self.N_b_in = np.zeros(N_GAS)
        self.R_gas_d = np.zeros(N_GAS)
        self.R_gas_b = np.zeros(N_GAS)
        self.N_d = np.zeros(N_GAS)
        self.N_b = np.zeros(N_GAS)
        self._thermo_cache_valid = True
        self._hydro_cache_valid = True


def test_oneshot_absorbs_positive_res_sum_into_dense_only() -> None:
    cell = _FakeCell()
    j = GAS_SPECIES_INDEX["CO"]
    cell.R_gas_d[j] = 5.0
    cell.N_d[j] = 1.0
    cell.N_b[j] = 0.5
    # res_sum = 5 - 1.5 = 3.5
    added = oneshot_bed_dense_species_sum_absorb(cell, species=("CO", "H2"))
    assert added["CO"] == 3.5
    assert cell.N_d[j] == 4.5
    assert cell.N_b[j] == 0.5
    assert cell._thermo_cache_valid is False


def test_oneshot_skips_nonpositive_res_sum() -> None:
    cell = _FakeCell()
    j = GAS_SPECIES_INDEX["CO"]
    cell.N_d[j] = 10.0
    cell.R_gas_d[j] = 1.0
    added = oneshot_bed_dense_species_sum_absorb(cell, species=("CO",))
    assert added["CO"] == 0.0
    assert cell.N_d[j] == 10.0


def test_oneshot_respects_fraction() -> None:
    cell = _FakeCell()
    j = GAS_SPECIES_INDEX["CO"]
    cell.R_gas_d[j] = 5.0
    cell.N_d[j] = 1.0
    cell.N_b[j] = 0.5
    added = oneshot_bed_dense_species_sum_absorb(cell, species=("CO",), fraction=0.1)
    assert abs(added["CO"] - 0.35) < 1e-12
    assert abs(cell.N_d[j] - 1.35) < 1e-12
