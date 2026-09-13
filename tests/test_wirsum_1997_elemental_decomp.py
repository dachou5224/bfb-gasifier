"""P5 续：H/O 开口来自裸 R12 对消，不是热解 allocator。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.wirsum_1997_elemental_decomp import (
    build_p5_elemental_decomp_payload,
    extract_r12_extent,
    subtract_r12_pair,
)


def test_extract_r12_extent_matches_exact_pair() -> None:
    gas = np.zeros(N_GAS, dtype=np.float64)
    xi = 1.0e13
    gas[GAS_SPECIES_INDEX["H2"]] = -xi
    gas[GAS_SPECIES_INDEX["O2"]] = -0.5 * xi
    gas[GAS_SPECIES_INDEX["H2O"]] = xi
    got = extract_r12_extent(gas)
    assert got["extent_mol_s"] == pytest.approx(xi)
    assert got["h2_over_h2o"] == pytest.approx(-1.0)
    assert got["o2_over_h2o"] == pytest.approx(-0.5)
    leftover = subtract_r12_pair(gas, got["extent_mol_s"])
    assert float(np.max(np.abs(leftover))) == pytest.approx(0.0, abs=1e-6)


def test_huge_r12_pair_leaves_only_float_residue_in_h_o() -> None:
    from src.core.elemental_ledger import gas_element_molar_rates

    gas = np.zeros(N_GAS, dtype=np.float64)
    xi = 5.64e13
    gas[GAS_SPECIES_INDEX["H2"]] = -xi
    gas[GAS_SPECIES_INDEX["O2"]] = -0.5 * xi
    gas[GAS_SPECIES_INDEX["H2O"]] = xi
    atoms = gas_element_molar_rates(gas, fuel_type="coal")
    assert abs(atoms["H"]) / xi < 1e-12
    assert abs(atoms["O"]) / xi < 1e-12
    assert abs(atoms["C"]) == pytest.approx(0.0, abs=0.0)


def test_lu_x0_ho_gap_is_not_the_allocator_an_eq_b() -> None:
    payload = build_p5_elemental_decomp_payload()
    assert payload["max_r12_extent_mol_s"] > 1.0e12
    bed0, bed1, bed2 = payload["per_cell"]
    for row in payload["per_cell"]:
        assert row["r12_looks_like_pair"] is True
        assert row["pyrolysis_cache_used"] is False
        assert abs(float(row["total"]["C"] or 0.0)) < 1e-9

    assert bed0["has_solid_vm_or_moisture"] is True
    assert abs(float(bed0["gas_only"]["H"] or 0.0)) > 10.0
    assert abs(float(bed0["solid_only"]["H"] or 0.0)) > 10.0
    assert float(bed0["gas_solid_H_mismatch_over_solid_H"] or 1.0) < 1e-3

    assert bed1["has_solid_vm_or_moisture"] is False
    assert bed2["has_solid_vm_or_moisture"] is False
    assert abs(float(bed1["total"]["H"] or 0.0)) < 0.02
    assert abs(float(bed1["total"]["O"] or 0.0)) < 0.02


def test_decomp_does_not_touch_experiment_co() -> None:
    dumped = str(build_p5_elemental_decomp_payload())
    assert "experiment" not in dumped.lower()
    assert "y_CO_target" not in dumped
