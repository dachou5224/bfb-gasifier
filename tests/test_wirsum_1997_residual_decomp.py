"""R12 Startwert 之后 packed ‖F‖₂ 由 bed0 能量残差撑住，不是 R5。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.wirsum_1997_residual_decomp import (
    build_residual_decomp_payload,
    extract_r5_extent,
)


def test_extract_r5_extent_matches_exact_pair() -> None:
    gas = np.zeros(N_GAS, dtype=np.float64)
    xi = 2.0e4
    gas[GAS_SPECIES_INDEX["CO"]] = -xi
    gas[GAS_SPECIES_INDEX["O2"]] = -0.5 * xi
    gas[GAS_SPECIES_INDEX["CO2"]] = xi
    got = extract_r5_extent(gas)
    assert got["extent_mol_s"] == pytest.approx(xi)
    assert got["co_over_co2"] == pytest.approx(-1.0)
    assert got["o2_over_co2"] == pytest.approx(-0.5)


def test_after_r12_startwert_raw_l2_is_bed0_energy() -> None:
    payload = build_residual_decomp_payload()
    before = payload["before"]
    after = payload["after"]
    assert before["dominant_family"] == "gas"
    assert before["max_r12_extent_mol_s"] > 1.0e12
    assert abs(float(before["family_l2"]["energy"] or 0.0)) > 1.0e6

    assert after["dominant_family"] == "energy"
    assert after["linf_family"] == "energy"
    assert after["linf_label"] == "c0.energy"
    assert after["linf_over_l2"] > 0.99
    assert after["max_r12_extent_mol_s"] < 1.0e2
    assert float(after["family_l2"]["gas"] or 0.0) < 1.0e6
    assert abs(float(after["family_l2"]["energy"] or 0.0)) > 1.0e6
    bed0 = after["per_cell"][0]
    assert bed0["r5_looks_like_pair"] is False
    assert abs(float(bed0["energy_residual_W"] or 0.0)) > 1.0e6


def test_residual_decomp_does_not_project_or_use_experiment_co() -> None:
    import src.solvers.wirsum_1997_residual_decomp as mod

    source = inspect.getsource(mod)
    assert "same_phase_r5_startwert" in source
    assert "residual_clip" in source
    dumped = str(build_residual_decomp_payload())
    assert "experiment" not in dumped.lower()
    assert "y_CO_target" not in dumped
