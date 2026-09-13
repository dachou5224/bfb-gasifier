"""同相 R12 Startwert：只改 x0 holdup，不进 residual。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.wirsum_1997_startwert_r12 import (
    apply_same_phase_r12_startwert,
    apply_same_phase_r12_startwert_to_cells,
    build_startwert_r12_ab_payload,
    same_phase_r12_overlap_mol_s,
)


def _vec(**kwargs: float) -> np.ndarray:
    n = np.zeros(N_GAS, dtype=np.float64)
    for name, value in kwargs.items():
        n[GAS_SPECIES_INDEX[name]] = value
    return n


def test_r12_startwert_consumes_limiting_h2() -> None:
    n = _vec(H2=1.0, O2=1.0, H2O=0.2)
    diag = apply_same_phase_r12_startwert(n)
    assert diag["xi12_mol_s"] == pytest.approx(1.0)
    assert n[GAS_SPECIES_INDEX["H2"]] == pytest.approx(0.0)
    assert n[GAS_SPECIES_INDEX["O2"]] == pytest.approx(0.5)
    assert n[GAS_SPECIES_INDEX["H2O"]] == pytest.approx(1.2)
    assert same_phase_r12_overlap_mol_s(n) == pytest.approx(0.0)


def test_r12_startwert_consumes_limiting_o2() -> None:
    n = _vec(H2=4.0, O2=1.0, H2O=0.0)
    apply_same_phase_r12_startwert(n)
    assert n[GAS_SPECIES_INDEX["H2"]] == pytest.approx(2.0)
    assert n[GAS_SPECIES_INDEX["O2"]] == pytest.approx(0.0)
    assert n[GAS_SPECIES_INDEX["H2O"]] == pytest.approx(2.0)
    assert same_phase_r12_overlap_mol_s(n) == pytest.approx(0.0)


def test_r12_startwert_conserves_h_and_o_atoms() -> None:
    n = _vec(H2=1.25, O2=5.94, H2O=3.0, CO=2.0, CH4=0.4)
    h_before = 2.0 * n[GAS_SPECIES_INDEX["H2"]] + 2.0 * n[GAS_SPECIES_INDEX["H2O"]]
    o_before = 2.0 * n[GAS_SPECIES_INDEX["O2"]] + n[GAS_SPECIES_INDEX["H2O"]]
    apply_same_phase_r12_startwert(n)
    h_after = 2.0 * n[GAS_SPECIES_INDEX["H2"]] + 2.0 * n[GAS_SPECIES_INDEX["H2O"]]
    o_after = 2.0 * n[GAS_SPECIES_INDEX["O2"]] + n[GAS_SPECIES_INDEX["H2O"]]
    assert h_after == pytest.approx(h_before)
    assert o_after == pytest.approx(o_before)
    assert n[GAS_SPECIES_INDEX["CO"]] == pytest.approx(2.0)
    assert n[GAS_SPECIES_INDEX["CH4"]] == pytest.approx(0.4)


def test_startwert_is_not_residual_clip_or_engineering_fastox() -> None:
    import src.solvers.wirsum_1997_startwert_r12 as mod

    source = inspect.getsource(mod)
    assert "apply_same_phase_r12_startwert" in source
    assert "project_bed_holdup_fast_oxidation_closure" in source
    assert "nr_line_search_fast_oxidation_projection_thesis" in source
    apply_src = inspect.getsource(apply_same_phase_r12_startwert)
    assert "global_residual" not in apply_src
    assert "clip" not in apply_src
    assert "rate_R5" not in apply_src
    cells_src = inspect.getsource(apply_same_phase_r12_startwert_to_cells)
    assert "global_residual" not in cells_src


def test_lu_x0_same_phase_r12_startwert_drops_overlap_and_extent() -> None:
    payload = build_startwert_r12_ab_payload()
    before = payload["before"]
    after = payload["after"]
    assert before["max_xi12_possible_mol_s"] > 0.1
    assert before["max_r12_extent_mol_s"] > 1.0e12
    assert after["max_xi12_possible_mol_s"] == pytest.approx(0.0, abs=1e-12)
    assert after["max_r12_extent_mol_s"] < 1.0e6
    assert after["residual_l2"] < 0.1 * before["residual_l2"]
    assert payload["name"] == "same_phase_r12_startwert"
    assert "residual_clip" in payload["not"]


def test_startwert_ab_does_not_touch_experiment_co() -> None:
    dumped = str(build_startwert_r12_ab_payload())
    assert "experiment" not in dumped.lower()
    assert "y_CO_target" not in dumped
