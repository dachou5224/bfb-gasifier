"""``M = m_in / K`` Startwert：只填 K>0 的 char/ash，不把 VM 写入 Eq.2.7。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.solvers.wirsum_1997_holdup_startwert import (
    apply_km_holdup_startwert_to_cells,
    build_holdup_startwert_ab_payload,
    km_support_holdup,
)


def test_km_support_inverts_positive_k_and_leaves_k0_zero() -> None:
    m_in = np.array([[0.32, 0.37], [0.08, 0.0]], dtype=np.float64)
    K = np.array([[0.08, 0.0], [0.08, 0.0]], dtype=np.float64)
    M = km_support_holdup(m_in, K)
    assert M[0, 0] == pytest.approx(0.32 / 0.08)
    assert M[1, 0] == pytest.approx(0.08 / 0.08)
    assert M[0, 1] == pytest.approx(0.0)
    assert M[1, 1] == pytest.approx(0.0)


def test_holdup_startwert_is_not_residual_anchor() -> None:
    import src.solvers.wirsum_1997_holdup_startwert as mod

    source = inspect.getsource(mod)
    assert "holdup_k0_anchor" in source
    assert "source_minus_M" in source
    apply_src = inspect.getsource(apply_km_holdup_startwert_to_cells)
    assert "global_residual" not in apply_src


def test_lu_km_holdup_startwert_closes_char_ash_energy_outlet() -> None:
    payload = build_holdup_startwert_ab_payload()
    before = payload["before"]["per_cell"][0]
    after = payload["after"]["per_cell"][0]
    assert before["holdup_zero"] is True
    assert before["H_solid_out_zero"] is True
    assert abs(float(before["residual_W"] or 0.0)) > 1.0e7

    assert after["holdup_zero"] is False
    assert after["H_solid_out_zero"] is False
    assert float(after["m_holdup_kg"]["char"] or 0.0) > 1.0
    assert float(after["m_holdup_kg"]["ash"] or 0.0) > 0.1
    assert float(after["m_holdup_kg"]["vm"] or 0.0) == pytest.approx(0.0, abs=1e-12)
    assert after["K_auf_max_1_per_s"]["vm"] == pytest.approx(0.0, abs=0.0)
    assert abs(float(after["residual_W"] or 0.0)) < abs(float(before["residual_W"] or 0.0))
    char_gap = float(after["H_zu_by_comp_W"]["char"] or 0.0) - float(
        after["H_solid_out_by_comp_W"]["char"] or 0.0
    )
    assert abs(char_gap) < 2.0e6


def test_holdup_ab_does_not_touch_experiment_co() -> None:
    dumped = str(build_holdup_startwert_ab_payload())
    assert "experiment" not in dumped.lower()
    assert "y_CO_target" not in dumped
