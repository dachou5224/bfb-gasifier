"""阶段 F：关回流的床–自由空间界面。Wirsum exact 与 Hamel 扩展分名。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.physics.wirsum_1997_freeboard import (
    FREEBOARD_MODEL_HAMEL_EXTENDED,
    FREEBOARD_MODEL_WIRSUM_1997_EXACT,
)
from src.solvers.wirsum_1997_bed_freeboard import (
    OUTER_CRITERIA,
    BedFreeboardGeometry,
    build_stage_f_payload,
    evaluate_interface,
    gas_interface_gap,
    outer_criterion_coverage,
    route_total_gas,
    solid_edge_gaps,
)


def test_gas_interface_ignores_phase_split() -> None:
    total = np.array([1.0, 2.0, 0.5], dtype=np.float64)
    a = route_total_gas(0.2 * total, 0.8 * total)
    b = route_total_gas(0.9 * total, 0.1 * total)
    np.testing.assert_allclose(a, b)
    np.testing.assert_allclose(a, total)


def test_hamel_extended_ejects_more_than_wirsum_roof_only() -> None:
    geo = BedFreeboardGeometry()
    N = np.array([1.0, 0.0], dtype=np.float64)
    wir = evaluate_interface(
        geo=geo, N_d=N, N_b=N, m_char=0.32, m_ash=0.08,
        model=FREEBOARD_MODEL_WIRSUM_1997_EXACT, recycle=False,
    )
    ham = evaluate_interface(
        geo=geo, N_d=N, N_b=N, m_char=0.32, m_ash=0.08,
        model=FREEBOARD_MODEL_HAMEL_EXTENDED, recycle=False,
    )
    assert float(np.sum(ham.m_eject)) > float(np.sum(wir.m_eject))
    assert wir.model == FREEBOARD_MODEL_WIRSUM_1997_EXACT
    assert ham.model == FREEBOARD_MODEL_HAMEL_EXTENDED


def test_recycle_closed_has_zero_recycle_and_closed_edges() -> None:
    s = evaluate_interface(
        geo=BedFreeboardGeometry(),
        N_d=np.array([0.2, 0.8]),
        N_b=np.array([0.6, 0.2]),
        m_char=0.32,
        m_ash=0.08,
        model=FREEBOARD_MODEL_HAMEL_EXTENDED,
        recycle=False,
    )
    assert gas_interface_gap(s) == pytest.approx(0.0)
    gaps = solid_edge_gaps(s)
    assert gaps["eject_minus_return_leave_kg_s"] == pytest.approx(0.0)
    assert gaps["leave_minus_catch_fly_kg_s"] == pytest.approx(0.0)
    assert gaps["recycle_kg_s"] == pytest.approx(0.0)


def test_recycle_open_returns_cyclone_catch() -> None:
    kwargs = dict(
        geo=BedFreeboardGeometry(d_p=2.0e-4),
        N_d=np.array([0.2, 0.8]),
        N_b=np.array([0.6, 0.2]),
        m_char=0.32,
        m_ash=0.08,
        model=FREEBOARD_MODEL_HAMEL_EXTENDED,
        eta_cyclone=0.95,
    )
    off = evaluate_interface(**kwargs, recycle=False)
    on = evaluate_interface(**kwargs, recycle=True)
    assert off.fate == "reaches_next_cell"
    assert float(np.sum(off.m_leave)) > 0.0
    np.testing.assert_allclose(on.m_recycle, 0.95 * on.m_leave)
    np.testing.assert_allclose(off.m_recycle, 0.0)
    np.testing.assert_allclose(on.m_flyash, 0.05 * on.m_leave)


def test_outer_criteria_do_not_claim_temperature_or_energy() -> None:
    cov = outer_criterion_coverage(recycle=False)
    assert set(cov) == set(OUTER_CRITERIA)
    assert cov["temperature"] == "not_covered"
    assert cov["energy"] == "not_covered"
    assert cov["pyrolysis_release"] == "not_covered"
    assert cov["total_gas"] == "covered"
    assert cov["bed_surface_fallback_entrainment"] == "covered"
    assert cov["hydro_exchange"] == "frozen_context"
    assert cov["element_balance"] == "not_covered"


def test_payload_keeps_models_and_skips_experiment_co() -> None:
    payload = build_stage_f_payload()
    assert payload["solved"] is False
    assert payload["eject_hamel_minus_wirsum_kg_s"] > 0.0
    coarse = payload["recycle_closed_coarse"]
    assert coarse["wirsum_1997_exact"]["fate"] == "returns"
    assert coarse["hamel_extended"]["fate"] == "returns"
    fine_on = payload["recycle_open_fine_hamel"]
    assert fine_on["fate"] == "reaches_next_cell"
    assert float(np.sum(fine_on["m_recycle_kg_s"])) > 0.0
    dumped = str(payload)
    assert "y_CO_target" not in dumped
    assert "experiment" not in dumped.lower()


def test_module_does_not_import_production_freeboard() -> None:
    import src.solvers.wirsum_1997_bed_freeboard as mod

    src = inspect.getsource(mod)
    assert "from src.core.freeboard" not in src
    assert "import simulate" not in src
    assert "fastox" not in src
