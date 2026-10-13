"""阶段 F5：床–自由空间与旋风 cell 共用 ``N_fb``、``m_leave``。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.solvers.wirsum_1997_bed_fb_cyclone import (
    OUTER_CRITERIA,
    BedFreeboardGeometry,
    build_stage_f5_payload,
    couple,
    default_inlets,
    interface_to_cyclone,
    outer_criterion_coverage_coupled,
)
from src.solvers.wirsum_1997_bed_freeboard import route_total_gas


def test_coupler_does_not_cut_eta() -> None:
    import src.solvers.wirsum_1997_bed_fb_cyclone as mod

    src = inspect.getsource(mod)
    assert "cyclone_split(" not in src
    assert "evaluate_interface(" not in src
    assert "from src.physics.wirsum_1997_cyclone" not in src
    assert "from src.solvers.wirsum_1997_cyclone_cells import" in src


def test_handoff_uses_same_N_fb_and_leave() -> None:
    inn = default_inlets()
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    out = couple(geo=geo, inn=inn, recycle=False)
    expect_N = float(np.sum(route_total_gas(inn.N_d, inn.N_b)))
    face = interface_to_cyclone(geo=geo, inn=inn)
    assert out["face"]["N_fb"] == pytest.approx(expect_N)
    assert out["N_fb_handoff_gap"] == pytest.approx(0.0)
    assert out["m_leave_handoff_gap"] == pytest.approx(0.0)
    assert float(out["face"]["m_leave"][0]) == pytest.approx(float(face["m_leave"][0]))
    assert float(out["face"]["m_leave"][0]) > 0.0
    assert out["face"]["eject_minus_return_leave"] == pytest.approx(0.0)


def test_coarse_sends_zero_solids_to_cyclone() -> None:
    out = couple(geo=BedFreeboardGeometry(), inn=default_inlets(), recycle=False)
    assert out["face"]["fate"] == "returns"
    assert float(out["face"]["m_leave"][0]) == pytest.approx(0.0)
    assert out["cyc"]["m_underflow_kg_s"] == pytest.approx(0.0)
    assert out["cyc"]["physical_convergence"] is True


def test_recycle_flag_only_changes_bed_rez() -> None:
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    inn = default_inlets()
    off = couple(geo=geo, inn=inn, recycle=False)
    on = couple(geo=geo, inn=inn, recycle=True)
    assert off["cyc"]["m_rez_kg_s"] == pytest.approx(0.0)
    assert on["cyc"]["m_rez_kg_s"] == pytest.approx(on["cyc"]["m_underflow_kg_s"])
    assert off["face"]["N_fb"] == pytest.approx(on["face"]["N_fb"])
    assert off["face"]["m_leave"][0] == pytest.approx(on["face"]["m_leave"][0])
    assert on["cyc"]["dF_bed_dMret"] > 0.0
    assert off["cyc"]["dF_bed_dMret"] == pytest.approx(0.0)


def test_coverage_still_not_furnace() -> None:
    cov = outer_criterion_coverage_coupled()
    assert set(cov) == set(OUTER_CRITERIA)
    assert cov["total_gas"] == "covered"
    assert cov["energy"] == "not_covered"
    assert cov["pyrolysis_release"] == "not_covered"
    assert cov["element_balance"] == "not_covered"


def test_payload_skips_experiment() -> None:
    payload = build_stage_f5_payload()
    assert payload["solved"] is False
    assert payload["fine_recycle_closed"]["N_fb_handoff_gap"] == pytest.approx(0.0)
    assert payload["fine_recycle_closed"]["m_leave_handoff_gap"] == pytest.approx(0.0)
    assert payload["coarse_recycle_closed"]["m_leave_char_kg_s"] == pytest.approx(0.0)
    dumped = str(payload)
    assert "y_CO_target" not in dumped
    assert "engineering_lapple_1951" in payload["not"]
