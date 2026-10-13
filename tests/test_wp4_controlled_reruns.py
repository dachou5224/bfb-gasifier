"""工作包 4 的记录：参数分支和 Wesseling 条件性预测。不重解。"""

import json
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_DIR = _ROOT / "docs/audit_reports/plan_ab_diagnosis_2026-09-28"


def test_fuel_plus_5pct_reaches_the_target_and_air_stops_at_the_half_step() -> None:
    fuel = json.loads((_DIR / "wp4_parameter_branches_fuel_rest.json").read_text(encoding="utf-8"))
    air = json.loads((_DIR / "wp4_parameter_branches_air_rest.json").read_text(encoding="utf-8"))
    fuel_row = fuel["parameter_branch"][0]
    air_row = air["parameter_branch"][0]
    assert fuel_row["class"] == "parameter_branch"
    assert fuel_row["reached_target"] is True
    assert fuel_row["path"][-1]["passed"] is True
    assert fuel_row["path"][-1]["reason"] == "block_gates"
    assert fuel_row["path"][-1]["gate_ratio"]["gas"] < 1.0
    assert air_row["reached_target"] is False
    assert air_row["certified_factor"] == pytest.approx(0.975)
    assert air_row["path"][-1]["T_bed"] == pytest.approx(air_row["path"][0]["T_bed"])
    assert air_row["path"][-1]["gate_ratio"]["gas"] > 1.0


def test_wesseling1_unconverged_state_is_not_reported_as_a_prediction() -> None:
    report = json.loads((_DIR / "wp4_wesseling1.json").read_text(encoding="utf-8"))
    assert report["conditional"] is True
    assert report["passed"] is False
    assert report["prediction"] is None
    assert report["initial_w_undecided"] == [1.0, 1.0]
    assert report["assumptions"]["secondary_air_placed_at_xi"] is False
    assert report["assumptions"]["steam_status"] == "assumed_zero_because_missing"
