"""工作包 D 的现场记录。不重解。"""

import json
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "docs/audit_reports/wp1_wp4_review_2026-09-28"


def test_p3_live_recert_passed_and_branches_stopped_on_the_first_explicit_step() -> None:
    p3 = json.loads((_DIR / "wpd_p3.json").read_text(encoding="utf-8"))
    fuel = json.loads((_DIR / "wpd_fuel.json").read_text(encoding="utf-8"))
    air = json.loads((_DIR / "wpd_air.json").read_text(encoding="utf-8"))
    assert p3["passed"] is True
    assert p3["gas_l2"] < 1.0e-6
    assert p3["transport_closed"] is True
    assert fuel["reached_target"] is False
    assert fuel["path"][0]["factor"] == 1.025
    assert fuel["path"][0]["reason"] == "merit_reject"
    assert fuel["stop_is_not_a_0p01_threshold"] is True
    assert air["reached_target"] is False
    assert [step["factor"] for step in air["path"]] == [0.975]
    assert air["path"][0]["reason"] == "merit_reject"


def test_r_series_keeps_the_near_root_direction_and_the_trajectory_end_states() -> None:
    rdir = _DIR.parent / "review_abcd_2026-09-28"
    t50 = json.loads((rdir / "r2_tplus50.json").read_text(encoding="utf-8"))
    traces = json.loads((rdir / "r1_wesseling_traces.json").read_text(encoding="utf-8"))
    assert t50["passed"] is False
    assert t50["stats"]["direction"] == "frozen_down_composite_up"
    assert t50["stats"]["checkpoint"] == "least_squares_x"
    assert t50["stats"]["j_budget"] == 4.0
    assert len(t50["x"]) == 50
    assert traces["continue_budget_s"] == 15.0
    assert traces["w_undecided"][1] == 1.0
    assert traces["reason_counts"]["both_window_exhausted"] > 0
    assert "dt_mismatch" not in traces["reason_counts"]
    assert traces["reason_counts"]["furnace_exit"] > 0


def test_wesseling_has_no_exit_prediction_while_transport_is_open() -> None:
    row = json.loads((_DIR / "wpd_wesseling.json").read_text(encoding="utf-8"))
    assert row["transport_closed"] is False
    assert row["prediction"] is None
    assert row["reason"] == "transport_not_closed"
    assert row["wesseling2_solved"] is False
    assert "steam_kg_h" in row["wesseling2_missing"]
    assert "primary_O2_Nm3_h" in row["wesseling2_missing"]
    assert row["trajectory_continue_s"] == 15.0
    assert row["blocks"]["passed"] is False
    assert row["input_path"]["target_reached"] is False
    assert row["input_path"]["intermediate_name"] == "not CASE_HTW_WESSELING_1"
    assert row["input_path"]["cold_start_of_target_run"] is False
