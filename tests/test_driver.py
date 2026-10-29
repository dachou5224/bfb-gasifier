"""求解驱动。完整工况标 slow，不进默认集。"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from hamel.cases.balance import BalanceGap
from hamel.cases.config import SolverConfig
from hamel.kinetics.common import KineticsInputError
from hamel.solver.campaign import _process, execute, overlay, run_campaign, run_specs, write_report
from hamel.solver.driver import SolveResult, solve_aligned
from hamel.solver.gates import OuterMonitor, ResidualBlocks, arrays_allowed, balance_met, inner_met, outer_met
from hamel.solver.newton import damped_newton, normalized_merit

_CASE = Path(__file__).resolve().parents[1] / "data" / "cases" / "vtt_pressurized.json"


def _linear(config=None):
    config = config or SolverConfig()
    target = np.array([1.0, 2.0, 3.0])
    blocks = ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 10.0, 10.0, 10.0)

    def evaluate(trial):
        return np.asarray(trial, dtype=float) - target

    result = damped_newton(
        np.array([0.2, 0.2, 0.2]), evaluate, np.ones(3),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        lambda trial: bool(np.all(trial > 0.0)),
        lambda residual: normalized_merit(residual, blocks), config,
    )
    return result


def test_damped_newton_converges_on_a_linear_system():
    result = _linear()
    assert result.status == "CONVERGED"
    assert result.iterations == 1
    assert result.x == pytest.approx([1.0, 2.0, 3.0])


def test_damping_fails_when_no_trial_is_feasible():
    """差分只动一个分量，牛顿试探同时动两个分量，后者全部拒绝。"""
    config = SolverConfig()
    start = np.array([1.0, 1.0])
    blocks = ResidualBlocks(np.array([0, 1]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)

    def feasible(trial):
        delta = np.asarray(trial, dtype=float) - start
        changed = np.flatnonzero(np.abs(delta) > 0.0)
        return bool(changed.size == 0 or (changed.size == 1 and abs(float(delta[changed[0]])) < 1.0e-6))

    result = damped_newton(
        start, lambda trial: trial - np.array([0.2, 0.2]), np.ones(2),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        feasible, lambda residual: normalized_merit(residual, blocks), config,
    )
    assert result.status == "DAMPING_FAILED"
    assert result.iterations == 1


def test_rate_exception_shrinks_the_step_and_still_converges():
    """求解手段：ADR 0022。越界异常只否决这一步。"""
    config = SolverConfig()
    blocks = ResidualBlocks(np.array([0]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)
    raised = {"count": 0}

    def evaluate(trial):
        value = float(trial[0])
        if value > 2.0:
            raised["count"] += 1
            raise KineticsInputError("速率式超出原文适用范围")
        return np.array([value * value - 1.0])

    result = damped_newton(
        np.array([0.1]), evaluate, np.ones(1),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        lambda trial: bool(np.all(np.isfinite(trial))),
        lambda residual: normalized_merit(residual, blocks), config,
    )
    assert raised["count"] >= 1
    assert result.status == "CONVERGED"
    assert result.x[0] == pytest.approx(1.0)


def test_outer_criteria_use_plant_scales_and_accept_paired_zeros():
    config = SolverConfig()
    previous = OuterMonitor(np.array([1000.0]), np.array([0.0]), np.array([2.0]), np.array([0.0]))
    current = OuterMonitor(np.array([1000.05]), np.array([0.0]), np.array([2.0]), np.array([0.0]))
    assert outer_met(
        current, previous, m_dot_ref=1.0, n_dot_ref=1.0, m_tot=4.0,
        temperature_tolerance=config.outer_temperature_tolerance_k, relative_tolerance=config.outer_relative_tolerance,
    )
    moved = OuterMonitor(np.array([1000.05]), np.array([0.0]), np.array([2.0]), np.array([1.0]))
    assert not outer_met(
        moved, previous, m_dot_ref=1.0, n_dot_ref=1.0, m_tot=4.0,
        temperature_tolerance=config.outer_temperature_tolerance_k, relative_tolerance=config.outer_relative_tolerance,
    )
    assert outer_met(
        moved, previous, m_dot_ref=1.0, n_dot_ref=1.0e5, m_tot=4.0,
        temperature_tolerance=config.outer_temperature_tolerance_k, relative_tolerance=config.outer_relative_tolerance,
    )


def test_balance_and_nonnegative_checks():
    closed = balance_met(BalanceGap(0.0, 0.0, 0.0, 5.0), (2.0, 3.0, 4.0), 5.0, 1.0e4, 1.0e-3, 5.0e-3)
    assert closed["ok"]
    carbon = balance_met(BalanceGap(0.01, 0.0, 0.0, 5.0), (1.0, 3.0, 4.0), 5.0, 1.0e4, 1.0e-3, 5.0e-3)
    assert not carbon["ok"]
    assert arrays_allowed(([1.0, 0.0],), ([0.2],), [1100.0], [0.0])
    assert not arrays_allowed(([-1e-9],), ([0.2],), [1100.0], [0.0])
    assert not arrays_allowed(([1.0],), ([0.2],), [0.0], [0.0])
    assert not arrays_allowed(([1.0],), ([0.2],), [1100.0], [-0.1])


def test_outer_loop_recomputes_and_passes_both_gates():
    config = SolverConfig()
    target = np.array([1.0, 0.5, 1100.0])
    blocks = ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 1.0, 1.0, 1.0e3)
    seen = {"consistent": 0}

    def refresh(x, iteration, consistency):
        if x is not None:
            target[:] = np.asarray(x, dtype=float)
        vector = np.array([1.2, 0.6, 1110.0]) if x is None else np.asarray(x, dtype=float)
        return type("Opened", (), {"x": vector, "info": {"round": iteration, "consistency": consistency}})()

    def consistent(values, iteration):
        seen["consistent"] += 1
        residual = np.asarray(values, dtype=float) - target
        return type("Checked", (), {"residual": residual, "blocks": blocks, "info": {"round": iteration}})()

    result = solve_aligned(
        refresh, lambda trial, opened: np.asarray(trial, dtype=float) - target,
        lambda opened: (blocks, np.ones(3), lambda trial: bool(np.all(trial > 0.0)), (1.0, 1.0, 1.0e3, 1.0)),
        lambda opened, values: OuterMonitor(values[2:3], values[1:2], values[1:2], values[0:1]),
        consistent, config,
    )
    assert result.status == "CONVERGED"
    assert result.rounds == 2
    assert seen["consistent"] == 1
    assert result.x == pytest.approx([1.0, 0.5, 1100.0])


def test_run_table_changes_one_item_and_records_unwired_rows():
    specs = run_specs()
    # ADR 0030 的敏感性运行由 16 次增为 17 次（ADR 0048 第 3 项加 S17）。
    assert [item.name for item in specs] == ["baseline", *[f"S{index}" for index in range(1, 18)]]
    assert len({item.change for item in specs}) == 18
    assert next(item for item in specs if item.name == "S1").wired
    assert [item.name for item in specs if not item.wired] == []
    assert all(item.wired for item in specs)


def test_overlay_changes_only_the_named_value(tmp_path):
    original = json.loads(_CASE.read_text(encoding="utf-8"))
    specs = {item.name: item for item in run_specs()}
    dest = tmp_path / "s2.json"
    assert overlay(_CASE, dest, specs["S2"]) is False
    patched = json.loads(dest.read_text(encoding="utf-8"))
    assert _field(patched, "distributor_orifice_count") == 95.6
    assert _field(patched, "fuel_size_class_diameters") == _field(original, "fuel_size_class_diameters")
    fine = tmp_path / "s7.json"
    overlay(_CASE, fine, specs["S7"])
    diameters = _field(json.loads(fine.read_text(encoding="utf-8")), "fuel_size_class_diameters")
    assert diameters[0] == pytest.approx(3.28e-5)
    assert diameters[1:] == _field(original, "fuel_size_class_diameters")[1:]
    assert overlay(_CASE, tmp_path / "s1.json", specs["S1"]) is True


def test_toy_report_covers_the_required_sections(tmp_path):
    solved = _toy(SolverConfig())
    record = {
        "name": "baseline", "change": "不改冻结项", "source": "ADR 0030", "status": solved.status,
        "reason": solved.reason, "rounds": solved.rounds, "iterations": solved.iterations,
        "ratios": {"gas": 0.0, "solid": 0.0, "energy": 0.0},
        "balance": {"ok": True, "ratios": {"C": 0.0, "H": 0.0, "O": 0.0, "energy": 0.0}, "k_w": 0.0, "heat_w": 0.0},
        "nonnegative": True,
        "measurements": {"T_FRB_K": 1100.0, "X_C": 50.0, "CO": 0.1, "CO2": 0.1, "H2": 0.1, "CH4": 0.01, "H2O": 0.1},
        "diagnostics": {"rounds": solved.diagnostics["rounds"]},
    }
    payload = write_report(tmp_path, [record, {
        "name": "S1", "change": "k_W = 0", "source": "ADR 0029", "status": "NOT_APPLIED",
        "reason": "玩具不跑装置", "measurements": None, "balance": None, "nonnegative": None,
        "ratios": {}, "diagnostics": {}, "rounds": 0, "iterations": 0,
    }], title="玩具工况", toy=True)
    text = (tmp_path / "driver_report.md").read_text(encoding="utf-8")
    assert "玩具工况" in text and "数值加密" in text and "F24" in text and "F25" in text and "F26" in text
    assert "各验收量共享壁面散热这一校准输入" in text
    assert "全自由空间闭合" in text and "数值积分与闭式" in text and "超过容差一半" in text
    assert "T_FRB_K" in text and "TIMEOUT" not in payload["runs"][0]["status"]
    assert payload["runs"][0]["acceptance"]["premises_met"] is True
    assert payload["differences"][0]["status"] == "NOT_APPLIED"


def _diagnostics():
    return {"rounds": [{"info": {
        "a_free_per_s": np.ones((2, 2)),
        "drain_char_kg_s": np.array([0.25]),
        "hover": np.array([[True, False]]),
        "k_w": np.float64(3.5),
        "bed_intervals": np.int64(64),
    }}]}


def test_execute_report_keeps_numpy_array_shapes(tmp_path, monkeypatch):
    """单元素数组保持列表。失败和收敛都要能写进 JSON。"""
    queued = [
        SolveResult("NOT_CONVERGED", "内层未收敛", 1, 4, diagnostics=_diagnostics()),
        SolveResult(
            "CONVERGED", "", 2, 5, {"gas": np.float64(0.0)},
            {"ok": True, "ratios": {"C": np.float64(0.0)}}, True,
            {"T_FRB_K": np.float64(1100.0), "X_C": np.float64(50.0), "CO": np.float64(0.1), "CO2": np.float64(0.1), "H2": np.float64(0.1), "CH4": np.float64(0.01), "H2O": np.float64(0.1)},
            _diagnostics(),
        ),
    ]

    def fake_solve(path, case_id="1", *, adiabatic=False, deadline=None, config=None):
        return queued.pop(0)

    monkeypatch.setattr("hamel.solver.campaign.solve_case", fake_solve)
    failed = execute("baseline", "case.json", False)
    solved = execute("S1", "case.json", True)
    for record in (failed, solved):
        info = record["diagnostics"]["rounds"][0]["info"]
        assert info["a_free_per_s"] == [[1.0, 1.0], [1.0, 1.0]]
        assert info["drain_char_kg_s"] == [0.25]
        assert type(info["drain_char_kg_s"][0]) is float
        assert info["hover"] == [[True, False]]
        assert type(info["k_w"]) is float and info["k_w"] == 3.5
        assert type(info["bed_intervals"]) is int and info["bed_intervals"] == 64
    assert failed["status"] == "NOT_CONVERGED"
    assert solved["status"] == "CONVERGED"
    assert type(solved["ratios"]["gas"]) is float
    write_report(tmp_path, [failed, solved], title="诊断序列化", toy=True)
    saved = json.loads((tmp_path / "driver_report.json").read_text(encoding="utf-8"))
    assert [item["status"] for item in saved["runs"]] == ["NOT_CONVERGED", "CONVERGED"]
    assert saved["runs"][1]["diagnostics"]["rounds"][0]["info"]["a_free_per_s"] == [[1.0, 1.0], [1.0, 1.0]]
    assert saved["runs"][0]["diagnostics"]["rounds"][0]["info"]["drain_char_kg_s"] == [0.25]


def _toy(config):
    target = np.array([1.0, 0.5, 1100.0])
    blocks = ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 1.0, 1.0, 1.0e3)

    def refresh(x, iteration, consistency):
        if x is not None:
            target[:] = np.asarray(x, dtype=float)
        vector = np.array([1.2, 0.6, 1110.0]) if x is None else np.asarray(x, dtype=float)
        return type("Opened", (), {"x": vector, "info": {"round": iteration}})()

    return solve_aligned(
        refresh, lambda trial, opened: np.asarray(trial, dtype=float) - target,
        lambda opened: (blocks, np.ones(3), lambda trial: bool(np.all(trial > 0.0)), (1.0, 1.0, 1.0e3, 1.0)),
        lambda opened, values: OuterMonitor(np.asarray(values)[2:3], np.asarray(values)[1:2], np.asarray(values)[1:2], np.asarray(values)[0:1]),
        lambda values, iteration: type("Checked", (), {"residual": np.asarray(values) - target, "blocks": blocks, "info": {}})(),
        config,
    )


def test_newton_journal_does_not_change_the_accepted_step(tmp_path):
    """记录只追加。有无 journal 时 x 与残差逐位相同，失败点写出 npz。"""
    plain = _linear()
    config = SolverConfig()
    target = np.array([1.0, 2.0, 3.0])
    blocks = ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 10.0, 10.0, 10.0)
    journal = {"iterations": [], "blocks": blocks, "names": ["gas", "solid", "energy"], "bottom_norm": 0.0}
    traced = damped_newton(
        np.array([0.2, 0.2, 0.2]), lambda trial: np.asarray(trial, dtype=float) - target, np.ones(3),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        lambda trial: bool(np.all(trial > 0.0)),
        lambda residual: normalized_merit(residual, blocks), config, journal=journal,
    )
    assert traced.status == plain.status == "CONVERGED"
    assert np.array_equal(traced.x, plain.x)
    assert np.array_equal(traced.residual, plain.residual)
    assert journal["iterations"][0]["trials"][0]["accepted"] is True
    assert journal["iterations"][0]["blocks"]["bottom"] == 0.0

    start = np.array([1.0, 1.0])
    fail_blocks = ResidualBlocks(np.array([0, 1]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)

    def feasible(trial):
        delta = np.asarray(trial, dtype=float) - start
        changed = np.flatnonzero(np.abs(delta) > 0.0)
        return bool(changed.size == 0 or (changed.size == 1 and abs(float(delta[changed[0]])) < 1.0e-6))

    def evaluate(trial):
        return np.asarray(trial, dtype=float) - np.array([0.2, 0.2])

    silent = damped_newton(
        start, evaluate, np.ones(2), lambda residual: inner_met(residual, fail_blocks, config.residual_tolerance),
        feasible, lambda residual: normalized_merit(residual, fail_blocks), config,
    )
    failure = {"iterations": [], "blocks": fail_blocks, "names": ["n0", "n1"], "failure_path": str(tmp_path / "failure.npz"),
               "explain": lambda trial: "可行条件：不能同时改两个分量"}
    logged = damped_newton(
        start, evaluate, np.ones(2), lambda residual: inner_met(residual, fail_blocks, config.residual_tolerance),
        feasible, lambda residual: normalized_merit(residual, fail_blocks), config, journal=failure,
    )
    assert silent.status == logged.status == "DAMPING_FAILED"
    assert np.array_equal(silent.x, logged.x)
    assert np.array_equal(silent.residual, logged.residual)
    trials = failure["iterations"][-1]["trials"]
    assert len(trials) == config.damping_index_limit + 1
    assert all(item["status"] == "reject" and item["reason"] == "可行条件：不能同时改两个分量" for item in trials)
    saved = np.load(tmp_path / "failure.npz")
    assert np.array_equal(saved["x"], logged.x)
    assert np.array_equal(saved["F"], logged.residual)
    assert saved["J"].shape == (2, 2) and saved["dx"].shape == (2,)


def test_shrink_exception_is_named_in_the_trial_record():
    config = SolverConfig()
    blocks = ResidualBlocks(np.array([0]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)
    journal = {"iterations": [], "blocks": blocks, "names": ["x"]}

    def evaluate(trial):
        value = float(trial[0])
        if value > 2.0:
            raise KineticsInputError("速率式超出原文适用范围")
        return np.array([value * value - 1.0])

    result = damped_newton(
        np.array([0.1]), evaluate, np.ones(1),
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        lambda trial: bool(np.all(np.isfinite(trial))),
        lambda residual: normalized_merit(residual, blocks), config, journal=journal,
    )
    assert result.status == "CONVERGED"
    reasons = [item["reason"] for row in journal["iterations"] for item in row["trials"]]
    assert any(reason.startswith("KineticsInputError:") for reason in reasons)


def test_wall_clock_kill_is_recorded_as_timeout():
    assert _process([sys.executable, "-c", "import time; time.sleep(30)"], 0.2) is None


@pytest.mark.slow
def test_case1_campaign_stops_at_the_wall_clock(tmp_path):
    """第 1 关全套运行。默认集不跑。每次运行到 wall_timeout_s 记 TIMEOUT。"""
    records = run_campaign(_CASE, tmp_path)
    assert [item["name"] for item in records] == [item.name for item in run_specs()]
    assert all(item["status"] in {"CONVERGED", "NOT_CONVERGED", "DAMPING_FAILED", "JACOBIAN_MISMATCH", "TIMEOUT", "DOMAIN_STOP", "NOT_APPLIED", "ERROR"} for item in records)


def _field(doc, name):
    found = []

    def walk(node):
        if isinstance(node, dict):
            child = node.get(name)
            if isinstance(child, dict) and "value" in child:
                found.append(child["value"])
            for item in node.values():
                walk(item)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(doc)
    assert len(found) == 1
    return found[0]
