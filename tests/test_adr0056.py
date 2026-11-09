"""ADR 0056：起步飞灰份额与（第 5 项的）O₂ 幂变量投影。"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.solver import driver
from hamel.solver.driver import PHI_INIT_CAP, solve_aligned
from hamel.solver.gates import OuterMonitor, ResidualBlocks
from hamel.solver.initial import InitialStateError
from hamel.solver.newton import NewtonResult, IDENTITY
from hamel.solver.precalc import PrecalcError, PrecalcResult, SolidContext, _drain_terms
from test_assembly import make_frozen, make_state

_CASE = Path("data/cases/vtt_pressurized.json")


def test_phi_init_cap_is_defined_once_and_clips_at_s18():
    """0,3 只在 PHI_INIT_CAP 一处；φ ≤ 0,3 时 φ_init = φ。求解手段：ADR 0056。"""
    assert PHI_INIT_CAP == 0.3
    assert min(1.0, PHI_INIT_CAP) == 0.3
    assert min(0.3, PHI_INIT_CAP) == 0.3
    assert min(0.2, PHI_INIT_CAP) == 0.2


def test_precalculate_override_writes_phi_before_any_physics(monkeypatch):
    import hamel.solver.precalc as precalc

    case = SimpleNamespace(fly_ash_fraction=1.0)
    monkeypatch.setattr(precalc, "_load", lambda *a, **k: case)
    monkeypatch.setattr(precalc, "_equilibrium", lambda *a, **k: (_ for _ in ()).throw(PrecalcError("after phi")))
    with pytest.raises(PrecalcError, match="after phi"):
        precalc._precalculate("x", "1", None, None, 0.0, 0.0, False, 0, False, None, 0.0, 0.0, None, None, 0.3)
    assert case.fly_ash_fraction == 0.3


def test_drain_and_fly_ash_target_follow_phi():
    """排灰与飞灰目标随 φ 切换。ADR 0055 第 4 项；由本卡的 φ 覆盖传入。"""
    released = SimpleNamespace(
        fuel_up_kg_s=np.zeros((1, 1)), undried_exit_kg_s=np.zeros(1), dried_exit_kg_s=np.zeros(1),
    )
    state = SimpleNamespace(bed_char_kg=np.array([[0.01]]))
    physical = SimpleNamespace(moisture_mass_fraction=0.0, ash_mass_fraction_wf=1.0e-3)
    inert, masses = np.array([0.38]), np.array([0.4])

    def terms(phi):
        case = SimpleNamespace(physical=physical, fuel_kg_s=1.0, n_free=0, fly_ash_fraction=phi)
        return _drain_terms(case, state, released, inert, masses)

    start, switched = terms(0.3), terms(1.0)
    assert start["target"] == pytest.approx(0.3 * start["available"])
    assert start["ash"] == pytest.approx(0.7 * start["available"])
    assert switched["target"] == pytest.approx(switched["available"])
    assert switched["ash"] == 0.0 and np.all(switched["char"] == 0.0)


def _linear_aligned(config, **kwargs):
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

    class Rule:
        transform = IDENTITY
        final = None

        def trial(self, trial):
            return bool(np.all(trial > 0.0))

    result = solve_aligned(
        refresh, lambda trial, opened: np.asarray(trial, dtype=float) - target,
        lambda opened: (blocks, np.ones(3), Rule(), (1.0, 1.0, 1.0e3, 1.0)),
        lambda opened, values: OuterMonitor(values[2:3], values[1:2], values[1:2], values[0:1]),
        consistent, config, **kwargs,
    )
    return result, seen


def test_outer_gate_waits_until_round_three_when_phi_switches():
    """φ_init ≠ φ 时第 2 轮不判外层门，最早第 3 轮。"""
    result, seen = _linear_aligned(SolverConfig(outer_round_limit=2), gate_from=3)
    assert result.status == "NOT_CONVERGED"
    assert result.reason == "外层超过 outer_round_limit"
    assert result.rounds == 2 and seen["consistent"] == 0
    result, seen = _linear_aligned(SolverConfig(), gate_from=3)
    assert result.status == "CONVERGED" and result.rounds == 3 and seen["consistent"] == 1


def test_outer_gate_timing_unchanged_when_phi_init_equals_phi():
    """φ_init = φ 时判定时机与现行相同（第 2 轮）。"""
    result, seen = _linear_aligned(SolverConfig())
    assert result.status == "CONVERGED" and result.rounds == 2 and seen["consistent"] == 1


def test_switch_precalc_failure_is_labelled():
    def refresh(x, iteration, consistency):
        if iteration == 2:
            raise PrecalcError("预计算失败")
        return type("Opened", (), {"x": np.array([1.0]), "info": {}})()

    config = SolverConfig(outer_round_limit=3)
    blocks = ResidualBlocks(np.array([0]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)

    class Rule:
        transform = IDENTITY
        final = None

        def trial(self, trial):
            return True

    result = solve_aligned(
        refresh, lambda trial, opened: np.asarray(trial, dtype=float) - 1.0,
        lambda opened: (blocks, np.ones(1), Rule(), (1.0, 1.0, 1.0, 1.0)),
        lambda opened, values: OuterMonitor(np.ones(1), np.ones(1), np.ones(1), np.ones(1)),
        lambda values, iteration: None, config, switch_at=2,
    )
    assert result.status == "DOMAIN_STOP" and result.rounds == 2
    assert "切换时" in result.reason


def test_round_one_inner_failure_is_labelled_before_the_switch(monkeypatch):
    def fake_newton(x0, *args, **kwargs):
        x = np.asarray(x0, dtype=float)
        return NewtonResult("DAMPING_FAILED", "没有使归一残差严格下降的可行步", x, x, 1, 1.0)

    monkeypatch.setattr(driver, "damped_newton", fake_newton)
    config = SolverConfig(outer_round_limit=3)
    blocks = ResidualBlocks(np.array([0]), np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)

    class Rule:
        transform = IDENTITY
        final = None

        def trial(self, trial):
            return True

    result = solve_aligned(
        lambda x, i, c: type("Opened", (), {"x": np.array([1.0]), "info": {}})(),
        lambda trial, opened: np.asarray(trial, dtype=float) - 1.0,
        lambda opened: (blocks, np.ones(1), Rule(), (1.0, 1.0, 1.0, 1.0)),
        lambda opened, values: OuterMonitor(np.ones(1), np.ones(1), np.ones(1), np.ones(1)),
        lambda values, iteration: None, config, switch_at=2,
    )
    assert result.status == "DAMPING_FAILED" and result.rounds == 1
    assert "切换前" in result.reason


def _stub(phi):
    available = 1.0e-4
    frozen = make_frozen(
        fly_ash_fraction=phi, ash_target_kg_s=phi * available,
        drain_ash_kg_s=(1.0 - phi) * available, bed_ash_fine_share=0.1,
    )
    state = make_state(frozen)
    return PrecalcResult(
        frozen, state, 971.0, {}, "B", 1e-4, 1e-4, 0.0, 0.0, (), 0.0, 1073.0, "B",
        (1, 1e-4, "B", 0.0), (), (), 8, (), {}, SolidContext(SimpleNamespace(), SimpleNamespace()),
    )


def _patch_inner(monkeypatch, newton=None):
    def default_newton(x0, *args, **kwargs):
        x = np.asarray(x0, dtype=float)
        journal = kwargs.get("journal")
        if journal is not None:
            journal.setdefault("iterations", []).append({"iteration": 1, "trials": [], "merit": 0.0})
        return NewtonResult("CONVERGED", "", x, np.zeros_like(x), 1, 0.0)

    monkeypatch.setattr(driver, "damped_newton", newton or default_newton)


def test_phi_one_uses_init_then_switches(monkeypatch):
    """φ = 1：初值构造与第 1 轮冻结量 φ = 0,3，第 2 轮起为 1；排灰与飞灰目标随 φ 切换。"""
    seen = []

    def fake_precalculate(*args, **kwargs):
        phi = float(kwargs["fly_ash_fraction"])
        seen.append(phi)
        return _stub(phi)

    def fake_construct(first, precalc_fn, m_ref, **kwargs):
        nxt = precalc_fn(first.state, (first.frozen.drain_ash_kg_s, first.frozen.drain_char_kg_s))
        records = kwargs.get("records")
        if records is not None:
            records.append({"loop": 1, "phi": nxt.frozen.fly_ash_fraction})
        return nxt, records if records is not None else [], ""

    _patch_inner(monkeypatch)
    monkeypatch.setattr(driver, "precalculate", fake_precalculate)
    monkeypatch.setattr(driver, "construct", fake_construct)
    result = driver.solve_case(str(_CASE), "1", config=SolverConfig(outer_round_limit=2, nu_max=1))
    assert seen[:2] == [0.3, 0.3] and seen[2] == 1.0
    assert result.diagnostics["initial"]["phi_init"] == 0.3
    assert result.diagnostics["initial"]["phi"] == 1.0
    round1, round2 = result.diagnostics["rounds"][0], result.diagnostics["rounds"][1]
    assert round1["phi_init"] == 0.3 and round1["phi"] == 1.0
    assert round1["info"]["fly_ash_fraction"] == 0.3
    assert round1["info"]["ash_target_kg_s"] == pytest.approx(0.3e-4)
    assert round1["info"]["drain_ash_kg_s"] == pytest.approx(0.7e-4)
    assert round2["switched"] is True
    assert round2["info"]["fly_ash_fraction"] == 1.0
    assert round2["info"]["ash_target_kg_s"] == pytest.approx(1.0e-4)
    assert round2["info"]["drain_ash_kg_s"] == 0.0
    assert result.status == "NOT_CONVERGED" and result.rounds == 2


def test_phi_below_cap_does_not_switch(monkeypatch, tmp_path):
    """φ = 0,2 时 φ_init = 0,2，切换逻辑不起作用。"""
    doc = json.loads(_CASE.read_text(encoding="utf-8"))
    doc["cases"]["1"]["fly_ash_fraction"]["value"] = 0.2
    path = tmp_path / "phi02.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    seen = []

    def fake_precalculate(*args, **kwargs):
        phi = float(kwargs["fly_ash_fraction"])
        seen.append(phi)
        return _stub(phi)

    def fake_construct(first, precalc_fn, m_ref, **kwargs):
        return first, kwargs.get("records") or [], ""

    _patch_inner(monkeypatch)
    monkeypatch.setattr(driver, "precalculate", fake_precalculate)
    monkeypatch.setattr(driver, "construct", fake_construct)
    result = driver.solve_case(str(path), "1", config=SolverConfig(outer_round_limit=2, nu_max=1))
    assert seen and all(phi == 0.2 for phi in seen)
    assert result.diagnostics["initial"]["phi_init"] == 0.2
    assert result.diagnostics["initial"]["phi"] == 0.2
    assert "switched" not in (result.diagnostics.get("rounds") or [{}])[-1]


def test_infeasible_at_phi_init_stops_without_lowering(monkeypatch):
    seen = []

    def fake_precalculate(*args, **kwargs):
        seen.append(float(kwargs["fly_ash_fraction"]))
        return _stub(float(kwargs["fly_ash_fraction"]))

    def fake_construct(*args, **kwargs):
        raise InitialStateError("φ_init 下固相不可行")

    _patch_inner(monkeypatch)
    monkeypatch.setattr(driver, "precalculate", fake_precalculate)
    monkeypatch.setattr(driver, "construct", fake_construct)
    result = driver.solve_case(str(_CASE), "1", config=SolverConfig(outer_round_limit=2, nu_max=1))
    assert result.status == "DOMAIN_STOP"
    assert "切换前" in result.reason
    assert seen == [0.3]
    assert result.diagnostics["initial"]["phi_init"] == 0.3


def test_switch_precalc_failure_on_the_driver_path(monkeypatch):
    def fake_precalculate(*args, **kwargs):
        phi = float(kwargs["fly_ash_fraction"])
        if phi == 1.0:
            raise PrecalcError("切换预计算失败")
        return _stub(phi)

    def fake_construct(first, precalc_fn, m_ref, **kwargs):
        return first, kwargs.get("records") or [], ""

    _patch_inner(monkeypatch)
    monkeypatch.setattr(driver, "precalculate", fake_precalculate)
    monkeypatch.setattr(driver, "construct", fake_construct)
    result = driver.solve_case(str(_CASE), "1", config=SolverConfig(outer_round_limit=3, nu_max=1))
    assert result.status == "DOMAIN_STOP" and result.rounds == 2
    assert "切换时" in result.reason
