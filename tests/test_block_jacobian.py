"""分组差分雅可比。求解手段：ADR 0046。"""
import pickle
from pathlib import Path

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.solver.driver import _feasible
from hamel.solver.gates import OuterMonitor, ResidualBlocks
from hamel.solver.newton import _SQRT_EPS, _grouped_jacobian, _jacobian, damped_newton
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import jacobian_layout, unknown_slots
from test_assembly import make_frozen


def _toy(layout, hidden=0.0):
    """节点邻域内的非线性玩具残差。hidden 是未登记的非相邻耦合。"""
    nodes = layout[0]
    by_node = {int(node): np.flatnonzero(nodes == node) for node in np.unique(nodes)}

    def evaluate(x):
        out = 0.4 * x + 0.05 * x * x + 1.0e-3
        for node, index in by_node.items():
            for other in (node - 1, node + 1):
                neighbour = by_node.get(other)
                if neighbour is not None:
                    for position, column in enumerate(index):
                        if position < neighbour.size:
                            out[column] += 0.07 * x[neighbour[position]]
        if hidden:
            out[0] += hidden * x[-1]
        return out

    return evaluate


def _case(frozen):
    slots = unknown_slots(frozen)
    layout = jacobian_layout(slots, frozen)
    x = np.linspace(0.5, 1.5, len(slots))
    return slots, layout, x, np.linspace(1.0, 2.0, len(slots))


def test_grouped_matches_columnwise_on_a_toy_problem():
    frozen = make_frozen()
    slots, layout, x, scales = _case(frozen)
    evaluate = _toy(layout)
    residual = evaluate(x)
    grouped, stats = _grouped_jacobian(x, residual, evaluate, scales, lambda trial: True, layout)
    reference = _jacobian(x, residual, evaluate, scales, lambda trial: True)
    assert np.array_equal(grouped, reference)
    assert stats["fallbacks"] == 0 and stats["separate"] == 0
    assert 0 < stats["groups"] < x.size


def test_unregistered_non_adjacent_coupling_stops_with_the_rows_and_columns():
    frozen = make_frozen()
    slots, layout, x, scales = _case(frozen)
    evaluate = _toy(layout, hidden=0.2)
    journal = {"names": [slot.name for slot in slots]}
    result = damped_newton(
        x, evaluate, scales, lambda residual: False, lambda trial: True,
        lambda residual: float(np.linalg.norm(residual)), SolverConfig(),
        journal=journal, jacobian_layout=layout, verify_jacobian=True,
    )
    assert result.status == "JACOBIAN_MISMATCH"
    assert journal["jacobian_mismatch"] == {"rows": [0], "columns": [x.size - 1]}


def test_infeasible_group_falls_back_to_single_columns():
    frozen = make_frozen()
    slots = unknown_slots(frozen)
    layout = jacobian_layout(slots, frozen)
    x, scales = np.ones(len(slots)), np.ones(len(slots))
    nodes, within, _ = layout
    first = int(np.flatnonzero(nodes == 0)[0])
    third = int(np.flatnonzero((nodes == 3) & (within == within[first]))[0])
    threshold = 2.0 + 1.5 * _SQRT_EPS
    feasible = lambda trial: float(trial[first] + trial[third]) <= threshold
    evaluate = _toy(layout)
    residual = evaluate(x)
    grouped, stats = _grouped_jacobian(x, residual, evaluate, scales, feasible, layout)
    assert stats["fallbacks"] >= 1
    assert np.array_equal(grouped, _jacobian(x, residual, evaluate, scales, feasible))


def test_registered_non_adjacent_columns_are_separate_groups():
    frozen = make_frozen(hover_closure="B")
    slots, layout, x, scales = _case(frozen)
    nodes = layout[0]
    base = _toy(layout)

    def evaluate(values):
        out = base(values)
        out[int(np.flatnonzero(nodes == 1)[0])] += 0.3 * values[int(np.flatnonzero(nodes == 3)[-1])]
        return out

    residual = evaluate(x)
    grouped, stats = _grouped_jacobian(x, residual, evaluate, scales, lambda trial: True, layout)
    reference = _jacobian(x, residual, evaluate, scales, lambda trial: True)
    assert np.array_equal(grouped, reference)
    assert stats["separate"] == int(layout[2].sum()) > 0


def _six_nodes(hidden=0.0, starts=None):
    """6 节点、每节点 1 未知量。hidden 是未登记的非相邻耦合 F_0 += hidden·x_5。"""
    nodes = np.arange(6)
    layout = (nodes, np.zeros(6, dtype=int), np.zeros(6, dtype=bool))
    blocks = ResidualBlocks(nodes, np.array([], dtype=int), np.array([], dtype=int), 1.0, 1.0, 1.0)
    journal = {"iterations": [], "blocks": blocks, "names": [f"x{i}" for i in range(6)], "jacobian_layout": layout}
    plan = starts or {1: np.zeros(6), 2: np.ones(6)}

    def evaluate(trial, _opened=None):
        out = np.asarray(trial, dtype=float).copy()
        out[0] += hidden * float(out[5])
        return out

    opened = lambda x, _i, _c: type("Opened", (), {"x": np.asarray(x, dtype=float).copy(), "info": {}})()
    checked = lambda values, _i: type("Checked", (), {"residual": evaluate(values), "blocks": blocks, "info": {}})()
    return layout, journal, plan, evaluate, opened, checked, blocks


def _aligned(plan, opened, evaluate, checked, blocks, config, journal):
    from hamel.solver.driver import solve_aligned
    return solve_aligned(
        lambda x, iteration, consistency: opened(plan[iteration], iteration, consistency), evaluate,
        lambda _opened: (blocks, np.ones(6), lambda trial: bool(np.all(np.isfinite(trial))), (1.0, 1.0, 1.0, 1.0)),
        lambda _opened, values: OuterMonitor(np.ones(1), np.ones(1), np.ones(1), np.asarray(values)[:1]),
        checked, config, journal=journal,
    )


def test_verification_waits_for_the_first_constructed_jacobian():
    """外层第 1 轮 0 次内层迭代不消耗核对。ADR 0046 第 5 项。"""
    _layout, journal, plan, evaluate, opened, checked, blocks = _six_nodes(hidden=0.2)
    result = _aligned(plan, opened, evaluate, checked, blocks, SolverConfig(), journal)
    assert result.status == "JACOBIAN_MISMATCH"
    assert result.rounds == 2
    assert journal["jacobian_mismatch"] == {"rows": [0], "columns": [5]}
    assert "jacobian_verified" not in journal


def test_verification_runs_once_per_run(monkeypatch):
    from hamel.solver import newton
    calls = {"full": 0}
    reference = newton._jacobian

    def counted(x, residual, evaluate, scales, feasible, *, blocks=None, columns=None):
        calls["full"] += columns is None
        return reference(x, residual, evaluate, scales, feasible, blocks=blocks, columns=columns)

    monkeypatch.setattr(newton, "_jacobian", counted)
    _layout, journal, plan, evaluate, opened, checked, blocks = _six_nodes(starts={1: np.ones(6), 2: np.ones(6)})
    result = _aligned(plan, opened, evaluate, checked, blocks, SolverConfig(), journal)
    assert result.status == "CONVERGED" and result.rounds == 2
    assert journal["jacobian_verified"] is True
    assert calls["full"] == 1


@pytest.mark.slow
def test_case1_saved_state_grouped_matches_columnwise():
    root = Path(__file__).resolve().parents[1]
    folder = root / "reports/t26c-diagnosis/baseline"
    with (folder / "frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    with np.load(folder / "failure.npz", allow_pickle=True) as saved:
        x, scales = np.array(saved["x"]), np.array(saved["scales"])
    slots = unknown_slots(frozen)
    layout = jacobian_layout(slots, frozen)
    evaluate = lambda trial: evaluate_residual(np.asarray(trial, dtype=float), frozen).values
    residual = evaluate(x)
    feasible = lambda trial: _feasible(trial, frozen)
    grouped, stats = _grouped_jacobian(x, residual, evaluate, scales, feasible, layout)
    reference = _jacobian(x, residual, evaluate, scales, feasible)
    assert np.array_equal(grouped, reference)
    assert stats["groups"] < x.size
