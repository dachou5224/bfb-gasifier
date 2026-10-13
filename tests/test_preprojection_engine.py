"""OPT-002 preprojection_engine 单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.preprojection_engine import (
    LsqPreprojectionSpec,
    PreprojectionBcContext,
    run_improvement_gated_lsq,
    run_lsq_preprojection,
    run_relax_schedule,
    scaled_rms,
)


def test_scaled_rms_basic():
    v = np.array([3.0, 4.0], dtype=np.float64)
    assert scaled_rms(v, 5.0) == pytest.approx(np.sqrt(0.5))


def test_run_improvement_gated_lsq_accepts_on_improvement():
    state = {"x": 2.0}
    calls = {"restore": 0, "sync": 0, "finalize": 0}

    def restore() -> None:
        calls["restore"] += 1
        state["x"] = 2.0

    def set_values(values: np.ndarray) -> None:
        state["x"] = float(values[0])

    def objective(values: np.ndarray) -> np.ndarray:
        return np.array([float(values[0])], dtype=np.float64)

    bc = PreprojectionBcContext(
        cell_index=0,
        apply_bc_fn=lambda: calls.__setitem__("finalize", calls["finalize"] + 1),
    )

    ok = run_improvement_gated_lsq(
        spec=LsqPreprojectionSpec(
            x0=np.array([2.0]),
            lower=np.array([0.0]),
            upper=np.array([10.0]),
            max_nfev=8,
        ),
        set_values=set_values,
        objective=objective,
        initial_metric=2.0,
        final_metric=lambda: abs(state["x"]),
        improvement_factor=0.75,
        restore=restore,
        sync_full_residuals=lambda: calls.__setitem__("sync", calls["sync"] + 1),
        bc=bc,
    )
    assert ok is True
    assert state["x"] < 0.75 * 2.0
    assert calls["sync"] == 1
    assert calls["restore"] == 0


def test_run_improvement_gated_lsq_rejects_and_restores():
    state = {"x": 5.0}
    restored = {"flag": False}

    def restore() -> None:
        state["x"] = 5.0
        restored["flag"] = True

    def set_values(values: np.ndarray) -> None:
        state["x"] = float(values[0])

    def objective(values: np.ndarray) -> np.ndarray:
        return np.array([float(values[0]) + 1.0], dtype=np.float64)

    sync_count = {"n": 0}
    ok = run_improvement_gated_lsq(
        spec=LsqPreprojectionSpec(
            x0=np.array([5.0]),
            lower=np.array([4.0]),
            upper=np.array([6.0]),
            max_nfev=4,
        ),
        set_values=set_values,
        objective=objective,
        initial_metric=1.0,
        final_metric=lambda: 1.0,
        improvement_factor=0.5,
        restore=restore,
        sync_full_residuals=lambda: sync_count.__setitem__("n", sync_count["n"] + 1),
        bc=None,
    )
    assert ok is False
    assert restored["flag"] is True
    assert sync_count["n"] == 1


def test_run_improvement_gated_lsq_lsq_exception_restores():
    restored = {"flag": False}

    def restore() -> None:
        restored["flag"] = True

    def set_values(_values: np.ndarray) -> None:
        raise RuntimeError("bad set")

    def objective(_values: np.ndarray) -> np.ndarray:
        return np.array([0.0])

    ok = run_improvement_gated_lsq(
        spec=LsqPreprojectionSpec(
            x0=np.array([1.0]),
            lower=np.array([0.0]),
            upper=np.array([2.0]),
            max_nfev=2,
        ),
        set_values=set_values,
        objective=objective,
        initial_metric=1.0,
        final_metric=lambda: 0.0,
        improvement_factor=0.75,
        restore=restore,
        sync_full_residuals=lambda: None,
        bc=None,
    )
    assert ok is False
    assert restored["flag"] is True


def test_run_relax_schedule_accepts_first_success():
    state = {"x": 0}
    log: list[float] = []

    def restore() -> None:
        state["x"] = 0

    def try_trial(omega: float) -> bool:
        log.append(omega)
        state["x"] = int(omega * 10)
        return omega <= 0.5

    ok, value = run_relax_schedule(
        [1.0, 0.5, 0.25],
        restore=restore,
        prepare_trial=lambda: None,
        try_trial=try_trial,
    )
    assert ok is True
    assert value == 0.5
    assert log == [1.0, 0.5]


def test_run_relax_schedule_failure_restores():
    restored = {"n": 0}

    def restore() -> None:
        restored["n"] += 1

    ok, value = run_relax_schedule(
        [1.0, 0.5],
        restore=restore,
        try_trial=lambda _omega: False,
        finalize_failure=lambda: restored.__setitem__("n", restored["n"] + 10),
    )
    assert ok is False
    assert value == 0.0
    # 每个 candidate 前 restore 一次 + 失败后再 restore 一次
    assert restored["n"] == 13


def test_run_lsq_preprojection_custom_accept():
    state = {"x": 5.0}

    def restore() -> None:
        state["x"] = 5.0

    def apply_solution(raw: np.ndarray) -> None:
        state["x"] = float(raw[0])

    def objective(x: np.ndarray) -> np.ndarray:
        return np.array([float(x[0])], dtype=np.float64)

    ok = run_lsq_preprojection(
        spec=LsqPreprojectionSpec(
            x0=np.array([5.0]),
            lower=np.array([0.0]),
            upper=np.array([10.0]),
            max_nfev=8,
        ),
        objective=objective,
        apply_solution=apply_solution,
        restore=restore,
        accept_predicate=lambda: state["x"] < 1.0,
        sync_full_residuals=lambda: None,
    )
    assert ok is True
    assert state["x"] < 1.0
