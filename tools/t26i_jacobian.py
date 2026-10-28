"""分组差分雅可比与逐列差分的逐位比较。只报告，不跑牛顿。求解手段：ADR 0046。

初值点的冻结量与状态由第 1 轮预计算生成后存于 reports/t26i-jacobian/initial/。
"""
from __future__ import annotations

import json
import pickle
import time
from pathlib import Path

import numpy as np

from hamel.cases.assemble import assemble_case
from hamel.solver.driver import _feasible, _references, _scales
from hamel.solver.newton import _grouped_jacobian, _jacobian
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import jacobian_layout, unknown_slots

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "t26i-jacobian"
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_S13 = _ROOT / "reports" / "t26c-diagnosis" / "S13"


def _measure(frozen, x) -> dict:
    n_ref, _m, _q = _references(assemble_case(_CASE, "1").physical)
    scales = _scales(frozen, n_ref, float(frozen.geometric_mass_kg.sum()))
    layout = jacobian_layout(unknown_slots(frozen), frozen)
    calls = {"n": 0}

    def evaluate(trial):
        calls["n"] += 1
        return evaluate_residual(np.asarray(trial, dtype=float), frozen).values

    residual = evaluate(x)
    feasible = lambda trial: _feasible(trial, frozen)
    started, offset = time.perf_counter(), calls["n"]
    grouped, stats = _grouped_jacobian(x, residual, evaluate, scales, feasible, layout)
    grouped_seconds, grouped_calls = time.perf_counter() - started, calls["n"] - offset
    started, offset = time.perf_counter(), calls["n"]
    reference = _jacobian(x, residual, evaluate, scales, feasible)
    column_seconds, column_calls = time.perf_counter() - started, calls["n"] - offset
    return {
        "unknowns": int(x.size), "groups": int(stats["groups"]), "separate": int(stats["separate"]),
        "fallbacks": int(stats["fallbacks"]), "grouped_calls": int(grouped_calls), "column_calls": int(column_calls),
        "grouped_seconds": grouped_seconds, "column_seconds": column_seconds,
        "bitwise_equal": bool(np.array_equal(grouped, reference)),
        "max_abs_difference": float(np.max(np.abs(grouped - reference))),
    }


def run() -> None:
    with (_OUT / "initial" / "frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    with np.load(_OUT / "initial" / "state.npz") as saved:
        x = np.array(saved["x"])
    with (_S13 / "frozen.pkl").open("rb") as handle:
        s13_frozen = pickle.load(handle)
    with np.load(_S13 / "failure.npz", allow_pickle=True) as saved:
        s13_x = np.array(saved["x"])
    records = {"initial": _measure(frozen, x), "S13": _measure(s13_frozen, s13_x)}
    (_OUT / "report.json").write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(records, ensure_ascii=False))


if __name__ == "__main__":
    run()
