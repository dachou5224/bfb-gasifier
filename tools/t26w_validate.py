"""T26w：主线基准重跑与取证运行（`.scratch/t26v/evidence.json`）的对照。

不打任何 Corella 补丁：950 °C 已从速率式移到终解检查。跨环境比较用 STATUS“已决定”第 2 条
（`|a − b| ≤ 1e-6·max(|a|, |b|)`），不要求逐位、不要求轮数相同。

``python -m tools.t26w_validate [--limit-s …]``
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hamel.solver.driver as driver
from hamel.cases.config import SolverConfig

_ROOT = Path(__file__).resolve().parents[1]
_CASE, _EVIDENCE = _ROOT / "data/cases/vtt_pressurized.json", _ROOT / ".scratch/t26v/evidence.json"
_OUT = _ROOT / ".scratch/t26w/baseline.json"
_ROUNDS, _LIMIT, _RTOL = 26, 8.0 * 3600.0, 1.0e-6


def _relative(left, right) -> float:
    """跨环境判据：相对差；两侧都极小时看绝对差。STATUS“已决定”第 2 条。"""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    assert left.shape == right.shape, (left.shape, right.shape)
    difference = np.abs(left - right)
    scale = np.maximum(np.abs(left), np.abs(right))
    tiny = scale < 1.0e-12
    return float(np.max(np.where(tiny, difference / 1.0e-15, difference / np.maximum(scale, 1.0e-300))))


def _peak_k(row) -> float:
    """取证文件与本次记录的键名不同：取证用 ``bed_peak_c``，本次用 ``peak_k``。"""
    return float(row["peak_k"]) if "peak_k" in row else float(row["bed_peak_c"]) + 273.15


def compare(record, evidence) -> dict:
    rows = evidence["rounds"]
    out = {
        "final_bed_relative": _relative(record["final_bed_k"], rows[-1]["bed_temperature_k"]),
        "trial_bottom_c": None if record["trial_bottom_peak_k"] is None else record["trial_bottom_peak_k"] - 273.15,
        "evidence_trial_bottom_c": evidence.get("trial_bottom_max_c"),
        "evidence_rounds": evidence.get("rounds_done"), "evidence_iterations": evidence.get("iterations"),
    }
    mine = record["bed_peak_per_round"]
    paired, worst = [], 0.0
    for index, row in enumerate(mine[:len(rows)]):
        value = _relative(_peak_k(row), _peak_k(rows[index]))
        worst = max(worst, value)
        paired.append({"round": row["round"], "peak_c": _peak_k(row) - 273.15,
                       "evidence_peak_c": _peak_k(rows[index]) - 273.15, "relative": value})
    out.update(per_round=paired, per_round_worst_relative=worst, rounds_covered=len(paired),
               over_950c_rounds=[row["round"] for row in mine if _peak_k(row) > 950.0 + 273.15])
    out["passed"] = (record["status"] == "CONVERGED" and out["final_bed_relative"] <= _RTOL
                     and worst <= _RTOL and len(paired) == len(rows))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit-s", type=float, default=_LIMIT)
    args = parser.parse_args(argv)
    started, last = time.monotonic(), {"bed_k": None}
    peaks: list = []
    real_monitor = driver._monitor

    def monitor(state, frozen, assembly=None):
        last["bed_k"] = np.asarray(state.bed_temperature_k, dtype=float).tolist()
        peaks.append(max(last["bed_k"]))
        print(f"[round {len(peaks)}] bed peak {peaks[-1] - 273.15:.4f} °C "
              f"t={time.monotonic() - started:.0f} s", flush=True)
        return real_monitor(state, frozen, assembly)

    driver._monitor = monitor
    try:
        result = driver.solve_case(str(_CASE), "1", config=SolverConfig(outer_round_limit=_ROUNDS),
                                   deadline=lambda: time.monotonic() - started >= args.limit_s)
    finally:
        driver._monitor = real_monitor
    diagnostics = result.diagnostics or {}
    record = {"status": result.status, "reason": result.reason, "rounds": result.rounds,
              "iterations": result.iterations, "elapsed_s": time.monotonic() - started,
              "final_bed_k": last["bed_k"],
              "bed_peak_per_round": diagnostics.get("bed_peak_per_round") or [],
              "trial_bottom_peak_k": diagnostics.get("trial_bottom_peak_k"),
              "bed_temperature_peak": diagnostics.get("bed_temperature_peak")}
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    # 先落盘再比较：比较出错时原始结果仍保留（T26w 的第一次运行因此在最后一步丢了结果）。
    _OUT.write_text(json.dumps(record, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    if _EVIDENCE.exists() and record["final_bed_k"] is not None:
        record["comparison"] = compare(record, json.loads(_EVIDENCE.read_text(encoding="utf-8")))
        _OUT.write_text(json.dumps(record, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    print(json.dumps({k: record.get(k) for k in ("status", "reason", "rounds", "iterations", "elapsed_s")},
                     ensure_ascii=False), flush=True)
    comparison = record.get("comparison")
    if comparison is not None:
        print(json.dumps({k: comparison[k] for k in
                          ("final_bed_relative", "per_round_worst_relative", "rounds_covered",
                           "over_950c_rounds", "passed")}, ensure_ascii=False), flush=True)
    return 0 if comparison and comparison["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
