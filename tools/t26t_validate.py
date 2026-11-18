"""T26t 的 A、C、D 项入口；插值、比较与外层包装复用 `tools/t26s_scan.py`。Corella 只在本工具外推，不记录出口组成、T̄_FRB、X_C。"""
from __future__ import annotations

import argparse
import dataclasses
import json
import pickle
from pathlib import Path

import numpy as np

import t26s_scan as scan
import hamel.solver.driver as driver
from hamel.cases.config import SolverConfig

_ROOT = Path(__file__).resolve().parents[1]
_SCRATCH = _ROOT / ".scratch/t26t"
_CASE = _ROOT / "data/cases/vtt_pressurized.json"
_BASE = "b87fafb48da9f2a236745cbf5800dd0fb414d993"
scan._SCRATCH, scan._BASE = _SCRATCH, _BASE


def _write(name, payload):
    (_SCRATCH / name).write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=str), encoding="utf-8")


def _solver(saved):
    """旧进程 pickle 的 SolverConfig 没有 ADR 0058 的新字段，按字段重建。"""
    return SolverConfig(**{item.name: getattr(saved, item.name) for item in dataclasses.fields(saved) if item.name in {f.name for f in dataclasses.fields(SolverConfig)}})


_CLASS_FIELDS = ("a_free_per_s", "a_hover_per_s", "b_free_per_s", "b_hover_per_s", "c_up_per_s", "c_down_per_s", "theta_up", "theta_down", "theta", "drift_hover_per_s", "drift_free_per_s", "migration_mu_per_s", "ejection_to_free", "ejection_to_hover", "fuel_feed_kg_s", "drain_char_kg_s", "char_gen_free_kg_s", "char_gen_hover_kg_s", "ash_gen_free_kg_s", "ash_gen_hover_kg_s", "char_gen_bed_kg_s", "ash_gen_bed_kg_s", "fuel_bed_kg", "fuel_up_bed_kg_s", "fuel_down_bed_kg_s", "fuel_exit_bed_kg_s", "fuel_up_free_kg_s", "fuel_down_free_kg_s", "u_r_bed_m_s", "u_r_freeboard_m_s", "u_r_hover_m_s", "hover", "free_absent")
# 只依赖该粒级自身轨迹的量（记账与悬停窗口）。生成项、床区流与排灰经床区总固体流与其它粒级耦合，不属于此类。
_OWN_FIELDS = ("a_free_per_s", "b_free_per_s", "c_up_per_s", "c_down_per_s", "theta", "theta_up", "theta_down", "u_r_freeboard_m_s", "drift_hover_per_s", "drift_free_per_s", "migration_mu_per_s", "hover", "free_absent")


def _class_arrays(values, klass, count):
    return {name: np.asarray(values[name])[..., klass] for name in _CLASS_FIELDS if name in values}


def _relative(left, right):
    left, right = np.asarray(left), np.asarray(right)
    if left.dtype == bool or right.dtype == bool:
        return 0.0 if np.array_equal(left, right) else 1.0
    left, right = left.astype(float), right.astype(float)
    scale = max(float(np.max(np.abs(left), initial=0.0)), float(np.max(np.abs(right), initial=0.0)))
    return 0.0 if scale == 0.0 else float(np.max(np.abs(left - right), initial=0.0) / scale)


def endpoints():
    saved = {label: scan._load(f"precalc_args_{label}.pkl") for label in ("r17", "r18")}
    scan._corella(True)
    report = {"baseline_sha": _BASE, "discrete_differences": scan._call(saved["r17"], saved["r18"], 0.0)[2]}
    try:
        for label in ("r17", "r18"):
            kwargs = dict(saved[label]["kwargs"])
            kwargs["solver"] = _solver(kwargs["solver"])
            result, books, elapsed = scan._run(saved[label]["args"], kwargs)
            np.savez(_SCRATCH / f"frozen_{label}.npz", **scan._numeric(result.frozen))
            report[label] = {"w_f": float(result.frozen.bed_ash_fine_share), "class2_bottom_m_hat": books[2]["book"].cells[0].m_hat, "elapsed_s": elapsed}
        report["check"] = {}
        for name, ess, label in (("s0", 0.0, "r17"), ("s1", 1.0, "r18")):
            args, kwargs, _ = scan._call(saved["r17"], saved["r18"], ess)
            result, _books, _elapsed = scan._run(args, kwargs)
            report["check"][name] = scan._compare(scan._numeric(result.frozen), _SCRATCH / f"frozen_{label}.npz")
    finally:
        scan._corella(False)
    _write("endpoints.json", report)
    if any(report["check"][n][key] for n in report["check"] for key in ("missing", "extra", "different")):
        raise SystemExit("端点复现有不等，停下")

def initial(name):
    scan._corella(True)
    try:
        result, books, elapsed = scan._run(
            (str(_CASE), "1", None),
            {"solver": SolverConfig(), "iteration": 1, "consistency": False, "fly_ash_fraction": 0.3},
            trace=True,
        )
        np.savez(_SCRATCH / f"initial_{name}.npz", **scan._numeric(result.frozen))
        _write(f"initial_{name}.json", {"elapsed_s": elapsed, "w_f": float(result.frozen.bed_ash_fine_share),
                                        "boundary": (result.reported.get("trajectory") or {}).get("boundary"),
                                        "classes": scan._cells(books, range(len(books)))})
    finally:
        scan._corella(False)


def capture(limit_s):
    scan._hamel_clean()
    summary = scan._solve(18, limit_s, {17, 18})
    _write("capture.json", dict(summary, baseline_sha=_BASE))
    saved = summary.get("saved") or {}
    if saved.get(17) != 1 or saved.get(18) != 1:
        raise SystemExit(f"第 18 轮之前停止或预计算次数异常：{summary.get('status')} {summary.get('reason')} saved={saved}")
    endpoints()
    initial("old")


def compare():
    saved = scan._load("precalc_args_r17.pkl")
    kwargs = dict(saved["kwargs"])
    kwargs["solver"] = _solver(kwargs["solver"])
    scan._corella(True)
    try:
        result, _books, elapsed = scan._run(saved["args"], kwargs)
    finally:
        scan._corella(False)
    got, snap = scan._numeric(result.frozen), np.load(_SCRATCH / "frozen_r17.npz")
    uniform = [item[1] == 0 for item in ((result.reported.get("trajectory") or {}).get("boundary") or [])]
    count = got["diameter_m"].size
    report = {"baseline_sha": _BASE, "elapsed_s": elapsed, "uniform": uniform, "classes": {}}
    for klass in range(count):
        rows = {name: _relative(value, snap[name][..., klass]) for name, value in _class_arrays(got, klass, count).items()
                if name in snap.files}
        own = max((value for name, value in rows.items() if name in _OWN_FIELDS), default=0.0)
        report["classes"][klass] = {"uniform": bool(uniform[klass]) if klass < len(uniform) else None,
                                    "own_relative": own, "max_relative": max(rows.values(), default=0.0), "fields": rows}
    _write("compare.json", report)
    bad = [k for k, row in report["classes"].items() if row["uniform"] and row["own_relative"] > 1.0e-12]
    if bad:
        raise SystemExit(f"第 C.2 项不满足：粒级 {bad}")


def _align(rows):
    """按 (d, φ_S, 起点高度) 分组，组内按初速排序取秩：两点间同一轨迹位的归宿可直接对照。"""
    groups = {}
    for row in rows:
        groups.setdefault((round(row[5], 12), round(row[6], 6), round(row[0], 12)), []).append((row[1], f"{row[2]}|{row[3]}|{row[4]}"))
    return {key + (rank,): sig for key, items in groups.items() for rank, (_u, sig) in enumerate(sorted(items))}


def _frozen_points(saved17, saved18):
    points = []
    for index in range(11):
        args, kwargs, _ = scan._call(saved17, saved18, 0.1 * index)
        kwargs = dict(kwargs)
        kwargs["solver"] = _solver(kwargs["solver"])
        result, books, elapsed = scan._run(args, kwargs, samples=8, trace=True)
        traces = {str(k): _align(books[k]["traces"]) for k in (2, 3)}
        points.append({"s": 0.1 * index, "frozen": scan._numeric(result.frozen), "elapsed_s": elapsed,
                       "classes": scan._cells(books, range(len(books))), "traces": traces})
        print(f"scan s={0.1 * index:.1f} {elapsed:.1f}s", flush=True)
    return points




def _jumps(points):
    rows = {}
    for name in points[0]["frozen"]:
        arrays = [np.asarray(point["frozen"][name], dtype=float).ravel() for point in points]
        scale = max(float(np.max(np.abs(array), initial=0.0)) for array in arrays)
        if scale == 0.0:
            continue
        steps = [float(np.max(np.abs(arrays[i + 1] - arrays[i]), initial=0.0) / scale) for i in range(len(arrays) - 1)]
        rows[name] = {"max": max(steps), "absolute": max(steps) * scale, "scale": scale, "at": int(np.argmax(steps))}
    return sorted(rows.items(), key=lambda item: -item[1]["max"])

def scan_twostate():
    saved17, saved18 = scan._load("precalc_args_r17.pkl"), scan._load("precalc_args_r18.pkl")
    scan._corella(True)
    try:
        points = _frozen_points(saved17, saved18)
    finally:
        scan._corella(False)
    coefficients = {}
    for name in ("theta_up", "theta_down", "m_hat", "a_f", "b_f", "c_up", "c_down"):
        values = [[cell.get(name) for item in point["classes"] for cell in item["cells"] if cell.get(name) is not None]
                  for point in points]
        scale = max((abs(value) for row in values for value in row), default=0.0)
        steps = [max((abs(a - b) for a, b in zip(values[i], values[i + 1])), default=0.0) / scale if scale else 0.0
                 for i in range(len(values) - 1)]
        coefficients[name] = {"max": max(steps), "absolute": max(steps) * scale, "at": int(np.argmax(steps))}
    changes = []
    for index in range(10):
        for klass in ("2", "3"):
            left, right = points[index]["traces"][klass], points[index + 1]["traces"][klass]
            changes += [{"s": (0.1 * index, 0.1 * (index + 1)), "class": klass, "key": list(key),
                         "old": left[key], "new": right[key]} for key in sorted(set(left) & set(right)) if left[key] != right[key]]
    _write("scan.json", {"baseline_sha": _BASE, "points": len(points), "jumps": _jumps(points),
                         "coefficients": coefficients, "signature_changes": changes})

def outer(limit_s):
    store, calls = {}, []
    real_solve, real_precalc = driver.solve_case, driver.precalculate

    def solve(*args, **kwargs):
        result = real_solve(*args, **kwargs)
        store["result"] = result
        return result

    def record(*args, **kwargs):
        calls.append((args, kwargs))
        return real_precalc(*args, **kwargs)

    driver.solve_case, driver.precalculate = solve, record
    try:
        summary = scan._solve(26, limit_s, set())
    finally:
        driver.solve_case, driver.precalculate = real_solve, real_precalc
    result = store.get("result")
    if result is not None:
        summary["result"] = {"status": result.status, "reason": result.reason, "rounds": result.rounds,
                             "iterations": result.iterations, "ratios": result.ratios, "balance": result.balance,
                             "nonnegative": result.nonnegative,
                             "bottom_relative": (result.diagnostics or {}).get("bottom_relative"),
                             "negative_floor": (result.diagnostics or {}).get("negative_floor"),
                             "final": (result.diagnostics or {}).get("final")}
        # 外层门通过的那一轮 2(b) 复核的入参（置零前的状态）。
        candidates = [(args, kwargs) for args, kwargs in calls
                      if kwargs.get("iteration") == result.rounds and kwargs.get("consistency")]
        if candidates:
            with (_SCRATCH / "final_args.pkl").open("wb") as handle:
                pickle.dump({"args": candidates[0][0], "kwargs": candidates[0][1]}, handle, protocol=pickle.HIGHEST_PROTOCOL)
            summary["final_args_saved"] = len(candidates)
    _write("outer.json", dict(summary, baseline_sha=_BASE))

def nodebits():
    saved = scan._load("final_args.pkl")
    report = {"baseline_sha": _BASE}
    scan._corella(True)
    try:
        for bisections in (6, 8):
            kwargs = dict(saved["kwargs"])
            kwargs["solver"] = dataclasses.replace(kwargs["solver"], freeboard_boundary_bisections=bisections)
            result, _books, elapsed = scan._run(saved["args"], kwargs)
            report[str(bisections)] = {"frozen": scan._numeric(result.frozen), "elapsed_s": elapsed}
    finally:
        scan._corella(False)
    rows = {name: _relative(report["6"]["frozen"][name], report["8"]["frozen"][name]) for name in report["6"]["frozen"]}
    _write("nodebits.json", {"baseline_sha": _BASE, "max": max(rows.values(), default=0.0), "fields": rows, "elapsed_s": {k: report[k]["elapsed_s"] for k in ("6", "8")}})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "endpoints", "compare", "initial", "scan", "outer", "nodebits"))
    parser.add_argument("--limit-s", type=float, default=None)
    args = parser.parse_args(argv)
    _SCRATCH.mkdir(parents=True, exist_ok=True)
    if args.mode in ("capture", "outer"):
        hours = 2.5 if args.mode == "capture" else 8.0
        (capture if args.mode == "capture" else outer)(hours * 3600.0 if args.limit_s is None else args.limit_s)
    else:
        {"endpoints": endpoints, "compare": compare, "initial": lambda: initial("new"),
         "scan": scan_twostate, "nodebits": nodebits}[args.mode]()


if __name__ == "__main__":
    main()
