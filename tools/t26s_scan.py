"""T26s 端点、扫描与外层。插值 (1−s)·x8+s·x9：ReactorState 全部数组，以及 reaction_solid_kg_s、
previous_drain_kg_s、previous_fuel_return_kg_s、previous_release 的数值字段、redirected_return_kg_s、
redirected_fuel_out_kg_s、previous_fine_share。其余离散量固定取第 8 轮。形状不同或非有限即停止。
不记录出口组成、T̄_FRB、X_C。Corella 只在本工具外推。
"""
from __future__ import annotations

import argparse, dataclasses, json, pickle, signal, subprocess, time
from pathlib import Path
import numpy as np
import hamel.freeboard.accounting as accounting, hamel.kinetics.homogeneous as homogeneous
import hamel.solver.driver as driver, hamel.solver.precalc as precalc
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.solver.gates import _abs_max
from hamel.solver.newton import SolveTimeout
from hamel.solver.vector import ReactorState

_ROOT = Path(__file__).resolve().parents[1]
_SCRATCH = _ROOT / ".scratch/t26s"
_CASE = _ROOT / "data/cases/vtt_pressurized.json"
_BASE = "4ab8555895bf566e5a159d693268ecb2f9446ca7"
_NUMERIC = {"reaction_solid_kg_s", "previous_drain_kg_s", "previous_fuel_return_kg_s", "redirected_return_kg_s", "redirected_fuel_out_kg_s", "previous_fine_share"}

class _Timeout(BaseException):
    """墙钟到点。CHARTER §7。"""

def _write(name, payload):
    (_SCRATCH / name).write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=str), encoding="utf-8")

def _load(name):
    with (_SCRATCH / name).open("rb") as handle:
        return pickle.load(handle)

def _blend(left, right, s, label):
    if left is None and right is None:
        return None
    if isinstance(left, dict) or isinstance(right, dict):
        if set(left) != set(right):
            raise SystemExit(f"字典键不同：{label}")
        return {key: _blend(left[key], right[key], s, f"{label}.{key}") for key in left}
    if dataclasses.is_dataclass(left) and not isinstance(left, type):
        if type(left) is not type(right):
            raise SystemExit(f"类型不同：{label}")
        values = {item.name: _blend(getattr(left, item.name), getattr(right, item.name), s, f"{label}.{item.name}") for item in dataclasses.fields(left)}
        return type(left)(**values)
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        first, second = np.asarray(left, float), np.asarray(right, float)
        if first.shape != second.shape or not np.isfinite(first).all() or not np.isfinite(second).all():
            raise SystemExit(f"形状不同或非有限：{label}")
        return (1.0 - s) * first + s * second
    if isinstance(left, (bool, np.bool_)) or isinstance(right, (bool, np.bool_)):
        raise SystemExit(f"布尔值不参与插值：{label}")
    try:
        first, second = float(left), float(right)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"不能插值：{label}") from exc
    if not np.isfinite([first, second]).all():
        raise SystemExit(f"非有限值：{label}")
    return (1.0 - s) * first + s * second

def _same(left, right):
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        return np.array_equal(np.asarray(left), np.asarray(right))
    if dataclasses.is_dataclass(left) and not isinstance(left, type):
        return type(left) is type(right) and all(_same(getattr(left, item.name), getattr(right, item.name)) for item in dataclasses.fields(left))
    return left == right

def _call(saved8, saved9, s):
    a8, k8, a9, k9 = saved8["args"], saved8["kwargs"], saved9["args"], saved9["kwargs"]
    if set(k8) != set(k9):
        raise SystemExit(f"关键字集合不同：{sorted(set(k8) ^ set(k9))}")
    state = ReactorState(**{item.name: _blend(getattr(a8[2], item.name), getattr(a9[2], item.name), s, item.name) for item in dataclasses.fields(ReactorState)})
    kwargs, discrete = {}, []
    for key, value in k8.items():
        numeric = key in _NUMERIC or key == "previous_release"
        kwargs[key] = _blend(value, k9[key], s, key) if numeric else value
        if not numeric and not _same(value, k9[key]):
            discrete.append(key)
    for index, name in ((0, "path"), (1, "case_id")):
        if not _same(a8[index], a9[index]):
            discrete.append(name)
    return (a8[0], a8[1], state), kwargs, discrete

def _numeric(frozen):
    values = {}
    for item in dataclasses.fields(frozen):
        value = getattr(frozen, item.name)
        if isinstance(value, np.ndarray):
            values[item.name] = np.array(value)
        elif isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, (bool, np.bool_)):
            values[item.name] = np.asarray(value)
    return values

def _compare(got, path):
    snap = np.load(path)
    return {
        "missing": sorted(set(snap.files) - set(got)),
        "extra": sorted(set(got) - set(snap.files)),
        "different": sorted(name for name in snap.files if name in got and not np.array_equal(got[name], snap[name])),
    }

def _config(saved, u_min):
    if not np.isfinite(u_min) or float(u_min) < 0.0:
        raise SystemExit("u_min 必须是有限非负数")
    names = {item.name for item in dataclasses.fields(SolverConfig)}
    data = {item.name: getattr(saved, item.name) for item in dataclasses.fields(saved) if item.name in names and item.name != "freeboard_residence_floor_m_s"}
    if "freeboard_residence_floor_m_s" in names:
        data["freeboard_residence_floor_m_s"] = float(u_min)
    return SolverConfig(**data)

def _run(args, kwargs, *, samples=None, trace=False):
    books, box = [], []
    real_acc, real_tr, real_samples = precalc.account_freeboard, accounting.trace_until_fate, precalc.within_class_samples

    def account(*pos, **kw):
        box.clear()
        book = real_acc(*pos, **kw)
        use_bed, kept = None, []
        if trace:
            shares = np.asarray(pos[8], float)
            rows = kw.get("samples") or ((float(pos[4]), float(pos[6]), 1.0),)
            cursor = 0
            bed_res, bed_out = np.zeros(len(book.cells)), np.zeros(len(book.cells))
            for _d, _phi, weight in rows:
                for index in range(shares.size):
                    rec = box[cursor]
                    cursor += 1
                    scale = float(weight) * float(shares[index])
                    bed_res += scale * rec[5]
                    bed_out += scale * (rec[6] + rec[7])
                cursor += len(book.cells)
            use_bed = (bed_res > 0.0) | (bed_out > 0.0)
            if samples is None or len(books) in (2, 3):
                kept = [(rec[0], rec[1], rec[2], rec[3], rec[4], rec[8], rec[9], rec[10]) for rec in box]
        books.append({"book": book, "use_bed": use_bed, "traces": kept})
        return book

    def follow(height, speed, *rest, **_kw):
        fate = real_tr(height, speed, *rest, **_kw)
        if trace:
            box.append((float(height), float(speed), fate.fate, fate.entered_cell, fate.entry_kind, np.asarray(fate.residence_s, float), np.asarray(fate.upward_free, float), np.asarray(fate.downward_free, float), float(rest[4]), float(rest[6]), tuple(map(float, fate.residence_s))))
        return fate

    def sampled(lower, upper, diameters=8, sphericities=8, sigma=accounting.SPHERICITY_SIGMA):
        count = samples if samples is not None else diameters
        return real_samples(lower, upper, diameters=count, sphericities=count, sigma=sigma)

    precalc.account_freeboard, accounting.trace_until_fate = account, follow
    if samples is not None:
        precalc.within_class_samples = sampled
    started = time.perf_counter()
    try:
        result = precalc.precalculate(*args, **kwargs)
    finally:
        precalc.account_freeboard, accounting.trace_until_fate, precalc.within_class_samples = real_acc, real_tr, real_samples
    return result, books, time.perf_counter() - started

def _cells(books, classes):
    rows = []
    for klass in classes:
        item, cells = books[klass], []
        for index, cell in enumerate(item["book"].cells):
            cells.append({"cell": index, "theta_up": cell.theta_up, "theta_down": cell.theta_down, "m_hat": cell.m_hat, "a_f": cell.a_f, "b_f": cell.b_f, "c_up": cell.c_up, "c_down": cell.c_down, "free_state": cell.free_state, "use_bed": None if item["use_bed"] is None else bool(item["use_bed"][index]), "windows": [window.fate for window in cell.windows]})
        rows.append({"class": klass, "capture": list(item["book"].capture_cells), "cells": cells})
    return rows

def _solve(rounds, limit_s, save):
    original, started, seen = homogeneous._CORELLA_HIGH_K, time.monotonic(), {}
    homogeneous._CORELLA_HIGH_K = float("inf")
    real, monitor, rows, monitors = driver.precalculate, driver._monitor, [], []
    folder = _SCRATCH / ("capture" if save else "outer")
    folder.mkdir(parents=True, exist_ok=True)

    def precalculate(*args, **kwargs):
        iteration, consistency = kwargs.get("iteration"), bool(kwargs.get("consistency", False))
        print(f"precalc round={iteration} consistency={consistency}", flush=True)
        if iteration in save and not consistency:
            seen[iteration] = seen.get(iteration, 0) + 1
            if seen[iteration] == 1:
                (_SCRATCH / f"precalc_args_r{iteration}.pkl").write_bytes(pickle.dumps({"args": args, "kwargs": kwargs}, protocol=pickle.HIGHEST_PROTOCOL))
                print(f"saved round {iteration}", flush=True)
        return real(*args, **kwargs)

    def watching(state, frozen, assembly=None):
        value = monitor(state, frozen, assembly)
        monitors.append((value, float(np.sum(frozen.geometric_mass_kg))))
        top = float(np.sum(frozen.a_free_per_s[0] * state.freeboard_ash_free_kg[0]))
        if frozen.hover_closure == "C" and assembly is not None:
            top += float(assembly.redirect_ash_kg_s)
        inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
        rows.append({"phi": float(frozen.fly_ash_fraction), "w_f": float(frozen.bed_ash_fine_share), "E": top, "ash_target": float(frozen.ash_target_kg_s), "ash_relative": value.ash_relative, "bed_char_kg": float(np.sum(state.bed_char_kg)), "min_inert_kg": float(np.min(inert)), "bed_temperature_k": np.asarray(state.bed_temperature_k, float).tolist()})
        _write("solve_progress.json", {"rounds": len(rows), "w_f": rows[-1]["w_f"], "inner_pending": True})
        print(f"round {len(rows)} w_f={rows[-1]['w_f']}", flush=True)
        return value

    driver.precalculate, driver._monitor = precalculate, watching
    signal.signal(signal.SIGALRM, lambda *_a: (_ for _ in ()).throw(_Timeout()))
    signal.setitimer(signal.ITIMER_REAL, float(limit_s))
    summary = {"saved": {}}
    try:
        result = driver.solve_case(str(_CASE), "1", config=SolverConfig(outer_round_limit=int(rounds)), deadline=lambda: time.monotonic() - started >= limit_s, artifact_dir=str(folder))
        summary.update(status=result.status, reason=result.reason, rounds=result.rounds, iterations=result.iterations, saved=dict(seen))
        if result.ratios:
            summary["block_ratios"] = {key: float(value) for key, value in result.ratios.items()}
        summary["consistent_ash_relative"] = (result.diagnostics or {}).get("consistent_ash_relative")
        initial = (result.diagnostics or {}).get("initial") or {}
        n_ref, m_ref, _q = driver._references(assemble_case(str(_CASE), "1").physical)
        gate_from = 3 if initial.get("phi_init") != initial.get("phi") else 2
        for index, item in enumerate((result.diagnostics or {}).get("rounds") or []):
            if index >= len(rows):
                break
            newton = item.get("newton") or []
            rows[index].update(round=item.get("round"), inner=item.get("inner"), inner_reason=item.get("reason"), inner_iterations=len(newton), projected=int(sum(int(entry.get("projected", 0)) for entry in newton)))
            if newton and isinstance(newton[-1], dict) and "final_check" in newton[-1]:
                rows[index]["final_check"] = newton[-1]["final_check"]
            if index:
                current, m_tot = monitors[index]
                previous = monitors[index - 1][0]
                rows[index]["gate"] = {"dT": _abs_max(current.temperature_k, previous.temperature_k), "solid": _abs_max(current.solid_flow_kg_s, previous.solid_flow_kg_s) / m_ref, "inventory": _abs_max(current.inventory_kg, previous.inventory_kg) / m_tot, "gas": _abs_max(current.gas_flow_mol_s, previous.gas_flow_mol_s) / n_ref, "checked": item.get("round", 0) >= gate_from}
        failure = folder / "failure.npz"
        if failure.exists() and (result.status == "DAMPING_FAILED" or "雅可比" in str(result.reason)):
            data = np.load(failure, allow_pickle=True)
            names = list(map(str, data["names"])) if "names" in data.files else []

            def named(axis):
                return [{"index": int(i), "name": names[i] if i < len(names) else str(i), "y": float(data["y"][i]), "x": float(data["x"][i])} for i in np.where(~np.any(data["J"], axis=axis))[0]]

            summary["failure"] = {"zero_columns": named(0), "zero_rows": named(1)}
    except (_Timeout, SolveTimeout) as exc:
        summary.update(status="TIMEOUT", reason=str(exc) or "运行超过墙钟上限", rounds=len(rows), saved=dict(seen))
    except Exception as exc:  # noqa: BLE001
        summary.update(status=type(exc).__name__, reason=str(exc)[:800], rounds=len(rows), saved=dict(seen))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        driver.precalculate, driver._monitor, homogeneous._CORELLA_HIGH_K = real, monitor, original
    summary["elapsed_s"], summary["rounds_detail"] = time.monotonic() - started, rows
    return summary

def _corella(enter):
    if enter:
        _corella.previous = homogeneous._CORELLA_HIGH_K
        homogeneous._CORELLA_HIGH_K = float("inf")
    else:
        homogeneous._CORELLA_HIGH_K = _corella.previous

def _hamel_clean():
    diff = subprocess.check_output(["git", "diff", "--stat", _BASE, "--", "hamel"], cwd=_ROOT, text=True)
    if diff.strip():
        raise SystemExit(f"hamel/ 已偏离基准 {_BASE}：{diff}")

def endpoints():
    saved8, saved9 = _load("precalc_args_r8.pkl"), _load("precalc_args_r9.pkl")
    _corella(True)
    report = {"baseline_sha": _BASE, "discrete_differences": _call(saved8, saved9, 0.0)[2]}
    try:
        snaps = {}
        for label, saved in (("r8", saved8), ("r9", saved9)):
            result, books, elapsed = _run(saved["args"], saved["kwargs"])
            np.savez(_SCRATCH / f"frozen_{label}.npz", **_numeric(result.frozen))
            snaps[label] = {"w_f": float(result.frozen.bed_ash_fine_share), "class2_bottom_m_hat": books[2]["book"].cells[0].m_hat, "elapsed_s": elapsed}
        report["endpoints"] = snaps
        report["commander"] = {"w_f_s0": 1.8791e-3, "w_f_s1": 1.9310e-3, "class2_bottom_m_hat_s0": 3.25}
        report["endpoint_check"] = {}
        for name, ess, label in (("s0", 0.0, "r8"), ("s1", 1.0, "r9")):
            args, kwargs, _discrete = _call(saved8, saved9, ess)
            result, _books, _elapsed = _run(args, kwargs)
            report["endpoint_check"][name] = _compare(_numeric(result.frozen), _SCRATCH / f"frozen_{label}.npz")
    finally:
        _corella(False)
    report["stop"] = any(report["endpoint_check"][name][key] for name in report["endpoint_check"] for key in ("missing", "extra", "different"))
    _write("endpoints.json", report)
    if report["stop"]:
        raise SystemExit("端点复现有不等，停下")
    return report

def capture(limit_s):
    _SCRATCH.mkdir(parents=True, exist_ok=True)
    _hamel_clean()
    summary = _solve(9, limit_s, {8, 9})
    summary["baseline_sha"] = _BASE
    _write("capture.json", summary)
    saved = summary.get("saved") or {}
    if saved.get(8) != 1 or saved.get(9) != 1:
        raise SystemExit(f"第 9 轮之前停止或预计算次数异常：{summary.get('status')} {summary.get('reason')} saved={saved}")
    return endpoints()

def compare():
    _corella(True)
    report = {}
    try:
        for label in ("r8", "r9"):
            saved = _load(f"precalc_args_{label}.pkl")
            kwargs = dict(saved["kwargs"])
            kwargs["solver"] = _config(kwargs["solver"], 0.0)
            result, _books, elapsed = _run(saved["args"], kwargs)
            report[label] = {**_compare(_numeric(result.frozen), _SCRATCH / f"frozen_{label}.npz"), "elapsed_s": elapsed}
    finally:
        _corella(False)
    _write("compare.json", report)
    if any(report[label][key] for label in report for key in ("missing", "extra", "different")):
        raise SystemExit("第 C 项有不等")

def initial():
    _corella(True)
    packed = {}
    try:
        for u_min in (0.0, 1.0e-2):
            result, books, elapsed = _run((str(_CASE), "1", None), {"solver": _config(SolverConfig(), u_min), "iteration": 1, "consistency": False, "fly_ash_fraction": 0.3}, trace=True)
            packed[u_min] = (books, elapsed, float(result.frozen.bed_ash_fine_share))
    finally:
        _corella(False)
    base, new = packed[0.0], packed[1.0e-2]
    shortened = []
    for klass, (left, right) in enumerate(zip(base[0], new[0])):
        same = len(left["traces"]) == len(right["traces"])
        shortened.append({"class": klass, "traces": len(left["traces"]), "shortened": None if not same else sum(a[7] != b[7] for a, b in zip(left["traces"], right["traces"]))})
    _write("initial.json", {"phi_init": 0.3, "note": "第一次预计算与驱动一致，φ_init = min(φ, 0.3)", "shortened": shortened, "elapsed_s": {"0": base[1], "0.01": new[1]}, "w_f": {"0": base[2], "0.01": new[2]}, "classes_0": _cells(base[0], range(len(base[0]))), "classes_0.01": _cells(new[0], range(len(new[0])))})

def scan():
    saved8, saved9 = _load("precalc_args_r8.pkl"), _load("precalc_args_r9.pkl")
    _corella(True)
    path = _SCRATCH / "scan.jsonl"
    try:
        for u_min in (0.0, 1.0e-2):
            for count in (8, 16):
                points = [index / 10.0 for index in range(11)]
                if count == 8:
                    points += [index / 10.0 + 0.05 for index in range(10)]
                for ess in points:
                    args, kwargs, _discrete = _call(saved8, saved9, ess)
                    kwargs["solver"] = _config(kwargs["solver"], u_min)
                    result, books, elapsed = _run(args, kwargs, samples=count, trace=True)
                    traces = {str(klass): [{"h": rec[0], "u": rec[1], "fate": rec[2], "entered": rec[3], "kind": rec[4], "d": rec[5], "phi": rec[6]} for rec in books[klass]["traces"]] for klass in (2, 3)}
                    row = {"u_min": u_min, "n": count, "s": ess, "w_f": float(result.frozen.bed_ash_fine_share), "elapsed_s": elapsed, "classes": _cells(books, (1, 2, 3, 4)), "traces": traces}
                    with path.open("a", encoding="utf-8") as handle:
                        handle.write(json.dumps(row, ensure_ascii=False, default=str) + "\n")
                    print(f"scan u={u_min} n={count} s={ess} {elapsed:.1f}s", flush=True)
    finally:
        _corella(False)

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "endpoints", "compare", "initial", "scan", "outer"))
    parser.add_argument("--limit-s", type=float, default=None)
    args = parser.parse_args(argv)
    _SCRATCH.mkdir(parents=True, exist_ok=True)
    if args.mode == "capture":
        capture(4.0 * 3600.0 if args.limit_s is None else args.limit_s)
    elif args.mode == "endpoints":
        endpoints()
    elif args.mode == "compare":
        compare()
    elif args.mode == "initial":
        initial()
    elif args.mode == "scan":
        scan()
    else:
        summary = _solve(26, 6.0 * 3600.0 if args.limit_s is None else args.limit_s, set())
        summary["baseline_sha"] = _BASE
        _write("outer.json", summary)
        if summary.get("failure") or summary.get("status") in {"DAMPING_FAILED", "TIMEOUT"} or "雅可比" in str(summary.get("reason")):
            raise SystemExit(f"外层停止：{summary.get('status')} {summary.get('reason')}")

if __name__ == "__main__":
    main()
