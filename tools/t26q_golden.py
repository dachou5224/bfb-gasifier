"""T26q 基准快照：``hamel/`` 行数清理前后的逐位比较（工况 1）。

``python -m tools.t26q_golden capture --out .scratch/t26q/golden.npz``
``python -m tools.t26q_golden compare --golden .scratch/t26q/golden.npz``

全程 φ = 0,3（工况副本，不做 φ 切换）；另记 φ = 1 时初值构造被拒的完整 diagnostics。
外层只进第 1 轮、内层记前 3 次迭代，用包装函数在 3 次后中止，不改正式路径。
Corella 外推只在本工具里设置，与 T26o、T26p 相同。比较用 ``np.array_equal``。
"""
from __future__ import annotations

import argparse
import json
import re
import signal
import subprocess
import time
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

import hamel.kinetics.homogeneous as homogeneous
import hamel.solver.driver as driver
import hamel.solver.newton as newton
from hamel.cases.config import SolverConfig
from hamel.solver.initial import InitialStateError
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_PHI = 0.3
_SOLID = ("bed_char_kg", "freeboard_char_free_kg", "freeboard_char_hover_kg",
          "freeboard_ash_free_kg", "freeboard_ash_hover_kg")
# ADR 0052 遗留的成功交替报告项；B2 允许删除，快照不比较它们。
_DROP = {"step_history", "volume_ratio", "min_stock_kg"}
# ADR 0056 第 5 项写入牛顿记录；φ = 0,3 第 1 轮不投影，值为 0。比较时排除，不改记录内容。
_EXCLUDE = (".projected", ".overshoot")


class _Stop(BaseException):
    """内层记满 3 次迭代后中止，不改变正式路径。"""


class _Timeout(BaseException):
    """硬超时。CHARTER §7。"""


def _sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=_ROOT, text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _phi_copy() -> Path:
    folder = _ROOT / ".scratch" / "t26q"
    folder.mkdir(parents=True, exist_ok=True)
    doc = json.loads(_CASE.read_text(encoding="utf-8"))
    doc["cases"]["1"]["fly_ash_fraction"]["value"] = _PHI
    path = folder / "case_phi03.json"
    path.write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
    return path


def _flatten(prefix: str, value, out: dict) -> None:
    """把嵌套结果压成 npz 的扁平数组。键按字典序，比较不依赖插入顺序。"""
    if isinstance(value, dict):
        for key in sorted(value, key=str):
            _flatten(f"{prefix}.{key}", value[key], out)
    elif isinstance(value, (list, tuple)) and not all(
        np.isscalar(item) and not isinstance(item, (str, bytes, dict)) for item in value
    ):
        for index, item in enumerate(value):
            _flatten(f"{prefix}[{index}]", item, out)
    elif value is None:
        out[prefix] = np.array("<none>")
    elif isinstance(value, (list, tuple)):
        out[prefix] = np.asarray(value)
    elif callable(value):
        out[prefix] = np.array("<callable>")
    elif isinstance(value, (np.ndarray, np.generic)) or np.isscalar(value):
        out[prefix] = np.asarray(value)
    else:
        out[prefix] = np.array(str(value))


def _run(limit_s: float, phi03: bool) -> dict:
    """φ = 0,3 的完整快照，或 φ = 1 的初值构造拒绝终态。"""
    path = _phi_copy() if phi03 else _CASE
    box = {"pre": [], "construct": None, "diagnostics": None, "first": None}
    notes: list = []
    originals = (driver.precalculate, driver.construct, newton._note)

    def precalculate(*args, **kwargs):
        result = originals[0](*args, **kwargs)
        if not phi03:
            return result
        if box["first"] is None:
            box["first"] = result
        previous = args[2] if len(args) > 2 else kwargs.get("previous")
        entry = {
            "w_f": result.frozen.bed_ash_fine_share, "phi": result.frozen.fly_ash_fraction,
            "drain_ash_kg_s": result.frozen.drain_ash_kg_s, "drain_char_kg_s": result.frozen.drain_char_kg_s,
            "ash_balance": result.reported.get("ash_balance"),
        }
        if previous is not None:
            entry["state"] = {name: getattr(previous, name) for name in _SOLID}
        box["pre"].append(entry)
        return result

    def construct(*args, **kwargs):
        try:
            box["construct"] = originals[1](*args, **kwargs)
            return box["construct"]
        except InitialStateError as exc:
            box["diagnostics"], box["exception"] = exc.diagnostics, f"{type(exc).__name__}: {exc}"
            raise

    def note(journal, iteration, merit, residual, delta, trials, scales, jacobian_ok=True):
        originals[2](journal, iteration, merit, residual, delta, trials, scales, jacobian_ok)
        if not phi03 or journal is None or journal.get("blocks") is None:
            return
        entry = journal["iterations"][-1]
        notes.append({
            "iteration": iteration, "merit": merit,
            "delta": None if delta is None else np.asarray(delta, dtype=float),
            "blocks": entry.get("blocks"), "step_max": (entry.get("step") or {}).get("max"),
            "trials": entry.get("trials"), "jacobian_ok": jacobian_ok,
        })
        if len(notes) >= 3:
            raise _Stop()

    driver.precalculate, driver.construct, newton._note = precalculate, construct, note
    started, summary = time.monotonic(), {}
    try:
        result = driver.solve_case(str(path), "1", config=replace(SolverConfig(), outer_round_limit=1),
                                   deadline=lambda: time.monotonic() - started >= limit_s)
        summary.update(status=result.status, reason=result.reason,
                       rounds=result.rounds, iterations=result.iterations)
    except _Stop:
        summary.update(status="STOPPED_AFTER_THREE_INNER_ITERATIONS", reason="包装函数中止")
    except Exception as exc:  # noqa: BLE001 - _Timeout 继续上抛
        summary.update(status=type(exc).__name__, reason=str(exc)[:800])
    finally:
        driver.precalculate, driver.construct, newton._note = originals
    if not phi03:
        return {"diagnostics": box["diagnostics"], "exception": box.get("exception"),
                "raised": box.get("raised") or (None if box["diagnostics"] else summary["status"])}
    data = {"summary": summary, "pre": box["pre"], "notes": tuple(notes)}
    if box["construct"] is not None:
        data["construct"] = {
            "note": box["construct"][2],
            "records": [{k: v for k, v in row.items() if k not in _DROP} for row in box["construct"][1]],
        }
    if box["first"] is not None:
        data["first"] = {field.name: getattr(box["first"].frozen, field.name)
                         for field in fields(box["first"].frozen)}
        data["first_residual"] = evaluate_residual(
            pack_state(box["first"].state, box["first"].frozen), box["first"].frozen,
        ).values
    return data


def _expire(signum, frame):
    raise _Timeout()


def _snapshot(limit_s: float, rejection: bool = True) -> dict:
    signal.signal(signal.SIGALRM, _expire)
    signal.setitimer(signal.ITIMER_REAL, float(limit_s))
    original = homogeneous._CORELLA_HIGH_K
    homogeneous._CORELLA_HIGH_K = float("inf")
    try:
        phi = _run(limit_s, True)
        extra = _run(limit_s, False) if rejection else {}
    finally:
        homogeneous._CORELLA_HIGH_K = original
        signal.setitimer(signal.ITIMER_REAL, 0)
    out: dict = {}
    _flatten("data.phi03", phi, out)
    if extra:
        _flatten("data.rejection", extra, out)
    return {re.sub(r"[^0-9A-Za-z_.-]+", "_", key): np.asarray(value) for key, value in out.items()}


def capture(path: Path, limit_s: float, base_sha: str) -> int:
    try:
        data = _snapshot(limit_s)
    except _Timeout:
        print(json.dumps({"out": str(path), "status": "TIMEOUT", "limit_s": limit_s}, ensure_ascii=False))
        return 2
    data.update({"meta.base_sha": np.array(base_sha), "meta.phi": np.array(_PHI),
                 "meta.case": np.array("CASE_VTT_PRESSURIZED_1")})
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)
    print(json.dumps({"out": str(path), "base_sha": base_sha, "arrays": len(data)}, ensure_ascii=False))
    return 0


def compare(golden: Path, limit_s: float, prefix: str = "data.") -> int:
    reference = np.load(golden, allow_pickle=False)
    try:
        current = _snapshot(limit_s, rejection=not prefix.startswith("data.phi03."))
    except _Timeout:
        print(json.dumps({"golden": str(golden), "status": "TIMEOUT"}, ensure_ascii=False))
        return 2

    def kept(key: str) -> bool:
        return not any(key.endswith(item) for item in _EXCLUDE)

    keys = sorted(key for key in reference.files if key.startswith(prefix) and kept(key))
    missing = [key for key in keys if key not in current]
    extra = sorted(key for key in current if key.startswith(prefix) and key not in reference.files and kept(key))
    excluded = sorted(key for key in current if key.startswith(prefix) and not kept(key))
    different = [key for key in keys if key in current and not np.array_equal(reference[key], current[key])]
    print(json.dumps({
        "golden": str(golden), "prefix": prefix, "arrays": len(keys),
        "equal": not (missing or extra or different),
        "missing": missing[:10], "extra": extra[:10], "different": different[:10],
        "different_count": len(different), "excluded": excluded[:10], "excluded_count": len(excluded),
        "base_sha": str(reference["meta.base_sha"]),
    }, ensure_ascii=False))
    return 0 if not (missing or extra or different) else 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("capture", "compare"))
    parser.add_argument("--out", default=str(_ROOT / ".scratch" / "t26q" / "golden.npz"))
    parser.add_argument("--golden", default=None)
    parser.add_argument("--limit-s", type=float, default=6.0 * 3600.0)
    parser.add_argument("--base-sha", default=None)
    parser.add_argument("--prefix", default="data.", help="只比较以此为前缀的数组（T26r：data.phi03.）")
    args = parser.parse_args(argv)
    if args.action == "capture":
        return capture(Path(args.out), args.limit_s, args.base_sha or _sha())
    return compare(Path(args.golden or args.out), args.limit_s, args.prefix)


if __name__ == "__main__":
    raise SystemExit(main())
