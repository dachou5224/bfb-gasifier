"""T26o 验证：ADR 0054 预先登记的两条判据（工况 1，外层 1 轮）。

只报告，不改 ``hamel/``。Corella 外推只在本工具里设置（``_CORELLA_HIGH_K = +inf``），
不进入正式路径（ADR 0010 第 9 条未批准外推）。

子命令
- frozen：不外推 Corella。判据：初值构造 4 次交替结束，A5 排灰灰分量 +0,153 g/s（±1 %）；
  第 1 轮内层只能停在最底床格 900 °C 处的 Corella 拒绝差分点（O4），不出现 O1、O2、O3、O5。
- extrapolated：临时外推 Corella。判据：第 1 轮内层 CONVERGED，迭代 ≤ 100，三块归一残差 ≤ 1·10⁻⁶，
  各格温度与 ``reports/research-solver-path/diagnostic/outer.json`` 的 rounds[0] 相差 < 1 K。

两次运行都记录每次迭代 y 上最大的步长分量与雅可比求解是否成功（ADR 0054 验证第 4 条）。
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.driver as driver
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.solver.gates import block_ratios
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import unpack_state

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "t26o-validation"
_BASELINE = _ROOT / "reports" / "research-solver-path" / "diagnostic" / "outer.json"
_TOLERANCE = 1.0e-6
_MAX_STEP = 1.0e6
# Corella 有效域拒绝的原文（hamel/kinetics/homogeneous.py，ADR 0010 第 9 条）。
_CORELLA_REJECTION = "Corella 参数只覆盖"
_BOTTOM_BOUNDARY_K = 900.0
_BOTTOM_TOLERANCE_K = 0.01


def run(extrapolate: bool, limit_s: float) -> dict:
    """跑一次工况 1 的外层 1 轮，收集初值构造与第 1 轮内层的记录。"""
    if extrapolate:
        import hamel.kinetics.homogeneous as homogeneous
        homogeneous._CORELLA_HIGH_K = float("inf")
    box: dict = {}
    original = driver.precalculate

    def precalculate(*args, **kwargs):
        result = original(*args, **kwargs)
        box["frozen"] = result.frozen
        return result

    driver.precalculate = precalculate
    started = time.monotonic()
    try:
        result = driver.solve_case(
            str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=1),
            deadline=lambda: time.monotonic() - started >= limit_s,
        )
        summary = {"status": result.status, "reason": result.reason, "rounds": result.rounds,
                   "iterations": result.iterations}
    except Exception as exc:  # noqa: BLE001 - 运行记录要留下任何异常
        result, summary = None, {"status": type(exc).__name__, "reason": str(exc)[:800]}
    finally:
        driver.precalculate = original
    summary["elapsed_s"] = time.monotonic() - started
    if result is None:
        return {"summary": summary, "extrapolated": extrapolate}
    summary["initial"] = _initial(result.diagnostics.get("initial") or {})
    summary["inner"] = _inner(result.diagnostics.get("rounds") or [])
    summary["final_check"] = result.diagnostics.get("final_check")
    if result.x is not None and "frozen" in box:
        values = np.asarray(result.x, dtype=float)
        frozen = box["frozen"]
        state = unpack_state(values, frozen)
        summary["temperature_c"] = {
            "bed": (np.asarray(state.bed_temperature_k) - 273.15).tolist(),
            "freeboard": (np.asarray(state.freeboard_temperature_k) - 273.15).tolist(),
        }
        # 收敛点上的残差（journal 里是每次迭代开始时的残差，不是最后一跳之后的）。
        blocks = driver._blocks(frozen, *driver._references(_physical()))
        summary["final_ratios"] = block_ratios(evaluate_residual(values, frozen).values, blocks)
    return {"summary": summary, "extrapolated": extrapolate, "x": None if result.x is None else np.asarray(result.x)}


def _physical():
    return assemble_case(str(_CASE), "1").physical


def _initial(initial: dict) -> dict:
    loops = initial.get("loops") or []
    last = loops[-1] if loops else {}
    # 停止的那次交替：pre_solve 是牛顿前的状态，rejected 是被第 3 层拒的收敛解。
    before = last.get("pre_solve") or {}
    return {
        "note": initial.get("note"), "loops": len(loops),
        "drain_ash_kg_s": before.get("drain_ash_kg_s", last.get("drain_ash_kg_s")),
        "iterations": [item.get("iterations") for item in loops],
        "bed_char_kg": before.get("bed_char_kg", last.get("bed_char_kg")),
        "solid_residual": last.get("solid_residual"),
        "rejected": last.get("rejected"),
    }


def _inner(rounds: list) -> dict:
    """第 1 轮内层：停止原因、每次迭代的最大步长与雅可比是否可解。"""
    if not rounds:
        return {"rows": 0}
    first = rounds[0]
    journal = first.get("newton") or []
    iterations = []
    for entry in journal:
        step = entry.get("step") or {}
        largest = step.get("max") or {}
        iterations.append({
            "iteration": entry.get("iteration"), "merit": entry.get("merit"),
            "blocks": entry.get("blocks"), "jacobian_ok": entry.get("jacobian_ok"),
            "delta_y": largest.get("raw"), "delta_y_name": largest.get("name"),
            "delta_y_scaled": largest.get("value"),
            "accepted_i": next((item.get("i") for item in entry.get("trials") or [] if item.get("accepted")), None),
            "rejections": _rejections(entry.get("trials") or []),
        })
    largest = max((abs(item["delta_y"]) for item in iterations if item["delta_y"] is not None), default=0.0)
    return {
        "status": first.get("inner"), "reason": first.get("reason"), "rows": len(rounds),
        "iterations": len(iterations), "final_blocks": iterations[-1]["blocks"] if iterations else None,
        "final_check": None,  # 逐次迭代的结果在 steps 里；收敛解的那一次见 diagnostics。
        "max_delta_y": largest, "all_jacobians_ok": all(item["jacobian_ok"] for item in iterations),
        "steps": iterations,
    }


def _rejections(trials: list) -> dict:
    counts: dict = {}
    for item in trials:
        if item.get("status") != "ok":
            key = str(item.get("reason", ""))[:60]
            counts[key] = counts.get(key, 0) + 1
    return counts


def compare(extrapolated: dict, frozen: dict) -> dict:
    """按 ADR 0054 预先登记的判据逐条给结论。"""
    baseline = json.loads(_BASELINE.read_text(encoding="utf-8"))["rounds"][0]
    initial = frozen["summary"].get("initial") or {}
    inner = frozen["summary"].get("inner") or {}
    cool = extrapolated["summary"].get("inner") or {}
    out: dict = {}
    drain = initial.get("drain_ash_kg_s")
    target = 0.153e-3
    out["initial_loops"] = {"value": initial.get("loops"), "target": 4,
                            "ok": initial.get("loops") == 4}
    out["initial_drain_ash_g_s"] = {
        "value": None if drain is None else drain * 1.0e3, "target": 0.153,
        "ok": drain is not None and abs(drain - target) <= 0.01 * target,
    }
    bed = (frozen["summary"].get("temperature_c") or {}).get("bed") or []
    out["frozen_stop"] = {"status": inner.get("status"), "reason": inner.get("reason"),
                          "bottom_bed_c": None if not bed else bed[-1]}
    for name in ("O1", "O2", "O3", "O5"):
        out[f"frozen_no_{name}"] = {"ok": not _has(inner, name)}
    out["frozen_o4"] = _is_o4(inner, bed)
    out["extrapolated_converged"] = {"status": cool.get("status"), "ok": cool.get("status") == "CONVERGED"}
    blocks = extrapolated["summary"].get("final_ratios") or {}
    out["extrapolated_residuals"] = {
        "value": {name: blocks.get(name) for name in ("gas", "solid", "energy")},
        "ok": bool(blocks) and max(float(blocks[name]) for name in ("gas", "solid", "energy")) <= _TOLERANCE,
    }
    out["extrapolated_iterations"] = {"value": cool.get("iterations"), "limit": 100,
                                      "ok": (cool.get("iterations") or 0) <= 100}
    temperature = extrapolated["summary"].get("temperature_c") or {}
    out["temperature_match"] = _temperature(temperature, baseline)
    # ADR 0054 验证第 4 条：Δy ≥ 1·10⁶ 即停止上报。冻结模型里“雅可比不可解”是判据 1 预期的 O4 停止，
    # 由 frozen_o4 单独判定；外推组内层收敛，逐次雅可比都必须可解。
    for side in ("frozen", "extrapolated"):
        inner_side = (frozen if side == "frozen" else extrapolated)["summary"].get("inner") or {}
        ok = (inner_side.get("max_delta_y") or 0.0) < _MAX_STEP
        if side == "extrapolated":
            ok = ok and inner_side.get("all_jacobians_ok") is True
        out[f"{side}_steps"] = {
            "max_delta_y": inner_side.get("max_delta_y"), "jacobians_ok": inner_side.get("all_jacobians_ok"),
            "ok": ok,
        }
    out["all_ok"] = all(item.get("ok") for key, item in out.items() if isinstance(item, dict) and "ok" in item)
    return out


def _temperature(temperature: dict, baseline: dict) -> dict:
    differences = {}
    for name, key in (("bed", "bed_temperature_c"), ("freeboard", "freeboard_temperature_c")):
        ours = temperature.get(name)
        theirs = baseline.get(key)
        if not ours or not theirs or len(ours) != len(theirs):
            return {"ok": False, "differences": differences, "reason": "温度形状不一致"}
        differences[name] = [float(a - b) for a, b in zip(ours, theirs)]
    largest = max(abs(value) for values in differences.values() for value in values)
    return {"ok": largest < 1.0, "largest_k": largest, "differences": differences}


def _has(inner: dict, code: str) -> bool:
    """按诊断底稿的障碍定义判定停止原因是否属于 O1、O2、O3、O5。"""
    reason = str(inner.get("reason") or "")
    rejections = " ".join(
        str(item.get("rejections") or {}) for item in (inner.get("steps") or [])
    )
    if code == "O1":
        return inner.get("iterations") == 0
    if code == "O2":
        return "O2" in reason or "O₂" in reason
    if code == "O3":
        return "DAMPING_FAILED" == inner.get("status") and "m < 0" in rejections
    if code == "O5":
        return "DAMPING_FAILED" == inner.get("status") and "H2" in rejections
    return False


def _is_o4(inner: dict, bed_celsius: list) -> dict:
    """判据 1 的停止点：Corella 有效域在差分点上被拒，且最底床格在 900 °C（PR #152 审查 P2）。

    两条都要满足：(a) 停止原因是 Corella 范围拒绝的原文；(b) 最底床格温度与 900 °C
    相差 < 0.01 K。只凭“雅可比不能求解”这个前缀不算 O4 —— 任何差分点不可评价都带它。
    """
    reason = str(inner.get("reason") or "")
    bottom = None if not bed_celsius else float(bed_celsius[-1])
    corella = _CORELLA_REJECTION in reason
    near = bottom is not None and abs(bottom - _BOTTOM_BOUNDARY_K) < _BOTTOM_TOLERANCE_K
    return {
        "ok": bool(corella and near), "reason": reason,
        "corella_rejection": corella, "difference_point": "差分点不可评价" in reason,
        "bottom_bed_c": bottom, "bottom_within_tolerance_k": _BOTTOM_TOLERANCE_K,
        "bottom_within_tolerance": near,
    }


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit-s", type=float, default=7200.0)
    parser.add_argument("--out", default=str(_OUT))
    parser.add_argument("--only", choices=("frozen", "extrapolated"))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    frozen, extrapolated = None, None
    if args.only in (None, "frozen"):
        frozen = run(False, args.limit_s)
        _write(out / "frozen.json", frozen)
        print("frozen", frozen["summary"]["status"], frozen["summary"].get("inner", {}).get("status"))
    if args.only in (None, "extrapolated"):
        extrapolated = run(True, args.limit_s)
        _write(out / "extrapolated.json", extrapolated)
        print("extrapolated", extrapolated["summary"]["status"], extrapolated["summary"].get("inner", {}).get("status"))
    if frozen is not None and extrapolated is not None:
        verdict = compare(extrapolated, frozen)
        _write(out / "verdict.json", verdict)
        print(json.dumps({key: value for key, value in verdict.items() if isinstance(value, dict)},
                         ensure_ascii=False, indent=1))


def _write(path: Path, payload: dict) -> None:
    values = payload.pop("x", None)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=_plain), encoding="utf-8")
    if values is not None:
        np.save(path.with_suffix(".npy"), np.asarray(values, dtype=float))


def _plain(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return str(value)


if __name__ == "__main__":
    main()
