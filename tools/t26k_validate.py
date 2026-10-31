"""T26k 验证：按 ADR 0051 实现粒级内漂移后重做初值构造，正常结束时跑基准外层第 1 轮的内层。

只报告，不整批重跑。报告第 1 次预计算的 δ、μ 与逐 φ_p 的去向，初值构造每一步的存量、
体积比、回代残差与界面流，与 T26j 预启动点（闭合 A）同一阶段的悬停存量与体积比对照。
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.solver.driver import solve_case
from hamel.solver.newton import SolveTimeout

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "t26k-validation"
_OLD = _ROOT / "reports" / "t26j-validation" / "validation.json"


def main(out: str | None = None) -> None:
    folder = Path(out) if out else _OUT
    folder.mkdir(parents=True, exist_ok=True)
    config = replace(SolverConfig(), outer_round_limit=1)
    started = time.monotonic()
    deadline = started + float(config.wall_timeout_s)
    try:
        result = solve_case(
            str(_CASE), "1", deadline=lambda: time.monotonic() >= deadline, config=config,
            artifact_dir=str(folder),
        )
        status, outcome = result.status, result
    except SolveTimeout:
        status, outcome = "TIMEOUT", None
    elapsed = time.monotonic() - started
    payload = {"elapsed_s": elapsed, "status": status, "t26j": _baseline()}
    payload.update({"reason": None, "rounds": None, "iterations": None, "initial": None,
                    "rounds_detail": None, "hover_drift": None} if outcome is None else {
        "reason": outcome.reason, "rounds": outcome.rounds, "iterations": outcome.iterations,
        "initial": _plain(outcome.diagnostics.get("initial")),
        "rounds_detail": _plain(outcome.diagnostics.get("rounds")),
        "hover_drift": _plain(outcome.diagnostics.get("hover_drift")),
    })
    (folder / "validation.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    (folder / "validation.md").write_text(_markdown(payload), encoding="utf-8")
    print(f"{payload['status']} 轮 {payload['rounds']} 内层 {payload['iterations']} 墙钟 {elapsed:.0f} s")


def _baseline():
    """T26j（闭合 A）预启动点：悬停存量与体积比，用于同一阶段对照。"""
    if not _OLD.exists():
        return None
    data = json.loads(_OLD.read_text(encoding="utf-8"))
    pick = [item for item in (data.get("initial") or {}).get("loops") or []
            if str(item.get("stage", "")).startswith("预启动自由空间焦炭块")]
    if not pick:
        return None
    return {"status": data.get("status"), "reason": data.get("reason"),
            "char_hover_kg": pick[0].get("char_hover_kg"), "volume_ratio": pick[0].get("volume_ratio")}


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def _markdown(payload: dict) -> str:
    lines = ["# T26k 验证：粒级内漂移（ADR 0051）、初值构造与基准外层第 1 轮内层", ""]
    lines.append(f"- 结论：`{payload['status']}`，外层 {payload['rounds']} 轮，内层 {payload['iterations']} 次，"
                 f"墙钟 {payload['elapsed_s']:.0f} s")
    if payload.get("reason"):
        lines.append(f"- 停止原因：{payload['reason']}")
    lines += ["", *_drift_lines(payload.get("hover_drift")), *_initial_lines(payload),
              *_baseline_lines(payload), *_inner_lines(payload)]
    return "\n".join(lines).rstrip("\n") + "\n"


def _drift_lines(drift) -> list:
    """第 1 次预计算的 δ^h、δ^f、μ 与逐 φ_p 的去向。ADR 0051 第 2 项。"""
    lines = ["## 第 1 次预计算：悬停窗口与漂移因子", ""]
    if not drift:
        return lines + ["（没有取到，运行未产生预计算记录）", ""]
    for item in drift.get("classes", []):
        span = f"{item['lower_m'] * 1e3:.4g}–{item['upper_m'] * 1e3:.4g} mm，ρ_k = {item['rho_log']:.4g}"
        active = [cell for cell in item["cells"]
                  if cell["drift_hover"] or cell["drift_free"] or cell["migration_mu"]]
        if not active:
            lines.append(f"- 粒级 {item['class']}（{span}）：没有悬停窗口")
            continue
        lines.append(f"- 粒级 {item['class']}（{span}）：")
        for cell in active:
            verdicts = {}
            for window in cell["windows"]:
                verdicts[window["fate"]] = verdicts.get(window["fate"], 0) + 1
            spread = ", ".join(f"{name} {count}" for name, count in sorted(verdicts.items()))
            lines.append(f"  - 几何格 {cell['cell']}：δ^h = {cell['drift_hover']:.4g}，"
                         f"δ^f = {cell['drift_free']:.4g}，μ = {cell['migration_mu']:.4g}；φ_p 去向 {spread}")
            detail = [
                f"φ={window['phi']:.2f}: {window['fate']} {window['low_m'] * 1e3:.4g}–{window['high_m'] * 1e3:.4g} mm"
                f"（L = {window['length']:.4g}，ω = {window['weight']:.4g}，出口因子 {window['exit']:.4g}）"
                for window in cell["windows"] if window["weight"] > 0.0
            ]
            if detail:
                lines.append(f"    - {'；'.join(detail)}")
    return [*lines, ""]


def _initial_lines(payload) -> list:
    lines = ["## 初值构造", ""]
    initial = payload.get("initial") or {}
    lines.append(f"- 交替说明：{initial.get('note') or '—'}")
    for record in initial.get("loops") or []:
        if record.get("stage"):
            lines.append(f"- {record['stage']}：体积比（自下而上）{_fmt(record.get('volume_ratio'), 3)}")
            lines.append(f"  - 各格悬停焦炭 kg：{_fmt_table(record.get('char_hover_kg'))}")
            lines.append(f"  - 各格自由焦炭 kg：{_fmt_table(record.get('char_free_kg'))}")
            if record.get("freeboard_residual") is not None:
                lines.append(f"  - 自由空间固相门 {record['freeboard_residual']:.3e}，"
                             f"床区固相行（未求解，只记录）{record['bed_residual']:.3e}")
            continue
        lines.append(f"- 交替 {record.get('loop')}：固相门 {_fmt(record.get('solid_residual'))}，"
                     f"最小 m_C {_fmt(record.get('minimum_char_kg'))} kg，"
                     f"床区焦炭 {_fmt(record.get('bed_char_kg'))} kg，相对变化 {_fmt_change(record.get('change'))}")
        for name, key in (("上行界面流", "up_kg_s"), ("下行界面流", "down_kg_s"), ("床面回落", "return_kg_s")):
            if record.get(key):
                lines.append(f"  - {name} kg/s：{_fmt(record[key], 4)}")
        if record.get("note"):
            lines.append(f"  - 备注：{record['note']}")
    return [*lines, ""]


def _baseline_lines(payload) -> list:
    old = payload.get("t26j")
    lines = ["## 与 T26j 预启动点（闭合 A）对照", ""]
    if not old:
        return lines + ["（本地没有 T26j 的验证结果）", ""]
    loops = [item for item in (payload.get("initial") or {}).get("loops") or []
             if str(item.get("stage", "")).startswith("预启动自由空间焦炭块")]
    lines.append(f"- T26j：`{old['status']}`；停止原因 {old['reason']}")
    lines.append(f"- T26j 悬停焦炭 kg（自上而下）：{_fmt_table(old.get('char_hover_kg'))}")
    lines.append(f"- T26j 体积比（自下而上）：{_fmt(old.get('volume_ratio'), 3)}")
    if not loops:
        return lines + ["- 本次没有走到预启动焦炭块，无法同阶段对照", ""]
    lines.append(f"- 本次同一阶段悬停焦炭 kg（自上而下）：{_fmt_table(loops[0].get('char_hover_kg'))}")
    lines.append(f"- 本次同一阶段体积比（自下而上）：{_fmt(loops[0].get('volume_ratio'), 3)}")
    total_new = float(np.sum(np.asarray(loops[0].get("char_hover_kg"), dtype=float)))
    total_old = float(np.sum(np.asarray(old.get("char_hover_kg"), dtype=float)))
    lines.append(f"- 悬停焦炭合计：本次 {total_new:.4g} kg，T26j {total_old:.4g} kg（比 {total_new / total_old:.3g}）")
    ratio_new = max(float(item) for item in loops[0].get("volume_ratio") or [0.0])
    ratio_old = max(float(item) for item in old.get("volume_ratio") or [0.0])
    lines.append(f"- 最大体积比：本次 {ratio_new:.4g}，T26j {ratio_old:.4g}（(1−ε_mf)·V 的倍数，上限 1）")
    return [*lines, ""]


def _inner_lines(payload) -> list:
    lines = ["## 基准外层第 1 轮的内层牛顿", ""]
    detail = payload.get("rounds_detail") or []
    if not detail:
        return lines + ["（未运行：初值构造未正常结束）", ""]
    for row in detail:
        lines.append(f"- 外层第 {row['round']} 轮内层：`{row['inner']}`，{row.get('reason') or '—'}")
        for item in row.get("newton") or []:
            blocks = item.get("blocks") or {}
            lines.append(f"  - 迭代 {item['iteration']}：归一残差 {item['merit']:.6g}，"
                         f"气相 {blocks.get('gas', 0.0):.4g}／固相 {blocks.get('solid', 0.0):.4g}／"
                         f"能量 {blocks.get('energy', 0.0):.4g}")
            for trial in item.get("trials") or []:
                lines.append(f"    - 试探 {trial['i']}（1/{2 ** trial['i']:g} 步）：{trial['status']}"
                             f" {trial.get('reason') or ''}")
    return [*lines, ""]


def _fmt(values, digits: int = 4) -> str:
    """标量（固相门、最小 m_C、床区焦炭）与序列都按同一写法输出。"""
    if values is None:
        return "—"
    if np.isscalar(values) or not hasattr(values, "__len__"):
        return f"{float(values):.{digits}g}"
    return ", ".join(f"{float(value):.{digits}g}" for value in values)


def _fmt_table(rows) -> str:
    if not rows:
        return "—"
    return "；".join("[" + ", ".join(f"{float(value):.4g}" for value in row) + "]" for row in rows)


def _fmt_change(change) -> str:
    if not change:
        return "—"
    return ", ".join(f"{name} {float(value):.2e}" for name, value in change.items())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
