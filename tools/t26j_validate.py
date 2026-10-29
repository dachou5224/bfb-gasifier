"""T26j 验证：重做初值构造，再从该初值跑基准外层第 1 轮的内层牛顿。

按 T26g 任务卡 C 部分的判据（ADR 0045、ADR 0049）：报告每一步的存量、体积比与停止原因；
初值构造正常结束时才进入基准第 1 轮内层。只报告，不整批重跑。
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
_OUT = _ROOT / "reports" / "t26j-validation"


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
        status = result.status
    except SolveTimeout:
        result, status = None, "TIMEOUT"
    elapsed = time.monotonic() - started
    payload = {
        "elapsed_s": elapsed,
        "status": status,
        "reason": None if result is None else result.reason,
        "rounds": None if result is None else result.rounds,
        "iterations": None if result is None else result.iterations,
        "ratios": None if result is None else result.ratios,
        "initial": None if result is None else result.diagnostics.get("initial"),
        "rounds_detail": None if result is None else _plain(result.diagnostics.get("rounds")),
        "geldart": None if result is None else _plain(result.diagnostics.get("geldart")),
    }
    (folder / "validation.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    (folder / "validation.md").write_text(_markdown(payload), encoding="utf-8")
    print(f"{payload['status']} 轮 {payload['rounds']} 内层 {payload['iterations']} 墙钟 {elapsed:.0f} s")


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
    lines = ["# T26j 验证：初值构造与基准外层第 1 轮内层", ""]
    lines.append(f"- 结论：`{payload['status']}`，外层 {payload['rounds']} 轮，内层 {payload['iterations']} 次，墙钟 {payload['elapsed_s']:.0f} s")
    if payload["reason"]:
        lines.append(f"- 停止原因：{payload['reason']}")
    lines.append("")
    initial = payload.get("initial") or {}
    lines.append("## 初值构造")
    lines.append("")
    lines.append(f"- 交替说明：{initial.get('note') or '—'}")
    for record in initial.get("loops") or []:
        if "stage" in record:
            lines.append(f"- 预启动：自由空间固相门 {record['freeboard_residual']:.3e}，床区固相行（未求解）{record['bed_residual']:.3e}")
            lines.append(f"  - 各格自由焦炭 kg：{_fmt(record['char_free_kg'])}")
            lines.append(f"  - 各格悬停焦炭 kg：{_fmt(record['char_hover_kg'])}")
            lines.append(f"  - 体积比（自上而下）：{_fmt(record['volume_ratio'])}")
            lines.append(f"  - 界面流前：{_fmt(record['interface_before']['up'])} → 后：{_fmt(record['interface_after']['up'])}（上行 kg/s）")
            lines.append(f"  - 回落前：{_fmt(record['interface_before']['return'])} → 后：{_fmt(record['interface_after']['return'])}")
            continue
        lines.append(
            f"- 交替 {record['loop']}：固相门 {_fmt_scalar(record['solid_residual'])}，"
            f"最小 m_C {_fmt_scalar(record['minimum_char_kg'])} kg，床区焦炭 {_fmt_scalar(record['bed_char_kg'])} kg，"
            f"相对变化 {_fmt_change(record.get('change'))}"
        )
        if record.get("note"):
            lines.append(f"  - 备注：{record['note']}")
    lines.append("")
    return "\n".join(lines) + "\n"


def _fmt(values) -> str:
    if values is None:
        return "—"
    return ", ".join(f"{float(value):.4g}" for value in values)


def _fmt_scalar(value) -> str:
    return "—" if value is None else f"{float(value):.4g}"


def _fmt_change(change) -> str:
    if not change:
        return "—"
    return ", ".join(f"{name} {float(value):.2e}" for name, value in change.items())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
