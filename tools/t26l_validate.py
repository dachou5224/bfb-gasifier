"""T26l 验证：按 ADR 0052 重做初值构造，正常结束后跑基准外层第 1 轮的内层。只报告。

判据预先登记：初值构造正常结束，且基准外层第 1 轮内层以 CONVERGED 终止。不做整批重跑。
报告每次交替的牛顿迭代数与残差历史、床区焦炭与各格惰性物、排灰、带出灰、回落、Ṁ_ab
（重算值与冻结值）、递推最底格闭合差、自由空间各格体积比与悬停存量，并与研究底稿
reports/research-initial-solids/study.json 的 1073 K 诊断结果对照。
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
_OUT = _ROOT / "reports" / "t26l-validation"
_RESEARCH = _ROOT / "reports" / "research-initial-solids" / "study.json"


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
    payload = {
        "elapsed_s": elapsed, "status": status, "outer_round_limit": 1,
        "reason": None if outcome is None else outcome.reason,
        "rounds": None if outcome is None else outcome.rounds,
        "iterations": None if outcome is None else outcome.iterations,
        "initial": None if outcome is None else _plain(outcome.diagnostics.get("initial")),
        "rounds_detail": None if outcome is None else _plain(outcome.diagnostics.get("rounds")),
        "research": _research(),
    }
    payload["inner_status"] = _inner_status(payload)
    payload["criteria"] = _criteria(payload)
    (folder / "validation.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    (folder / "validation.md").write_text(_markdown(payload), encoding="utf-8")
    print(
        f"{payload['status']} 轮 {payload['rounds']} 内层 {payload['iterations']} "
        f"外层第 1 轮 {payload['inner_status']} 墙钟 {elapsed:.0f} s"
    )


def _inner_status(payload) -> str | None:
    rows = payload.get("rounds_detail") or []
    return rows[0].get("inner") if rows else None


def _criteria(payload) -> dict:
    initial = payload.get("initial") or {}
    ended = payload.get("initial") is not None and not str(initial.get("note", "")).startswith("停止")
    return {
        "initial_normal_end": bool(ended),
        "first_round_inner_converged": payload.get("inner_status") == "CONVERGED",
    }


def _research():
    """研究底稿 1073 K、重算 Ṁ_ab 组的诊断结果。"""
    if not _RESEARCH.exists():
        return None
    data = json.loads(_RESEARCH.read_text(encoding="utf-8"))
    for run in data.get("runs") or []:
        if run.get("temperature_k") == 1073.0 and run.get("recompute_down"):
            return run
    return None


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
    lines = ["# T26l 验证：初值构造的固相牛顿（ADR 0052）与基准外层第 1 轮内层", ""]
    lines += ["## 判据", "", "- 初值构造正常结束。", "- 基准外层第 1 轮内层以 `CONVERGED` 终止。", ""]
    criteria = payload.get("criteria") or {}
    lines += ["## 结论", ""]
    lines.append(f"- `{payload['status']}`，外层 {payload['rounds']} 轮，内层 {payload['iterations']} 次，"
                 f"墙钟 {payload['elapsed_s']:.0f} s")
    lines.append(f"- 初值构造正常结束：{'是' if criteria.get('initial_normal_end') else '否'}；"
                 f"外层第 1 轮内层：`{payload.get('inner_status')}`")
    if payload.get("reason"):
        lines.append(f"- 停止原因：{payload['reason']}")
    lines += ["", *_initial_lines(payload), *_research_lines(payload), *_inner_lines(payload)]
    return "\n".join(lines).rstrip("\n") + "\n"


def _initial_lines(payload) -> list:
    lines = ["## 初值构造", ""]
    initial = payload.get("initial") or {}
    if payload.get("initial") is None:
        return lines + ["（没有取到，运行未进入初值构造）", ""]
    lines.append(f"- 交替说明：{initial.get('note') or '—'}")
    for record in initial.get("loops") or []:
        if record.get("stopped"):
            lines += _stop_lines(record)
            continue
        history = record.get("residual_history") or []
        lines.append(
            f"- 交替 {record.get('loop')}：牛顿 {record.get('iterations')} 次，"
            f"固相门 {_fmt(record.get('solid_residual'))}，残差历史 "
            f"{' → '.join(f'{float(value):.3g}' for value in history) if history else '—'}"
        )
        if record.get("bed_char_kg") is not None:
            lines.append(f"  - 床区焦炭 {_fmt(record.get('bed_char_kg'))} kg；各床格惰性物 kg："
                         f"{_fmt(record.get('bed_inert_kg'))}")
            lines.append(f"  - 排灰：灰 {_fmt(record.get('drain_ash_kg_s'))} kg/s，焦炭 {_fmt(record.get('drain_char_kg_s'))} kg/s；"
                         f"带出灰 {_fmt(record.get('elutriated_ash_kg_s'))} kg/s；床面回落 {_fmt(record.get('return_kg_s'))} kg/s")
            lines.append(f"  - Ṁ_ab 冻结值 kg/s：{_fmt(record.get('frozen_down_kg_s'))}；重算值 kg/s：{_fmt(record.get('down_kg_s'))}；"
                         f"递推最底格闭合差 {_fmt(record.get('bottom_closure_kg_s'))} kg/s")
            lines.append(f"  - 悬停焦炭 {_fmt(record.get('hover_char_kg'))} kg；体积比（自下而上）：{_fmt(record.get('volume_ratio'), 3)}")
            lines.append(f"  - 界面相对变化：{_fmt_change(record.get('change'))}")
    return [*lines, ""]


def _stop_lines(record: dict) -> list:
    """牛顿停止的那次交替：起点、残差历史与阻尼试探。ADR 0052 第 3 项。"""
    lines = [f"- 交替 {record.get('loop')}：牛顿停止 —— {record.get('reason')}"]
    lines.append(f"  - 牛顿状态 `{record.get('status')}`，迭代 {record.get('iterations')} 次，"
                 f"ν_max = {record.get('nu_max')}；起点固相未知量中为 0 的有 {record.get('zero_stocks')} 个")
    if record.get("bed_char_kg") is not None:
        lines.append(f"  - 起点的预计算冻结量：床区焦炭 {_fmt(record.get('bed_char_kg'))} kg；各床格惰性物 kg："
                     f"{_fmt(record.get('bed_inert_kg'))}")
        lines.append(f"  - 起点的排灰：灰 {_fmt(record.get('drain_ash_kg_s'))} kg/s，焦炭 "
                     f"{_fmt(record.get('drain_char_kg_s'))} kg/s；Ṁ_ab 冻结值 kg/s：{_fmt(record.get('frozen_down_kg_s'))}")
    history = record.get("merit_history") or []
    if history:
        lines.append(f"  - 残差历史（归一欧氏范数）：{' → '.join(f'{float(value):.4g}' for value in history)}")
    iterations = ((record.get("journal") or {}).get("iterations")) or []
    if iterations:
        last = iterations[-1]
        trials = last.get("trials") or []
        reasons = {}
        for trial in trials:
            key = f"{trial.get('status')} {trial.get('reason') or ''}".strip()
            reasons[key] = reasons.get(key, 0) + 1
        spread = "；".join(f"{name} {count} 次" for name, count in reasons.items())
        lines.append(f"  - 第 {last.get('iteration')} 次迭代的阻尼：{len(trials)} 个试探点（i = 0…{max(len(trials) - 1, 0)}），{spread}")
        step = last.get("step") or {}
        largest = step.get("max") or {}
        if largest:
            lines.append(f"  - 牛顿步：归一范数 {_fmt(step.get('norm'))}，最大分量 {largest.get('name')} = "
                         f"{_fmt(largest.get('value'))}（× 尺度），原始 {_fmt(largest.get('raw'))} kg")
        if len(trials) == 31 and len(reasons) == 1:
            lines.append("  - 起点按 ADR 0045 第 3 项第 1 步把自由空间存量取 0；牛顿步在这些 0 存量上为负时，"
                         "任何 i ≤ 30 的缩步仍为负，ADR 0017 第 7 项的可行条件无法满足，按 ADR 0052 第 3 项记为阻尼失败。")
    return lines


def _research_lines(payload) -> list:
    lines = ["## 与研究底稿 1073 K 诊断对照", ""]
    old = payload.get("research")
    if not old:
        return lines + ["（本地没有研究底稿的结果）", ""]
    lines.append(f"- 研究底稿：温度 1073 K、重算 Ṁ_ab；收敛 {old.get('converged')}，"
                 f"残差历史 {' → '.join(f'{float(v):.3g}' for v in old.get('residual_history') or [])}")
    lines.append(f"- 研究底稿床区焦炭 {_fmt(old.get('bed_char_kg'))} kg；各床格惰性物 kg：{_fmt(old.get('bed_inert_by_cell_kg'))}")
    lines.append(f"- 研究底稿排灰：灰 {_fmt(old.get('drain_ash_kg_s'))} kg/s，焦炭 {_fmt(_sum(old.get('drain_char_kg_s')))} kg/s；"
                 f"带出灰 {_fmt(old.get('elutriated_ash_kg_s'))} kg/s；悬停焦炭 {_fmt(old.get('hover_char_kg'))} kg")
    flows = old.get("flows") or {}
    lines.append(f"- 研究底稿 Ṁ_ab 冻结值 kg/s：{_fmt(flows.get('frozen_down'))}；重算值 kg/s：{_fmt(flows.get('down'))}；"
                 f"闭合差 {_fmt(flows.get('closure'))} kg/s；最大体积比 {_fmt(_max(old.get('volume_ratio')), 3)}")
    latest = _latest(payload)
    if not latest:
        lines.append("- 本次没有完成一次交替（固相牛顿在第一次交替即停止），没有同一阶段的存量与排灰可比。")
    if latest:
        lines.append("- 本次同一阶段（最后一次交替）：")
        lines.append(f"  - 床区焦炭 {_fmt(latest.get('bed_char_kg'))} kg；各床格惰性物 kg：{_fmt(latest.get('bed_inert_kg'))}")
        lines.append(f"  - 排灰：灰 {_fmt(latest.get('drain_ash_kg_s'))} kg/s，焦炭 {_fmt(_sum(latest.get('drain_char_kg_s')))} kg/s；"
                     f"带出灰 {_fmt(latest.get('elutriated_ash_kg_s'))} kg/s；悬停焦炭 {_fmt(latest.get('hover_char_kg'))} kg")
        lines.append(f"  - Ṁ_ab 冻结值 kg/s：{_fmt(latest.get('frozen_down_kg_s'))}；重算值 kg/s：{_fmt(latest.get('down_kg_s'))}；"
                     f"闭合差 {_fmt(latest.get('bottom_closure_kg_s'))} kg/s；最大体积比 {_fmt(_max(latest.get('volume_ratio')), 3)}")
    return [*lines, ""]


def _latest(payload):
    loops = [
        item for item in (payload.get("initial") or {}).get("loops") or []
        if item.get("bed_char_kg") is not None and not item.get("stopped")
    ]
    return loops[-1] if loops else None


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
    return [*lines, ""]


def _sum(values):
    if not values:
        return None
    return float(np.sum(np.asarray(values, dtype=float)))


def _max(values):
    if not values:
        return None
    return float(np.max(np.asarray(values, dtype=float)))


def _fmt(values, digits: int = 4) -> str:
    """标量（固相门、床区焦炭）与序列按同一写法输出。"""
    if values is None:
        return "—"
    if np.isscalar(values) or not hasattr(values, "__len__"):
        return f"{float(values):.{digits}g}"
    return ", ".join(f"{float(value):.{digits}g}" for value in values)


def _fmt_change(change) -> str:
    if not change:
        return "—"
    return ", ".join(f"{name} {float(value):.2e}" for name, value in change.items())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
