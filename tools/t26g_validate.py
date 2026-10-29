"""T26g 验证：初值构造的每次循环，以及从该初值运行基准外层第 1 轮的内层牛顿。只报告。

求解手段：ADR 0045；解释性选择：ADR 0044。判据见任务卡第 10、11 项。
"""
from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.properties.data import GAS_SPECIES
from hamel.solver.campaign import _process, config_for
from hamel.solver.driver import SolveResult, _references, solve_case
from hamel.solver.initial import DomainRejected, InitialStateError, NegativeStock, construct
from hamel.solver.newton import SolveTimeout
from hamel.solver.precalc import precalculate
from hamel.solver.sensitivity import section_volume_m3
from tools.driver_report import convert

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "t26g-validation"
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_TOLERANCE = 1.0e-2


def _relative(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(_ROOT))
    except ValueError:
        return str(path)


def _num(value):
    return "—" if value is None else f"{float(value):.6e}"


def _cell_volume(frozen, cell: int) -> float:
    height = float(frozen.freeboard_height_m)
    z0 = (frozen.n_freeboard - 1 - cell) * height
    if math.isclose(float(frozen.bed_diameter_m), float(frozen.freeboard_diameter_m)):
        return math.pi * (0.5 * float(frozen.freeboard_diameter_m)) ** 2 * height
    return section_volume_m3(z0, z0 + height, frozen.bed_diameter_m, frozen.freeboard_diameter_m)


def _inventory(state, frozen) -> dict:
    """最终固相存量、最小 m_C 与各格体积比。ADR 0045 第 3、4 项。"""
    arrays = {
        "bed": np.asarray(state.bed_char_kg, dtype=float),
        "free": np.asarray(state.freeboard_char_free_kg, dtype=float),
        "hover": np.asarray(state.freeboard_char_hover_kg, dtype=float),
    }
    position = min(
        ((name, array, np.unravel_index(int(np.argmin(array)), array.shape)) for name, array in arrays.items() if array.size),
        key=lambda item: float(item[1][item[2]]),
    )
    cells = []
    for cell in range(frozen.n_freeboard):
        solid = float(
            state.freeboard_char_free_kg[cell].sum() + state.freeboard_char_hover_kg[cell].sum()
            + state.freeboard_ash_free_kg[cell].sum() + state.freeboard_ash_hover_kg[cell].sum()
        )
        volume = solid / float(frozen.rho_p_kg_m3)
        limit = (1.0 - float(frozen.epsilon_mf)) * _cell_volume(frozen, cell)
        cells.append({
            "cell": cell, "geometric": frozen.n_freeboard - 1 - cell, "solid_kg": solid, "volume_m3": volume,
            "limit_m3": limit, "ratio": None if limit <= 0.0 else volume / limit,
        })
    return {
        "bed_char_kg": float(state.bed_char_kg.sum()),
        "free_char_kg": float(state.freeboard_char_free_kg.sum()),
        "hover_char_kg": float(state.freeboard_char_hover_kg.sum()),
        "free_ash_kg": float(state.freeboard_ash_free_kg.sum()),
        "hover_ash_kg": float(state.freeboard_ash_hover_kg.sum()),
        "minimum_char": {"state": position[0], "m_C_kg": float(position[1][position[2]]),
                         "position": [int(item) for item in position[2]]},
        "cells": cells,
        "over_limit": [row["cell"] for row in cells if row["ratio"] is not None and row["ratio"] > 1.0],
    }


def _tar(state) -> list:
    index = {name: position for position, name in enumerate(GAS_SPECIES)}
    return [
        {"species": name, "bed_mol_s": float(state.bed_gas_mol_s[:, :, index[name]].sum()),
         "freeboard_mol_s": float(state.freeboard_gas_mol_s[:, index[name]].sum())}
        for name in ("C6H6", "C10H8")
    ]


def _precalc_closure():
    """与驱动的第一次外层同一条预计算链：排灰与燃料回落逐次回馈。"""
    holder = {"drain": 0.0, "fuel": 0.0}

    def precalc(state):
        result = precalculate(
            _CASE, "1", state, previous_drain_kg_s=holder["drain"],
            previous_fuel_return_kg_s=holder["fuel"], adiabatic=False, iteration=1,
        )
        holder["drain"] = float(result.frozen.drain_ash_kg_s + result.frozen.drain_char_kg_s.sum())
        holder["fuel"] = float(result.frozen.fuel_down_free_kg_s[-1].sum())
        return result

    return precalc


def initial() -> None:
    _n, m_ref, _q = _references(assemble_case(_CASE, "1").physical)
    started = time.perf_counter()
    first = precalculate(_CASE, "1", None, adiabatic=False, iteration=1)
    body = {"tar": _tar(first.state), "t_init_k": float(first.t_init_k), "t_ad_k": float(first.t_ad_k)}
    loops: list = []
    try:
        result, loops, note = construct(
            first, _precalc_closure(), m_ref, gate=float(SolverConfig().residual_tolerance),
            tolerance=_TOLERANCE, records=loops,
        )
    except (DomainRejected, NegativeStock, InitialStateError) as exc:
        body.update(status="STOPPED", reason=str(exc), loops=loops, note="初值构造未正常结束")
        if isinstance(exc, (DomainRejected, NegativeStock)):
            body.update(phase=exc.phase, rejected=_inventory(exc.state, first.frozen))
        if isinstance(exc, NegativeStock):
            body.update(condition=exc.condition, char_residual_before=exc.residual,
                        negative=sorted(exc.items, key=lambda item: item["solved_kg"]))
    else:
        body.update(status="OK", loops=loops, note=note, final=_inventory(result.state, result.frozen))
    body["seconds"] = time.perf_counter() - started
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "initial.json").write_text(json.dumps(convert(body), ensure_ascii=False, indent=2), encoding="utf-8")
    _write()
    print(json.dumps({"status": body["status"], "loops": len(loops), "seconds": body["seconds"]}, ensure_ascii=False))


def _execute_round() -> dict:
    config = replace(config_for("baseline"), outer_round_limit=1)
    folder = _OUT / "round1"
    folder.mkdir(parents=True, exist_ok=True)
    started, end = time.perf_counter(), time.monotonic() + config.wall_timeout_s
    try:
        result = solve_case(_CASE, adiabatic=False, deadline=lambda: time.monotonic() >= end, config=config, artifact_dir=folder)
    except SolveTimeout:
        result = SolveResult("TIMEOUT", "运行超过墙钟上限")
    record = convert({
        "name": "baseline", "change": "不改冻结项", "status": result.status, "reason": result.reason,
        "rounds": result.rounds, "iterations": result.iterations, "diagnostics": result.diagnostics,
        "wall_clock_s": time.perf_counter() - started, "outer_round_limit": 1,
    })
    rounds = (record.get("diagnostics") or {}).get("rounds") or []
    record["inner_status"] = rounds[0].get("inner") if rounds else None
    return record


def _spawn(command: str, timeout: float, folder: Path) -> None:
    """受监督的子进程，带硬超时。到点杀整个进程组，记 TIMEOUT。"""
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = _process([sys.executable, "-m", "tools.t26g_validate", command], timeout)
    if completed is None:
        record = {"status": "TIMEOUT", "reason": "墙钟上限内没有返回", "inner_status": None}
    elif completed.returncode != 0:
        log = folder / "stderr.log"
        log.write_text(completed.stderr or "", encoding="utf-8")
        record = {"status": "ERROR", "reason": f"子进程退出码 {completed.returncode}", "inner_status": None,
                  "stderr_log": _relative(log)}
    else:
        record = json.loads(completed.stdout)
        if completed.stderr:
            (folder / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if command == "solve-round":
        (folder / "record.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
        _write()
    kept = {key: record[key] for key in ("status", "inner_status", "reason", "iterations") if key in record}
    kept["wall_clock_s"] = record.get("wall_clock_s", time.perf_counter() - started)
    print(json.dumps(kept, ensure_ascii=False))


def _read(name: str):
    path = _OUT / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def _write() -> None:
    body, record = _read("initial.json"), _read("round1/record.json")
    lines = [
        "# T26g 验证", "",
        "求解手段：ADR 0045；解释性选择：ADR 0044。判据：初值构造正常结束，且基准外层第 1 轮内层 CONVERGED。",
        "相对变化以 ADR 0012 的固定参考尺度 ṁ_ref（燃料干基质量流）作分母；上限 1·10⁻²，至多 5 次。",
    ]
    if body:
        lines += ["", "## 焦油初值", "", "| 组分 | 床区合计 mol/s | 自由空间合计 mol/s |", "|---|---|---|"]
        lines += [f"| {row['species']} | {_num(row['bed_mol_s'])} | {_num(row['freeboard_mol_s'])} |" for row in body["tar"]]
        lines += ["", "表内是各格之和，每格相同（工况 1 的燃料加入格在床底，16 格都加焦油）。"]
        lines += ["", f"T_ad = {_num(body['t_ad_k'])} K，T_init = {_num(body['t_init_k'])} K。"
                       f"构造状态 {body['status']}，墙钟 {body['seconds']:.1f} s。"]
        if body.get("reason"):
            lines += ["", f"停止原因：{body['reason']}", f"停止阶段：{body.get('phase') or '—'}。"]
        if body.get("negative"):
            lines += [
                "", f"焦炭块（或灰块）的解出现负存量，共 {len(body['negative'])} 项；条件数 "
                f"{_num(body.get('condition'))}，块残差（解前）{_num(body.get('char_residual_before'))}。"
                "按 ADR 0045 第 4 项停止，不作近似、不截断存量。最负的 12 项：", "",
                "| 变量 | 初始 kg | 解 kg |", "|---|---|---|",
            ]
            lines += [f"| {row['name']} | {_num(row['initial_kg'])} | {_num(row['solved_kg'])} |" for row in body["negative"][:12]]
        if body["loops"]:
            lines += ["", "## 每次循环", "",
                      "| 次 | 焦炭块 ‖r‖∞/ṁ_ref | 回代 ‖r‖∞/ṁ_ref | 最小 m_C kg | 界面流相对变化 max | 备注 |",
                      "|---|---|---|---|---|---|"]
            for row in body["loops"]:
                change = row.get("change") or {}
                lines.append(
                    f"| {row['loop']} | {_num(row.get('char_residual'))} | {_num(row.get('solid_residual'))} "
                    f"| {_num(row.get('minimum_char_kg'))} | {_num(max(change.values()) if change else None)} "
                    f"| {row.get('note', '')} {row.get('domain_rejected', '')} |"
                )
            detail = [row for row in body["loops"] if row.get("up_kg_s")]
            if detail:
                lines += ["", "| 次 | 床上行 kg/s（衡算 0…N−1） | 床下行 kg/s | 自由空间回落 kg/s | 相对变化 |", "|---|---|---|---|---|"]
                for row in detail:
                    join = lambda key: "，".join(f"{value:.6e}" for value in row[key])
                    change = "，".join(f"{name}={_num(value)}" for name, value in row["change"].items())
                    lines.append(f"| {row['loop']} | {join('up_kg_s')} | {join('down_kg_s')} | {join('return_kg_s')} | {change} |")
        final = body.get("final") or body.get("rejected")
        if final:
            lines += ["", "## 固相存量" if body.get("final") else "## 停止处的固相存量", "",
                      f"床区焦炭 {_num(final['bed_char_kg'])} kg，自由空间焦炭（自由 {_num(final['free_char_kg'])}，"
                      f"悬停 {_num(final['hover_char_kg'])}）kg，灰（自由 {_num(final['free_ash_kg'])}，"
                      f"悬停 {_num(final['hover_ash_kg'])}）kg。",
                      f"最小 m_C = {_num(final['minimum_char']['m_C_kg'])} kg（{final['minimum_char']['state']}，"
                      f"位置 {final['minimum_char']['position']}）。", "",
                      "| 衡算下标 | 几何下标 | 固体 kg | 固体 m³ | 上限 m³ | 比值 |", "|---|---|---|---|---|---|"]
            lines += [
                f"| {row['cell']} | {row['geometric']} | {row['solid_kg']:.6e} | {row['volume_m3']:.6e} "
                f"| {row['limit_m3']:.6e} | {_num(row['ratio'])} |" for row in final["cells"]
            ]
            over = final.get("over_limit") or []
            lines += ["", "超过上限的格：" + ("无" if not over else "，".join(str(item) for item in over)) + "。"]
            if body.get("final") and over:
                lines.append("回代通过但固体体积超出 ADR 0026 第 5 项的上限；按该项报告并停止。")
    lines += ["", "## 外层第 1 轮", ""]
    if record:
        verdict = "第 1 轮内层以 CONVERGED 终止。" if record.get("inner_status") == "CONVERGED" else "未通过。停止，不重跑，不另试其它手段。"
        lines += [
            f"诊断记录在 `reports/t26g-validation/round1/record.json`。求解状态 {record.get('status')}。"
            f"原因：{record.get('reason') or '无'}。内层状态 {record.get('inner_status')}。"
            f"外层 {record.get('rounds')} 轮，内层迭代 {record.get('iterations')} 次，墙钟 {record.get('wall_clock_s')} s。",
            verdict,
        ]
    elif body:
        lines.append("初值构造未正常结束，按任务卡第 11 项不运行第 1 轮。" if body["status"] != "OK" else "尚未运行。")
    (_OUT / "validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv) -> None:
    calls = {
        "write": _write, "initial": initial,
        "solve-round": lambda: sys.stdout.write(json.dumps(_execute_round(), ensure_ascii=False)),
        "run-initial": lambda: _spawn("initial", 7200, _OUT / "initial-logs"),
        "round": lambda: _spawn("solve-round", SolverConfig().wall_timeout_s, _OUT / "round1"),
    }
    command = argv[1] if len(argv) > 1 else "write"
    if command not in calls:
        raise SystemExit("用法：initial | run-initial | solve-round | round | write")
    calls[command]()


if __name__ == "__main__":
    main(sys.argv)
