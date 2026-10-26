"""T26d 验证。失败点重算雅可比与牛顿步，并跑基准外层第 1 轮。只报告。"""
from __future__ import annotations

import json
import pickle
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.solver.campaign import _process, config_for
from hamel.solver.driver import SolveResult, _feasible, solve_case
from hamel.solver.gates import ResidualBlocks
from hamel.solver.newton import SolveTimeout, _SQRT_EPS, _difference_basis, _jacobian, _newton_step, normalized_merit
from hamel.solver.residual import evaluate_residual
from tools.driver_report import convert
from tools.t26c_diagnosis import _direction, _direction_sentence, _list, _num, _probe

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "t26d-validation"
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_FAILURE = _ROOT / "reports" / "t26c-diagnosis" / "baseline"
_SMALL = 1.0e-8
_FACTORS = (1.0, 1.0e2, 1.0e4)
_HOVER = (
    "freeboard[b=7,geo=0].char.hover.k3",
    "freeboard[b=7,geo=0].ash.hover.k3",
    "freeboard[b=6,geo=1].char.hover.k3",
    "freeboard[b=6,geo=1].ash.hover.k3",
)


def _load_point():
    with np.load(_FAILURE / "failure.npz", allow_pickle=True) as saved:
        payload = {key: np.array(saved[key]) for key in saved.files}
    with (_FAILURE / "frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    blocks = ResidualBlocks(
        payload["gas"], payload["solid"], payload["energy"],
        float(payload["n_dot_ref"]), float(payload["m_dot_ref"]), float(payload["q_dot_ref"]),
    )
    return payload, frozen, blocks, [str(item) for item in payload["names"].tolist()]


def _evaluate(frozen, calls):
    def evaluate(trial):
        calls["n"] += 1
        if calls["n"] % 50 == 0:
            print(f"残差评价 {calls['n']}", file=sys.stderr, flush=True)
        return evaluate_residual(np.asarray(trial, dtype=float), frozen).values
    return evaluate


def _rel(column, base):
    if column is None or base is None:
        return None
    denom = float(np.linalg.norm(base))
    if denom == 0.0:
        return 0.0 if float(np.linalg.norm(column)) == 0.0 else None
    return float(np.linalg.norm(np.asarray(column) - np.asarray(base)) / denom)


def _hover_rows(frozen, x, residual, scales, matrix, saved, names) -> list:
    rows = []
    for name in _HOVER:
        index = names.index(name)
        step = _SQRT_EPS * _difference_basis(x[index], scales[index])
        base = matrix[:, index]
        relatives, status = {}, {}
        for factor in _FACTORS:
            point = np.array(x, dtype=float, copy=True)
            point[index] += step * factor
            values, state, reason = _probe(frozen, point)
            column = None if state != "ok" else (values - residual) / (step * factor)
            relatives[str(factor)] = _rel(column, base)
            status[str(factor)] = state if state == "ok" else f"{state}:{reason}"
        order = np.argsort(-np.abs(base))[:5]
        rows.append({
            "name": name, "x": float(x[index]), "scale": float(scales[index]), "step": float(step),
            "old_step": float(_SQRT_EPS * max(abs(float(x[index])), float(scales[index]))),
            "column_norm": float(np.linalg.norm(base)),
            "saved_column_norm": float(np.linalg.norm(saved[:, index])),
            "relative_to_saved": _rel(base, saved[:, index]),
            "relative": relatives, "status": status, "column": base,
            "largest": [{"name": names[int(item)], "value": float(base[int(item)])} for item in order],
        })
    return rows

def _small_descent(direction) -> dict:
    merit = direction["merit"]
    ok, fell = [], []
    for row in direction["rows"]:
        if row["h"] > _SMALL or row["status"] != "ok" or row["merit"] is None:
            continue
        ok.append(row["h"])
        if row["merit"] < merit:
            fell.append(row["h"])
    return {"ok": ok, "decreased": fell, "descends": bool(fell)}

def point() -> dict:
    payload, frozen, blocks, names = _load_point()
    x, scales = payload["x"], payload["scales"]
    calls = {"n": 0}
    evaluate = _evaluate(frozen, calls)
    started = time.perf_counter()
    fresh = evaluate(x)
    matrix = _jacobian(x, fresh, evaluate, scales, lambda trial, frozen=frozen: _feasible(trial, frozen))
    delta = _newton_step(matrix, fresh, x, scales, blocks)
    direction = _direction(frozen, blocks, x, fresh, matrix, delta)
    hover = _hover_rows(frozen, x, fresh, scales, matrix, payload["J"], names)
    body = {
        "repeated_residual_norm": normalized_merit(fresh - payload["F"], blocks),
        "residual_calls": calls["n"], "seconds": time.perf_counter() - started,
        "direction": direction, "small_step": _small_descent(direction),
        "step_norm": float(np.linalg.norm(delta)), "saved_step_norm": float(np.linalg.norm(payload["dx"])),
        "hover": [{key: value for key, value in row.items() if key != "column"} for row in hover],
    }
    _OUT.mkdir(parents=True, exist_ok=True)
    np.savez(
        _OUT / "point.npz", x=x, F=fresh, J=matrix, dx=delta, scales=scales,
        hover_names=np.array([row["name"] for row in hover]),
        hover_columns=np.column_stack([row["column"] for row in hover]),
    )
    (_OUT / "point.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    (_OUT / "validation.md").write_text(_markdown(body, None), encoding="utf-8")
    print(json.dumps({"descends": body["small_step"]["descends"], "decreased": body["small_step"]["decreased"],
                      "seconds": body["seconds"], "calls": calls["n"]}, ensure_ascii=False), flush=True)
    return body


def _execute_round() -> dict:
    config = replace(config_for("baseline"), outer_round_limit=1)
    folder = _OUT / "round1"
    folder.mkdir(parents=True, exist_ok=True)
    end = time.monotonic() + config.wall_timeout_s
    started = time.perf_counter()
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

def run_round() -> dict:
    folder = _OUT / "round1"
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = _process([sys.executable, "-m", "tools.t26d_validate", "solve-round"], SolverConfig().wall_timeout_s)
    if completed is None:
        record = {"name": "baseline", "status": "TIMEOUT", "reason": "墙钟上限内没有返回", "inner_status": None,
                  "rounds": 0, "iterations": 0, "diagnostics": {}, "wall_clock_s": time.perf_counter() - started}
    elif completed.returncode != 0:
        (folder / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")
        record = {"name": "baseline", "status": "ERROR", "reason": f"子进程退出码 {completed.returncode}",
                  "inner_status": None, "rounds": 0, "iterations": 0, "diagnostics": {},
                  "stderr_log": "reports/t26d-validation/round1/stderr.log",
                  "wall_clock_s": time.perf_counter() - started}
    else:
        record = json.loads(completed.stdout)
    (folder / "record.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    body = json.loads((_OUT / "point.json").read_text(encoding="utf-8"))
    (_OUT / "validation.md").write_text(_markdown(body, record), encoding="utf-8")
    print(json.dumps({"status": record["status"], "inner_status": record.get("inner_status"),
                      "reason": record.get("reason"), "iterations": record.get("iterations"),
                      "wall_clock_s": record.get("wall_clock_s")}, ensure_ascii=False), flush=True)
    return record


def _markdown(body, record) -> str:
    direction, small = body["direction"], body["small_step"]
    lines = [
        "# T26d 验证", "",
        "失败点取 `reports/t26c-diagnosis/baseline/failure.npz` 与同目录 `frozen.pkl`。雅可比与牛顿步按 ADR 0041 重算。方向检查的 h 与 T26c 第 6 项相同。",
        f"小步长指 h ≤ {_SMALL:.0e}。这些 h 里至少有一个可行且 g(h) < g(0)，才记为小步长上的下降方向。", "",
        f"重算 F 与保存 F 的归一差 = {_num(body['repeated_residual_norm'])}。残差评价 {body['residual_calls']} 次，用时 {body['seconds']:.1f} s。",
        f"新牛顿步范数 {_num(body['step_norm'])}，保存的牛顿步范数 {_num(body['saved_step_norm'])}。", "",
        "## 方向", "",
        f"g(0) = {_num(direction['merit'])}。g'(0) = {_num(direction['slope'])}。", "",
        "| h | 状态 | g(h) | ‖JΔx − [F(x+hΔx)−F(x)]/h‖ / ‖JΔx‖ | 原因 |",
        "|---|---|---|---|---|",
    ]
    for row in direction["rows"]:
        lines.append(f"| {row['h']:.0e} | {row['status']} | {_num(row['merit'])} | {_num(row['relative'])} | {row['reason']} |")
    decreased = [row["h"] for row in direction["rows"] if row["merit"] is not None and row["merit"] < direction["merit"]]
    agreed = [row["h"] for row in direction["rows"] if row["relative"] is not None and row["relative"] < 1.0e-2]
    lines += [
        "", f"全部已评价 h 中 g(h) < g(0) 的步长：{_list(decreased)}。相对差小于 1e-02 的步长：{_list(agreed)}。",
        _direction_sentence(direction, decreased, agreed),
        f"小步长上状态为 ok 的 h：{_list(small['ok'])}。其中 g(h) < g(0) 的 h：{_list(small['decreased'])}。",
        f"小步长上是下降方向：{small['descends']}。", "",
        "## 悬停焦炭与悬停灰", "",
        "列是失败点上实际存在的悬停未知量。新步长为 √ε·|x|。旧步长为 √ε·max(|x|, s)。",
        "1×、10²×、10⁴× 是新步长的倍数，相对差以新步长的差分列为分母。完整列在 `reports/t26d-validation/point.npz` 的 `hover_columns`。", "",
        "| 列 | x | 新步长 | 旧步长 | 新列范数 | 与保存列的相对差 | 1× | 10²× | 10⁴× |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in body["hover"]:
        relative = row["relative"]
        lines.append(
            f"| {row['name']} | {row['x']:.6e} | {row['step']:.6e} | {row['old_step']:.6e} | {row['column_norm']:.6e} | "
            f"{_num(row['relative_to_saved'])} | {_num(relative.get('1.0'))} | {_num(relative.get('100.0'))} | {_num(relative.get('10000.0'))} |"
        )
    lines.append("")
    for row in body["hover"]:
        text = "；".join(f"{item['name']} = {item['value']:.6e}" for item in row["largest"])
        lines.append(f"- {row['name']} 绝对值最大的分量：{text}")
    lines += ["", "## 外层第 1 轮", ""]
    damping = False
    if record is None:
        lines.append("尚未运行。第 1 轮还没有结果，不进入重跑。")
    else:
        damping = record.get("inner_status") == "DAMPING_FAILED" or record.get("status") == "DAMPING_FAILED"
        lines.append("`outer_round_limit = 1` 只用于取出第 1 轮，不改默认 100。诊断记录在 `reports/t26d-validation/round1/record.json`。")
        lines.append(
            f"求解状态 {record.get('status')}。原因：{record.get('reason') or '无'}。内层状态 {record.get('inner_status')}。"
            f"外层 {record.get('rounds')} 轮，内层迭代 {record.get('iterations')} 次，墙钟 {record.get('wall_clock_s')} s。"
        )
        if not small["descends"] or damping:
            reason = "小步长上方向不是下降方向。" if not small["descends"] else ""
            reason += "第 1 轮以内层 DAMPING_FAILED 终止。" if damping else ""
            lines.append("停止，不重跑。" + reason)
        elif record.get("status") in {"ERROR", "TIMEOUT"}:
            lines.append(f"第 1 轮没有正常返回（{record.get('status')}），不进入重跑。")
        elif record.get("inner_status") != "CONVERGED":
            lines.append(f"停止，不重跑。第 1 轮终止状态为 {record.get('status')}：{record.get('reason') or '无'}。")
        else:
            lines.append("方向检查与第 1 轮都没有触发停止条件。")
    lines.append("")
    return "\n".join(lines)

def main(argv) -> None:
    command = argv[1]
    if command == "point":
        point()
    elif command == "solve-round":
        sys.stdout.write(json.dumps(_execute_round(), ensure_ascii=False))
    elif command == "round":
        run_round()
    else:
        raise SystemExit("用法：point | round | solve-round")


if __name__ == "__main__":
    main(sys.argv)
