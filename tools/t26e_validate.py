"""T26e 验证。失败点方向、初值雅可比，以及基准外层第 1 轮。只报告。"""
from __future__ import annotations
import json
import pickle
import sys
import time
from dataclasses import replace
from pathlib import Path
import numpy as np
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.solver.campaign import _process, config_for
from hamel.solver.driver import SolveResult, _blocks, _feasible, _references, _scales, solve_case
from hamel.solver.gates import ResidualBlocks
from hamel.solver.newton import SolveTimeout, normalized_merit
from hamel.solver.precalc import precalculate
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots
from tools.driver_report import convert
from tools.t26c_diagnosis import _list, _num
_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "t26e-validation"
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_FAILURE = _ROOT / "reports" / "t26c-diagnosis" / "baseline"
_LOG = ("gas_outflow", "char", "ash")
def _evaluate(frozen, calls):
    def evaluate(trial):
        calls["n"] += 1
        if calls["n"] % 50 == 0:
            print(f"残差评价 {calls['n']}", file=sys.stderr, flush=True)
        return evaluate_residual(np.asarray(trial, dtype=float), frozen).values
    return evaluate
def _names(frozen):
    slots = unknown_slots(frozen)
    return [slot.name for slot in slots], [slot.kind for slot in slots]
def _column(x, residual, evaluate, scales, feasible, blocks, names, kinds):
    """对数变量与逐列放大已撤回。测试替换本函数，不再调用求解器。"""
    del x, residual, evaluate, scales, feasible, blocks, names, kinds
    raise RuntimeError("对数变量与逐列放大步长已由 ADR 0043 撤回")
def _direction(frozen, blocks, x, residual, matrix, delta, mask):
    del frozen, blocks, x, residual, matrix, delta, mask
    raise RuntimeError("对数变量与逐列放大步长已由 ADR 0043 撤回")
def point() -> None:
    with np.load(_FAILURE / "failure.npz", allow_pickle=True) as saved:
        payload = {key: np.array(saved[key]) for key in saved.files}
    with (_FAILURE / "frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    names, kinds = _names(frozen)
    if names != [str(item) for item in payload["names"].tolist()]:
        raise SystemExit("失败点变量名与冻结量不一致")
    blocks = ResidualBlocks(
        payload["gas"], payload["solid"], payload["energy"],
        float(payload["n_dot_ref"]), float(payload["m_dot_ref"]), float(payload["q_dot_ref"]),
    )
    x, scales = payload["x"], payload["scales"]
    calls, started = {"n": 0}, time.perf_counter()
    evaluate = _evaluate(frozen, calls)
    fresh = evaluate(x)
    feasible = lambda trial, frozen=frozen: _feasible(trial, frozen)
    mask, matrix, delta, unresolved, enlarged = _column(x, fresh, evaluate, scales, feasible, blocks, names, kinds)
    direction = None if delta is None else _direction(frozen, blocks, x, fresh, matrix, delta, mask)
    body = {"repeated_residual_norm": normalized_merit(fresh - payload["F"], blocks), "residual_calls": calls["n"],
            "seconds": time.perf_counter() - started, "direction": direction,
            "step_norm": None if delta is None else float(np.linalg.norm(delta)), "saved_step_norm": float(np.linalg.norm(payload["dx"])),
            "log_count": int(mask.sum()), "unresolved": unresolved, "enlarged": enlarged}
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "point.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    _write()
    print(json.dumps({"step_norm": body["step_norm"], "unresolved": len(unresolved), "seconds": body["seconds"]}, ensure_ascii=False), flush=True)
def initial() -> None:
    assembled = assemble_case(_CASE, "1")
    n_ref, m_ref, q_ref = _references(assembled.physical)
    prepared = precalculate(_CASE, "1", None, adiabatic=False, iteration=1, solver=SolverConfig())
    frozen, names, kinds = prepared.frozen, *_names(prepared.frozen)
    x = pack_state(prepared.state, frozen)
    scales = _scales(frozen, n_ref, float(frozen.geometric_mass_kg.sum()))
    blocks = _blocks(frozen, n_ref, m_ref, q_ref)
    calls, started = {"n": 0}, time.perf_counter()
    evaluate = _evaluate(frozen, calls)
    residual = evaluate(x)
    mask, matrix, _delta, unresolved, enlarged = _column(
        x, residual, evaluate, scales, lambda trial, frozen=frozen: _feasible(trial, frozen), blocks, names, kinds,
    )
    singular = np.linalg.svd(matrix, compute_uv=False)
    smallest, largest = float(singular[-1]), float(singular[0])
    body = {"unknowns": int(x.size), "log_count": int(mask.sum()),
            "zero_kept_linear": int(sum(kind in _LOG and float(x[index]) == 0.0 for index, kind in enumerate(kinds))),
            "singular_min": smallest, "singular_max": largest, "condition": None if smallest == 0.0 else largest / smallest,
            "enlarged": enlarged, "unresolved": unresolved, "residual_calls": calls["n"],
            "seconds": time.perf_counter() - started, "merit": normalized_merit(residual, blocks)}
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "initial.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    _write()
    print(json.dumps({"singular_min": smallest, "condition": body["condition"], "enlarged": len(enlarged), "unresolved": len(unresolved), "seconds": body["seconds"]}, ensure_ascii=False), flush=True)
def _execute_round() -> dict:
    config = replace(config_for("baseline"), outer_round_limit=1)
    folder = _OUT / "round1"
    folder.mkdir(parents=True, exist_ok=True)
    started, end = time.perf_counter(), time.monotonic() + config.wall_timeout_s
    try:
        result = solve_case(_CASE, adiabatic=False, deadline=lambda: time.monotonic() >= end, config=config, artifact_dir=folder)
    except SolveTimeout:
        result = SolveResult("TIMEOUT", "运行超过墙钟上限")
    record = convert({"name": "baseline", "change": "不改冻结项", "status": result.status, "reason": result.reason, "rounds": result.rounds, "iterations": result.iterations, "diagnostics": result.diagnostics, "wall_clock_s": time.perf_counter() - started, "outer_round_limit": 1})
    rounds = (record.get("diagnostics") or {}).get("rounds") or []
    record["inner_status"] = rounds[0].get("inner") if rounds else None
    return record
def _spawn(command, timeout, folder) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = _process([sys.executable, "-m", "tools.t26e_validate", command], timeout)
    if completed is None:
        record = {"status": "TIMEOUT", "reason": "墙钟上限内没有返回", "inner_status": None,
                  "wall_clock_s": time.perf_counter() - started}
    elif completed.returncode != 0:
        (folder / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")
        record = {"status": "ERROR", "reason": f"子进程退出码 {completed.returncode}", "inner_status": None,
                  "stderr_log": str((folder / "stderr.log").relative_to(_ROOT))}
    else:
        record = json.loads(completed.stdout)
        if completed.stderr:
            (folder / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if command == "solve-round":
        (folder / "record.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
        _write()
    kept = {key: record[key] for key in ("status", "inner_status", "reason", "drop", "max_h", "singular_min", "condition", "iterations", "step_norm", "unresolved") if key in record}
    kept["wall_clock_s"] = record.get("wall_clock_s", time.perf_counter() - started)
    print(json.dumps(kept, ensure_ascii=False), flush=True)
def _read(name):
    path = _OUT / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
def _write() -> None:
    point_body, initial_body, record = _read("point.json"), _read("initial.json"), _read("round1/record.json")
    lines = ["# T26e 验证", "", "求解手段：ADR 0042。失败点取 `reports/t26c-diagnosis/baseline/`。探针为 y + h·Δy，F 在 x = exp(y) 处评价。"]
    if point_body:
        direction = point_body["direction"]
        lines += ["", f"重算 F 与保存 F 的归一差 = {_num(point_body['repeated_residual_norm'])}。雅可比残差评价 {point_body['residual_calls']} 次，连同方向探针用时 {point_body['seconds']:.1f} s。对数变量 {point_body['log_count']} 个。Δy 范数 {_num(point_body['step_norm'])}，保存的 x 空间步范数 {_num(point_body['saved_step_norm'])}。被放大 {len(point_body['enlarged'])} 列，到上限仍不可分辨 {len(point_body['unresolved'])} 列。", "", "## 失败点方向", ""]
        if direction is None:
            lines.append("未形成牛顿步，无方向探针。")
        else:
            fell = [row["h"] for row in direction["rows"] if row["merit"] is not None and row["merit"] < direction["merit"]]
            lines += [f"g(0) = {_num(direction['merit'])}。g'(0) = {_num(direction['slope'])}。可行的最大 h = {_num(direction['max_h'])}。该步 g = {_num(direction['max_h_merit'])}，下降量 g(0)−g(h) = {_num(direction['drop'])}。", "", "| h | 状态 | g(h) | ‖JΔy − [F(y+hΔy)−F]/h‖ / ‖JΔy‖ | 原因 |", "|---|---|---|---|---|"]
            lines += [f"| {row['h']:.0e} | {row['status']} | {_num(row['merit'])} | {_num(row['relative'])} | {row['reason']} |" for row in direction["rows"]]
            lines += ["", f"g(h) < g(0) 的步长：{_list(fell)}。"]
    lines += ["", "## 初值点"]
    if initial_body:
        lines += [
            "", f"未知量 {initial_body['unknowns']} 个，对数变量 {initial_body['log_count']} 个，值为 0 而保持线性的流量或存量 {initial_body['zero_kept_linear']} 个。归一残差 {_num(initial_body['merit'])}。雅可比 {initial_body['seconds']:.1f} s，残差评价 {initial_body['residual_calls']} 次。",
            f"最小奇异值 {_num(initial_body['singular_min'])}，最大奇异值 {_num(initial_body['singular_max'])}，条件数 {_num(initial_body['condition'])}。被放大 {len(initial_body['enlarged'])} 列，到上限仍不可分辨 {len(initial_body['unresolved'])} 列。",
            "", "| 列 | x | 起始步长 | 最终步长 | 上限 |", "|---|---|---|---|---|",
        ]
        lines += [f"| {row['name']} | {row['x']:.6e} | {row['start']:.6e} | {row['step']:.6e} | {row['cap']:.6e} |" for row in initial_body["enlarged"]]
        if initial_body["unresolved"]:
            lines += ["", "| 不可分辨列 | x | 最终步长 |", "|---|---|---|"]
            lines += [f"| {row['name']} | {row['x']:.6e} | {row['step']:.6e} |" for row in initial_body["unresolved"]]
    lines += ["", "## 外层第 1 轮", ""]
    if record:
        verdict = "第 1 轮内层以 CONVERGED 终止。" if record.get("inner_status") == "CONVERGED" else "未通过。停止，不重跑，不另试其它手段。"
        lines += [f"诊断记录在 `reports/t26e-validation/round1/record.json`。求解状态 {record.get('status')}。原因：{record.get('reason') or '无'}。内层状态 {record.get('inner_status')}。外层 {record.get('rounds')} 轮，内层迭代 {record.get('iterations')} 次，墙钟 {record.get('wall_clock_s')} s。", verdict]
    elif point_body or initial_body:
        lines.append("尚未运行。")
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
def main(argv) -> None:
    command = argv[1]
    calls = {
        "point": point, "initial": initial,
        "solve-round": lambda: sys.stdout.write(json.dumps(_execute_round(), ensure_ascii=False)),
        "round": lambda: _spawn("solve-round", SolverConfig().wall_timeout_s, _OUT / "round1"),
        "run-point": lambda: _spawn("point", 7200, _OUT / "point-logs"),
        "run-initial": lambda: _spawn("initial", 7200, _OUT / "initial-logs"),
    }
    if command not in calls:
        raise SystemExit("用法：point | initial | round | solve-round | run-point | run-initial")
    calls[command]()
if __name__ == "__main__":
    main(sys.argv)
