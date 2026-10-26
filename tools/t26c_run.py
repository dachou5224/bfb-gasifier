"""重跑基准到失败点，并对 S13、S14 做外层第 1 轮冒烟。求解在子进程里跑，到 wall_timeout_s 时终止整个进程组。"""
from __future__ import annotations
import json
import sys
import time
from dataclasses import replace
from pathlib import Path
from hamel.solver.campaign import _process, config_for
from hamel.solver.driver import SolveResult, solve_case
from hamel.solver.newton import SolveTimeout
from tools.driver_report import convert
from tools.t26c_diagnosis import diagnose
_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "t26c-diagnosis"
def _relative(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(_ROOT))
    except ValueError:
        return str(path)
def _change(name: str) -> str:
    if name == "S13":
        return "悬停样本计为床面回落"
    if name == "S14":
        return "悬停样本计为出口带出"
    return "不改冻结项"
def _worker_command(name: str, outer_limit: int | None) -> list[str]:
    command = [sys.executable, "-m", "tools.t26c_run", "solve", name]
    if outer_limit is not None:
        command.append(str(outer_limit))
    return command
def _execute(name: str, outer_limit: int | None) -> dict:
    config = config_for(name)
    if outer_limit is not None:
        config = replace(config, outer_round_limit=outer_limit)
    folder = _OUT / name
    folder.mkdir(parents=True, exist_ok=True)
    end = time.monotonic() + config.wall_timeout_s
    started = time.perf_counter()
    try:
        result = solve_case(
            _CASE, adiabatic=name == "S1", deadline=lambda: time.monotonic() >= end,
            config=config, artifact_dir=folder,
        )
    except SolveTimeout:
        result = SolveResult("TIMEOUT", "运行超过墙钟上限")
    record = convert({
        "name": name, "change": _change(name), "status": result.status, "reason": result.reason,
        "rounds": result.rounds, "iterations": result.iterations, "diagnostics": result.diagnostics,
        "wall_clock_s": time.perf_counter() - started, "smoke_outer_round_limit": outer_limit,
    })
    diagnostics = record["diagnostics"]
    for key in ("failure_npz", "frozen_path"):
        if diagnostics.get(key):
            diagnostics[key] = _relative(diagnostics[key])
    return record
def _blank(name: str, outer_limit, status: str, reason: str, started: float) -> dict:
    return {
        "name": name, "change": _change(name), "status": status, "reason": reason,
        "rounds": 0, "iterations": 0, "diagnostics": {},
        "wall_clock_s": time.perf_counter() - started, "smoke_outer_round_limit": outer_limit,
    }
def run(name: str, outer_limit: int | None) -> dict:
    config = config_for(name)
    if outer_limit is not None:
        config = replace(config, outer_round_limit=outer_limit)
    folder = _OUT / name
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = _process(_worker_command(name, outer_limit), config.wall_timeout_s)
    if completed is None:
        record = _blank(name, outer_limit, "TIMEOUT", "墙钟上限内没有返回", started)
    elif completed.returncode != 0:
        (folder / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")
        record = _blank(name, outer_limit, "ERROR", f"子进程退出码 {completed.returncode}", started)
        record["stderr_log"] = _relative(folder / "stderr.log")
    else:
        record = json.loads(completed.stdout)
    (folder / "record.json").write_text(json.dumps(record, ensure_ascii=False), encoding="utf-8")
    if name == "baseline" and record["status"] == "DAMPING_FAILED":
        diagnose(folder / "record.json", _OUT)
    return record
def main(argv) -> None:
    command = argv[1]
    if command == "solve":
        outer = None if len(argv) < 4 else int(argv[3])
        sys.stdout.write(json.dumps(_execute(argv[2], outer), ensure_ascii=False))
        return
    if command == "baseline":
        record = run("baseline", None)
    elif command == "smoke":
        record = run(argv[2], 1)
    elif command == "diagnose":
        diagnose(_OUT / "baseline" / "record.json", _OUT)
        print(json.dumps({"wrote": str(_OUT / "diagnosis.md")}, ensure_ascii=False))
        return
    else:
        raise SystemExit("用法：baseline | smoke S13|S14 | diagnose")
    print(json.dumps({"name": record["name"], "status": record["status"], "reason": record["reason"], "rounds": record["rounds"], "iterations": record["iterations"], "wall_clock_s": record["wall_clock_s"]}, ensure_ascii=False))
if __name__ == "__main__":
    main(sys.argv)
