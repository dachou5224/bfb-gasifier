"""ADR 0030 的基准与 S1–S16。不能通过现有接口表达的项记 NOT_APPLIED。"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, replace
from pathlib import Path

from hamel.cases.config import SolverConfig
from hamel.solver.driver import SolveResult, solve_case
from hamel.solver.newton import SolveTimeout
from tools.driver_report import convert, write_report

_ARITHMETIC_M = (6.45e-5, 1.75e-4, 3.75e-4, 6.05e-4, 8.55e-4, 1.5e-3, 2.575e-3)
_ROWS = (
    ("baseline", "不改冻结项", "ADR 0030", True, ""),
    ("S1", "k_W = 0，初值不扣除壁面散热", "ADR 0030 F19；ADR 0029", True, ""),
    ("S2", "第二块布风板 n_0 = 95.6", "ADR 0019", True, ""),
    ("S3", "ε_P = 0.79", "ADR 0017 第 17 项", True, ""),
    ("S4", "ε_P = 0.91", "ADR 0017 第 17 项", True, ""),
    ("S5", "D_g 改为 O2–N2 的 Fuller 二元值", "ADR 0021 第 4 项", True, ""),
    ("S6", "七个粒级改用 ADR 0023 的算术平均径", "ADR 0023；ADR 0030", True, ""),
    ("S7", "最细档代表径改为 0.0328 mm", "ADR 0023", True, ""),
    ("S8", "最细档代表径改为 0.0753 mm", "ADR 0023", True, ""),
    ("S9", "ε_mf = 0.40，φ_S = 0.6", "ADR 0024", True, ""),
    ("S10", "ε_mf = 0.40，φ_S = 0.9", "ADR 0024", True, ""),
    ("S11", "ε_mf = 0.55，φ_S = 0.6", "ADR 0024", True, ""),
    ("S12", "ε_mf = 0.55，φ_S = 0.9", "ADR 0024", True, ""),
    ("S13", "悬停样本计为床面回落", "ADR 0026 第 6 项；ADR 0038", True, ""),
    ("S14", "悬停样本计为出口带出", "ADR 0026 第 6 项；ADR 0038", True, ""),
    ("S15", "Serio 第二对参数", "ADR 0013 第 4 项", True, ""),
    ("S16", "cell 数加倍", "ADR 0013 第 7 项", True, ""),
)

_BED_PAIRS = {"S9": (0.40, 0.6), "S10": (0.40, 0.9), "S11": (0.55, 0.6), "S12": (0.55, 0.9)}


@dataclass(frozen=True)
class RunSpec:
    name: str
    change: str
    source: str
    wired: bool
    reason: str = ""


def run_specs() -> tuple:
    order = ["baseline", *[f"S{index}" for index in range(1, 17)]]
    by_name = {row[0]: row for row in _ROWS}
    return tuple(RunSpec(*by_name[name]) for name in order)


def overlay(source, dest, spec: RunSpec) -> bool:
    """把这一项写进工况副本。返回是否绝热。不改原文件。"""
    doc = json.loads(Path(source).read_text(encoding="utf-8"))
    if spec.name == "S1":
        Path(dest).write_text(json.dumps(doc), encoding="utf-8")
        return True
    if spec.name == "S2":
        _assign(doc, "distributor_orifice_count", 95.6)
    elif spec.name == "S6":
        _assign(doc, "fuel_size_class_diameters", list(_ARITHMETIC_M))
    elif spec.name in {"S7", "S8"}:
        diameters = list(_locate(doc, "fuel_size_class_diameters")["value"])
        diameters[0] = 3.28e-5 if spec.name == "S7" else 7.53e-5
        _assign(doc, "fuel_size_class_diameters", diameters)
    elif spec.name in _BED_PAIRS:
        voidage, sphericity = _BED_PAIRS[spec.name]
        _assign(doc, "bed_minimum_fluidization_voidage", voidage)
        _assign(doc, "bed_sphericity", sphericity)
    elif spec.name == "S16":
        _assign(doc, "cell_heights", [0.075] * 16 + [0.14375] * 16)
    elif spec.name not in {"S3", "S4", "S5", "S13", "S14", "S15"}:
        raise ValueError(spec.name)
    Path(dest).write_text(json.dumps(doc), encoding="utf-8")
    return False


def run_campaign(source, report_dir, *, timeout_s: float | None = None) -> list:
    """依次跑基准和 S1–S16。超时记 TIMEOUT。未接线的项不假装跑过。"""
    timeout = float(SolverConfig().wall_timeout_s if timeout_s is None else timeout_s)
    records = []
    with tempfile.TemporaryDirectory() as directory:
        for spec in run_specs():
            if not spec.wired:
                records.append(_absent(spec))
                continue
            dest = Path(directory) / f"{spec.name}.json"
            adiabatic = overlay(source, dest, spec) if spec.name != "baseline" else False
            path = source if spec.name == "baseline" else dest
            records.append(_spawn(spec, path, adiabatic, timeout))
    write_report(report_dir, records, title="第 1 关运行", toy=False)
    return records


def _spawn(spec, path, adiabatic, timeout) -> dict:
    completed = _process(
        [sys.executable, "-m", "hamel.solver.campaign", spec.name, str(path), "1" if adiabatic else "0"],
        timeout,
    )
    if completed is None:
        return _absent(spec, "TIMEOUT", "墙钟上限内没有返回")
    if completed.returncode != 0:
        return _absent(spec, "NOT_CONVERGED", (completed.stderr or "")[-500:])
    return json.loads(completed.stdout)


def _process(args, timeout):
    proc = subprocess.Popen(
        args, cwd=_root(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True,
    )
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        return None
    return subprocess.CompletedProcess(args, proc.returncode, out, err)


def execute(name: str, path: str, adiabatic: bool) -> dict:
    spec = next(item for item in run_specs() if item.name == name)
    end = time.monotonic() + SolverConfig().wall_timeout_s
    try:
        result = solve_case(path, adiabatic=adiabatic, deadline=lambda: time.monotonic() >= end, config=config_for(name))
    except SolveTimeout:
        result = SolveResult("TIMEOUT", "运行超过墙钟上限")
    return convert({
        "name": spec.name, "change": spec.change, "source": spec.source, "status": result.status,
        "reason": result.reason, "rounds": result.rounds, "iterations": result.iterations,
        "ratios": result.ratios, "balance": result.balance, "nonnegative": result.nonnegative,
        "measurements": result.measurements, "diagnostics": result.diagnostics,
    })


def _absent(spec: RunSpec, status: str = "NOT_APPLIED", reason: str | None = None) -> dict:
    return {
        "name": spec.name, "change": spec.change, "source": spec.source, "status": status,
        "reason": spec.reason if reason is None else reason, "measurements": None, "balance": None,
        "nonnegative": None, "ratios": {}, "diagnostics": {}, "rounds": 0, "iterations": 0,
    }


def config_for(name: str) -> SolverConfig:
    """敏感性只改这份副本。默认字段等于冻结值。"""
    base = SolverConfig()
    if name == "S3":
        return replace(base, particle_emissivity=0.79)
    if name == "S4":
        return replace(base, particle_emissivity=0.91)
    if name == "S5":
        return replace(base, interphase_diffusivity="fuller_o2_n2")
    if name == "S13":
        return replace(base, hover_closure="B")
    if name == "S14":
        return replace(base, hover_closure="C")
    if name == "S15":
        return replace(base, serio_frequency_per_s=9.96e2, serio_activation_j_mol=64.5e3)
    return base


def _assign(doc, name, value) -> None:
    _locate(doc, name)["value"] = value


def _locate(doc, name):
    found = []

    def walk(node) -> None:
        if isinstance(node, dict):
            child = node.get(name)
            if isinstance(child, dict) and "value" in child:
                found.append(child)
            for item in node.values():
                walk(item)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(doc)
    if len(found) != 1:
        raise ValueError(f"{name} 出现 {len(found)} 次")
    return found[0]


def _root() -> Path:
    return Path(__file__).resolve().parents[2]


def main(argv) -> None:
    json.dump(execute(argv[1], argv[2], argv[3] == "1"), sys.stdout, ensure_ascii=False)


if __name__ == "__main__":
    main(sys.argv)
