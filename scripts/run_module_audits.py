#!/usr/bin/env python3
"""按模块执行 pytest + 审计脚本。

目标：把“模块 -> 测试 -> 审计 -> 排障入口”收敛到单一入口，
用于分阶段链路中的 module-audit gate。
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter


@dataclass(frozen=True)
class ModuleAudit:
    module_id: str
    scope: str
    tests: tuple[str, ...]
    audits: tuple[str, ...]
    pass_rule: str
    triage_hint: str


MODULE_AUDITS: tuple[ModuleAudit, ...] = (
    ModuleAudit(
        module_id="physics",
        scope="最小流化、气泡动力学、相间传质",
        tests=(
            "python3 -m pytest tests/test_bubble_dynamics.py tests/test_mass_transfer.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_hydrodynamics_consistency_lu.py",
        ),
        pass_rule="pytest 全通过，且 hydrodynamics consistency 审计退出码为 0",
        triage_hint="src/physics/*, src/core/cell_hydrodynamics.py",
    ),
    ModuleAudit(
        module_id="kinetics",
        scope="均相/异相反应、tar 氧化与修正",
        tests=(
            "python3 -m pytest tests/test_tar_reactions.py tests/test_cell_kinetics.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_r10_pressure_correction.py",
            "python3 scripts/verify_r10_correction.py",
        ),
        pass_rule="pytest 全通过，R10 压力口径审计与当前实现验证均退出码为 0",
        triage_hint="src/kinetics/*, src/core/cell_kinetics.py",
    ),
    ModuleAudit(
        module_id="thermal",
        scope="干燥/热解与热模型",
        tests=(
            "python3 -m pytest tests/test_thermal_models.py tests/test_cell_pyrolysis.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_zone_profile_lu.py",
        ),
        pass_rule="pytest 全通过，zone profile 审计退出码为 0",
        triage_hint="src/thermal/*, src/core/cell_pyrolysis.py",
    ),
    ModuleAudit(
        module_id="thermodynamics",
        scope="平衡常数、Gibbs 最小化、微量组分",
        tests=(
            "python3 -m pytest tests/test_thermodynamics.py tests/test_gibbs_minimizer.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_global_nr_profile_lu.py",
        ),
        pass_rule="pytest 全通过，global NR profile 审计退出码为 0",
        triage_hint="src/thermodynamics/*, src/core/cell_minor_species.py",
    ),
    ModuleAudit(
        module_id="core_cell",
        scope="cell 守恒、固体流、相间交换闭合",
        tests=(
            "python3 -m pytest tests/test_cell_balances.py tests/test_conservation_priority.py tests/test_reactor_solid_streams.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_phase_exchange_conservation.py",
        ),
        pass_rule="pytest 全通过，phase exchange conservation 审计退出码为 0",
        triage_hint="src/core/cell*.py, src/core/connectivity*.py",
    ),
    ModuleAudit(
        module_id="solvers",
        scope="outer/inner 收敛语义与 NR 初始化策略",
        tests=(
            "python3 -m pytest tests/test_outer_loop_convergence.py "
            "tests/test_global_nr_solver.py::test_global_nr_outer_loop_forces_vorabrechnung_refresh "
            "tests/test_global_nr_solver.py::test_global_nr_stops_when_outer_is_matched_and_inner_stalls "
            "-q --tb=no",
        ),
        audits=(
            "python3 scripts/benchmark_nr_init_gibbs_bootstrap.py",
        ),
        pass_rule="pytest 全通过，NR init benchmark 脚本退出码为 0",
        triage_hint="src/solvers/*, src/workflow/steps/nr_inner_step.py",
    ),
    ModuleAudit(
        module_id="workflow",
        scope="流程图映射与输入口径一致性",
        tests=(
            "python3 -m pytest tests/test_flowchart_mapping_contract.py tests/test_connectivity_graph.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_flowchart_mapping.py",
            "python3 scripts/_audit_inputs.py",
        ),
        pass_rule="pytest 全通过，flowchart mapping 与输入审计均退出码为 0",
        triage_hint="src/workflow/*, docs/gasifier_model_flowchart_trilingual.mmd",
    ),
    ModuleAudit(
        module_id="freeboard",
        scope="freeboard 轨迹与反应集",
        tests=(
            "python3 -m pytest tests/test_freeboard_analytical_trajectory.py -q --tb=no",
        ),
        audits=(
            "python3 scripts/audit_freeboard_reaction_sets_lu.py",
        ),
        pass_rule="pytest 全通过，freeboard reaction-set 审计退出码为 0",
        triage_hint="src/core/freeboard_*.py, src/physics/freeboard.py",
    ),
)


def _pick_modules(module_ids: list[str]) -> list[ModuleAudit]:
    if not module_ids:
        return list(MODULE_AUDITS)
    selected: list[ModuleAudit] = []
    unknown = []
    index = {m.module_id: m for m in MODULE_AUDITS}
    for module_id in module_ids:
        module_id = module_id.strip().lower()
        if not module_id:
            continue
        mod = index.get(module_id)
        if mod is None:
            unknown.append(module_id)
            continue
        selected.append(mod)
    if unknown:
        raise ValueError(f"Unknown module ids: {unknown}")
    return selected


def _run_cmd(cmd: str, *, cwd: Path, timeout_sec: int, dry_run: bool) -> tuple[int, float]:
    print(f"$ {cmd}")
    if dry_run:
        return 0, 0.0
    t0 = perf_counter()
    try:
        ret = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            timeout=max(int(timeout_sec), 1),
        )
    except subprocess.TimeoutExpired:
        return 124, perf_counter() - t0
    return int(ret.returncode), perf_counter() - t0


def main() -> int:
    parser = argparse.ArgumentParser(description="Run per-module tests and audits")
    parser.add_argument(
        "--module",
        action="append",
        default=[],
        help="Module id to run; repeat for multiple (default: all)",
    )
    parser.add_argument("--list", action="store_true", help="Print module matrix and exit")
    parser.add_argument("--skip-tests", action="store_true", help="Skip pytest commands")
    parser.add_argument("--skip-audits", action="store_true", help="Skip audit scripts")
    parser.add_argument("--dry-run", action="store_true", help="Print commands only")
    parser.add_argument("--json", action="store_true", help="Print JSON summary at end")
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=180,
        help="Per-command timeout in seconds (default: 180)",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent

    if args.list:
        for mod in MODULE_AUDITS:
            print(f"[{mod.module_id}] {mod.scope}")
            print(f"  tests : {len(mod.tests)} command(s)")
            print(f"  audits: {len(mod.audits)} command(s)")
            print(f"  pass  : {mod.pass_rule}")
            print(f"  triage: {mod.triage_hint}")
        return 0

    modules = _pick_modules(args.module)
    summary: dict[str, dict[str, object]] = {}

    for mod in modules:
        print(f"\n=== [module:{mod.module_id}] {mod.scope} ===")
        module_ok = True
        module_elapsed = 0.0
        commands_total = 0

        if not args.skip_tests:
            for cmd in mod.tests:
                commands_total += 1
                code, elapsed = _run_cmd(
                    cmd,
                    cwd=repo_root,
                    timeout_sec=args.timeout_sec,
                    dry_run=args.dry_run,
                )
                module_elapsed += elapsed
                if code != 0:
                    module_ok = False
                    print(
                        f"[FAIL] module={mod.module_id} stage=tests code={code} "
                        f"triage={mod.triage_hint}"
                    )
                    break

        if module_ok and not args.skip_audits:
            for cmd in mod.audits:
                commands_total += 1
                code, elapsed = _run_cmd(
                    cmd,
                    cwd=repo_root,
                    timeout_sec=args.timeout_sec,
                    dry_run=args.dry_run,
                )
                module_elapsed += elapsed
                if code != 0:
                    module_ok = False
                    print(
                        f"[FAIL] module={mod.module_id} stage=audits code={code} "
                        f"triage={mod.triage_hint}"
                    )
                    break

        summary[mod.module_id] = {
            "scope": mod.scope,
            "passed": bool(module_ok),
            "commands": int(commands_total),
            "elapsed_s": float(module_elapsed),
            "pass_rule": mod.pass_rule,
            "triage_hint": mod.triage_hint,
        }
        if not module_ok:
            if args.json:
                print(json.dumps({"passed": False, "modules": summary}, ensure_ascii=False, indent=2))
            return 1

    print("\nAll selected module audits passed.")
    if args.json:
        print(json.dumps({"passed": True, "modules": summary}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
