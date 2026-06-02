#!/usr/bin/env python3
"""BFB 大模型分阶段测试编排器（独立模块 -> 耦合模块 -> 全模型）。

设计目标：
1) 先测可独立验证模块，快速定位单点问题；
2) 再测依赖耦合与求解器策略，验证 Check1/Check2 链路；
3) 最后测全模型回归（LU / 文献回归 / phase gates）。
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Stage:
    id: str
    name: str
    commands: tuple[str, ...]


STAGES: tuple[Stage, ...] = (
    Stage(
        id="00_sanity",
        name="数量级与接线门禁",
        commands=(
            "python3 tests/sanity_checks.py",
        ),
    ),
    Stage(
        id="10_independent_modules",
        name="独立模块测试（无/弱耦合）",
        commands=(
            "python3 -m pytest tests/test_bubble_dynamics.py tests/test_mass_transfer.py tests/test_tar_reactions.py tests/test_thermal_models.py tests/test_thermodynamics.py tests/test_gibbs_minimizer.py tests/test_feed_inlet.py tests/test_correlation_priority.py -q --tb=no",
        ),
    ),
    Stage(
        id="15_module_audits",
        name="模块审计矩阵（逐模块 pytest + audit）",
        commands=(
            "python3 scripts/run_module_audits.py",
        ),
    ),
    Stage(
        id="20_coupled_cell",
        name="单元耦合测试（cell 级）",
        commands=(
            "python3 -m pytest tests/test_cell_balances.py tests/test_cell_kinetics.py tests/test_cell_pyrolysis.py tests/test_cell_solver_acceptance.py tests/test_reactor_heat_loss_distribution.py tests/test_reactor_solid_streams.py -q --tb=no",
        ),
    ),
    Stage(
        id="30_solver_policy",
        name="求解器策略与收敛路径",
        commands=(
            "python3 -m pytest "
            "tests/test_global_nr_solver.py::test_resolve_nr_init_strategy_prefers_vorabrechnung_by_default "
            "tests/test_global_nr_solver.py::test_thesis_mode_forces_global_nr_even_when_solver_arg_is_gauss_seidel "
            "tests/test_global_nr_solver.py::test_gauss_seidel_solver_is_disabled_unless_legacy_flag_enabled "
            "tests/test_global_nr_solver.py::test_gauss_seidel_solver_allowed_with_legacy_flag "
            "tests/test_global_nr_solver.py::test_non_thesis_mode_still_disables_gs_warmup_when_legacy_gs_is_off "
            "tests/test_solver_sequence.py::test_refresh_cell_vorabrechnung_runs_hydrodynamics_then_cache_update "
            "tests/test_solver_sequence.py::test_refresh_cell_vorabrechnung_force_invalidates_before_recompute "
            "-q --tb=no",
        ),
    ),
    Stage(
        id="35_convergence_ladder",
        name="收敛阶梯（cell → inner NR → outer Abgleich）",
        commands=(
            "python3 -m pytest "
            "tests/test_outer_loop_convergence.py "
            "tests/test_global_nr_solver.py::test_global_nr_outer_loop_forces_vorabrechnung_refresh "
            "tests/test_global_nr_solver.py::test_global_nr_stops_when_outer_is_matched_and_inner_stalls "
            "-q --tb=no",
            "python3 scripts/audit_global_nr_profile_lu.py --strict",
        ),
    ),
    Stage(
        id="40_phase_gates",
        name="Phase Gate 守恒与流程门禁",
        commands=(
            "python3 -m pytest "
            "tests/test_phase_gate_00.py tests/test_phase_gate_10.py tests/test_phase_gate_20.py "
            "tests/test_phase_gate_30.py tests/test_phase_gate_40.py "
            "-m 'not slow' -q --tb=no",
        ),
    ),
    Stage(
        id="50_whole_model",
        name="全模型回归（global_nr 主线）",
        commands=(
            "python3 scripts/audit_global_nr_profile_lu.py --strict",
            "python3 -m pytest tests/test_table2_LU_global_nr.py tests/test_literature_regression_priority.py tests/test_phase1.py -q --tb=no",
        ),
    ),
)

WHOLE_MODEL_PROMOTION_PREREQS: tuple[str, ...] = (
    "00_sanity",
    "10_independent_modules",
    "15_module_audits",
    "20_coupled_cell",
    "30_solver_policy",
    "35_convergence_ladder",
    "40_phase_gates",
)


def _iter_selected(
    *, from_stage: str | None, only_stage: str | None
) -> Iterable[Stage]:
    if only_stage is not None:
        picked = [s for s in STAGES if s.id == only_stage]
        if not picked:
            raise ValueError(f"Unknown --only-stage: {only_stage}")
        return picked
    if from_stage is None:
        return STAGES
    start_idx = next((i for i, s in enumerate(STAGES) if s.id == from_stage), None)
    if start_idx is None:
        raise ValueError(f"Unknown --from-stage: {from_stage}")
    return STAGES[start_idx:]


def _missing_promotion_prereqs(
    current_stage: str,
    passed_stage_ids: set[str],
) -> list[str]:
    if current_stage != "50_whole_model":
        return []
    return [sid for sid in WHOLE_MODEL_PROMOTION_PREREQS if sid not in passed_stage_ids]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run scientific step-by-step test sequence")
    parser.add_argument("--from-stage", help="Start from a stage id")
    parser.add_argument("--only-stage", help="Run only one stage id")
    parser.add_argument("--dry-run", action="store_true", help="Print commands only")
    parser.add_argument(
        "--include-slow-lu",
        action="store_true",
        help="Additionally run slow LU test (tests/test_table2_LU.py)",
    )
    parser.add_argument(
        "--include-heavy-audits",
        action="store_true",
        help="Include long-running audits (global NR profile + script phase gates)",
    )
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=180,
        help="Per-command timeout in seconds (default: 180)",
    )
    parser.add_argument(
        "--show-deprecation-warnings",
        action="store_true",
        help="Do not suppress DeprecationWarning during sequence execution",
    )
    parser.add_argument(
        "--allow-direct-whole-model",
        action="store_true",
        help="Bypass promotion gate and allow running whole-model stage without prior stages",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    selected = list(_iter_selected(from_stage=args.from_stage, only_stage=args.only_stage))
    child_env = os.environ.copy()
    if not args.show_deprecation_warnings:
        child_env.setdefault("PYTHONWARNINGS", "ignore::DeprecationWarning")
    passed_stage_ids: set[str] = set()

    for stage in selected:
        missing_prereqs = _missing_promotion_prereqs(stage.id, passed_stage_ids)
        if missing_prereqs:
            msg = (
                "[PROMOTION-GATE] 50_whole_model requires all earlier gates passed: "
                f"missing={missing_prereqs}. "
                "Run full sequence or pass --allow-direct-whole-model for manual override."
            )
            if args.allow_direct_whole_model:
                print(f"{msg} [OVERRIDDEN]")
            elif args.dry_run:
                print(f"{msg} [DRY-RUN ONLY]")
            else:
                print(msg)
                return 2

        print(f"\n=== [{stage.id}] {stage.name} ===")
        for cmd in stage.commands:
            print(f"$ {cmd}")
            if args.dry_run:
                continue
            try:
                ret = subprocess.run(
                    cmd,
                    shell=True,
                    cwd=repo_root,
                    timeout=max(int(args.timeout_sec), 1),
                    env=child_env,
                )
            except subprocess.TimeoutExpired:
                print(
                    "[TIMEOUT] "
                    f"stage={stage.id} cmd={cmd} exceeded {int(args.timeout_sec)}s; "
                    "treat as potential script problem and stop."
                )
                return 124
            if ret.returncode != 0:
                print(f"[FAIL] stage={stage.id} cmd={cmd}")
                return ret.returncode
        if stage.id == "30_solver_policy" and args.include_heavy_audits:
            cmd = "python3 scripts/audit_global_nr_profile_lu.py"
            print(f"$ {cmd}")
            if not args.dry_run:
                try:
                    ret = subprocess.run(
                        cmd,
                        shell=True,
                        cwd=repo_root,
                        timeout=max(int(args.timeout_sec), 1),
                        env=child_env,
                    )
                except subprocess.TimeoutExpired:
                    print(
                        "[TIMEOUT] "
                        f"stage={stage.id} cmd={cmd} exceeded {int(args.timeout_sec)}s; stop."
                    )
                    return 124
                if ret.returncode != 0:
                    print(f"[FAIL] stage={stage.id} cmd={cmd}")
                    return ret.returncode
        if stage.id == "40_phase_gates" and args.include_heavy_audits:
            heavy_cmds = (
                "python3 scripts/run_phase_gate_00.py --strict",
                "python3 scripts/run_phase_gate_10.py",
                "python3 scripts/run_phase_gate_20.py",
                "python3 scripts/run_phase_gate_30.py",
                "python3 scripts/run_phase_gate_40.py",
            )
            for cmd in heavy_cmds:
                print(f"$ {cmd}")
                if args.dry_run:
                    continue
                try:
                    ret = subprocess.run(
                        cmd,
                        shell=True,
                        cwd=repo_root,
                        timeout=max(int(args.timeout_sec), 1),
                        env=child_env,
                    )
                except subprocess.TimeoutExpired:
                    print(
                        "[TIMEOUT] "
                        f"stage={stage.id} cmd={cmd} exceeded {int(args.timeout_sec)}s; stop."
                    )
                    return 124
                if ret.returncode != 0:
                    print(f"[FAIL] stage={stage.id} cmd={cmd}")
                    return ret.returncode
        passed_stage_ids.add(stage.id)

    if args.include_slow_lu:
        cmd = "python3 -m pytest tests/test_table2_LU.py -q --tb=no -m slow"
        print("\n=== [60_slow_lu] 慢速文献工况全量校验 ===")
        print(f"$ {cmd}")
        if not args.dry_run:
            try:
                ret = subprocess.run(
                    cmd,
                    shell=True,
                    cwd=repo_root,
                    timeout=max(int(args.timeout_sec), 1),
                    env=child_env,
                )
            except subprocess.TimeoutExpired:
                print(
                    "[TIMEOUT] stage=60_slow_lu "
                    f"cmd exceeded {int(args.timeout_sec)}s; stop."
                )
                return 124
            if ret.returncode != 0:
                print("[FAIL] stage=60_slow_lu")
                return ret.returncode

    print("\nAll selected stages passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
