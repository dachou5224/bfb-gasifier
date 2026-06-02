#!/usr/bin/env python3
"""Run slimming-plan validation tiers and emit a machine-readable report.

Tier A: slimming-plan audit tests (fast, encodes pre-cleanup assumptions)
Tier B: sanity + non-slow pytest (broad regression baseline)
Tier C: run_test_sequence stages (optional, can be slow / known-broken)
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


@dataclass
class StepResult:
    name: str
    command: str
    exit_code: int
    elapsed_s: float
    note: str = ""


@dataclass
class ValidationReport:
    tier_a_slimming_audit: StepResult | None = None
    tier_b_sanity: StepResult | None = None
    tier_b_pytest_not_slow: StepResult | None = None
    tier_c_stages: list[StepResult] = field(default_factory=list)
    summary: dict[str, object] = field(default_factory=dict)


def _run(cmd: list[str], *, name: str, timeout_sec: int, note: str = "") -> StepResult:
    print(f"\n=== {name} ===")
    print("$", " ".join(cmd))
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            timeout=max(int(timeout_sec), 1),
            text=True,
        )
        code = int(proc.returncode)
    except subprocess.TimeoutExpired:
        code = 124
        note = (note + " TIMEOUT").strip()
    elapsed = time.perf_counter() - t0
    status = "PASS" if code == 0 else "FAIL"
    print(f"[{status}] exit={code} elapsed={elapsed:.1f}s {note}".strip())
    return StepResult(name=name, command=" ".join(cmd), exit_code=code, elapsed_s=elapsed, note=note)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate codebase slimming plan assumptions")
    parser.add_argument("--skip-pytest-full", action="store_true", help="Skip tier B full pytest")
    parser.add_argument("--run-stages", action="store_true", help="Run tier C run_test_sequence stages")
    parser.add_argument(
        "--stage-timeout-sec",
        type=int,
        default=300,
        help="Timeout per stage command (default: 300)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "slimming_plan_validation_latest.json",
        help="JSON report output path",
    )
    args = parser.parse_args()

    report = ValidationReport()

    report.tier_a_slimming_audit = _run(
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/test_codebase_slimming_plan.py",
            "-q",
            "--tb=short",
            "-m",
            "slimming_audit",
        ],
        name="tier_a_slimming_audit",
        timeout_sec=600,
    )

    report.tier_b_sanity = _run(
        [sys.executable, "tests/sanity_checks.py"],
        name="tier_b_sanity",
        timeout_sec=120,
    )

    if not args.skip_pytest_full:
        report.tier_b_pytest_not_slow = _run(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests/",
                "-q",
                "--tb=no",
                "-m",
                "not slow",
                "--ignore=tests/test_global_nr_solver.py",
            ],
            name="tier_b_pytest_not_slow_excluding_global_nr",
            timeout_sec=3600,
        )
        report.summary["tier_b_global_nr_solver_note"] = (
            "Run tests/test_global_nr_solver.py separately; full-suite merge may segfault in scipy path."
        )

    if args.run_stages:
        stage_ids = [
            "00_sanity",
            "10_independent_modules",
            "15_module_audits",
            "20_coupled_cell",
            "30_solver_policy",
            "35_convergence_ladder",
            "40_phase_gates",
            "50_whole_model",
        ]
        for stage_id in stage_ids:
            note = ""
            if stage_id in {"30_solver_policy", "35_convergence_ladder"}:
                note = "expected_fail_before_P0_gate_fix"
            elif stage_id == "50_whole_model":
                note = "expected_promotion_gate_when_stages_not_chained"
            report.tier_c_stages.append(
                _run(
                    [
                        sys.executable,
                        "scripts/run_test_sequence.py",
                        "--only-stage",
                        stage_id,
                        f"--timeout-sec={int(args.stage_timeout_sec)}",
                    ],
                    name=f"tier_c_stage_{stage_id}",
                    timeout_sec=max(int(args.stage_timeout_sec) * 4, 300),
                    note=note,
                )
            )

    expected_stage_failures = {
        "tier_c_stage_30_solver_policy",
        "tier_c_stage_35_convergence_ladder",
        "tier_c_stage_50_whole_model",
    }
    tier_a_ok = report.tier_a_slimming_audit.exit_code == 0
    tier_b_ok = (
        report.tier_b_sanity.exit_code == 0
        and (
            report.tier_b_pytest_not_slow is None
            or report.tier_b_pytest_not_slow.exit_code == 0
        )
    )
    unexpected_stage_failures = [
        step.name
        for step in report.tier_c_stages
        if step.exit_code != 0 and step.name not in expected_stage_failures
    ]
    expected_stage_failures_hit = [
        step.name
        for step in report.tier_c_stages
        if step.name in expected_stage_failures and step.exit_code != 0
    ]

    report.summary = {
        "slimming_plan_audit_confirmed": tier_a_ok,
        "active_core_regression_ok": tier_b_ok,
        "expected_gate_failures_before_p0": expected_stage_failures_hit,
        "unexpected_stage_failures": unexpected_stage_failures,
        "recommendation": (
            "Proceed with P0 gate/index fixes, then rerun tier C."
            if tier_a_ok and tier_b_ok
            else "Do not slim yet; fix failing tier A/B checks first."
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "summary": report.summary,
        "tier_a_slimming_audit": asdict(report.tier_a_slimming_audit),
        "tier_b_sanity": asdict(report.tier_b_sanity),
        "tier_b_pytest_not_slow": (
            None
            if report.tier_b_pytest_not_slow is None
            else asdict(report.tier_b_pytest_not_slow)
        ),
        "tier_c_stages": [asdict(step) for step in report.tier_c_stages],
    }
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nReport written to {args.output}")
    print(json.dumps(report.summary, ensure_ascii=False, indent=2))

    if not tier_a_ok:
        return 1
    if not tier_b_ok:
        return 2
    if unexpected_stage_failures:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
