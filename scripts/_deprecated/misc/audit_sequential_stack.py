#!/usr/bin/env python3
"""顺序化审计入口：子模块 -> 集成 cell -> 整炉 reactor。

默认遇到失败即停止；可用 ``--continue-on-failure`` 继续跑完后续阶段。
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


STAGES: list[tuple[str, list[list[str]]]] = [
    (
        "submodules",
        [
            [sys.executable, "tests/sanity_checks.py"],
            [sys.executable, "-m", "pytest", "-q", "tests/test_cell_balances.py", "tests/test_cell_pyrolysis.py", "tests/test_thermal_models.py"],
            [sys.executable, "scripts/audit_cell_submodels.py"],
            [sys.executable, "scripts/audit_species_index_alignment.py"],
        ],
    ),
    (
        "integrated-cell",
        [
            [sys.executable, "scripts/audit_cell_conservation.py"],
            [sys.executable, "scripts/audit_phase_exchange_conservation.py"],
            [sys.executable, "scripts/audit_single_cell_stiff_solver.py"],
        ],
    ),
    (
        "whole-reactor",
        [
            [sys.executable, "-m", "pytest", "-q", "tests/test_table2_LU.py"],
            [sys.executable, "scripts/audit_phase1_htw_lu.py"],
            [sys.executable, "scripts/audit_zone_profile_lu.py"],
            [sys.executable, "scripts/audit_drying_pyrolysis_extent_lu.py"],
            [sys.executable, "scripts/audit_temperature_profile_lu.py"],
        ],
    ),
]


def run_cmd(cmd: list[str]) -> int:
    print(f"$ {' '.join(cmd)}")
    completed = subprocess.run(cmd, cwd=REPO_ROOT)
    return int(completed.returncode)


def main() -> int:
    ap = argparse.ArgumentParser(description="按固定顺序运行 BFB 审计栈")
    ap.add_argument(
        "--continue-on-failure",
        action="store_true",
        help="某阶段失败后继续执行后续阶段",
    )
    args = ap.parse_args()

    print("=" * 88)
    print("BFB audit stack")
    print("Order: submodules -> integrated cell -> whole reactor")
    print("=" * 88)

    overall_rc = 0
    for stage_name, commands in STAGES:
        print(f"\n[stage] {stage_name}")
        print("-" * 88)
        for cmd in commands:
            rc = run_cmd(cmd)
            if rc != 0:
                overall_rc = rc
                print(f"[stage:{stage_name}] failed with rc={rc}")
                if not args.continue_on_failure:
                    return rc
                break
        else:
            print(f"[stage:{stage_name}] ok")

    if overall_rc == 0:
        print("\nAll audit stages completed.")
    return overall_rc


if __name__ == "__main__":
    raise SystemExit(main())
