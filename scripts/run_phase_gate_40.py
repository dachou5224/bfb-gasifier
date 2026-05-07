#!/usr/bin/env python3
"""Run the phase-4 axial extent and recycle gate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.phase_gate_checks import gate_40_axial_extent_recycle


def main() -> int:
    parser = argparse.ArgumentParser(description="Run phase-4 axial extent/recycle gate")
    parser.add_argument("--json", action="store_true", help="Print JSON only")
    args = parser.parse_args()

    result = gate_40_axial_extent_recycle()
    payload = result.to_dict()

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(f"{payload['gate_id']} | phase={payload['phase']} | passed={payload['passed']}")
        extent = payload["metrics"]["extent_audit"]
        rows = extent.get("rows", [])
        top_suffix = extent.get("char_dominant_from_cell")
        top_suffix_len = 0
        if top_suffix is not None and rows:
            top_suffix_len = max(len(rows) - int(top_suffix), 0)
        print(
            "  vm_done_cell = "
            f"{extent.get('vm_done_cell')} | moisture_done_cell = {extent.get('moisture_done_cell')} | "
            f"char_dominant_from = {extent.get('char_dominant_from_cell')} (suffix_len={top_suffix_len})"
        )
        print(
            "  rebound_release(vm/moist) = "
            f"{extent.get('vm_release_rebound', float('nan')):.3e}/{extent.get('moist_release_rebound', float('nan')):.3e} | "
            "rebound_stock(vm/moist) = "
            f"{extent.get('vm_rebound', float('nan')):.3e}/{extent.get('moist_rebound', float('nan')):.3e} | "
            f"converged={extent.get('converged')} n_iter={extent.get('n_iter')}"
        )
        if payload["hard_failures"]:
            print("  HARD FAILURES:")
            for item in payload["hard_failures"]:
                print(f"    - {item}")

    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
