#!/usr/bin/env python3
"""Run the phase-2 elemental closure gate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.phase_gate_checks import gate_20_elemental_closure


def main() -> int:
    parser = argparse.ArgumentParser(description="Run phase-2 elemental closure gate")
    parser.add_argument("--json", action="store_true", help="Print JSON only")
    args = parser.parse_args()

    result = gate_20_elemental_closure()
    payload = result.to_dict()

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(f"{payload['gate_id']} | phase={payload['phase']} | passed={payload['passed']}")
        numbering = payload["metrics"]["reaction"].get("reaction_numbering", {})
        if numbering:
            print(f"  reaction_label_authority = {numbering.get('authority', 'unknown')}")
        print(
            "  pyrolysis_max_abs = "
            f"{payload['metrics']['pyrolysis']['max_abs_closure']:.3e} | "
            "reaction_max_abs = "
            f"{payload['metrics']['reaction']['max_abs_closure']:.3e}"
        )
        thesis_view = payload["metrics"]["reaction"].get("net_molar_gas_source_by_thesis", {})
        if thesis_view:
            print(
                "  thesis net sources (R7/R8/R9/R11) = "
                f"{thesis_view.get('R7', 0.0):.3e}, "
                f"{thesis_view.get('R8', 0.0):.3e}, "
                f"{thesis_view.get('R9', 0.0):.3e}, "
                f"{thesis_view.get('R11', 0.0):.3e}"
            )
        if payload["hard_failures"]:
            print("  HARD FAILURES:")
            for item in payload["hard_failures"]:
                print(f"    - {item}")

    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
