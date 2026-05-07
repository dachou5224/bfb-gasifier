#!/usr/bin/env python3
"""Run the phase-1 structural parity gate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.phase_gate_checks import gate_10_structural_parity


def main() -> int:
    parser = argparse.ArgumentParser(description="Run phase-1 structural parity gate")
    parser.add_argument("--json", action="store_true", help="Print JSON only")
    args = parser.parse_args()

    result = gate_10_structural_parity()
    payload = result.to_dict()

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(f"{payload['gate_id']} | phase={payload['phase']} | passed={payload['passed']}")
        print(
            "  hard_failures = "
            f"{len(payload['hard_failures'])} | structural_mutators = "
            f"{len(payload['metrics'].get('structural_contract', {}).get('mutators', []))}"
        )
        if payload["hard_failures"]:
            print("  HARD FAILURES:")
            for item in payload["hard_failures"]:
                print(f"    - {item}")

    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
