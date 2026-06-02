#!/usr/bin/env python3
"""Run the phase-0 baseline freeze gate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.phase_gate_checks import gate_00_baseline_frozen


def main() -> int:
    parser = argparse.ArgumentParser(description="Run phase-0 baseline freeze gate")
    parser.add_argument("--json", action="store_true", help="Print JSON only")
    parser.add_argument("--strict", action="store_true", help="Treat warnings as failures")
    args = parser.parse_args()

    result = gate_00_baseline_frozen(strict=args.strict)
    payload = result.to_dict()

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(f"{payload['gate_id']} | phase={payload['phase']} | passed={payload['passed']}")
        print(
            "  sanity_passed = "
            f"{payload['metrics'].get('sanity', {}).get('passed')} | "
            f"hard_failures = {len(payload['hard_failures'])} | warnings = {len(payload['warnings'])}"
        )
        if payload["hard_failures"]:
            print("  HARD FAILURES:")
            for item in payload["hard_failures"]:
                print(f"    - {item}")
        if payload["warnings"]:
            print("  WARNINGS:")
            for item in payload["warnings"]:
                print(f"    - {item}")
        reactor = payload["metrics"].get("current_reactor_snapshot", {})
        if reactor:
            print(
                "  LU baseline: "
                f"T_exit={reactor.get('T_exit_K', float('nan')):.2f} K, "
                f"carbon_conv={reactor.get('carbon_conv', float('nan')):.4f}, "
                f"runtime={reactor.get('runtime_s', float('nan')):.2f} s"
            )

    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
