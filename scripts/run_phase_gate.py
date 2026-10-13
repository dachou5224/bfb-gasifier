#!/usr/bin/env python3
"""统一 Phase Gate runner：gate 00 / 10 / 20 / 30 / 40。

Usage
-----
    python3 scripts/run_phase_gate.py 20
    python3 scripts/run_phase_gate.py 20 --json
    python3 scripts/run_phase_gate.py all
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Callable, Dict, List, Optional, Tuple

import _lu_harness  # noqa: F401 — bootstrap repo + scripts on sys.path

from src.core.phase_gate_checks import (
    gate_00_baseline_frozen,
    gate_10_structural_parity,
    gate_20_elemental_closure,
    gate_30_phase_partition_closure,
    gate_40_axial_extent_recycle,
)

GateFn = Callable[..., Any]
PrintFn = Callable[[Dict[str, Any]], None]

GATE_ORDER = ("00", "10", "20", "30", "40")


def _print_gate_00(payload: Dict[str, Any]) -> None:
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


def _print_gate_10(payload: Dict[str, Any]) -> None:
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


def _print_gate_20(payload: Dict[str, Any]) -> None:
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


def _print_gate_30(payload: Dict[str, Any]) -> None:
    print(f"{payload['gate_id']} | phase={payload['phase']} | passed={payload['passed']}")
    ledger = payload["metrics"]["phase_partition_ledger"]
    numbering = ledger.get("reaction_numbering", {})
    if numbering:
        print(f"  reaction_label_authority = {numbering.get('authority', 'unknown')}")
    print(
        "  N_ex_sum = "
        f"{ledger['N_ex_sum']:.3e} | "
        "bubble O2 raw/lim/supply = "
        f"{ledger['o2_demand_bubble_raw']:.3e}/{ledger['o2_demand_bubble_limited']:.3e}/{ledger['o2_supply_bubble']:.3e}"
    )
    thesis_view = ledger.get("net_molar_gas_source_by_thesis", {})
    if thesis_view:
        print(
            "  thesis net sources (R7/R8/R9/R11) = "
            f"{thesis_view.get('R7', 0.0):.3e}, "
            f"{thesis_view.get('R8', 0.0):.3e}, "
            f"{thesis_view.get('R9', 0.0):.3e}, "
            f"{thesis_view.get('R11', 0.0):.3e}"
        )
    print(
        "  bubble slack / dense shortfall / transfer potential = "
        f"{ledger['o2_slack_bubble']:.3e} / {ledger['o2_shortfall_dense']:.3e} / "
        f"{ledger['o2_transfer_potential_bd']:.3e}"
    )
    if payload["hard_failures"]:
        print("  HARD FAILURES:")
        for item in payload["hard_failures"]:
            print(f"    - {item}")


def _print_gate_40(payload: Dict[str, Any]) -> None:
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


GATE_REGISTRY: Dict[str, Tuple[GateFn, PrintFn]] = {
    "00": (gate_00_baseline_frozen, _print_gate_00),
    "10": (gate_10_structural_parity, _print_gate_10),
    "20": (gate_20_elemental_closure, _print_gate_20),
    "30": (gate_30_phase_partition_closure, _print_gate_30),
    "40": (gate_40_axial_extent_recycle, _print_gate_40),
}


def run_gate(gate_id: str, *, json_only: bool = False, strict: bool = False) -> int:
    if gate_id not in GATE_REGISTRY:
        raise ValueError(f"Unknown gate {gate_id!r}; expected one of {GATE_ORDER}")

    gate_fn, print_fn = GATE_REGISTRY[gate_id]
    if gate_id == "00":
        result = gate_fn(strict=strict)
    else:
        result = gate_fn()
    payload = result.to_dict()

    if json_only:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print_fn(payload)

    return 0 if result.passed else 1


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gate",
        choices=[*GATE_ORDER, "all"],
        help="Gate id (00–40) or all sequential",
    )
    parser.add_argument("--json", action="store_true", help="Print JSON only")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Gate 00 only: treat warnings as failures",
    )
    args = parser.parse_args(argv)

    if args.gate == "all":
        exit_code = 0
        for gate_id in GATE_ORDER:
            if not args.json:
                print(f"=== gate {gate_id} ===", flush=True)
            code = run_gate(
                gate_id,
                json_only=args.json,
                strict=bool(args.strict and gate_id == "00"),
            )
            exit_code = max(exit_code, code)
        return exit_code

    return run_gate(str(args.gate), json_only=args.json, strict=args.strict)


if __name__ == "__main__":
    raise SystemExit(main())
