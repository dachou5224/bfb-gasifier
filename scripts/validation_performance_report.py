#!/usr/bin/env python3
"""生成 LU 基线的验证性能报告（稳定性 + 精度 + 耗时）。"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from time import perf_counter
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tests.validation_case_utils import validation_numeric_tolerances


def _load_phase1_run_audit():
    mod_path = REPO_ROOT / "scripts" / "audit_phase1_htw_lu.py"
    spec = importlib.util.spec_from_file_location("audit_phase1_htw_lu", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module spec: {mod_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod.run_audit


def _max_species_rel_err(species_err: dict[str, float]) -> float:
    if not species_err:
        return 0.0
    return max(float(v) for v in species_err.values())


def build_report(*, strict: bool) -> tuple[dict[str, Any], bool]:
    t0 = perf_counter()
    run_audit = _load_phase1_run_audit()
    out, pass_validation, warnings = run_audit(strict=False)
    elapsed = perf_counter() - t0
    tol = validation_numeric_tolerances()

    temp_rel_err = float(out.get("err_T_vs_json_best") or 0.0)
    species_err = dict(out.get("species_relative_error") or {})
    carbon_rel = out.get("carbon_relative_error")
    carbon_rel = None if carbon_rel is None else float(carbon_rel)
    converged = bool(out.get("converged", False))
    n_iter = int(out.get("n_iter", 0) or 0)

    max_species_rel = _max_species_rel_err(species_err)
    thresholds = {
        "rtol_T": float(tol["rtol_T"]),
        "rtol_species_main": float(tol["rtol_CO_CO2_H2"]),
        "rtol_CH4": float(tol["rtol_CH4"]),
        "rtol_carbon_conv": float(tol["rtol_carbon_conv"]),
    }

    score = {
        "stability": {
            "converged": converged,
            "n_iter": n_iter,
        },
        "accuracy": {
            "temp_rel_err": temp_rel_err,
            "max_species_rel_err": max_species_rel,
            "carbon_rel_err": carbon_rel,
            "validation_pass": bool(pass_validation),
        },
        "runtime": {
            "wall_time_s": float(elapsed),
        },
    }

    report = {
        "case_key": out.get("case_key"),
        "solver": out.get("solver"),
        "strict_mode": bool(strict),
        "thresholds": thresholds,
        "score": score,
        "warnings": warnings,
        "raw": out,
    }

    strict_ok = True
    if strict:
        strict_ok = bool(pass_validation and converged)
    return report, strict_ok


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate validation performance report")
    parser.add_argument(
        "--output",
        default=str(REPO_ROOT / "data" / "validation_performance_latest.json"),
        help="Output JSON path (default: data/validation_performance_latest.json)",
    )
    parser.add_argument("--strict", action="store_true", help="Fail when validation does not pass or solver not converged")
    parser.add_argument("--print-json", action="store_true", help="Print report JSON to stdout")
    args = parser.parse_args()

    report, ok = build_report(strict=args.strict)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("Validation performance report generated")
    print(f"  output: {output_path}")
    print(f"  converged: {report['score']['stability']['converged']}  n_iter: {report['score']['stability']['n_iter']}")
    print(
        "  errors: "
        f"T_rel={report['score']['accuracy']['temp_rel_err']:.4f}, "
        f"species_max_rel={report['score']['accuracy']['max_species_rel_err']:.4f}, "
        f"carbon_rel={report['score']['accuracy']['carbon_rel_err']}"
    )
    print(f"  wall_time_s: {report['score']['runtime']['wall_time_s']:.2f}")
    print(f"  validation_pass: {report['score']['accuracy']['validation_pass']}")

    if args.print_json:
        print(json.dumps(report, ensure_ascii=False, indent=2))

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
