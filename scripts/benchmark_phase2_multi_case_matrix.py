#!/usr/bin/env python3
"""Phase2 多 case 验收矩阵：零 diagnostic scale，三分支趋势 baseline。

矩阵定义见 ``docs/validation_cases_chemistry_map.md`` §7：
  - CASE_HTW_WESSELING_1   (Air+Steam 基准)
  - CASE_HTW_WESSELING_2   (O₂+Steam 对照)
  - CASE_VTT_PRESSURIZED_PEAT_6 (高水分泥炭)

Usage
-----
    python3 scripts/benchmark_phase2_multi_case_matrix.py
    python3 scripts/benchmark_phase2_multi_case_matrix.py --cases CASE_HTW_WESSELING_2
    python3 scripts/benchmark_phase2_multi_case_matrix.py --unity-scales --max-global-iter 80
"""

from __future__ import annotations

import argparse
from time import perf_counter
from typing import Any

import numpy as np

from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    apply_baseline_kinetics,
    build_phase2_reactor,
    ensure_repo_path,
    phase2_solve_kwargs,
    write_json,
)

ensure_repo_path()

from tests.validation.case_lu_loader import (
    MULTI_CASE_MATRIX_KEYS,
    json_numeric_or_none,
    load_validation_case_flat,
    load_validation_case_node,
    validation_numeric_tolerances,
)

MATRIX_SLOTS: dict[str, dict[str, str]] = {
    "CASE_HTW_WESSELING_1": {
        "slot": "htw_air_steam_baseline",
        "branch": "Air+Steam",
        "focus": "CO/CO₂/H₂ ±15%；轴向 T、ξ≈0.6 峰",
    },
    "CASE_HTW_WESSELING_2": {
        "slot": "htw_o2_steam",
        "branch": "O2+Steam",
        "focus": "绝对 CO/H₂ ~3× LU；CO/(CO+CO₂) 趋势不下降",
    },
    "CASE_VTT_PRESSURIZED_PEAT_6": {
        "slot": "vtt_peat_high_moisture",
        "branch": "Air+Steam",
        "focus": "低床温、高 CO₂",
    },
}

TREND_BANDS: dict[str, tuple[float, float]] = {
    "w2_co_abs_vs_lu": (2.0, 4.0),
    "w2_co_ratio_vs_lu": (0.85, 1.15),
    "peat_co2_abs_vs_lu": (1.2, 2.5),
    "peat_t_exit_vs_lu_K": (-150.0, -20.0),
}


def _co_co2_ratio(dry: dict[str, float]) -> float:
    co = float(dry.get("CO", 0.0))
    co2 = float(dry.get("CO2", 0.0))
    denom = max(co + co2, 1e-12)
    return co / denom


def _ref_dry_outputs(node: dict) -> dict[str, float | None]:
    outs = node.get("outputs", {})
    dry = outs.get("exit_gas_dry_mol_frac", {})
    ref: dict[str, float | None] = {}
    for sp in ("CO", "CO2", "H2", "CH4"):
        ref[sp] = json_numeric_or_none(dry.get(sp))
    ref["T_exit_K"] = json_numeric_or_none(outs.get("exit_temperature_K"))
    ref["co_co2_ratio"] = (
        _co_co2_ratio({k: float(v) for k, v in ref.items() if k in ("CO", "CO2") and v is not None})
        if ref.get("CO") is not None and ref.get("CO2") is not None
        else None
    )
    return ref


def _relative_error(sim: float, ref: float) -> float:
    return abs(sim - ref) / max(abs(ref), 1e-12)


def _run_case(
    case_key: str,
    *,
    max_global_iter: int,
    ladder_level: str,
    unity_scales: bool,
) -> dict[str, Any]:
    node = load_validation_case_node(case_key)
    case = load_validation_case_flat(case_key)
    ref = _ref_dry_outputs(node)
    slot = MATRIX_SLOTS.get(case_key, {"slot": case_key, "branch": "unknown", "focus": ""})

    cfg_overrides: dict[str, object] | None = None
    if unity_scales:
        cfg_overrides = {}

    t0 = perf_counter()
    reactor = build_phase2_reactor(
        case=case,
        ladder_level=ladder_level,
        cfg_overrides=cfg_overrides,
    )
    if unity_scales:
        apply_baseline_kinetics(reactor.config)

    solve_kwargs = phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False)
    out = reactor.solve(**solve_kwargs)
    wall_s = perf_counter() - t0

    dry = out.get("exit_gas_dry") or {}
    sim = {
        "CO": float(dry.get("CO", 0.0)),
        "CO2": float(dry.get("CO2", 0.0)),
        "H2": float(dry.get("H2", 0.0)),
        "CH4": float(dry.get("CH4", 0.0)),
        "T_exit_K": float(out["T_profile"][-1]),
        "co_co2_ratio": _co_co2_ratio(dry),
    }

    tol = validation_numeric_tolerances()
    errors: dict[str, float | None] = {}
    within_hamel: dict[str, bool | None] = {}
    for sp in ("CO", "CO2", "H2", "CH4"):
        ref_v = ref.get(sp)
        if ref_v is None:
            errors[sp] = None
            within_hamel[sp] = None
        else:
            errors[sp] = _relative_error(sim[sp], float(ref_v))
            within_hamel[sp] = errors[sp] <= tol["rtol_CO_CO2_H2" if sp != "CH4" else "rtol_CH4"]

    ref_t = ref.get("T_exit_K")
    if ref_t is None:
        errors["T_exit_K"] = None
        within_hamel["T_exit_K"] = None
    else:
        errors["T_exit_K"] = _relative_error(sim["T_exit_K"], float(ref_t))
        within_hamel["T_exit_K"] = errors["T_exit_K"] <= tol["rtol_T"]

    return {
        "case_key": case_key,
        "slot": slot["slot"],
        "branch": slot["branch"],
        "focus": slot["focus"],
        "primary_agent": case.get("primary_agent"),
        "fuel_type": case.get("fuel_type"),
        "ER": case.get("ER"),
        "P_MPa": float(case["P"]) / 1e6,
        "unity_scales": bool(unity_scales),
        "ladder_level": ladder_level,
        "wall_time_s": float(wall_s),
        "converged_fully": bool(out.get("converged_fully", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", np.nan) or np.nan),
        "nr_outer_iters": int(out.get("nr_outer_iters", 0)),
        "bed0_T_K": float(reactor.cells[0].T),
        "reference": ref,
        "simulated": sim,
        "relative_error": errors,
        "within_hamel_band": within_hamel,
        "cfg_scales": {
            "r2_scale": float(reactor.config.r2_scale),
            "r4_scale": float(reactor.config.r4_scale),
            "r5_scale": float(reactor.config.r5_scale),
            "r8_scale": float(reactor.config.r8_scale),
        },
    }


def _evaluate_trends(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    lu = rows.get("CASE_HTW_WESSELING_1")
    w2 = rows.get("CASE_HTW_WESSELING_2")
    peat = rows.get("CASE_VTT_PRESSURIZED_PEAT_6")
    trends: dict[str, Any] = {}

    def _trend_item(
        name: str,
        value: float,
        *,
        valid: bool,
        note: str = "",
    ) -> dict[str, Any]:
        band = TREND_BANDS[name]
        return {
            "value": float(value),
            "band": list(band),
            "pass": bool(valid and band[0] <= value <= band[1]),
            "valid": bool(valid),
            "note": note,
        }

    if lu and w2:
        lu_co = lu["simulated"]["CO"]
        w2_co = w2["simulated"]["CO"]
        lu_ratio = lu["simulated"]["co_co2_ratio"]
        w2_ratio = w2["simulated"]["co_co2_ratio"]
        valid = bool(lu["converged_fully"] and w2["converged_fully"])
        note = "" if valid else "至少一方未 converged_fully，趋势仅诊断"
        co_mult = w2_co / max(lu_co, 1e-12)
        ratio_mult = w2_ratio / max(lu_ratio, 1e-12)
        trends["w2_co_abs_vs_lu"] = _trend_item("w2_co_abs_vs_lu", co_mult, valid=valid, note=note)
        trends["w2_co_ratio_vs_lu"] = _trend_item("w2_co_ratio_vs_lu", ratio_mult, valid=valid, note=note)

    if lu and peat:
        co2_mult = peat["simulated"]["CO2"] / max(lu["simulated"]["CO2"], 1e-12)
        t_delta = peat["simulated"]["T_exit_K"] - lu["simulated"]["T_exit_K"]
        valid = bool(lu["converged_fully"] and peat["converged_fully"])
        note = "" if valid else "Peat 未 converged_fully，趋势仅诊断"
        trends["peat_co2_abs_vs_lu"] = _trend_item("peat_co2_abs_vs_lu", co2_mult, valid=valid, note=note)
        trends["peat_t_exit_vs_lu_K"] = _trend_item("peat_t_exit_vs_lu_K", t_delta, valid=valid, note=note)

    return trends


def _print_summary(payload: dict[str, Any]) -> None:
    print("=== Phase2 multi-case matrix ===", flush=True)
    for key, row in payload["cases"].items():
        conv = "OK" if row["converged_fully"] else "FAIL"
        sim = row["simulated"]
        ref = row["reference"]
        print(
            f"{key} [{row['branch']}] conv={conv} rms={row['rms_scaled_final']:.4f} "
            f"T={sim['T_exit_K']:.0f}K (ref {ref.get('T_exit_K')}) | "
            f"CO={sim['CO']:.3f} CO2={sim['CO2']:.3f} ratio={sim['co_co2_ratio']:.3f} | "
            f"({row['wall_time_s']:.0f}s)",
            flush=True,
        )
    trends = payload.get("trends", {})
    if trends:
        print("--- trend checks ---", flush=True)
        for name, item in trends.items():
            valid = item.get("valid", True)
            status = "PASS" if item.get("pass") else ("DIAG" if not valid else "FAIL")
            note = item.get("note", "")
            suffix = f" ({note})" if note else ""
            print(f"  {name}: {item['value']:.3f} band={item['band']} -> {status}{suffix}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2 multi-case validation matrix")
    parser.add_argument(
        "--cases",
        nargs="*",
        default=list(MULTI_CASE_MATRIX_KEYS),
        help="validation case keys (default: full §7 matrix)",
    )
    parser.add_argument("--ladder-level", default="p2", choices=["baseline", "p0", "p1", "p2"])
    parser.add_argument(
        "--unity-scales",
        action="store_true",
        help="reset r2/r4/r5/r6/r7/r8_scale to 1.0 (diagnostic unity baseline)",
    )
    add_max_global_iter_arg(parser, default=150)
    add_json_output_arg(parser, "phase2_multi_case_matrix_latest.json")
    args = parser.parse_args()

    rows: dict[str, dict[str, Any]] = {}
    for case_key in args.cases:
        print(f"Running {case_key} ...", flush=True)
        rows[case_key] = _run_case(
            case_key,
            max_global_iter=int(args.max_global_iter),
            ladder_level=str(args.ladder_level),
            unity_scales=bool(args.unity_scales),
        )

    payload = {
        "matrix": "validation_cases_chemistry_map §7",
        "cases_requested": list(args.cases),
        "ladder_level": str(args.ladder_level),
        "max_global_iter": int(args.max_global_iter),
        "unity_scales": bool(args.unity_scales),
        "cases": rows,
        "trends": _evaluate_trends(rows),
    }
    write_json(args.output, payload)
    _print_summary(payload)
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
