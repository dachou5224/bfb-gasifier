#!/usr/bin/env python3
"""全量模拟 merit/RMS 追踪：Init → Outer(refresh/pre-inner) → Inner NR。

用法：
  python scripts/audit_full_sim_merit_trace.py
  python scripts/audit_full_sim_merit_trace.py --max-global-iter 40 --output data/phase2_full_sim_merit_trace.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

import importlib.util

_breakdown_path = Path(__file__).resolve().parent / "audit_init_residual_breakdown.py"
_breakdown_spec = importlib.util.spec_from_file_location(
    "audit_init_residual_breakdown",
    _breakdown_path,
)
_breakdown_mod = importlib.util.module_from_spec(_breakdown_spec)
assert _breakdown_spec.loader is not None
_breakdown_spec.loader.exec_module(_breakdown_mod)
_collect_init_residual_breakdown = _breakdown_mod._collect_init_residual_breakdown


def _finite(v: float) -> bool:
    return math.isfinite(float(v))


def _merit_driver(*, rms: float, split: float, energy: float) -> str:
    merit = max(rms, split, energy)
    eps = 1.0e-12
    if abs(merit - split) <= eps:
        return "split"
    if abs(merit - energy) <= eps:
        return "energy"
    if abs(merit - rms) <= eps:
        return "rms"
    return "mixed"


def _pick(d: dict[str, Any], key: str, default: float = 0.0) -> float:
    try:
        return float(d.get(key, default))
    except (TypeError, ValueError):
        return float(default)


def _summarize_init(reactor: Any) -> dict[str, Any]:
    from src.core.connectivity import _measure_nr_init_scaled_metrics

    cells = reactor._solver_cells_for_nr()
    metrics = _measure_nr_init_scaled_metrics(
        cells,
        reactor.config,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
    )
    breakdown = _collect_init_residual_breakdown(reactor)
    abgleich = getattr(reactor, "_vorab_a_tier_init_abgleich_diag", None)
    joint = getattr(reactor, "_bottom_joint_preprojection_run_diag", None)
    return {
        "scaled_metrics": metrics,
        "residual_breakdown": {
            k: breakdown[k]
            for k in (
                "rms_scaled",
                "merit_scaled",
                "gas_rms",
                "solid_rms",
                "energy_rms",
                "gas_phase_split_rms",
                "component_max_rms",
                "top_families_by_rms",
                "top_equations",
            )
            if k in breakdown
        },
        "init_abgleich_diag": abgleich,
        "bottom_joint_preprojection_diag": joint,
        "merit_driver": _merit_driver(
            rms=_pick(metrics, "rms_scaled"),
            split=_pick(metrics, "gas_phase_split_rms"),
            energy=_pick(metrics, "energy_rms"),
        ),
    }


def _inner_peak(history: list[dict[str, Any]]) -> dict[str, Any]:
    if not history:
        return {"enabled": False}
    merit_vals = [_pick(h, "merit") for h in history]
    idx = int(max(range(len(merit_vals)), key=lambda i: merit_vals[i]))
    row = history[idx]
    return {
        "enabled": True,
        "peak_iter": int(row.get("iter", idx)),
        "peak_merit": _pick(row, "merit"),
        "peak_split_rms": _pick(row, "split_rms"),
        "peak_energy_rms": _pick(row, "energy_rms"),
        "peak_rms_scaled": _pick(row, "rms_scaled"),
        "t_max_step_K_at_peak": _pick(row, "t_max_step_K", float("nan")),
        "driver_at_peak": _merit_driver(
            rms=_pick(row, "rms_scaled"),
            split=_pick(row, "split_rms"),
            energy=_pick(row, "energy_rms"),
        ),
    }


def _build_outer_trace(outer_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    trace: list[dict[str, Any]] = []
    prev_merit = float("nan")
    prev_split = float("nan")
    for row in outer_rows:
        oi = int(row.get("outer_iter", len(trace) + 1))
        rms = _pick(row, "inner_rms_scaled_final")
        split = _pick(row, "gas_phase_split_rms")
        energy = _pick(row, "outer_energy_rms_final")
        if not _finite(energy):
            energy = _pick(row, "inner_rms_scaled_final")
        composite = _pick(row, "outer_merit_composite")
        if not _finite(composite):
            composite = max(rms, split, energy)
        driver = _merit_driver(rms=rms, split=split, energy=energy)
        pre_split_before = _pick(row, "preinner_split_rms_before")
        pre_split_after = _pick(row, "preinner_split_rms_after")
        pre_merit_before = _pick(row, "preinner_merit_before")
        pre_merit_after = _pick(row, "preinner_merit_after")
        dT_refresh = _pick(row, "dT_refresh")
        inner_hist = list(row.get("inner_group_history") or [])
        inner_peak = _inner_peak(inner_hist)
        reasons: list[str] = []
        if _finite(prev_merit) and composite > prev_merit + 0.01:
            if _finite(prev_split) and split > prev_split + 0.02:
                reasons.append("outer_end_split_up")
            if pre_split_before > split + 0.02:
                reasons.append("preinner_refresh_split_shock")
            if dT_refresh > 15.0:
                reasons.append("vorab_refresh_dT")
            if bool(row.get("refresh_skipped_for_check1_gate")):
                reasons.append("refresh_skipped_check1_continuation")
            if pre_merit_after > pre_merit_before + 0.005:
                reasons.append("preinner_merit_rejected_or_holdup")
            elif pre_merit_after < pre_merit_before - 0.005:
                reasons.append("preinner_merit_improved")
            if not bool(row.get("preinner_global_gate_pass", True)):
                reasons.append("preinner_global_gate_rollback")
            if inner_peak.get("enabled") and inner_peak.get("driver_at_peak") == "split":
                reasons.append("inner_chased_split_via_T")
            if not reasons:
                reasons.append("composite_merit_up_untagged")
        elif _finite(prev_merit) and composite < prev_merit - 0.01:
            reasons.append("composite_merit_down")
        trace.append(
            {
                "outer_iter": oi,
                "composite_merit": composite,
                "inner_rms": rms,
                "split_rms": split,
                "energy_rms": energy,
                "merit_driver": driver,
                "dT_refresh_K": dT_refresh,
                "dT_outer_K": _pick(row, "dT_outer"),
                "preinner_split_before": pre_split_before,
                "preinner_split_after": pre_split_after,
                "preinner_merit_delta": pre_merit_after - pre_merit_before,
                "preinner_global_gate_pass": bool(row.get("preinner_global_gate_pass", False)),
                "preinner_split_shock": bool(row.get("preinner_split_shock", False)),
                "preinner_bed0_accepted": bool(row.get("preinner_bed0_accepted", False)),
                "inner_n_iter": int(row.get("inner_n_iter", 0)),
                "inner_converged": bool(row.get("inner_converged", False)),
                "inner_peak": inner_peak,
                "delta_composite_vs_prev": (
                    composite - prev_merit if _finite(prev_merit) else float("nan")
                ),
                "likely_reasons": reasons,
            }
        )
        prev_merit = composite
        prev_split = split
    return trace


def _detect_drift_events(trace: list[dict[str, Any]], *, jump_tol: float = 0.015) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    for row in trace:
        delta = _pick(row, "delta_composite_vs_prev")
        if _finite(delta) and delta >= jump_tol:
            events.append(
                {
                    "outer_iter": int(row["outer_iter"]),
                    "delta_composite": delta,
                    "composite_merit": _pick(row, "composite_merit"),
                    "merit_driver": row.get("merit_driver"),
                    "likely_reasons": list(row.get("likely_reasons") or []),
                    "preinner_split_before": _pick(row, "preinner_split_before"),
                    "split_rms_after": _pick(row, "split_rms"),
                }
            )
    return events


def run_trace(*, max_global_iter: int, ladder_level: str) -> dict[str, Any]:
    # Init 诊断：独立 fresh reactor，与 solve 内 init 同配置、不污染求解轨迹
    init_reactor = build_phase2_reactor(ladder_level=ladder_level)
    init_result = run_phase2_init(init_reactor, init_strategy="vorabrechnung")
    init_summary = _summarize_init(init_reactor)

    reactor = build_phase2_reactor(ladder_level=ladder_level)
    reactor.config.nr_inner_record_group_history_thesis = True
    solve_kwargs = phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False)
    result = reactor.solve(**solve_kwargs)

    outer_rows = list(result.get("nr_outer_history") or [])
    trace = _build_outer_trace(outer_rows)
    drift_events = _detect_drift_events(trace)

    solve_summary = {
        "converged": bool(result.get("converged", False)),
        "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
        "nr_merit_early_stop": bool(result.get("nr_merit_early_stop", False)),
        "nr_merit_early_stop_reason": str(result.get("nr_merit_early_stop_reason", "")),
        "nr_inner_budget_total": int(result.get("nr_inner_budget_total", 0)),
        "nr_inner_budget_used": int(result.get("nr_inner_budget_used", 0)),
        "rms_scaled_final": _pick(result, "rms_scaled_final"),
        "rms_scaled_energy_final": _pick(result, "rms_scaled_energy_final"),
        "rms_scaled_gas_phase_split_final": _pick(result, "rms_scaled_gas_phase_split_final"),
    }
    final_driver = _merit_driver(
        rms=solve_summary["rms_scaled_final"],
        split=solve_summary["rms_scaled_gas_phase_split_final"],
        energy=solve_summary["rms_scaled_energy_final"],
    )

    return {
        "max_global_iter": int(max_global_iter),
        "ladder_level": ladder_level,
        "init_timing_s": {
            "nr_init_s_total": float(getattr(init_result, "nr_init_s_total", float("nan"))),
            "nr_vorabrechnung_s": float(getattr(init_result, "nr_vorabrechnung_s", float("nan"))),
        },
        "init": init_summary,
        "solve_summary": solve_summary,
        "final_merit_driver": final_driver,
        "outer_trace": trace,
        "drift_events": drift_events,
        "outer_rows": outer_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=40)
    parser.add_argument("--ladder-level", default="p2")
    parser.add_argument("--output", type=Path, default=data_path("phase2_full_sim_merit_trace.json"))
    args = parser.parse_args()

    payload = run_trace(
        max_global_iter=int(args.max_global_iter),
        ladder_level=str(args.ladder_level),
    )
    write_json(args.output, payload)

    init_m = payload["init"]["scaled_metrics"]
    drift = payload["drift_events"]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "init_merit": init_m.get("merit_scaled"),
                "init_split": init_m.get("gas_phase_split_rms"),
                "init_energy": init_m.get("energy_rms"),
                "init_driver": payload["init"]["merit_driver"],
                "final_merit": payload["solve_summary"]["rms_scaled_final"],
                "final_split": payload["solve_summary"]["rms_scaled_gas_phase_split_final"],
                "final_driver": payload["final_merit_driver"],
                "outers": payload["solve_summary"]["nr_outer_iters"],
                "early_stop": payload["solve_summary"]["nr_merit_early_stop"],
                "inner_budget": f"{payload['solve_summary']['nr_inner_budget_used']}/{payload['solve_summary']['nr_inner_budget_total']}",
                "drift_events": drift[:8],
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
