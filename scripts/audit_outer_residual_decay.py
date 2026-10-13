#!/usr/bin/env python3
"""Outer/inner 残差衰减审计：识别平台项与 pre-inner 有效性。

用法：
  python scripts/audit_outer_residual_decay.py --max-global-iter 20
  python scripts/audit_outer_residual_decay.py --input-json data/some_run.json
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

from src.core.preprojection_engine import preproj_cell_energy_residual_abs
from src.core.cell import S_CHAR, S_ASH


def _plateau_slope(values: list[float], *, tail: int = 3) -> float | None:
    if len(values) < tail + 1:
        return None
    tail_vals = values[-tail:]
    if not all(math.isfinite(v) for v in tail_vals):
        return None
    return float(tail_vals[-1] - tail_vals[0])


def _decay_ratio(first: float, last: float) -> float | None:
    if not math.isfinite(first) or not math.isfinite(last) or first <= 1.0e-12:
        return None
    return float(last / first)


def _analyze_inner_history(history: list[dict[str, Any]]) -> dict[str, Any]:
    if not history:
        return {"enabled": False}
    keys = ("merit", "rms_scaled", "gas_rms", "solid_rms", "energy_rms", "split_rms")
    series = {k: [float(row.get(k, float("nan"))) for row in history] for k in keys}
    out: dict[str, Any] = {"enabled": True, "n_iter": len(history), "series": series}
    for key, vals in series.items():
        finite = [v for v in vals if math.isfinite(v)]
        if not finite:
            continue
        out[f"{key}_first"] = finite[0]
        out[f"{key}_last"] = finite[-1]
        out[f"{key}_decay_ratio"] = _decay_ratio(finite[0], finite[-1])
        out[f"{key}_tail_slope"] = _plateau_slope(finite)
    return out


def _analyze_outer_history(outer_history: list[dict[str, Any]]) -> dict[str, Any]:
    if not outer_history:
        return {"enabled": False}

    def _col(name: str) -> list[float]:
        return [float(row.get(name, float("nan"))) for row in outer_history]

    inner_rms = _col("inner_rms_scaled_final")
    split_rms = _col("gas_phase_split_rms")
    preinner_delta_merit = [
        float(row.get("preinner_merit_after", 0.0)) - float(row.get("preinner_merit_before", 0.0))
        for row in outer_history
    ]
    skipped = sum(1 for row in outer_history if bool(row.get("preinner_skipped_all_tiers", False)))
    reconcile = sum(1 for row in outer_history if bool(row.get("preinner_holdup_reconciled", False)))

    diagnosis: list[str] = []
    if _decay_ratio(inner_rms[0], inner_rms[-1]) is not None and _decay_ratio(inner_rms[0], inner_rms[-1]) > 0.9:
        diagnosis.append("outer_inner_rms_plateau: inner NR 末态 RMS 下降不足，可能为结构型 BC/merit 冲突")
    if split_rms and max(split_rms) > 0.05:
        diagnosis.append("split_elevated: phase-split 持续高于 0.05，优先 T2 split 而非 solid LSQ")
    if preinner_delta_merit and sum(abs(d) for d in preinner_delta_merit) / max(len(preinner_delta_merit), 1) < 1e-3:
        diagnosis.append("preinner_low_impact: pre-inner merit 增量均值接近 0，考虑 merit-driven 跳过或换 tier")
    if reconcile == 0 and any(
        float(row.get("preinner_run_holdup_reconcile", 0.0)) for row in outer_history
    ):
        diagnosis.append("holdup_reconcile_noop: T1 已调度但未改状态，检查 flux 阈值或 bed 过滤")
    if any(float(row.get("preinner_holdup_flux_mismatch_before", 0.0)) > 0.15 for row in outer_history):
        diagnosis.append("holdup_flux_mismatch: bed1+ auf_in/out_auf 偏离，优先 T1 reconcile")
    if skipped > 0:
        diagnosis.append(f"preinner_skipped_tiers: {skipped}/{len(outer_history)} outer 轮跳过全部 tier")

    return {
        "enabled": True,
        "n_outer": len(outer_history),
        "inner_rms_first": inner_rms[0],
        "inner_rms_last": inner_rms[-1],
        "inner_rms_decay_ratio": _decay_ratio(inner_rms[0], inner_rms[-1]),
        "split_rms_max": max(split_rms) if split_rms else None,
        "preinner_merit_delta_mean": sum(preinner_delta_merit) / max(len(preinner_delta_merit), 1),
        "preinner_skipped_outer_count": skipped,
        "preinner_holdup_reconcile_count": reconcile,
        "diagnosis": diagnosis,
        "outer_rows": outer_history,
    }


def _bed_energy_hotspots(reactor: Any) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for i, cell in enumerate(reactor.cells):
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        holdup = float(cell.m_solid.sum())
        outflow = float(cell._solid_outflow_rates()[:, (S_CHAR, S_ASH)].sum())
        auf_in = float(cell.m_solid_auf_in[:, (S_CHAR, S_ASH)].sum())
        rows.append(
            {
                "bed": int(i),
                "T_K": float(cell.T),
                "energy_residual_MW": float(preproj_cell_energy_residual_abs(cell) / 1e6),
                "solid_holdup_kg": holdup,
                "solid_out_kg_s": outflow,
                "solid_auf_in_kg_s": auf_in,
                "auf_in_over_outflow": auf_in / max(outflow, 1e-12),
            }
        )
    rows.sort(key=lambda r: abs(float(r["energy_residual_MW"])), reverse=True)
    return rows


def run_audit(
    *,
    max_global_iter: int,
    ladder_level: str,
    record_inner: bool,
    cfg_overrides: dict[str, object] | None = None,
) -> dict[str, Any]:
    reactor = build_phase2_reactor(ladder_level=ladder_level, cfg_overrides=cfg_overrides)
    reactor.config.nr_inner_record_group_history_thesis = bool(record_inner)
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    from src.core.connectivity import _measure_nr_init_scaled_metrics

    init_metrics = _measure_nr_init_scaled_metrics(
        reactor._solver_cells_for_nr(),
        reactor.config,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
    )

    solve_kwargs = phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False)
    result = reactor.solve(**solve_kwargs)
    outer_history = list(result.get("nr_outer_history") or [])

    last_inner: list[dict[str, Any]] = []
    if outer_history:
        last_inner = list(outer_history[-1].get("inner_group_history") or [])

    payload: dict[str, Any] = {
        "max_global_iter": int(max_global_iter),
        "ladder_level": ladder_level,
        "record_inner_group_history": bool(record_inner),
        "init_metrics": init_metrics,
        "solve_summary": {
            "converged": bool(result.get("converged", False)),
            "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
            "nr_merit_early_stop": bool(result.get("nr_merit_early_stop", False)),
            "nr_merit_early_stop_reason": str(result.get("nr_merit_early_stop_reason", "")),
            "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
            "rms_scaled_energy_final": float(result.get("rms_scaled_energy_final", float("nan"))),
            "rms_scaled_gas_phase_split_final": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        },
        "outer_decay": _analyze_outer_history(outer_history),
        "last_outer_inner_decay": _analyze_inner_history(last_inner),
        "bed_energy_hotspots_final": _bed_energy_hotspots(reactor),
    }
    if cfg_overrides:
        payload["cfg_overrides"] = dict(cfg_overrides)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=20)
    parser.add_argument("--ladder-level", default="p2")
    parser.add_argument("--record-inner", action="store_true", help="开启 inner group history")
    parser.add_argument(
        "--split-t-cap-k",
        type=float,
        default=None,
        help="覆盖 nr_inner_split_dominant_t_step_cap_K_thesis（0=关闭 cap）",
    )
    parser.add_argument("--bed1-split", action="store_true", help="启用 bed1 T2 phase-split")
    parser.add_argument(
        "--split-plateau-improvement",
        type=float,
        default=None,
        help="覆盖 nr_outer_preinner_split_plateau_bed0_improvement_thesis",
    )
    parser.add_argument(
        "--split-tier-min-rms",
        type=float,
        default=None,
        help="覆盖 nr_outer_preinner_split_tier_min_rms_thesis",
    )
    parser.add_argument(
        "--split-max-bed-index",
        type=int,
        default=None,
        help="覆盖 nr_outer_preinner_phase_split_max_bed_index_thesis（2=含 selective bed2）",
    )
    parser.add_argument(
        "--no-merit-early-stop",
        action="store_true",
        help="关闭 nr_outer_merit_early_stop_thesis",
    )
    parser.add_argument("--input-json", type=Path, default=None, help="仅分析已有 JSON（含 outer_rows）")
    parser.add_argument("--output", type=Path, default=data_path("phase2_outer_residual_decay_audit.json"))
    args = parser.parse_args()

    if args.input_json is not None:
        raw = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
        outer_rows = raw.get("outer_decay", {}).get("outer_rows") or raw.get("nr_outer_history") or raw.get("results")
        if isinstance(outer_rows, list) and outer_rows and "inner_rms_scaled_final" in outer_rows[0]:
            payload = {"outer_decay": _analyze_outer_history(outer_rows), "source": str(args.input_json)}
        else:
            print("无法从 JSON 解析 outer_history", file=sys.stderr)
            return 1
    else:
        cfg_overrides: dict[str, object] | None = None
        if args.split_t_cap_k is not None:
            cfg_overrides = {"nr_inner_split_dominant_t_step_cap_K_thesis": float(args.split_t_cap_k)}
        if args.bed1_split:
            cfg_overrides = dict(cfg_overrides or {})
            cfg_overrides["nr_outer_preinner_bed1_phase_split_thesis"] = True
            cfg_overrides["nr_outer_preinner_phase_split_max_bed_index_thesis"] = 1
        if args.split_plateau_improvement is not None:
            cfg_overrides = dict(cfg_overrides or {})
            cfg_overrides["nr_outer_preinner_split_plateau_bed0_improvement_thesis"] = float(
                args.split_plateau_improvement
            )
        if args.split_tier_min_rms is not None:
            cfg_overrides = dict(cfg_overrides or {})
            cfg_overrides["nr_outer_preinner_split_tier_min_rms_thesis"] = float(args.split_tier_min_rms)
        if args.split_max_bed_index is not None:
            cfg_overrides = dict(cfg_overrides or {})
            cfg_overrides["nr_outer_preinner_phase_split_max_bed_index_thesis"] = int(
                args.split_max_bed_index
            )
        if args.no_merit_early_stop:
            cfg_overrides = dict(cfg_overrides or {})
            cfg_overrides["nr_outer_merit_early_stop_thesis"] = False
        payload = run_audit(
            max_global_iter=int(args.max_global_iter),
            ladder_level=str(args.ladder_level),
            record_inner=bool(args.record_inner),
            cfg_overrides=cfg_overrides,
        )

    write_json(args.output, payload)
    print(json.dumps(
        {
            "output": str(args.output),
            "diagnosis": payload.get("outer_decay", {}).get("diagnosis", []),
            "solve_summary": payload.get("solve_summary"),
            "hotspots": payload.get("bed_energy_hotspots_final", [])[:3],
        },
        indent=2,
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
