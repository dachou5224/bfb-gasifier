#!/usr/bin/env python3
"""延迟 rate-gate 后压 rms：加长 outer / 组合热重启 A/B（handoff §5.16）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_rate_gate_rms_stack_ab.py
"""

from __future__ import annotations

import argparse
from typing import Any

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.species import GAS_SPECIES_INDEX

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}

PS = {
    "nr_outer_preinner_upper_co_local2x2_thesis": True,
    "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
    "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
    "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.3,
    "nr_outer_preinner_upper_dense_absorb_mode_thesis": "phase_share",
    "nr_outer_preinner_upper_dense_absorb_post_gate_thesis": True,
}

GATE_AFTER = {
    "nr_trace_o2_rate_gate_thesis": True,
    "nr_trace_o2_rate_gate_min_bed_index_thesis": 2,
    "nr_trace_o2_rate_gate_after_y_co_thesis": True,
    "nr_trace_o2_rate_gate_y_co_threshold_thesis": 0.14,
}

HOT = {
    "nr_outer_hot_restart_after_y_co_thesis": True,
    "nr_outer_hot_restart_y_co_threshold_thesis": 0.14,
}


def _bed2_R(reactor: Any) -> dict[str, float]:
    c = reactor.cells[2]
    c.residuals()
    out: dict[str, float] = {"T": float(c.T)}
    for sp in ("H2O", "H2", "O2", "CO"):
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"R_{sp}"] = float(c.R_gas_d[j] + c.R_gas_b[j])
    return out


def _run(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    hist = list(result.get("nr_outer_history") or [])
    gate_events = [
        {
            "outer": h.get("outer_iter"),
            "armed": h.get("trace_o2_rate_gate_armed"),
            "newly": h.get("trace_o2_rate_gate_newly_armed"),
            "y_co": h.get("trace_o2_rate_gate_y_co_top") or h.get("y_co_top"),
        }
        for h in hist
        if h.get("trace_o2_rate_gate_newly_armed") or h.get("trace_o2_rate_gate_armed")
    ]
    bed2 = _bed2_R(reactor)
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "converged": bool(result.get("converged", False)),
        "n_acc": sum(
            1 for x in (result.get("nr_accepted_lambda_history") or []) if x is not None
        ),
        "bed2": bed2,
        "gate_events": gate_events[:8],
        "hot_restart": dict(result.get("nr_outer_hot_restart") or {}),
        "overrides": dict(overrides),
        "max_global_iter": int(max_iter),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} split={out['split']:.4f} "
        f"E={out['energy']:.4f} bed2_RH2O={bed2['R_H2O']:+.2f} "
        f"n_acc={out['n_acc']} wall={wall_s:.0f}s",
        flush=True,
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="只跑 o8 四臂")
    args = ap.parse_args()
    arms = [
        _run("ps_o8", {}, 8),
        _run("gate_after_yco_o8", dict(GATE_AFTER), 8),
        _run("gate_after_yco_o16", dict(GATE_AFTER), 16),
        _run("hot_restart_o8", dict(HOT), 8),
        _run("gate_after_yco_plus_hot_o8", {**GATE_AFTER, **HOT}, 8),
    ]
    if not args.quick:
        arms.append(_run("gate_after_yco_plus_hot_o12", {**GATE_AFTER, **HOT}, 12))
    out = data_path("phase2_rate_gate_rms_stack_ab.json")
    write_json(
        out,
        {
            "note": "delay rate-gate + longer outer / hot_restart for rms; §5.16",
            "arms": arms,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
