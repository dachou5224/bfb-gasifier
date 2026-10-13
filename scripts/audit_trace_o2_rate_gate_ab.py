#!/usr/bin/env python3
"""微量 O₂ 速率门控 A/B（handoff §5.15）：保 CO 且压 bed2 翻转？

对比：ps 基线 / rate_gate bed2+ / rate_gate all / snap bed2+ / hot_restart。

用法：
  python3 scripts/audit_trace_o2_rate_gate_ab.py --max-global-iter 8
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


def _bed2_R(reactor: Any) -> dict[str, float]:
    c = reactor.cells[2]
    c.residuals()
    out: dict[str, float] = {"T": float(c.T)}
    for sp in ("H2O", "H2", "O2", "CO"):
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"R_{sp}"] = float(c.R_gas_d[j] + c.R_gas_b[j])
    j_o2 = int(GAS_SPECIES_INDEX["O2"])
    out["N_O2"] = float(c.N_d[j_o2] + c.N_b[j_o2])
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
    bed2 = _bed2_R(reactor)
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "converged": bool(result.get("converged", False)),
        "bed2": bed2,
        "overrides": dict(overrides),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} split={out['split']:.4f} "
        f"bed2_RH2O={bed2['R_H2O']:+.2f} N_O2={bed2['N_O2']:.2e} wall={wall_s:.0f}s"
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()
    arms = [
        _run("ps_f0.3", {}, args.max_global_iter),
        _run(
            "rate_gate_bed2plus_immediate",
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 2,
            },
            args.max_global_iter,
        ),
        _run(
            "rate_gate_bed2plus_after_yco",
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 2,
                "nr_trace_o2_rate_gate_after_y_co_thesis": True,
                "nr_trace_o2_rate_gate_y_co_threshold_thesis": 0.14,
            },
            args.max_global_iter,
        ),
        _run(
            "hot_restart",
            {
                "nr_outer_hot_restart_after_y_co_thesis": True,
                "nr_outer_hot_restart_y_co_threshold_thesis": 0.14,
            },
            args.max_global_iter,
        ),
        _run(
            "rate_gate_after_yco_plus_hot_restart",
            {
                "nr_trace_o2_rate_gate_thesis": True,
                "nr_trace_o2_rate_gate_min_bed_index_thesis": 2,
                "nr_trace_o2_rate_gate_after_y_co_thesis": True,
                "nr_trace_o2_rate_gate_y_co_threshold_thesis": 0.14,
                "nr_outer_hot_restart_after_y_co_thesis": True,
                "nr_outer_hot_restart_y_co_threshold_thesis": 0.14,
            },
            args.max_global_iter,
        ),
    ]
    out = data_path("phase2_trace_o2_rate_gate_ab_o8.json")
    write_json(
        out,
        {
            "note": "trace O2 rate-gate A/B vs snap/hot_restart; isomorphic ps_f0.3",
            "max_global_iter": int(args.max_global_iter),
            "arms": arms,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
