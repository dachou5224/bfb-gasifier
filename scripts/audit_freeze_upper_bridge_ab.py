#!/usr/bin/env python3
"""同一次 solve 内冻结 upper CO 桥 A/B（handoff §5.13）。

对比 y_CO 阈值冻结 vs local2 接受次数冻结；记录 bed2 R_H2O/H2 与冻结时点。

用法：
  python3 scripts/audit_freeze_upper_bridge_ab.py --max-global-iter 8
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
    out: dict[str, float] = {}
    for sp in ("H2O", "H2", "O2", "CO"):
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"R_{sp}"] = float(c.R_gas_d[j] + c.R_gas_b[j])
    out["T"] = float(c.T)
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
    freeze_events = [
        {
            "outer": h.get("outer_iter"),
            "frozen": h.get("upper_co_bridge_frozen"),
            "y_co_top": h.get("y_co_top"),
            "local2": h.get("upper_co_local2x2_accepted"),
            "absorb": h.get("upper_dense_absorb_accepted"),
        }
        for h in hist
        if h.get("upper_co_bridge_frozen")
        or h.get("upper_co_local2x2_enabled")
        or h.get("upper_dense_absorb_enabled")
    ]
    first_freeze = next((e for e in freeze_events if e.get("frozen")), None)
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
        "local2_n_accepted": sum(
            1 for h in hist if h.get("upper_co_local2x2_accepted")
        ),
        "local2_n_enabled": sum(1 for h in hist if h.get("upper_co_local2x2_enabled")),
        "first_freeze_outer": None if first_freeze is None else first_freeze.get("outer"),
        "freeze_events": freeze_events,
        "bed2": bed2,
        "overrides": dict(overrides),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} split={out['split']:.4f} "
        f"bed2_RH2O={bed2['R_H2O']:+.2f} freeze@{out['first_freeze_outer']} "
        f"local2={out['local2_n_accepted']}/{out['local2_n_enabled']} wall={wall_s:.0f}s"
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()
    arms = [
        _run("ps_f0.3", {}, args.max_global_iter),
        _run(
            "freeze_y0.14",
            {
                "nr_outer_freeze_upper_co_bridge_after_y_co_thesis": True,
                "nr_outer_freeze_upper_co_bridge_y_co_threshold_thesis": 0.14,
            },
            args.max_global_iter,
        ),
        _run(
            "freeze_y0.10",
            {
                "nr_outer_freeze_upper_co_bridge_after_y_co_thesis": True,
                "nr_outer_freeze_upper_co_bridge_y_co_threshold_thesis": 0.10,
            },
            args.max_global_iter,
        ),
        _run(
            "freeze_after_local2_1",
            {"nr_outer_freeze_upper_co_bridge_after_local2_accepts_thesis": 1},
            args.max_global_iter,
        ),
        _run(
            "freeze_after_local2_2",
            {"nr_outer_freeze_upper_co_bridge_after_local2_accepts_thesis": 2},
            args.max_global_iter,
        ),
    ]
    out = data_path("phase2_freeze_upper_bridge_ab_o8.json")
    write_json(
        out,
        {
            "note": "in-solve freeze upper CO bridge A/B; bed2 R + freeze timeline",
            "max_global_iter": int(args.max_global_iter),
            "arms": arms,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
