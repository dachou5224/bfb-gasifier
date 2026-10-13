#!/usr/bin/env python3
"""上段 CO 局部 2×2 Pre-Inner A/B（handoff §5.7）。

用法：
  python3 scripts/audit_upper_co_local2x2_ab.py --max-global-iter 8
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

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}


def _run(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    counts = dict(result.get("nr_counts") or {})
    hist = list(result.get("nr_outer_history") or [])
    hits = [
        {
            "outer": h.get("outer_iter"),
            "accepted": h.get("upper_co_local2x2_accepted"),
            "beds": h.get("upper_co_local2x2_beds"),
        }
        for h in hist
        if h.get("upper_co_local2x2_enabled")
    ]
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "ls_fail": int(counts.get("line_search_failures", -1)),
        "n_acc": sum(
            1 for x in (result.get("nr_accepted_lambda_history") or []) if x is not None
        ),
        "local2_hits": hits,
        "local2_n_accepted": sum(1 for h in hits if h.get("accepted")),
        "overrides": dict(overrides),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} "
        f"split={out['split']:.4f} E={out['energy']:.4f} "
        f"ls_fail={out['ls_fail']} n_acc={out['n_acc']} "
        f"local2_acc={out['local2_n_accepted']}/{len(hits)} wall={wall_s:.0f}s"
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()
    arms = [
        _run("baseline", {}, args.max_global_iter),
        _run(
            "local2x2_post_gate",
            {"nr_outer_preinner_upper_co_local2x2_thesis": True},
            args.max_global_iter,
        ),
        _run(
            "local2x2_bundle",
            {
                "nr_outer_preinner_upper_co_local2x2_thesis": True,
                "nr_outer_preinner_upper_co_local2x2_post_gate_thesis": False,
            },
            args.max_global_iter,
        ),
        _run(
            "local2x2_post+soft_CO",
            {
                "nr_outer_preinner_upper_co_local2x2_thesis": True,
                "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
                "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
                "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.1,
                "nr_outer_preinner_upper_dense_absorb_skip_local_gate_thesis": True,
            },
            args.max_global_iter,
        ),
    ]
    out = data_path("phase2_upper_co_local2x2_ab_o8.json")
    write_json(
        out,
        {
            "note": "upper CO local 2x2 preinner A/B; isomorphic o8",
            "max_global_iter": int(args.max_global_iter),
            "arms": arms,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
