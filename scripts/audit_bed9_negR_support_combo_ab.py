#!/usr/bin/env python3
"""负 R support vs 床顶 solid LSQ + 延迟焓透传组合 A/B（handoff §5.22）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed9_negR_support_combo_ab.py
"""

from __future__ import annotations

from typing import Any

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from audit_bed9_holdup_closure_ab import _chain_table, _metrics
from audit_rate_gate_rms_stack_ab import GATE_AFTER, HOT, ISOMORPH, PS

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

STACK = {**GATE_AFTER, **HOT}

BED_TOP_LSQ_AFTER = {
    "nr_outer_preinner_bed_top_solid_holdup_thesis": True,
    "nr_outer_preinner_bed_top_solid_holdup_after_y_co_thesis": True,
    "nr_outer_preinner_bed_top_solid_holdup_y_co_threshold_thesis": 0.14,
}

ENERGY_PT_AFTER = {
    "nr_bed_top_solid_energy_passthrough_thesis": True,
    "nr_bed_top_solid_energy_passthrough_after_y_co_thesis": True,
    "nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis": 0.14,
}

NEG_R = {"nr_holdup_support_include_negative_R_thesis": True}


def _run(label: str, overrides: dict[str, object], max_iter: int = 8) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(STACK)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    met = _metrics(reactor, result)
    chain = _chain_table(reactor)
    print(
        f"[{label}] CO={met['exit_CO']:.4f} rms={met['rms']:.4f} "
        f"E={met['groups']['energy_rms']:.4f} "
        f"Fhat9={met['Fhat_T_bed9']:+.3f} res_char9={met['bed9']['res_char']:+.3e} "
        f"over_sup={met['bed9']['over_support_char']:+.3e} "
        f"over_m*={met['bed9']['over_mstar_char']:+.3e} wall={wall_s:.0f}s",
        flush=True,
    )
    print(
        f"  top3: "
        + ", ".join(
            f"c{t['cell']}:{t.get('species') or t.get('kind')}={t['F_hat']:+.3f}"
            for t in met["topk"][:3]
        ),
        flush=True,
    )
    return {
        "label": label,
        "wall_s": float(wall_s),
        "metrics": met,
        "chain_tail": chain[7:],
        "overrides": dict(overrides),
        "passthrough_armed": bool(
            getattr(reactor.cells[9], "_solid_energy_passthrough", False)
        ),
    }


def main() -> None:
    arms = [
        _run("gate_hot_base", {}),
        _run("negR_support", NEG_R),
        _run("bed_top_lsq_after_y_co", BED_TOP_LSQ_AFTER),
        _run("energy_pt_after_y_co", ENERGY_PT_AFTER),
        _run("lsq_plus_energy_pt_after_y_co", {**BED_TOP_LSQ_AFTER, **ENERGY_PT_AFTER}),
        _run(
            "negR_plus_lsq_after_y_co",
            {**NEG_R, **BED_TOP_LSQ_AFTER},
        ),
        _run(
            "full_combo_after_y_co",
            {**NEG_R, **BED_TOP_LSQ_AFTER, **ENERGY_PT_AFTER},
        ),
    ]
    out = data_path("phase2_bed9_negR_support_combo_ab_o8.json")
    write_json(
        out,
        {
            "note": "neg-R holdup support + bed-top LSQ/energy passthrough combos; §5.22",
            "arms": arms,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
