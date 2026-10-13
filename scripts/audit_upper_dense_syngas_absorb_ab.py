#!/usr/bin/env python3
"""上段密相 CO/H2 总量吸收 A/B（handoff §5.5）。

对照臂：
1. baseline：同构栈（bed1 preinner，max_bed=默认 0）
2. absorb：``nr_outer_preinner_upper_dense_syngas_absorb_thesis=True``
3. maxbed9：现有冻总量 LSQ 扩到 bed9（对照：不能关 res_sum）

用法：
  python3 scripts/audit_upper_dense_syngas_absorb_ab.py
  python3 scripts/audit_upper_dense_syngas_absorb_ab.py --max-global-iter 8
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


def _summarize(result: dict[str, Any], reactor: Any) -> dict[str, Any]:
    exit_gas = dict(result.get("bed_exit_gas") or result.get("exit_gas") or {})
    counts = dict(result.get("nr_counts") or {})
    hist = list(result.get("nr_outer_history") or [])
    absorb_hits = []
    for h in hist:
        if h.get("upper_dense_absorb_accepted") or h.get("upper_dense_absorb_enabled"):
            absorb_hits.append(
                {
                    "outer": h.get("outer_iter"),
                    "accepted": h.get("upper_dense_absorb_accepted"),
                    "beds": h.get("upper_dense_absorb_beds"),
                    "added_CO": h.get("upper_dense_absorb_added_CO_mol_s"),
                    "added_H2": h.get("upper_dense_absorb_added_H2_mol_s"),
                }
            )
    # 末态上段 CO res_sum
    j = int(GAS_SPECIES_INDEX["CO"])
    upper_co = []
    for i, cell in enumerate(reactor.cells):
        if i < 6 or str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        cell.residuals()
        in_d = float(cell.N_zu_d[j] + cell.N_rez_d[j] + cell.N_d_in[j])
        in_b = float(cell.N_zu_b[j] + cell.N_rez_b[j] + cell.N_b_in[j])
        res_sum = in_d + in_b + float(cell.R_gas_d[j] + cell.R_gas_b[j]) - float(
            cell.N_d[j] + cell.N_b[j]
        )
        res_split = (
            in_d
            + float(cell.R_gas_d[j])
            - float(cell.N_d[j])
            + float(cell.N_ex[j])
        ) - (
            in_b
            + float(cell.R_gas_b[j])
            - float(cell.N_b[j])
            - float(cell.N_ex[j])
        )
        upper_co.append(
            {
                "cell": i,
                "res_sum": res_sum,
                "res_split": res_split,
                "N_d": float(cell.N_d[j]),
                "N_b": float(cell.N_b[j]),
                "R_d": float(cell.R_gas_d[j]),
            }
        )
    return {
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "ls_fail": int(counts.get("line_search_failures", -1)),
        "n_acc": sum(1 for x in (result.get("nr_accepted_lambda_history") or []) if x is not None),
        "absorb_hits": absorb_hits,
        "upper_CO_res": upper_co,
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
    out = _summarize(result, reactor)
    out["label"] = label
    out["wall_s"] = float(wall_s)
    out["overrides"] = dict(overrides)
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} "
        f"split={out['split']:.4f} E={out['energy']:.4f} "
        f"ls_fail={out['ls_fail']} n_acc={out['n_acc']} wall={wall_s:.0f}s"
    )
    if out["upper_CO_res"]:
        rs = [r["res_sum"] for r in out["upper_CO_res"]]
        sp = [r["res_split"] for r in out["upper_CO_res"]]
        print(
            f"  upper CO |res_sum| mean={sum(abs(x) for x in rs)/len(rs):.3f} "
            f"|res_split| mean={sum(abs(x) for x in sp)/len(sp):.3f}"
        )
    if out["absorb_hits"]:
        print(f"  absorb_hits={out['absorb_hits'][:4]}...")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()

    arms = [
        ("baseline", {}),
        (
            "absorb",
            {"nr_outer_preinner_upper_dense_syngas_absorb_thesis": True},
        ),
        (
            "absorb_soft_CO",
            {
                "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
                "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
                "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.1,
            },
        ),
        (
            "maxbed9",
            {"nr_outer_preinner_phase_split_max_bed_index_thesis": 9},
        ),
        (
            "absorb+maxbed9",
            {
                "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
                "nr_outer_preinner_phase_split_max_bed_index_thesis": 9,
            },
        ),
    ]
    rows = [_run(label, ov, args.max_global_iter) for label, ov in arms]
    payload = {
        "note": "upper dense CO/H2 absorb A/B; isomorphic o8",
        "max_global_iter": int(args.max_global_iter),
        "arms": rows,
    }
    out = data_path("phase2_upper_dense_syngas_absorb_ab_o8.json")
    write_json(out, payload)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
