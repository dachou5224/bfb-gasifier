#!/usr/bin/env python3
"""local2(bed6) → bed2 H2O/H2 翻转：耦合路径解剖（handoff §5.13）。

假设：
  H1 回流/BC 直耦：仅 oneshot 改 bed6 CO 并 refresh BC 后，bed2 R_H2O 立即翻转
  H2 仅 NR 轨迹：oneshot 后 bed2 不变，多 outer 后才翻转
  H3 轴向下传：改 bed6 不影响 bed2 入口，只改下游/回流

用法：
  python3 scripts/audit_bed6_local2_to_bed2_coupling.py --max-global-iter 4
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

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

from src.core.connectivity.preprojection import run_outer_preinner_upper_co_local2x2
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


def _bed_gas_snapshot(cell: Any, species: tuple[str, ...] = ("H2O", "H2", "O2", "CO", "CO2", "N2")) -> dict[str, Any]:
    cell.residuals()
    out: dict[str, Any] = {
        "T": float(getattr(cell, "T", float("nan"))),
    }
    for sp in species:
        j = int(GAS_SPECIES_INDEX[sp])
        Nd = float(cell.N_d[j])
        Nb = float(cell.N_b[j])
        tot = Nd + Nb
        out[sp] = {
            "Nd": Nd,
            "Nb": Nb,
            "tot": tot,
            "Rd": float(cell.R_gas_d[j]),
            "Rb": float(cell.R_gas_b[j]),
            "Rnet": float(cell.R_gas_d[j] + cell.R_gas_b[j]),
            "N_in": float(cell.N_d_in[j] + cell.N_b_in[j]),
            "N_zu": float(cell.N_zu_d[j] + cell.N_zu_b[j]),
            "N_rez": float(cell.N_rez_d[j] + cell.N_rez_b[j]),
            "N_ex": float(cell.N_ex[j]),
        }
    return out


def _snap_pair(reactor: Any) -> dict[str, Any]:
    return {
        "bed2": _bed_gas_snapshot(reactor.cells[2]),
        "bed6": _bed_gas_snapshot(reactor.cells[6]),
        "bed1": _bed_gas_snapshot(reactor.cells[1]),
        "bed0": _bed_gas_snapshot(reactor.cells[0]),
    }


def _delta_R(a: dict[str, Any], b: dict[str, Any], bed: str, sp: str) -> float:
    return float(b[bed][sp]["Rnet"] - a[bed][sp]["Rnet"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=4)
    args = ap.parse_args()

    # 1) Baseline plateau (no local2) — short
    r0 = build_phase2_reactor(ladder_level="p2", cfg_overrides=dict(ISOMORPH))
    run_phase2_init(r0)
    res0, w0 = timed_phase2_solve(r0, **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter)))
    snap0 = _snap_pair(r0)
    print(
        f"[baseline] CO={float((res0.get('bed_exit_gas') or {}).get('CO', float('nan'))):.4f} "
        f"bed2 R_H2O={snap0['bed2']['H2O']['Rnet']:+.3f} R_H2={snap0['bed2']['H2']['Rnet']:+.3f}"
    )

    # 2) Same init, oneshot local2 restricted to single beds
    trials = []
    for min_bed, max_bed, label in (
        (6, 6, "oneshot_bed6_only"),
        (7, 9, "oneshot_bed7_9"),
        (6, 9, "oneshot_bed6_9"),
    ):
        r = build_phase2_reactor(ladder_level="p2", cfg_overrides=dict(ISOMORPH))
        run_phase2_init(r)
        # march to same short plateau without local2
        timed_phase2_solve(r, **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter)))
        before = _snap_pair(r)
        cfg = r.config
        cfg.nr_outer_preinner_upper_co_local2x2_thesis = True
        cfg.nr_outer_preinner_upper_co_local2x2_min_bed_index_thesis = int(min_bed)
        cfg.nr_outer_preinner_upper_co_local2x2_fraction_thesis = 1.0
        # Temporarily mask beds outside [min_bed, max_bed] by freezing their N after local2
        # Easier: patch dof by setting min_bed and zeroing other beds' CO step via restore
        apply = r._apply_all_bc_for_nr
        # snapshot all bed N
        snaps = [
            (i, np.asarray(c.N_d).copy(), np.asarray(c.N_b).copy())
            for i, c in enumerate(r.cells)
            if str(getattr(c, "cell_type", "bed")) == "bed"
        ]
        diag = run_outer_preinner_upper_co_local2x2(r.cells, cfg, apply_bc_fn=apply)
        # restore beds outside range
        j = int(GAS_SPECIES_INDEX["CO"])
        for i, nd, nb in snaps:
            if i < min_bed or i > max_bed:
                r.cells[i].N_d[:] = nd
                r.cells[i].N_b[:] = nb
                r.cells[i]._thermo_cache_valid = False
                r.cells[i]._hydro_cache_valid = False
        apply()
        after = _snap_pair(r)
        trial = {
            "label": label,
            "min_bed": min_bed,
            "max_bed": max_bed,
            "local2_accepted": bool(diag.get("accepted")),
            "beds_touched": list(diag.get("beds_touched") or []),
            "before": before,
            "after": after,
            "dR_bed2_H2O": _delta_R(before, after, "bed2", "H2O"),
            "dR_bed2_H2": _delta_R(before, after, "bed2", "H2"),
            "dR_bed2_O2": _delta_R(before, after, "bed2", "O2"),
            "dN_bed6_CO": float(after["bed6"]["CO"]["tot"] - before["bed6"]["CO"]["tot"]),
            "dN_bed2_CO": float(after["bed2"]["CO"]["tot"] - before["bed2"]["CO"]["tot"]),
            "d_rez_bed0_CO": float(after["bed0"]["CO"]["N_rez"] - before["bed0"]["CO"]["N_rez"]),
            "d_rez_bed2_H2O": float(after["bed2"]["H2O"]["N_rez"] - before["bed2"]["H2O"]["N_rez"]),
            "d_in_bed2_H2O": float(after["bed2"]["H2O"]["N_in"] - before["bed2"]["H2O"]["N_in"]),
            "dT_bed2": float(after["bed2"]["T"] - before["bed2"]["T"]),
        }
        trials.append(trial)
        print(
            f"[{label}] acc={trial['local2_accepted']} dN6_CO={trial['dN_bed6_CO']:+.3f} "
            f"dR2_H2O={trial['dR_bed2_H2O']:+.3f} dR2_H2={trial['dR_bed2_H2']:+.3f} "
            f"d_rez0_CO={trial['d_rez_bed0_CO']:+.3f} d_in2_H2O={trial['d_in_bed2_H2O']:+.3f} "
            f"dT2={trial['dT_bed2']:+.2f}"
        )

    # 3) Full local2 stack few outers — confirm flip appears in NR path
    rL = build_phase2_reactor(
        ladder_level="p2",
        cfg_overrides={**ISOMORPH, "nr_outer_preinner_upper_co_local2x2_thesis": True},
    )
    run_phase2_init(rL)
    resL, wL = timed_phase2_solve(rL, **phase2_solve_kwargs(max_global_iter=8))
    snapL = _snap_pair(rL)
    print(
        f"[local2_o8] CO={float((resL.get('bed_exit_gas') or {}).get('CO', float('nan'))):.4f} "
        f"bed2 R_H2O={snapL['bed2']['H2O']['Rnet']:+.3f} R_H2={snapL['bed2']['H2']['Rnet']:+.3f}"
    )

    # Verdict
    h1 = any(abs(t["dR_bed2_H2O"]) > 5.0 for t in trials if t["label"] == "oneshot_bed6_only")
    h1_79 = any(abs(t["dR_bed2_H2O"]) > 5.0 for t in trials if t["label"] == "oneshot_bed7_9")
    verdict = {
        "H1_bed6_oneshot_flips_bed2": bool(h1),
        "H1b_bed7_9_oneshot_flips_bed2": bool(h1_79),
        "H2_needs_NR_if_oneshot_small": bool(
            (not h1)
            and abs(float(snapL["bed2"]["H2O"]["Rnet"] - snap0["bed2"]["H2O"]["Rnet"])) > 5.0
        ),
    }
    print("verdict", verdict)

    out = data_path("phase2_bed6_local2_to_bed2_coupling_o4.json")
    write_json(
        out,
        {
            "note": "bed6 local2 → bed2 H2O/H2 coupling; oneshot vs NR",
            "baseline": {
                "exit_CO": float((res0.get("bed_exit_gas") or {}).get("CO", float("nan"))),
                "snap": snap0,
                "wall_s": float(w0),
            },
            "local2_o8": {
                "exit_CO": float((resL.get("bed_exit_gas") or {}).get("CO", float("nan"))),
                "snap": snapL,
                "wall_s": float(wL),
            },
            "oneshot_trials": trials,
            "verdict": verdict,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
