#!/usr/bin/env python3
"""解剖 Pre-Inner phase-split LSQ 爆炸：哪张床 / 哪个方程把 max_abs 打爆（handoff §5.9）。

在 local2x2 post-gate 栈上，逐 bed 记录 split LSQ 前后的全局 rms/max_abs，
并用 F_hat argmax 定位主导残差（cell / 物种 / 相）。

用法：
  python3 scripts/audit_preinner_split_lsq_detonation.py --max-global-iter 6
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

from src.core import connectivity as conn
from src.core.connectivity import preprojection as pp
from src.core.connectivity.joint_gate import _measure_nr_init_scaled_metrics
from src.core.species import GAS_SPECIES
from src.solvers.global_nr_solver import build_equation_scales, global_residual, pack_reactor
from src.solvers.nr_indexing import gas_mode_for_nr, n_gas_var, n_var

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
    "nr_outer_preinner_upper_co_local2x2_thesis": True,
}


def _locate_max_abs(cells: list[Any], cfg: Any, apply_bc_fn: Any) -> dict[str, Any]:
    if apply_bc_fn is not None:
        apply_bc_fn()
    ref_gas = max(
        float(getattr(cfg, "O2_feed", 0.0) + getattr(cfg, "H2O_feed", 0.0) + getattr(cfg, "N2_feed", 0.0)),
        1.0,
    )
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    x = pack_reactor(cells)
    F = np.asarray(global_residual(x, cells, apply_bc_fn or (lambda: None)), dtype=np.float64)
    scale = np.asarray(build_equation_scales(cells, ref_gas, ref_solid, ref_energy), dtype=np.float64)
    F_hat = F / np.maximum(scale, 1e-30)
    k = int(np.argmax(np.abs(F_hat)))
    # map flat index → cell / var
    offset = 0
    loc: dict[str, Any] = {
        "flat_k": k,
        "F_hat": float(F_hat[k]),
        "F_raw": float(F[k]),
        "scale": float(scale[k]),
        "rms": float(np.sqrt(np.mean(F_hat**2))),
        "max_abs": float(np.max(np.abs(F_hat))),
    }
    for i, cell in enumerate(cells):
        nv = int(n_var(cell))
        if offset <= k < offset + nv:
            local = k - offset
            ng = int(n_gas_var(cell))
            loc["cell"] = i
            loc["cell_type"] = str(getattr(cell, "cell_type", "?"))
            loc["local_idx"] = local
            loc["n_gas"] = ng
            if local < ng and gas_mode_for_nr(cell) == "two_phase":
                half = ng // 2
                phase = "dense" if local < half else "bubble"
                gi = local if local < half else local - half
                # gas unknown order follows main NR gas index
                from src.solvers.nr_indexing import main_nr_gas_idx_for_cell

                gidx = main_nr_gas_idx_for_cell(cell)
                sp = str(GAS_SPECIES[int(gidx[gi])]) if gi < len(gidx) else f"g{gi}"
                loc["phase"] = phase
                loc["species"] = sp
                loc["kind"] = "gas"
            elif local == ng and str(getattr(cell, "cell_type", "bed")) == "bed":
                loc["kind"] = "energy_or_first_solid"
            else:
                loc["kind"] = "solid_or_other"
            break
        offset += nv
    # top-5
    top = np.argsort(-np.abs(F_hat))[:5]
    loc["top5"] = []
    for t in top:
        # reuse mapping lightly
        off = 0
        entry: dict[str, Any] = {"flat_k": int(t), "F_hat": float(F_hat[t])}
        for i, cell in enumerate(cells):
            nv = int(n_var(cell))
            if off <= int(t) < off + nv:
                local = int(t) - off
                ng = int(n_gas_var(cell))
                entry["cell"] = i
                if local < ng and gas_mode_for_nr(cell) == "two_phase":
                    half = ng // 2
                    phase = "dense" if local < half else "bubble"
                    gi = local if local < half else local - half
                    from src.solvers.nr_indexing import main_nr_gas_idx_for_cell

                    gidx = main_nr_gas_idx_for_cell(cell)
                    sp = str(GAS_SPECIES[int(gidx[gi])]) if gi < len(gidx) else f"g{gi}"
                    entry["phase"] = phase
                    entry["species"] = sp
                else:
                    entry["kind"] = "non_gas"
                    entry["local_idx"] = local
                break
            off += nv
        loc["top5"].append(entry)
    return loc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=6)
    args = ap.parse_args()

    events: list[dict[str, Any]] = []
    bridge_n = {"n": 0}
    _orig_split = pp.run_outer_preinner_gas_phase_split_preprojection
    _orig_bed = pp.preproject_bed_major_gas_phase_split_to_exchange_closure
    _orig_bridge = pp.run_outer_preinner_start_value_bridge

    ctx: dict[str, Any] = {"apply": None, "cfg": None, "cells": None, "bridge": 0}

    def w_bridge(cells, cfg, **kwargs):
        bridge_n["n"] += 1
        ctx["cells"] = cells
        ctx["cfg"] = cfg
        ctx["apply"] = kwargs.get("apply_bc_fn")
        ctx["bridge"] = bridge_n["n"]
        return _orig_bridge(cells, cfg, **kwargs)

    def w_bed(*a, **k):
        cells = ctx["cells"]
        cfg = ctx["cfg"]
        apply = ctx["apply"]
        bed_index = int(k.get("cell_index", a[2] if len(a) > 2 else -1))
        before = _locate_max_abs(cells, cfg, apply)
        ok = bool(_orig_bed(*a, **k))
        after = _locate_max_abs(cells, cfg, apply)
        detonated = float(after["max_abs"]) > max(float(before["max_abs"]) * 5.0, 1.0)
        ev = {
            "bridge": int(ctx["bridge"]),
            "bed": bed_index,
            "lsq_ok": ok,
            "detonated": detonated,
            "before": {
                "rms": before["rms"],
                "max_abs": before["max_abs"],
                "argmax": {
                    kk: before[kk]
                    for kk in ("cell", "phase", "species", "kind", "F_hat")
                    if kk in before
                },
            },
            "after": {
                "rms": after["rms"],
                "max_abs": after["max_abs"],
                "argmax": {
                    kk: after[kk]
                    for kk in ("cell", "phase", "species", "kind", "F_hat")
                    if kk in after
                },
                "top5": after["top5"],
            },
            "d_rms": float(after["rms"] - before["rms"]),
            "d_max": float(after["max_abs"] - before["max_abs"]),
        }
        events.append(ev)
        if detonated or abs(ev["d_max"]) > 0.5:
            am = ev["after"]["argmax"]
            print(
                f"  bridge#{ev['bridge']} bed{bed_index} ok={ok} "
                f"max {before['max_abs']:.3e}→{after['max_abs']:.3e} "
                f"argmax→ cell{am.get('cell')} {am.get('phase')}/{am.get('species')} "
                f"F={am.get('F_hat', float('nan')):.3e} detonated={detonated}"
            )
        return ok

    def w_split(*a, **k):
        print(f"[split] bridge#{ctx['bridge']}")
        return _orig_split(*a, **k)

    pp.run_outer_preinner_start_value_bridge = w_bridge  # type: ignore[assignment]
    conn.run_outer_preinner_start_value_bridge = w_bridge  # type: ignore[assignment]
    pp.preproject_bed_major_gas_phase_split_to_exchange_closure = w_bed  # type: ignore[assignment]
    pp.run_outer_preinner_gas_phase_split_preprojection = w_split  # type: ignore[assignment]

    try:
        reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=dict(ISOMORPH))
        run_phase2_init(reactor)
        result, wall_s = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter))
        )
    finally:
        pp.run_outer_preinner_start_value_bridge = _orig_bridge  # type: ignore[assignment]
        conn.run_outer_preinner_start_value_bridge = _orig_bridge  # type: ignore[assignment]
        pp.preproject_bed_major_gas_phase_split_to_exchange_closure = _orig_bed  # type: ignore[assignment]
        pp.run_outer_preinner_gas_phase_split_preprojection = _orig_split  # type: ignore[assignment]

    assert isinstance(result, dict)
    dets = [e for e in events if e["detonated"]]
    # summary by bed
    by_bed: dict[str, int] = {}
    by_sp: dict[str, int] = {}
    for e in dets:
        by_bed[str(e["bed"])] = by_bed.get(str(e["bed"]), 0) + 1
        sp = str((e["after"].get("argmax") or {}).get("species", "?"))
        by_sp[sp] = by_sp.get(sp, 0) + 1
    out = {
        "note": "preinner split LSQ detonation bed×species; local2x2 post-gate stack",
        "max_global_iter": int(args.max_global_iter),
        "wall_s": float(wall_s),
        "exit_CO": float((result.get("bed_exit_gas") or {}).get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "n_events": len(events),
        "n_detonations": len(dets),
        "detonation_by_bed": by_bed,
        "detonation_argmax_species": by_sp,
        "events": events,
        "first_detonation": dets[0] if dets else None,
    }
    path = data_path("phase2_preinner_split_lsq_detonation_o6.json")
    write_json(path, out)
    print(
        f"done CO={out['exit_CO']:.4f} dets={len(dets)} by_bed={by_bed} by_sp={by_sp} "
        f"wall={wall_s:.0f}s → {path}"
    )


if __name__ == "__main__":
    main()
