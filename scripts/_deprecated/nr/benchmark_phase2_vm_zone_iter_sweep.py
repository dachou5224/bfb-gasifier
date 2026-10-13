#!/usr/bin/env python3
"""Phase2 LU：VM 脱挥发分区 (0 vs 3) × max_global_iter 收敛对比。"""

from __future__ import annotations

from _lu_harness import (
    LADDER_P2,
    add_max_global_iter_arg,
    apply_config_overrides,
    build_phase2_config,
    data_path,
    phase2_solve_kwargs,
    repo_root,
    run_phase2_init,
    write_json,
)

REPO_ROOT = repo_root()

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter

from src.core.connectivity import _measure_nr_init_scaled_metrics
from src.core.reactor import Reactor
REF_DRY = {"CO": 0.13, "CO2": 0.11}


def _run_case(*, vm_zone: int, max_global_iter: int, r5_scale: float | None = None) -> dict[str, object]:
    cfg = build_phase2_config()
    cfg.vorab_vm_devolatilization_zone_cells_thesis = int(vm_zone)
    if r5_scale is not None:
        cfg.r5_scale = float(r5_scale)
    apply_config_overrides(cfg, LADDER_P2)

    reactor = Reactor(cfg)
    t0 = perf_counter()
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    init_s = perf_counter() - t0
    solver_cells = reactor._solver_cells_for_nr()
    init_metrics = _measure_nr_init_scaled_metrics(
        solver_cells, cfg, apply_bc_fn=reactor._apply_all_bc_for_nr
    )
    abgleich = getattr(reactor, "_vorab_a_tier_init_abgleich_diag", None) or {}
    abgleich_final = abgleich.get("final") or {}

    t1 = perf_counter()
    result = reactor._solve_global_nr(
        max_iter=int(max_global_iter),
        tol=float(phase2_solve_kwargs().get("tol_global", 1.0)),
        verbose=False,
        init_strategy="vorabrechnung",
        jacobian_strategy=phase2_solve_kwargs().get("nr_jacobian_strategy"),
        skip_precalc=True,
    )
    solve_s = perf_counter() - t1

    bed_exit = result.get("bed_exit_gas_dry") or {}
    exit_dry = result.get("exit_gas_dry") or {}
    bed_indices = [
        i for i, c in enumerate(reactor.cells) if str(getattr(c, "cell_type", "bed")) == "bed"
    ]
    bottom_profile = []
    from src.core.cell_thermo_helpers import gas_mole_fractions
    from src.core.species import GAS_SPECIES_INDEX as idx

    for i in bed_indices[:6]:
        c = reactor.cells[i]
        y = gas_mole_fractions(c.N_b, c.N_d, "mixed")
        bottom_profile.append(
            {
                "cell": i,
                "T_K": float(c.T),
                "y_CO": float(y[idx["CO"]]),
                "y_CO2": float(y[idx["CO2"]]),
                "y_O2": float(y[idx["O2"]]),
            }
        )

    return {
        "vm_zone_cells": int(vm_zone),
        "r5_scale": float(cfg.r5_scale),
        "max_global_iter": int(max_global_iter),
        "init_wall_s": float(init_s),
        "solve_wall_s": float(solve_s),
        "init_rms_scaled": float(init_metrics["rms_scaled"]),
        "init_abgleich_aligned": bool(abgleich_final.get("aligned", False)),
        "init_product_holdup_rel_max": float(abgleich_final.get("product_source_holdup_rel_max", 0.0) or 0.0),
        "init_solid_reactive_rel_max": float(abgleich_final.get("solid_reactive_budget_rel_max", 0.0) or 0.0),
        "converged": bool(result.get("converged", False)),
        "converged_fully": bool(result.get("converged_fully", False)),
        "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
        "rms_scaled_final": float(result.get("rms_scaled_final", 0.0) or 0.0),
        "rms_scaled_gas_phase_split_final": float(result.get("rms_scaled_gas_phase_split_final", 0.0) or 0.0),
        "bed_exit_gas_dry": {k: float(bed_exit.get(k, 0.0)) for k in REF_DRY},
        "reactor_exit_gas_dry": {k: float(exit_dry.get(k, 0.0)) for k in REF_DRY},
        "ref_dry": dict(REF_DRY),
        "bottom_profile": bottom_profile,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=150)
    parser.add_argument("--r5-scale", type=float, default=None, help="覆盖 cfg.r5_scale（R5=0 对照）")
    parser.add_argument(
        "--zones",
        type=int,
        nargs="+",
        default=[0, 3],
        help="vorab_vm_devolatilization_zone_cells_thesis 取值列表",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_vm_zone_iter150_sweep.json",
    )
    args = parser.parse_args()

    results = [
        _run_case(
            vm_zone=int(z),
            max_global_iter=int(args.max_global_iter),
            r5_scale=args.r5_scale,
        )
        for z in args.zones
    ]
    payload = {
        "max_global_iter": int(args.max_global_iter),
        "r5_scale": float(args.r5_scale) if args.r5_scale is not None else None,
        "zones": [int(z) for z in args.zones],
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
