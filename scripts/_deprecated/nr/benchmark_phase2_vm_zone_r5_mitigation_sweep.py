#!/usr/bin/env python3
"""VM zone=3：R5 缓解策略扫描（全局 r5_scale / fixed50 / bed1 局部 R5）。"""

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
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

from src.core.connectivity import _measure_nr_init_scaled_metrics
from src.core.reactor import Reactor
REF_DRY = {"CO": 0.13, "CO2": 0.11}


@dataclass(frozen=True)
class MitigationCase:
    label: str
    r5_scale: float = 0.75
    bed_r5_multipliers: tuple[float, ...] = ()
    fixed50: bool = False
    product_holdup_closure: bool = True
    init_abgleich_bridge: bool = True


DEFAULT_CASES: tuple[MitigationCase, ...] = (
    MitigationCase("baseline_r5_075"),
    MitigationCase("global_r5_050", r5_scale=0.50),
    MitigationCase("global_r5_025", r5_scale=0.25),
    MitigationCase("bed1_r5_x025", bed_r5_multipliers=(1.0, 0.25, 1.0)),
    MitigationCase("bed1_r5_x010", bed_r5_multipliers=(1.0, 0.10, 1.0)),
    MitigationCase("bed1_r5_x000", bed_r5_multipliers=(1.0, 0.0, 1.0)),
    MitigationCase("fixed50_r5_075", fixed50=True),
    MitigationCase("fixed50_r5_050", fixed50=True, r5_scale=0.50),
    MitigationCase("fixed50_r5_025", fixed50=True, r5_scale=0.25),
    MitigationCase("fixed50_bed1_x025", fixed50=True, bed_r5_multipliers=(1.0, 0.25, 1.0)),
    MitigationCase("fixed50_bed1_x000", fixed50=True, bed_r5_multipliers=(1.0, 0.0, 1.0)),
)


def _apply_bed_r5_multipliers(reactor: Reactor, multipliers: tuple[float, ...]) -> list[float]:
    """按 bed 序号缩放 cell.r5_scale（VM 分区 bed1 局部限 R5）。"""
    applied: list[float] = []
    bed_rank = 0
    for cell in reactor.cells:
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        if bed_rank >= len(multipliers):
            break
        factor = float(max(multipliers[bed_rank], 0.0))
        cell.r5_scale = float(max(cell.r5_scale * factor, 0.0))
        applied.append(float(cell.r5_scale))
        bed_rank += 1
    return applied


def _apply_fixed50(cfg) -> None:
    cfg.gas_inlet_split_strategy = "fixed"
    cfg.gas_inlet_dense_frac = 0.50
    cfg.vorab_staged_primary_gas_dense_frac_thesis = 0.50


def _run_mitigation_case(case: MitigationCase, *, vm_zone: int, max_global_iter: int) -> dict[str, object]:
    cfg = build_phase2_config()
    cfg.vorab_vm_devolatilization_zone_cells_thesis = int(vm_zone)
    cfg.r5_scale = float(case.r5_scale)
    cfg.vorab_a_tier_product_holdup_source_closure_thesis = bool(case.product_holdup_closure)
    cfg.vorab_a_tier_init_abgleich_bridge_thesis = bool(case.init_abgleich_bridge)
    if case.fixed50:
        _apply_fixed50(cfg)
    for key, value in {
        "nr_line_search_gas_phase_split_merit_thesis": True,
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1180.0,
        "freeboard_trajectory_coeff_model": "exact_hamel",
    }.items():
        setattr(cfg, key, value)

    reactor = Reactor(cfg)
    bed_r5_applied = _apply_bed_r5_multipliers(reactor, case.bed_r5_multipliers)

    t0 = perf_counter()
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    init_s = perf_counter() - t0
    init_metrics = _measure_nr_init_scaled_metrics(
        reactor._solver_cells_for_nr(),
        cfg,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
    )
    abgleich_final = (getattr(reactor, "_vorab_a_tier_init_abgleich_diag", None) or {}).get("final") or {}

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
    from src.core.cell_thermo_helpers import gas_mole_fractions
    from src.core.species import GAS_SPECIES_INDEX as idx

    bed_indices = [i for i, c in enumerate(reactor.cells) if str(getattr(c, "cell_type", "bed")) == "bed"]
    bottom_profile = []
    for i in bed_indices[:3]:
        c = reactor.cells[i]
        y = gas_mole_fractions(c.N_b, c.N_d, "mixed")
        bottom_profile.append(
            {
                "cell": i,
                "T_K": float(c.T),
                "r5_scale": float(c.r5_scale),
                "y_CO": float(y[idx["CO"]]),
                "y_CO2": float(y[idx["CO2"]]),
                "y_O2": float(y[idx["O2"]]),
            }
        )

    co_err = abs(float(bed_exit.get("CO", 0.0)) - REF_DRY["CO"]) / REF_DRY["CO"]
    co2_err = abs(float(bed_exit.get("CO2", 0.0)) - REF_DRY["CO2"]) / REF_DRY["CO2"]
    return {
        "label": case.label,
        "vm_zone_cells": int(vm_zone),
        "r5_scale_cfg": float(case.r5_scale),
        "bed_r5_multipliers": list(case.bed_r5_multipliers),
        "bed_r5_applied": bed_r5_applied,
        "fixed50": bool(case.fixed50),
        "product_holdup_closure": bool(case.product_holdup_closure),
        "init_abgleich_bridge": bool(case.init_abgleich_bridge),
        "max_global_iter": int(max_global_iter),
        "init_wall_s": float(init_s),
        "solve_wall_s": float(solve_s),
        "init_rms_scaled": float(init_metrics["rms_scaled"]),
        "init_abgleich_aligned": bool(abgleich_final.get("aligned", False)),
        "converged_fully": bool(result.get("converged_fully", False)),
        "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
        "rms_scaled_final": float(result.get("rms_scaled_final", 0.0) or 0.0),
        "bed_exit_gas_dry": {k: float(bed_exit.get(k, 0.0)) for k in REF_DRY},
        "bed_exit_co_co2_l1_err": float(co_err + co2_err),
        "bottom_profile": bottom_profile,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vm-zone", type=int, default=3)
    parser.add_argument("--max-global-iter", type=int, default=40)
    parser.add_argument(
        "--labels",
        nargs="*",
        default=None,
        help="只跑指定 label（默认全表）",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_vm_zone_r5_mitigation_sweep.json",
    )
    args = parser.parse_args()

    label_set = set(args.labels) if args.labels else None
    cases = [c for c in DEFAULT_CASES if label_set is None or c.label in label_set]
    if not cases:
        raise SystemExit("No matching mitigation cases")

    results = [_run_mitigation_case(c, vm_zone=int(args.vm_zone), max_global_iter=int(args.max_global_iter)) for c in cases]
    results.sort(key=lambda r: (float(r["rms_scaled_final"]), float(r["bed_exit_co_co2_l1_err"])))

    payload = {
        "vm_zone_cells": int(args.vm_zone),
        "max_global_iter": int(args.max_global_iter),
        "results": results,
        "best_by_rms": results[0]["label"] if results else None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
