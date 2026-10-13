#!/usr/bin/env python3
"""zone=3（Kap.2.1）下 force bed0 T0 越过点火阈，对照 zone=0。

保留 extent limiter；仅收窄 bed0 温度围栏并设 T=T0。

Usage
-----
    python3 scripts/audit_vm_zone_force_T0_lu.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

import src.workflow.steps.init_precalc_step as ips
from _lu_harness import (
    REF_CO,
    REF_CO2,
    add_json_output_arg,
    add_max_global_iter_arg,
    bed_cell_indices,
    build_phase2_reactor,
    co_co2_score,
    phase2_solve_kwargs,
    timed_phase2_solve,
    write_json,
)
from src.core.species import GAS_SPECIES_INDEX as idx


def _y(cell: Any) -> dict[str, float]:
    y = cell._mole_fractions("combined")
    return {sp: float(y[j]) for sp, j in idx.items()}


def _summarize(reactor: Any, result: dict[str, Any], *, wall_s: float, label: str) -> dict[str, Any]:
    beds = bed_cell_indices(reactor)
    axial = []
    for i in beds:
        c = reactor.cells[i]
        y = _y(c)
        axial.append(
            {
                "i": i,
                "T_K": float(c.T),
                "y_O2": y.get("O2", 0.0),
                "y_CO": y.get("CO", 0.0),
                "y_CO2": y.get("CO2", 0.0),
                "y_H2": y.get("H2", 0.0),
                "y_CH4": y.get("CH4", 0.0),
            }
        )
    T_prof = [r["T_K"] for r in axial]
    exit_dry = {k: float(v) for k, v in (result.get("exit_gas_dry") or {}).items()}
    gas0 = np.asarray(reactor.cells[beds[0]]._vm_gas_source_cache, dtype=np.float64)
    vm_shares = []
    for i in beds[:4]:
        g = np.asarray(reactor.cells[i]._vm_gas_source_cache, dtype=np.float64)
        vm_shares.append({"i": i, "vm_gas_src_mol_s": float(np.sum(g))})
    row = {
        "label": label,
        "wall_s": wall_s,
        "vm_zone": int(getattr(reactor.config, "vorab_vm_devolatilization_zone_cells_thesis", 0)),
        "converged": bool(result.get("converged")),
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "T0_K": float(T_prof[0]) if T_prof else float("nan"),
        "T_bed_max_K": float(np.max(T_prof)) if T_prof else float("nan"),
        "T_bed_mean_K": float(np.mean(T_prof)) if T_prof else float("nan"),
        "T_bed_profile_K": T_prof,
        "bed_top_y_O2": float(axial[-1]["y_O2"]) if axial else float("nan"),
        "bed0_y_O2": float(axial[0]["y_O2"]) if axial else float("nan"),
        "exit_dry": exit_dry,
        "co_co2_score_vs_table2": co_co2_score(
            float(exit_dry.get("CO", float("nan"))),
            float(exit_dry.get("CO2", float("nan"))),
            ref_co=REF_CO,
            ref_co2=REF_CO2,
        ),
        "vm_sources": vm_shares,
        "bed0_vm_CO_H2_CH4": {
            "CO": float(gas0[idx["CO"]]),
            "H2": float(gas0[idx["H2"]]),
            "CH4": float(gas0[idx["CH4"]]),
        },
        "axial_bed": axial,
        "bed0_fence": [
            float(getattr(reactor.cells[beds[0]], "_nr_temperature_min_K", float("nan"))),
            float(getattr(reactor.cells[beds[0]], "_nr_temperature_max_K", float("nan"))),
        ]
        if beds
        else [float("nan"), float("nan")],
    }
    print(
        f"{label}: zone={row['vm_zone']} conv={row['converged_fully']} "
        f"rms={row['rms_scaled_final']:.3e} T0={row['T0_K']:.1f} "
        f"fence={row['bed0_fence']} bedTopO2={row['bed_top_y_O2']:.4f} "
        f"ex CO/CO2/H2/O2="
        f"{exit_dry.get('CO', float('nan')):.4f}/"
        f"{exit_dry.get('CO2', float('nan')):.4f}/"
        f"{exit_dry.get('H2', float('nan')):.4f}/"
        f"{exit_dry.get('O2', float('nan')):.2e} "
        f"score={row['co_co2_score_vs_table2']:.3f} wall={wall_s:.1f}s",
        flush=True,
    )
    return row


def run_force(
    *,
    vm_zone: int,
    T0_target: float,
    half: float,
    max_global_iter: int,
    force: bool,
) -> dict[str, Any]:
    _orig = ips._install_temperature_fence_from_vorabrechnung

    def _install(reactor, T_est):
        _orig(reactor, T_est)
        if not force:
            return
        beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if not beds:
            return
        b0 = beds[0]
        b0._nr_temperature_min_K = float(T0_target - half)
        b0._nr_temperature_max_K = float(T0_target + half)
        b0.T = float(T0_target)
        b0._vorab_temperature_estimate_K = float(T0_target)

    label = (
        f"zone{vm_zone}_force_T0={T0_target:.0f}_pm{half:.0f}"
        if force
        else f"zone{vm_zone}_default"
    )
    ips._install_temperature_fence_from_vorabrechnung = _install
    try:
        r = build_phase2_reactor(
            ladder_level="p2",
            builder_kwargs={"enable_vorab_march_staged": True},
            cfg_overrides={
                "vorab_vm_devolatilization_zone_cells_thesis": int(vm_zone),
                "vorab_vm_devolatilization_zone_profile_thesis": "linear",
            },
        )
        out, wall = timed_phase2_solve(
            r, **phase2_solve_kwargs(max_global_iter=int(max_global_iter))
        )
        return _summarize(r, out, wall_s=wall, label=label)
    finally:
        ips._install_temperature_fence_from_vorabrechnung = _orig


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(p, default=8)
    add_json_output_arg(p, "phase2_vm_zone_force_T0_o8.json")
    p.add_argument("--targets", type=float, nargs="+", default=[1100.0])
    p.add_argument("--half", type=float, default=50.0)
    p.add_argument("--zones", type=str, default="3,0")
    p.add_argument("--include-default", action="store_true", help="also run no-force cold seed")
    args = p.parse_args()
    zones = [int(x.strip()) for x in str(args.zones).split(",") if x.strip()]

    rows: list[dict[str, Any]] = []
    if args.include_default:
        for z in zones:
            print(f"=== zone={z} default (no force) ===", flush=True)
            rows.append(
                run_force(
                    vm_zone=z,
                    T0_target=900.0,
                    half=50.0,
                    max_global_iter=args.max_global_iter,
                    force=False,
                )
            )
    for z in zones:
        for T0 in args.targets:
            print(f"=== zone={z} force_T0={T0:.0f} ±{args.half:.0f} ===", flush=True)
            rows.append(
                run_force(
                    vm_zone=z,
                    T0_target=float(T0),
                    half=float(args.half),
                    max_global_iter=args.max_global_iter,
                    force=True,
                )
            )

    # delta zone3_force - zone0_force for same T0
    deltas = []
    by = {r["label"]: r for r in rows}
    for T0 in args.targets:
        a = by.get(f"zone3_force_T0={T0:.0f}_pm{args.half:.0f}")
        b = by.get(f"zone0_force_T0={T0:.0f}_pm{args.half:.0f}")
        if a and b:
            ea, eb = a["exit_dry"], b["exit_dry"]
            d = {
                "T0_target": float(T0),
                "d_T0_K": a["T0_K"] - b["T0_K"],
                "d_bed_top_y_O2": a["bed_top_y_O2"] - b["bed_top_y_O2"],
                "d_rms": a["rms_scaled_final"] - b["rms_scaled_final"],
                "d_exit_CO": ea.get("CO", 0.0) - eb.get("CO", 0.0),
                "d_exit_CO2": ea.get("CO2", 0.0) - eb.get("CO2", 0.0),
                "d_exit_H2": ea.get("H2", 0.0) - eb.get("H2", 0.0),
                "d_exit_O2": ea.get("O2", 0.0) - eb.get("O2", 0.0),
                "d_score": a["co_co2_score_vs_table2"] - b["co_co2_score_vs_table2"],
            }
            deltas.append(d)
            print(
                f"delta(z3-z0)@force{T0:.0f}: dT0={d['d_T0_K']:+.1f} "
                f"dBedO2={d['d_bed_top_y_O2']:+.4f} dCO={d['d_exit_CO']:+.4f} "
                f"dCO2={d['d_exit_CO2']:+.4f} dH2={d['d_exit_H2']:+.4f} "
                f"dO2={d['d_exit_O2']:+.2e} d_score={d['d_score']:+.3f}",
                flush=True,
            )

    payload = {
        "policy": "keep lim ON; force bed0 T seed+fence; A/B vm_zone 3 vs 0",
        "max_global_iter": int(args.max_global_iter),
        "half_K": float(args.half),
        "cases": rows,
        "delta_zone3_minus_zone0_at_force": deltas,
    }
    out = write_json(args.output, payload)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
