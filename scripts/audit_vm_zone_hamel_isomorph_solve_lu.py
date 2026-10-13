#!/usr/bin/env python3
"""Phase2：Hamel Kap.2.1 多格 VM 释放（zone=3）vs legacy bed0 瞬时（zone=0）求解对照。

Usage
-----
    python3 scripts/audit_vm_zone_hamel_isomorph_solve_lu.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

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


def _vm_sources(reactor: Any) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for i in bed_cell_indices(reactor)[:6]:
        c = reactor.cells[i]
        gas = np.asarray(getattr(c, "_vm_gas_source_cache", np.zeros(1)), dtype=np.float64)
        rows.append(
            {
                "i": float(i),
                "vm_gas_src_mol_s": float(np.sum(gas)),
                "vm_CO": float(gas[idx["CO"]]) if gas.size > idx["CO"] else 0.0,
                "vm_H2": float(gas[idx["H2"]]) if gas.size > idx["H2"] else 0.0,
                "vm_CH4": float(gas[idx["CH4"]]) if gas.size > idx["CH4"] else 0.0,
                "vm_CO2": float(gas[idx["CO2"]]) if gas.size > idx["CO2"] else 0.0,
                "cache_valid": float(bool(getattr(c, "_vm_cache_valid", False))),
            }
        )
    return rows


def _summarize(reactor: Any, result: dict[str, Any], *, wall_s: float, label: str) -> dict[str, Any]:
    beds = bed_cell_indices(reactor)
    axial = []
    for i in beds:
        c = reactor.cells[i]
        y = _y(c)
        axial.append(
            {
                "i": i,
                "z_m": float(c.geo.h_center),
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
    bed_exit = {k: float(v) for k, v in (result.get("bed_exit_gas_dry") or {}).items()}
    fb = list(getattr(reactor, "freeboard_cells", []) or [])
    score = co_co2_score(
        float(exit_dry.get("CO", float("nan"))),
        float(exit_dry.get("CO2", float("nan"))),
    )
    row = {
        "label": label,
        "wall_s": wall_s,
        "vm_zone": int(getattr(reactor.config, "vorab_vm_devolatilization_zone_cells_thesis", 0)),
        "vm_profile": str(
            getattr(reactor.config, "vorab_vm_devolatilization_zone_profile_thesis", "")
        ),
        "r5_scale": float(reactor.config.r5_scale),
        "converged": bool(result.get("converged")),
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "T0_K": float(T_prof[0]) if T_prof else float("nan"),
        "T_bed_mean_K": float(np.mean(T_prof)) if T_prof else float("nan"),
        "T_bed_max_K": float(np.max(T_prof)) if T_prof else float("nan"),
        "T_bed_profile_K": T_prof,
        "T_fb_profile_K": [float(c.T) for c in fb],
        "bed_top_y_O2": float(axial[-1]["y_O2"]) if axial else float("nan"),
        "exit_dry": exit_dry,
        "bed_exit_dry": bed_exit,
        "co_co2_score_vs_table2": score,
        "ref_table2": {"CO": REF_CO, "CO2": REF_CO2},
        "vm_sources_init_or_final": _vm_sources(reactor),
        "axial_bed": axial,
    }
    print(
        f"{label}: zone={row['vm_zone']} conv={row['converged_fully']} "
        f"rms={row['rms_scaled_final']:.3e} n_iter={row['n_iter']} "
        f"T0={row['T0_K']:.1f} Tmax={row['T_bed_max_K']:.1f} "
        f"bedTopO2={row['bed_top_y_O2']:.4f} "
        f"ex CO/CO2/H2/CH4/O2="
        f"{exit_dry.get('CO', float('nan')):.4f}/"
        f"{exit_dry.get('CO2', float('nan')):.4f}/"
        f"{exit_dry.get('H2', float('nan')):.4f}/"
        f"{exit_dry.get('CH4', float('nan')):.5f}/"
        f"{exit_dry.get('O2', float('nan')):.2e} "
        f"score={score:.3f} wall={wall_s:.1f}s",
        flush=True,
    )
    return row


def run_case(*, vm_zone: int, max_global_iter: int) -> dict[str, Any]:
    label = f"vm_zone_{vm_zone}"
    reactor = build_phase2_reactor(
        ladder_level="p2",
        builder_kwargs={"enable_vorab_march_staged": True},
        cfg_overrides={
            "vorab_vm_devolatilization_zone_cells_thesis": int(vm_zone),
            "vorab_vm_devolatilization_zone_profile_thesis": "linear",
        },
    )
    kw = phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False)
    result, wall_s = timed_phase2_solve(reactor, **kw)
    return _summarize(reactor, result, wall_s=wall_s, label=label)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(p, default=8)
    add_json_output_arg(p, "phase2_vm_zone_hamel_isomorph_o8.json")
    p.add_argument(
        "--zones",
        type=str,
        default="3,0",
        help="comma-separated vm_zone values (default: 3,0)",
    )
    args = p.parse_args()
    zones = [int(x.strip()) for x in str(args.zones).split(",") if x.strip()]
    rows = []
    for z in zones:
        rows.append(run_case(vm_zone=z, max_global_iter=int(args.max_global_iter)))
    by_label = {r["label"]: r for r in rows}
    delta = None
    if "vm_zone_3" in by_label and "vm_zone_0" in by_label:
        a, b = by_label["vm_zone_3"], by_label["vm_zone_0"]
        ea, eb = a["exit_dry"], b["exit_dry"]
        delta = {
            "d_T0_K": a["T0_K"] - b["T0_K"],
            "d_T_bed_max_K": a["T_bed_max_K"] - b["T_bed_max_K"],
            "d_bed_top_y_O2": a["bed_top_y_O2"] - b["bed_top_y_O2"],
            "d_rms": a["rms_scaled_final"] - b["rms_scaled_final"],
            "d_exit_CO": ea.get("CO", 0.0) - eb.get("CO", 0.0),
            "d_exit_CO2": ea.get("CO2", 0.0) - eb.get("CO2", 0.0),
            "d_exit_H2": ea.get("H2", 0.0) - eb.get("H2", 0.0),
            "d_exit_CH4": ea.get("CH4", 0.0) - eb.get("CH4", 0.0),
            "d_exit_O2": ea.get("O2", 0.0) - eb.get("O2", 0.0),
            "d_co_co2_score": a["co_co2_score_vs_table2"] - b["co_co2_score_vs_table2"],
        }
        print(
            "delta(zone3-zone0): "
            f"dT0={delta['d_T0_K']:+.1f} dTmax={delta['d_T_bed_max_K']:+.1f} "
            f"dBedO2={delta['d_bed_top_y_O2']:+.4f} d_rms={delta['d_rms']:+.3e} "
            f"dCO={delta['d_exit_CO']:+.4f} dCO2={delta['d_exit_CO2']:+.4f} "
            f"dH2={delta['d_exit_H2']:+.4f} dCH4={delta['d_exit_CH4']:+.5f}",
            flush=True,
        )
    payload = {
        "max_global_iter": int(args.max_global_iter),
        "ladder": "p2",
        "note": (
            "Hamel Kap.2.1 multi-cell VM release (zone=3 linear) vs legacy "
            "bed0 dump (zone=0). Production knobs kept (r5_scale=1.0, R12 off)."
        ),
        "cases": rows,
        "delta_zone3_minus_zone0": delta,
    }
    out = write_json(args.output, payload)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
