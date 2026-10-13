#!/usr/bin/env python3
"""保留 extent limiter，强制 bed0 温度种子越过点火阈。

做法：在 Vorab 装温度围栏后，把 bed0 的
``(_nr_temperature_min_K, _nr_temperature_max_K)`` 收窄到 ``T0±half``，
并把 ``bed0.T`` 设为 ``T0``。动力学 knobs 保持生产口径（R12 关、lim 开）。

Usage
-----
    python3 scripts/audit_hamel_force_T0_keep_lim_lu.py --max-global-iter 8
"""

from __future__ import annotations

from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    bed_cell_indices,
    build_phase2_reactor,
    phase2_solve_kwargs,
    timed_phase2_solve,
    write_json,
)

import argparse
from typing import Any

import numpy as np

import src.workflow.steps.init_precalc_step as ips
from src.core.species import GAS_SPECIES_INDEX as idx


def _y(cell: Any) -> dict[str, float]:
    y = cell._mole_fractions("combined")
    return {sp: float(y[j]) for sp, j in idx.items()}


def _smooth(T: list[float]) -> dict[str, float]:
    arr = np.asarray(T, dtype=float)
    if arr.size < 2:
        return {"T_std_K": float("nan"), "T_max_abs_dT_K": float("nan"), "T_range_K": float("nan")}
    return {
        "T_std_K": float(np.std(arr)),
        "T_max_abs_dT_K": float(np.max(np.abs(np.diff(arr)))),
        "T_range_K": float(np.max(arr) - np.min(arr)),
    }


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
                "fence": [
                    float(getattr(c, "_nr_temperature_min_K", float("nan"))),
                    float(getattr(c, "_nr_temperature_max_K", float("nan"))),
                ],
            }
        )
    T_prof = [r["T_K"] for r in axial]
    exit_dry = {k: float(v) for k, v in (result.get("exit_gas_dry") or {}).items()}
    bed_exit = {k: float(v) for k, v in (result.get("bed_exit_gas_dry") or {}).items()}
    fb = list(getattr(reactor, "freeboard_cells", []) or [])
    row = {
        "label": label,
        "wall_s": wall_s,
        "converged": bool(result.get("converged")),
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "T0_K": float(T_prof[0]) if T_prof else float("nan"),
        "T_bed_mean_K": float(np.mean(T_prof)) if T_prof else float("nan"),
        "T_bed_profile_K": T_prof,
        "T_fb_profile_K": [float(c.T) for c in fb],
        "bed_top_y_O2": float(axial[-1]["y_O2"]) if axial else float("nan"),
        "exit_dry": exit_dry,
        "bed_exit_dry": bed_exit,
        "bed0_fence": axial[0]["fence"] if axial else [float("nan"), float("nan")],
        "smoothness": _smooth(T_prof),
        "axial_bed": axial,
        "cfg_note": {
            "enable_r12": bool(reactor.config.enable_r12),
            "r5_scale": float(reactor.config.r5_scale),
            "extent_limiters": "ON (production default)",
            "seed_bottom_K": float(
                getattr(reactor.config, "vorab_nr_temperature_seed_bottom_K_thesis", float("nan"))
            ),
            "cold_anchor": bool(
                getattr(reactor.config, "vorab_nr_temperature_seed_bed0_cold_anchor_thesis", False)
            ),
        },
    }
    sm = row["smoothness"]
    print(
        f"{label}: conv={row['converged_fully']} rms={row['rms_scaled_final']:.3e} "
        f"T0={row['T0_K']:.1f} fence={row['bed0_fence']} "
        f"dTmax={sm['T_max_abs_dT_K']:.1f} bedO2={row['bed_top_y_O2']:.4f} "
        f"ex CO/H2/O2="
        f"{exit_dry.get('CO', float('nan')):.4f}/"
        f"{exit_dry.get('H2', float('nan')):.4f}/"
        f"{exit_dry.get('O2', float('nan')):.4f} "
        f"wall={wall_s:.1f}s",
        flush=True,
    )
    return row


def run_default(*, max_global_iter: int) -> dict[str, Any]:
    r = build_phase2_reactor(ladder_level="p2")
    out, wall = timed_phase2_solve(
        r, **phase2_solve_kwargs(max_global_iter=max_global_iter)
    )
    return _summarize(r, out, wall_s=wall, label="P0_default")


def run_force_T0(T0_target: float, *, half: float, max_global_iter: int) -> dict[str, Any]:
    _orig = ips._install_temperature_fence_from_vorabrechnung

    def _install(reactor, T_est):
        _orig(reactor, T_est)
        beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if not beds:
            return
        b0 = beds[0]
        b0._nr_temperature_min_K = float(T0_target - half)
        b0._nr_temperature_max_K = float(T0_target + half)
        b0.T = float(T0_target)
        b0._vorab_temperature_estimate_K = float(T0_target)

    ips._install_temperature_fence_from_vorabrechnung = _install
    try:
        r = build_phase2_reactor(ladder_level="p2")
        out, wall = timed_phase2_solve(
            r, **phase2_solve_kwargs(max_global_iter=max_global_iter)
        )
        return _summarize(
            r, out, wall_s=wall, label=f"force_T0={T0_target:.0f}_pm{half:.0f}"
        )
    finally:
        ips._install_temperature_fence_from_vorabrechnung = _orig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(parser, default=8)
    add_json_output_arg(parser, "phase2_hamel_force_T0_keep_lim_o8.json")
    parser.add_argument(
        "--targets",
        type=float,
        nargs="+",
        default=[1100.0, 1150.0],
        help="bed0 强制目标温度 [K]",
    )
    parser.add_argument("--half", type=float, default=50.0, help="围栏半宽 [K]")
    parser.add_argument("--skip-default", action="store_true")
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []
    if not args.skip_default:
        print("=== P0_default (lim ON, no force) ===", flush=True)
        rows.append(run_default(max_global_iter=args.max_global_iter))
    for T0 in args.targets:
        print(f"=== force_T0={T0:.0f} ±{args.half:.0f} (lim ON) ===", flush=True)
        rows.append(
            run_force_T0(T0, half=args.half, max_global_iter=args.max_global_iter)
        )

    base = next((r for r in rows if r["label"] == "P0_default"), None)
    deltas = []
    if base is not None:
        for r in rows:
            if r is base:
                continue
            deltas.append(
                {
                    "label": r["label"],
                    "dT0_K": r["T0_K"] - base["T0_K"],
                    "d_bed_top_y_O2": r["bed_top_y_O2"] - base["bed_top_y_O2"],
                    "d_exit_CO": r["exit_dry"].get("CO", 0.0) - base["exit_dry"].get("CO", 0.0),
                    "d_exit_H2": r["exit_dry"].get("H2", 0.0) - base["exit_dry"].get("H2", 0.0),
                    "d_exit_O2": r["exit_dry"].get("O2", 0.0) - base["exit_dry"].get("O2", 0.0),
                    "d_rms": r["rms_scaled_final"] - base["rms_scaled_final"],
                    "d_T_max_abs_dT": (
                        r["smoothness"]["T_max_abs_dT_K"] - base["smoothness"]["T_max_abs_dT_K"]
                    ),
                }
            )
        print("=== delta vs P0 ===", flush=True)
        for d in deltas:
            print(
                f"  {d['label']}: dT0={d['dT0_K']:+.1f} d_bedO2={d['d_bed_top_y_O2']:+.4f} "
                f"d_exCO={d['d_exit_CO']:+.4f} d_exH2={d['d_exit_H2']:+.4f} "
                f"d_exO2={d['d_exit_O2']:+.4f} d_rms={d['d_rms']:+.3e} "
                f"d_dTmax={d['d_T_max_abs_dT']:+.1f}",
                flush=True,
            )

    payload = {
        "policy": "keep extent limiters ON; force bed0 T seed+fence past ignition",
        "max_global_iter": args.max_global_iter,
        "half_K": args.half,
        "cases": rows,
        "deltas_vs_P0": deltas,
    }
    write_json(args.output, payload)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
