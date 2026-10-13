#!/usr/bin/env python3
"""bed0 多 T 种子扫描 + 热枝偏好选解（对齐气流床 cell0 策略）。

对每个 T0 种子：放宽 bed0 围栏（种子在区间内），设 T=T0，跑 Phase2 oN；
选解优先 ignited（T0≥阈）、再收敛、再低 rms、再高出口 CO。

Usage
-----
    PYTHONPATH=.:scripts python3 scripts/audit_bed0_multistart_hot_branch_lu.py \\
        --max-global-iter 8 --profile front_loaded
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

IGNITED_T_K = 1200.0  # [K] 与气流床 _select_best_multistart 同量级


def _y(cell: Any) -> dict[str, float]:
    y = cell._mole_fractions("combined")
    return {sp: float(y[j]) for sp, j in idx.items()}


def _summarize(
    reactor: Any,
    result: dict[str, Any],
    *,
    wall_s: float,
    label: str,
    T0_seed: float,
) -> dict[str, Any]:
    beds = bed_cell_indices(reactor)
    axial = []
    o2_cons: list[float] = []
    for i in beds:
        c = reactor.cells[i]
        y = _y(c)
        nin = float(
            (c.N_b_in + c.N_d_in + c.N_zu_b + c.N_zu_d + c.N_rez_b + c.N_rez_d)[idx["O2"]]
        )
        nout = float((c.N_b + c.N_d)[idx["O2"]])
        o2_cons.append(nin - nout)
        axial.append(
            {
                "i": i,
                "T_K": float(c.T),
                "y_O2": y.get("O2", 0.0),
                "y_CO": y.get("CO", 0.0),
                "y_CO2": y.get("CO2", 0.0),
                "y_H2": y.get("H2", 0.0),
                "E_MW": float(c.calc_energy_balance()) / 1e6,
            }
        )
    exit_dry = {k: float(v) for k, v in (result.get("exit_gas_dry") or {}).items()}
    T0 = float(axial[0]["T_K"]) if axial else float("nan")
    row = {
        "label": label,
        "T0_seed_K": float(T0_seed),
        "wall_s": wall_s,
        "profile": str(
            getattr(reactor.config, "vorab_vm_devolatilization_zone_profile_thesis", "")
        ),
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "T0_K": T0,
        "T1_K": float(axial[1]["T_K"]) if len(axial) > 1 else float("nan"),
        "T_bed_max_K": float(max(a["T_K"] for a in axial)) if axial else float("nan"),
        "ignited": bool(T0 >= IGNITED_T_K),
        "bed0_dO2": float(o2_cons[0]) if o2_cons else float("nan"),
        "bed_top_y_O2": float(axial[-1]["y_O2"]) if axial else float("nan"),
        "E0_MW": float(axial[0]["E_MW"]) if axial else float("nan"),
        "exit_dry": exit_dry,
        "exit_CO": float(exit_dry.get("CO", 0.0)),
        "exit_CO2": float(exit_dry.get("CO2", 0.0)),
        "co_co2_score_vs_table2": co_co2_score(
            float(exit_dry.get("CO", float("nan"))),
            float(exit_dry.get("CO2", float("nan"))),
            ref_co=REF_CO,
            ref_co2=REF_CO2,
        ),
        "bed0_fence": [
            float(getattr(reactor.cells[beds[0]], "_nr_temperature_min_K", float("nan"))),
            float(getattr(reactor.cells[beds[0]], "_nr_temperature_max_K", float("nan"))),
        ]
        if beds
        else [float("nan"), float("nan")],
        "axial_bed": axial,
    }
    print(
        f"{label}: seed={T0_seed:.0f} conv={row['converged_fully']} "
        f"ign={row['ignited']} rms={row['rms_scaled_final']:.3e} "
        f"T0={row['T0_K']:.0f} dO2={row['bed0_dO2']:.2f} "
        f"topO2={row['bed_top_y_O2']:.4f} CO={row['exit_CO']:.4f} "
        f"E0={row['E0_MW']:+.3f}MW score={row['co_co2_score_vs_table2']:.3f}",
        flush=True,
    )
    return row


def _select_best(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    """对齐气流床：优先自持热枝，再收敛，再低 rms，再高 CO。

    自持热枝 = T0≥点火阈 **且** 未贴死在围栏下沿（否则只是 force 假热枝）。
    """

    def _self_sustaining(row: dict[str, Any]) -> bool:
        """热枝自持：点火温度 + 未贴下沿 + 能量不强烈求冷 + 有 O₂ 消耗。"""
        if not row.get("ignited"):
            return False
        fence = row.get("bed0_fence") or [float("nan"), float("nan")]
        t_lo = float(fence[0])
        if np.isfinite(t_lo) and float(row["T0_K"]) < t_lo + 25.0:
            return False
        if float(row.get("E0_MW", 0.0)) < -1.0:
            return False
        if float(row.get("bed0_dO2", 0.0)) < 0.5:
            return False
        if float(row.get("rms_scaled_final", 1e9)) > 0.08:
            return False
        return True

    best: dict[str, Any] | None = None
    for row in rows:
        if best is None:
            best = row
            continue
        r_ss, b_ss = _self_sustaining(row), _self_sustaining(best)
        if r_ss and not b_ss:
            best = row
            continue
        if r_ss != b_ss:
            continue
        # 无自持时：收敛优先，再 rms / CO / T0
        if row["converged_fully"] and not best["converged_fully"]:
            best = row
            continue
        if row["converged_fully"] != best["converged_fully"]:
            continue
        r_rms, b_rms = row["rms_scaled_final"], best["rms_scaled_final"]
        if r_rms < b_rms * 0.95:
            best = row
            continue
        if abs(r_rms - b_rms) <= max(b_rms * 0.1, 1e-4):
            if row["exit_CO"] > best["exit_CO"] + 1e-4:
                best = row
            elif abs(row["exit_CO"] - best["exit_CO"]) <= 1e-4 and row["T0_K"] > best["T0_K"]:
                best = row
        elif r_rms < b_rms:
            best = row
    return best


def run_seed(
    *,
    profile: str,
    T0_seed: float,
    half_lo: float,
    half_hi: float,
    max_global_iter: int,
    disable_cold_anchor: bool,
) -> dict[str, Any]:
    _orig = ips._install_temperature_fence_from_vorabrechnung

    def _install(reactor, T_est):
        _orig(reactor, T_est)
        beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if not beds:
            return
        b0 = beds[0]
        # 软围栏：种子在区间内，上界放宽以允许向热枝走
        t_lo = float(max(700.0, T0_seed - half_lo))
        t_hi = float(min(2500.0, max(T0_seed + half_hi, 1400.0)))
        b0._nr_temperature_min_K = t_lo
        b0._nr_temperature_max_K = t_hi
        b0.T = float(T0_seed)
        b0._vorab_temperature_estimate_K = float(T0_seed)

    overrides: dict[str, object] = {
        "vorab_vm_devolatilization_zone_cells_thesis": 3,
        "vorab_vm_devolatilization_zone_profile_thesis": profile,
    }
    if disable_cold_anchor:
        overrides["vorab_nr_temperature_seed_bed0_cold_anchor_thesis"] = False
        overrides["vorab_nr_temperature_seed_bottom_K_thesis"] = float(T0_seed)

    label = f"{profile}_seed{T0_seed:.0f}"
    ips._install_temperature_fence_from_vorabrechnung = _install
    try:
        r = build_phase2_reactor(
            ladder_level="p2",
            builder_kwargs={"enable_vorab_march_staged": True},
            cfg_overrides=overrides,
        )
        out, wall = timed_phase2_solve(
            r, **phase2_solve_kwargs(max_global_iter=int(max_global_iter))
        )
        return _summarize(r, out, wall_s=wall, label=label, T0_seed=T0_seed)
    finally:
        ips._install_temperature_fence_from_vorabrechnung = _orig


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(p, default=8)
    add_json_output_arg(p, "phase2_bed0_multistart_hot_branch_o8.json")
    p.add_argument(
        "--profile",
        type=str,
        default="front_loaded",
        choices=("front_loaded", "linear"),
    )
    p.add_argument(
        "--seeds",
        type=float,
        nargs="+",
        default=[750.0, 950.0, 1100.0, 1250.0, 1400.0, 1600.0],
    )
    p.add_argument("--half-lo", type=float, default=200.0, help="fence lower = seed - half_lo")
    p.add_argument("--half-hi", type=float, default=400.0, help="fence upper ≥ seed + half_hi")
    p.add_argument(
        "--keep-cold-anchor",
        action="store_true",
        help="保留 bed0 cold-anchor（默认关闭以便热种子生效）",
    )
    args = p.parse_args()

    rows: list[dict[str, Any]] = []
    for T0 in args.seeds:
        print(
            f"=== {args.profile} seed T0={T0:.0f} "
            f"(fence≈[{max(700, T0 - args.half_lo):.0f}, "
            f"{min(2500, max(T0 + args.half_hi, 1400)):.0f}]) ===",
            flush=True,
        )
        rows.append(
            run_seed(
                profile=str(args.profile),
                T0_seed=float(T0),
                half_lo=float(args.half_lo),
                half_hi=float(args.half_hi),
                max_global_iter=int(args.max_global_iter),
                disable_cold_anchor=not bool(args.keep_cold_anchor),
            )
        )

    best = _select_best(rows)
    if best is not None:
        print(
            f"SELECTED: {best['label']} ign={best['ignited']} "
            f"conv={best['converged_fully']} rms={best['rms_scaled_final']:.4f} "
            f"T0={best['T0_K']:.0f} CO={best['exit_CO']:.4f}",
            flush=True,
        )

    payload = {
        "policy": (
            "bed0 multi-start T seeds; soft fence; prefer ignited "
            f"(T0>={IGNITED_T_K:.0f}K) then converged then rms then CO"
        ),
        "ignited_T_K": IGNITED_T_K,
        "profile": args.profile,
        "max_global_iter": int(args.max_global_iter),
        "seeds_K": [float(x) for x in args.seeds],
        "half_lo_K": float(args.half_lo),
        "half_hi_K": float(args.half_hi),
        "cold_anchor_disabled": not bool(args.keep_cold_anchor),
        "cases": rows,
        "selected": best,
    }
    write_json(args.output, payload)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
