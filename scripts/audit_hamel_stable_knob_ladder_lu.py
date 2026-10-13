#!/usr/bin/env python3
"""可收敛优先的单因子梯子（强制 R12 关闭）。

相对生产口径，逐项只动一个（或可加性组合）knob，观察：
  - 收敛（converged_fully / rms / n_iter）
  - T0、床温轴向平滑度
  - bed-top / 出口 O₂ 与还原气

Cases
-----
  P0   生产
  S1   仅 r*_scale=1
  H0   仅入口 hydro 分流
  L0   仅关 extent limiter
  S1H  scales=1 + hydro
  S1L  scales=1 + 关 lim
  S1HL scales=1 + hydro + 关 lim

Usage
-----
    python3 scripts/audit_hamel_stable_knob_ladder_lu.py --max-global-iter 8
"""

from __future__ import annotations

from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    apply_baseline_kinetics,
    bed_cell_indices,
    build_phase2_reactor,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)

import argparse
import json
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np

from src.core import cell_kinetics as ck
from src.core.connectivity.gas_inlet import resolve_bottom_gas_inlet_dense_fraction
from src.core.species import GAS_SPECIES_INDEX as idx


@contextmanager
def _extent_limiters(enabled: bool):
    prev = ck._EXTENT_LIMITERS_ENABLED_OVERRIDE
    ck._EXTENT_LIMITERS_ENABLED_OVERRIDE = bool(enabled)
    try:
        yield
    finally:
        ck._EXTENT_LIMITERS_ENABLED_OVERRIDE = prev


def _y(cell: Any) -> dict[str, float]:
    y = cell._mole_fractions("combined")
    return {sp: float(y[j]) for sp, j in idx.items()}


def _axial_smoothness(T: list[float]) -> dict[str, float]:
    arr = np.asarray(T, dtype=float)
    if arr.size < 2:
        return {"T_std_K": float("nan"), "T_max_abs_dT_K": float("nan"), "T_range_K": float("nan")}
    dT = np.diff(arr)
    return {
        "T_std_K": float(np.std(arr)),
        "T_max_abs_dT_K": float(np.max(np.abs(dT))),
        "T_range_K": float(np.max(arr) - np.min(arr)),
    }


def _snapshot(
    reactor: Any,
    result: dict[str, Any],
    *,
    wall_s: float,
    extent_limiters_enabled: bool,
) -> dict[str, Any]:
    beds = bed_cell_indices(reactor)
    axial: list[dict[str, Any]] = []
    for i in beds:
        cell = reactor.cells[i]
        y = _y(cell)
        axial.append(
            {
                "i": i,
                "z_m": float(cell.geo.h_center),
                "T_K": float(cell.T),
                "y_O2": y.get("O2", 0.0),
                "y_CO": y.get("CO", 0.0),
                "y_CO2": y.get("CO2", 0.0),
                "y_H2": y.get("H2", 0.0),
                "y_CH4": y.get("CH4", 0.0),
            }
        )
    T_prof = [r["T_K"] for r in axial]
    smooth = _axial_smoothness(T_prof)
    exit_dry = result.get("exit_gas_dry") or {}
    c0 = reactor.cells[beds[0]]
    dense_frac = float(resolve_bottom_gas_inlet_dense_fraction(c0, reactor.config))
    o2_d = float(c0.N_zu_d[idx["O2"]])
    o2_b = float(c0.N_zu_b[idx["O2"]])
    return {
        "wall_s": wall_s,
        "converged": bool(result.get("converged")),
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "T0_K": float(T_prof[0]) if T_prof else float("nan"),
        "T_bed_mean_K": float(np.mean(T_prof)) if T_prof else float("nan"),
        "T_bed_profile_K": T_prof,
        "bed_top_y_O2": float(axial[-1]["y_O2"]) if axial else float("nan"),
        "exit_dry": {k: float(v) for k, v in exit_dry.items()},
        "inlet": {
            "strategy": str(reactor.config.gas_inlet_split_strategy),
            "cfg_dense_frac": float(reactor.config.gas_inlet_dense_frac),
            "resolved_dense_frac": dense_frac,
            "N_zu_O2_dense": o2_d,
            "N_zu_O2_bubble": o2_b,
            "O2_dense_share": o2_d / max(o2_d + o2_b, 1e-30),
        },
        "smoothness": smooth,
        "axial_bed": axial,
        "cfg": {
            "enable_r12": bool(reactor.config.enable_r12),
            "r4_scale": float(reactor.config.r4_scale),
            "r5_scale": float(reactor.config.r5_scale),
            "r7_scale": float(reactor.config.r7_scale),
            "extent_limiters_enabled": bool(extent_limiters_enabled),
        },
    }


def _push_cells(reactor: Any, *, scales_one: bool) -> None:
    # 本梯子强制 R12 关闭
    reactor.config.enable_r12 = False
    if scales_one:
        apply_baseline_kinetics(reactor.config)
    for cell in list(reactor.cells) + list(getattr(reactor, "freeboard_cells", []) or []):
        cell.enable_r12 = False
        if scales_one:
            for s in ("r2", "r4", "r5", "r6", "r7", "r8"):
                setattr(cell, f"{s}_scale", 1.0)
            cell.r5_bubble_scale = 1.0


def _run_case(
    name: str,
    *,
    max_global_iter: int,
    scales_one: bool = False,
    limiters_on: bool = True,
    hydro_inlet: bool = False,
) -> dict[str, Any]:
    overrides: dict[str, object] = {"enable_r12": False}
    if hydro_inlet:
        overrides.update(
            {
                "gas_inlet_split_strategy": "precalc_hydrodynamic_flux",
                "gas_inlet_dense_frac": 0.20,
                "r5_bubble_scale_bed0": 1.0,
            }
        )
    reactor = build_phase2_reactor(
        ladder_level="p2",
        cfg_overrides=overrides,
        apply_kinetics=scales_one,
    )
    _push_cells(reactor, scales_one=scales_one)
    run_phase2_init(reactor)
    with _extent_limiters(enabled=limiters_on):
        result, wall_s = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=max_global_iter)
        )
    snap = _snapshot(
        reactor, result, wall_s=wall_s, extent_limiters_enabled=limiters_on
    )
    snap["name"] = name
    snap["knobs"] = {
        "enable_r12": False,
        "scales_one": scales_one,
        "limiters_on": limiters_on,
        "hydro_inlet": hydro_inlet,
    }
    return snap


CASES: list[tuple[str, dict[str, Any]]] = [
    # P0：Hamel 同构生产默认（裸动力学，limiters off）
    ("P0_production", dict(scales_one=False, limiters_on=False, hydro_inlet=False)),
    # 工程对照：旧生产限幅路径
    ("P0_eng_limiters", dict(scales_one=False, limiters_on=True, hydro_inlet=False)),
    ("S1_scales1", dict(scales_one=True, limiters_on=True, hydro_inlet=False)),
    ("H0_hydro", dict(scales_one=False, limiters_on=True, hydro_inlet=True)),
    ("L0_nolim", dict(scales_one=False, limiters_on=False, hydro_inlet=False)),
    ("S1H_scales1_hydro", dict(scales_one=True, limiters_on=True, hydro_inlet=True)),
    ("S1L_scales1_nolim", dict(scales_one=True, limiters_on=False, hydro_inlet=False)),
    ("S1HL_scales1_hydro_nolim", dict(scales_one=True, limiters_on=False, hydro_inlet=True)),
]


def _print(c: dict[str, Any]) -> None:
    sm = c["smoothness"]
    print(
        f"  conv={c['converged_fully']} rms={c['rms_scaled_final']:.3e} n={c['n_iter']} "
        f"T0={c['T0_K']:.1f} dTmax={sm['T_max_abs_dT_K']:.1f} Trange={sm['T_range_K']:.1f} "
        f"bedO2={c['bed_top_y_O2']:.4f} "
        f"exO2={c['exit_dry'].get('O2', float('nan')):.4f} "
        f"exCO={c['exit_dry'].get('CO', float('nan')):.4f} "
        f"exH2={c['exit_dry'].get('H2', float('nan')):.4f} "
        f"f_dense={c['inlet']['O2_dense_share']:.3f} wall={c['wall_s']:.1f}s",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(parser, default=8)
    add_json_output_arg(parser, "phase2_hamel_stable_knob_ladder_o8.json")
    parser.add_argument("--only", nargs="+", default=None)
    parser.add_argument("--reuse", type=Path, default=None)
    args = parser.parse_args()

    reuse: dict[str, Any] = {}
    if args.reuse and args.reuse.is_file():
        prev = json.loads(args.reuse.read_text(encoding="utf-8"))
        for c in prev.get("cases", []):
            reuse[str(c["name"])] = c

    wanted = list(args.only) if args.only else [n for n, _ in CASES]
    step_map = dict(CASES)
    cases: list[dict[str, Any]] = []
    for name in wanted:
        if name not in step_map:
            raise SystemExit(f"unknown {name}; choose from {list(step_map)}")
        if name in reuse:
            print(f"=== {name} (reused) ===", flush=True)
            cases.append(reuse[name])
            _print(reuse[name])
            continue
        print(f"=== {name} ===", flush=True)
        c = _run_case(name, max_global_iter=args.max_global_iter, **step_map[name])
        cases.append(c)
        _print(c)

    by = {c["name"]: c for c in cases}
    base = by.get("P0_production")
    ranking = []
    if base is not None:
        for n, c in by.items():
            if n == "P0_production":
                continue
            # 收敛优先分：未满收敛重罚；平滑差重罚；再看点火/耗氧
            conv_pen = 0.0 if c["converged_fully"] else 200.0
            # rms 对数惩罚（生产也未必 fully，用相对）
            rms = max(float(c["rms_scaled_final"]), 1e-12)
            rms0 = max(float(base["rms_scaled_final"]), 1e-12)
            rms_pen = 20.0 * max(np.log10(rms / rms0), 0.0)
            smooth_pen = 0.15 * float(c["smoothness"]["T_max_abs_dT_K"])
            dT0 = float(c["T0_K"] - base["T0_K"])
            dO2 = float(c["bed_top_y_O2"] - base["bed_top_y_O2"])
            dCO = float(c["exit_dry"].get("CO", 0.0) - base["exit_dry"].get("CO", 0.0))
            score = dT0 - 400.0 * dO2 + 50.0 * dCO - conv_pen - rms_pen - smooth_pen
            ranking.append(
                {
                    "name": n,
                    "score": score,
                    "dT0_K": dT0,
                    "d_bed_top_y_O2": dO2,
                    "d_exit_CO": dCO,
                    "d_exit_H2": float(
                        c["exit_dry"].get("H2", 0.0) - base["exit_dry"].get("H2", 0.0)
                    ),
                    "converged_fully": c["converged_fully"],
                    "rms_scaled_final": c["rms_scaled_final"],
                    "T_max_abs_dT_K": c["smoothness"]["T_max_abs_dT_K"],
                }
            )
        ranking.sort(key=lambda r: r["score"], reverse=True)

    summary = {
        "policy": "R12 always off; single-factor / additive knobs; prefer convergence",
        "max_global_iter": args.max_global_iter,
        "cases": cases,
        "ranking_vs_P0": ranking,
    }
    if ranking:
        print("=== ranking vs P0 (higher better; penalize nonconv/rms/jagged T) ===", flush=True)
        for r in ranking:
            print(
                f"  {r['name']}: score={r['score']:.1f} "
                f"conv={r['converged_fully']} rms={r['rms_scaled_final']:.2e} "
                f"dT0={r['dT0_K']:+.1f} d_bedO2={r['d_bed_top_y_O2']:+.4f} "
                f"d_exCO={r['d_exit_CO']:+.4f} dTmax={r['T_max_abs_dT_K']:.1f}",
                flush=True,
            )
    write_json(args.output, summary)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
