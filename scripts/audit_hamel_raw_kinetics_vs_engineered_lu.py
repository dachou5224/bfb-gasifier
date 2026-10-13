#!/usr/bin/env python3
"""Hamel 原版动力学 knobs 分层 A/B：温度与 O₂。

从生产口径逐项打开：
  A  ``r*_scale=1``（仍关 R12、仍开 limiter、仍 fixed 入口）
  B  A + 关 extent limiter
  C  A + ``enable_r12=True``（仍开 limiter）
  D  A + 关 limiter + R12
  E  D + 入口改流体力学闭包（= 完整 hamel_raw）

Usage
-----
    python3 scripts/audit_hamel_raw_kinetics_vs_engineered_lu.py --ladder
    python3 scripts/audit_hamel_raw_kinetics_vs_engineered_lu.py --ladder --max-global-iter 8
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
from src.core.species import GAS_SPECIES_INDEX as idx


@contextmanager
def _extent_limiters(enabled: bool):
    prev = ck._EXTENT_LIMITERS_ENABLED_OVERRIDE
    ck._EXTENT_LIMITERS_ENABLED_OVERRIDE = bool(enabled)
    try:
        yield
    finally:
        ck._EXTENT_LIMITERS_ENABLED_OVERRIDE = prev


def _y_wet(cell: Any) -> dict[str, float]:
    y = cell._mole_fractions("combined")
    return {sp: float(y[j]) for sp, j in idx.items()}


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
        y = _y_wet(cell)
        N = np.asarray(cell.N_b, dtype=float) + np.asarray(cell.N_d, dtype=float)
        axial.append(
            {
                "i": i,
                "cell_type": "bed",
                "T_K": float(cell.T),
                "y_O2": y.get("O2", 0.0),
                "y_CO": y.get("CO", 0.0),
                "y_CO2": y.get("CO2", 0.0),
                "y_H2": y.get("H2", 0.0),
                "y_CH4": y.get("CH4", 0.0),
                "N_O2_mol": float(N[idx["O2"]]),
            }
        )
    for j, cell in enumerate(getattr(reactor, "freeboard_cells", []) or []):
        y = _y_wet(cell)
        axial.append(
            {
                "i": j,
                "cell_type": "freeboard",
                "T_K": float(cell.T),
                "y_O2": y.get("O2", 0.0),
                "y_CO": y.get("CO", 0.0),
                "y_CO2": y.get("CO2", 0.0),
                "y_H2": y.get("H2", 0.0),
                "y_CH4": y.get("CH4", 0.0),
            }
        )

    T_prof = list(result.get("bed_T_profile") or [reactor.cells[i].T for i in beds])
    exit_dry = result.get("exit_gas_dry") or {}
    bed_top = axial[beds[-1]] if beds else {}
    T0 = float(T_prof[0]) if T_prof else float("nan")
    return {
        "wall_s": wall_s,
        "converged": bool(result.get("converged")),
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "T0_K": T0,
        "T_bed_mean_K": float(np.mean(T_prof)) if T_prof else float("nan"),
        "T_bed_max_K": float(np.max(T_prof)) if T_prof else float("nan"),
        "T_bed_profile_K": [float(t) for t in T_prof],
        "bed_top_y_O2": float(bed_top.get("y_O2", float("nan"))),
        "exit_dry": {k: float(v) for k, v in exit_dry.items()} if exit_dry else {},
        "O2_feed_mol_s": float(getattr(reactor.config, "O2_feed", 0.0)),
        "axial": axial,
        "cfg": {
            "enable_r12": bool(reactor.config.enable_r12),
            "r2_scale": float(reactor.config.r2_scale),
            "r4_scale": float(reactor.config.r4_scale),
            "r5_scale": float(reactor.config.r5_scale),
            "r6_scale": float(reactor.config.r6_scale),
            "r7_scale": float(reactor.config.r7_scale),
            "r8_scale": float(reactor.config.r8_scale),
            "gas_inlet_split_strategy": str(reactor.config.gas_inlet_split_strategy),
            "gas_inlet_dense_frac": float(reactor.config.gas_inlet_dense_frac),
            "extent_limiters_enabled": bool(extent_limiters_enabled),
            "cold_anchor": bool(
                getattr(reactor.config, "vorab_nr_temperature_seed_bed0_cold_anchor_thesis", False)
            ),
            "seed_bottom_K": float(
                getattr(reactor.config, "vorab_nr_temperature_seed_bottom_K_thesis", float("nan"))
            ),
        },
    }


def _push_scales(reactor: Any, *, scales_one: bool, enable_r12: bool) -> None:
    if scales_one:
        apply_baseline_kinetics(reactor.config)
    reactor.config.enable_r12 = bool(enable_r12)
    for cell in list(reactor.cells) + list(getattr(reactor, "freeboard_cells", []) or []):
        cell.enable_r12 = bool(enable_r12)
        if scales_one:
            for s in ("r2", "r4", "r5", "r6", "r7", "r8"):
                setattr(cell, f"{s}_scale", 1.0)
            cell.r5_bubble_scale = 1.0


def _run_case(
    name: str,
    *,
    max_global_iter: int,
    scales_one: bool = False,
    enable_r12: bool = False,
    limiters_on: bool = True,
    hydro_inlet: bool = False,
) -> dict[str, Any]:
    overrides: dict[str, object] = {"enable_r12": bool(enable_r12)}
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
    _push_scales(reactor, scales_one=scales_one, enable_r12=enable_r12)

    run_phase2_init(reactor)
    solve_kw = phase2_solve_kwargs(max_global_iter=max_global_iter)
    with _extent_limiters(enabled=limiters_on):
        result, wall_s = timed_phase2_solve(reactor, **solve_kw)
    snap = _snapshot(
        reactor, result, wall_s=wall_s, extent_limiters_enabled=limiters_on
    )
    snap["name"] = name
    snap["knobs"] = {
        "scales_one": scales_one,
        "enable_r12": enable_r12,
        "limiters_on": limiters_on,
        "hydro_inlet": hydro_inlet,
    }
    return snap


def _print_case(c: dict[str, Any]) -> None:
    print(
        f"  T0={c['T0_K']:.1f} Tmean={c['T_bed_mean_K']:.1f} "
        f"bed_top_yO2={c['bed_top_y_O2']:.4f} "
        f"exit_O2={c['exit_dry'].get('O2', float('nan')):.5f} "
        f"exit_CO={c['exit_dry'].get('CO', float('nan')):.5f} "
        f"exit_H2={c['exit_dry'].get('H2', float('nan')):.5f} "
        f"conv={c['converged_fully']} wall={c['wall_s']:.1f}s",
        flush=True,
    )


def _metrics(c: dict[str, Any]) -> dict[str, float]:
    return {
        "T0_K": float(c["T0_K"]),
        "T_bed_mean_K": float(c["T_bed_mean_K"]),
        "bed_top_y_O2": float(c["bed_top_y_O2"]),
        "exit_O2": float(c["exit_dry"].get("O2", 0.0)),
        "exit_CO": float(c["exit_dry"].get("CO", 0.0)),
        "exit_H2": float(c["exit_dry"].get("H2", 0.0)),
        "exit_CH4": float(c["exit_dry"].get("CH4", 0.0)),
    }


def _delta(a: dict[str, Any], b: dict[str, Any]) -> dict[str, float]:
    ma, mb = _metrics(a), _metrics(b)
    return {f"d_{k}": mb[k] - ma[k] for k in ma}


LADDER_STEPS: list[tuple[str, dict[str, Any]]] = [
    ("production", dict(scales_one=False, enable_r12=False, limiters_on=True, hydro_inlet=False)),
    ("A_scales1", dict(scales_one=True, enable_r12=False, limiters_on=True, hydro_inlet=False)),
    ("B_scales1_nolim", dict(scales_one=True, enable_r12=False, limiters_on=False, hydro_inlet=False)),
    ("C_scales1_r12", dict(scales_one=True, enable_r12=True, limiters_on=True, hydro_inlet=False)),
    ("D_scales1_nolim_r12", dict(scales_one=True, enable_r12=True, limiters_on=False, hydro_inlet=False)),
    ("E_hamel_raw", dict(scales_one=True, enable_r12=True, limiters_on=False, hydro_inlet=True)),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(parser, default=8)
    add_json_output_arg(parser, "phase2_hamel_kinetics_knob_ladder_o8.json")
    parser.add_argument("--ladder", action="store_true", help="跑完整分层梯子")
    parser.add_argument(
        "--reuse",
        type=Path,
        default=None,
        help="复用已有 JSON 中同名 case（跳过重跑）",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        default=None,
        help="只跑指定 case 名，例如 A_scales1 C_scales1_r12",
    )
    args = parser.parse_args()

    reuse: dict[str, Any] = {}
    if args.reuse and args.reuse.is_file():
        prev = json.loads(args.reuse.read_text(encoding="utf-8"))
        for c in prev.get("cases", []):
            reuse[str(c["name"])] = c
        print(f"reuse {len(reuse)} cases from {args.reuse}", flush=True)

    wanted = list(args.only) if args.only else [n for n, _ in LADDER_STEPS]
    if not args.ladder and args.only is None:
        # 兼容旧入口：production + E
        wanted = ["production", "E_hamel_raw"]

    cases: list[dict[str, Any]] = []
    step_map = dict(LADDER_STEPS)
    for name in wanted:
        if name not in step_map:
            raise SystemExit(f"unknown case {name!r}; choose from {list(step_map)}")
        if name in reuse:
            print(f"=== {name} (reused) ===", flush=True)
            cases.append(reuse[name])
            _print_case(reuse[name])
            continue
        print(f"=== {name} ===", flush=True)
        c = _run_case(name, max_global_iter=args.max_global_iter, **step_map[name])
        cases.append(c)
        _print_case(c)

    by_name = {c["name"]: c for c in cases}
    summary: dict[str, Any] = {"cases": cases, "metrics": {n: _metrics(c) for n, c in by_name.items()}}

    # 相邻增量（若梯子齐全）
    order = [n for n, _ in LADDER_STEPS if n in by_name]
    increments: dict[str, Any] = {}
    for a, b in zip(order[:-1], order[1:]):
        key = f"{a}→{b}"
        increments[key] = _delta(by_name[a], by_name[b])
    summary["increments"] = increments
    if increments:
        print("=== increments ===", flush=True)
        for key, d in increments.items():
            print(
                f"  {key}: dT0={d['d_T0_K']:+.1f} d_bedO2={d['d_bed_top_y_O2']:+.4f} "
                f"d_exO2={d['d_exit_O2']:+.4f} d_exCO={d['d_exit_CO']:+.4f} d_exH2={d['d_exit_H2']:+.5f}",
                flush=True,
            )

    if "production" in by_name and "E_hamel_raw" in by_name:
        summary["delta_total"] = _delta(by_name["production"], by_name["E_hamel_raw"])

    # 因果排序：相对 production 的 |dT0| + 耗氧改善
    if "production" in by_name:
        base = by_name["production"]
        ranking = []
        for n, c in by_name.items():
            if n == "production":
                continue
            d = _delta(base, c)
            ranking.append(
                {
                    "name": n,
                    "dT0_K": d["d_T0_K"],
                    "d_bed_top_y_O2": d["d_bed_top_y_O2"],
                    "d_exit_O2": d["d_exit_O2"],
                    "d_exit_CO": d["d_exit_CO"],
                    "d_exit_H2": d["d_exit_H2"],
                    "score_ignition": d["d_T0_K"] - 500.0 * d["d_bed_top_y_O2"],
                }
            )
        ranking.sort(key=lambda r: r["score_ignition"], reverse=True)
        summary["ranking_vs_production"] = ranking
        print("=== ranking vs production (higher = hotter + less bed O2) ===", flush=True)
        for r in ranking:
            print(
                f"  {r['name']}: score={r['score_ignition']:.1f} "
                f"dT0={r['dT0_K']:+.1f} d_bedO2={r['d_bed_top_y_O2']:+.4f}",
                flush=True,
            )

    write_json(args.output, summary)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
