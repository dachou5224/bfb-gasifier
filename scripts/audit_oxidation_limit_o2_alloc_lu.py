#!/usr/bin/env python3
"""氧化限速 vs O₂ 工程分配：是否限制点火/收敛（zone=3 生产基线）。

单因子：
  P0     生产（r5=0.25, lim ON, dense=0.5, staged ON）
  R5     仅 r5_scale=1
  L0     仅关 extent limiter
  D20    dense_frac=0.20（更多 O₂ 进气泡相）
  D80    dense_frac=0.80（更多 O₂ 进密相）
  ST0    关闭 staged primary（全部一次风进 bed0）
  R5L    r5=1 + 关 lim
  D80L   dense=0.80 + 关 lim

Usage
-----
    python3 scripts/audit_oxidation_limit_o2_alloc_lu.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from typing import Any

import numpy as np

from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    bed_cell_indices,
    build_phase2_reactor,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
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


def _probe_bed0_o2_limit_factors(reactor: Any) -> dict[str, float]:
    """在当前状态算一次残差，读 bed0 氧化 extent 限幅诊断。"""
    beds = bed_cell_indices(reactor)
    c0 = reactor.cells[beds[0]]
    c0.residuals()
    # cell_kinetics 把限幅写在最近一次 kinetics 诊断属性上（若存在）
    keys = [
        "limit_factor_o2",
        "limit_factor_o2_dense",
        "limit_factor_o2_bubble",
        "o2_budget_dense",
        "o2_budget_bubble",
        "dense_o2_demand",
        "bubble_o2_demand",
        "r5d",
        "r5d_lim",
        "r6d",
        "r6d_lim",
    ]
    out: dict[str, float] = {}
    diag = getattr(c0, "_last_kinetics_diag", None) or getattr(c0, "_kinetics_diag", None)
    if isinstance(diag, dict):
        for k in keys:
            if k in diag:
                out[k] = float(diag[k])
    # 回退：从模块级最近一次（若 cell 无 diag）
    for k in keys:
        if k not in out and hasattr(ck, f"_LAST_{k.upper()}"):
            out[k] = float(getattr(ck, f"_LAST_{k.upper()}"))
    # 直接从 holdup/进料估 O₂ 密相份额
    o2_d = float(c0.N_zu_d[idx["O2"]])
    o2_b = float(c0.N_zu_b[idx["O2"]])
    out["N_zu_O2_dense"] = o2_d
    out["N_zu_O2_bubble"] = o2_b
    out["O2_dense_share"] = o2_d / max(o2_d + o2_b, 1e-30)
    out["resolved_dense_frac"] = float(
        resolve_bottom_gas_inlet_dense_fraction(c0, reactor.config)
    )
    out["T0_K"] = float(c0.T)
    out["y_O2"] = float(_y(c0).get("O2", 0.0))
    # 若 kinetics 把 lf 挂在 cell
    for attr in (
        "limit_factor_o2",
        "limit_factor_o2_dense",
        "limit_factor_o2_bubble",
        "_limit_factor_o2",
        "_limit_factor_o2_dense",
    ):
        if hasattr(c0, attr):
            try:
                out[attr.lstrip("_")] = float(getattr(c0, attr))
            except Exception:
                pass
    return out


def _summarize(
    reactor: Any,
    result: dict[str, Any],
    *,
    wall_s: float,
    name: str,
    knobs: dict[str, Any],
    init_probe: dict[str, float],
    final_probe: dict[str, float],
) -> dict[str, Any]:
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
                "y_H2": y.get("H2", 0.0),
                "y_CO2": y.get("CO2", 0.0),
            }
        )
    T = [r["T_K"] for r in axial]
    exit_dry = {k: float(v) for k, v in (result.get("exit_gas_dry") or {}).items()}
    row = {
        "name": name,
        "knobs": knobs,
        "wall_s": wall_s,
        "converged_fully": bool(result.get("converged_fully")),
        "n_iter": result.get("n_iter"),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "T0_K": float(T[0]) if T else float("nan"),
        "T_bed_max_K": float(np.max(T)) if T else float("nan"),
        "T_bed_profile_K": T,
        "bed_top_y_O2": float(axial[-1]["y_O2"]) if axial else float("nan"),
        "exit_dry": exit_dry,
        "vm_zone": int(getattr(reactor.config, "vorab_vm_devolatilization_zone_cells_thesis", 0)),
        "init_o2_probe": init_probe,
        "final_o2_probe": final_probe,
        "axial_bed": axial,
    }
    hot = row["T0_K"] >= 1100.0 and row["bed_top_y_O2"] < 1e-3
    print(
        f"{name}: hot={hot} conv={row['converged_fully']} rms={row['rms_scaled_final']:.3e} "
        f"T0={row['T0_K']:.1f} Tmax={row['T_bed_max_K']:.1f} bedO2={row['bed_top_y_O2']:.4f} "
        f"exCO/H2/O2={exit_dry.get('CO', float('nan')):.4f}/"
        f"{exit_dry.get('H2', float('nan')):.4f}/{exit_dry.get('O2', float('nan')):.3e} "
        f"init_lfO2={init_probe.get('limit_factor_o2', float('nan'))} "
        f"dense={init_probe.get('O2_dense_share', float('nan')):.2f} "
        f"wall={wall_s:.1f}s",
        flush=True,
    )
    return row


def _run(
    name: str,
    *,
    max_global_iter: int,
    r5_scale: float | None = None,
    limiters_on: bool = True,
    dense_frac: float | None = None,
    staged_primary: bool | None = None,
) -> dict[str, Any]:
    overrides: dict[str, object] = {"enable_r12": False}
    if r5_scale is not None:
        overrides["r5_scale"] = float(r5_scale)
    if dense_frac is not None:
        overrides["gas_inlet_split_strategy"] = "fixed"
        overrides["gas_inlet_dense_frac"] = float(dense_frac)
        overrides["vorab_staged_primary_gas_dense_frac_thesis"] = float(dense_frac)
    if staged_primary is not None:
        overrides["vorab_staged_primary_gas_injection_thesis"] = bool(staged_primary)
        if not staged_primary:
            overrides["vorab_a_tier_bed_axial_primary_profile_thesis"] = False

    reactor = build_phase2_reactor(
        ladder_level="p2",
        builder_kwargs={"enable_vorab_march_staged": True},
        cfg_overrides=overrides,
    )
    # 同步到 cell
    for cell in list(reactor.cells) + list(getattr(reactor, "freeboard_cells", []) or []):
        cell.enable_r12 = False
        if r5_scale is not None:
            cell.r5_scale = float(r5_scale)

    run_phase2_init(reactor)
    with _extent_limiters(enabled=limiters_on):
        init_probe = _probe_bed0_o2_limit_factors(reactor)
        result, wall_s = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=int(max_global_iter))
        )
        final_probe = _probe_bed0_o2_limit_factors(reactor)

    knobs = {
        "r5_scale": float(reactor.config.r5_scale),
        "limiters_on": bool(limiters_on),
        "dense_frac": float(reactor.config.gas_inlet_dense_frac),
        "staged_primary": bool(
            getattr(reactor.config, "vorab_staged_primary_gas_injection_thesis", False)
        ),
        "vm_zone": int(getattr(reactor.config, "vorab_vm_devolatilization_zone_cells_thesis", 0)),
    }
    return _summarize(
        reactor,
        result,
        wall_s=wall_s,
        name=name,
        knobs=knobs,
        init_probe=init_probe,
        final_probe=final_probe,
    )


CASES: list[tuple[str, dict[str, Any]]] = [
    ("P0_production", dict()),
    ("R5_scale1", dict(r5_scale=1.0)),
    ("L0_nolim", dict(limiters_on=False)),
    ("D20_more_bubble_O2", dict(dense_frac=0.20)),
    ("D80_more_dense_O2", dict(dense_frac=0.80)),
    ("ST0_no_staged", dict(staged_primary=False)),
    ("R5L_scale1_nolim", dict(r5_scale=1.0, limiters_on=False)),
]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(p, default=8)
    add_json_output_arg(p, "phase2_oxidation_limit_o2_alloc_o8.json")
    p.add_argument("--only", nargs="+", default=None)
    args = p.parse_args()
    wanted = list(args.only) if args.only else [n for n, _ in CASES]
    step = dict(CASES)
    rows = []
    for name in wanted:
        if name not in step:
            raise SystemExit(f"unknown {name}; choose from {list(step)}")
        print(f"=== {name} ===", flush=True)
        rows.append(_run(name, max_global_iter=args.max_global_iter, **step[name]))

    base = next((r for r in rows if r["name"] == "P0_production"), rows[0] if rows else None)
    deltas = []
    if base is not None:
        for r in rows:
            if r is base:
                continue
            ea, eb = r["exit_dry"], base["exit_dry"]
            deltas.append(
                {
                    "name": r["name"],
                    "d_T0_K": r["T0_K"] - base["T0_K"],
                    "d_bed_top_y_O2": r["bed_top_y_O2"] - base["bed_top_y_O2"],
                    "d_rms": r["rms_scaled_final"] - base["rms_scaled_final"],
                    "d_exit_CO": ea.get("CO", 0.0) - eb.get("CO", 0.0),
                    "d_exit_H2": ea.get("H2", 0.0) - eb.get("H2", 0.0),
                    "d_exit_O2": ea.get("O2", 0.0) - eb.get("O2", 0.0),
                    "jumped_hot": bool(r["T0_K"] >= 1100.0 and r["bed_top_y_O2"] < 1e-3),
                }
            )
            print(
                f"  delta vs P0 {r['name']}: dT0={deltas[-1]['d_T0_K']:+.1f} "
                f"dBedO2={deltas[-1]['d_bed_top_y_O2']:+.4f} "
                f"dCO={deltas[-1]['d_exit_CO']:+.4f} hot={deltas[-1]['jumped_hot']}",
                flush=True,
            )

    payload = {
        "question": "Do oxidation rate limits or engineered O2 allocation limit ignition/convergence?",
        "baseline": "Phase2 p2, vm_zone=3, R12 off, cold-anchor",
        "max_global_iter": int(args.max_global_iter),
        "cases": rows,
        "deltas_vs_P0": deltas,
    }
    out = write_json(args.output, payload)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
