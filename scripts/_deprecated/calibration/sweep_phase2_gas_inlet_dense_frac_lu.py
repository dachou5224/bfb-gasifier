#!/usr/bin/env python3
"""Phase2 LU：床底 O₂ 密相/气泡 split 敏感性。

扫描 ``gas_inlet_split_strategy`` / ``gas_inlet_dense_frac``，报告出口与 bed0–2
T、CO/CO₂、R5b/R5d 分反应预算。

Usage
-----
    PYTHONUNBUFFERED=1 python3 scripts/sweep_phase2_gas_inlet_dense_frac_lu.py
    PYTHONUNBUFFERED=1 python3 scripts/sweep_phase2_gas_inlet_dense_frac_lu.py \\
        --dense-fracs 0.35 0.50 0.65 0.80 --max-global-iter 20 --kinetics phase2
"""

from __future__ import annotations

from _lu_harness import (
    REF_CO,
    REF_CO2,
    add_max_global_iter_arg,
    build_phase2_config,
    data_path,
    phase2_solve_kwargs,
    repo_root,
    write_json,
)

REPO_ROOT = repo_root()

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

from audit_bed_co_co2_budget_lu import _analyze_cell_co_co2
from src.core.connectivity import resolve_bottom_gas_inlet_dense_fraction
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx


@dataclass
class BedCellRow:
    cell_index: int
    xi_over_H: float
    T_K: float
    y_CO: float
    y_CO2: float
    y_O2: float
    R5b_CO: float
    R5d_CO: float
    R1_CO2: float
    co_net: float
    co2_net: float


@dataclass
class Row:
    label: str
    gas_inlet_split_strategy: str
    gas_inlet_dense_frac_cfg: float
    dense_frac_resolved_bed0: float
    staged_dense_frac_cfg: float
    r2_scale: float
    r4_scale: float
    r5_scale: float
    r8_scale: float
    converged_fully: bool
    rms_scaled_final: float
    co: float
    co2: float
    bed_co: float
    bed_co2: float
    h2: float
    ch4: float
    t_exit_K: float
    score_co_co2: float
    bed0_N_zu_d_O2_frac: float
    bed0_N_zu_b_O2_frac: float
    bed0_y_O2_bubble: float
    bed0_y_O2_dense: float
    bed_cells: list[dict[str, float | int]]
    wall_time_s: float


def _score_co_co2(co: float, co2: float) -> float:
    return abs(co - REF_CO) / REF_CO + abs(co2 - REF_CO2) / REF_CO2


def _apply_kinetics(cfg, mode: str) -> None:
    """phase2：保留 validation 标定倍率；baseline：全部 1.0。"""
    if str(mode).strip().lower() == "baseline":
        cfg.r2_scale = 1.0
        cfg.r4_scale = 1.0
        cfg.r5_scale = 1.0
        cfg.r6_scale = 1.0
        cfg.r7_scale = 1.0
        cfg.r8_scale = 1.0


def _bed_cell_rows(reactor: Reactor, *, n_cells: int = 3) -> list[BedCellRow]:
    H_bed = float(reactor.config.H_bed)
    out: list[BedCellRow] = []
    for i in range(min(n_cells, len(reactor.cells))):
        cell = reactor.cells[i]
        row = _analyze_cell_co_co2(cell)
        co, co2 = row["co_terms"], row["co2_terms"]
        out.append(
            BedCellRow(
                cell_index=i,
                xi_over_H=float(cell.geo.h_center / max(H_bed, 1e-12)),
                T_K=float(row["T_K"]),
                y_CO=float(row["wet_y"]["CO"]),
                y_CO2=float(row["wet_y"]["CO2"]),
                y_O2=float(row["wet_y"]["O2"]),
                R5b_CO=float(co["R5b"]),
                R5d_CO=float(co["R5d"]),
                R1_CO2=float(co2["R1"]),
                co_net=float(row["co_net"]),
                co2_net=float(row["co2_net"]),
            )
        )
    return out


def run_one(
    *,
    label: str,
    split_strategy: str,
    dense_frac_cfg: float,
    max_global_iter: int,
    kinetics: str,
    sync_staged_dense: bool,
) -> Row:
    cfg = build_phase2_config()
    cfg.vorab_major_gibbs_x0 = False
    _apply_kinetics(cfg, kinetics)
    cfg.gas_inlet_split_strategy = str(split_strategy)
    cfg.gas_inlet_dense_frac = float(dense_frac_cfg)
    if sync_staged_dense and split_strategy in {"fixed", "legacy_fixed"}:
        cfg.vorab_staged_primary_gas_dense_frac_thesis = float(dense_frac_cfg)

    solve_kwargs = phase2_solve_kwargs()
    solve_kwargs["max_global_iter"] = int(max_global_iter)
    solve_kwargs["verbose"] = False

    t0 = perf_counter()
    reactor = Reactor(cfg)
    out = reactor.solve(**solve_kwargs)
    wall_s = perf_counter() - t0

    bed0 = reactor.cells[0]
    o2_feed = float(cfg.O2_feed)
    n_zu_d = float(bed0.N_zu_d[idx["O2"]])
    n_zu_b = float(bed0.N_zu_b[idx["O2"]])
    dense_resolved = float(resolve_bottom_gas_inlet_dense_fraction(bed0, cfg))
    y_b = bed0._mole_fractions("b")
    y_d = bed0._mole_fractions("d")

    y = out.get("exit_gas_dry") or {}
    bed = out.get("bed_exit_gas_dry") or {}
    co = float(y.get("CO", 0.0))
    co2 = float(y.get("CO2", 0.0))
    return Row(
        label=label,
        gas_inlet_split_strategy=str(split_strategy),
        gas_inlet_dense_frac_cfg=float(dense_frac_cfg),
        dense_frac_resolved_bed0=dense_resolved,
        staged_dense_frac_cfg=float(getattr(cfg, "vorab_staged_primary_gas_dense_frac_thesis", 0.70)),
        r2_scale=float(cfg.r2_scale),
        r4_scale=float(cfg.r4_scale),
        r5_scale=float(cfg.r5_scale),
        r8_scale=float(cfg.r8_scale),
        converged_fully=bool(out.get("converged_fully", False)),
        rms_scaled_final=float(out.get("rms_scaled_final", 0.0) or 0.0),
        co=co,
        co2=co2,
        bed_co=float(bed.get("CO", 0.0)),
        bed_co2=float(bed.get("CO2", 0.0)),
        h2=float(y.get("H2", 0.0)),
        ch4=float(y.get("CH4", 0.0)),
        t_exit_K=float(out["T_profile"][-1]),
        score_co_co2=_score_co_co2(co, co2),
        bed0_N_zu_d_O2_frac=n_zu_d / max(o2_feed, 1e-12),
        bed0_N_zu_b_O2_frac=n_zu_b / max(o2_feed, 1e-12),
        bed0_y_O2_bubble=float(y_b[idx["O2"]]),
        bed0_y_O2_dense=float(y_d[idx["O2"]]),
        bed_cells=[asdict(b) for b in _bed_cell_rows(reactor, n_cells=3)],
        wall_time_s=float(wall_s),
    )


def _print_row(row: Row) -> None:
    conv = "OK" if row.converged_fully else "FAIL"
    print(
        f"[{row.label}] strat={row.gas_inlet_split_strategy} cfg_dense={row.gas_inlet_dense_frac_cfg:.2f} "
        f"resolved={row.dense_frac_resolved_bed0:.3f} staged={row.staged_dense_frac_cfg:.2f} | "
        f"conv={conv} rms={row.rms_scaled_final:.4f} | "
        f"zu_O2 d/b={row.bed0_N_zu_d_O2_frac:.2f}/{row.bed0_N_zu_b_O2_frac:.2f} | "
        f"CO={row.co:.4f} CO2={row.co2:.4f} bedCO={row.bed_co:.4f} bedCO2={row.bed_co2:.4f} | "
        f"score={row.score_co_co2:.3f} ({row.wall_time_s:.0f}s)",
        flush=True,
    )
    for bc in row.bed_cells:
        print(
            f"    cell{bc['cell_index']} xi={bc['xi_over_H']:.2f} T={bc['T_K']:.0f}K "
            f"yCO={bc['y_CO']:.3f} yCO2={bc['y_CO2']:.3f} yO2={bc['y_O2']:.4f} "
            f"R5b={bc['R5b_CO']:.2f} R5d={bc['R5d_CO']:.3f} co_net={bc['co_net']:.2f}",
            flush=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dense-fracs",
        type=float,
        nargs="+",
        default=[0.35, 0.50, 0.65, 0.80],
        help="fixed 策略下的 gas_inlet_dense_frac 扫描点",
    )
    parser.add_argument(
        "--kinetics",
        choices=("phase2", "baseline"),
        default="phase2",
        help="phase2=保留 validation 标定 r4/r5/r7；baseline=全部 1.0",
    )
    parser.add_argument(
        "--sync-staged-dense",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="fixed 扫描时同步 vorab_staged_primary_gas_dense_frac_thesis",
    )
    parser.add_argument(
        "--include-precalc",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否包含 precalc_hydrodynamic_flux 基线（当前 Phase2 默认）",
    )
    parser.add_argument("--max-global-iter", type=int, default=25)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_gas_inlet_dense_frac_sweep_latest.json",
    )
    args = parser.parse_args()

    scenarios: list[tuple[str, str, float]] = []
    if args.include_precalc:
        scenarios.append(
            (
                "precalc_baseline",
                "precalc_hydrodynamic_flux",
                float(build_phase2_config().gas_inlet_dense_frac),
            )
        )
    for frac in args.dense_fracs:
        scenarios.append((f"fixed_{frac:.2f}", "fixed", float(frac)))

    rows: list[Row] = []
    print(
        f"Phase2 LU gas inlet dense/bubble sweep | kinetics={args.kinetics} "
        f"r4/r5/r7 per config unless baseline",
        flush=True,
    )
    print(f"ref CO={REF_CO} CO2={REF_CO2} | {len(scenarios)} scenarios", flush=True)
    print("=" * 120, flush=True)

    for label, strategy, dense_cfg in scenarios:
        row = run_one(
            label=label,
            split_strategy=strategy,
            dense_frac_cfg=dense_cfg,
            max_global_iter=int(args.max_global_iter),
            kinetics=str(args.kinetics),
            sync_staged_dense=bool(args.sync_staged_dense),
        )
        rows.append(row)
        _print_row(row)

    converged = [r for r in rows if r.converged_fully]
    best = min(converged, key=lambda r: r.score_co_co2) if converged else None

    payload = {
        "note": f"kinetics={args.kinetics}; sync_staged_dense={args.sync_staged_dense}",
        "ref_co": REF_CO,
        "ref_co2": REF_CO2,
        "scenarios": [asdict(r) for r in rows],
        "best_converged": asdict(best) if best is not None else None,
        "best_any_rms": asdict(min(rows, key=lambda r: r.rms_scaled_final)) if rows else None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("=" * 120, flush=True)
    print(f"Wrote {args.output}", flush=True)
    if best is not None:
        print(
            f"Best converged: {best.label} score={best.score_co_co2:.3f} "
            f"CO={best.co:.4f} CO2={best.co2:.4f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
