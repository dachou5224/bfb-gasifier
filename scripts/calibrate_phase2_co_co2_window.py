#!/usr/bin/env python3
"""Phase2 LU：CO/CO2 干基组分局部标定（甲烷暂不纳入目标）。

固定 Phase2 freeboard global-NR 口径，在 ``r4_scale`` / ``r5_scale`` / ``r8_scale``
小窗口内搜索，目标仅最小化 CO、CO2 相对 JSON 参考偏差，并要求 ``converged_fully``。

Usage
-----
    python3 scripts/calibrate_phase2_co_co2_window.py --mode sweep1d --axis r8
    python3 scripts/calibrate_phase2_co_co2_window.py --mode grid_r2_r4_r5
"""

from __future__ import annotations

from _lu_harness import (
    REF_CO,
    REF_CO2,
    add_max_global_iter_arg,
    build_phase2_reactor,
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

@dataclass
class Row:
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
    wall_time_s: float


def _score_co_co2(co: float, co2: float) -> float:
    return abs(co - REF_CO) / REF_CO + abs(co2 - REF_CO2) / REF_CO2


def run_one(
    *,
    r2_scale: float,
    r4_scale: float,
    r5_scale: float,
    r8_scale: float,
    max_global_iter: int,
) -> Row:
    cfg_overrides = {
        "vorab_major_gibbs_x0": False,
        "r2_scale": float(max(r2_scale, 0.0)),
        "r4_scale": float(max(r4_scale, 0.0)),
        "r5_scale": float(max(r5_scale, 0.0)),
        "r8_scale": float(max(r8_scale, 0.0)),
    }

    solve_kwargs = phase2_solve_kwargs()
    solve_kwargs["max_global_iter"] = int(max_global_iter)
    solve_kwargs["verbose"] = False

    t0 = perf_counter()
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg_overrides)
    out = reactor.solve(**solve_kwargs)
    wall_s = perf_counter() - t0

    y = out.get("exit_gas_dry") or {}
    bed = out.get("bed_exit_gas_dry") or {}
    co = float(y.get("CO", 0.0))
    co2 = float(y.get("CO2", 0.0))
    return Row(
        r2_scale=float(r2_scale),
        r4_scale=float(r4_scale),
        r5_scale=float(r5_scale),
        r8_scale=float(r8_scale),
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
        wall_time_s=float(wall_s),
    )


def _print_row(prefix: str, row: Row) -> None:
    conv = "OK" if row.converged_fully else "FAIL"
    print(
        f"{prefix} r2={row.r2_scale:.2f} r4={row.r4_scale:.3f} r5={row.r5_scale:.3f} r8={row.r8_scale:.3f} | "
        f"conv={conv} rms={row.rms_scaled_final:.4f} T={row.t_exit_K:.1f} | "
        f"CO={row.co:.4f} CO2={row.co2:.4f} bedCO={row.bed_co:.4f} bedCO2={row.bed_co2:.4f} | "
        f"H2={row.h2:.4f} CH4={row.ch4:.4f} | score={row.score_co_co2:.4f} ({row.wall_time_s:.0f}s)",
        flush=True,
    )


def sweep1d(
    axis: str,
    *,
    base_r2: float,
    base_r4: float,
    base_r5: float,
    base_r8: float,
    max_global_iter: int,
) -> list[Row]:
    grids = {
        "r2": [1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0],
        "r4": [0.50, 0.75, 1.00, 1.25, 1.50, 2.00],
        "r5": [0.25, 0.35, 0.50, 0.65, 0.75, 1.00],
        "r8": [0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00],
    }
    if axis not in grids:
        raise ValueError(f"Unknown axis {axis!r}")

    rows: list[Row] = []
    n = len(grids[axis])
    for i, value in enumerate(grids[axis], start=1):
        r2, r4, r5, r8 = base_r2, base_r4, base_r5, base_r8
        if axis == "r2":
            r2 = float(value)
        elif axis == "r4":
            r4 = float(value)
        elif axis == "r5":
            r5 = float(value)
        else:
            r8 = float(value)
        row = run_one(
            r2_scale=r2,
            r4_scale=r4,
            r5_scale=r5,
            r8_scale=r8,
            max_global_iter=max_global_iter,
        )
        rows.append(row)
        _print_row(f"[{i}/{n}]", row)
    return rows


def grid_r2_r5(
    *,
    r2_grid: list[float],
    r5_grid: list[float],
    r4_scale: float,
    r8_scale: float,
    max_global_iter: int,
) -> list[Row]:
    rows: list[Row] = []
    total = len(r2_grid) * len(r5_grid)
    n = 0
    for r2 in r2_grid:
        for r5 in r5_grid:
            n += 1
            row = run_one(
                r2_scale=float(r2),
                r4_scale=float(r4_scale),
                r5_scale=float(r5),
                r8_scale=float(r8_scale),
                max_global_iter=max_global_iter,
            )
            rows.append(row)
            _print_row(f"[{n}/{total}]", row)
    return rows


def grid_r2_r4_r5(
    *,
    r2_grid: list[float],
    r4_grid: list[float],
    r5_grid: list[float],
    r8_scale: float,
    max_global_iter: int,
) -> list[Row]:
    rows: list[Row] = []
    total = len(r2_grid) * len(r4_grid) * len(r5_grid)
    n = 0
    for r2 in r2_grid:
        for r4 in r4_grid:
            for r5 in r5_grid:
                n += 1
                row = run_one(
                    r2_scale=float(r2),
                    r4_scale=float(r4),
                    r5_scale=float(r5),
                    r8_scale=float(r8_scale),
                    max_global_iter=max_global_iter,
                )
                rows.append(row)
                _print_row(f"[{n}/{total}]", row)
    return rows


def grid_search(
    *,
    base_r2: float,
    r4_grid: list[float],
    r5_grid: list[float],
    r8_grid: list[float],
    max_global_iter: int,
) -> list[Row]:
    rows: list[Row] = []
    total = len(r4_grid) * len(r5_grid) * len(r8_grid)
    n = 0
    for r4 in r4_grid:
        for r5 in r5_grid:
            for r8 in r8_grid:
                n += 1
                row = run_one(
                    r2_scale=float(base_r2),
                    r4_scale=float(r4),
                    r5_scale=float(r5),
                    r8_scale=float(r8),
                    max_global_iter=max_global_iter,
                )
                rows.append(row)
                _print_row(f"[{n}/{total}]", row)
    return rows


def _best_converged(rows: list[Row]) -> Row | None:
    ok = [r for r in rows if r.converged_fully]
    if not ok:
        return None
    return min(ok, key=lambda r: r.score_co_co2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("sweep1d", "grid", "grid_r2_r5", "grid_r2_r4_r5"),
        default="sweep1d",
    )
    parser.add_argument("--axis", choices=("r2", "r4", "r5", "r8"), default="r5")
    parser.add_argument("--max-global-iter", type=int, default=20)
    parser.add_argument("--base-r2", type=float, default=1.00)
    parser.add_argument("--base-r4", type=float, default=0.50)
    parser.add_argument("--base-r5", type=float, default=0.25)
    parser.add_argument("--base-r8", type=float, default=1.00)
    parser.add_argument("--grid-r2", type=str, default="1,2,3,5,8,12")
    parser.add_argument("--grid-r4", type=str, default="0.25,0.50,1.00,1.50,2.00")
    parser.add_argument("--grid-r5", type=str, default="0.30,0.45,0.60,0.75")
    parser.add_argument("--grid-r8", type=str, default="0.50,1.00,1.50,2.00")
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_co_co2_calibration_latest.json",
    )
    args = parser.parse_args()

    def _parse_grid(text: str) -> list[float]:
        return [float(x.strip()) for x in text.split(",") if x.strip()]

    print("=" * 120)
    print("Phase2 CO/CO2 calibration (CH4 excluded from objective)")
    print(f"ref: CO={REF_CO:.3f}, CO2={REF_CO2:.3f}; max_global_iter={args.max_global_iter}")
    print("=" * 120)

    if args.mode == "sweep1d":
        rows = sweep1d(
            args.axis,
            base_r2=args.base_r2,
            base_r4=args.base_r4,
            base_r5=args.base_r5,
            base_r8=args.base_r8,
            max_global_iter=args.max_global_iter,
        )
    elif args.mode == "grid_r2_r5":
        rows = grid_r2_r5(
            r2_grid=_parse_grid(args.grid_r2),
            r5_grid=_parse_grid(args.grid_r5),
            r4_scale=args.base_r4,
            r8_scale=args.base_r8,
            max_global_iter=args.max_global_iter,
        )
    elif args.mode == "grid_r2_r4_r5":
        rows = grid_r2_r4_r5(
            r2_grid=_parse_grid(args.grid_r2),
            r4_grid=_parse_grid(args.grid_r4),
            r5_grid=_parse_grid(args.grid_r5),
            r8_scale=args.base_r8,
            max_global_iter=args.max_global_iter,
        )
    else:
        rows = grid_search(
            base_r2=args.base_r2,
            r4_grid=_parse_grid(args.grid_r4),
            r5_grid=_parse_grid(args.grid_r5),
            r8_grid=_parse_grid(args.grid_r8),
            max_global_iter=args.max_global_iter,
        )

    best = _best_converged(rows)
    payload = {
        "ref": {"CO": REF_CO, "CO2": REF_CO2},
        "mode": args.mode,
        "axis": args.axis if args.mode == "sweep1d" else None,
        "max_global_iter": int(args.max_global_iter),
        "rows": [asdict(r) for r in rows],
        "best_converged": asdict(best) if best is not None else None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("-" * 120)
    if best is None:
        print("无 converged_fully 点；请放宽网格或提高 max_global_iter。")
        return 1
    print(
        f"BEST converged: r2={best.r2_scale:.2f} r4={best.r4_scale:.3f} r5={best.r5_scale:.3f} r8={best.r8_scale:.3f} | "
        f"CO={best.co:.4f} CO2={best.co2:.4f} bedCO={best.bed_co:.4f} bedCO2={best.bed_co2:.4f} score={best.score_co_co2:.4f}"
    )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
