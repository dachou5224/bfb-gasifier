#!/usr/bin/env python3
"""Phase2 LU：inlet_fraction vs hamel_eq27_explicit 热损模式对比探针。

先用 ``inlet_fraction`` 跑短 NR，按 cell 级 H_in 估算 implied Q_total，
再以等效 ``Q_wall_W_total`` 跑显式模式，对比出口温度 / RMS / 能量残差。

用法::

    python3 scripts/audit_energy_heat_loss_mode_compare.py
    python3 scripts/audit_energy_heat_loss_mode_compare.py --max-global-iter 8 --json

Source: ``docs/audit_reports/hamel_parity_energy_2026-07-06.md``
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from typing import Any

import _lu_harness  # noqa: F401
from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    build_phase2_reactor,
    data_path,
    ensure_repo_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    HEAT_LOSS_MODE_INLET_FRACTION,
    resolve_inlet_fraction_calibration_W,
    summarize_reactor_inlet_fraction_heat_loss_W,
)
from src.core.reactor import resolve_freeboard_heat_loss_bundle

ensure_repo_path()


@dataclass
class RunSummary:
    label: str
    heat_loss_mode: str
    heat_loss_frac: float
    Q_wall_W_total: float
    Q_heat_exchanger_W_total: float
    exit_T_K: float
    bed_top_T_K: float
    rms_scaled_final: float
    converged: bool
    global_iters: int
    wall_time_s: float
    implied_heat_loss_W: dict[str, float] | None
    implied_heat_loss_positive_W: dict[str, float] | None
    calibration_mode: str
    mean_abs_energy_residual_W: float
    max_abs_energy_residual_W: float


def _mean_energy_residuals(reactor) -> tuple[float, float]:
    vals = [abs(float(c.calc_energy_balance())) for c in reactor.cells]
    if reactor.freeboard_cells:
        vals.extend(abs(float(c.calc_energy_balance())) for c in reactor.freeboard_cells)
    if not vals:
        return 0.0, 0.0
    return float(sum(vals) / len(vals)), float(max(vals))


def _run_case(
    *,
    label: str,
    heat_loss_mode: str,
    heat_loss_frac: float,
    Q_wall_W_total: float,
    Q_heat_exchanger_W_total: float,
    max_global_iter: int,
    n_cells: int,
    n_freeboard_cells: int,
    init_only: bool = False,
    cfg_extra: dict[str, object] | None = None,
) -> RunSummary:
    overrides = {
        "n_cells": n_cells,
        "energy_balance_heat_loss_mode": heat_loss_mode,
        "heat_loss_frac": heat_loss_frac,
        "Q_wall_W_total": Q_wall_W_total,
        "Q_heat_exchanger_W_total": Q_heat_exchanger_W_total,
    }
    if cfg_extra:
        overrides.update(cfg_extra)
    reactor = build_phase2_reactor(
        builder_kwargs={"n_freeboard_cells": n_freeboard_cells},
        cfg_overrides=overrides,
        ladder_level="p2",
    )
    run_phase2_init(reactor)
    if init_only:
        solve_out = {"rms_scaled": 0.0, "converged": False, "global_iters": 0}
        wall_s = 0.0
    else:
        solve_out, wall_s = timed_phase2_solve(
            reactor,
            **phase2_solve_kwargs(max_global_iter=max_global_iter, tol_global=1.0),
        )
    implied_positive = (
        summarize_reactor_inlet_fraction_heat_loss_W(reactor, use_abs_H_in=False)
        if heat_loss_mode == HEAT_LOSS_MODE_INLET_FRACTION
        else None
    )
    implied = None
    calibration_mode = ""
    if heat_loss_mode == HEAT_LOSS_MODE_INLET_FRACTION:
        implied, calibration_mode = resolve_inlet_fraction_calibration_W(reactor)
    mean_e, max_e = _mean_energy_residuals(reactor)
    top = reactor.cells[-1]
    exit_cell = reactor.freeboard_cells[-1] if reactor.freeboard_cells else top
    return RunSummary(
        label=label,
        heat_loss_mode=heat_loss_mode,
        heat_loss_frac=float(heat_loss_frac),
        Q_wall_W_total=float(Q_wall_W_total),
        Q_heat_exchanger_W_total=float(Q_heat_exchanger_W_total),
        exit_T_K=float(exit_cell.T),
        bed_top_T_K=float(top.T),
        rms_scaled_final=float(solve_out.get("rms_scaled", solve_out.get("rms", 0.0)) or 0.0),
        converged=bool(solve_out.get("converged", False)),
        global_iters=int(solve_out.get("global_iters", solve_out.get("iterations", 0)) or 0),
        wall_time_s=float(wall_s),
        implied_heat_loss_W=implied,
        implied_heat_loss_positive_W=implied_positive,
        calibration_mode=calibration_mode,
        mean_abs_energy_residual_W=float(mean_e),
        max_abs_energy_residual_W=float(max_e),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(parser, default=8)
    add_json_output_arg(parser, "phase2_energy_heat_loss_mode_compare_latest.json")
    parser.add_argument("--n-cells", type=int, default=5)
    parser.add_argument("--n-freeboard-cells", type=int, default=4)
    parser.add_argument(
        "--init-only",
        action="store_true",
        help="仅 Vorab 初值后对比 implied Q 与能量残差，不跑 global NR（快探针）",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="缩小网格并关闭 staged Vorab / side blocks（开发快探针）",
    )
    args = parser.parse_args()

    n_cells = 3 if args.fast else int(args.n_cells)
    n_freeboard = 2 if args.fast else int(args.n_freeboard_cells)
    fast_overrides: dict[str, object] = {}
    if args.fast:
        fast_overrides = {
            "explicit_side_blocks_enabled": False,
            "vorab_staged_primary_gas_injection_thesis": False,
            "vorab_a_tier_bed_axial_primary_profile_thesis": False,
            "freeboard_age_quadrature_bins": 2,
        }

    frac_run = _run_case(
        label="inlet_fraction_baseline",
        heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION,
        heat_loss_frac=0.1,
        Q_wall_W_total=0.0,
        Q_heat_exchanger_W_total=0.0,
        max_global_iter=int(args.max_global_iter),
        n_cells=n_cells,
        n_freeboard_cells=n_freeboard,
        init_only=bool(args.init_only),
        cfg_extra=fast_overrides,
    )
    assert frac_run.implied_heat_loss_W is not None
    q_total_signed = float(frac_run.implied_heat_loss_W["total_W"])
    q_total = abs(q_total_signed)
    explicit_run = _run_case(
        label="hamel_eq27_explicit_calibrated",
        heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        heat_loss_frac=0.0,
        Q_wall_W_total=q_total,
        Q_heat_exchanger_W_total=0.0,
        max_global_iter=int(args.max_global_iter),
        n_cells=n_cells,
        n_freeboard_cells=n_freeboard,
        init_only=bool(args.init_only),
        cfg_extra=fast_overrides,
    )
    bundle = resolve_freeboard_heat_loss_bundle(
        build_phase2_reactor(
            builder_kwargs={"n_freeboard_cells": n_freeboard},
            cfg_overrides={"n_cells": n_cells, "Q_wall_W_total": q_total, **fast_overrides},
            ladder_level="p2",
        ).config
    )
    payload: dict[str, Any] = {
        "probe": {
            "n_cells": n_cells,
            "n_freeboard_cells": n_freeboard,
            "max_global_iter": int(args.max_global_iter),
            "init_only": bool(args.init_only),
            "fast": bool(args.fast),
            "calibration_mode": frac_run.calibration_mode,
            "calibrated_Q_wall_W_total": q_total,
            "calibrated_Q_wall_W_total_signed": q_total_signed,
            "inlet_fraction_positive_total_W": float(
                (frac_run.implied_heat_loss_positive_W or {}).get("total_W", 0.0)
            ),
            "freeboard_Q_wall_W_from_bundle": float(bundle.Q_wall_W),
        },
        "inlet_fraction": asdict(frac_run),
        "explicit_calibrated": asdict(explicit_run),
        "delta": {
            "exit_T_K": float(explicit_run.exit_T_K - frac_run.exit_T_K),
            "bed_top_T_K": float(explicit_run.bed_top_T_K - frac_run.bed_top_T_K),
            "rms_scaled_final": float(explicit_run.rms_scaled_final - frac_run.rms_scaled_final),
        },
    }
    out_path = write_json(args.output, payload)
    print(f"Wrote {out_path}")
    print(f"implied Q_total signed = {q_total_signed/1e3:.1f} kW  |  |Q| = {q_total/1e3:.1f} kW  mode={frac_run.calibration_mode}")
    print(
        f"inlet_fraction: exit_T={frac_run.exit_T_K:.1f} K  rms={frac_run.rms_scaled_final:.4g}  "
        f"energy|res|_mean={frac_run.mean_abs_energy_residual_W:.2e} W"
    )
    print(
        f"explicit:       exit_T={explicit_run.exit_T_K:.1f} K  rms={explicit_run.rms_scaled_final:.4g}  "
        f"energy|res|_mean={explicit_run.mean_abs_energy_residual_W:.2e} W"
    )
    print(f"Δexit_T = {payload['delta']['exit_T_K']:+.1f} K")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
