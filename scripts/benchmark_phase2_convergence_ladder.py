#!/usr/bin/env python3
"""Phase2 LU 收敛阶梯基准：baseline → P0 → P1 → P2。

固定 ``PHASE2_HTW_LU_FREEBOARD`` 口径，逐级打开收敛改进开关并量化 init/final 残差。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

from _lu_harness import (
    LADDER_LEVELS,
    add_json_output_arg,
    add_max_global_iter_arg,
    build_phase2_reactor,
    ensure_repo_path,
    phase2_solve_kwargs,
    run_phase2_init,
)

ensure_repo_path()

from src.core.connectivity import _measure_nr_init_scaled_metrics
from src.core.reactor import Reactor
from src.solvers.global_nr_solver import (
    _residual_gas_phase_metrics,
    build_equation_scales,
    global_residual,
    pack_reactor,
)


def _init_split_metrics(reactor: Reactor, cells: list) -> dict[str, float]:
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    gas_phase = _residual_gas_phase_metrics(F_hat, cells)
    return {
        "rms_scaled_gas_phase_split": float(gas_phase["gas_phase_split_rms"]),
        "max_abs_scaled_gas_phase_split": float(gas_phase["gas_phase_split_max_abs"]),
    }


def run_level(level: str, *, max_global_iter: int) -> dict[str, object]:
    if level not in LADDER_LEVELS:
        raise ValueError(f"Unknown level {level!r}; expected one of {sorted(LADDER_LEVELS)}")

    reactor = build_phase2_reactor(ladder_level=level)
    init_strategy = str(phase2_solve_kwargs().get("nr_init_strategy", "vorabrechnung"))
    run_phase2_init(reactor, init_strategy=init_strategy)
    solver_cells = reactor._solver_cells_for_nr()
    cfg = reactor.config
    init_metrics = _measure_nr_init_scaled_metrics(
        solver_cells, cfg, apply_bc_fn=reactor._apply_all_bc_for_nr
    )
    init_split = _init_split_metrics(reactor, solver_cells)

    t0 = perf_counter()
    solve_kwargs = phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False)
    result = reactor._solve_global_nr(
        max_iter=int(max_global_iter),
        tol=float(solve_kwargs.get("tol_global", 1.0)),
        verbose=False,
        init_strategy=init_strategy,
        jacobian_strategy=solve_kwargs.get("nr_jacobian_strategy"),
        skip_precalc=True,
    )
    wall_s = perf_counter() - t0

    joint_diag = getattr(reactor, "_bottom_joint_preprojection_run_diag", {}) or {}
    zone_diag = getattr(reactor, "_bottom_zone_solid_energy_preprojection_diag", {}) or {}
    counts = result.get("counts") or {}

    return {
        "level": level,
        "max_global_iter": int(max_global_iter),
        "wall_time_s": float(wall_s),
        "init_rms_scaled": float(init_metrics["rms_scaled"]),
        "init_max_abs_scaled": float(init_metrics["max_abs_scaled"]),
        "init_solid_rms": float(init_metrics.get("solid_rms", np.nan)) if "solid_rms" in init_metrics else float("nan"),
        "init_energy_rms": float(init_metrics.get("energy_rms", np.nan)) if "energy_rms" in init_metrics else float("nan"),
        "init_split_max": float(init_split["max_abs_scaled_gas_phase_split"]),
        "zone_bed1_energy_t_preprojection": bool(zone_diag.get("bed1_energy_temperature_preprojection_accepted", False)),
        "zone_solid_preprojections_accepted": float(zone_diag.get("solid_preprojections_accepted", 0.0)),
        "outer_preinner_last_split_rms_after": float(
            (getattr(reactor, "_outer_preinner_last_diag", {}) or {}).get("split_rms_after", 0.0)
        ),
        "outer_preinner_last_solid_rms_after": float(
            (getattr(reactor, "_outer_preinner_last_diag", {}) or {}).get("solid_rms_after", 0.0)
        ),
        "outer_preinner_last_energy_rms_after": float(
            (getattr(reactor, "_outer_preinner_last_diag", {}) or {}).get("energy_rms_after", 0.0)
        ),
        "outer_preinner_last_bed0_accepted": bool(
            (getattr(reactor, "_outer_preinner_last_diag", {}) or {}).get("bed0_accepted", False)
        ),
        "joint_accepted": bool(joint_diag.get("joint", False)),
        "joint_global_gate_pass": bool(joint_diag.get("joint_global_gate_pass", False)),
        "joint_gas_relax": float(joint_diag.get("joint_gas_relax", 0.0)),
        "joint_merit_scaled_before": float(joint_diag.get("joint_merit_scaled_before", 0.0)),
        "joint_merit_scaled_after": float(joint_diag.get("joint_merit_scaled_after", 0.0)),
        "bed0_T_K": float(reactor.cells[0].T),
        "converged": bool(result.get("converged", False)),
        "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
        "nr_outer_max": int(result.get("nr_outer_max", 0)),
        "n_iter": int(result.get("n_iter", 0)),
        "rms_scaled_final": float(result.get("rms_scaled_final", np.nan)),
        "rms_scaled_solid_final": float(result.get("rms_scaled_solid_final", np.nan)),
        "rms_scaled_energy_final": float(result.get("rms_scaled_energy_final", np.nan)),
        "rms_scaled_gas_phase_split_final": float(result.get("rms_scaled_gas_phase_split_final", np.nan)),
        "max_abs_scaled_final": float(result.get("max_abs_scaled_final", np.nan)),
        "max_abs_scaled_gas_phase_split_final": float(result.get("max_abs_scaled_gas_phase_split_final", np.nan)),
        "line_search_failures": int(counts.get("line_search_failures", 0)),
        "hydro_refresh_on_accept": int(
            (result.get("nr_counts") or result.get("counts") or {}).get(
                "hydrodynamics_refresh_on_accept",
                0,
            )
        ),
        "nr_include_pyrolysis_gas_in_unknowns": bool(
            getattr(cfg, "nr_include_pyrolysis_gas_species_in_unknowns_thesis", False)
        ),
        "overrides": dict(LADDER_LEVELS[level]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--level",
        choices=[*LADDER_LEVELS.keys(), "all"],
        default="all",
        help="Run one ladder level or all sequentially",
    )
    add_max_global_iter_arg(parser, default=6)
    add_json_output_arg(parser, "phase2_convergence_ladder_latest.json")
    args = parser.parse_args()

    levels = list(LADDER_LEVELS.keys()) if args.level == "all" else [str(args.level)]
    if args.level == "all":
        import subprocess

        results: list[dict[str, object]] = []
        for level in levels:
            out_path = args.output.with_name(f"{args.output.stem}_{level}{args.output.suffix}")
            cmd = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--level",
                level,
                "--max-global-iter",
                str(int(args.max_global_iter)),
                "--output",
                str(out_path),
            ]
            proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"Ladder level {level!r} failed:\n{proc.stderr or proc.stdout}")
            results.append(json.loads(out_path.read_text(encoding="utf-8"))["results"][0])
    else:
        results = [run_level(levels[0], max_global_iter=int(args.max_global_iter))]

    payload = {
        "max_global_iter": int(args.max_global_iter),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
