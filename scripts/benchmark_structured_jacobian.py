#!/usr/bin/env python3
"""Small benchmark for structured Jacobian assembly/solve vs dense FD baseline."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter
from typing import Any

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.core.reactor import Reactor
from src.solvers.global_nr_solver import build_equation_scales, build_jacobian_fd, global_residual, pack_reactor
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config
from tests.test_global_nr_solver import _build_initialized_thesis_freeboard_reactor


def _run_case(
    *,
    label: str,
    n_bed: int,
    n_freeboard: int,
    max_global_iter: int,
    tol_global: float,
    jacobian_strategy: str,
    linear_backend: str,
) -> dict[str, Any]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = int(max(n_bed, 1))
    cfg.n_freeboard_cells = int(max(n_freeboard, 0))
    cfg.vorab_major_gibbs_x0 = False

    reactor = Reactor(cfg)
    t0 = perf_counter()
    out = reactor.solve(
        solver="global_nr",
        max_global_iter=int(max(max_global_iter, 1)),
        tol_global=float(tol_global),
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=str(jacobian_strategy),
        nr_linear_solver_backend=str(linear_backend),
    )
    wall_s = perf_counter() - t0
    timing = out.get("nr_timing", {}) or {}
    return {
        "label": str(label),
        "wall_s": float(wall_s),
        "nr_jacobian_strategy": out.get("nr_jacobian_strategy"),
        "nr_linear_solver_backend": out.get("nr_linear_solver_backend_last") or out.get("nr_linear_solver_backend"),
        "jacobian_build_s": float(timing.get("jacobian_build_s", 0.0) or 0.0),
        "linear_solve_s": float(timing.get("linear_solve_s", 0.0) or 0.0),
        "band_lu_s": float(timing.get("band_lu_s", 0.0) or 0.0),
        "side_update_s": float(timing.get("side_update_s", 0.0) or 0.0),
        "nr_structure_validation_ok": out.get("nr_structure_validation_ok"),
        "nr_side_element_count": out.get("nr_side_element_count"),
        "nr_linear_fallback_used": bool(out.get("nr_linear_fallback_used", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", float("inf"))),
        "T_exit_K": float(out["T_profile"][-1]),
        "carbon_conv": float(out.get("carbon_conv", 0.0) or 0.0),
    }


def _run_initialized_jacobian_case(
    *,
    label: str,
    n_bed: int,
    n_freeboard: int,
    jacobian_strategy: str,
) -> dict[str, Any]:
    reactor = _build_initialized_thesis_freeboard_reactor(n_bed=int(max(n_bed, 1)), n_freeboard=int(max(n_freeboard, 0)))
    solver_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(solver_cells)
    F0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    cfg = reactor.config
    eq_scale = build_equation_scales(
        solver_cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )

    t0 = perf_counter()
    J, meta = build_jacobian_fd(
        x,
        F0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy=str(jacobian_strategy),
    )
    wall_s = perf_counter() - t0
    structure = meta.get("jacobian_structure") or {}
    return {
        "label": str(label),
        "wall_s": float(wall_s),
        "jacobian_strategy": str(meta.get("strategy")),
        "residual_cell_calls": int(meta.get("residual_cell_calls", 0)),
        "nnz": int(J.nnz),
        "zero_cols": int(meta.get("zero_cols", 0)),
        "zero_rows": int(meta.get("zero_rows", 0)),
        "side_element_count": structure.get("side_element_count"),
        "structure_validation_ok": structure.get("structure_validation_ok"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark structured Jacobian path against dense_fd baseline")
    parser.add_argument(
        "--mode",
        choices=("solve", "initialized-jacobian"),
        default="solve",
        help="benchmark full Reactor.solve or just initialized Jacobian assembly",
    )
    parser.add_argument("--n-bed", type=int, default=3, help="number of bed cells")
    parser.add_argument("--n-freeboard", type=int, default=2, help="number of freeboard cells")
    parser.add_argument("--max-global-iter", type=int, default=2, help="global NR iterations")
    parser.add_argument("--tol-global", type=float, default=1.0, help="global NR tolerance")
    parser.add_argument("--output-json", type=Path, default=None, help="optional JSON output path")
    args = parser.parse_args()

    if args.mode == "solve":
        structured = _run_case(
            label="structured",
            n_bed=args.n_bed,
            n_freeboard=args.n_freeboard,
            max_global_iter=args.max_global_iter,
            tol_global=args.tol_global,
            jacobian_strategy="band_plus_side_elements_structured",
            linear_backend="structured_direct",
        )
        dense = _run_case(
            label="dense_baseline",
            n_bed=args.n_bed,
            n_freeboard=args.n_freeboard,
            max_global_iter=args.max_global_iter,
            tol_global=args.tol_global,
            jacobian_strategy="dense_fd",
            linear_backend="sparse_direct_fallback",
        )
    else:
        structured = _run_initialized_jacobian_case(
            label="structured",
            n_bed=args.n_bed,
            n_freeboard=args.n_freeboard,
            jacobian_strategy="band_plus_side_elements_structured",
        )
        dense = _run_initialized_jacobian_case(
            label="dense_baseline",
            n_bed=args.n_bed,
            n_freeboard=args.n_freeboard,
            jacobian_strategy="dense_fd",
        )
    speedup = (dense["wall_s"] - structured["wall_s"]) / max(dense["wall_s"], 1e-12)

    payload = {
        "config": {
            "mode": str(args.mode),
            "n_bed": int(args.n_bed),
            "n_freeboard": int(args.n_freeboard),
            "max_global_iter": int(args.max_global_iter),
            "tol_global": float(args.tol_global),
        },
        "structured": structured,
        "dense_baseline": dense,
        "wall_speedup_fraction": float(speedup),
    }

    print("=== Structured Jacobian Benchmark ===")
    print(json.dumps(payload, ensure_ascii=False, indent=2))

    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"saved: {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
