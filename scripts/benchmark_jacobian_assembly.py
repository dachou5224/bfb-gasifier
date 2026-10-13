#!/usr/bin/env python3
"""Benchmark Jacobian assembly: dense_fd vs structured vs hybrid (LU / VTT).

Usage:
    python3 scripts/benchmark_jacobian_assembly.py
    python3 scripts/benchmark_jacobian_assembly.py --case vtt --skip-dense
    python3 scripts/benchmark_jacobian_assembly.py --output-json data/jacobian_assembly_benchmark_latest.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter
from typing import Any

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np

from src.solvers.global_nr_solver import (
    build_equation_scales,
    build_jacobian_fd,
    global_residual,
    pack_reactor,
    resolve_hybrid_jacobian_strategy,
)
from tests.global_nr_fixtures import (
    build_initialized_lu_reactor,
    build_initialized_vtt_minimal_loop_reactor,
)


def _assemble_once(
    *,
    label: str,
    strategy: str,
    reactor,
    solver_cells,
    connectivity_graph=None,
) -> dict[str, Any]:
    cfg = reactor.config
    x = pack_reactor(solver_cells)
    f0 = global_residual(x, solver_cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        solver_cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=cfg.fuel_feed,
        ref_energy_W=cfg.fuel_feed * 20e6,
    )
    t0 = perf_counter()
    jacobian, meta = build_jacobian_fd(
        x,
        f0,
        solver_cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy=strategy,
        connectivity_graph=connectivity_graph,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        affected_residual_cells_fn=reactor._affected_nr_residual_cells,
    )
    wall_s = perf_counter() - t0
    structure = meta.get("jacobian_structure") or {}
    return {
        "label": label,
        "strategy": str(meta.get("strategy")),
        "wall_s": float(wall_s),
        "residual_cell_calls": int(meta.get("residual_cell_calls", 0)),
        "nnz": int(jacobian.nnz),
        "zero_cols": int(meta.get("zero_cols", 0)),
        "zero_rows": int(meta.get("zero_rows", 0)),
        "unexpected_block_count": int(structure.get("unexpected_block_count", -1)),
        "structure_validation_ok": bool(structure.get("structure_validation_ok", False)),
    }


def _max_abs_diff(
    build_reactor,
    strategy_a: str,
    strategy_b: str,
    *,
    solver_cells_fn,
    connectivity_graph=None,
) -> float:
    def _one(strategy: str):
        reactor = build_reactor()
        cells = solver_cells_fn(reactor)
        cfg = reactor.config
        x = pack_reactor(cells)
        f0 = global_residual(x, cells, reactor._apply_all_bc_for_nr)
        eq_scale = build_equation_scales(
            cells,
            ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
            ref_solid_kg_s=cfg.fuel_feed,
            ref_energy_W=cfg.fuel_feed * 20e6,
        )
        jacobian, _ = build_jacobian_fd(
            x,
            f0,
            cells,
            reactor._apply_all_bc_for_nr,
            eq_scale,
            strategy=strategy,
            connectivity_graph=connectivity_graph,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
            affected_residual_cells_fn=reactor._affected_nr_residual_cells,
        )
        return jacobian

    ja = _one(strategy_a)
    jb = _one(strategy_b)
    return float(np.max(np.abs((ja - jb).toarray())))


def benchmark_lu(*, skip_dense: bool) -> dict[str, Any]:
    reactor = build_initialized_lu_reactor(n_cells=3)
    cells = reactor.cells
    structured = "block_tridiag_structured"
    hybrid = resolve_hybrid_jacobian_strategy(structured)
    rows = [
        _assemble_once(label="structured", strategy=structured, reactor=reactor, solver_cells=cells),
        _assemble_once(label="hybrid", strategy=hybrid, reactor=reactor, solver_cells=cells),
    ]
    if not skip_dense:
        rows.append(_assemble_once(label="dense", strategy="dense_fd", reactor=reactor, solver_cells=cells))
    parity = {
        "hybrid_vs_structured_max_abs": _max_abs_diff(
            lambda: build_initialized_lu_reactor(n_cells=3),
            hybrid,
            structured,
            solver_cells_fn=lambda r: r.cells,
        ),
    }
    if not skip_dense:
        parity["structured_vs_dense_max_abs"] = _max_abs_diff(
            lambda: build_initialized_lu_reactor(n_cells=3),
            structured,
            "dense_fd",
            solver_cells_fn=lambda r: r.cells,
        )
    return {"case": "lu_3bed", "runs": rows, "parity": parity}


def benchmark_vtt(*, skip_dense: bool) -> dict[str, Any]:
    reactor = build_initialized_vtt_minimal_loop_reactor()
    cells = reactor._solver_cells_for_nr()
    graph = reactor.connectivity_graph
    structured = "band_plus_side_elements_structured"
    hybrid = resolve_hybrid_jacobian_strategy(structured)
    rows = [
        _assemble_once(
            label="structured",
            strategy=structured,
            reactor=reactor,
            solver_cells=cells,
            connectivity_graph=graph,
        ),
        _assemble_once(
            label="hybrid",
            strategy=hybrid,
            reactor=reactor,
            solver_cells=cells,
            connectivity_graph=graph,
        ),
    ]
    if not skip_dense:
        rows.append(
            _assemble_once(
                label="dense",
                strategy="dense_fd",
                reactor=reactor,
                solver_cells=cells,
                connectivity_graph=graph,
            )
        )
    parity = {
        "hybrid_vs_structured_max_abs": _max_abs_diff(
            build_initialized_vtt_minimal_loop_reactor,
            hybrid,
            structured,
            solver_cells_fn=lambda r: r._solver_cells_for_nr(),
            connectivity_graph=graph,
        ),
    }
    if not skip_dense:
        parity["structured_vs_dense_max_abs"] = _max_abs_diff(
            build_initialized_vtt_minimal_loop_reactor,
            structured,
            "dense_fd",
            solver_cells_fn=lambda r: r._solver_cells_for_nr(),
            connectivity_graph=graph,
        )
    return {"case": "vtt_minimal_loop", "runs": rows, "parity": parity}


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Jacobian assembly strategies")
    parser.add_argument(
        "--case",
        choices=("lu", "vtt", "both"),
        default="both",
        help="which initialized reactor graph to benchmark",
    )
    parser.add_argument(
        "--skip-dense",
        action="store_true",
        help="skip dense_fd (slow on VTT)",
    )
    parser.add_argument("--output-json", type=Path, default=None)
    args = parser.parse_args()

    payload: dict[str, Any] = {"cases": []}
    if args.case in ("lu", "both"):
        payload["cases"].append(benchmark_lu(skip_dense=bool(args.skip_dense)))
    if args.case in ("vtt", "both"):
        payload["cases"].append(benchmark_vtt(skip_dense=bool(args.skip_dense)))

    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"saved: {args.output_json}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
