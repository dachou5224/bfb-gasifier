#!/usr/bin/env python3
"""VTT 单回路全量 NR 残差审计：全局 / 逐 cell / 主组分 phase-split 分解。

用法::

    python scripts/vtt_full_sim_residual_audit.py
    python scripts/vtt_full_sim_residual_audit.py --max-iter 40 --json-out data/vtt_residual_audit.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from _lu_harness import repo_root, run_phase2_init

REPO_ROOT = repo_root()

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES
from src.solvers.global_nr_solver import (
    _project_cell_residual_to_main_nr,
    _residual_gas_phase_metrics,
    _residual_group_metrics,
    build_equation_scales,
    global_residual,
    pack_reactor,
)
from src.solvers.nr_indexing import (
    main_nr_gas_idx_for_cell,
    main_nr_gas_phase_split_positions_for_cell,
    n_gas_var,
)
from tests.validation_case_utils import (
    VTT_MINIMAL_LOOP_SOLVE_KWARGS,
    build_vtt_minimal_loop_nr_acceptance_config,
)


def _scaled_F(reactor: Reactor, cells: list) -> np.ndarray:
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    return np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)


def _species_split_rows(reactor: Reactor, cell) -> list[dict[str, float | str]]:
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    res = np.asarray(_project_cell_residual_to_main_nr(cell), dtype=np.float64)
    sc = build_equation_scales([cell], ref_gas, ref_solid, ref_energy)
    hat = res / sc
    gas_idx = main_nr_gas_idx_for_cell(cell)
    pos = main_nr_gas_phase_split_positions_for_cell(cell)
    half = n_gas_var(cell) // 2
    rows: list[dict[str, float | str]] = []
    for p in pos:
        sp = GAS_SPECIES[int(gas_idx[p])]
        dense = float(hat[p])
        bubble = float(hat[half + p])
        rows.append(
            {
                "species": sp,
                "dense_minus_bubble_abs": abs(dense - bubble),
                "dense_scaled": dense,
                "bubble_scaled": bubble,
            }
        )
    rows.sort(key=lambda r: float(r["dense_minus_bubble_abs"]), reverse=True)
    return rows


def _audit(reactor: Reactor, cells: list, label: str) -> dict:
    F_hat = _scaled_F(reactor, cells)
    gm = _residual_group_metrics(F_hat, cells)
    gp = _residual_gas_phase_metrics(F_hat, cells)
    bed_splits = {}
    for i, cell in enumerate(cells):
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        bed_splits[f"bed{i}"] = _species_split_rows(reactor, cell)[:7]
    return {
        "label": label,
        "global": {**gm, **gp, "rms_all": float(np.sqrt(np.mean(F_hat**2)))},
        "bed_major_phase_split": bed_splits,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="VTT full NR residual audit")
    parser.add_argument("--max-iter", type=int, default=int(VTT_MINIMAL_LOOP_SOLVE_KWARGS["max_iter"]))
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    cfg = build_vtt_minimal_loop_nr_acceptance_config()
    reactor = Reactor(cfg)
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    cells = reactor._solver_cells_for_nr()
    init_audit = _audit(reactor, cells, "init")

    result = reactor._solve_global_nr(
        skip_precalc=True,
        max_iter=int(args.max_iter),
        tol=float(VTT_MINIMAL_LOOP_SOLVE_KWARGS["tol"]),
        init_strategy=str(VTT_MINIMAL_LOOP_SOLVE_KWARGS["init_strategy"]),
    )
    cells = reactor._solver_cells_for_nr()
    final_audit = _audit(reactor, cells, "final")

    payload = {
        "solve": {
            "converged": bool(result.get("converged")),
            "converged_fully": bool(result.get("converged_fully", result.get("converged"))),
            "rms_scaled_final": float(result.get("rms_scaled_final", np.nan)),
            "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
        },
        "init": init_audit,
        "final": final_audit,
    }

    print("=== VTT full sim residual audit ===")
    print(
        f"converged={payload['solve']['converged']} "
        f"rms={payload['solve']['rms_scaled_final']:.5f} "
        f"outers={payload['solve']['nr_outer_iters']}"
    )
    print(
        f"init  split_rms={init_audit['global']['gas_phase_split_rms']:.4f} "
        f"gas_rms={init_audit['global']['gas_rms']:.4f}"
    )
    print(
        f"final split_rms={final_audit['global']['gas_phase_split_rms']:.5f} "
        f"gas={final_audit['global']['gas_rms']:.5f} "
        f"solid={final_audit['global']['solid_rms']:.5f} "
        f"energy={final_audit['global']['energy_rms']:.5f}"
    )
    for bed, rows in final_audit["bed_major_phase_split"].items():
        top = rows[0] if rows else None
        if top:
            print(
                f"  {bed} worst split: {top['species']} "
                f"|d-b|={float(top['dense_minus_bubble_abs']):.4e}"
            )

    if args.json_out:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
