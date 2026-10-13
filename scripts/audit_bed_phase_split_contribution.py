#!/usr/bin/env python3
"""逐 bed 主组分 phase-split 贡献审计（对齐 global NR split merit）。

用法：
  python scripts/audit_bed_phase_split_contribution.py
  python scripts/audit_bed_phase_split_contribution.py --max-global-iter 20 --output data/phase2_bed_split_contribution.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity.preprojection import (
    _bed_major_gas_phase_split_rms,
    preproject_bed_major_gas_phase_split_to_exchange_closure,
)
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
    n_var,
)


def _scaled_F_hat(reactor: Any, cells: list[Any]) -> np.ndarray:
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    return np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)


def _cell_nr_split_terms(F_hat: np.ndarray, cell: Any, *, offset: int) -> np.ndarray:
    from src.solvers.nr_indexing import gas_mode_for_nr

    gas_mode = gas_mode_for_nr(cell)
    ng = n_gas_var(cell)
    if gas_mode != "two_phase" or ng <= 0:
        return np.zeros(0, dtype=np.float64)
    half = ng // 2
    dense = F_hat[offset : offset + half]
    bubble = F_hat[offset + half : offset + ng]
    split_pos = main_nr_gas_phase_split_positions_for_cell(cell)
    if split_pos.size == 0:
        return np.zeros(0, dtype=np.float64)
    return np.asarray(dense[split_pos] - bubble[split_pos], dtype=np.float64)


def _species_split_rows(reactor: Any, cell: Any, F_hat: np.ndarray, *, offset: int) -> list[dict[str, float | str]]:
    from src.solvers.nr_indexing import gas_mode_for_nr

    gas_mode = gas_mode_for_nr(cell)
    ng = n_gas_var(cell)
    if gas_mode != "two_phase" or ng <= 0:
        return []
    half = ng // 2
    dense = F_hat[offset : offset + half]
    bubble = F_hat[offset + half : offset + ng]
    gas_idx = main_nr_gas_idx_for_cell(cell)
    pos = main_nr_gas_phase_split_positions_for_cell(cell)
    rows: list[dict[str, float | str]] = []
    for p in pos:
        sp = GAS_SPECIES[int(gas_idx[p])]
        d_val = float(dense[p])
        b_val = float(bubble[p])
        rows.append(
            {
                "species": sp,
                "dense_minus_bubble": float(d_val - b_val),
                "abs_split": float(abs(d_val - b_val)),
            }
        )
    rows.sort(key=lambda r: float(r["abs_split"]), reverse=True)
    return rows


def _audit_cells(reactor: Any, cells: list[Any], *, label: str) -> dict[str, Any]:
    reactor._apply_all_bc_for_nr()
    F_hat = _scaled_F_hat(reactor, cells)
    global_gas = _residual_gas_phase_metrics(F_hat, cells)
    global_group = _residual_group_metrics(F_hat, cells)

    offset = 0
    bed_rows: list[dict[str, Any]] = []
    split_sq_sum = 0.0
    split_term_count = 0
    for i, cell in enumerate(cells):
        nv = n_var(cell)
        terms = _cell_nr_split_terms(F_hat, cell, offset=offset)
        offset += nv
        if terms.size == 0:
            continue
        cell_rms = float(np.sqrt(np.mean(terms**2)))
        cell_max = float(np.max(np.abs(terms)))
        split_sq_sum += float(np.sum(terms**2))
        split_term_count += int(terms.size)
        preproj_rms = float(_bed_major_gas_phase_split_rms(cell, reaction_rate_multiplier=1.0))
        bed_rows.append(
            {
                "cell_index": int(i),
                "cell_type": str(getattr(cell, "cell_type", "bed")),
                "T_K": float(getattr(cell, "T", float("nan"))),
                "nr_split_rms": cell_rms,
                "nr_split_max_abs": cell_max,
                "preproj_major_split_rms": preproj_rms,
                "n_split_terms": int(terms.size),
                "top_species": _species_split_rows(reactor, cell, F_hat, offset=offset - nv)[:5],
            }
        )

    bed_rows.sort(key=lambda r: float(r["nr_split_rms"]), reverse=True)
    for row in bed_rows:
        n = int(row["n_split_terms"])
        row["split_variance_fraction"] = float(
            (float(row["nr_split_rms"]) ** 2) * n / max(split_sq_sum, 1.0e-30)
        )

    return {
        "label": label,
        "global": {
            **global_group,
            **global_gas,
            "merit_like": float(
                max(
                    global_group["gas_rms"],
                    global_group["solid_rms"],
                    global_group["energy_rms"],
                    global_gas["gas_phase_split_rms"],
                )
            ),
        },
        "bed_rows": bed_rows,
        "split_term_count": int(split_term_count),
    }


def _probe_single_bed_t2(
    reactor: Any,
    cells: list[Any],
    *,
    bed_index: int,
    improvement_factor: float,
) -> dict[str, float | bool]:
    from src.core.connectivity import _measure_nr_init_scaled_metrics

    cfg = reactor.config
    metrics_before = _measure_nr_init_scaled_metrics(
        cells, cfg, apply_bc_fn=reactor._apply_all_bc_for_nr
    )
    accepted = bool(
        preproject_bed_major_gas_phase_split_to_exchange_closure(
            cells,
            cfg,
            cell_index=int(bed_index),
            improvement_factor=float(improvement_factor),
            preserve_primary_inlet_totals=(int(bed_index) == 0),
            apply_bc_fn=reactor._apply_all_bc_for_nr,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        )
    )
    metrics_after = _measure_nr_init_scaled_metrics(
        cells, cfg, apply_bc_fn=reactor._apply_all_bc_for_nr
    )
    return {
        "bed_index": float(bed_index),
        "accepted": accepted,
        "split_rms_before": float(metrics_before["gas_phase_split_rms"]),
        "split_rms_after": float(metrics_after["gas_phase_split_rms"]),
        "merit_before": float(metrics_before["merit_scaled"]),
        "merit_after": float(metrics_after["merit_scaled"]),
        "improvement_factor": float(improvement_factor),
    }


def run_audit(*, max_global_iter: int, ladder_level: str) -> dict[str, Any]:
    reactor = build_phase2_reactor(ladder_level=ladder_level)
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    cells_init = reactor._solver_cells_for_nr()
    init_audit = _audit_cells(reactor, cells_init, label="init")

    solve_kwargs = phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False)
    result = reactor.solve(**solve_kwargs)
    cells_final = reactor._solver_cells_for_nr()
    final_audit = _audit_cells(reactor, cells_final, label="final")

    t2_probes: list[dict[str, float | bool]] = []
    top_beds = [int(r["cell_index"]) for r in final_audit["bed_rows"][:4]]
    for bed_index in sorted(set(top_beds)):
        t2_probes.append(
            _probe_single_bed_t2(reactor, cells_final, bed_index=bed_index, improvement_factor=0.85)
        )

    return {
        "max_global_iter": int(max_global_iter),
        "ladder_level": ladder_level,
        "solve_summary": {
            "converged": bool(result.get("converged", False)),
            "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
            "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
            "rms_scaled_gas_phase_split_final": float(
                result.get("rms_scaled_gas_phase_split_final", float("nan"))
            ),
        },
        "init": init_audit,
        "final": final_audit,
        "single_bed_t2_probes_final_state": t2_probes,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=20)
    parser.add_argument("--ladder-level", default="p2")
    parser.add_argument(
        "--output",
        type=Path,
        default=data_path("phase2_bed_phase_split_contribution.json"),
    )
    args = parser.parse_args()

    payload = run_audit(
        max_global_iter=int(args.max_global_iter),
        ladder_level=str(args.ladder_level),
    )
    write_json(args.output, payload)

    final = payload["final"]
    print(json.dumps(
        {
            "output": str(args.output),
            "split_rms_final": final["global"]["gas_phase_split_rms"],
            "top_beds": [
                {
                    "cell": r["cell_index"],
                    "type": r["cell_type"],
                    "split_rms": round(float(r["nr_split_rms"]), 5),
                    "var_frac": round(float(r["split_variance_fraction"]), 3),
                    "worst_species": (r["top_species"][0]["species"] if r["top_species"] else None),
                }
                for r in final["bed_rows"][:6]
            ],
            "t2_probes": payload["single_bed_t2_probes_final_state"],
        },
        indent=2,
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
