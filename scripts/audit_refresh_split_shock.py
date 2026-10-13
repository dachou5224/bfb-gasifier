#!/usr/bin/env python3
"""第一次 Vorab refresh 前后 phase-split 分解审计。

定位 init merit=0.104 → refresh 后 split≈0.223 的 cell/bed/物种来源。

用法：
  python scripts/audit_refresh_split_shock.py
  python scripts/audit_refresh_split_shock.py --output data/phase2_refresh_split_shock.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from _lu_harness import build_phase2_reactor, data_path, run_phase2_init, write_json

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity import (
    _measure_nr_init_scaled_metrics,
    run_outer_preinner_start_value_bridge,
)
from src.core.connectivity.preprojection import _bed_major_gas_phase_split_rms
from src.core.species import GAS_SPECIES
from src.solvers.global_nr_solver import (
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


def _species_split_rows(cell: Any, F_hat: np.ndarray, *, offset: int) -> list[dict[str, float | str]]:
    from src.solvers.nr_indexing import gas_mode_for_nr

    if gas_mode_for_nr(cell) != "two_phase":
        return []
    ng = n_gas_var(cell)
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
                "dense": d_val,
                "bubble": b_val,
                "dense_minus_bubble": float(d_val - b_val),
                "abs_split": float(abs(d_val - b_val)),
            }
        )
    rows.sort(key=lambda r: float(r["abs_split"]), reverse=True)
    return rows


def _audit_snapshot(reactor: Any, cells: list[Any], *, label: str) -> dict[str, Any]:
    reactor._apply_all_bc_for_nr()
    metrics = _measure_nr_init_scaled_metrics(
        cells, reactor.config, apply_bc_fn=reactor._apply_all_bc_for_nr
    )
    F_hat = _scaled_F_hat(reactor, cells)
    global_gas = _residual_gas_phase_metrics(F_hat, cells)
    global_group = _residual_group_metrics(F_hat, cells)

    offset = 0
    cell_rows: list[dict[str, Any]] = []
    split_sq_sum = 0.0
    for i, cell in enumerate(cells):
        nv = n_var(cell)
        block_offset = offset
        offset += nv
        from src.solvers.nr_indexing import gas_mode_for_nr

        if gas_mode_for_nr(cell) != "two_phase":
            continue
        ng = n_gas_var(cell)
        half = ng // 2
        dense = F_hat[block_offset : block_offset + half]
        bubble = F_hat[block_offset + half : block_offset + ng]
        pos = main_nr_gas_phase_split_positions_for_cell(cell)
        if pos.size == 0:
            continue
        terms = np.asarray(dense[pos] - bubble[pos], dtype=np.float64)
        cell_rms = float(np.sqrt(np.mean(terms**2)))
        split_sq_sum += float(np.sum(terms**2))
        cell_rows.append(
            {
                "cell_index": int(i),
                "cell_type": str(getattr(cell, "cell_type", "bed")),
                "T_K": float(getattr(cell, "T", float("nan"))),
                "split_rms": cell_rms,
                "split_max_abs": float(np.max(np.abs(terms))),
                "preproj_major_split_rms": float(
                    _bed_major_gas_phase_split_rms(cell, reaction_rate_multiplier=1.0)
                ),
                "n_split_terms": int(terms.size),
                "top_species": _species_split_rows(cell, F_hat, offset=block_offset)[:4],
            }
        )

    cell_rows.sort(key=lambda r: float(r["split_rms"]), reverse=True)
    for row in cell_rows:
        n = int(row["n_split_terms"])
        row["split_variance_fraction"] = float(
            (float(row["split_rms"]) ** 2) * n / max(split_sq_sum, 1.0e-30)
        )

    bed_only = [r for r in cell_rows if r["cell_type"] == "bed"]
    bed_agg: dict[str, float] = {}
    for name, indices in (
        ("bed0", [0]),
        ("bed1", list(range(1, 10))),
        ("bed2", list(range(10, 20))),
    ):
        subset = [r for r in bed_only if int(r["cell_index"]) in indices]
        if not subset:
            continue
        sq = sum(float(r["split_rms"]) ** 2 * int(r["n_split_terms"]) for r in subset)
        n = sum(int(r["n_split_terms"]) for r in subset)
        bed_agg[name] = float(np.sqrt(sq / max(n, 1)))

    return {
        "label": label,
        "metrics": metrics,
        "global_group": global_group,
        "global_gas": global_gas,
        "bed_aggregate_split_rms": bed_agg,
        "cell_rows": cell_rows,
        "vorab_signature": reactor._vorabrechnung_snapshot_signature(),
    }


def _delta_table(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    bm = before["metrics"]
    am = after["metrics"]
    out: dict[str, Any] = {
        "merit_delta": float(am["merit_scaled"] - bm["merit_scaled"]),
        "split_delta": float(am["gas_phase_split_rms"] - bm["gas_phase_split_rms"]),
        "energy_delta": float(am["energy_rms"] - bm["energy_rms"]),
        "rms_delta": float(am["rms_scaled"] - bm["rms_scaled"]),
        "bed_aggregate_delta": {
            k: float(after["bed_aggregate_split_rms"].get(k, 0.0) - before["bed_aggregate_split_rms"].get(k, 0.0))
            for k in sorted(set(before["bed_aggregate_split_rms"]) | set(after["bed_aggregate_split_rms"]))
        },
        "top_cell_deltas": [],
    }
    before_by_cell = {int(r["cell_index"]): r for r in before["cell_rows"]}
    deltas: list[dict[str, Any]] = []
    for row in after["cell_rows"]:
        ci = int(row["cell_index"])
        prev = before_by_cell.get(ci)
        if prev is None:
            continue
        d_split = float(row["split_rms"] - prev["split_rms"])
        if abs(d_split) < 1.0e-6:
            continue
        deltas.append(
            {
                "cell_index": ci,
                "cell_type": row["cell_type"],
                "split_rms_before": float(prev["split_rms"]),
                "split_rms_after": float(row["split_rms"]),
                "split_rms_delta": d_split,
                "worst_species_before": (prev["top_species"][0] if prev["top_species"] else None),
                "worst_species_after": (row["top_species"][0] if row["top_species"] else None),
            }
        )
    deltas.sort(key=lambda r: abs(float(r["split_rms_delta"])), reverse=True)
    out["top_cell_deltas"] = deltas[:12]
    return out


def run_audit(*, ladder_level: str) -> dict[str, Any]:
    reactor = build_phase2_reactor(ladder_level=ladder_level)
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    cells = reactor._solver_cells_for_nr()

    init_snap = _audit_snapshot(reactor, cells, label="init")

    reactor._refresh_vorabrechnung_sources_for_nr(True)
    refresh_snap = _audit_snapshot(reactor, cells, label="after_first_refresh")

    preinner_diag = run_outer_preinner_start_value_bridge(
        cells,
        reactor.config,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        top_above_cell=reactor.freeboard_cells[0] if reactor.freeboard_cells else None,
    )
    preinner_snap = _audit_snapshot(reactor, cells, label="after_preinner_bridge")

    return {
        "ladder_level": ladder_level,
        "init": init_snap,
        "after_first_refresh": refresh_snap,
        "after_preinner_bridge": preinner_snap,
        "init_to_refresh_delta": _delta_table(init_snap, refresh_snap),
        "refresh_to_preinner_delta": _delta_table(refresh_snap, preinner_snap),
        "preinner_bridge_diag": preinner_diag,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ladder-level", default="p2")
    parser.add_argument("--output", type=Path, default=data_path("phase2_refresh_split_shock.json"))
    args = parser.parse_args()

    payload = run_audit(ladder_level=str(args.ladder_level))
    write_json(args.output, payload)

    d = payload["init_to_refresh_delta"]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "init_merit": payload["init"]["metrics"]["merit_scaled"],
                "refresh_merit": payload["after_first_refresh"]["metrics"]["merit_scaled"],
                "refresh_split": payload["after_first_refresh"]["metrics"]["gas_phase_split_rms"],
                "merit_delta": d["merit_delta"],
                "split_delta": d["split_delta"],
                "bed_agg_refresh": payload["after_first_refresh"]["bed_aggregate_split_rms"],
                "top_cell_deltas": d["top_cell_deltas"][:6],
                "preinner": {
                    "global_gate_pass": payload["preinner_bridge_diag"].get("global_gate_pass"),
                    "split_before": payload["preinner_bridge_diag"].get("split_rms_before"),
                    "split_after": payload["preinner_bridge_diag"].get("split_rms_after"),
                    "bed0_accepted": payload["preinner_bridge_diag"].get("bed0_accepted"),
                },
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
