#!/usr/bin/env python3
"""bed1 固相输运链快照：init 末 → outer1 refresh → pre-inner 前后。

用法：
  python scripts/audit_bed1_solid_transport_chain.py
  python scripts/audit_bed1_solid_transport_chain.py --output data/bed1_chain.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from _lu_harness import build_phase2_reactor, data_path, run_phase2_init, write_json

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.cell import S_ASH, S_CHAR
from src.core.connectivity import (
    max_upper_bed_solid_auf_flux_ratio_error,
    run_outer_preinner_start_value_bridge,
)
from src.core.preprojection_engine import preproj_cell_energy_residual_abs


def _bed_transport_row(reactor: Any, bed_index: int = 1) -> dict[str, float]:
    bed_cells = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    if bed_index >= len(bed_cells):
        return {"bed": float(bed_index), "missing": 1.0}
    cell = bed_cells[bed_index]
    below = bed_cells[bed_index - 1] if bed_index > 0 else None
    auf_in = float(cell.m_solid_auf_in[:, (S_CHAR, S_ASH)].sum())
    out_auf = float(
        (cell.K_solid_auf[:, (S_CHAR, S_ASH)] * cell.m_solid[:, (S_CHAR, S_ASH)]).sum()
    )
    out_total = float(cell._solid_outflow_rates()[:, (S_CHAR, S_ASH)].sum())
    below_out = (
        float(
            (below.K_solid_auf[:, (S_CHAR, S_ASH)] * below.m_solid[:, (S_CHAR, S_ASH)]).sum()
        )
        if below is not None
        else 0.0
    )
    return {
        "bed": float(bed_index),
        "T_K": float(cell.T),
        "energy_residual_MW": float(preproj_cell_energy_residual_abs(cell) / 1e6),
        "char_ash_holdup_kg": float(cell.m_solid[:, (S_CHAR, S_ASH)].sum()),
        "solid_auf_in_kg_s": auf_in,
        "solid_out_auf_kg_s": out_auf,
        "solid_out_total_kg_s": out_total,
        "auf_in_over_out_auf": auf_in / max(out_auf, 1e-12),
        "below_out_auf_kg_s": below_out,
        "below_matches_auf_in": float(abs(below_out - auf_in) <= max(1e-9, 1e-6 * max(below_out, auf_in, 1.0))),
        "flux_ratio_error": float(max_upper_bed_solid_auf_flux_ratio_error(bed_cells, min_bed_index=1)),
    }


def run_probe(*, output: Path | None = None) -> dict[str, Any]:
    reactor = build_phase2_reactor(ladder_level="p2")
    run_phase2_init(reactor, init_strategy="vorabrechnung")

    stages: dict[str, dict[str, float]] = {}
    stages["init_end"] = _bed_transport_row(reactor, bed_index=1)

    reactor._refresh_vorabrechnung_sources_for_nr(force=True)
    stages["outer1_post_refresh"] = _bed_transport_row(reactor, bed_index=1)

    solver_cells = reactor._solver_cells_for_nr()
    pre_diag = run_outer_preinner_start_value_bridge(
        solver_cells,
        reactor.config,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        top_above_cell=reactor.freeboard_cells[0] if reactor.freeboard_cells else None,
    )
    stages["outer1_post_preinner"] = _bed_transport_row(reactor, bed_index=1)

    payload: dict[str, Any] = {
        "stages": stages,
        "preinner_diag": {
            "holdup_reconciled": bool(pre_diag.get("holdup_reconciled", False)),
            "holdup_flux_mismatch_before": float(pre_diag.get("holdup_flux_mismatch_before", 0.0)),
            "split_rms_before": float(pre_diag.get("split_rms_before", 0.0)),
            "split_rms_after": float(pre_diag.get("split_rms_after", 0.0)),
            "merit_before": float(pre_diag.get("merit_before", 0.0)),
            "merit_after": float(pre_diag.get("merit_after", 0.0)),
            "tier_plan": pre_diag.get("tier_plan"),
        },
    }
    if output is not None:
        write_json(output, payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=data_path("bed1_solid_transport_chain.json"))
    args = parser.parse_args()
    payload = run_probe(output=args.output)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
