"""同一雅可比上重解不同信赖域，再用完整 P 看预测和实际。不沿路径求根。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    material_gate_ratios,
    progress_phi,
    step_improves_progress,
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_X = Path(__file__).resolve().parents[1] / "t1_t4_review_2026-09-30" / "u3_path.json"
_OUT = Path(__file__).resolve().parent / "w1_fixed.json"


def _plain(obj):
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def main() -> None:
    saved = json.loads(_X.read_text())["s_0_16_controlled"]
    x = np.asarray(saved["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.16)
    q = outer_P_named(x, nq)
    f0 = coupled_residual_named(x, q)
    stats: dict = {}
    rows = []
    for radius in (0.5, 0.05, 0.005, 0.0005):
        proposed = _e8_log_scaled_correct(
            q, x, max_j=4, log_span=radius, stats=stats,
        )
        pred = np.asarray(stats["predicted_raw"], dtype=np.float64)
        qy = outer_P_named(proposed, nq)
        fy = coupled_residual_named(proposed, qy)
        phi0 = progress_phi(f0)
        row = {
            "radius": radius,
            "J_calls": stats.get("J_calls"),
            "fd_residual_calls": stats.get("fd_residual_calls"),
            "chain_residual_calls": stats.get("chain_residual_calls"),
            "bed_face_calls": stats.get("bed_face_calls"),
            "linear_status": stats.get("linear_status"),
            "active_bound_count": stats.get("active_bound_count"),
            "predicted_phi": progress_phi(pred),
            "actual_phi": progress_phi(fy),
            "predicted_drop": phi0 - progress_phi(pred),
            "actual_drop": phi0 - progress_phi(fy),
            "predicted_ratios": material_gate_ratios(pred),
            "actual_ratios": material_gate_ratios(fy),
            "energy_W": fy[-2:].tolist(),
            "accepted_by_progress": step_improves_progress(f0, fy),
            "transport_closed": bool(np.all(np.asarray(qy.hydro.w_undecided) <= 1.0e-12)),
        }
        drop = row["predicted_drop"]
        row["rho"] = None if drop <= 1.0e-12 else row["actual_drop"] / drop
        rows.append(row)
        print(radius, row["predicted_ratios"], row["actual_ratios"], row["rho"], row["accepted_by_progress"], flush=True)
    out = {
        "state": "U3 s=0.16, not the unsaved V2 stop",
        "entrance_ratios": material_gate_ratios(f0),
        "entrance_phi": progress_phi(f0),
        "rows": rows,
    }
    _OUT.write_text(json.dumps(_plain(out)), encoding="utf-8")


if __name__ == "__main__":
    main()
