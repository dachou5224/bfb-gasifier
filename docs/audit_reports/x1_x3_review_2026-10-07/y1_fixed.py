"""s=0.48 失败点上，门尺度子问题是否给出与预测一致的下降。不跑长路径。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    progress_phi,
    solve_named_controlled,
    step_improves_progress,
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_ROOT = Path(__file__).resolve().parents[1] / "w1_w3_review_2026-10-07"
_OUT = Path(__file__).resolve().parent / "y1_fixed.json"


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
    failed = json.loads((_ROOT / "x2_s048.json").read_text())
    passed = json.loads((_ROOT / "x2_s032.json").read_text())
    x = np.asarray(failed["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.48)
    q = outer_P_named(x, nq)
    f0 = coupled_residual_named(x, q)
    stats: dict = {}
    proposed = _e8_log_scaled_correct(
        q, x, log_span=0.001953125, ratio_span=0.001953125, stats=stats,
    )
    pred = np.asarray(stats["predicted_raw"], dtype=np.float64)
    qy = outer_P_named(proposed, nq)
    fy = coupled_residual_named(proposed, qy)
    x32 = np.asarray(passed["x"], dtype=np.float64)
    q32 = outer_P_named(x32, q_on_wesseling1_path(0.32))
    f32 = coupled_residual_named(x32, q32)
    cert32 = certify_controlled_state(x32, f32, q32.hydro)
    led32 = external_imbalance(x32, q32)
    recover_x = np.asarray(
        json.loads(
            (Path(__file__).resolve().parents[1] / "v1_v3_review_2026-10-05" / "w2_from_008.json").read_text()
        )["x"],
        dtype=np.float64,
    )
    recovered = solve_named_controlled(
        q_on_wesseling1_path(0.16), recover_x, deadline_s=30.0,
    )
    z = np.asarray(recovered["x"], dtype=np.float64)
    qr = outer_P_named(z, q_on_wesseling1_path(0.16))
    fr = coupled_residual_named(z, qr)
    certr = certify_controlled_state(z, fr, qr.hydro)
    out = {
        "s048": {
            "phi0": progress_phi(f0),
            "predicted_phi": progress_phi(pred),
            "actual_phi": progress_phi(fy),
            "energy_W": fy[-2:].tolist(),
            "accepted_by_progress": step_improves_progress(f0, fy),
            "prediction_gap": abs(progress_phi(pred) - progress_phi(fy)),
            "objective": stats.get("objective"),
            "linear_status": stats.get("linear_status"),
        },
        "s032_still_passes": bool(cert32["passed"]),
        "s032_cho": {k: float(led32["gap_mol_s"][k]) for k in ("C", "H", "O")},
        "s016_ratio_recovery_passes": bool(certr["passed"]),
        "s016_elapsed_s": recovered["elapsed_s"],
    }
    print(json.dumps(_plain(out), indent=2), flush=True)
    _OUT.write_text(json.dumps(_plain(out)), encoding="utf-8")


if __name__ == "__main__":
    main()
