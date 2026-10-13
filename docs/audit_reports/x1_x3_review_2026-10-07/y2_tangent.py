"""从 s=0.32 的根做雅可比切线预测，按预测门比值选一个参数步。"""

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
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_PREV = Path(__file__).resolve().parents[1] / "w1_w3_review_2026-10-07" / "x2_s032.json"
_OUT = Path(__file__).resolve().parent / "y2_tangent.json"


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


def _residual(s: float, x: np.ndarray) -> np.ndarray:
    q = outer_P_named(x, q_on_wesseling1_path(s))
    return coupled_residual_named(x, q)


def main() -> None:
    prev = json.loads(_PREV.read_text())
    if not prev["joint_pass"]:
        raise SystemExit("s=0.32 is not a joint pass")
    s0 = 0.32
    x = np.asarray(prev["x"], dtype=np.float64)
    h = 1.0e-4
    dF = (_residual(s0 + h, x) - _residual(s0 - h, x)) / (2.0 * h)
    q = outer_P_named(x, q_on_wesseling1_path(s0))
    stats: dict = {}
    rows = []
    chosen = None
    for ds in (0.16, 0.08, 0.04, 0.02, 0.01, 0.005, 0.002, 0.001):
        proposed = _e8_log_scaled_correct(
            q, x, log_span=0.5, ratio_span=1.0, stats=stats,
            residual_offset=dF * ds,
        )
        pred = np.asarray(stats["predicted_raw"], dtype=np.float64)
        row = {
            "ds": ds,
            "s": s0 + ds,
            "predicted_phi": progress_phi(pred),
            "predicted_energy_W": pred[-2:].tolist(),
            "v_abs_max": stats.get("v_abs_max"),
        }
        rows.append(row)
        print(ds, row["predicted_phi"], row["predicted_energy_W"], flush=True)
        energy_ok = max(abs(float(v)) for v in pred[-2:]) <= 1.0
        if chosen is None and energy_ok and row["predicted_phi"] <= 20.0:
            chosen = (ds, proposed.copy(), row)
    report = {
        "predictor": "jacobian_tangent",
        "s0": s0,
        "difference_h": h,
        "candidates": rows,
    }
    if chosen is None:
        report["taken"] = None
        print("no ds keeps predicted energy within 1 W and phi<=20", flush=True)
    else:
        ds, x_pred, row = chosen
        s = s0 + ds
        out = solve_named_controlled(q_on_wesseling1_path(s), np.asarray(x_pred), deadline_s=150.0)
        z = np.asarray(out["x"], dtype=np.float64)
        qq = outer_P_named(z, q_on_wesseling1_path(s))
        f = coupled_residual_named(z, qq)
        cert = certify_controlled_state(z, f, qq.hydro)
        led = external_imbalance(z, qq)
        cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
        report["taken"] = {
            "ds": ds,
            "s": s,
            "predicted_phi": row["predicted_phi"],
            "reason": out["reason"],
            "elapsed_s": out["elapsed_s"],
            "actual_phi": progress_phi(f),
            "passed_block_gates": bool(cert["passed"]),
            "transport_closed": bool(cert["transport_closed"]),
            "cho_gap_mol_s": cho,
            "joint_pass": bool(
                cert["passed"] and cert["transport_closed"] and max(abs(v) for v in cho.values()) < 1.0e-4
            ),
            "gate_ratio": cert["gate_ratio"],
            "calls": {
                "F": out["stats"].get("F_calls"),
                "J": out["stats"].get("J_calls"),
                "P": out["stats"].get("P_calls"),
            },
            "x": z.tolist(),
            "F": np.asarray(f, dtype=np.float64).tolist(),
        }
        print("taken", ds, report["taken"]["joint_pass"], report["taken"]["reason"], round(report["taken"]["elapsed_s"], 1), flush=True)
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")


if __name__ == "__main__":
    main()
