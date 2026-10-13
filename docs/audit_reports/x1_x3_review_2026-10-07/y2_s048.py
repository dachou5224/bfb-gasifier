"""从 s=0.48 的失败检查点做一次有预算恢复。不改参数步长。"""

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
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_PREV = Path(__file__).resolve().parents[1] / "w1_w3_review_2026-10-07" / "x2_s048.json"
_OUT = Path(__file__).resolve().parent / "y2_s048.json"


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
    prev = json.loads(_PREV.read_text())
    x = np.asarray(prev["x"], dtype=np.float64)
    out = solve_named_controlled(q_on_wesseling1_path(0.48), x, deadline_s=120.0)
    z = np.asarray(out["x"], dtype=np.float64)
    q = outer_P_named(z, q_on_wesseling1_path(0.48))
    f = coupled_residual_named(z, q)
    cert = certify_controlled_state(z, f, q.hydro)
    led = external_imbalance(z, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    report = {
        "s": 0.48,
        "start": "x2_s048 accepted state",
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "phi": progress_phi(f),
        "passed_block_gates": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "cho_gap_mol_s": cho,
        "joint_pass": bool(
            cert["passed"] and cert["transport_closed"] and max(abs(v) for v in cho.values()) < 1.0e-4
        ),
        "gate_ratio": cert["gate_ratio"],
        "history": out["stats"].get("history"),
        "calls": {
            "F": out["stats"].get("F_calls"),
            "J": out["stats"].get("J_calls"),
            "P": out["stats"].get("P_calls"),
            "direction": out["stats"].get("direction"),
        },
        "x": z.tolist(),
        "F": np.asarray(f, dtype=np.float64).tolist(),
    }
    print(report["joint_pass"], report["reason"], round(report["elapsed_s"], 1), report["phi"], report["gate_ratio"], flush=True)
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")


if __name__ == "__main__":
    main()
