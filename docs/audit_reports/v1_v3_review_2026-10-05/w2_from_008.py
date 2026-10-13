"""暖启动已过门之后，从 s=0.08 的联合通过点再解一次 s=0.16。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    solve_named_controlled,
)
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_OLD = Path(__file__).resolve().parents[1] / "s1_s3_review_2026-09-30" / "t3_corrected_path.json"
_OUT = Path(__file__).resolve().parent / "w2_from_008.json"


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
    old = json.loads(_OLD.read_text())
    saved = [row for row in old["toward_022"] if abs(row["s"] - 0.08) < 1.0e-12 and row["joint_pass"]][-1]
    x = np.asarray(saved["x"], dtype=np.float64)
    out = solve_named_controlled(q_on_wesseling1_path(0.16), x, deadline_s=300.0)
    z = np.asarray(out["x"], dtype=np.float64)
    q = outer_P_named(z, q_on_wesseling1_path(0.16))
    f = coupled_residual_named(z, q)
    cert = certify_controlled_state(z, f, q.hydro)
    led = external_imbalance(z, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    cho_ok = max(abs(v) for v in cho.values()) < 1.0e-4
    report = {
        "start": "joint-pass s=0.08, target s=0.16",
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "passed_block_gates": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "cho_gap_mol_s": cho,
        "cho_gap_compatible": cho_ok,
        "joint_pass": bool(cert["passed"] and cert["transport_closed"] and cho_ok),
        "gate_ratio": cert["gate_ratio"],
        "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
        "T_bed": float(z[-2]),
        "T_fb": float(z[-1]),
        "history": out["stats"].get("history"),
        "calls": {
            "F": out["stats"].get("F_calls"),
            "J": out["stats"].get("J_calls"),
            "P": out["stats"].get("P_calls"),
            "outers": out["stats"].get("outers"),
            "direction": out["stats"].get("direction"),
            "log_span": out["stats"].get("log_span"),
        },
        "x": z.tolist(),
        "F": np.asarray(f, dtype=np.float64).tolist(),
    }
    print(report["joint_pass"], report["reason"], round(report["elapsed_s"], 1), report["gate_ratio"], flush=True)
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")


if __name__ == "__main__":
    main()
