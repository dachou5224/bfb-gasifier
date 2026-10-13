"""V1 固定点通过后：重验 s=0.08，再从 U3 的 s=0.16 失败点做一次不超过 320 s 的求根。"""

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
_U3 = Path(__file__).resolve().parents[1] / "t1_t4_review_2026-09-30" / "u3_path.json"
_OUT = Path(__file__).resolve().parent / "v2_once.json"


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


def _assess(s: float, x: np.ndarray) -> dict:
    nq = q_on_wesseling1_path(s)
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, f, q.hydro)
    led = external_imbalance(x, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    cho_ok = max(abs(v) for v in cho.values()) < 1.0e-4
    return {
        "s": s,
        "passed_block_gates": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "cho_gap_mol_s": cho,
        "cho_gap_compatible": cho_ok,
        "joint_pass": bool(cert["passed"] and cert["transport_closed"] and cho_ok),
        "gate_ratio": cert["gate_ratio"],
        "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
    }


def main() -> None:
    old = json.loads(_OLD.read_text())
    saved08 = [row for row in old["toward_022"] if abs(row["s"] - 0.08) < 1.0e-12 and row["joint_pass"]][-1]
    recheck = _assess(0.08, np.asarray(saved08["x"], dtype=np.float64))
    print("recheck", recheck["joint_pass"], recheck["gate_ratio"], flush=True)
    u3 = json.loads(_U3.read_text())["s_0_16_controlled"]
    report = {"s_0_08_recheck": recheck, "start": "U3 s=0.16 failure, not the s=0.08 state"}
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")
    if not recheck["joint_pass"]:
        print("s=0.08 no longer jointly passes", flush=True)
        return
    x = np.asarray(u3["x"], dtype=np.float64)
    out = solve_named_controlled(q_on_wesseling1_path(0.16), x, deadline_s=320.0)
    z = np.asarray(out["x"], dtype=np.float64)
    solved = _assess(0.16, z)
    solved["reason"] = out["reason"]
    solved["elapsed_s"] = out["elapsed_s"]
    solved["merit"] = out["merit"]
    stats = out["stats"]
    solved["calls"] = {
        "F": stats.get("F_calls"),
        "J": stats.get("J_calls"),
        "P": stats.get("P_calls"),
        "outers": stats.get("outers"),
        "direction": stats.get("direction"),
        "accepted_alpha": stats.get("accepted_alpha"),
        "material_ratios": stats.get("material_ratios"),
    }
    print("s=0.16", solved["joint_pass"], solved["reason"], round(solved["elapsed_s"], 1), solved["gate_ratio"], flush=True)
    report["s_0_16"] = solved
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")


if __name__ == "__main__":
    main()
