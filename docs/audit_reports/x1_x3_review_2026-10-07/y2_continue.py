"""从切线得到的 s=0.325 继续。每步重算切线，取预测能量仍在 1 W 内的最大步长。"""

from __future__ import annotations

import hashlib
import json
import time
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

_START = Path(__file__).resolve().parent / "y2_tangent.json"
_OUT = Path(__file__).resolve().parent / "y2_continue.json"
_ROOT = Path(__file__).resolve().parents[3]
_SOURCES = (
    "src/solvers/wirsum_1997_named_f.py",
    "src/solvers/wirsum_1997_named_f_e8.py",
    "src/solvers/wirsum_1997_named_f_control.py",
    "src/solvers/wesseling_case_contract.py",
    "src/solvers/wesseling_two_zone.py",
    "data/validation_cases.json",
)
_CANDIDATES = (0.02, 0.01, 0.005, 0.002)


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


def _source() -> str:
    digest = hashlib.sha256()
    for rel in _SOURCES:
        digest.update(rel.encode())
        digest.update((_ROOT / rel).read_bytes())
    return digest.hexdigest()


def _residual(s: float, x: np.ndarray) -> np.ndarray:
    return coupled_residual_named(x, outer_P_named(x, q_on_wesseling1_path(s)))


def _inputs(s: float) -> dict:
    q = q_on_wesseling1_path(s)
    feed = q.energy.feed
    return {
        "s": s,
        "steam_kg_h_assumed": 0.0,
        "Q_W": float(q.Q_W),
        "h_f_dry_J_kg": float(feed.h_f_dry_J_kg),
        "conditional_case": True,
    }


def main() -> None:
    start = json.loads(_START.read_text())["taken"]
    if not start["joint_pass"]:
        raise SystemExit("s=0.325 tangent step is not a joint pass")
    s = float(start["s"])
    x = np.asarray(start["x"], dtype=np.float64)
    source = _source()
    points = []
    failure = None
    t0 = time.monotonic()
    budget = 600.0
    while s < 1.0 - 1.0e-12 and time.monotonic() - t0 < budget:
        h = 1.0e-4
        dF = (_residual(s + h, x) - _residual(s - h, x)) / (2.0 * h)
        q = outer_P_named(x, q_on_wesseling1_path(s))
        stats: dict = {}
        chosen = None
        scan = []
        for ds in _CANDIDATES:
            if s + ds > 1.0:
                ds = 1.0 - s
            proposed = _e8_log_scaled_correct(
                q, x, log_span=0.5, ratio_span=1.0, stats=stats, residual_offset=dF * ds,
            )
            pred = np.asarray(stats["predicted_raw"], dtype=np.float64)
            row = {
                "ds": ds,
                "predicted_phi": progress_phi(pred),
                "predicted_energy_W": pred[-2:].tolist(),
            }
            scan.append(row)
            energy_ok = max(abs(float(v)) for v in pred[-2:]) <= 1.0
            if chosen is None and energy_ok and row["predicted_phi"] <= 20.0:
                chosen = (ds, np.asarray(proposed, dtype=np.float64).copy(), row)
            if s + ds >= 1.0:
                break
        if chosen is None:
            failure = {"s": s, "reason": "no_predicted_step", "scan": scan}
            print("stop", s, "no predicted step", flush=True)
            break
        ds, x_pred, row = chosen
        target = s + ds
        t_step = time.monotonic()
        out = solve_named_controlled(q_on_wesseling1_path(target), x_pred, deadline_s=60.0)
        z = np.asarray(out["x"], dtype=np.float64)
        qq = outer_P_named(z, q_on_wesseling1_path(target))
        f = coupled_residual_named(z, qq)
        cert = certify_controlled_state(z, f, qq.hydro)
        led = external_imbalance(z, qq)
        cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
        joint = bool(
            cert["passed"] and cert["transport_closed"] and max(abs(v) for v in cho.values()) < 1.0e-4
        )
        point = {
            "s": target,
            "ds": ds,
            "predictor": "jacobian_tangent",
            "predicted_phi": row["predicted_phi"],
            "predicted_energy_W": row["predicted_energy_W"],
            "reason": out["reason"],
            "elapsed_s": time.monotonic() - t_step,
            "actual_phi": progress_phi(f),
            "passed_block_gates": bool(cert["passed"]),
            "transport_closed": bool(cert["transport_closed"]),
            "cho_gap_mol_s": cho,
            "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
            "joint_pass": joint,
            "gate_ratio": cert["gate_ratio"],
            "inputs": _inputs(target),
            "source_sha256": source,
            "calls": {
                "F": out["stats"].get("F_calls"),
                "J": out["stats"].get("J_calls"),
                "P": out["stats"].get("P_calls"),
            },
            "x": z.tolist(),
            "F": np.asarray(f, dtype=np.float64).tolist(),
        }
        print(target, joint, out["reason"], round(point["elapsed_s"], 1), ds, flush=True)
        if not joint:
            failure = point
            break
        points.append(point)
        s = target
        x = z
        _OUT.write_text(json.dumps(_plain({
            "source_sha256": source,
            "wall_s": time.monotonic() - t0,
            "last_joint_pass_s": s,
            "n_steps": len(points),
            "points": points,
            "failure": failure,
            "reached_s1": False,
        })), encoding="utf-8")
    report = {
        "source_sha256": source,
        "wall_s": time.monotonic() - t0,
        "last_joint_pass_s": s if points else float(start["s"]),
        "n_steps": len(points),
        "points": points,
        "failure": failure,
        "reached_s1": bool(points) and abs(points[-1]["s"] - 1.0) <= 1.0e-12,
    }
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")
    print("end", report["last_joint_pass_s"], report["n_steps"], round(report["wall_s"], 1), flush=True)


if __name__ == "__main__":
    main()
