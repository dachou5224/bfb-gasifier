"""从已认证的 s=0.08 沿同一路径加倍前进。失败则步长减半，小于 0.01 停止。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    solve_named_controlled,
)
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_PATH = Path(__file__).resolve().parent / "s3_path_step.json"


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
    doc = json.loads(_PATH.read_text())
    rows = list(doc.get("continued", []))
    passed_rows = [row for row in rows if row.get("passed")]
    if passed_rows:
        current = float(passed_rows[-1]["s"])
        x = np.asarray(passed_rows[-1]["x"], dtype=np.float64)
        step = 0.08
    else:
        current = 0.08
        step = 0.08
        x = np.asarray(doc["next2"]["x"], dtype=np.float64)
    t0 = time.monotonic()
    while current < 1.0 - 1.0e-12 and step >= 0.01 and time.monotonic() - t0 < 900.0:
        trial = min(1.0, current + step)
        nq = q_on_wesseling1_path(trial)
        out = solve_named_controlled(nq, x, deadline_s=320.0)
        z = np.asarray(out["x"], dtype=np.float64)
        q = outer_P_named(z, nq)
        f = coupled_residual_named(z, q)
        cert = certify_controlled_state(z, f, q.hydro)
        passed = bool(cert["passed"])
        led = external_imbalance(z, q)
        row = {
            "name": "not CASE_HTW_WESSELING_1" if trial < 1.0 - 1.0e-12 else "CASE_HTW_WESSELING_1",
            "s": trial,
            "target_reached": bool(passed and trial >= 1.0 - 1.0e-12),
            "passed": passed,
            "solver_passed": bool(out["passed"]),
            "reason": "block_gates" if passed else out["reason"],
            "elapsed_s": out["elapsed_s"],
            "gate_ratio": cert["gate_ratio"],
            "T_bed": float(z[-2]),
            "T_fb": float(z[-1]),
            "ledger_gap_mol_s": {k: float(v) for k, v in led["gap_mol_s"].items()},
            "ledger_enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
            "x": z.tolist(),
            "F": np.asarray(f, dtype=np.float64).tolist(),
        }
        rows.append(row)
        doc["continued"] = rows
        doc["certified_s"] = current if not passed else trial
        _PATH.write_text(json.dumps(_plain(doc)), encoding="utf-8")
        print(trial, passed, row["reason"], round(out["elapsed_s"], 1), cert["gate_ratio"], flush=True)
        if passed:
            x = z
            current = trial
            step = min(1.0 - current, step * 2.0)
            if current >= 1.0 - 1.0e-12:
                break
        else:
            step *= 0.5
    print("certified_s", current, "next_step", step, flush=True)


if __name__ == "__main__":
    main()
