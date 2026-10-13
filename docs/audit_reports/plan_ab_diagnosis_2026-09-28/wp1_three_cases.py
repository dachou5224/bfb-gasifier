"""工作包 1 三例。基准只给初值。通过只认刷新 P 之后的分块门。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_named_f import default_x0_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import edges_from_x

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp1_three_cases.json"
_ENDPOINT = _ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz"


def _cold() -> np.ndarray:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = default_x0_named(nq)
    if float(x[42]) <= 0.0:
        x[42] = 1.2e-5
    if float(x[43]) <= 0.0:
        x[43] = 1.8e-5
    q = outer_P_named(x, nq)
    ed = edges_from_x(x[:39], q.energy.species)
    tau = np.maximum(np.asarray(q.hydro.tau_s, dtype=np.float64), 1.0e-6)
    x[46:48] = np.maximum(np.asarray(ed["m_eject"], dtype=np.float64) * tau, 1.0e-12)
    return x


def _row(name: str, out: dict) -> dict:
    x = np.asarray(out["x"], dtype=np.float64)
    return {
        "name": name,
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "merit": out["merit"],
        "stats": out["stats"],
        "blocks": out["blocks"],
        "gate_ratio": out["gate_ratio"],
        "opened_zeros": out["opened_zeros"],
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
        "worst_index": int(np.argmax(np.abs(out["F"]))),
        "worst_value": float(out["F"][int(np.argmax(np.abs(out["F"])))]),
    }


def main() -> None:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    base = np.load(_ENDPOINT)["x"].copy()
    hot = base.copy()
    hot[-2:] += 50.0
    half = base.copy()
    pos = half[:-2] > 0.0
    half[:-2][pos] *= 0.5
    cases = [
        ("T_plus_50K", hot, 120.0),
        ("positive_half", half, 120.0),
        ("cold_start", _cold(), 120.0),
    ]
    rows = []
    for name, x0, budget in cases:
        out = solve_named_controlled(nq, x0, deadline_s=budget)
        row = _row(name, out)
        rows.append(row)
        print(name, row["passed"], row["reason"], round(row["elapsed_s"], 1), flush=True)
    _OUT.write_text(json.dumps({"cases": rows}, indent=2), encoding="utf-8")
    print("wrote", _OUT, flush=True)


if __name__ == "__main__":
    main()
