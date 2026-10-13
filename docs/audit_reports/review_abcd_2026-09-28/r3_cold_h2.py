"""不读端点的冷初值。检查气泡氢是否还能被锁在零。"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import edges_from_x
from src.solvers.wirsum_1997_named_f import default_x0_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

_OUT = Path(__file__).resolve().parent / "r3_cold_h2.json"


def main() -> None:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = default_x0_named(nq)
    if float(x[42]) <= 0.0:
        x[42] = 1.2e-5
    if float(x[43]) <= 0.0:
        x[43] = 1.8e-5
    q = outer_P_named(x, nq)
    ed = edges_from_x(x[:39], q.energy.species)
    tau = np.asarray(q.hydro.tau_s, dtype=np.float64)
    und = np.asarray(q.hydro.w_undecided, dtype=np.float64)
    tau = np.where(und <= 1.0e-12, np.maximum(tau, 1.0e-6), 0.0)
    x[46:48] = np.maximum(np.asarray(ed["m_eject"], dtype=np.float64) * tau, 0.0)
    out = solve_named_controlled(nq, x, deadline_s=90.0)
    got = np.asarray(out["x"], dtype=np.float64)
    row = {
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "bubble_h2": float(got[9]),
        "linear_zero_ids": out["stats"].get("linear_zero_ids", []),
        "h2_locked_at_zero": float(got[9]) == 0.0,
        "T_bed": float(got[-2]),
        "T_fb": float(got[-1]),
        "gate_ratio": out["gate_ratio"],
        "stats": {k: v for k, v in out["stats"].items() if k != "linear_zero_ids"},
        "x": got.tolist(),
    }
    _OUT.write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(row["passed"], row["reason"], round(row["elapsed_s"], 1), "h2", row["bubble_h2"], flush=True)


if __name__ == "__main__":
    main()
