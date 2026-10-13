"""温度 +50 K。自适应雅可比预算 4→8→16，保存完整状态。"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "docs/audit_reports/review_abcd_2026-09-28/r2_tplus50.json"
_X = np.load(_ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]


def main() -> None:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x0 = np.asarray(_X, dtype=np.float64).copy()
    x0[-2:] += 50.0
    out = solve_named_controlled(nq, x0, deadline_s=120.0)
    x = np.asarray(out["x"], dtype=np.float64)
    row = {
        "name": "T_plus_50K",
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "stats": out["stats"],
        "gate_ratio": out["gate_ratio"],
        "blocks": out["blocks"],
        "transport_closed": out["transport_closed"],
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
        "x": x.tolist(),
        "F": np.asarray(out["F"], dtype=np.float64).tolist(),
    }
    _OUT.write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(row["passed"], row["reason"], round(row["elapsed_s"], 1), row["stats"].get("direction"), row["stats"].get("j_budget"), flush=True)


if __name__ == "__main__":
    main()
