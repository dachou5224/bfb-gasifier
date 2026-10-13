"""固定 s=0.16 失败点上只跑一次原校正器。不启动长路径。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import coupled_residual_named, inner_hydro_frozen_residual, outer_P_named
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, residual_block_merit
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_X = Path(__file__).resolve().parents[1] / "s1_s3_review_2026-09-30" / "t3_corrected_path.json"
_OUT = Path(__file__).resolve().parent / "u1_s016_once.json"


def main() -> None:
    doc = json.loads(_X.read_text())
    row = [item for item in doc["toward_022"] if abs(item["s"] - 0.16) < 1.0e-12 and not item["passed_block_gates"]][-1]
    x = np.asarray(row["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.16)
    q = outer_P_named(x, nq)
    entrance = residual_block_merit(inner_hydro_frozen_residual(x, q))
    stats: dict = {}
    got = _e8_log_scaled_correct(
        q, x, max_nfev=8, max_j=4, ftol=1.0e-4, temperature_factor=1.5, log_span=4.0, stats=stats,
    )
    inner = residual_block_merit(inner_hydro_frozen_residual(got, q))
    full_q = outer_P_named(got, nq)
    full = residual_block_merit(coupled_residual_named(got, full_q))
    out = {
        "s": 0.16,
        "inner_status": stats.get("inner_status"),
        "checkpoint": stats.get("checkpoint"),
        "F_calls": stats.get("F_calls"),
        "J_calls": stats.get("J_calls"),
        "inner_merit_entrance": entrance,
        "inner_merit_best": stats.get("inner_merit_best"),
        "inner_merit_returned": inner,
        "full_merit_returned": full,
        "returned_inner_not_worse": bool(inner <= entrance + 1.0e-9),
        "outer_gate_unchanged": "full P still required; this record is not an accepted root",
    }
    _OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
