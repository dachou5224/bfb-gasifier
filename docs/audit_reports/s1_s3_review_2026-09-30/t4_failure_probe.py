"""修正后 s=0.16 失败点的分块尺度和焦油列。一次内层校正，不重跑 320 s。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import _block_report
from src.solvers.wirsum_1997_named_f_e8 import (
    _e8_log_scaled_correct,
    residual_block_merit,
    residual_block_scale,
)
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_PATH = Path(__file__).resolve().parent / "t3_corrected_path.json"
_OUT = Path(__file__).resolve().parent / "t4_failure_probe.json"


def _tar_ids(n: int) -> list[int]:
    return [i for i in range(40, 46) if i < n]


def main() -> None:
    doc = json.loads(_PATH.read_text())
    row = [item for item in doc["toward_022"] if abs(item["s"] - 0.16) < 1.0e-12 and not item["passed_block_gates"]][-1]
    x = np.asarray(row["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.16)
    q = outer_P_named(x, nq)
    f0 = coupled_residual_named(x, q)
    scale = residual_block_scale(int(f0.size))
    blocks0 = _block_report(f0)
    tar_ids = _tar_ids(int(x.size))
    columns = []
    for i in tar_ids + [0, 5, 12, 21]:
        if x[i] <= 0.0:
            continue
        xp = x.copy()
        xm = x.copy()
        xp[i] = x[i] * np.exp(1.0e-4)
        xm[i] = x[i] * np.exp(-1.0e-4)
        col = (coupled_residual_named(xp, q) - coupled_residual_named(xm, q)) / (2.0e-4 * scale)
        columns.append({
            "index": int(i),
            "x": float(x[i]),
            "F": float(f0[i]),
            "scaled_column_norm": float(np.linalg.norm(col)),
            "own_row": float(col[i]),
        })
    stats: dict = {}
    proposed = _e8_log_scaled_correct(
        q, x, max_nfev=8, max_j=4, ftol=1.0e-4, temperature_factor=1.5, log_span=4.0, stats=stats,
    )
    f_frozen = coupled_residual_named(proposed, q)
    q2 = outer_P_named(proposed, nq)
    f_new = coupled_residual_named(proposed, q2)
    out = {
        "s": 0.16,
        "blocks0": blocks0,
        "merit0": residual_block_merit(f0),
        "tar_rows": [float(f0[i]) for i in tar_ids],
        "tar_x": [float(x[i]) for i in tar_ids],
        "columns": columns,
        "corrector": {
            "inner_status": stats.get("inner_status"),
            "F_calls": stats.get("F_calls"),
            "J_calls": stats.get("J_calls"),
            "frozen_merit": residual_block_merit(f_frozen),
            "composite_merit": residual_block_merit(f_new),
            "frozen_blocks": _block_report(f_frozen),
            "composite_blocks": _block_report(f_new),
        },
    }
    _OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("blocks0", "merit0", "tar_rows", "columns", "corrector")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
