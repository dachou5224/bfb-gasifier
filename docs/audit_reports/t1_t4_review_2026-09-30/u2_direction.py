"""s=0.16 上，行变换后的一次校正方向对完整残差是否下降。不跑长路径。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import (
    coupled_residual_named,
    inner_hydro_frozen_residual,
    outer_P_named,
)
from src.solvers.wirsum_1997_named_f_control import _block_report
from src.solvers.wirsum_1997_named_f_e8 import (
    _e8_log_scaled_correct,
    residual_block_merit,
)
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_X = Path(__file__).resolve().parents[1] / "s1_s3_review_2026-09-30" / "t3_corrected_path.json"
_OUT = Path(__file__).resolve().parent / "u2_direction.json"


def main() -> None:
    doc = json.loads(_X.read_text())
    saved = [item for item in doc["toward_022"] if abs(item["s"] - 0.16) < 1.0e-12 and not item["passed_block_gates"]][-1]
    x = np.asarray(saved["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.16)
    q = outer_P_named(x, nq)
    stats: dict = {}
    proposed = _e8_log_scaled_correct(
        q, x, max_nfev=8, max_j=4, ftol=1.0e-4, temperature_factor=1.5, log_span=4.0, stats=stats,
    )
    step = proposed - x
    rows = []
    for alpha in (0.0, 0.001, 0.01, 0.1, 1.0):
        y = x + alpha * step
        f_inner = inner_hydro_frozen_residual(y, q)
        qy = outer_P_named(y, nq)
        f_full = coupled_residual_named(y, qy)
        led = external_imbalance(y, qy)
        rows.append({
            "alpha": alpha,
            "inner_merit": residual_block_merit(f_inner),
            "full_merit": residual_block_merit(f_full),
            "blocks": _block_report(f_full),
            "cho": {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")},
            "transport_closed": bool(np.all(np.asarray(qy.hydro.w_undecided) <= 1.0e-12)),
        })
        print(alpha, rows[-1]["inner_merit"], rows[-1]["full_merit"], rows[-1]["blocks"]["tar_l2"], flush=True)
    base = rows[0]["full_merit"]
    small = [r for r in rows if 0.0 < r["alpha"] <= 0.1]
    descent = [r for r in small if r["full_merit"] < base]
    out = {
        "s": 0.16,
        "inner_status": stats.get("inner_status"),
        "F_calls": stats.get("F_calls"),
        "J_calls": stats.get("J_calls"),
        "rows": rows,
        "small_step_full_descent": bool(descent),
        "best_small_alpha": None if not descent else min(descent, key=lambda r: r["full_merit"])["alpha"],
    }
    _OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("descent", out["small_step_full_descent"], out["best_small_alpha"], flush=True)


if __name__ == "__main__":
    main()
