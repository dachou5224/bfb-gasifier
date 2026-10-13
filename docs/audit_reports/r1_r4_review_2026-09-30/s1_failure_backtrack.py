"""在保存的温度 +50 K 失败点上，检查同一方向的分数步能否降低复合残差。"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_control import composite_backtrack, residual_block_merit
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q
from src.solvers.wirsum_1997_named_f import coupled_residual_named

_ROOT = Path(__file__).resolve().parents[3]
_X = json.loads(
    (_ROOT / "docs/audit_reports/review_abcd_2026-09-28/r2_tplus50.json").read_text()
)["x"]
_OUT = Path(__file__).resolve().parent / "s1_failure_backtrack.json"


def main() -> None:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = np.asarray(_X, dtype=np.float64)
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q, rate_scale=1.0)
    best = residual_block_merit(f)
    stats: dict = {"P_calls": 1.0, "F_calls": 1.0}
    proposed = _e8_log_scaled_correct(
        q, x, max_nfev=8, max_j=4, ftol=1.0e-4, temperature_factor=1.5,
        log_span=4.0, stats=stats,
    )
    chosen, rows = composite_backtrack(nq, x, proposed, best, stats)
    report = {
        "merit_at_zero": best,
        "accepted_alpha": None if chosen is None else chosen["alpha"],
        "accepted_merit": None if chosen is None else chosen["merit"],
        "rows": rows,
        "checkpoint": stats.get("checkpoint"),
        "inner_status": stats.get("inner_status"),
    }
    _OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(report["accepted_alpha"], report["accepted_merit"], best, flush=True)


if __name__ == "__main__":
    main()
