"""Wesseling 1 各抛射轨迹的末态。不加大续积分预算，不求根。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_control import transport_is_closed
from src.solvers.wirsum_1997_named_f_p2 import finite_hydro_from_x
from src.solvers.wesseling_two_zone import (
    build_wesseling1_two_zone_q,
    wesseling1_physical_x0,
)

_OUT = Path(__file__).resolve().parent / "r1_wesseling_traces.json"


def main() -> None:
    nq = build_wesseling1_two_zone_q()
    x = wesseling1_physical_x0(nq)
    q = outer_P_named(x, nq)
    traces: list[dict] = []
    finite_hydro_from_x(
        x, nq.reference,
        T_bed_K=float(x[-2]), T_fb_K=float(x[-1]), N_fb=x[16:23],
        launch_traces=traces,
    )
    reasons: dict[str, int] = {}
    for row in traces:
        reasons[row["reason"]] = reasons.get(row["reason"], 0) + 1
    report = {
        "continue_budget_s": float(nq.reference.trajectory_continue_s),
        "transport_closed": transport_is_closed(q.hydro),
        "w_undecided": np.asarray(q.hydro.w_undecided, dtype=float).tolist(),
        "n_traces": len(traces),
        "reason_counts": reasons,
        "traces": traces,
    }
    _OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(report["transport_closed"], report["w_undecided"], reasons, flush=True)


if __name__ == "__main__":
    main()
