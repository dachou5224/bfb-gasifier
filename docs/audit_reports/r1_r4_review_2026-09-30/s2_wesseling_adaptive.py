"""固定 Wesseling 气场上的轨迹验收。不求根，不把 15 s 续窗加长成停留时间。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.physics.wirsum_1997_freeboard import classify_refined_pair
from src.solvers.wirsum_1997_named_f_control import transport_is_closed
from src.solvers.wirsum_1997_named_f_p2 import LAUNCH_DT_S, finite_hydro_from_x
from src.solvers.wesseling_two_zone import (
    build_wesseling1_two_zone_q,
    wesseling1_physical_x0,
)

_OUT = Path(__file__).resolve().parent / "s2_wesseling_adaptive.json"
_RUNTIME_CAP_S = 20.0  # [s] 运行预算，不是停留时间


def _weights(hydro) -> dict:
    return {
        "w_leave": np.asarray(hydro.w_leave, dtype=float).tolist(),
        "w_return": np.asarray(hydro.w_return, dtype=float).tolist(),
        "w_undecided": np.asarray(hydro.w_undecided, dtype=float).tolist(),
        "tau_s": np.asarray(hydro.tau_s, dtype=float).tolist(),
        "occupancy_s": np.asarray(hydro.occupancy_s, dtype=float).tolist(),
    }


def _key(row: dict) -> tuple[float, float]:
    return (round(float(row["d_p_m"]), 12), round(float(row["u_launch_m_s"]), 8))


def _run(x, ref, dt: float, runtime_cap_s: float, *, one_cell: bool = False) -> dict:
    traces: list[dict] = []
    snap, _geo, _k = finite_hydro_from_x(
        x, ref,
        T_bed_K=float(x[-2]), T_fb_K=float(x[-1]), N_fb=x[16:23],
        dt=float(dt),
        H_cell=float(ref.H_freeboard) if one_cell else None,
        launch_traces=traces,
        runtime_cap_s=float(runtime_cap_s),
    )
    reasons: dict[str, int] = {}
    for row in traces:
        reasons[row["reason"]] = reasons.get(row["reason"], 0) + 1
    return {
        "dt_s": float(dt),
        "runtime_cap_s": float(runtime_cap_s),
        "one_cell": bool(one_cell),
        "transport_closed": transport_is_closed(snap),
        "weights": _weights(snap),
        "reason_counts": reasons,
        "n_traces": len(traces),
        "traces": traces,
    }


def _compare(coarse: dict, fine: dict) -> dict:
    left = {_key(row): row for row in coarse["traces"]}
    right = {_key(row): row for row in fine["traces"]}
    kinds: dict[str, int] = {}
    rows = []
    for key in sorted(left):
        a = left[key]
        b = right[key]
        kind = classify_refined_pair(
            (a["fate"], a["t_s"], a["height_m"], a["u_m_s"]),
            (b["fate"], b["t_s"], b["height_m"], b["u_m_s"]),
        )
        kinds[kind] = kinds.get(kind, 0) + 1
        rows.append({
            "d_p_m": key[0],
            "u_launch_m_s": key[1],
            "weight": a.get("weight"),
            "kind": kind,
            "reason_dt": a["reason"],
            "reason_dt_half": b["reason"],
            "t_s": a["t_s"],
            "t_half_s": b["t_s"],
        })
    return {"kind_counts": kinds, "pairs": rows}


def main() -> None:
    nq = build_wesseling1_two_zone_q()
    x = wesseling1_physical_x0(nq)
    ref = nq.reference
    fixed_two = _run(x, ref, LAUNCH_DT_S, 0.0)
    fixed_one = _run(x, ref, LAUNCH_DT_S, 0.0, one_cell=True)
    coarse_rows = [row for row in fixed_two["traces"] if abs(float(row["d_p_m"]) - 3.0e-4) < 1.0e-12]
    adaptive_dt = _run(x, ref, LAUNCH_DT_S, _RUNTIME_CAP_S)
    adaptive_half = _run(x, ref, 0.5 * LAUNCH_DT_S, _RUNTIME_CAP_S)
    compared = _compare(adaptive_dt, adaptive_half)
    report = {
        "continue_budget_s": float(ref.trajectory_continue_s),
        "runtime_cap_s": _RUNTIME_CAP_S,
        "fixed_window": {
            "two_cell": {k: fixed_two[k] for k in ("transport_closed", "weights", "reason_counts", "n_traces")},
            "one_cell": {k: fixed_one[k] for k in ("transport_closed", "weights", "reason_counts", "n_traces")},
            "coarse_d_p_3e-4": [
                {k: row[k] for k in ("fate", "reason", "t_s", "cell", "height_m", "u_m_s", "weight") if k in row}
                for row in coarse_rows
            ],
        },
        "adaptive": {
            "dt": {k: adaptive_dt[k] for k in ("transport_closed", "weights", "reason_counts", "n_traces")},
            "dt_half": {k: adaptive_half[k] for k in ("transport_closed", "weights", "reason_counts", "n_traces")},
            "step_agreement": compared,
        },
    }
    _OUT.write_text(json.dumps(report), encoding="utf-8")
    print("fixed two", fixed_two["reason_counts"], fixed_two["weights"]["w_undecided"], flush=True)
    print("fixed one", fixed_one["weights"]["w_undecided"], flush=True)
    print("adaptive", adaptive_dt["reason_counts"], adaptive_dt["transport_closed"], flush=True)
    print("half", adaptive_half["reason_counts"], adaptive_half["transport_closed"], flush=True)
    print("agree", compared["kind_counts"], flush=True)


if __name__ == "__main__":
    main()
