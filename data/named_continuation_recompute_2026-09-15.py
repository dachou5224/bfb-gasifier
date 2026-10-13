"""Walk the named continuation chain once. Not a production entry."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import freeze_stage_q, solve_inner_named
from src.solvers.wirsum_1997_named_f_continuation import certify_named_full_rate_endpoint
from src.solvers.wirsum_1997_named_f_e3 import freeze_e3_q, solve_named_e3_from_r234
from src.solvers.wirsum_1997_named_f_e4 import freeze_e4_q, solve_named_e4_from_e3
from src.solvers.wirsum_1997_named_f_e5 import freeze_e5_q, solve_named_e5_from_e4
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q, solve_named_e6_from_e5
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q, solve_named_e6_then_r10
from src.solvers.wirsum_1997_named_f_r12 import freeze_r12_only_q, solve_named_r12_only
from src.solvers.wirsum_1997_named_f_r12_r5 import freeze_r12_r5_q, solve_named_r12_then_r5
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import (
    freeze_r12_r5_r8_q,
    solve_named_r12_r5_then_r8,
)
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import (
    freeze_r12_r5_r8_r6_q,
    solve_named_r12_r5_r8_then_r6,
)
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import (
    freeze_r12_r5_r8_r6_r7_q,
    solve_named_r12_r5_r8_r6_then_r7,
)
from src.solvers.wirsum_1997_named_f_r234 import freeze_homog_r234_q, solve_named_homog_then_r234
from src.solvers.wirsum_1997_stiffness_full_gas import IDX

OUT = Path("data/named_continuation_recompute_2026-09-15.json")


def _hist(sol: dict) -> list[dict]:
    rows = []
    for row in sol.get("continuation") or []:
        item = {}
        for k, v in row.items():
            if isinstance(v, (str, bool)) or v is None:
                item[k] = v
            elif isinstance(v, (int, float, np.integer, np.floating)):
                item[k] = float(v)
        rows.append(item)
    return rows


def summarize(name: str, sol: dict, nq, coupling: str, wall_s: float) -> dict:
    x = np.asarray(sol["x"], dtype=np.float64)
    cert = certify_named_full_rate_endpoint(x, nq, coupling=coupling)
    cert_out = {k: v for k, v in cert.items() if k != "F"}
    n = int(x.size)
    out = {
        "name": name,
        "f_name": sol.get("f_name"),
        "enabled_rxns": list(sol.get("enabled_rxns") or []),
        "n_unknowns": n,
        "T_K": float(x[-1]),
        "residual_l2": float(sol.get("residual_l2", np.linalg.norm(sol["F"]))),
        "physical_convergence": bool(sol.get("physical_convergence", False)),
        "solver_stop_reason": sol.get("solver_stop_reason", sol.get("failure_reason")),
        "advance_mode": sol.get("advance_mode"),
        "coupling_mode": sol.get("coupling_mode", coupling),
        "branch_tracked": sol.get("branch_tracked"),
        "endpoint_certified": bool(cert_out["endpoint_certified"]),
        "cert": cert_out,
        "N_d_H2": float(x[IDX["H2"]]) if n > IDX["H2"] else None,
        "N_d_O2": float(x[IDX["O2"]]) if n > IDX["O2"] else None,
        "N_d_CO": float(x[IDX["CO"]]) if n > IDX["CO"] else None,
        "N_d_CH4": float(x[IDX["CH4"]]) if n > IDX["CH4"] else None,
        "N_b_O2": float(x[7 + IDX["O2"]]) if n > 7 + IDX["O2"] else None,
        "N_tar": float(x[40] + x[41]) if n >= 43 else None,
        "wall_s": float(wall_s),
        "continuation": _hist(sol),
    }
    print(
        f"{name:12s} T={out['T_K']:8.2f} K  "
        f"|F|={out['residual_l2']:.3e}  "
        f"stop={out['solver_stop_reason']}  "
        f"phys={out['physical_convergence']}  "
        f"cert={out['endpoint_certified']}  "
        f"{wall_s:.1f}s",
        flush=True,
    )
    return out


def timed(fn):
    t0 = time.perf_counter()
    sol = fn()
    return sol, time.perf_counter() - t0


def main() -> None:
    rows = []
    t_all = time.perf_counter()

    e0, dt = timed(lambda: solve_inner_named(freeze_stage_q(0), max_iter=40))
    rows.append(summarize("F_E0", e0, freeze_stage_q(0), "frozen_q", dt))

    r12, dt = timed(lambda: solve_named_r12_only(max_iter=20))
    rows.append(summarize("R12", r12, freeze_r12_only_q(), "frozen_q", dt))

    r5, dt = timed(
        lambda: solve_named_r12_then_r5(max_iter=20, apply_startwert=True, r12=r12)
    )
    rows.append(summarize("R12+R5", r5, freeze_r12_r5_q(), "frozen_q", dt))

    r8, dt = timed(lambda: solve_named_r12_r5_then_r8(max_iter=20, r5=r5))
    rows.append(summarize("R12+R5+R8", r8, freeze_r12_r5_r8_q(), "frozen_q", dt))

    r6, dt = timed(
        lambda: solve_named_r12_r5_r8_then_r6(
            max_iter=20, apply_startwert=True, r8=r8
        )
    )
    rows.append(summarize("R12+R5+R8+R6", r6, freeze_r12_r5_r8_r6_q(), "frozen_q", dt))

    r7, dt = timed(
        lambda: solve_named_r12_r5_r8_r6_then_r7(max_iter=20, r6=r6)
    )
    rows.append(
        summarize("homog+R7", r7, freeze_r12_r5_r8_r6_r7_q(), "frozen_q", dt)
    )

    r234, dt = timed(lambda: solve_named_homog_then_r234(max_iter=20, homog=r7))
    rows.append(summarize("R2-R4", r234, freeze_homog_r234_q(), "frozen_q", dt))

    e3, dt = timed(lambda: solve_named_e3_from_r234(max_iter=20, r234=r234))
    rows.append(summarize("F_E3", e3, freeze_e3_q(), "frozen_q", dt))

    e4, dt = timed(lambda: solve_named_e4_from_e3(max_iter=20, e3=e3))
    rows.append(summarize("F_E4", e4, freeze_e4_q(), "composite_P", dt))

    e5, dt = timed(lambda: solve_named_e5_from_e4(max_iter=20, e4=e4))
    rows.append(summarize("F_E5", e5, freeze_e5_q(Q_W=0.0), "composite_P", dt))

    e6, dt = timed(lambda: solve_named_e6_from_e5(max_iter=20, e5=e5))
    rows.append(summarize("F_E6", e6, freeze_e6_q(Q_W=0.0), "composite_P", dt))

    r10, dt = timed(
        lambda: solve_named_e6_then_r10(
            Q_W=0.0, max_iter=20, apply_startwert=False, e6=e6
        )
    )
    rows.append(summarize("F_E6+R10", r10, freeze_e6_r10_q(Q_W=0.0), "composite_P", dt))

    payload = {
        "date": "2026-09-15",
        "note": (
            "Named-core chain after continuation review §7. "
            "R11 not connected. Not a furnace / HTW claim."
        ),
        "total_wall_s": time.perf_counter() - t_all,
        "r11_tar_conversion": False,
        "production_hamel_outer_refresh": False,
        "stages": rows,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {OUT}  total {payload['total_wall_s']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
