"""U3 新失败点上的一次有界方向。不跑 320 s。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    material_gate_ratios,
    step_reduces_material_violation,
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, residual_block_merit
from src.solvers.wirsum_1997_named_f_p2 import bed_face_exchange
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_X = Path(__file__).resolve().parents[1] / "t1_t4_review_2026-09-30" / "u3_path.json"
_OUT = Path(__file__).resolve().parent / "v1_fixed.json"


def main() -> None:
    saved = json.loads(_X.read_text())["s_0_16_controlled"]
    x = np.asarray(saved["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.16)
    q = outer_P_named(x, nq)
    face_k, face_b, *_ = bed_face_exchange(x, q.reference, T_bed_K=float(x[-2]))
    stats: dict = {}
    proposed = _e8_log_scaled_correct(
        q, x, max_nfev=8, max_j=4, ftol=1.0e-4,
        temperature_factor=1.5, log_span=0.5, stats=stats,
    )
    f0 = coupled_residual_named(x, q)
    positive = (x > 0.0) & (proposed > 0.0)
    rows = []
    for alpha in (0.0, 0.001, 0.01, 0.1, 1.0):
        y = x.copy()
        y[positive] = x[positive] * (proposed[positive] / x[positive]) ** alpha
        y[~positive] = x[~positive] + alpha * (proposed[~positive] - x[~positive])
        qy = outer_P_named(y, nq)
        fy = coupled_residual_named(y, qy)
        rows.append({
            "alpha": alpha,
            "merit": residual_block_merit(fy),
            "material_ratios": material_gate_ratios(fy),
            "energy_W": fy[-2:].tolist(),
            "material_improved": alpha == 0.0 or step_reduces_material_violation(f0, fy),
            "transport_closed": bool(np.all(np.asarray(qy.hydro.w_undecided) <= 1.0e-12)),
        })
        print(alpha, rows[-1]["material_ratios"], rows[-1]["energy_W"], rows[-1]["merit"], flush=True)
    small = [r for r in rows if 0.0 < r["alpha"] <= 0.1 and r["material_improved"]]
    out = {
        "state": "U3 s=0.16 controlled, not a root",
        "K_eject_matches_bed_face": bool(np.allclose(face_k, q.energy.species.K_eject)),
        "K_bd_matches_bed_face": bool(np.isclose(face_b, q.energy.species.geo_ex.K_bd)),
        "inner_status": stats.get("inner_status"),
        "F_calls": stats.get("F_calls"),
        "J_calls": stats.get("J_calls"),
        "optimizer_entrance": stats.get("inner_merit_entrance"),
        "optimizer_best": stats.get("inner_merit_best"),
        "rows": rows,
        "small_step_material_improved": bool(small),
    }
    _OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("improved", out["small_step_material_improved"], "K", out["K_eject_matches_bed_face"], out["K_bd_matches_bed_face"], flush=True)


if __name__ == "__main__":
    main()
