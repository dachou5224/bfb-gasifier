"""从保存的 s=0.700 trial 做一次联合校正。不写认证文件。"""
import json
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.solvers.wirsum_1997_named_f import char_face_flows_kg_s, coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state, solve_named_controlled
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_SOLID = (14, 15, 37, 38, 39, 46, 47)
_HERE = Path(__file__).resolve().parent
_TRIAL = _HERE.parents[0] / "equation_loop_2026-10-08" / "advance_0.700.json"


def _picture(x, q, f):
    _eject, _back, leave = char_face_flows_kg_s(x, q)
    k = q.energy.species.k
    return {
        "solid_rows_kg_s": {str(i): float(f[i]) for i in _SOLID},
        "cyclone_mass_kg": [float(x[37]), float(x[38])],
        "cyclone_inflow_kg_s": float(np.sum(leave)),
        "recycle_kg_s": float(k.K_ab) * float(x[37]),
        "energy_rows_W": [float(f[-2]), float(f[-1])],
        "w_leave": [float(v) for v in q.hydro.w_leave],
        "tau_s": [float(v) for v in q.hydro.tau_s],
    }


def main() -> None:
    saved = json.loads(_TRIAL.read_text())
    x0 = np.asarray(saved["trial_x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.7)
    q0 = outer_P_named(x0, nq)
    f0 = coupled_residual_named(x0, q0)
    out = solve_named_controlled(
        nq, x0, deadline_s=60.0, temperature_radius_scale=4.0,
        settling_margin_eta=None, energy_guard=False,
    )
    x = np.asarray(out["x"], dtype=np.float64)
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, f, q.hydro)
    led = external_imbalance(x, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    joint = bool(cert["passed"] and max(abs(v) for v in cho.values()) < 1.0e-4)
    row = {
        "source_sha256": named_kernel_source_sha256(),
        "temperature_radius_scale": 4.0,
        "settling_margin_eta": None,
        "energy_guard": False,
        "deadline_s": 60.0,
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "passed": bool(out["passed"]),
        "joint_pass": joint,
        "gate_ratio": cert["gate_ratio"],
        "cho_gap_mol_s": cho,
        "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
        "calls": {
            "F": out["stats"].get("F_calls"),
            "J": out["stats"].get("J_calls"),
            "P": out["stats"].get("P_calls"),
        },
        "entrance": _picture(x0, q0, f0),
        "exit": _picture(x, q, f),
        "history": out["stats"].get("history"),
        "x": x,
        "F": f,
    }
    text = json.dumps(row, indent=2, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else v)
    (_HERE / "w2_result.json").write_text(text)
    brief = {k: row[k] for k in row if k not in ("x", "F", "history")}
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
