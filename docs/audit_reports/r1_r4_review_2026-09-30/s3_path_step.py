"""Wesseling 1 输入路径的第一小步。中间点不叫 Wesseling，不过门不往前走。"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
from src.solvers.wesseling_two_zone import q_on_wesseling1_path, wesseling1_input_path

_ROOT = Path(__file__).resolve().parents[3]
_OUT = Path(__file__).resolve().parent / "s3_path_step.json"
_BASE = np.load(_ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
_S = 0.02


def _git() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=_ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _inputs(nq) -> dict:
    feed = nq.energy.feed
    ref = nq.reference
    return {
        "s": _S,
        "name": "not CASE_HTW_WESSELING_1",
        "O2_mol_s": float(nq.energy.species.N_in_b[5]),
        "N2_mol_s": float(nq.energy.species.N_in_b[6]),
        "fuel_C_mol_s": float(feed.nC_fuel),
        "fuel_H_mol_s": float(feed.nH_fuel),
        "fuel_O_mol_s": float(feed.nO_fuel),
        "moisture_kg_s": float(feed.m_moist_kg_s),
        "h_f_dry_J_kg": float(feed.h_f_dry_J_kg),
        "area_m2": float(ref.area),
        "H_bed_m": float(ref.H_bed),
        "H_freeboard_m": float(ref.H_freeboard),
        "P_Pa": float(ref.P),
        "T_in_K": float(nq.energy.T_in),
        "Q_W": float(nq.Q_W),
        "d_p_m": [float(v) for v in nq.energy.species.d_p],
        "trajectory_continue_s": float(ref.trajectory_continue_s),
        "fuel_size_applied_to_bins": False,
        "fuel_heating_value": None,
    }


def main() -> None:
    nq = q_on_wesseling1_path(_S)
    x0 = np.asarray(_BASE, dtype=np.float64).copy()
    out = solve_named_controlled(nq, x0, deadline_s=320.0)
    x = np.asarray(out["x"], dtype=np.float64)
    q = outer_P_named(x, nq)
    ledger = external_imbalance(x, q)
    row = {
        "name": "not CASE_HTW_WESSELING_1",
        "s": _S,
        "target_reached": False,
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "gate_ratio": out["gate_ratio"],
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
        "inputs": _inputs(nq),
        "input_path": wesseling1_input_path(),
        "git": _git(),
        "x": x.tolist(),
        "F": np.asarray(out["F"], dtype=np.float64).tolist(),
        "ledger": ledger,
    }
    def _plain(obj):
        if isinstance(obj, dict):
            return {str(k): _plain(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_plain(v) for v in obj]
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        return obj

    _OUT.write_text(json.dumps(_plain(row)), encoding="utf-8")
    print(row["passed"], row["reason"], round(row["elapsed_s"], 1), row["gate_ratio"], flush=True)


if __name__ == "__main__":
    main()
