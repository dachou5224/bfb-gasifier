"""修正输运后的一个中间点和旧 0.22 附近。旧 x 只作初猜。不把 git HEAD 当成源码身份。"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    solve_named_controlled,
)
from src.solvers.wirsum_1997_named_f_p2 import finite_hydro_from_x
from src.solvers.wesseling_two_zone import (
    build_wesseling1_two_zone_q,
    q_on_wesseling1_path,
    wesseling1_physical_x0,
)

_ROOT = Path(__file__).resolve().parents[3]
_OUT = Path(__file__).resolve().parent / "t3_corrected_path.json"
_P3 = _ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz"
_OLD = _ROOT / "docs/audit_reports/r1_r4_review_2026-09-30/s3_path_step.json"
_SOURCES = (
    "src/physics/wirsum_1997_freeboard.py",
    "src/solvers/wirsum_1997_named_f.py",
    "src/solvers/wirsum_1997_named_f_e8.py",
    "src/solvers/wirsum_1997_named_f_control.py",
    "src/solvers/wirsum_1997_named_f_p2.py",
    "src/solvers/wirsum_1997_coupled_f_two_bin_pyro.py",
    "src/solvers/external_state_ledger.py",
    "src/solvers/wesseling_two_zone.py",
)


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


def _source_hash() -> dict:
    digest = hashlib.sha256()
    for rel in _SOURCES:
        digest.update(rel.encode())
        digest.update(Path(_ROOT / rel).read_bytes())
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=_ROOT, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--", *_SOURCES], cwd=_ROOT, text=True)
    return {
        "git_head": head,
        "git_head_is_not_the_source_identity": True,
        "source_sha256": digest.hexdigest(),
        "sources": list(_SOURCES),
        "dirty": dirty.strip().splitlines(),
    }


def _inputs(nq, s: float) -> dict:
    feed = nq.energy.feed
    ref = nq.reference
    return {
        "s": s,
        "name": "not CASE_HTW_WESSELING_1",
        "O2_mol_s": float(nq.energy.species.N_in_b[5]),
        "fuel_C_mol_s": float(feed.nC_fuel),
        "moisture_kg_s": float(feed.m_moist_kg_s),
        "h_f_dry_J_kg": float(feed.h_f_dry_J_kg),
        "area_m2": float(ref.area),
        "H_bed_m": float(ref.H_bed),
        "H_freeboard_m": float(ref.H_freeboard),
        "trajectory_continue_s": float(ref.trajectory_continue_s),
        "trajectory_extra_physical_s": float(ref.trajectory_extra_physical_s),
        "trajectory_max_work_steps": int(ref.trajectory_max_work_steps),
        "trajectory_wall_s": float(ref.trajectory_wall_s),
        "fuel_heating_value": None,
        "fuel_size_applied_to_bins": False,
    }


def _solve(s: float, x0: np.ndarray, deadline_s: float) -> dict:
    nq = q_on_wesseling1_path(s)
    out = solve_named_controlled(nq, x0, deadline_s=deadline_s)
    z = np.asarray(out["x"], dtype=np.float64)
    q = outer_P_named(z, nq)
    f = coupled_residual_named(z, q)
    cert = certify_controlled_state(z, f, q.hydro)
    led = external_imbalance(z, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    cho_ok = max(abs(v) for v in cho.values()) < 1.0e-4
    row = {
        "name": "not CASE_HTW_WESSELING_1",
        "s": s,
        "target_reached": False,
        "passed_block_gates": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "cho_gap_mol_s": cho,
        "cho_gap_compatible": cho_ok,
        "joint_pass": bool(cert["passed"] and cert["transport_closed"] and cho_ok),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "gate_ratio": cert["gate_ratio"],
        "T_bed": float(z[-2]),
        "T_fb": float(z[-1]),
        "inputs": _inputs(nq, s),
        "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
        "x": z.tolist(),
        "F": np.asarray(f, dtype=np.float64).tolist(),
    }
    print(s, row["joint_pass"], row["passed_block_gates"], row["transport_closed"], cho, round(out["elapsed_s"], 1), flush=True)
    return row


def main() -> None:
    nq = build_wesseling1_two_zone_q()
    x_phys = wesseling1_physical_x0(nq)
    snap, _geo, _k = finite_hydro_from_x(
        x_phys, nq.reference, T_bed_K=float(x_phys[-2]), T_fb_K=float(x_phys[-1]), N_fb=x_phys[16:23],
    )
    fixed = {
        "w_undecided": np.asarray(snap.w_undecided, dtype=float).tolist(),
        "w_leave": np.asarray(snap.w_leave, dtype=float).tolist(),
        "extra_physical_s": float(nq.reference.trajectory_extra_physical_s),
        "max_work_steps": int(nq.reference.trajectory_max_work_steps),
        "wall_s": float(nq.reference.trajectory_wall_s),
    }
    print("fixed", fixed["w_undecided"], flush=True)
    p3 = np.load(_P3)["x"]
    old = json.loads(_OLD.read_text())
    old_rows = [row for row in old["continued"] if row.get("passed") and abs(float(row["s"]) - 0.22) < 1.0e-12]
    old_x = np.asarray(old_rows[-1]["x"], dtype=np.float64)
    report = {
        "source": _source_hash(),
        "fixed_field": fixed,
        "s_0_02_from_p3": _solve(0.02, p3, 200.0),
        "s_0_22_from_old_guess": _solve(0.22, old_x, 320.0),
    }
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")


if __name__ == "__main__":
    main()
