"""空气 -5% 与燃料 +5%。从 P3 端点延拓到目标输入，通过只认分块门。

同一参数的三例和 E7 链没有新的求解改动，本脚本不重跑。
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f_control import (
    advance_parameter_factor,
    solve_named_controlled,
)
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_pyrolysis_ledger import dry_mass_from_cho_atoms

_ROOT = Path(__file__).resolve().parents[3]
_ENDPOINT = _ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz"
_OUT = _ROOT / "docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp4_parameter_branches.json"
_ATTEMPT_S = 120.0
_WALL_S = 300.0


def _base():
    return replace(freeze_e8_q(Q_W=0.0), stage=11)


def _air(factor: float):
    nq = _base()
    sp = nq.energy.species
    return replace(
        nq,
        energy=replace(
            nq.energy,
            species=replace(sp, N_in_b=np.asarray(sp.N_in_b, dtype=np.float64) * factor),
        ),
    )


def _fuel(factor: float):
    nq = _base()
    feed = nq.energy.feed
    scaled = replace(
        feed,
        nC_fuel=float(feed.nC_fuel) * factor,
        nH_fuel=float(feed.nH_fuel) * factor,
        nO_fuel=float(feed.nO_fuel) * factor,
        nN_fuel=float(feed.nN_fuel) * factor,
        nS_fuel=float(feed.nS_fuel) * factor,
        nC_fixed=float(feed.nC_fixed) * factor,
        nC_vm=float(feed.nC_vm) * factor,
        m_moist_kg_s=float(feed.m_moist_kg_s) * factor,
    )
    scaled = replace(scaled, m_dry_kg_s=dry_mass_from_cho_atoms(scaled))
    return replace(nq, energy=replace(nq.energy, feed=scaled))


def _record(factor: float, out: dict) -> dict:
    x = np.asarray(out["x"], dtype=np.float64)
    blocks = out["blocks"] or {}
    return {
        "factor": factor,
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "merit": out["merit"],
        "gate_ratio": out["gate_ratio"],
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
        "gas_l2": None if "gas_l2" not in blocks else float(blocks["gas_l2"]),
        "worst_index": int(np.argmax(np.abs(out["F"]))),
        "worst_value": float(out["F"][int(np.argmax(np.abs(out["F"])))]),
    }


def continue_branch(name: str, make, target: float, x0: np.ndarray) -> dict:
    current = 1.0
    step = abs(target - 1.0)
    x = np.asarray(x0, dtype=np.float64).copy()
    path = []
    t0 = time.monotonic()
    stop = False
    while abs(current - target) > 1.0e-9 and not stop:
        remain = _WALL_S - (time.monotonic() - t0)
        if remain < 20.0:
            break
        direction = 1.0 if target > current else -1.0
        trial = current + direction * min(step, abs(target - current))
        out = solve_named_controlled(make(trial), x, deadline_s=min(_ATTEMPT_S, remain))
        row = _record(trial, out)
        path.append(row)
        print(name, round(trial, 4), row["passed"], row["reason"], round(row["elapsed_s"], 1), flush=True)
        if out["passed"]:
            x = np.asarray(out["x"], dtype=np.float64)
        current, step, stop = advance_parameter_factor(
            current, target, step, passed=bool(out["passed"])
        )
    return {
        "name": name,
        "class": "parameter_branch",
        "target_factor": target,
        "certified_factor": current,
        "reached_target": abs(current - target) <= 1.0e-9,
        "elapsed_s": time.monotonic() - t0,
        "path": path,
        "T_bed": None if not path else path[-1]["T_bed"],
        "T_fb": None if not path else path[-1]["T_fb"],
    }


def finish_remaining(name: str, make, mid: float, target: float) -> dict:
    """从端点重解已认证的半步，再用完整时限走向目标因子。"""
    x = np.load(_ENDPOINT)["x"]
    path = []
    out_mid = solve_named_controlled(make(mid), x, deadline_s=150.0)
    path.append(_record(mid, out_mid))
    print(name, mid, path[-1]["passed"], path[-1]["reason"], round(path[-1]["elapsed_s"], 1), flush=True)
    current = mid if out_mid["passed"] else 1.0
    if out_mid["passed"]:
        out = solve_named_controlled(make(target), out_mid["x"], deadline_s=160.0)
        path.append(_record(target, out))
        print(name, target, path[-1]["passed"], path[-1]["reason"], round(path[-1]["elapsed_s"], 1), flush=True)
        if out["passed"]:
            current = target
    return {
        "name": name,
        "class": "parameter_branch",
        "target_factor": target,
        "certified_factor": current,
        "reached_target": abs(current - target) <= 1.0e-9,
        "path": path,
        "T_bed": path[-1]["T_bed"],
        "T_fb": path[-1]["T_fb"],
    }


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("air-rest", "fuel-rest"):
        if which == "air-rest":
            row = finish_remaining("air_minus_5pct", _air, 0.975, 0.95)
            out = _OUT.with_name("wp4_parameter_branches_air_rest.json")
        else:
            row = finish_remaining("fuel_plus_5pct", _fuel, 1.025, 1.05)
            out = _OUT.with_name("wp4_parameter_branches_fuel_rest.json")
        out.write_text(json.dumps({"parameter_branch": [row]}, indent=2), encoding="utf-8")
        print("wrote", out, row["reached_target"], flush=True)
        return
    out = _OUT if which == "both" else _OUT.with_name(f"wp4_parameter_branches_{which}.json")
    x0 = np.load(_ENDPOINT)["x"]
    jobs = []
    if which in ("both", "air"):
        jobs.append(("air_minus_5pct", _air, 0.95))
    if which in ("both", "fuel"):
        jobs.append(("fuel_plus_5pct", _fuel, 1.05))
    rows = []
    for name, make, target in jobs:
        rows.append(continue_branch(name, make, target, x0))
        payload = {
            "same_parameter_initial_value": "cited wp1_three_cases.json; not rerun",
            "e7_chain": {
                "status": "not_rerun",
                "reason": "solve_named_e8_from_e7 unchanged; prior timeout remains hydro 0.925 fate 0",
            },
            "cold_start": "cited wp1 cold_start; not rerun",
            "parameter_branch": rows,
        }
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
