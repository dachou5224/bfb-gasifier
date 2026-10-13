"""工作包 D。现场重认证，显式步长，输运未闭合时不给出 Wesseling 出口预测。"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_control import (
    solve_named_controlled,
    transport_is_closed,
)
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_pyrolysis_ledger import dry_mass_from_cho_atoms
from src.solvers.wesseling_case_contract import map_wesseling_inputs
from src.solvers.wesseling_two_zone import (
    build_wesseling1_two_zone_q,
    wesseling1_input_path,
    wesseling1_physical_x0,
)

_ROOT = Path(__file__).resolve().parents[3]
_DIR = _ROOT / "docs/audit_reports/wp1_wp4_review_2026-09-28"
_ENDPOINT = _ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz"


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
    return replace(
        nq,
        energy=replace(nq.energy, feed=replace(scaled, m_dry_kg_s=dry_mass_from_cho_atoms(scaled))),
    )


def _blocks(x, nq) -> dict:
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q, rate_scale=1.0)
    blocks = named_endpoint_blocks(f)
    return {
        "passed": bool(blocks["passed"]),
        "transport_closed": transport_is_closed(q.hydro),
        "gas_l2": float(blocks["gas_l2"]),
        "solid_l2": float(blocks["solid_l2"]),
        "tar_l2": float(blocks["tar_l2"]),
        "energy_W": float(blocks["energy_W"]),
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
    }


def recert_p3() -> dict:
    x = np.load(_ENDPOINT)["x"]
    row = _blocks(x, _base())
    row["name"] = "P3"
    row["x_source"] = str(_ENDPOINT)
    return row


def _step(name: str, make, factor: float, x: np.ndarray) -> tuple[dict, np.ndarray]:
    out = solve_named_controlled(make(factor), x, deadline_s=120.0)
    got = np.asarray(out["x"], dtype=np.float64)
    row = {
        "name": name,
        "factor": factor,
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "transport_closed": bool(out["transport_closed"]),
        "gate_ratio": out["gate_ratio"],
        "T_bed": float(got[-2]),
        "T_fb": float(got[-1]),
        "x": got.tolist(),
    }
    print(name, factor, row["passed"], row["reason"], round(row["elapsed_s"], 1), flush=True)
    return row, got if out["passed"] else x


def branch(name: str, make, factors: list[float]) -> dict:
    x = np.load(_ENDPOINT)["x"]
    path = []
    for factor in factors:
        row, x = _step(name, make, factor, x)
        path.append(row)
        if not row["passed"]:
            break
    target = factors[-1]
    reached = bool(path) and path[-1]["factor"] == target and path[-1]["passed"]
    return {
        "name": name,
        "factors_requested": factors,
        "reached_target": reached,
        "stop_is_not_a_0p01_threshold": True,
        "path": path,
    }


def wesseling() -> dict:
    nq = build_wesseling1_two_zone_q()
    x = wesseling1_physical_x0(nq)
    q = outer_P_named(x, nq)
    closed = transport_is_closed(q.hydro)
    book = external_imbalance(x, q)
    blocks = named_endpoint_blocks(coupled_residual_named(x, q, rate_scale=1.0))
    two = map_wesseling_inputs("CASE_HTW_WESSELING_2")
    missing = [row["name"] for row in two["fields"] if row["status"] == "missing"]
    report = {
        "transport_closed": closed,
        "w_undecided": np.asarray(q.hydro.w_undecided, dtype=float).tolist(),
        "prediction": None,
        "reason": "transport_not_closed" if not closed else "not_solved",
        "trajectory_continue_s": float(nq.reference.trajectory_continue_s),
        "blocks": {
            "passed": bool(blocks["passed"]),
            "gas_l2": float(blocks["gas_l2"]),
            "solid_l2": float(blocks["solid_l2"]),
            "tar_l2": float(blocks["tar_l2"]),
            "energy_W": float(blocks["energy_W"]),
        },
        "gap_mol_s": book["gap_mol_s"],
        "enthalpy_gap_W": book["enthalpy_W"]["gap"],
        "input_path": wesseling1_input_path(),
        "full_element_certification": False,
        "wesseling2_missing": missing,
        "wesseling2_solved": False,
    }
    if closed:
        out = solve_named_controlled(nq, x, deadline_s=120.0)
        report["reason"] = out["reason"]
        report["passed"] = bool(out["passed"])
        if out["passed"]:
            y = np.asarray(out["x"], dtype=np.float64)
            report["prediction"] = {"T_bed": float(y[-2]), "T_fb": float(y[-1])}
    return report


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "p3"):
        row = recert_p3()
        path = _DIR / "wpd_p3.json"
        path.write_text(json.dumps(row, indent=2), encoding="utf-8")
        print("p3", row["passed"], flush=True)
    if which in ("all", "fuel"):
        row = branch("fuel_plus_5pct", _fuel, [1.025, 1.05])
        (_DIR / "wpd_fuel.json").write_text(json.dumps(row, indent=2), encoding="utf-8")
    if which in ("all", "air"):
        row = branch("air_minus_5pct", _air, [0.975, 0.95])
        (_DIR / "wpd_air.json").write_text(json.dumps(row, indent=2), encoding="utf-8")
    if which in ("all", "wesseling"):
        row = wesseling()
        (_DIR / "wpd_wesseling.json").write_text(json.dumps(row, indent=2), encoding="utf-8")
        print("wesseling", row["reason"], row["transport_closed"], flush=True)


if __name__ == "__main__":
    main()
