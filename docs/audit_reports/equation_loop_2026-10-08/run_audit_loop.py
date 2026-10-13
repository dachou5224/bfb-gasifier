"""方程审计闭环：baseline / audit / replay / advance。不覆盖 aa1_aa2.json。"""

from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import edges_from_x
from src.solvers.wirsum_1997_named_f import (
    coupled_residual_named,
    outer_P_named,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    solve_named_controlled,
    state_in_temporary_domain,
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.vorabrechnung.vorab_gibbs_seed import pyrolysis_gibbs_account
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_DIR = Path(__file__).resolve().parent
_ROOT = _DIR.parents[2]
_PREV = _ROOT / "docs/audit_reports/z1_z3_review_2026-10-08/aa1_aa2.json"
_PREDICT = _ROOT / "docs/audit_reports/aa1_aa3_review_2026-10-08/predict_evidence.json"
_RUNS = _DIR / "runs.jsonl"
_MAX_SOLVES = 6


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


def _atomic(path: Path, obj) -> None:
    text = json.dumps(_plain(obj), indent=2)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)
    json.loads(path.read_text(encoding="utf-8"))


def _append(row: dict) -> None:
    line = json.dumps(_plain(row), ensure_ascii=False)
    with _RUNS.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def _versions() -> dict:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256": named_kernel_source_sha256(),
    }


def _solve_count() -> int:
    if not _RUNS.is_file():
        return 0
    count = 0
    for line in _RUNS.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        if json.loads(line).get("counts_as_coupled_solve"):
            count += 1
    return count


def _cold(s: float, x: np.ndarray) -> dict:
    q = outer_P_named(x, q_on_wesseling1_path(s))
    f = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, f, q.hydro)
    led = external_imbalance(x, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    cho_ok = max(abs(v) for v in cho.values()) < 1.0e-4
    u = unpack_named(x, 11)
    ed = edges_from_x(x[:39], q.energy.species)
    n_d = float(np.sum(u["N_d"]))
    tar_d = max(float(u["N_tar1"]), 0.0) + max(float(u["N_tar2"]), 0.0)
    return {
        "s": float(s),
        "passed_block_gates": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "finite": bool(cert["finite"]),
        "in_domain": bool(cert["in_domain"]),
        "joint_pass": bool(cert["passed"] and cho_ok),
        "gate_ratio": cert["gate_ratio"],
        "energy_rows_W": cert["energy_rows_W"],
        "cho_gap_mol_s": cho,
        "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
        "T_bed_K": float(u["T"]),
        "T_fb_K": float(u["T_fb"]),
        "y_tar_dense": tar_d / (n_d + tar_d) if n_d + tar_d > 0.0 else 0.0,
        "m_leave_kg_s": [float(v) for v in ed["m_leave"]],
        "m_rez_kg_s": float(ed["m_rez"]),
        "m_return_kg_s": [float(v) for v in ed["m_return"]],
        "F": np.asarray(f, dtype=np.float64),
        "x": np.asarray(x, dtype=np.float64),
    }


def _saved_points() -> list[dict]:
    blob = json.loads(_PREV.read_text(encoding="utf-8"))
    return [row for row in blob["aa2"]["points"] if row.get("accepted")]


def baseline() -> None:
    rows = []
    for point in _saved_points():
        s = float(point["s"])
        x = np.asarray(point["x"], dtype=np.float64)
        cold = _cold(s, x)
        saved_f = np.asarray(point["F"], dtype=np.float64)
        delta = np.asarray(cold["F"], dtype=np.float64) - saved_f
        cold["max_abs_F_minus_saved"] = float(np.max(np.abs(delta)))
        cold["saved_F_matches"] = bool(np.allclose(cold["F"], saved_f, rtol=0.0, atol=0.0))
        rows.append({k: v for k, v in cold.items() if k not in ("x", "F")})
        _append({
            "run_id": f"baseline-{s:.3f}",
            "mode": "baseline",
            "counts_as_coupled_solve": False,
            "s": s,
            "joint_pass": cold["joint_pass"],
            "max_abs_F_minus_saved": cold["max_abs_F_minus_saved"],
            **_versions(),
        })
    report = {"parent": str(_PREV), **_versions(), "points": rows}
    _atomic(_DIR / "baseline.json", report)
    if rows and rows[-1]["joint_pass"]:
        last = _saved_points()[-1]
        cold = _cold(float(last["s"]), np.asarray(last["x"], dtype=np.float64))
        _atomic(_DIR / "certified.json", {
            "s": cold["s"],
            "x": cold["x"],
            "F": cold["F"],
            "joint_pass": cold["joint_pass"],
            "gate_ratio": cold["gate_ratio"],
            "cho_gap_mol_s": cold["cho_gap_mol_s"],
            "energy_rows_W": cold["energy_rows_W"],
            **_versions(),
        })
    print(json.dumps(_plain(report), indent=2))


def _correct(s: float, x0: np.ndarray, budget_s: float, scale: float, settling_margin_eta: float | None) -> dict:
    t0 = time.monotonic()
    out = solve_named_controlled(
        q_on_wesseling1_path(s), x0, deadline_s=budget_s,
        temperature_radius_scale=scale,
        settling_margin_eta=settling_margin_eta,
    )
    trial = np.asarray(out["x"], dtype=np.float64)
    cold = _cold(s, trial)
    return {
        "reason": out["reason"],
        "elapsed_s": time.monotonic() - t0,
        "solver_elapsed_s": out["elapsed_s"],
        "calls": {
            "F": out["stats"].get("F_calls"),
            "J": out["stats"].get("J_calls"),
            "P": out["stats"].get("P_calls"),
        },
        "history": out["stats"].get("history"),
        "cold": cold,
        "trial_x": trial,
    }


def replay(s: float, budget_s: float, scale: float) -> None:
    if _solve_count() >= _MAX_SOLVES:
        raise SystemExit("coupled solve budget is exhausted")
    blob = json.loads(_PREDICT.read_text(encoding="utf-8"))
    row = next(item for item in blob if abs(float(item["s"]) - s) < 1.0e-9)
    x0 = np.asarray(row["predictor_x"], dtype=np.float64)
    if not state_in_temporary_domain(x0):
        raise SystemExit("predictor_x is outside the temporary domain")
    result = _correct(s, x0, budget_s, scale, None)
    record = {
        "run_id": f"replay-{s:.3f}-scale-{scale:g}",
        "mode": "replay",
        "counts_as_coupled_solve": True,
        "parent": str(_PREDICT),
        "s": s,
        "budget_s": budget_s,
        "temperature_radius_scale": scale,
        "reason": result["reason"],
        "elapsed_s": result["elapsed_s"],
        "calls": result["calls"],
        "joint_pass": result["cold"]["joint_pass"],
        "gate_ratio": result["cold"]["gate_ratio"],
        "energy_rows_W": result["cold"]["energy_rows_W"],
        "cho_gap_mol_s": result["cold"]["cho_gap_mol_s"],
        "history": result["history"],
        "trial_unavailable": False,
        **_versions(),
    }
    _append(record)
    _atomic(_DIR / f"replay_{s:.3f}_scale_{scale:g}.json", {
        **record,
        "x0": x0,
        "trial_x": result["trial_x"],
        "F": result["cold"]["F"],
    })
    print(json.dumps(_plain({k: record[k] for k in record if k != "history"}), indent=2))


def _predict(s: float, x: np.ndarray, ds: float, deadline_at: float, scale: float, settling_margin_eta: float | None) -> dict:
    h = 1.0e-4
    def residual(step: float) -> np.ndarray:
        return coupled_residual_named(x, outer_P_named(x, q_on_wesseling1_path(step)))
    dF = (residual(s + h) - residual(s - h)) / (2.0 * h)
    q = outer_P_named(x, q_on_wesseling1_path(s))
    remain = deadline_at - time.monotonic()
    stats: dict = {}
    if remain <= 0.0:
        return {"ok": False, "reason": "budget_exhausted"}
    proposed = _e8_log_scaled_correct(
        q, x, deadline_s=remain, log_span=0.5, ratio_span=1.0, stats=stats,
        residual_offset=dF * ds, temperature_radius_scale=scale,
        couple_trajectory_jacobian=True,
        settling_margin_eta=settling_margin_eta,
    )
    z = np.asarray(proposed, dtype=np.float64)
    return {
        "ok": bool(stats.get("subproblem_solved")) and state_in_temporary_domain(z),
        "subproblem_solved": bool(stats.get("subproblem_solved")),
        "x": z,
        "reason": None if stats.get("subproblem_solved") else "subproblem_not_solved",
    }


def advance(
    max_accepted: int,
    budget_s: float,
    scale: float,
    ds: float,
    *,
    fresh_budget: bool,
    settling_margin_eta: float | None = None,
) -> None:
    cert_path = _DIR / "certified.json"
    if not cert_path.is_file():
        raise SystemExit("no certified.json; run baseline first")
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    if not cert.get("joint_pass"):
        raise SystemExit("latest checkpoint is not a joint pass")
    s = float(cert["s"])
    x = np.asarray(cert["x"], dtype=np.float64)
    if float(ds) <= 0.0:
        raise SystemExit("ds must be positive")
    accepted = 0
    t0 = time.monotonic()
    while accepted < max_accepted and time.monotonic() - t0 < budget_s and s < 1.0 - 1.0e-12:
        if not fresh_budget and _solve_count() >= _MAX_SOLVES:
            _append({"run_id": "advance-stop", "mode": "advance", "reason": "solve_budget", "counts_as_coupled_solve": False})
            break
        deadline_at = t0 + budget_s
        remain = deadline_at - time.monotonic()
        step_budget = min(60.0, remain)
        gibbs_before = dict(pyrolysis_gibbs_account)
        if step_budget <= 1.0:
            break
        pred = _predict(s, x, ds, time.monotonic() + step_budget, scale, settling_margin_eta)
        if not pred["ok"]:
            _append({
                "run_id": f"advance-{s + ds:.3f}",
                "mode": "advance",
                "counts_as_coupled_solve": True,
                "s": s + ds,
                "accepted": False,
                "reason": pred.get("reason"),
                "rolled_back_to": s,
                "pyrolysis_gibbs": {
                    key: int(pyrolysis_gibbs_account[key] - gibbs_before[key])
                    for key in gibbs_before
                },
                **_versions(),
            })
            break
        result = _correct(
            s + ds, pred["x"], min(60.0, max(0.0, deadline_at - time.monotonic())), scale,
            settling_margin_eta,
        )
        row = {
            "run_id": f"advance-{s + ds:.3f}",
            "mode": "advance",
            "counts_as_coupled_solve": True,
            "s": s + ds,
            "ds": float(ds),
            "parent_s": s,
            "temperature_radius_scale": scale,
            "settling_margin_eta": settling_margin_eta,
            "reason": result["reason"],
            "elapsed_s": result["elapsed_s"],
            "calls": result["calls"],
            "pyrolysis_gibbs": {
                key: int(pyrolysis_gibbs_account[key] - gibbs_before[key])
                for key in gibbs_before
            },
            "joint_pass": result["cold"]["joint_pass"],
            "gate_ratio": result["cold"]["gate_ratio"],
            "history": result["history"],
            "trial_unavailable": False,
            **_versions(),
        }
        _append(row)
        payload = {**row, "x0": x, "prediction_x": pred["x"], "trial_x": result["trial_x"], "F": result["cold"]["F"]}
        _atomic(_DIR / f"advance_{s + ds:.3f}.json", payload)
        if result["cold"]["joint_pass"]:
            s = s + ds
            x = np.asarray(result["cold"]["x"], dtype=np.float64)
            accepted += 1
            _atomic(cert_path, {
                "s": s,
                "x": x,
                "F": result["cold"]["F"],
                "joint_pass": True,
                "gate_ratio": result["cold"]["gate_ratio"],
                "cho_gap_mol_s": result["cold"]["cho_gap_mol_s"],
                "energy_rows_W": result["cold"]["energy_rows_W"],
                "ds": float(ds),
                "temperature_radius_scale": scale,
                "settling_margin_eta": settling_margin_eta,
                **_versions(),
            })
        else:
            row["rolled_back_to"] = s
            break
    print(json.dumps({"accepted": accepted, "s": s, "wall_s": time.monotonic() - t0}))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("baseline", "audit", "replay", "advance"))
    parser.add_argument("--hypothesis", default="")
    parser.add_argument("--s", type=float, default=0.410)
    parser.add_argument("--budget-s", type=float, default=60.0)
    parser.add_argument("--max-accepted", type=int, default=2)
    parser.add_argument("--temperature-scale", type=float, default=1.0)
    parser.add_argument("--ds", type=float, default=None)
    parser.add_argument("--settling-margin-eta", type=float, default=None)
    parser.add_argument("--fresh-budget", action="store_true")
    args = parser.parse_args()
    if args.mode == "baseline":
        baseline()
    elif args.mode == "audit":
        print("structural audits live in tests/test_named_f_equation_invariants.py")
    elif args.mode == "replay":
        replay(args.s, min(60.0, args.budget_s), args.temperature_scale)
    else:
        if args.ds is None or args.ds <= 0.0:
            raise SystemExit("advance requires --ds")
        cap_n = args.max_accepted if args.fresh_budget else min(2, args.max_accepted)
        cap_s = args.budget_s if args.fresh_budget else min(120.0, args.budget_s)
        advance(
            cap_n, cap_s, args.temperature_scale, args.ds,
            fresh_budget=args.fresh_budget,
            settling_margin_eta=args.settling_margin_eta,
        )


if __name__ == "__main__":
    main()
