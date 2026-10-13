"""AA1：s=0.40 预测一次，只在 s=0.402 校正。通过后才跑 AA2 的短批次。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    progress_phi,
    solve_named_controlled,
    state_in_temporary_domain,
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_ROOT = Path(__file__).resolve().parents[3]
_OUT = Path(__file__).resolve().parent / "aa1_aa2.json"
_SAVED = _ROOT / "docs/audit_reports/x1_x3_review_2026-10-07/y2_continue.json"
_CORRECT_S = 60.0
_BATCH_S = 300.0
_MAX_STEPS = 5


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


def _full(s: float, x: np.ndarray) -> np.ndarray:
    return coupled_residual_named(x, outer_P_named(x, q_on_wesseling1_path(s)))


def _ledger(x: np.ndarray, q) -> dict:
    led = external_imbalance(x, q)
    cho = {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")}
    return {"cho_gap_mol_s": cho, "enthalpy_gap_W": float(led["enthalpy_W"]["gap"])}


def _cold(s: float, x: np.ndarray) -> dict:
    q = outer_P_named(x, q_on_wesseling1_path(s))
    f = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, f, q.hydro)
    led = _ledger(x, q)
    cho_ok = max(abs(v) for v in led["cho_gap_mol_s"].values()) < 1.0e-4
    return {
        "passed_block_gates": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "finite": bool(cert["finite"]),
        "in_domain": bool(cert["in_domain"]),
        "gate_ratio": cert["gate_ratio"],
        "energy_rows_W": cert["energy_rows_W"],
        "phi": progress_phi(f),
        "joint_pass": bool(cert["passed"] and cert["transport_closed"] and cho_ok),
        **led,
        "x": np.asarray(x, dtype=np.float64),
        "F": np.asarray(f, dtype=np.float64),
    }


def _predict(s: float, x: np.ndarray, ds: float) -> dict:
    h = 1.0e-4
    dF = (_full(s + h, x) - _full(s - h, x)) / (2.0 * h)
    q = outer_P_named(x, q_on_wesseling1_path(s))
    stats: dict = {}
    t0 = time.monotonic()
    proposed = _e8_log_scaled_correct(
        q, x, log_span=0.5, ratio_span=1.0, stats=stats, residual_offset=dF * ds,
    )
    z = np.asarray(proposed, dtype=np.float64)
    q_t = outer_P_named(z, q_on_wesseling1_path(s + ds))
    f_t = coupled_residual_named(z, q_t)
    finite_f = bool(np.all(np.isfinite(f_t)))
    return {
        "elapsed_s": time.monotonic() - t0,
        "x": z,
        "F": f_t,
        "phi": progress_phi(f_t) if finite_f else None,
        "gate_ratio": certify_controlled_state(z, f_t, q_t.hydro)["gate_ratio"] if finite_f else None,
        "finite": bool(np.all(np.isfinite(z)) and finite_f),
        "in_domain": bool(state_in_temporary_domain(z)),
        "subproblem_solved": bool(stats.get("subproblem_solved")),
        "linear_status": stats.get("linear_status"),
        "solver_status": stats.get("solver_status"),
        "solver_message": stats.get("solver_message"),
        "calls": {
            "F": stats.get("F_calls"),
            "J": stats.get("J_calls"),
            "P": 2.0,
        },
    }


def _calls(stats: dict) -> dict:
    return {"F": stats.get("F_calls"), "J": stats.get("J_calls"), "P": stats.get("P_calls")}


def _diagnose(reason: str, history: list, cold: dict) -> str:
    if reason == "deadline":
        return "budget_exhausted"
    if reason == "subproblem_not_solved":
        return "subproblem_not_solved"
    if reason in ("model_mismatch", "constraint_stall"):
        return "model_mismatch" if reason == "model_mismatch" else "corrector_stall"
    if reason == "no_descent_direction":
        return "corrector_stall"
    if cold["passed_block_gates"] and not cold["joint_pass"]:
        return "external_ledger"
    ratios = cold["gate_ratio"]
    names = ("gas", "solid", "tar")
    dominant = max(names, key=lambda name: float(ratios[name]))
    energy = max(abs(float(v)) for v in ratios["energy_rows"])
    if energy >= float(ratios[dominant]):
        return "energy"
    return dominant


def _gas_mismatch(history: list) -> bool:
    hits = 0
    for row in history:
        if not row.get("accepted"):
            continue
        predicted = row.get("predicted_material_ratios")
        actual = row.get("material_ratios")
        if not predicted or not actual:
            continue
        if float(actual[0]) > 10.0 * max(float(predicted[0]), 1.0e-6):
            hits += 1
    return hits >= 2


def _jv_gap(s: float, x0: np.ndarray, x1: np.ndarray) -> dict:
    direction = np.asarray(x1, dtype=np.float64) - np.asarray(x0, dtype=np.float64)
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        return {"skipped": True}
    unit = direction / norm
    eps = 1.0e-4
    q = outer_P_named(x0, q_on_wesseling1_path(s))

    def full(step: float) -> np.ndarray:
        return coupled_residual_named(x0 + step * unit, outer_P_named(x0 + step * unit, q_on_wesseling1_path(s)))

    def frozen(step: float) -> np.ndarray:
        from src.solvers.wirsum_1997_named_f import inner_hydro_frozen_residual
        return inner_hydro_frozen_residual(x0 + step * unit, q)

    d_full = (full(eps) - full(-eps)) / (2.0 * eps)
    d_frozen = (frozen(eps) - frozen(-eps)) / (2.0 * eps)
    scale = max(float(np.linalg.norm(d_full)), 1.0e-30)
    return {
        "skipped": False,
        "relative_l2": float(np.linalg.norm(d_full - d_frozen) / scale),
        "full_l2": float(np.linalg.norm(d_full)),
    }


def _correct(s: float, x: np.ndarray, budget_s: float) -> dict:
    t0 = time.monotonic()
    out = solve_named_controlled(q_on_wesseling1_path(s), x, deadline_s=budget_s)
    cold = _cold(s, np.asarray(out["x"], dtype=np.float64))
    history = out["stats"].get("history", [])
    return {
        "reason": out["reason"],
        "elapsed_s": time.monotonic() - t0,
        "solver_elapsed_s": out["elapsed_s"],
        "calls": _calls(out["stats"]),
        "history": history,
        "temperature_bound_active": any(bool(row.get("temperature_bound_active")) for row in history),
        "diagnosis": _diagnose(out["reason"], history, cold),
        "cold": cold,
    }


def _next_ds(ds: float, corrected: dict) -> float:
    history = corrected["history"]
    span = float(history[-1]["log_span"]) if history else 0.5
    outers = len(history)
    gas = float(corrected["cold"]["gate_ratio"]["gas"])
    expensive = outers > 4 or span < 0.49 or gas > 0.5 or corrected["elapsed_s"] > 20.0
    active_temp = bool(corrected["temperature_bound_active"])
    if expensive or active_temp:
        return ds
    return min(0.02, ds * 2.0)


def main() -> None:
    saved = json.loads(_SAVED.read_text())["points"][-1]
    s = float(saved["s"])
    x = np.asarray(saved["x"], dtype=np.float64)
    base = _cold(s, x)
    if not base["joint_pass"]:
        raise SystemExit("saved s=0.40 is not a joint pass")
    ds = 0.002
    prediction = _predict(s, x, ds)
    entrance_ok = bool(prediction["finite"] and prediction["in_domain"])
    report: dict = {
        "source_sha256": named_kernel_source_sha256(),
        "s0": s,
        "target_s": s + ds,
        "prediction_updates_certified_checkpoint": False,
        "screen_1W_phi20_required": False,
        "base_joint_pass": True,
        "prediction": {k: v for k, v in prediction.items() if k not in ("x", "F")},
        "prediction_x": prediction["x"],
        "entrance_ok": entrance_ok,
    }
    if not entrance_ok or not prediction["subproblem_solved"]:
        report["aa1"] = {
            "joint_pass": False,
            "certified": False,
            "diagnosis": "predictor_inadmissible" if not entrance_ok else "subproblem_not_solved",
        }
        report["aa2"] = {"started": False}
        _OUT.write_text(json.dumps(_plain(report), indent=2), encoding="utf-8")
        print(json.dumps(_plain({k: v for k, v in report.items() if k != "prediction_x"}), indent=2))
        return
    corrected = _correct(s + ds, prediction["x"], _CORRECT_S)
    aa1 = {
        "joint_pass": corrected["cold"]["joint_pass"],
        "certified": corrected["cold"]["joint_pass"],
        "reason": corrected["reason"],
        "diagnosis": corrected["diagnosis"],
        "elapsed_s": corrected["elapsed_s"],
        "calls": corrected["calls"],
        "gate_ratio": corrected["cold"]["gate_ratio"],
        "energy_rows_W": corrected["cold"]["energy_rows_W"],
        "cho_gap_mol_s": corrected["cold"]["cho_gap_mol_s"],
        "enthalpy_gap_W": corrected["cold"]["enthalpy_gap_W"],
        "phi": corrected["cold"]["phi"],
        "temperature_bound_active": corrected["temperature_bound_active"],
        "history": corrected["history"],
        "gas_model_mismatch": _gas_mismatch(corrected["history"]),
    }
    if aa1["gas_model_mismatch"]:
        aa1["jv"] = _jv_gap(s + ds, prediction["x"], corrected["cold"]["x"])
    report["aa1"] = aa1
    if aa1["certified"]:
        report["certified_x"] = corrected["cold"]["x"]
        report["certified_F"] = corrected["cold"]["F"]
    if not aa1["joint_pass"]:
        report["trial_x"] = corrected["cold"]["x"]
        report["trial_F"] = corrected["cold"]["F"]
        report["aa2"] = {"started": False}
        _OUT.write_text(json.dumps(_plain(report), indent=2), encoding="utf-8")
        print(json.dumps(_plain({k: report[k] for k in report if k not in ("prediction_x", "trial_x", "trial_F", "certified_x", "certified_F")}), indent=2))
        return
    points = []
    cursor_s = s + ds
    cursor_x = np.asarray(corrected["cold"]["x"], dtype=np.float64)
    t_batch = time.monotonic()
    while len(points) < _MAX_STEPS and time.monotonic() - t_batch < _BATCH_S and cursor_s < 1.0 - 1.0e-12:
        remain = _BATCH_S - (time.monotonic() - t_batch)
        step_ds = ds if cursor_s + ds <= 1.0 else 1.0 - cursor_s
        pred = _predict(cursor_s, cursor_x, step_ds)
        if not (pred["finite"] and pred["in_domain"] and pred["subproblem_solved"]):
            points.append({
                "s": cursor_s + step_ds,
                "ds": step_ds,
                "accepted": False,
                "diagnosis": "predictor_inadmissible" if not (pred["finite"] and pred["in_domain"]) else "subproblem_not_solved",
                "prediction_phi": pred["phi"],
            })
            break
        budget = min(_CORRECT_S, remain)
        if budget <= 1.0:
            points.append({"s": cursor_s + step_ds, "accepted": False, "diagnosis": "budget_exhausted"})
            break
        step = _correct(cursor_s + step_ds, pred["x"], budget)
        row = {
            "s": cursor_s + step_ds,
            "ds": step_ds,
            "accepted": step["cold"]["joint_pass"],
            "reason": step["reason"],
            "diagnosis": step["diagnosis"],
            "elapsed_s": step["elapsed_s"],
            "calls": step["calls"],
            "gate_ratio": step["cold"]["gate_ratio"],
            "phi": step["cold"]["phi"],
            "temperature_bound_active": step["temperature_bound_active"],
            "prediction_phi": pred["phi"],
        }
        if step["cold"]["joint_pass"]:
            cursor_s = cursor_s + step_ds
            cursor_x = np.asarray(step["cold"]["x"], dtype=np.float64)
            row["x"] = cursor_x
            row["F"] = step["cold"]["F"]
            ds = _next_ds(step_ds, step)
            points.append(row)
        else:
            row["rolled_back_to"] = cursor_s
            points.append(row)
            break
    accepted = [row for row in points if row.get("accepted")]
    report["aa2"] = {
        "started": True,
        "n_accepted": len(accepted),
        "n_attempted": len(points),
        "wall_s": time.monotonic() - t_batch,
        "last_certified_s": cursor_s,
        "net_advance": cursor_s - (s + ds if False else s + 0.002),
        "points": points,
    }
    if accepted:
        report["certified_x"] = cursor_x
        report["certified_s"] = cursor_s
    _OUT.write_text(json.dumps(_plain(report), indent=2), encoding="utf-8")
    printable = {k: report[k] for k in report if not str(k).endswith("_x") and not str(k).endswith("_F")}
    print(json.dumps(_plain(printable), indent=2), flush=True)


if __name__ == "__main__":
    main()
