"""Z1/Z2 固定点验收。不缩小参数步长，不调用受控延拓。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import (
    named_kernel_source_files,
    named_kernel_source_sha256,
)
from src.solvers.wirsum_1997_named_f import (
    apply_tar_reaction_transport_rows,
    coupled_residual_named,
    inner_hydro_frozen_residual,
    linear_tar_row_map,
    outer_P_named,
)
from src.solvers.wirsum_1997_named_f_control import (
    certify_controlled_state,
    progress_phi,
)
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_ROOT = Path(__file__).resolve().parents[3]
_OUT = Path(__file__).resolve().parent / "z1_z2.json"


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


def _frozen(s: float, x: np.ndarray) -> np.ndarray:
    return inner_hydro_frozen_residual(x, q_on_wesseling1_path(s))


def _quality(stats: dict) -> dict:
    keys = (
        "linear_status",
        "subproblem_solved",
        "solver_success",
        "solver_status",
        "solver_message",
        "solver_nit",
        "max_nd_violation",
        "projected_gradient_norm",
        "inexact_descent",
        "active_bound_count",
        "v_abs_max",
        "objective",
        "fallback_success",
        "fallback_status",
        "fallback_message",
        "fallback_feasible",
    )
    return {key: stats.get(key) for key in keys}


def main() -> None:
    saved = json.loads(
        (_ROOT / "docs/audit_reports/x1_x3_review_2026-10-07/y2_continue.json").read_text()
    )["points"][-1]
    s = float(saved["s"])
    x = np.asarray(saved["x"], dtype=np.float64)
    q = outer_P_named(x, q_on_wesseling1_path(s))
    f = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, f, q.hydro)
    led = external_imbalance(x, q)
    derivatives = {}
    for h in (1.0e-3, 1.0e-4, 1.0e-5):
        derivatives[h] = (_full(s + h, x) - _full(s - h, x)) / (2.0 * h)
    reference = derivatives[1.0e-4]
    ref_norm = float(np.linalg.norm(reference))
    derivative_gap = {
        str(h): float(np.linalg.norm(row - reference) / ref_norm)
        for h, row in derivatives.items()
    }
    h_ref = 1.0e-4
    frozen = (_frozen(s + h_ref, x) - _frozen(s - h_ref, x)) / (2.0 * h_ref)
    linearization = []
    for ds in (2.0e-3, 1.0e-3, 5.0e-4):
        actual = _full(s + ds, x)
        remainder = actual - f - reference * ds
        linearization.append({
            "ds": ds,
            "remainder_l2": float(np.linalg.norm(remainder)),
            "remainder_over_ds2": float(np.linalg.norm(remainder) / ds ** 2),
            "energy_remainder_W": remainder[-2:].tolist(),
        })
    offset = reference * 2.0e-3
    bad = apply_tar_reaction_transport_rows(x, q, offset)
    good = linear_tar_row_map(offset)
    stats: dict = {}
    proposed = _e8_log_scaled_correct(
        q, x, log_span=0.5, ratio_span=1.0, stats=stats, residual_offset=offset,
    )
    pred = np.asarray(stats["predicted_raw"], dtype=np.float64)
    z = np.asarray(proposed, dtype=np.float64)
    f_same = _full(s, z)
    f_shift = _full(s + 2.0e-3, z)
    failed = json.loads(
        (_ROOT / "docs/audit_reports/x1_x3_review_2026-10-07/y2_s048.json").read_text()
    )
    x48 = np.asarray(failed["x"], dtype=np.float64)
    q48 = outer_P_named(x48, q_on_wesseling1_path(0.48))
    stats48: dict = {}
    _e8_log_scaled_correct(q48, x48, log_span=0.5, ratio_span=1.0, stats=stats48)
    pred48 = np.asarray(stats48["predicted_raw"], dtype=np.float64)
    report = {
        "continuation": False,
        "new_root": False,
        "source_sha256": named_kernel_source_sha256(),
        "source_file_count": len(named_kernel_source_files()),
        "s040": {
            "passed_block_gates": bool(cert["passed"]),
            "transport_closed": bool(cert["transport_closed"]),
            "gate_ratio": cert["gate_ratio"],
            "cho_gap_mol_s": {k: float(led["gap_mol_s"][k]) for k in ("C", "H", "O")},
            "enthalpy_gap_W": float(led["enthalpy_W"]["gap"]),
        },
        "derivative_relative_to_h_1e-4": derivative_gap,
        "full_minus_frozen_hydro_relative": float(np.linalg.norm(reference - frozen) / ref_norm),
        "linearization_at_fixed_x": linearization,
        "offset_transform_error_l2": float(np.linalg.norm(bad - good)),
        "ds_0_002": {
            "subproblem": _quality(stats),
            "predicted_phi": progress_phi(pred),
            "predicted_energy_W": pred[-2:].tolist(),
            "passes_1W_and_phi20_screen": bool(
                max(abs(float(v)) for v in pred[-2:]) <= 1.0 and progress_phi(pred) <= 20.0
            ),
            "screen_is_not_a_branch_test": True,
            "max_state_change": float(np.max(np.abs(z - x))),
            "same_s_phi": progress_phi(f_same),
            "same_s_energy_W": f_same[-2:].tolist(),
            "shifted_s_phi": progress_phi(f_shift),
            "shifted_s_energy_W": f_shift[-2:].tolist(),
            "probe_is_not_a_root": True,
        },
        "s048_local_subproblem": {
            "subproblem": _quality(stats48),
            "predicted_phi": progress_phi(pred48),
            "predicted_energy_W": pred48[-2:].tolist(),
            "entrance_phi": progress_phi(coupled_residual_named(x48, q48)),
        },
    }
    _OUT.write_text(json.dumps(_plain(report), indent=2), encoding="utf-8")
    print(json.dumps(_plain(report), indent=2), flush=True)


if __name__ == "__main__":
    main()
