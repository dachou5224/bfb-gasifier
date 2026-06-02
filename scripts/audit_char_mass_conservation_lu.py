#!/usr/bin/env python3
"""Audit LU case char mass conservation before any kinetics tuning.

This script separates the Hamel Eq.2-6 solid terms for char:
external/recycle/axial inflows, reaction sink, size-class migration, outflow,
and the residual actually solved by the cell.  It is intentionally a ledger,
not a calibration script.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.core.cell import S_ASH, S_CHAR
from src.core.reactor import Reactor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_case_LU,
)

_SOLID_COMPONENT_LABELS = {
    S_CHAR: "char",
    S_ASH: "ash",
}


def _sum_char(arr: np.ndarray, comp: int = S_CHAR) -> float:
    return float(np.sum(np.asarray(arr, dtype=np.float64)[:, comp]))


def _positive_char(arr: np.ndarray, comp: int = S_CHAR) -> float:
    return float(np.sum(np.maximum(np.asarray(arr, dtype=np.float64)[:, comp], 0.0)))


def _negative_char(arr: np.ndarray, comp: int = S_CHAR) -> float:
    return float(np.sum(np.minimum(np.asarray(arr, dtype=np.float64)[:, comp], 0.0)))


def _class_values(arr: np.ndarray, comp: int = S_CHAR) -> list[float]:
    return [float(x) for x in np.asarray(arr, dtype=np.float64)[:, comp]]


def _cell_row(label: str, cell) -> dict[str, float | int | str]:
    mig = cell._calc_size_migration()
    out = cell._solid_outflow_rates()
    res = cell.calc_solid_balance()
    return {
        "label": label,
        "solid_state_model": str(cell.solid_state_model),
        "fresh_zu_char_kg_s": _sum_char(cell.m_solid_zu),
        "rez_char_kg_s": _sum_char(cell.m_solid_rez),
        "neighbor_in_char_kg_s": _sum_char(cell.m_solid_in),
        "auf_in_char_kg_s": _sum_char(cell.m_solid_auf_in),
        "ab_in_char_kg_s": _sum_char(cell.m_solid_ab_in),
        "reaction_char_kg_s": _sum_char(cell.R_solid),
        "reaction_char_source_positive_kg_s": _positive_char(cell.R_solid),
        "reaction_char_sink_negative_kg_s": _negative_char(cell.R_solid),
        "size_migration_char_kg_s": _sum_char(mig),
        "out_char_kg_s": _sum_char(out),
        "residual_char_kg_s": _sum_char(res),
        "max_abs_class_residual_char_kg_s": float(np.max(np.abs(np.asarray(res)[:, S_CHAR]))),
        "holdup_or_stream_char": _sum_char(cell.m_solid),
        "up_out_char_kg_s": _sum_char(cell._solid_upflow_rates()),
        "down_out_char_kg_s": _sum_char(cell._solid_downflow_rates()),
        "ash_size_migration_kg_s": _sum_char(mig, S_ASH),
    }


def _bed_transport_profile(cells: list) -> list[dict]:
    """Vertical char migration profile through bed cells and size classes."""
    rows: list[dict] = []
    for i, cell in enumerate(cells):
        mig = cell._calc_size_migration()
        res = cell.calc_solid_balance()
        up = cell._solid_upflow_rates()
        down = cell._solid_downflow_rates()
        auf_in = np.asarray(cell.m_solid_auf_in, dtype=np.float64)
        ab_in = np.asarray(cell.m_solid_ab_in, dtype=np.float64)
        up_to_above = _sum_char(up)
        down_to_below = _sum_char(down)
        auf_from_below = _sum_char(auf_in)
        ab_from_above = _sum_char(ab_in)
        auf_in_gap = 0.0
        if i > 0:
            auf_in_gap = auf_from_below - _sum_char(cells[i - 1]._solid_upflow_rates())
        ab_in_gap = 0.0
        if i + 1 < len(cells):
            ab_in_gap = ab_from_above - _sum_char(cells[i + 1]._solid_downflow_rates())
        local_throughput_ref = max(
            auf_from_below
            + ab_from_above
            + _sum_char(cell.m_solid_zu)
            + _sum_char(cell.m_solid_rez)
            + abs(_sum_char(cell.R_solid))
            + up_to_above
            + down_to_below,
            1e-12,
        )
        rows.append(
            {
                "label": f"bed[{i}]",
                "height_center_m": float(getattr(cell.geo, "h_center", float("nan"))),
                "height_bottom_m": float(getattr(cell.geo, "h_center", 0.0) - 0.5 * getattr(cell.geo, "dh", 0.0)),
                "height_top_m": float(getattr(cell.geo, "h_center", 0.0) + 0.5 * getattr(cell.geo, "dh", 0.0)),
                "d_p_classes_m": [float(x) for x in np.asarray(cell.solid.d_p_classes, dtype=np.float64)],
                "K_auf_char_1_s": _class_values(cell.K_solid_auf),
                "K_ab_char_1_s": _class_values(cell.K_solid_ab),
                "char_holdup_by_class_kg": _class_values(cell.m_solid),
                "char_up_out_by_class_kg_s": _class_values(up),
                "char_down_out_by_class_kg_s": _class_values(down),
                "char_auf_in_by_class_kg_s": _class_values(auf_in),
                "char_ab_in_by_class_kg_s": _class_values(ab_in),
                "char_reaction_by_class_kg_s": _class_values(cell.R_solid),
                "char_size_migration_by_class_kg_s": _class_values(mig),
                "char_residual_by_class_kg_s": _class_values(res),
                "char_holdup_total_kg": _sum_char(cell.m_solid),
                "char_up_out_total_kg_s": up_to_above,
                "char_down_out_total_kg_s": down_to_below,
                "char_auf_in_total_kg_s": auf_from_below,
                "char_ab_in_total_kg_s": ab_from_above,
                "char_reaction_total_kg_s": _sum_char(cell.R_solid),
                "char_size_migration_total_kg_s": _sum_char(mig),
                "char_residual_total_kg_s": _sum_char(res),
                "char_residual_norm_local": float(abs(_sum_char(res)) / local_throughput_ref),
                "char_auf_in_gap_vs_below_kg_s": float(auf_in_gap),
                "char_ab_in_gap_vs_above_kg_s": float(ab_in_gap),
            }
        )
    return rows


def _bed_transport_summary(profile: list[dict]) -> dict[str, float | bool]:
    max_res_norm = max((float(row["char_residual_norm_local"]) for row in profile), default=0.0)
    max_abs_res = max((abs(float(row["char_residual_total_kg_s"])) for row in profile), default=0.0)
    max_auf_gap = max((abs(float(row["char_auf_in_gap_vs_below_kg_s"])) for row in profile), default=0.0)
    max_ab_gap = max((abs(float(row["char_ab_in_gap_vs_above_kg_s"])) for row in profile), default=0.0)
    max_migration = max((abs(float(row["char_size_migration_total_kg_s"])) for row in profile), default=0.0)
    return {
        "max_abs_cell_char_residual_kg_s": float(max_abs_res),
        "max_cell_char_residual_norm_local": float(max_res_norm),
        "max_abs_auf_in_gap_vs_below_kg_s": float(max_auf_gap),
        "max_abs_ab_in_gap_vs_above_kg_s": float(max_ab_gap),
        "max_abs_size_migration_char_kg_s": float(max_migration),
        "bed_transport_ok_10pct_local": bool(max_res_norm <= 0.10 and max_auf_gap <= 1e-9 and max_ab_gap <= 1e-9),
    }


def _component_balance_row(label: str, cell, comp: int, comp_name: str) -> dict[str, float | str]:
    mig = cell._calc_size_migration()
    up = cell._solid_upflow_rates()
    down = cell._solid_downflow_rates()
    out = up + down
    res = cell.calc_solid_balance()
    return {
        "label": label,
        "component": comp_name,
        "holdup_kg": _sum_char(cell.m_solid, comp),
        "fresh_zu_kg_s": _sum_char(cell.m_solid_zu, comp),
        "rez_kg_s": _sum_char(cell.m_solid_rez, comp),
        "neighbor_in_kg_s": _sum_char(cell.m_solid_in, comp),
        "auf_in_kg_s": _sum_char(cell.m_solid_auf_in, comp),
        "ab_in_kg_s": _sum_char(cell.m_solid_ab_in, comp),
        "reaction_kg_s": _sum_char(cell.R_solid, comp),
        "size_migration_kg_s": _sum_char(mig, comp),
        "up_out_kg_s": _sum_char(up, comp),
        "down_out_kg_s": _sum_char(down, comp),
        "out_kg_s": _sum_char(out, comp),
        "residual_kg_s": _sum_char(res, comp),
        "K_auf_mean_1_s": float(np.mean(np.asarray(cell.K_solid_auf, dtype=np.float64)[:, comp])),
        "K_ab_mean_1_s": float(np.mean(np.asarray(cell.K_solid_ab, dtype=np.float64)[:, comp])),
    }


def _bed_freeboard_interface_component_balance(reactor: Reactor) -> list[dict[str, float | str]]:
    """Component balance at the top bed/freeboard interface for char/ash."""
    if not reactor.cells:
        return []
    top = reactor.cells[-1]
    rows: list[dict[str, float | str]] = []
    for comp, name in _SOLID_COMPONENT_LABELS.items():
        row = _component_balance_row("bed_top", top, comp, name)
        if getattr(reactor, "freeboard_cells", None):
            fb0 = reactor.freeboard_cells[0]
            row["freeboard_bottom_return_kg_s"] = _sum_char(fb0._solid_downflow_rates(), comp)
            row["bed_top_ab_in_gap_vs_freeboard_return_kg_s"] = float(
                row["ab_in_kg_s"] - row["freeboard_bottom_return_kg_s"]
            )
        else:
            row["freeboard_bottom_return_kg_s"] = 0.0
            row["bed_top_ab_in_gap_vs_freeboard_return_kg_s"] = 0.0
        rows.append(row)
    if len(rows) == 2:
        total = {
            "label": "bed_top",
            "component": "char+ash",
        }
        for key in rows[0]:
            if key in {"label", "component"}:
                continue
            total[key] = float(rows[0][key]) + float(rows[1][key])
        rows.append(total)
    return rows


def _top_bed_solid_local_sensitivity(reactor: Reactor) -> list[dict[str, float | str]]:
    """Finite-difference top-bed solid residual sensitivity to local holdup."""
    if not reactor.cells:
        return []
    cell = reactor.cells[-1]
    rows: list[dict[str, float | str]] = []
    base_res = cell.calc_solid_balance()
    for comp, name in _SOLID_COMPONENT_LABELS.items():
        values = np.asarray(cell.m_solid, dtype=np.float64)[:, comp]
        if values.size == 0:
            continue
        class_idx = int(np.argmax(np.abs(np.asarray(base_res, dtype=np.float64)[:, comp])))
        base_m = float(values[class_idx])
        h = max(1e-6 * max(abs(base_m), 1.0), 1e-8)
        cell.m_solid[class_idx, comp] = base_m + h
        try:
            pert_res = cell.calc_solid_balance()
        finally:
            cell.m_solid[class_idx, comp] = base_m
        deriv = float((pert_res[class_idx, comp] - base_res[class_idx, comp]) / h)
        k_sum = float(cell.K_solid_auf[class_idx, comp] + cell.K_solid_ab[class_idx, comp])
        residual = float(base_res[class_idx, comp])
        newton_est = -residual / deriv if abs(deriv) > 1e-14 else float("nan")
        rows.append(
            {
                "component": name,
                "size_class": int(class_idx),
                "holdup_kg": base_m,
                "residual_kg_s": residual,
                "d_residual_d_holdup_1_s": deriv,
                "minus_K_sum_1_s": -k_sum,
                "local_newton_delta_holdup_kg": newton_est,
            }
        )
    return rows


def _freeboard_closure_transport_profile(result: dict) -> list[dict[str, float | int]]:
    """Read analytical freeboard solid transport from result-builder fields."""
    up = [float(v) for v in result.get("freeboard_entrained_char_kg_s", []) or []]
    down = [float(v) for v in result.get("freeboard_entrained_return_char_profile_kg_s", []) or []]
    hold = [float(v) for v in result.get("freeboard_solid_holdup_char_profile_kg", []) or []]
    hold_before = [float(v) for v in result.get("freeboard_solid_holdup_char_before_profile_kg", []) or []]
    sink = [float(v) for v in result.get("freeboard_explicit_char_sink_applied_profile_kg", []) or []]
    z = [float(v) for v in result.get("freeboard_axial_z_m", []) or []]
    xi = [float(v) for v in result.get("freeboard_axial_xi", []) or []]
    carry = [float(v) for v in result.get("freeboard_carry_ratio_profile", []) or []]
    n = max(len(up), len(down), len(hold), len(hold_before), len(sink), len(z), len(xi), len(carry))
    rows: list[dict[str, float | int]] = []
    for i in range(n):
        rows.append(
            {
                "segment": i,
                "z_m": z[i] if i < len(z) else float("nan"),
                "xi": xi[i] if i < len(xi) else float("nan"),
                "char_holdup_before_kg": hold_before[i] if i < len(hold_before) else 0.0,
                "char_holdup_kg": hold[i] if i < len(hold) else 0.0,
                "char_up_kg_s": up[i] if i < len(up) else 0.0,
                "char_return_kg_s": down[i] if i < len(down) else 0.0,
                "char_sink_applied_kg": sink[i] if i < len(sink) else 0.0,
                "carry_ratio": carry[i] if i < len(carry) else float("nan"),
            }
        )
    return rows


def _freeboard_closure_transport_summary(reactor: Reactor, result: dict) -> dict[str, float | bool]:
    profile = _freeboard_closure_transport_profile(result)
    current_bed_top_up = _sum_char(reactor.cells[-1]._solid_upflow_rates()) if reactor.cells else 0.0
    bed_top_up = float(
        result.get(
            "freeboard_bed_top_up_char_kg_s",
            current_bed_top_up,
        )
        or 0.0
    )
    exit_char = float(result.get("freeboard_entrained_exit_char_kg_s", 0.0) or 0.0)
    return_char = float(result.get("freeboard_entrained_return_char_kg_s", 0.0) or 0.0)
    cyclone_capture = float(result.get("freeboard_cyclone_capture_char_kg_s", 0.0) or 0.0)
    bottom_recycle = _sum_char(reactor.cells[0].m_solid_rez) if reactor.cells else 0.0
    closure_sink = float(sum(max(float(row["char_sink_applied_kg"]), 0.0) for row in profile))
    explicit_top_up = (
        _sum_char(reactor.freeboard_cells[-1]._solid_upflow_rates())
        if getattr(reactor, "freeboard_cells", None)
        else 0.0
    )
    explicit_bottom_return = (
        _sum_char(reactor.freeboard_cells[0]._solid_downflow_rates())
        if getattr(reactor, "freeboard_cells", None)
        else 0.0
    )
    cyclone_in = (
        _sum_char(reactor.cyclone_cell.m_solid_auf_in)
        if getattr(reactor, "cyclone_cell", None) is not None
        else 0.0
    )
    return_leg_in = (
        _sum_char(reactor.return_leg_cell.m_solid_ab_in)
        if getattr(reactor, "return_leg_cell", None) is not None
        else 0.0
    )
    closure_transport_gap = bed_top_up - return_char - exit_char - closure_sink
    return {
        "freeboard_active": bool(result.get("freeboard_active", False)),
        "bed_top_up_char_kg_s": float(bed_top_up),
        "current_bed_top_up_char_kg_s": float(current_bed_top_up),
        "gap_current_bed_top_vs_closure_bed_top_kg_s": float(current_bed_top_up - bed_top_up),
        "pre_result_refresh_bed_top_up_char_kg_s": float(
            result.get("freeboard_pre_result_refresh_bed_top_up_char_kg_s", bed_top_up) or 0.0
        ),
        "pre_result_refresh_gap_current_bed_top_kg_s": float(
            result.get("freeboard_pre_result_refresh_gap_current_bed_top_kg_s", 0.0) or 0.0
        ),
        "post_result_refresh_gap_current_bed_top_kg_s": float(
            result.get("freeboard_post_result_refresh_gap_current_bed_top_kg_s", current_bed_top_up - bed_top_up) or 0.0
        ),
        "freeboard_return_char_kg_s": float(return_char),
        "freeboard_exit_char_kg_s": float(exit_char),
        "freeboard_closure_sink_char_kg_equiv": float(closure_sink),
        "freeboard_closure_transport_gap_kg_s": float(closure_transport_gap),
        "cyclone_capture_char_kg_s": float(cyclone_capture),
        "system_escape_after_cyclone_char_kg_s": float(max(exit_char - cyclone_capture, 0.0)),
        "bottom_recycle_char_kg_s": float(bottom_recycle),
        "explicit_freeboard_top_up_char_kg_s": float(explicit_top_up),
        "explicit_freeboard_bottom_return_char_kg_s": float(explicit_bottom_return),
        "cyclone_in_char_kg_s": float(cyclone_in),
        "return_leg_in_char_kg_s": float(return_leg_in),
        "gap_result_exit_vs_explicit_top_kg_s": float(exit_char - explicit_top_up),
        "gap_result_return_vs_explicit_bottom_kg_s": float(return_char - explicit_bottom_return),
        "gap_cyclone_capture_vs_return_leg_in_kg_s": float(cyclone_capture - return_leg_in),
        "freeboard_profile_nonempty": bool(profile),
    }


def _build(mode: str) -> tuple[Reactor, dict]:
    case = load_case_LU()
    if mode == "phase2":
        return Reactor(build_phase2_htw_lu_freeboard_reactor_config(case)), dict(PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS)
    return Reactor(build_phase1_htw_lu_global_nr_reactor_config(case)), dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)


def summarize_reactor_char_ledger(reactor: Reactor, result: dict, *, mode: str, solve_kwargs: dict) -> dict:
    """Build a char ledger from an already solved reactor state."""
    bed_rows = [_cell_row(f"bed[{i}]", cell) for i, cell in enumerate(reactor.cells)]
    bed_transport_profile = _bed_transport_profile(list(reactor.cells))
    bed_transport_summary = _bed_transport_summary(bed_transport_profile)
    freeboard_rows = [_cell_row(f"freeboard[{i}]", cell) for i, cell in enumerate(getattr(reactor, "freeboard_cells", []))]
    side_rows = []
    for name in ("cyclone_cell", "return_leg_cell"):
        cell = getattr(reactor, name, None)
        if cell is not None:
            side_rows.append(_cell_row(name, cell))

    all_rows = bed_rows + freeboard_rows + side_rows
    top_char_residual_rows = sorted(
        all_rows,
        key=lambda row: abs(float(row["residual_char_kg_s"])),
        reverse=True,
    )[:5]
    fresh_char = float(sum(row["fresh_zu_char_kg_s"] for row in all_rows))
    reaction_char = float(sum(row["reaction_char_kg_s"] for row in all_rows))
    reaction_char_source = float(sum(row["reaction_char_source_positive_kg_s"] for row in all_rows))
    reaction_char_sink = float(sum(row["reaction_char_sink_negative_kg_s"] for row in all_rows))
    migration_char = float(sum(row["size_migration_char_kg_s"] for row in all_rows))
    residual_char = float(sum(row["residual_char_kg_s"] for row in all_rows))
    max_cell_res = float(max((abs(float(row["residual_char_kg_s"])) for row in all_rows), default=0.0))
    ash_migration = float(sum(row["ash_size_migration_kg_s"] for row in all_rows))

    top = reactor.cells[-1]
    top_up_char = _sum_char(top._solid_upflow_rates())
    bottom_recycle = _sum_char(reactor.cells[0].m_solid_rez)
    system_exit_proxy = max(top_up_char - bottom_recycle, 0.0)
    outer_history = list(result.get("nr_outer_history", []) or [])
    dT_refresh_vals = [
        float(item.get("dT_refresh", 0.0) or 0.0)
        for item in outer_history
        if isinstance(item, dict)
    ]
    outer_seed_reset_count = sum(
        1
        for item in outer_history
        if isinstance(item, dict) and bool(item.get("lambda_seed_reset_by_refresh", False))
    )

    return {
        "mode": mode,
        "solve_kwargs": solve_kwargs,
        "converged": bool(result.get("converged")),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
        "rms_scaled_solid_final": float(result.get("rms_scaled_solid_final", float("nan"))),
        "rms_scaled_component_max_final": float(result.get("rms_scaled_component_max_final", float("nan"))),
        "max_abs_scaled_final": float(result.get("max_abs_scaled_final", float("nan"))),
        "max_abs_scaled_solid_final": float(result.get("max_abs_scaled_solid_final", float("nan"))),
        "carbon_conv": float(result.get("carbon_conv", float("nan"))),
        "nr_outer_history": outer_history,
        "nr_outer_refresh_jump_summary": {
            "outer_iters": int(len(outer_history)),
            "max_dT_refresh_K": float(max(dT_refresh_vals)) if dT_refresh_vals else 0.0,
            "avg_dT_refresh_K": float(np.mean(dT_refresh_vals)) if dT_refresh_vals else 0.0,
            "lambda_seed_reset_count": int(outer_seed_reset_count),
        },
        "system_char_summary": {
            "fresh_char_kg_s": fresh_char,
            "reaction_char_kg_s": reaction_char,
            "reaction_char_source_positive_kg_s": reaction_char_source,
            "reaction_char_sink_negative_kg_s": reaction_char_sink,
            "size_migration_char_kg_s": migration_char,
            "size_migration_ash_kg_s": ash_migration,
            "residual_char_kg_s": residual_char,
            "max_abs_cell_char_residual_kg_s": max_cell_res,
            "top_up_char_kg_s": top_up_char,
            "bottom_recycle_char_kg_s": bottom_recycle,
            "system_exit_proxy_top_up_minus_recycle_kg_s": system_exit_proxy,
            "fresh_minus_reaction_minus_exit_proxy_kg_s": fresh_char + reaction_char - system_exit_proxy,
        },
        "top_char_residual_rows": [
            {
                "label": str(row["label"]),
                "solid_state_model": str(row["solid_state_model"]),
                "residual_char_kg_s": float(row["residual_char_kg_s"]),
                "fresh_zu_char_kg_s": float(row["fresh_zu_char_kg_s"]),
                "auf_in_char_kg_s": float(row["auf_in_char_kg_s"]),
                "ab_in_char_kg_s": float(row["ab_in_char_kg_s"]),
                "reaction_char_kg_s": float(row["reaction_char_kg_s"]),
                "reaction_char_source_positive_kg_s": float(row["reaction_char_source_positive_kg_s"]),
                "reaction_char_sink_negative_kg_s": float(row["reaction_char_sink_negative_kg_s"]),
                "out_char_kg_s": float(row["out_char_kg_s"]),
                "holdup_or_stream_char": float(row["holdup_or_stream_char"]),
            }
            for row in top_char_residual_rows
        ],
        "bed_rows": bed_rows,
        "bed_transport_profile": bed_transport_profile,
        "bed_transport_summary": bed_transport_summary,
        "bed_freeboard_interface_component_balance": _bed_freeboard_interface_component_balance(reactor),
        "top_bed_solid_local_sensitivity": _top_bed_solid_local_sensitivity(reactor),
        "freeboard_closure_transport_profile": _freeboard_closure_transport_profile(result),
        "freeboard_closure_transport_summary": _freeboard_closure_transport_summary(reactor, result),
        "freeboard_rows": freeboard_rows,
        "side_rows": side_rows,
    }


def run(mode: str, *, max_global_iter: int | None = None, tol_global: float | None = None) -> dict:
    reactor, solve_kwargs = _build(mode)
    if max_global_iter is not None:
        solve_kwargs["max_global_iter"] = int(max_global_iter)
    if tol_global is not None:
        solve_kwargs["tol_global"] = float(tol_global)
    result = reactor.solve(**solve_kwargs)
    payload = summarize_reactor_char_ledger(reactor, result, mode=mode, solve_kwargs=solve_kwargs)
    if mode == "phase2" and reactor._use_explicit_freeboard_solver_graph():
        refreshed = reactor._refresh_explicit_freeboard_transport_from_closure(preserve_gas_state=True)
        if refreshed is not None:
            refreshed_result = reactor._build_exit_summary()
            payload["refreshed_freeboard_closure_transport_profile"] = _freeboard_closure_transport_profile(
                refreshed_result
            )
            payload["refreshed_freeboard_closure_transport_summary"] = _freeboard_closure_transport_summary(
                reactor,
                refreshed_result,
            )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("phase1", "phase2"), default="phase1")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--max-global-iter", type=int, default=None)
    parser.add_argument("--tol-global", type=float, default=None)
    args = parser.parse_args()

    payload = run(args.mode, max_global_iter=args.max_global_iter, tol_global=args.tol_global)
    if args.json:
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0 if payload["converged"] else 1

    print(f"LU char mass conservation audit ({args.mode})")
    print("=" * 88)
    print(
        f"converged={payload['converged']} rms={payload['rms_scaled_final']:.4e} "
        f"carbon_conv={100.0 * payload['carbon_conv']:.2f}%"
    )
    print("system summary:")
    for key, value in payload["system_char_summary"].items():
        print(f"  {key:48s} {float(value): .6e}")
    print("\ncell rows:")
    cols = [
        "label",
        "fresh_zu_char_kg_s",
        "rez_char_kg_s",
        "auf_in_char_kg_s",
        "ab_in_char_kg_s",
        "reaction_char_kg_s",
        "reaction_char_source_positive_kg_s",
        "reaction_char_sink_negative_kg_s",
        "size_migration_char_kg_s",
        "out_char_kg_s",
        "residual_char_kg_s",
    ]
    print(" ".join(f"{c:>16s}" for c in cols))
    for row in payload["bed_rows"] + payload["freeboard_rows"] + payload["side_rows"]:
        vals = []
        for col in cols:
            vals.append(f"{str(row[col]):>16s}" if col == "label" else f"{float(row[col]):16.6e}")
        print(" ".join(vals))
    print("\nbed char vertical transport profile:")
    prof_cols = [
        "label",
        "height_center_m",
        "char_holdup_total_kg",
        "char_up_out_total_kg_s",
        "char_down_out_total_kg_s",
        "char_reaction_total_kg_s",
        "char_residual_total_kg_s",
    ]
    print(" ".join(f"{c:>18s}" for c in prof_cols))
    for row in payload["bed_transport_profile"]:
        vals = []
        for col in prof_cols:
            vals.append(f"{str(row[col]):>18s}" if col == "label" else f"{float(row[col]):18.6e}")
        print(" ".join(vals))
    fb_summary = payload.get("freeboard_closure_transport_summary", {})
    if fb_summary.get("freeboard_active"):
        print("\nfreeboard closure transport summary:")
        for key, value in fb_summary.items():
            if key == "freeboard_active" or key == "freeboard_profile_nonempty":
                print(f"  {key:48s} {bool(value)}")
            else:
                print(f"  {key:48s} {float(value): .6e}")
    interface_rows = payload.get("bed_freeboard_interface_component_balance", [])
    if interface_rows:
        print("\nbed/freeboard interface component balance:")
        interface_cols = [
            "component",
            "auf_in_kg_s",
            "ab_in_kg_s",
            "freeboard_bottom_return_kg_s",
            "bed_top_ab_in_gap_vs_freeboard_return_kg_s",
            "reaction_kg_s",
            "up_out_kg_s",
            "down_out_kg_s",
            "residual_kg_s",
        ]
        print(" ".join(f"{c:>24s}" for c in interface_cols))
        for row in interface_rows:
            vals = []
            for col in interface_cols:
                vals.append(f"{str(row[col]):>24s}" if col == "component" else f"{float(row[col]):24.6e}")
            print(" ".join(vals))
    sensitivity_rows = payload.get("top_bed_solid_local_sensitivity", [])
    if sensitivity_rows:
        print("\ntop-bed solid local sensitivity:")
        sens_cols = [
            "component",
            "size_class",
            "holdup_kg",
            "residual_kg_s",
            "d_residual_d_holdup_1_s",
            "minus_K_sum_1_s",
            "local_newton_delta_holdup_kg",
        ]
        print(" ".join(f"{c:>26s}" for c in sens_cols))
        for row in sensitivity_rows:
            vals = []
            for col in sens_cols:
                if col == "component":
                    vals.append(f"{str(row[col]):>26s}")
                elif col == "size_class":
                    vals.append(f"{int(row[col]):26d}")
                else:
                    vals.append(f"{float(row[col]):26.6e}")
            print(" ".join(vals))
    refreshed_fb_summary = payload.get("refreshed_freeboard_closure_transport_summary", {})
    if refreshed_fb_summary.get("freeboard_active"):
        print("\nrefreshed freeboard closure transport summary (solid-only, gas/T preserved):")
        for key, value in refreshed_fb_summary.items():
            if key == "freeboard_active" or key == "freeboard_profile_nonempty":
                print(f"  {key:48s} {bool(value)}")
            else:
                print(f"  {key:48s} {float(value): .6e}")
    return 0 if payload["converged"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
