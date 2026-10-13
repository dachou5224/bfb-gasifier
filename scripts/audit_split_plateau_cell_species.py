#!/usr/bin/env python3
"""平台期 cell×物种 phase-split 主导项审计（handoff P0）。

对齐 global NR split merit：``dense[pos]-bubble[pos]``（仅 MAJOR_NR_GAS_SPECIES）。
在同构限幅关 + bed1 preinner 栈上跑 o8，报告：
1. 末态全局 / 按床 / cell×物种方差贡献
2. 末态再走一步限幅 Newton + fastox 的 energy↓/split↑ 交易解剖

用法：
  python3 scripts/audit_split_plateau_cell_species.py
  python3 scripts/audit_split_plateau_cell_species.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.species import GAS_SPECIES
from src.solvers.global_nr_solver import (
    _clip_dx,
    _line_search_merit,
    _residual_gas_phase_metrics,
    _residual_group_metrics,
    _rms_norm,
    _solve_linear_step,
    build_equation_scales,
    build_jacobian_fd,
    global_residual,
    pack_reactor,
    unpack_reactor,
)
from src.solvers.nr_indexing import (
    gas_mode_for_nr,
    has_temperature_var,
    install_main_nr_gas_indexing_on_cells,
    main_nr_gas_idx_for_cell,
    main_nr_gas_phase_split_positions_for_cell,
    n_gas_var,
    n_var,
    temperature_offset,
)
from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

ISOMORPH_CFG = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}


def _scaled_state(reactor: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Any]]:
    cells = list(reactor._solver_cells_for_nr())
    install_main_nr_gas_indexing_on_cells(cells, reactor.config)
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    x = pack_reactor(cells)
    F = global_residual(x, cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    return x, F_hat, scale, cells


def _iter_split_entries(F_hat: np.ndarray, cells: list[Any]) -> list[dict[str, Any]]:
    """All (cell, species) major phase-split residual terms with abs and sq."""
    entries: list[dict[str, Any]] = []
    offset = 0
    for i, cell in enumerate(cells):
        nv = int(n_var(cell))
        ng = int(n_gas_var(cell))
        if gas_mode_for_nr(cell) == "two_phase" and ng > 0:
            half = ng // 2
            dense = F_hat[offset : offset + half]
            bubble = F_hat[offset + half : offset + ng]
            gas_idx = main_nr_gas_idx_for_cell(cell)
            for p in main_nr_gas_phase_split_positions_for_cell(cell):
                sp = str(GAS_SPECIES[int(gas_idx[int(p)])])
                d_val = float(dense[int(p)])
                b_val = float(bubble[int(p)])
                diff = d_val - b_val
                entries.append(
                    {
                        "cell": int(i),
                        "cell_type": str(getattr(cell, "cell_type", "bed")),
                        "species": sp,
                        "T_K": float(getattr(cell, "T", float("nan"))),
                        "dense": d_val,
                        "bubble": b_val,
                        "dense_minus_bubble": float(diff),
                        "abs_split": float(abs(diff)),
                        "sq": float(diff * diff),
                    }
                )
        offset += nv
    return entries


def _rank_entries(entries: list[dict[str, Any]], *, top_n: int = 20) -> list[dict[str, Any]]:
    total_sq = float(sum(float(e["sq"]) for e in entries))
    ranked = sorted(entries, key=lambda e: float(e["abs_split"]), reverse=True)
    out: list[dict[str, Any]] = []
    for e in ranked[: int(top_n)]:
        row = dict(e)
        row["variance_fraction"] = float(e["sq"] / max(total_sq, 1.0e-30))
        out.append(row)
    return out


def _aggregate_by_species(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    total_sq = float(sum(float(e["sq"]) for e in entries))
    buckets: dict[str, dict[str, float]] = {}
    for e in entries:
        sp = str(e["species"])
        b = buckets.setdefault(sp, {"sq": 0.0, "max_abs": 0.0, "n": 0.0})
        b["sq"] += float(e["sq"])
        b["max_abs"] = max(b["max_abs"], float(e["abs_split"]))
        b["n"] += 1.0
    rows = []
    for sp, b in buckets.items():
        rows.append(
            {
                "species": sp,
                "variance_fraction": float(b["sq"] / max(total_sq, 1.0e-30)),
                "rms": float(np.sqrt(b["sq"] / max(b["n"], 1.0))),
                "max_abs": float(b["max_abs"]),
                "n_terms": int(b["n"]),
            }
        )
    rows.sort(key=lambda r: float(r["variance_fraction"]), reverse=True)
    return rows


def _aggregate_by_cell(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    total_sq = float(sum(float(e["sq"]) for e in entries))
    buckets: dict[int, dict[str, Any]] = {}
    for e in entries:
        i = int(e["cell"])
        b = buckets.setdefault(
            i,
            {
                "cell": i,
                "cell_type": e["cell_type"],
                "T_K": e["T_K"],
                "sq": 0.0,
                "max_abs": 0.0,
                "n": 0.0,
                "worst_species": e["species"],
                "worst_abs": 0.0,
            },
        )
        b["sq"] += float(e["sq"])
        b["max_abs"] = max(float(b["max_abs"]), float(e["abs_split"]))
        b["n"] += 1.0
        if float(e["abs_split"]) >= float(b["worst_abs"]):
            b["worst_abs"] = float(e["abs_split"])
            b["worst_species"] = e["species"]
    rows = []
    for b in buckets.values():
        n = max(float(b["n"]), 1.0)
        rows.append(
            {
                "cell": int(b["cell"]),
                "cell_type": b["cell_type"],
                "T_K": float(b["T_K"]),
                "variance_fraction": float(b["sq"] / max(total_sq, 1.0e-30)),
                "rms": float(np.sqrt(float(b["sq"]) / n)),
                "max_abs": float(b["max_abs"]),
                "worst_species": b["worst_species"],
                "n_terms": int(b["n"]),
            }
        )
    rows.sort(key=lambda r: float(r["variance_fraction"]), reverse=True)
    return rows


def _metrics(F_hat: np.ndarray, cells: list[Any]) -> dict[str, float]:
    g = _residual_group_metrics(F_hat, cells)
    p = _residual_gas_phase_metrics(F_hat, cells)
    merit = _line_search_merit(
        F_hat, cells, use_gas_phase_split=True, use_energy=True
    )
    return {
        "rms": float(_rms_norm(F_hat)),
        "split": float(p["gas_phase_split_rms"]),
        "energy": float(g["energy_rms"]),
        "gas": float(g["gas_rms"]),
        "solid": float(g["solid_rms"]),
        "merit": float(merit),
        "driver": (
            "energy"
            if g["energy_rms"] >= max(p["gas_phase_split_rms"], _rms_norm(F_hat)) - 1e-15
            else (
                "split"
                if p["gas_phase_split_rms"] >= _rms_norm(F_hat) - 1e-15
                else "rms"
            )
        ),
    }


def _decompose_snapshot(F_hat: np.ndarray, cells: list[Any], *, top_n: int = 15) -> dict[str, Any]:
    entries = _iter_split_entries(F_hat, cells)
    return {
        "metrics": _metrics(F_hat, cells),
        "by_species": _aggregate_by_species(entries),
        "by_cell": _aggregate_by_cell(entries),
        "top_cell_species": _rank_entries(entries, top_n=top_n),
        "n_split_terms": int(len(entries)),
    }


def _apply_fastox_trial(cells: list[Any], x_trial: np.ndarray) -> np.ndarray:
    unpack_reactor(np.asarray(x_trial, dtype=np.float64), cells)
    fuel_type = str(getattr(cells[0], "fuel_type", "coal") if cells else "coal")
    project_bed_holdup_fast_oxidation_closure(
        cells,
        SimpleNamespace(
            vorab_transport_x0_fast_oxidation_closure_thesis=True,
            fuel_type=fuel_type,
        ),
        phase_update="local_participating",
    )
    return pack_reactor(cells)


def _energy_split_trade_probe(reactor: Any) -> dict[str, Any]:
    """From plateau state: one Newton step, find λ that lowers energy but raises merit."""
    cfg = reactor.config
    x, F_hat, scale, cells = _scaled_state(reactor)
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    F0 = F_hat * scale
    J, _meta = build_jacobian_fd(
        x,
        F0,
        cells,
        reactor._apply_all_bc_for_nr,
        scale,
        strategy="hybrid_band_plus_side_elements_structured",
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        affected_residual_cells_fn=reactor._affected_nr_residual_cells,
        connectivity_graph=reactor.connectivity_graph,
    )
    # Temperature DOFs for structured T rescale path
    t_dofs: list[int] = []
    off = 0
    for cell in cells:
        nv = n_var(cell)
        if has_temperature_var(cell):
            t_dofs.append(off + int(temperature_offset(cell)))
        off += nv
    dx, method, _reg, _diag = _solve_linear_step(
        J,
        -F_hat,
        structured_jacobian=None,
        linear_solver_backend="sparse_direct_fallback",
        t_dof_indices=np.asarray(t_dofs, dtype=np.intp) if t_dofs else None,
        t_max_step_K=600.0,
    )
    dx_clipped = _clip_dx(
        dx, cells, ref_gas, ref_solid, t_step_limit_K=350.0
    )
    base = _decompose_snapshot(F_hat, cells, top_n=12)
    base_m = base["metrics"]

    trials: list[dict[str, Any]] = []
    best_trade: dict[str, Any] | None = None
    lam = 1.0
    for _ in range(12):
        x_trial = np.asarray(x, dtype=np.float64) + lam * dx_clipped
        # nonneg / T bounds (mirror LS)
        off = 0
        for cell in cells:
            nv = n_var(cell)
            if has_temperature_var(cell):
                tl = temperature_offset(cell)
                x_trial[off : off + tl] = np.maximum(x_trial[off : off + tl], 0.0)
                t_min = float(getattr(cell, "_nr_temperature_min_K", 300.0))
                t_max = float(getattr(cell, "_nr_temperature_max_K", 2500.0))
                x_trial[off + tl] = float(np.clip(x_trial[off + tl], t_min, t_max))
            else:
                x_trial[off : off + nv] = np.maximum(x_trial[off : off + nv], 0.0)
            off += nv
        x_trial = _apply_fastox_trial(cells, x_trial)
        F_t = global_residual(x_trial, cells, reactor._apply_all_bc_for_nr)
        F_hat_t = np.asarray(F_t, dtype=np.float64) / scale
        snap = _decompose_snapshot(F_hat_t, cells, top_n=12)
        m = snap["metrics"]
        row = {
            "lam": float(lam),
            "metrics": m,
            "d_merit": float(m["merit"] - base_m["merit"]),
            "d_split": float(m["split"] - base_m["split"]),
            "d_energy": float(m["energy"] - base_m["energy"]),
            "top_cell_species": snap["top_cell_species"][:8],
            "by_species": snap["by_species"][:6],
        }
        trials.append(row)
        is_trade = (
            float(m["energy"]) < float(base_m["energy"])
            and float(m["merit"]) >= float(base_m["merit"])
            and float(m["split"]) > float(base_m["split"])
        )
        if is_trade and (
            best_trade is None or float(m["energy"]) < float(best_trade["metrics"]["energy"])
        ):
            best_trade = row
        lam *= 0.5
        if lam < 1.0 / 1024.0:
            break

    # restore plateau state after probes
    unpack_reactor(x, cells)
    reactor._apply_all_bc_for_nr()

    delta_leaders: list[dict[str, Any]] = []
    if best_trade is not None:
        unpack_reactor(x, cells)
        F_base = global_residual(x, cells, reactor._apply_all_bc_for_nr)
        F_hat_base = np.asarray(F_base, dtype=np.float64) / scale
        base_map = {
            (e["cell"], e["species"]): e for e in _iter_split_entries(F_hat_base, cells)
        }
        x_trial = np.asarray(x, dtype=np.float64) + float(best_trade["lam"]) * dx_clipped
        x_trial = _apply_fastox_trial(cells, x_trial)
        F_t = global_residual(x_trial, cells, reactor._apply_all_bc_for_nr)
        F_hat_t = np.asarray(F_t, dtype=np.float64) / scale
        for e in _iter_split_entries(F_hat_t, cells):
            key = (e["cell"], e["species"])
            b = base_map.get(key)
            b_abs = float(b["abs_split"]) if b else 0.0
            d_abs = float(e["abs_split"]) - b_abs
            if abs(d_abs) < 1e-6:
                continue
            delta_leaders.append(
                {
                    "cell": e["cell"],
                    "cell_type": e["cell_type"],
                    "species": e["species"],
                    "abs_base": b_abs,
                    "abs_trial": float(e["abs_split"]),
                    "d_abs": float(d_abs),
                    "dense_minus_bubble_trial": float(e["dense_minus_bubble"]),
                }
            )
        delta_leaders.sort(key=lambda r: float(r["d_abs"]), reverse=True)
        unpack_reactor(x, cells)
        reactor._apply_all_bc_for_nr()

    return {
        "linear_backend": str(method),
        "base": base,
        "n_trials": int(len(trials)),
        "best_energy_split_trade": best_trade,
        "delta_abs_split_leaders": delta_leaders[:15],
        "all_trials_metrics": [
            {
                "lam": t["lam"],
                "merit": t["metrics"]["merit"],
                "split": t["metrics"]["split"],
                "energy": t["metrics"]["energy"],
                "d_merit": t["d_merit"],
                "d_split": t["d_split"],
                "d_energy": t["d_energy"],
            }
            for t in trials
        ],
    }


def run_audit(*, max_global_iter: int) -> dict[str, Any]:
    reactor = build_phase2_reactor(
        ladder_level="p2",
        apply_kinetics=True,
        cfg_overrides=dict(ISOMORPH_CFG),
    )
    run_phase2_init(reactor)
    out, wall = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_global_iter))
    )
    _x, F_hat, _scale, cells = _scaled_state(reactor)
    final = _decompose_snapshot(F_hat, cells, top_n=20)
    trade = _energy_split_trade_probe(reactor)
    eg = out.get("exit_gas_dry") or {}
    return {
        "note": (
            "Plateau cell×species split leaders; isomorphic stack; "
            "trade probe = one clipped Newton + fastox λ ladder"
        ),
        "cfg": dict(ISOMORPH_CFG),
        "wall_s": float(wall),
        "solve": {
            "rms": float(out.get("rms_scaled_final") or 0.0),
            "split": float(out.get("rms_scaled_gas_phase_split_final") or 0.0),
            "energy": float(out.get("rms_scaled_energy_final") or 0.0),
            "CO": float(eg.get("CO") or 0.0),
            "H2": float(eg.get("H2") or 0.0),
            "ls_fail": int((out.get("nr_counts") or {}).get("line_search_failures") or 0),
            "n_acc": len(
                [x for x in (out.get("nr_accepted_lambda_history") or []) if x is not None]
            ),
        },
        "final": final,
        "energy_split_trade_probe": trade,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=data_path("phase2_split_plateau_cell_species_o8.json"),
    )
    args = parser.parse_args()
    payload = run_audit(max_global_iter=int(args.max_global_iter))
    write_json(args.output, payload)

    final = payload["final"]
    trade = payload["energy_split_trade_probe"]
    summary = {
        "output": str(args.output),
        "solve": payload["solve"],
        "final_metrics": final["metrics"],
        "top_species": final["by_species"][:6],
        "top_cells": [
            {
                "cell": r["cell"],
                "type": r["cell_type"],
                "var_frac": round(float(r["variance_fraction"]), 3),
                "worst": r["worst_species"],
                "rms": round(float(r["rms"]), 5),
            }
            for r in final["by_cell"][:6]
        ],
        "top_cell_species": [
            {
                "cell": r["cell"],
                "sp": r["species"],
                "abs": round(float(r["abs_split"]), 5),
                "var_frac": round(float(r["variance_fraction"]), 3),
            }
            for r in final["top_cell_species"][:10]
        ],
        "trade": None
        if trade.get("best_energy_split_trade") is None
        else {
            "lam": trade["best_energy_split_trade"]["lam"],
            "d_merit": round(float(trade["best_energy_split_trade"]["d_merit"]), 5),
            "d_split": round(float(trade["best_energy_split_trade"]["d_split"]), 5),
            "d_energy": round(float(trade["best_energy_split_trade"]["d_energy"]), 5),
            "delta_leaders": [
                {
                    "cell": r["cell"],
                    "sp": r["species"],
                    "d_abs": round(float(r["d_abs"]), 5),
                    "abs_trial": round(float(r["abs_trial"]), 5),
                }
                for r in trade.get("delta_abs_split_leaders", [])[:8]
            ],
        },
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
