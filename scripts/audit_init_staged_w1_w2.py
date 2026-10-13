#!/usr/bin/env python3
"""分阶段 init 审计：对比 W1/W2 在 init_precalc 各步后的 bed0 状态与 NR 初值残差。

Usage:
    python3 scripts/audit_init_staged_w1_w2.py
    python3 scripts/audit_init_staged_w1_w2.py --cases CASE_HTW_WESSELING_2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from _lu_harness import build_phase2_reactor, data_path, ensure_repo_path, write_json

ensure_repo_path()

from src.core.connectivity import cell_total_solid_holdup
from src.core.connectivity.solid_transport import bed_upflow_coefficient
from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX
from src.core.cell_hydrodynamics import _resolve_u_d_closure
from src.physics.phase_fractions import calc_visible_bubble_fraction
from src.solvers.global_nr_solver import (
    _line_search_merit,
    _residual_group_metrics,
    build_equation_scales,
    global_residual,
    pack_reactor,
)
from src.workflow.steps import init_precalc_step as init_step

STAGES: tuple[str, ...] = (
    "1_seed_T_and_upstream",
    "2_vorab_x0_snapshot",
    "3_holdup_chain_and_reactive",
    "4_bottom_zone_preprojections",
    "5_abgleich_joint_preprojection",
    "6_full_init",
)

SPECIES_TRACE = ("CO", "CO2", "H2", "CH4", "O2", "H2O")
N_B_HILLIGARDT = 2.7
F_W_WAKE = 0.25


def _bed_cells(reactor: Reactor) -> list[int]:
    return [
        i
        for i, c in enumerate(reactor.cells)
        if str(getattr(c, "cell_type", "bed")) == "bed"
    ]


def _mole_frac(cell, phase: str = "d") -> dict[str, float]:
    idx = GAS_SPECIES_INDEX
    y = cell._mole_fractions(phase)
    return {sp: float(y[idx[sp]]) for sp in SPECIES_TRACE if sp in idx}


def _pyro_source_mol_s(cell) -> dict[str, float]:
    cache = getattr(cell, "_vm_gas_source_cache", None)
    if cache is None:
        return {}
    arr = np.asarray(cache, dtype=np.float64).ravel()
    out: dict[str, float] = {}
    for sp in GAS_SPECIES:
        if sp not in GAS_SPECIES_INDEX:
            continue
        v = float(arr[GAS_SPECIES_INDEX[sp]])
        if abs(v) > 1e-6:
            out[sp] = v
    return out


def _init_merit(reactor: Reactor) -> dict[str, float]:
    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    reactor._apply_all_bc_for_nr()
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    group = _residual_group_metrics(F_hat, cells)
    merit = _line_search_merit(
        F_hat,
        cells,
        use_gas_phase_split=bool(getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)),
    )
    return {
        "merit_scaled": float(merit),
        "gas_rms": float(group["gas_rms"]),
        "solid_rms": float(group["solid_rms"]),
        "energy_rms": float(group["energy_rms"]),
        "max_abs_scaled": float(np.max(np.abs(F_hat))) if F_hat.size else 0.0,
    }


def _bed0_formula_chain(reactor: Reactor, *, after_bc: bool = False) -> dict[str, Any]:
    """逐步还原 bed0 水力学 → K_auf → holdup seed 估算（无 cap，纯公式）。"""
    beds = _bed_cells(reactor)
    if not beds:
        return {}
    if after_bc:
        reactor._apply_all_bc_for_nr()
    c0 = reactor.cells[beds[0]]
    cfg = reactor.config
    eps_mf = float(getattr(cfg, "eps_mf", 0.45))
    u_d_cl = str(getattr(cfg, "hydrodynamics_u_d_closure", "wein_1992_eq312"))
    area = float(np.pi * c0.geo.D_bed**2 / 4.0)
    n_tot = float((c0.N_d + c0.N_b).sum())
    u0_from_n = n_tot * Rg * float(c0.T) / float(c0.P) / max(area, 1e-12)
    u_d_wein = 1.45 * float(c0.u_mf)
    if float(c0.u0) > u_d_wein:
        u_d_branch = "wein_potential"
        u_d_expected = u_d_wein
    else:
        u_d_branch = "micro_bubble_0.999"
        u_d_expected = float(c0.u_mf + (float(c0.u0) - float(c0.u_mf)) * 0.999)
    denom_eps = float(c0.u_b + (N_B_HILLIGARDT - 1.0) * c0.u_d)
    eps_num = float(c0.u0 - c0.u_d)
    eps_raw = eps_num / max(denom_eps, 1e-12)
    eps_clip = float(calc_visible_bubble_fraction(float(c0.u0), float(c0.u_b), float(c0.u_d)))
    dh = float(c0.geo.dh)
    kauf_formula = F_W_WAKE * float(c0.eps_b) * float(c0.u_b) / max((1.0 - float(c0.eps_b)) * dh, 1e-12)
    geo = float(cell_total_solid_holdup(c0))
    m_hold = float(c0.m_solid.sum())
    zu = float(c0.m_solid_zu.sum())
    k_total = float(c0.K_solid_auf[0, 0] + c0.K_solid_ab[0, 0])
    closure_alts: dict[str, dict[str, float]] = {}
    for cl in ("current", "hilligardt_eq311", "wein_1992_eq312", "backsolve_visible_epsb"):
        ud = float(
            _resolve_u_d_closure(
                strategy=cl,
                u0=float(c0.u0),
                u_mf=float(c0.u_mf),
                u_b=float(c0.u_b),
                eps_mf=eps_mf,
            )
        )
        closure_alts[cl] = {
            "u_d_m_s": ud,
            "eps_b": float(calc_visible_bubble_fraction(float(c0.u0), float(c0.u_b), ud)),
        }
    return {
        "after_bc": bool(after_bc),
        "feed_mol_s": {
            "O2": float(cfg.O2_feed),
            "H2O": float(cfg.H2O_feed),
            "N2": float(cfg.N2_feed),
        },
        "u0_from_N_m_s": float(u0_from_n),
        "u0_stored_m_s": float(c0.u0),
        "u_mf_m_s": float(c0.u_mf),
        "u_b_m_s": float(c0.u_b),
        "u_d_stored_m_s": float(c0.u_d),
        "u_d_closure_config": u_d_cl,
        "u_d_wein_branch": u_d_branch,
        "u_d_wein_expected_m_s": float(u_d_expected),
        "eps_b_numerator_m_s": eps_num,
        "eps_b_denominator_m_s": denom_eps,
        "eps_b_raw": float(eps_raw),
        "eps_b_stored": float(c0.eps_b),
        "eps_b_after_clip_0p01_0p7": float(np.clip(eps_clip, 0.01, 0.7)),
        "K_auf_formula_s_inv": float(kauf_formula),
        "K_auf_stored_s_inv": float(c0.K_solid_auf[0, 0]),
        "K_ab_stored_s_inv": float(c0.K_solid_ab[0, 0]),
        "K_total_s_inv": k_total,
        "geo_holdup_kg": geo,
        "m_solid_kg": m_hold,
        "holdup_over_geo": float(m_hold / max(geo, 1e-12)),
        "zu_fresh_kg_s": zu,
        "seed_inflow_over_K_kg": float(zu / max(k_total, 1e-12)),
        "N_gas_mol_s": n_tot,
        "closure_sensitivity": closure_alts,
    }


def _bed0_snapshot(reactor: Reactor) -> dict[str, Any]:
    beds = _bed_cells(reactor)
    if not beds:
        return {}
    i0 = beds[0]
    c0 = reactor.cells[i0]
    geo = cell_total_solid_holdup(c0)
    m_hold = float(c0.m_solid.sum())
    return {
        "bed0_index": i0,
        "T_K": float(c0.T),
        "m_solid_kg_s": m_hold,
        "geo_holdup_kg": float(geo),
        "holdup_over_geo": float(m_hold / max(geo, 1e-12)),
        "solid_out_kg_s": float(c0._solid_outflow_rates().sum()),
        "ab_in_kg_s": float(c0.m_solid_ab_in.sum()),
        "zu_fresh_kg_s": float(c0.m_solid_zu.sum()),
        "N_gas_mol": float((c0.N_d + c0.N_b).sum()),
        "y_dense": _mole_frac(c0, "d"),
        "energy_residual_W": float(c0.calc_energy_balance()),
        "energy_scaled": float(c0.calc_energy_balance()) / max(abs(float(reactor.config.fuel_feed)) * 20.0e6, 1.0e6),
        "pyro_source_mol_s": _pyro_source_mol_s(c0),
        "u_mf_m_s": float(c0.u_mf),
        "eps_b": float(c0.eps_b),
    }


def _stage_snapshot(reactor: Reactor, stage: str, *, compute_merit: bool) -> dict[str, Any]:
    out: dict[str, Any] = {
        "stage": stage,
        "bed0": _bed0_snapshot(reactor) if _bed_cells(reactor) else {},
    }
    if compute_merit:
        out["merit"] = _init_merit(reactor)
    return out


def _run_staged_init_incremental_only(reactor: Reactor) -> list[dict[str, Any]]:
    """逐步执行 init_precalc 五阶段（阶段 5 即最终态）。"""
    snapshots: list[dict[str, Any]] = []
    state = init_step._resolve_init_strategy_and_seed_T(reactor, init_strategy="vorabrechnung")
    snapshots.append(_stage_snapshot(reactor, STAGES[0], compute_merit=False))
    state = init_step._run_vorabrechnung_snapshot(reactor, state)
    snap2 = _stage_snapshot(reactor, STAGES[1], compute_merit=False)
    snap2["bed0_formula_chain"] = _bed0_formula_chain(reactor, after_bc=False)
    snap2["bed0_formula_chain_pre_seed"] = _bed0_formula_chain(reactor, after_bc=True)
    snapshots.append(snap2)
    init_step._run_post_vorabrechnung_holdup_and_reactive(reactor, state)
    snapshots.append(_stage_snapshot(reactor, STAGES[2], compute_merit=False))
    init_step._run_bottom_zone_preprojections(reactor, state)
    snapshots.append(_stage_snapshot(reactor, STAGES[3], compute_merit=True))
    init_step._run_init_abgleich_and_joint_preprojection(reactor, state)
    snapshots.append(_stage_snapshot(reactor, STAGES[4], compute_merit=True))
    snapshots.append(_stage_snapshot(reactor, STAGES[5], compute_merit=True))
    return snapshots


def _delta_table(w1: list[dict], w2: list[dict]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s1, s2 in zip(w1, w2):
        stage = s1["stage"]
        m1, m2 = s1.get("merit") or {}, s2.get("merit") or {}
        b1, b2 = s1["bed0"], s2["bed0"]
        rows.append(
            {
                "stage": stage,
                "merit_W1": m1.get("merit_scaled"),
                "merit_W2": m2.get("merit_scaled"),
                "delta_merit": (
                    (m2["merit_scaled"] - m1["merit_scaled"])
                    if m1.get("merit_scaled") is not None and m2.get("merit_scaled") is not None
                    else None
                ),
                "energy_rms_W1": m1.get("energy_rms"),
                "energy_rms_W2": m2.get("energy_rms"),
                "bed0_m_solid_W1": b1.get("m_solid_kg_s"),
                "bed0_m_solid_W2": b2.get("m_solid_kg_s"),
                "bed0_holdup_over_geo_W2": b2.get("holdup_over_geo"),
                "bed0_E_scaled_W1": b1.get("energy_scaled"),
                "bed0_E_scaled_W2": b2.get("energy_scaled"),
            }
        )
    return rows


def audit_case(case_key: str) -> dict[str, Any]:
    reactor = build_phase2_reactor(
        validation_key=case_key,
        ladder_level="p2",
        apply_kinetics=True,
    )
    cfg = reactor.config
    stages = _run_staged_init_incremental_only(reactor)
    return {
        "case_key": case_key,
        "primary_agent": str(getattr(cfg, "primary_agent", "")),
        "ER": float(cfg.ER or 0.0),
        "fuel_kg_s": float(cfg.fuel_feed),
        "O2_feed_mol_s": float(cfg.O2_feed),
        "H2O_feed_mol_s": float(cfg.H2O_feed),
        "N2_feed_mol_s": float(cfg.N2_feed),
        "steam_to_o2_molar": float(cfg.steam_to_o2_molar),
        "stages": stages,
        "final_merit": stages[-1].get("merit"),
        "final_bed0": stages[-1]["bed0"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        nargs="*",
        default=["CASE_HTW_WESSELING_1", "CASE_HTW_WESSELING_2"],
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=data_path("init_staged_audit_w1_w2.json"),
    )
    args = parser.parse_args()

    results: dict[str, Any] = {}
    for key in args.cases:
        print(f"Auditing {key} ...")
        results[key] = audit_case(key)

    if "CASE_HTW_WESSELING_1" in results and "CASE_HTW_WESSELING_2" in results:
        results["w1_vs_w2_delta"] = _delta_table(
            results["CASE_HTW_WESSELING_1"]["stages"],
            results["CASE_HTW_WESSELING_2"]["stages"],
        )

    payload = {
        "_meta": {
            "description": "Init precalc 分阶段审计 — bed0 与 NR merit",
            "stages": list(STAGES),
            "init_pipeline": "src/workflow/steps/init_precalc_step.py",
        },
        "cases": results,
    }
    write_json(args.output, payload)

    if "CASE_HTW_WESSELING_2" in results:
        fc = results["CASE_HTW_WESSELING_2"]["stages"][1].get("bed0_formula_chain_pre_seed") or {}
        if fc:
            print("\n=== W2 bed0 公式链（阶段2末，BC 后、seed 前）===")
            print(
                f"  u0={fc['u0_stored_m_s']:.4f}  u_d(Wein)={fc['u_d_wein_expected_m_s']:.4f}  "
                f"Δu={fc['eps_b_numerator_m_s']:.4f} → eps_b={fc['eps_b_stored']:.4f}"
            )
            print(
                f"  K_auf={fc['K_auf_stored_s_inv']:.4f}  zu={fc['zu_fresh_kg_s']:.3f}  "
                f"seed≈{fc['seed_inflow_over_K_kg']:.1f} kg  (m_solid={fc['m_solid_kg']:.3f})"
            )

    if "w1_vs_w2_delta" in results:
        print("\n=== W1 vs W2 分阶段 merit / bed0 holdup ===")
        for row in results["w1_vs_w2_delta"]:
            mw1, mw2, dm = row["merit_W1"], row["merit_W2"], row["delta_merit"]
            if mw1 is not None and mw2 is not None and dm is not None:
                merit_str = f"merit W1={mw1:.4f}  W2={mw2:.4f}  Δ={dm:+.4f}"
            else:
                merit_str = "merit n/a"
            print(
                f"  {row['stage']:36s}  {merit_str}  |  "
                f"m_solid W2={row['bed0_m_solid_W2']:.1f}  hold/geo={row['bed0_holdup_over_geo_W2']:.2f}"
            )
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
