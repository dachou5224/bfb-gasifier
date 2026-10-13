#!/usr/bin/env python3
"""Phase2 LU：Char 还原（R2/R4）参数审计与 r2/r4 灵敏度。

优先级
------
1. R4（Boudouard）动力学常数与 intrinsic rate 审计（Hamel §5.1.4 / Weeda）
2. 收敛态床层 R2 vs R4 积分预算 + r4/r2 标定扫描
3. 气相/固相停留时间 exposure 对比（τ_gas vs τ_char proxy）

Usage
-----
    python3 scripts/audit_char_reduction_sensitivity_lu.py --phase all
    python3 scripts/audit_char_reduction_sensitivity_lu.py --phase sensitivity --r4-grid 0.5,1,5,10,20,50
"""

from __future__ import annotations

from _lu_harness import (
    add_max_global_iter_arg,
    build_phase2_config,
    data_path,
    load_case_lu,
    phase2_solve_kwargs,
    repo_root,
    write_json,
)

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

import numpy as np

REPO_ROOT = repo_root()

from src.core.cell import S_CHAR
from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx, gas_diffusivity_correlation
from src.kinetics import char_reactions as cr
from src.kinetics.char_reactions import (
    R2_E_Rg,
    R2_k0,
    R4_kb_E,
    R4_kb_k0,
    R4_kc_E,
    R4_kc_k0,
    R4_kf_E,
    R4_kf_k0,
    d_core_from_spm_char_conversion,
    rate_R2,
    rate_R4,
    rate_R4_effective,
)
# 复用床层 CO/CO2 审计
from audit_bed_co_co2_budget_lu import _analyze_cell_co_co2, _sum_terms

REF = {"CO": 0.13, "CO2": 0.11}


@dataclass
class SolveRow:
    r2_scale: float
    r4_scale: float
    r5_scale: float
    converged_fully: bool
    rms_scaled_final: float
    co: float
    co2: float
    h2: float
    bed_co: float
    bed_co2: float
    bed_R2_CO: float
    bed_R4_CO: float
    bed_R5_CO: float
    score_co_co2: float
    wall_time_s: float


def _score(co: float, co2: float) -> float:
    return abs(co - REF["CO"]) / REF["CO"] + abs(co2 - REF["CO2"]) / REF["CO2"]


def _intrinsic_audit(*, P: float, d_p: float) -> list[dict]:
    """Hamel/Weeda R2/R4 intrinsic rate 表（与 O2 无关的缺氧区代表工况）。"""
    rows: list[dict] = []
    y_h2o, y_co2, y_co = 0.23, 0.14, 0.03
    c_h2o = y_h2o * P / (Rg * 1100.0)
    p_co2 = y_co2 * P
    p_co = y_co * P
    D_g = gas_diffusivity_correlation(1100.0, P)
    for T in (1000.0, 1050.0, 1100.0, 1150.0):
        D_g_T = gas_diffusivity_correlation(T, P)
        c_h2o_T = y_h2o * P / (Rg * T)
        p_co2_T = y_co2 * P
        p_co_T = y_co * P
        for X in (0.3, 0.6, 0.9):
            dc = d_core_from_spm_char_conversion(X, d_p)
            r4_kin = rate_R4(T, p_co2_T, p_co_T)
            r4_eff = rate_R4_effective(T, p_co2_T, p_co_T, d_p, D_g_T, dc)
            r2_eff = rate_R2(T, c_h2o_T, d_p, D_g_T, dc)
            rows.append(
                {
                    "T_K": T,
                    "char_conversion": X,
                    "d_core_mm": float(dc * 1e3),
                    "y_H2O": y_h2o,
                    "y_CO2": y_co2,
                    "y_CO": y_co,
                    "r4_kin_mol_m2_s": float(r4_kin),
                    "r4_eff_mol_m2_s": float(r4_eff),
                    "r2_eff_mol_m2_s": float(r2_eff),
                    "r2_over_r4_eff": float(r2_eff / max(r4_eff, 1e-30)),
                    "area_7p3m2_R4_CO_mol_s": float(2.0 * r4_eff * 7.3),
                    "area_7p3m2_R2_CO_mol_s": float(r2_eff * 7.3),
                }
            )
    return rows


def _kinetics_constants_snapshot() -> dict:
    return {
        "R2": {"k0": R2_k0, "Ea_J_mol": R2_E_Rg * Rg, "ref": "Hamel (1999) Eq.5.24-5.25"},
        "R4_forward": {"k0": R4_kf_k0, "Ea_J_mol": R4_kf_E, "ref": "Hamel (1999) Eq.5.26 Weeda"},
        "R4_kb": {"k0": R4_kb_k0, "Ea_J_mol": R4_kb_E},
        "R4_kc": {"k0": R4_kc_k0, "Ea_J_mol": R4_kc_E},
        "D_A_ref_m2_s": float(cr.D_A_ref),
        "note": "R2/R3/R4 均按 Hamel §5.1.3–5.1.5 表面动力学；R4_effective 与 rate_R4 相同 (Eq.5.33 忽略粒径)",
    }


def _char_exposure_rows(reactor: Reactor) -> list[dict]:
    cfg = reactor.config
    rows = []
    for i, cell in enumerate(reactor.cells):
        cell.calc_hydrodynamics()
        m_char = float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0)))
        w_up = float(np.sum(np.maximum(cell._solid_upflow_rates()[:, S_CHAR], 0.0)))
        w_down = float(np.sum(np.maximum(cell._solid_downflow_rates()[:, S_CHAR], 0.0)))
        w_out = max(w_up + w_down, 1e-12)
        area = float(np.sum(cell._calc_char_surface_area_per_class()))
        x_char = float(cell._compute_char_conversion())
        tau_gas_s = float(cell.geo.dh / max(cell.u_mf, 1e-3))
        tau_char_s = float(m_char / w_out)
        y = cell._mole_fractions("combined")
        rows.append(
            {
                "cell_index": i,
                "xi_over_H": float(cell.geo.h_center / max(float(cfg.H_bed), 1e-12)),
                "T_K": float(cell.T),
                "y_O2": float(y[idx["O2"]]),
                "y_CO2": float(y[idx["CO2"]]),
                "y_CO": float(y[idx["CO"]]),
                "y_H2O": float(y[idx["H2O"]]),
                "m_char_kg": m_char,
                "char_area_m2": area,
                "char_conversion": x_char,
                "char_upflow_kg_s": w_up,
                "char_downflow_kg_s": w_down,
                "tau_gas_s": tau_gas_s,
                "tau_char_proxy_s": tau_char_s,
                "tau_char_over_tau_gas": float(tau_char_s / max(tau_gas_s, 1e-9)),
            }
        )
    return rows


def _bed_reduction_budget(reactor: Reactor) -> dict:
    cfg = reactor.config
    cell_rows = [_analyze_cell_co_co2(c) for c in reactor.cells]
    for i, row in enumerate(cell_rows):
        row["cell_index"] = i
        row["xi_over_H"] = float(reactor.cells[i].geo.h_center / max(float(cfg.H_bed), 1e-12))
    co_sum = _sum_terms(cell_rows, "co_terms")
    co2_sum = _sum_terms(cell_rows, "co2_terms")
    anoxic = [r for r in cell_rows if float(r["wet_y"]["O2"]) < 1e-4]
    anoxic_co = _sum_terms(anoxic, "co_terms") if anoxic else {}
    return {
        "bed_integrated_co": co_sum,
        "bed_integrated_co2": co2_sum,
        "anoxic_cells": len(anoxic),
        "anoxic_integrated_co": anoxic_co,
        "cells": cell_rows,
    }


def _run_case_fixed(
    *,
    r2_scale: float,
    r4_scale: float,
    r5_scale: float,
    max_global_iter: int,
) -> tuple[SolveRow, Reactor, dict]:
    cfg = build_phase2_config()
    cfg.vorab_major_gibbs_x0 = False
    cfg.r2_scale = float(max(r2_scale, 0.0))
    cfg.r4_scale = float(max(r4_scale, 0.0))
    cfg.r5_scale = float(max(r5_scale, 0.0))
    kw = phase2_solve_kwargs()
    kw["max_global_iter"] = int(max_global_iter)
    kw["verbose"] = False
    t0 = perf_counter()
    reactor = Reactor(cfg)
    out = reactor.solve(**kw)
    wall = perf_counter() - t0
    budget = _bed_reduction_budget(reactor)
    co_sum = budget["bed_integrated_co"]
    y = out.get("exit_gas_dry") or {}
    bed = out.get("bed_exit_gas_dry") or {}
    r5_co = float(co_sum.get("R5b", 0.0) + co_sum.get("R5d", 0.0))
    row = SolveRow(
        r2_scale=float(r2_scale),
        r4_scale=float(r4_scale),
        r5_scale=float(r5_scale),
        converged_fully=bool(out.get("converged_fully")),
        rms_scaled_final=float(out.get("rms_scaled_final", 0.0) or 0.0),
        co=float(y.get("CO", 0.0)),
        co2=float(y.get("CO2", 0.0)),
        h2=float(y.get("H2", 0.0)),
        bed_co=float(bed.get("CO", 0.0)),
        bed_co2=float(bed.get("CO2", 0.0)),
        bed_R2_CO=float(co_sum.get("R2", 0.0)),
        bed_R4_CO=float(co_sum.get("R4", 0.0)),
        bed_R5_CO=r5_co,
        score_co_co2=_score(float(y.get("CO", 0.0)), float(y.get("CO2", 0.0))),
        wall_time_s=float(wall),
    )
    return row, reactor, budget


def _parse_grid(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("intrinsic", "exposure", "sensitivity", "all"), default="all")
    parser.add_argument("--max-global-iter", type=int, default=25)
    parser.add_argument("--r5-scale", type=float, default=0.75)
    parser.add_argument("--r4-grid", type=str, default="0.5,1,2,5,10,20,50,100")
    parser.add_argument("--r2-grid", type=str, default="0.5,1,2,5,10")
    parser.add_argument("--baseline-r4", type=float, default=20.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_char_reduction_audit.json",
    )
    args = parser.parse_args()

    case = load_case_lu()
    P = float(case["P"])
    d_p = 0.5 * sum(case.get("particle_diameter_mm_range", (1.5, 3.0))) * 1.0e-3

    payload: dict = {"case_key": "CASE_HTW_WESSELING_1", "ref_dry": REF}

    if args.phase in ("intrinsic", "all"):
        payload["kinetics_constants"] = _kinetics_constants_snapshot()
        payload["intrinsic_rates"] = _intrinsic_audit(P=P, d_p=d_p)
        print("=" * 100)
        print("R2/R4 intrinsic rate audit (anoxic-like gas, per m² char)")
        print("=" * 100)
        print(f"{'T':>5} {'X':>4} {'r4_eff':>12} {'r2_eff':>12} {'R2/R4':>10} {'7.3m² R4→CO':>12} {'7.3m² R2→CO':>12}")
        for row in payload["intrinsic_rates"]:
            if row["T_K"] == 1100.0 and row["char_conversion"] in (0.3, 0.6, 0.9):
                print(
                    f"{row['T_K']:>5.0f} {row['char_conversion']:>4.1f} "
                    f"{row['r4_eff_mol_m2_s']:>12.3e} {row['r2_eff_mol_m2_s']:>12.3e} "
                    f"{row['r2_over_r4_eff']:>10.0f} {row['area_7p3m2_R4_CO_mol_s']:>12.4f} "
                    f"{row['area_7p3m2_R2_CO_mol_s']:>12.4f}"
                )

    exposure_reactor: Reactor | None = None
    exposure_budget: dict | None = None

    if args.phase in ("sensitivity", "all"):
        r4_grid = _parse_grid(args.r4_grid)
        sens_rows: list[dict] = []
        print("=" * 100)
        print(f"r4_scale sensitivity (r2=1, r5={args.r5_scale})")
        print("=" * 100)
        best_r4 = r4_grid[0]
        best_score = float("inf")
        for i, r4 in enumerate(r4_grid, start=1):
            row, reactor, budget = _run_case_fixed(
                r2_scale=1.0,
                r4_scale=r4,
                r5_scale=args.r5_scale,
                max_global_iter=args.max_global_iter,
            )
            sens_rows.append(asdict(row))
            print(
                f"[r4 {i}/{len(r4_grid)}] r4={r4:.1f} conv={row.converged_fully} "
                f"CO={row.co:.4f} CO2={row.co2:.4f} bed_R2={row.bed_R2_CO:.2f} bed_R4={row.bed_R4_CO:.3f} "
                f"score={row.score_co_co2:.3f} ({row.wall_time_s:.0f}s)"
            )
            if row.converged_fully and row.score_co_co2 < best_score:
                best_score = row.score_co_co2
                best_r4 = r4
                exposure_reactor = reactor
                exposure_budget = budget

        r2_grid = _parse_grid(args.r2_grid)
        r4_for_r2 = float(args.baseline_r4) if args.baseline_r4 > 0 else best_r4
        print("=" * 100)
        print(f"r2_scale sensitivity (r4={r4_for_r2}, r5={args.r5_scale})")
        print("=" * 100)
        for i, r2 in enumerate(r2_grid, start=1):
            row, reactor, budget = _run_case_fixed(
                r2_scale=r2,
                r4_scale=r4_for_r2,
                r5_scale=args.r5_scale,
                max_global_iter=args.max_global_iter,
            )
            sens_rows.append(asdict(row))
            print(
                f"[r2 {i}/{len(r2_grid)}] r2={r2:.1f} conv={row.converged_fully} "
                f"CO={row.co:.4f} CO2={row.co2:.4f} bed_R2={row.bed_R2_CO:.2f} bed_R4={row.bed_R4_CO:.3f} "
                f"score={row.score_co_co2:.3f} ({row.wall_time_s:.0f}s)"
            )
            if row.converged_fully and row.score_co_co2 < best_score:
                best_score = row.score_co_co2
                exposure_reactor = reactor
                exposure_budget = budget

        payload["sensitivity"] = {
            "r5_scale": float(args.r5_scale),
            "r4_grid": r4_grid,
            "r2_grid": r2_grid,
            "best_r4_from_r4_sweep": float(best_r4),
            "rows": sens_rows,
        }
        converged = [r for r in sens_rows if r["converged_fully"]]
        if converged:
            best = min(converged, key=lambda r: r["score_co_co2"])
            payload["best_converged"] = best

    if args.phase in ("exposure", "all"):
        if exposure_reactor is None:
            _, exposure_reactor, exposure_budget = _run_case_fixed(
                r2_scale=1.0,
                r4_scale=0.5,
                r5_scale=args.r5_scale,
                max_global_iter=args.max_global_iter,
            )
        payload["char_exposure"] = _char_exposure_rows(exposure_reactor)
        if exposure_budget is None:
            exposure_budget = _bed_reduction_budget(exposure_reactor)
        payload["bed_reduction_budget"] = {
            "bed_integrated_co": exposure_budget["bed_integrated_co"],
            "bed_integrated_co2": exposure_budget["bed_integrated_co2"],
            "anoxic_integrated_co": exposure_budget.get("anoxic_integrated_co", {}),
        }
        print("=" * 100)
        print("Char vs gas exposure (baseline or best sensitivity reactor)")
        print("=" * 100)
        print(f"{'cell':>4} {'xi/H':>5} {'yO2':>7} {'tau_g':>6} {'tau_c':>8} {'tc/tg':>7} {'R2co':>8} {'R4co':>8}")
        cells = exposure_budget.get("cells") or []
        exp = payload["char_exposure"]
        for crow, erow in zip(cells, exp):
            co_t = crow["co_terms"]
            print(
                f"{erow['cell_index']:>4d} {erow['xi_over_H']:>5.2f} {erow['y_O2']:>7.1e} "
                f"{erow['tau_gas_s']:>6.3f} {erow['tau_char_proxy_s']:>8.1f} {erow['tau_char_over_tau_gas']:>7.0f} "
                f"{co_t.get('R2', 0.0):>8.3f} {co_t.get('R4', 0.0):>8.4f}"
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
