#!/usr/bin/env python3
"""Phase2 LU 床层 CO/CO2 来源–去向审计。

逐 cell 按 ``build_reaction_sources`` 同口径重算 CO/CO2 各反应净贡献，
定位 bed-exit 组成偏差的主要驱动（R5 氧化 / R1 燃尽 / R4 Boudouard / R8 WGSR / 热解）。

Usage
-----
    python3 scripts/audit_bed_co_co2_budget_lu.py
    python3 scripts/audit_bed_co_co2_budget_lu.py --output data/phase2_bed_co_co2_audit.json
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
from _audit_common import compute_char_reaction_rates, prepare_cell_for_analysis

REPO_ROOT = repo_root()

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx, gas_diffusivity_correlation
from src.kinetics.gas_reactions import (
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R7,
    rate_R8,
)
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
REF_DRY = {"CO": 0.13, "CO2": 0.11}


def _co_co2_term_dict() -> dict[str, float]:
    return {
        "pyrolysis": 0.0,
        "R1": 0.0,
        "R2": 0.0,
        "R4": 0.0,
        "R5b": 0.0,
        "R5d": 0.0,
        "R7": 0.0,
        "R8": 0.0,
        "R10b": 0.0,
        "R10d": 0.0,
        "R11b": 0.0,
        "R11d": 0.0,
    }


def _analyze_cell_co_co2(cell) -> dict[str, object]:
    """按 cell_kinetics limiter 口径分解 CO/CO2 净源项 [mol/s]。"""
    prepare_cell_for_analysis(cell)

    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    gas_src_vm = cell._vm_gas_source_cache.copy()
    rho_cat = cell._catalyst_bulk_density()
    D_g = gas_diffusivity_correlation(T, P)
    n_ex = cell.N_ex

    r5b = cell.r5_scale * rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], C_b[idx["H2O"]], P, y_b, idx) * cell.V_b
    r5d = cell.r5_scale * rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * cell.V_d
    ext7 = cell.r7_scale * rate_R7(
        T, C_d[idx["CH4"]], C_d[idx["H2O"]], C_d[idx["CO"]], C_d[idx["H2"]],
        rho_ash_kg_m3=cell._ash_bulk_density(),
    ) * cell.V_d
    ext8 = cell.r8_scale * rate_R8(
        T,
        P,
        float(y_d[idx["CO"]]),
        float(y_d[idx["H2O"]]),
        float(y_d[idx["CO2"]]),
        float(y_d[idx["H2"]]),
        rho_ash_kg_m3=cell._ash_bulk_density(),
    ) * cell.V_d

    C_tar_b = max(float(C_b[idx["TAR1"]] + C_b[idx["TAR2"]]), 0.0)
    C_tar_d = max(float(C_d[idx["TAR1"]] + C_d[idx["TAR2"]]), 0.0)
    ext10b = rate_R10(T, C_tar_b, C_b[idx["O2"]], P, cell.fuel_type) * cell.V_b
    ext10d = rate_R10(T, C_tar_d, C_d[idx["O2"]], P, cell.fuel_type) * cell.V_d
    ext11b = rate_R11_bubble(T, C_tar_b) * cell.V_b
    ext11d = rate_R11_suspension(T, C_tar_d, rho_cat) * cell.V_d
    stoich_r10 = get_lumped_tar_stoichiometry("R10", cell.fuel_type)
    stoich_r11 = get_lumped_tar_stoichiometry("R11", cell.fuel_type)
    nu_o2_r10 = max(0.0, -float(stoich_r10.get("O2", 0.0)))

    char_rates = compute_char_reaction_rates(cell, D_g)
    r1 = char_rates["r1"]
    r2 = char_rates["r2"]
    r3 = char_rates["r3"]
    r4 = char_rates["r4"]
    alpha = char_rates["alpha"]
    total_area = char_rates["total_area"]

    o2_supply_bubble = max(
        float(cell.N_zu_b[idx["O2"]]) + float(cell.N_b_in[idx["O2"]])
        + max(-float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_b[idx["O2"]]),
        1e-6,
    )
    o2_supply_dense = max(
        float(cell.N_zu_d[idx["O2"]]) + float(cell.N_d_in[idx["O2"]])
        + max(float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_d[idx["O2"]]),
        1e-6,
    )
    # 工程对照镜像（非生产 SoT）：生产默认 extent_limiters_enabled_thesis=False。
    # R5: CO + 1/2 O2 → CO2；extent = mol CO/s → O₂ demand = 0.5·r5
    bubble_o2_demand = 0.5 * r5b + nu_o2_r10 * ext10b
    dense_o2_demand = 0.5 * r5d + alpha * r1 * total_area + nu_o2_r10 * ext10d
    limit_o2_b = min(1.0, 0.995 * o2_supply_bubble / max(bubble_o2_demand, 1e-9))
    limit_o2_d = min(1.0, 0.995 * o2_supply_dense / max(dense_o2_demand, 1e-9))

    co_supply_b = max(
        float(cell.N_zu_b[idx["CO"]]) + float(cell.N_b_in[idx["CO"]])
        + max(-float(n_ex[idx["CO"]]), 0.0) + float(cell.N_rez_b[idx["CO"]]),
        1e-6,
    )
    limit_co_b = min(1.0, 0.995 * co_supply_b / max(r5b, 1e-9))

    r5b_lim = r5b * min(limit_o2_b, limit_co_b)
    r5d_lim = r5d * limit_o2_d

    h2o_supply = max(
        float(cell.N_zu_d[idx["H2O"]]) + float(cell.N_d_in[idx["H2O"]])
        + float(cell.N_zu_b[idx["H2O"]]) + float(cell.N_b_in[idx["H2O"]])
        + float(cell.N_rez_d[idx["H2O"]]) + float(cell.N_rez_b[idx["H2O"]])
        + float(gas_src_vm[idx["H2O"]]),
        1e-6,
    )
    h2o_consume = r2 * total_area + ext7 + max(0.0, ext8)
    limit_h2o = min(1.0, 0.995 * h2o_supply / max(h2o_consume, 1e-9))
    r2_lim = r2 * limit_h2o
    ext7_lim = ext7 * limit_h2o
    ext8_lim = ext8 * limit_h2o if ext8 > 0.0 else ext8
    ext10b_lim = ext10b * limit_o2_b
    ext10d_lim = ext10d * limit_o2_d
    ext11b_lim = ext11b * limit_h2o
    ext11d_lim = ext11d * limit_h2o
    r1_lim = r1 * limit_o2_d
    r4_lim = r4

    co_terms = _co_co2_term_dict()
    co2_terms = _co_co2_term_dict()

    co_terms["pyrolysis"] += float(gas_src_vm[idx["CO"]])
    co2_terms["pyrolysis"] += float(gas_src_vm[idx["CO2"]])

    co_terms["R5b"] -= r5b_lim
    co2_terms["R5b"] += r5b_lim
    co_terms["R5d"] -= r5d_lim
    co2_terms["R5d"] += r5d_lim

    co_terms["R7"] += ext7_lim

    co_terms["R8"] -= ext8_lim
    co2_terms["R8"] += ext8_lim

    if total_area > 0.0:
        r1_ext = r1_lim * total_area
        r2_ext = r2_lim * total_area
        r4_ext = r4_lim * total_area
        r1_o2_ext = alpha * r1_ext
        co_terms["R1"] += 2.0 * (r1_ext - r1_o2_ext)
        co2_terms["R1"] += 2.0 * r1_o2_ext - r1_ext
        co_terms["R2"] += r2_ext
        co_terms["R4"] += 2.0 * r4_ext
        co2_terms["R4"] -= r4_ext

    for sp_key, nu in stoich_r10.items():
        if sp_key == "CO":
            co_terms["R10b"] += ext10b_lim * float(nu)
            co_terms["R10d"] += ext10d_lim * float(nu)
        elif sp_key == "CO2":
            co2_terms["R10b"] += ext10b_lim * float(nu)
            co2_terms["R10d"] += ext10d_lim * float(nu)
    for sp_key, nu in stoich_r11.items():
        if sp_key == "CO":
            co_terms["R11b"] += ext11b_lim * float(nu)
            co_terms["R11d"] += ext11d_lim * float(nu)
        elif sp_key == "CO2":
            co2_terms["R11b"] += ext11b_lim * float(nu)
            co2_terms["R11d"] += ext11d_lim * float(nu)

    wet = cell._mole_fractions("combined")
    return {
        "xi_over_H": float(cell.geo.h_center / max(cell.geo.h_center + cell.geo.dh, 1e-9)),
        "T_K": float(T),
        "wet_y": {sp: float(wet[idx[sp]]) for sp in ("CO", "CO2", "H2", "O2", "H2O")},
        "co_terms": co_terms,
        "co2_terms": co2_terms,
        "co_net": float(sum(co_terms.values())),
        "co2_net": float(sum(co2_terms.values())),
        "raw": {
            "r5b": float(r5b),
            "r5b_lim": float(r5b_lim),
            "r5d": float(r5d),
            "r5d_lim": float(r5d_lim),
            "ext8": float(ext8),
            "ext8_lim": float(ext8_lim),
            "r1_area": float(r1 * total_area),
            "r4_area": float(r4 * total_area),
            "gas_src_co": float(gas_src_vm[idx["CO"]]),
            "gas_src_co2": float(gas_src_vm[idx["CO2"]]),
        },
        "limiters": {
            "o2_bubble": float(limit_o2_b),
            "o2_dense": float(limit_o2_d),
            "co_bubble": float(limit_co_b),
            "h2o": float(limit_h2o),
        },
        "char_area_m2": float(total_area),
    }


def _sum_terms(rows: list[dict], key: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in rows:
        for term, val in row[key].items():
            out[term] = out.get(term, 0.0) + float(val)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_bed_co_co2_audit.json",
    )
    parser.add_argument("--max-global-iter", type=int, default=None, help="覆盖 PHASE2 默认 max_global_iter（默认 150）")
    parser.add_argument(
        "--vm-zone-cells",
        type=int,
        default=None,
        help="vorab_vm_devolatilization_zone_cells_thesis（默认 Phase2 配置值）",
    )
    parser.add_argument("--r2-scale", type=float, default=None)
    parser.add_argument("--r4-scale", type=float, default=None)
    parser.add_argument("--r5-scale", type=float, default=None)
    parser.add_argument("--r8-scale", type=float, default=None)
    args = parser.parse_args()

    case = load_case_lu()
    cfg = build_phase2_config(case)
    cfg.vorab_major_gibbs_x0 = False
    if args.vm_zone_cells is not None:
        cfg.vorab_vm_devolatilization_zone_cells_thesis = int(args.vm_zone_cells)
    if args.r2_scale is not None:
        cfg.r2_scale = float(args.r2_scale)
    if args.r4_scale is not None:
        cfg.r4_scale = float(args.r4_scale)
    if args.r5_scale is not None:
        cfg.r5_scale = float(args.r5_scale)
    if args.r8_scale is not None:
        cfg.r8_scale = float(args.r8_scale)
    reactor = Reactor(cfg)
    kw = phase2_solve_kwargs()
    if args.max_global_iter is not None:
        kw["max_global_iter"] = int(args.max_global_iter)
    kw["verbose"] = False
    result, monitor = solve_with_nr_monitor(reactor, kw, check_x0=False)

    rows = [_analyze_cell_co_co2(cell) for cell in reactor.cells]
    H_bed = float(cfg.H_bed)
    for i, row in enumerate(rows):
        row["cell_index"] = i
        row["xi_over_H"] = float(reactor.cells[i].geo.h_center / H_bed)

    co_sum = _sum_terms(rows, "co_terms")
    co2_sum = _sum_terms(rows, "co2_terms")
    bed_exit = result.get("bed_exit_gas_dry") or {}
    exit_dry = result.get("exit_gas_dry") or {}

    payload = {
        "case_key": "CASE_HTW_WESSELING_1",
        "vm_zone_cells": int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0)),
        "kinetics": {
            "r2_scale": cfg.r2_scale,
            "r4_scale": cfg.r4_scale,
            "r5_scale": cfg.r5_scale,
            "r8_scale": cfg.r8_scale,
        },
        "converged_fully": bool(result.get("converged_fully")),
        "rms_scaled_final": float(result.get("rms_scaled_final", 0.0) or 0.0),
        "ref_dry": REF_DRY,
        "bed_exit_gas_dry": {k: float(bed_exit.get(k, 0.0)) for k in REF_DRY},
        "reactor_exit_gas_dry": {k: float(exit_dry.get(k, 0.0)) for k in REF_DRY},
        "bed_integrated_co_terms_mol_s": co_sum,
        "bed_integrated_co2_terms_mol_s": co2_sum,
        "cells": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=" * 140)
    print("Phase2 LU bed CO/CO2 budget audit")
    print(
        f"kinetics: r2={cfg.r2_scale:.2f} r4={cfg.r4_scale:.2f} "
        f"r5={cfg.r5_scale:.2f} r8={cfg.r8_scale:.2f}"
    )
    print("=" * 140)
    print_nr_monitor(monitor, prefix="NR")
    print(
        f"conv={payload['converged_fully']} rms={payload['rms_scaled_final']:.4e} "
        f"bed_exit_dry CO={bed_exit.get('CO', 0):.4f} CO2={bed_exit.get('CO2', 0):.4f} "
        f"(ref {REF_DRY['CO']}/{REF_DRY['CO2']})"
    )
    print("-" * 140)
    print("Bed-integrated net terms [mol/s] (CO sink negative, CO2 source positive):")
    print(f"  {'term':>8} {'CO':>12} {'CO2':>12}")
    for term in sorted(set(co_sum) | set(co2_sum)):
        print(f"  {term:>8} {co_sum.get(term, 0.0):>12.4f} {co2_sum.get(term, 0.0):>12.4f}")
    print(f"  {'TOTAL':>8} {sum(co_sum.values()):>12.4f} {sum(co2_sum.values()):>12.4f}")

    print("-" * 140)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>7} {'yCO':>7} {'yCO2':>7} {'yO2':>7} "
        f"{'R5_net':>9} {'R1_CO':>9} {'R4_CO':>9} {'R8_CO':>9} {'pyr_CO':>9} "
        f"{'R5_CO2':>9} {'R1_CO2':>9} {'R8_CO2':>9} {'pyr_CO2':>9} {'co_net':>9}"
    )
    print("-" * 140)
    for row in rows:
        co, co2 = row["co_terms"], row["co2_terms"]
        r5_co = co["R5b"] + co["R5d"]
        r5_co2 = co2["R5b"] + co2["R5d"]
        print(
            f"{row['cell_index']:>4d} {row['xi_over_H']:>6.2f} {row['T_K']:>7.1f} "
            f"{row['wet_y']['CO']:>7.4f} {row['wet_y']['CO2']:>7.4f} {row['wet_y']['O2']:>7.4f} "
            f"{r5_co:>9.4f} {co['R1']:>9.4f} {co['R4']:>9.4f} {co['R8']:>9.4f} {co['pyrolysis']:>9.4f} "
            f"{r5_co2:>9.4f} {co2['R1']:>9.4f} {co2['R8']:>9.4f} {co2['pyrolysis']:>9.4f} {row['co_net']:>9.4f}"
        )

    print("-" * 140)
    top3 = sorted(rows, key=lambda r: abs(r["co_terms"]["R5b"] + r["co_terms"]["R5d"]))[-3:]
    print("Top-3 |R5| cells (CO oxidation hotspots):")
    for row in reversed(top3):
        i = row["cell_index"]
        r5 = row["co_terms"]["R5b"] + row["co_terms"]["R5d"]
        print(
            f"  cell {i} xi/H={row['xi_over_H']:.2f} T={row['T_K']:.0f}K "
            f"R5_co={r5:.4f} r5b_lim/raw={row['raw']['r5b_lim']:.4f}/{row['raw']['r5b']:.4f} "
            f"O2y={row['wet_y']['O2']:.4f} limit_o2_b={row['limiters']['o2_bubble']:.3f}"
        )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
