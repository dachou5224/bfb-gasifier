#!/usr/bin/env python3
"""Phase2 LU：均相区（气泡 vs 悬浮相）还原/氧化反应效应审计。

按 ``build_reaction_sources`` 同口径，逐 cell 拆分：
- 气泡相均相：R5b / R6b / R10b / R11b / R12b
- 悬浮相均相：R5d / R6d / R7 / R8 / R10d / R11d / R12d
- 异相 R1–R4 仅进入悬浮相

重点回答：气泡相是否存在 WGSR/重整等还原路径；各相对 CO/CO2 净贡献。

Usage
-----
    python3 scripts/audit_homogeneous_bubble_reduction_lu.py
"""

from __future__ import annotations

from _lu_harness import (
    add_max_global_iter_arg,
    build_phase2_config,
    data_path,
    phase2_solve_kwargs,
    repo_root,
    write_json,
)

REPO_ROOT = repo_root()

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10, rate_R11_bubble, rate_R11_suspension
from src.kinetics.gas_reactions import (
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R6,
    rate_R7,
    rate_R8,
    rate_R12,
)
from src.kinetics.char_reactions import (
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)
MODEL_MAP = {
    "bubble_homogeneous": ["R5b", "R6b", "R10b", "R11b", "R12b"],
    "dense_homogeneous": ["R5d", "R6d", "R7", "R8", "R10d", "R11d", "R12d", "R9"],
    "dense_heterogeneous": ["R1", "R2", "R3", "R4"],
    "reduction_like_bubble": ["R11b"],
    "reduction_like_dense": ["R7", "R8_reverse", "R11d"],
    "oxidation_bubble": ["R5b", "R6b", "R10b", "R12b"],
    "oxidation_dense": ["R5d", "R6d", "R10d", "R12d", "R1"],
    "note": "R8/R7 在源码中仅作用于 R_gas_d；气泡相无 WGSR/蒸汽重整。",
}


def _co_co2_from_stoich(stoich: dict[str, float], extent: float) -> tuple[float, float]:
    return float(stoich.get("CO", 0.0) * extent), float(stoich.get("CO2", 0.0) * extent)


def _analyze_cell_homogeneous(cell) -> dict:
    cell.calc_hydrodynamics()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau_val)
    cell.calc_exchange()

    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    y_c = cell._mole_fractions("combined")
    gas_src = cell._vm_gas_source_cache.copy()
    areas = cell._calc_char_surface_area_per_class()
    total_area = float(np.sum(areas))
    rho_cat = cell._catalyst_bulk_density()
    char_x = cell._compute_char_conversion()
    D_g = __import__("src.core.species", fromlist=["gas_diffusivity_correlation"]).gas_diffusivity_correlation(T, P)
    n_ex = cell.N_ex

    r5b = cell.r5_scale * rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], C_b[idx["H2O"]], P, y_b, idx) * cell.V_b
    r6b = cell.r6_scale * rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]]) * cell.V_b
    r12b = rate_R12(T, C_b[idx["H2"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b if cell.enable_r12 else 0.0
    r5d = cell.r5_scale * rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * cell.V_d
    r6d = cell.r6_scale * rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]]) * cell.V_d
    r12d = rate_R12(T, C_d[idx["H2"]], C_d[idx["O2"]], P, y_d, idx) * cell.V_d if cell.enable_r12 else 0.0
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

    c_tar_b = max(float(C_b[idx["TAR1"]] + C_b[idx["TAR2"]]), 0.0)
    c_tar_d = max(float(C_d[idx["TAR1"]] + C_d[idx["TAR2"]]), 0.0)
    ext10b = rate_R10(T, c_tar_b, C_b[idx["O2"]], P, cell.fuel_type) * cell.V_b
    ext10d = rate_R10(T, c_tar_d, C_d[idx["O2"]], P, cell.fuel_type) * cell.V_d
    ext11b = rate_R11_bubble(T, c_tar_b) * cell.V_b
    ext11d = rate_R11_suspension(T, c_tar_d, rho_cat) * cell.V_d
    sto10 = get_lumped_tar_stoichiometry("R10", cell.fuel_type)
    sto11 = get_lumped_tar_stoichiometry("R11", cell.fuel_type)

    # limiters (同 cell_kinetics)
    o2_supply_b = max(
        float(cell.N_zu_b[idx["O2"]]) + float(cell.N_b_in[idx["O2"]])
        + max(-float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_b[idx["O2"]]),
        1e-6,
    )
    o2_supply_d = max(
        float(cell.N_zu_d[idx["O2"]]) + float(cell.N_d_in[idx["O2"]])
        + max(float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_d[idx["O2"]]),
        1e-6,
    )
    nu_o2_10 = max(0.0, -float(sto10.get("O2", 0.0)))
    r1 = r2 = r3 = r4 = 0.0
    alpha = 1.0
    if total_area > 0.0:
        dc = d_core_from_spm_char_conversion(char_x, float(cell.solid.d_p))
        r1, alpha = rate_R1(T, C_d[idx["O2"]], float(cell.solid.d_p), D_g, dc, cell.fuel_type)
        r2 = rate_R2(T, C_d[idx["H2O"]], float(cell.solid.d_p), D_g, dc)
        r3 = rate_R3(T, C_d[idx["H2"]], float(cell.solid.d_p), D_g, dc)
        p_co2 = max(float(C_d[idx["CO2"]] * Rg * T), 0.0)
        p_co = max(float(C_d[idx["CO"]] * Rg * T), 0.0)
        r4 = cell.r4_scale * rate_R4_effective(T, p_co2, p_co, float(cell.solid.d_p), D_g, dc)

    # R5: CO + 1/2 O2 → CO2；extent = mol CO/s → O₂ demand = 0.5·r5
    bubble_o2_dem = 0.5 * r5b + 1.5 * r6b + 0.5 * r12b + nu_o2_10 * ext10b
    dense_o2_dem = r5d + 1.5 * r6d + 0.5 * r12d + alpha * r1 * total_area + nu_o2_10 * ext10d
    lim_o2_b = min(1.0, 0.995 * o2_supply_b / max(bubble_o2_dem, 1e-9))
    lim_o2_d = min(1.0, 0.995 * o2_supply_d / max(dense_o2_dem, 1e-9))
    co_supply_b = max(
        float(cell.N_zu_b[idx["CO"]]) + float(cell.N_b_in[idx["CO"]])
        + max(-float(n_ex[idx["CO"]]), 0.0) + float(cell.N_rez_b[idx["CO"]]),
        1e-6,
    )
    lim_co_b = min(1.0, 0.995 * co_supply_b / max(r5b, 1e-9))
    r5b_lim = r5b * min(lim_o2_b, lim_co_b)
    r6b_lim = r6b * lim_o2_b
    r12b_lim = r12b * lim_o2_b
    r5d_lim = r5d * lim_o2_d
    r6d_lim = r6d * lim_o2_d
    r12d_lim = r12d * lim_o2_d

    h2o_supply = max(
        float(cell.N_zu_d[idx["H2O"]]) + float(cell.N_d_in[idx["H2O"]])
        + float(cell.N_zu_b[idx["H2O"]]) + float(cell.N_b_in[idx["H2O"]])
        + float(cell.N_rez_d[idx["H2O"]]) + float(cell.N_rez_b[idx["H2O"]])
        + float(gas_src[idx["H2O"]]),
        1e-6,
    )
    h2o_consume = r2 * total_area + ext7 + max(0.0, ext8) + max(0.0, -float(sto11.get("H2O", 0.0))) * (ext11b + ext11d)
    lim_h2o = min(1.0, 0.995 * h2o_supply / max(h2o_consume, 1e-9))
    ext7_lim = ext7 * lim_h2o
    ext8_lim = ext8 * lim_h2o if ext8 > 0.0 else ext8
    ext10b_lim = ext10b * lim_o2_b
    ext10d_lim = ext10d * lim_o2_d
    ext11b_lim = ext11b * lim_h2o
    ext11d_lim = ext11d * lim_h2o
    r1_lim = r1 * lim_o2_d
    r2_lim = r2 * lim_h2o
    r4_lim = r4

    terms_b: dict[str, dict[str, float]] = {}
    terms_d: dict[str, dict[str, float]] = {}

    def _set(term: str, phase: str, dco: float, dco2: float) -> None:
        bucket = terms_b if phase == "b" else terms_d
        bucket[term] = {"CO": dco, "CO2": dco2}

    _set("R5", "b", -r5b_lim, r5b_lim)
    _set("R5", "d", -r5d_lim, r5d_lim)
    _set("R7", "d", ext7_lim, 0.0)
    _set("R8", "d", -ext8_lim, ext8_lim)
    co10b, co2_10b = _co_co2_from_stoich(sto10, ext10b_lim)
    co10d, co2_10d = _co_co2_from_stoich(sto10, ext10d_lim)
    _set("R10", "b", co10b, co2_10b)
    _set("R10", "d", co10d, co2_10d)
    co11b, co2_11b = _co_co2_from_stoich(sto11, ext11b_lim)
    co11d, co2_11d = _co_co2_from_stoich(sto11, ext11d_lim)
    _set("R11", "b", co11b, co2_11b)
    _set("R11", "d", co11d, co2_11d)

    if total_area > 0.0:
        r1_ext = r1_lim * total_area
        r2_ext = r2_lim * total_area
        r4_ext = r4_lim * total_area
        r1_o2 = alpha * r1_ext
        _set("R1", "d", 2.0 * (r1_ext - r1_o2), 2.0 * r1_o2 - r1_ext)
        _set("R2", "d", r2_ext, 0.0)
        _set("R4", "d", 2.0 * r4_ext, -r4_ext)

    _set("pyrolysis", "d", float(gas_src[idx["CO"]]), float(gas_src[idx["CO2"]]))

    def _sum_terms(tdict: dict) -> dict[str, float]:
        co = sum(v["CO"] for v in tdict.values())
        co2 = sum(v["CO2"] for v in tdict.values())
        return {"CO": co, "CO2": co2}

    return {
        "T_K": float(T),
        "V_b_m3": float(cell.V_b),
        "V_d_m3": float(cell.V_d),
        "V_b_frac": float(cell.V_b / max(cell.V_b + cell.V_d, 1e-12)),
        "eps_b": float(cell.eps_b),
        "K_bd": float(cell.K_bd),
        "y_O2_b": float(y_b[idx["O2"]]),
        "y_O2_d": float(y_d[idx["O2"]]),
        "y_CO_b": float(y_b[idx["CO"]]),
        "y_CO_d": float(y_d[idx["CO"]]),
        "y_CO2_b": float(y_b[idx["CO2"]]),
        "y_CO2_d": float(y_d[idx["CO2"]]),
        "y_CO_combined": float(y_c[idx["CO"]]),
        "y_CO2_combined": float(y_c[idx["CO2"]]),
        "N_ex_CO_mol_s": float(n_ex[idx["CO"]]),
        "N_ex_CO2_mol_s": float(n_ex[idx["CO2"]]),
        "bubble_terms": terms_b,
        "dense_terms": terms_d,
        "bubble_net": _sum_terms(terms_b),
        "dense_net": _sum_terms(terms_d),
        "raw": {
            "r5b_lim": float(r5b_lim),
            "r5d_lim": float(r5d_lim),
            "ext7_lim": float(ext7_lim),
            "ext8_lim": float(ext8_lim),
            "ext11b_lim": float(ext11b_lim),
            "ext11d_lim": float(ext11d_lim),
            "lim_o2_b": float(lim_o2_b),
            "lim_o2_d": float(lim_o2_d),
        },
    }


def _integrate(rows: list[dict], *, bubble: bool) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for row in rows:
        terms = row["bubble_terms"] if bubble else row["dense_terms"]
        for name, sp in terms.items():
            if name not in out:
                out[name] = {"CO": 0.0, "CO2": 0.0}
            out[name]["CO"] += float(sp["CO"])
            out[name]["CO2"] += float(sp["CO2"])
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_homogeneous_bubble_reduction_audit.json",
    )
    parser.add_argument("--max-global-iter", type=int, default=25)
    args = parser.parse_args()

    cfg = build_phase2_config()
    cfg.vorab_major_gibbs_x0 = False
    kw = phase2_solve_kwargs()
    kw["max_global_iter"] = int(args.max_global_iter)
    kw["verbose"] = False
    reactor = Reactor(cfg)
    result = reactor.solve(**kw)

    rows = []
    for i, cell in enumerate(reactor.cells):
        row = _analyze_cell_homogeneous(cell)
        row["cell_index"] = i
        row["xi_over_H"] = float(cell.geo.h_center / max(float(cfg.H_bed), 1e-12))
        rows.append(row)

    bub_int = _integrate(rows, bubble=True)
    den_int = _integrate(rows, bubble=False)
    anoxic = [r for r in rows if r["y_O2_b"] < 1e-4 and r["y_O2_d"] < 1e-4]
    bub_anoxic = _integrate(anoxic, bubble=True)
    den_anoxic = _integrate(anoxic, bubble=False)

    payload = {
        "case_key": "CASE_HTW_WESSELING_1",
        "model_map": MODEL_MAP,
        "converged_fully": bool(result.get("converged_fully")),
        "exit_gas_dry": result.get("exit_gas_dry"),
        "bed_integrated_bubble": bub_int,
        "bed_integrated_dense": den_int,
        "anoxic_integrated_bubble": bub_anoxic,
        "anoxic_integrated_dense": den_anoxic,
        "cells": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=" * 120)
    print("Phase2 homogeneous-phase (bubble vs dense) CO/CO2 audit")
    print("=" * 120)
    print(f"conv={payload['converged_fully']}  exit CO={payload['exit_gas_dry'].get('CO', 0):.4f} CO2={payload['exit_gas_dry'].get('CO2', 0):.4f}")
    print("\n[Model] R7/R8/WGSR/steam reforming → dense only; bubble has R5/R6/R10/R11/R12.")
    print("-" * 120)
    print("Bed-integrated CO/CO2 [mol/s] by phase and reaction:")
    print(f"  {'reaction':>6} {'bub CO':>10} {'bub CO2':>10} {'den CO':>10} {'den CO2':>10}")
    all_rxn = sorted(set(bub_int) | set(den_int))
    for rxn in all_rxn:
        b = bub_int.get(rxn, {"CO": 0.0, "CO2": 0.0})
        d = den_int.get(rxn, {"CO": 0.0, "CO2": 0.0})
        print(f"  {rxn:>6} {b['CO']:>10.3f} {b['CO2']:>10.3f} {d['CO']:>10.3f} {d['CO2']:>10.3f}")
    b_tot = {"CO": sum(v["CO"] for v in bub_int.values()), "CO2": sum(v["CO2"] for v in bub_int.values())}
    d_tot = {"CO": sum(v["CO"] for v in den_int.values()), "CO2": sum(v["CO2"] for v in den_int.values())}
    print(f"  {'TOTAL':>6} {b_tot['CO']:>10.3f} {b_tot['CO2']:>10.3f} {d_tot['CO']:>10.3f} {d_tot['CO2']:>10.3f}")

    print("-" * 120)
    print("Anoxic cells only (yO2_b,d < 1e-4):")
    for rxn in ("R5", "R7", "R8", "R11", "R2", "R4"):
        b = bub_anoxic.get(rxn, {"CO": 0.0, "CO2": 0.0})
        d = den_anoxic.get(rxn, {"CO": 0.0, "CO2": 0.0})
        if abs(b["CO"]) + abs(b["CO2"]) + abs(d["CO"]) + abs(d["CO2"]) > 1e-6:
            print(f"  {rxn:>4} bubble CO={b['CO']:+.4f} CO2={b['CO2']:+.4f}  |  dense CO={d['CO']:+.4f} CO2={d['CO2']:+.4f}")

    print("-" * 120)
    print(f"{'cell':>4} {'xi':>4} {'Vb%':>5} {'yCO_b':>7} {'yCO_d':>7} {'NexCO':>8} {'R5b_CO':>8} {'R5d_CO':>8} {'R8_CO':>8} {'R11b_CO':>8}")
    for row in rows:
        b = row["bubble_terms"]
        d = row["dense_terms"]
        print(
            f"{row['cell_index']:>4d} {row['xi_over_H']:>4.2f} {100*row['V_b_frac']:>5.1f} "
            f"{row['y_CO_b']:>7.4f} {row['y_CO_d']:>7.4f} {row['N_ex_CO_mol_s']:>8.3f} "
            f"{b.get('R5', {}).get('CO', 0.0):>8.3f} {d.get('R5', {}).get('CO', 0.0):>8.3f} "
            f"{d.get('R8', {}).get('CO', 0.0):>8.3f} {b.get('R11', {}).get('CO', 0.0):>8.3f}"
        )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
