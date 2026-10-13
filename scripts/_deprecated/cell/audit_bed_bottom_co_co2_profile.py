#!/usr/bin/env python3
"""审计床底各 cell 的 T、O2/CO/CO2 与分级进料 N_zu 关系（Phase2 LU staged_only）。"""

from __future__ import annotations

from _lu_harness import (
    REF_CO,
    REF_CO2,
    apply_baseline_kinetics,
    build_phase2_config,
    data_path,
    phase2_solve_kwargs,
    repo_root,
    write_json,
)

REPO_ROOT = repo_root()

import json
import sys
from pathlib import Path

import numpy as np

from src.core.cell_kinetics import build_reaction_sources
from src.core.cell_thermo_helpers import gas_mole_fractions
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
DRY_SPECIES = ("CO", "CO2", "H2", "CH4", "N2", "O2")


def _dry_y(N_b, N_d) -> dict[str, float]:
    y = gas_mole_fractions(N_b, N_d, "mixed")
    total = float(sum(y[idx[s]] for s in DRY_SPECIES))
    if total < 1e-15:
        return {s: 0.0 for s in DRY_SPECIES}
    return {s: float(y[idx[s]] / total) for s in DRY_SPECIES}


def _reaction_bundle(cell):
    cell.calc_hydrodynamics()
    cell.calc_exchange()
    thermo = cell._get_local_thermo_bundle()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    if cell._vm_cache_valid:
        gas_src_vm = cell._vm_gas_source_cache.copy()
        solid_sink_vm = cell._vm_solid_sink_cache.copy()
    else:
        gas_src_vm = cell._calc_drying_pyrolysis_gas_source(tau_val)
        solid_sink_vm = cell.R_solid.copy()
    areas = cell._calc_char_surface_area_per_class()
    gibbs_minor = None
    if cell.use_gibbs_minor:
        use_nitrogen = not bool(getattr(cell, "_nr_solves_nh3_transport", False))
        gibbs_minor = cell.calc_minor_species_gibbs(use_nitrogen=use_nitrogen)
    return build_reaction_sources(
        T=cell.T,
        P=cell.P,
        fuel_type=cell.fuel_type,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=np.asarray(thermo["C_b"], dtype=np.float64),
        C_d=np.asarray(thermo["C_d"], dtype=np.float64),
        y_b=np.asarray(thermo["y_b"], dtype=np.float64),
        y_d=np.asarray(thermo["y_d"], dtype=np.float64),
        gas_src_vm=gas_src_vm,
        solid_sink_vm=solid_sink_vm,
        areas=areas,
        solid_d_p=cell.solid.d_p,
        solid_d_p_classes=cell.solid.d_p_classes,
        D_g=float(thermo["D_g"]),
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        rho_ash_kg_m3=cell._ash_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=gibbs_minor,
        r4_scale=cell.r4_scale,
        r2_scale=cell.r2_scale,
        r5_scale=cell.r5_scale,
        r6_scale=cell.r6_scale,
        r7_scale=cell.r7_scale,
        r8_scale=cell.r8_scale,
        rate_multiplier=1.0,
        N_zu_d=cell.N_zu_d,
        N_d_in=cell.N_d_in,
        N_zu_b=cell.N_zu_b,
        N_b_in=cell.N_b_in,
        N_rez_d=cell.N_rez_d,
        N_rez_b=cell.N_rez_b,
        N_ex=cell.N_ex,
        solid_shape=cell.R_solid.shape,
        char_index=0,
    )


def _configure(label: str, cfg) -> None:
    for s in ("r2", "r4", "r5", "r6", "r7", "r8"):
        setattr(cfg, f"{s}_scale", 1.0)
    if label == "staged_only":
        cfg.vorab_staged_primary_gas_injection_thesis = True
        cfg.vorab_staged_primary_gas_bed_cells_thesis = 3
        cfg.vorab_a_tier_bed_axial_primary_profile_thesis = True
        cfg.vorab_bed_temperature_march_thesis = False
    elif label == "legacy":
        cfg.vorab_staged_primary_gas_injection_thesis = False
        cfg.vorab_a_tier_bed_axial_primary_profile_thesis = False
        cfg.vorab_bed_temperature_march_thesis = False


def audit_scenario(label: str, *, max_global_iter: int) -> dict:
    cfg = build_phase2_config(enable_vorab_march_staged=False)
    _configure(label, cfg)
    reactor = Reactor(cfg)
    kw = phase2_solve_kwargs()
    kw["max_global_iter"] = int(max_global_iter)
    kw["verbose"] = False
    out = reactor.solve(**kw)

    h_bed = float(cfg.H_bed)
    o2_feed = float(cfg.O2_feed)
    rows = []
    bed_idx = [i for i, c in enumerate(reactor.cells) if str(getattr(c, "cell_type", "bed")) == "bed"]
    for bi in bed_idx:
        cell = reactor.cells[bi]
        z = float(cell.geo.h_center)
        xi = z / h_bed if h_bed > 0 else 0.0
        y_d = _dry_y(cell.N_b, cell.N_d)
        y_d_b = _dry_y(np.zeros_like(cell.N_b), cell.N_d) if float(np.sum(cell.N_d)) > 0 else _dry_y(cell.N_b, cell.N_d)
        y_d_bub = _dry_y(cell.N_b, np.zeros_like(cell.N_d)) if float(np.sum(cell.N_b)) > 0 else _dry_y(cell.N_b, cell.N_d)
        o2_zu = float(cell.N_zu_d[idx["O2"]] + cell.N_zu_b[idx["O2"]])
        bundle = _reaction_bundle(cell)
        co_prod = float(bundle.R_gas_b[idx["CO"]] + bundle.R_gas_d[idx["CO"]])
        co2_prod = float(bundle.R_gas_b[idx["CO2"]] + bundle.R_gas_d[idx["CO2"]])
        co_co2 = y_d["CO"] + y_d["CO2"]
        co_share = y_d["CO"] / co_co2 if co_co2 > 1e-12 else float("nan")
        rows.append(
            {
                "bed_i": int(bi),
                "z_m": z,
                "xi": xi,
                "T_K": float(cell.T),
                "eps_b": float(cell.eps_b),
                "u_b_m_s": float(cell.u_b),
                "u_d_m_s": float(cell.u_d),
                "O2_zu_mol_s": o2_zu,
                "O2_zu_frac_of_feed": o2_zu / o2_feed if o2_feed > 0 else 0.0,
                "solid_feed": bi == bed_idx[0],
                "staged_injection": o2_zu > 1e-9,
                "y_dry_CO": y_d["CO"],
                "y_dry_CO2": y_d["CO2"],
                "y_dry_O2": y_d["O2"],
                "y_dry_H2": y_d["H2"],
                "CO_share_in_carbonox": co_share,
                "y_dense_CO": y_d_b["CO"],
                "y_dense_CO2": y_d_b["CO2"],
                "y_dense_O2": y_d_b["O2"],
                "y_bubble_CO": y_d_bub["CO"],
                "y_bubble_CO2": y_d_bub["CO2"],
                "y_bubble_O2": y_d_bub["O2"],
                "char_x": float(cell._compute_char_conversion()),
                "extent_R1_mol_s": float(bundle.extent_r1),
                "extent_R2_mol_s": float(bundle.extent_r2),
                "extent_R4_mol_s": float(bundle.extent_r4),
                "extent_R8_net": float(bundle.net_molar_gas_source_r8),
                "R_CO_mol_s": co_prod,
                "R_CO2_mol_s": co2_prod,
                "o2_budget_dense": float(bundle.o2_budget_dense),
                "o2_demand_dense": float(bundle.o2_demand_dense_limited),
                "o2_slack_dense": float(bundle.o2_slack_dense),
                "limit_o2_dense": float(bundle.limit_factor_o2_dense),
            }
        )

    exit_y = out.get("exit_gas_dry") or {}
    return {
        "label": label,
        "converged_fully": bool(out.get("converged_fully", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", float("nan"))),
        "co_exit": float(exit_y.get("CO", 0.0)),
        "co2_exit": float(exit_y.get("CO2", 0.0)),
        "O2_feed_mol_s": o2_feed,
        "bed_rows": rows,
    }


def main() -> None:
    max_iter = 25
    results = []
    for label in ("staged_only", "legacy"):
        print(f"=== solving {label} ===", flush=True)
        results.append(audit_scenario(label, max_global_iter=max_iter))

    out_path = REPO_ROOT / "data" / "phase2_bed_bottom_co_co2_audit.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    for res in results:
        print(
            f"\n{res['label']:12s} conv={res['converged_fully']} "
            f"CO={res['co_exit']:.4f} CO2={res['co2_exit']:.4f} "
            f"(ref {REF_CO}/{REF_CO2})",
            flush=True,
        )
        print(
            f"{'bed':>3} {'xi':>5} {'T_K':>7} {'O2_zu%':>7} "
            f"{'yCO':>6} {'yCO2':>6} {'yO2':>6} {'CO%':>5} "
            f"{'R1':>8} {'R2':>8} {'R4':>8} {'dCO':>8} {'dCO2':>8}",
            flush=True,
        )
        for row in res["bed_rows"][:6]:
            print(
                f"{row['bed_i']:3d} {row['xi']:5.2f} {row['T_K']:7.0f} "
                f"{100*row['O2_zu_frac_of_feed']:6.1f} "
                f"{row['y_dry_CO']:6.3f} {row['y_dry_CO2']:6.3f} {row['y_dry_O2']:6.3f} "
                f"{100*row['CO_share_in_carbonox']:4.0f} "
                f"{row['extent_R1_mol_s']:8.3f} {row['extent_R2_mol_s']:8.3f} "
                f"{row['extent_R4_mol_s']:8.3f} {row['R_CO_mol_s']:8.3f} {row['R_CO2_mol_s']:8.3f}",
                flush=True,
            )
    print(f"\nwrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
