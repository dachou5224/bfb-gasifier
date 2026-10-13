#!/usr/bin/env python3
"""R1 异相氧化偏慢诊断：速率分解 × 气泡/乳化相 O₂ 分配。

对 bed0–2：
1. y_O2 / C_O2 气泡 vs 乳化；O₂ 进料 dense 份额；N_ex；K_bd
2. R1 仅见 C_d[O2]；假想改用 C_b / 混合相时 extent 倍率（诊断用，非改模型）
3. 串联阻力：动力学 vs 气膜/灰层；φ_c
4. 与 R5b 的 O₂ 消耗对比；lf / slack

Usage
-----
    python3 scripts/audit_r1_o2_phase_bottleneck_lu.py
    python3 scripts/audit_r1_o2_phase_bottleneck_lu.py --max-global-iter 8 --seed-T0 1250
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    build_phase2_reactor,
    phase2_solve_kwargs,
    timed_phase2_solve,
    write_json,
)
import src.workflow.steps.init_precalc_step as ips
from src.core.cell_kinetics import build_reaction_sources
from src.core.constants import Rg
from src.core.cell import S_CHAR
from src.core.species import GAS_SPECIES_INDEX as idx, fuel_tar_target_hc_ratio
from src.kinetics.arrhenius import k_hobbs
from src.kinetics.char_reactions import (
    D_A_ref,
    D_d_A,
    _k_diff_sphere,
    d_core_from_spm_char_conversion,
    phi_c,
    r1_hamel_kinetic_constants,
    rate_R1,
)


def _r1_resistance_breakdown(
    T: float,
    C_O2: float,
    d_p: float,
    D_g: float,
    d_core: float,
    fuel: str = "coal",
) -> dict[str, float | str]:
    C_O2 = max(float(C_O2), 0.0)
    P_O2 = C_O2 * Rg * T
    d_c = max(min(float(d_core), float(d_p)), 1e-9)
    phi = float(phi_c(d_c, T))
    k0, e = r1_hamel_kinetic_constants("coal" if fuel == "coal" else "biomass")
    k_ch = float(k_hobbs(k0, e, T))
    RT = Rg * T
    k_d_g_p = _k_diff_sphere(D_g, d_p) / RT
    D_d_A_p = D_d_A(D_A_ref, 2.0, d_p, d_c) / RT
    res_gas = (1.0 / max(k_d_g_p, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_p, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_ext_p = 1.0 / max((1.0 / phi) * (res_gas + res_ash), 1e-30)
    r, alpha = rate_R1(T, C_O2, d_p, D_g, d_c, fuel)  # type: ignore[arg-type]
    r_kin_only = k_ch * P_O2
    r_diff_only = k_ext_p * P_O2
    ratio = float(r_kin_only / max(r_diff_only, 1e-30))
    if r_kin_only < 0.3 * r_diff_only:
        regime = "kinetic"
    elif r_diff_only < 0.3 * r_kin_only:
        regime = "diffusion"
    else:
        regime = "mixed"
    return {
        "phi_c": phi,
        "alpha": float(alpha),
        "k_ch": k_ch,
        "k_ext_p": float(k_ext_p),
        "res_gas": float(res_gas),
        "res_ash": float(res_ash),
        "P_O2_Pa": float(P_O2),
        "r_mol_m2_s": float(r),
        "r_kin_only_approx": float(r_kin_only),
        "r_diff_only_approx": float(r_diff_only),
        "kin_over_diff": ratio,
        "regime": regime,
    }


def _bundle_for_cell(cell: Any):
    thermo = cell._get_local_thermo_bundle()
    C_b = np.asarray(thermo["C_b"], dtype=np.float64)
    C_d = np.asarray(thermo["C_d"], dtype=np.float64)
    y_b = np.asarray(thermo["y_b"], dtype=np.float64)
    y_d = np.asarray(thermo["y_d"], dtype=np.float64)
    D_g = float(thermo["D_g"])
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    if cell._vm_cache_valid:
        gas_src_vm = cell._vm_gas_source_cache.copy()
        solid_sink_vm = cell._vm_solid_sink_cache.copy()
    else:
        gas_src_vm = cell._calc_drying_pyrolysis_gas_source(tau_val)
        solid_sink_vm = cell.R_solid.copy()
    areas = cell._calc_char_surface_area_per_class()
    target_tar_hc = None
    if cell.solid.C_dry > 1e-12 and cell.solid.H_dry >= 0.0:
        target_tar_hc = fuel_tar_target_hc_ratio(cell.solid.C_dry, cell.solid.H_dry)
    return build_reaction_sources(
        T=float(cell.T),
        P=float(cell.P),
        fuel_type=cell.fuel_type,
        target_tar_hc_ratio=target_tar_hc,
        V_b=float(cell.V_b),
        V_d=float(cell.V_d),
        C_b=C_b,
        C_d=C_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=gas_src_vm,
        solid_sink_vm=solid_sink_vm,
        areas=areas,
        solid_d_p=float(cell.solid.d_p),
        solid_d_p_classes=cell.solid.d_p_classes,
        D_g=D_g,
        char_conversion=float(cell._compute_char_conversion()),
        rho_cat=float(cell._catalyst_bulk_density()),
        rho_ash_kg_m3=float(cell._ash_bulk_density()),
        enable_r12=bool(cell.enable_r12),
        use_gibbs_minor=bool(cell.use_gibbs_minor),
        gibbs_minor_sources=None,
        r4_scale=float(cell.r4_scale),
        r2_scale=float(cell.r2_scale),
        r5_scale=float(cell.r5_scale),
        r5_bubble_scale=float(cell.r5_bubble_scale),
        r6_scale=float(cell.r6_scale),
        r7_scale=float(cell.r7_scale),
        r8_scale=float(cell.r8_scale),
        rate_multiplier=1.0,
        N_zu_d=cell.N_zu_d,
        N_d_in=cell.N_d_in,
        N_zu_b=cell.N_zu_b,
        N_b_in=cell.N_b_in,
        N_rez_d=cell.N_rez_d,
        N_rez_b=cell.N_rez_b,
        N_ex=cell.N_ex,
        solid_shape=cell.R_solid.shape,
        char_index=S_CHAR,
    ), thermo, areas, D_g


def _audit_cell(cell: Any, bed_rank: int) -> dict[str, Any]:
    cell.calc_hydrodynamics()
    tau = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau)
    cell.calc_exchange()
    bundle, thermo, areas, D_g = _bundle_for_cell(cell)
    y_b = np.asarray(thermo["y_b"], dtype=np.float64)
    y_d = np.asarray(thermo["y_d"], dtype=np.float64)
    C_b = np.asarray(thermo["C_b"], dtype=np.float64)
    C_d = np.asarray(thermo["C_d"], dtype=np.float64)
    T = float(cell.T)
    X = float(cell._compute_char_conversion())
    d_p = np.asarray(cell.solid.d_p_classes, dtype=np.float64)
    fuel = str(cell.fuel_type)
    areas_arr = np.asarray(areas, dtype=np.float64)

    r1_d = r1_b_hyp = r1_mix_hyp = 0.0
    C_mix = 0.5 * (float(C_b[idx["O2"]]) + float(C_d[idx["O2"]]))
    class_rows: list[dict[str, Any]] = []
    for k, (a, dp) in enumerate(zip(areas_arr, d_p)):
        if a <= 0.0 or dp <= 0.0:
            continue
        dc = d_core_from_spm_char_conversion(X, float(dp))
        rd, _ = rate_R1(T, float(C_d[idx["O2"]]), float(dp), D_g, dc, fuel)  # type: ignore[arg-type]
        rb, _ = rate_R1(T, float(C_b[idx["O2"]]), float(dp), D_g, dc, fuel)  # type: ignore[arg-type]
        rm, _ = rate_R1(T, C_mix, float(dp), D_g, dc, fuel)  # type: ignore[arg-type]
        r1_d += float(rd) * float(a)
        r1_b_hyp += float(rb) * float(a)
        r1_mix_hyp += float(rm) * float(a)
        bd = _r1_resistance_breakdown(T, float(C_d[idx["O2"]]), float(dp), D_g, dc, fuel)
        class_rows.append(
            {
                "class": k,
                "d_p_m": float(dp),
                "d_core_m": float(dc),
                "area_m2": float(a),
                "extent_d": float(rd) * float(a),
                **bd,
            }
        )
    top = max(class_rows, key=lambda r: r["area_m2"]) if class_rows else None
    j = idx["O2"]
    o2_zu_d = float(cell.N_zu_d[j])
    o2_zu_b = float(cell.N_zu_b[j])
    o2_feed = o2_zu_d + o2_zu_b
    return {
        "bed_rank": bed_rank,
        "T_K": T,
        "char_conversion_X": X,
        "eps_b": float(cell.eps_b),
        "V_b": float(cell.V_b),
        "V_d": float(cell.V_d),
        "K_bd_1_s": float(cell.K_bd),
        "D_g_m2_s": D_g,
        "y_O2_bubble": float(y_b[j]),
        "y_O2_dense": float(y_d[j]),
        "y_O2_ratio_b_over_d": float(y_b[j] / max(float(y_d[j]), 1e-30)),
        "C_O2_bubble": float(C_b[j]),
        "C_O2_dense": float(C_d[j]),
        "N_zu_b_O2": o2_zu_b,
        "N_zu_d_O2": o2_zu_d,
        "dense_frac_O2_feed": float(o2_zu_d / max(o2_feed, 1e-30)) if o2_feed > 0 else float("nan"),
        "N_ex_O2": float(cell.N_ex[j]),
        "N_b_O2": float(cell.N_b[j]),
        "N_d_O2": float(cell.N_d[j]),
        "A_char_total_m2": float(np.sum(areas_arr)),
        "R1_extent_manual": r1_d,
        "R1_extent_if_use_C_O2_bubble": r1_b_hyp,
        "R1_extent_if_use_C_O2_mix": r1_mix_hyp,
        "R1_boost_if_bubble_C_O2": float(r1_b_hyp / max(r1_d, 1e-30)),
        "lf_o2_b": float(bundle.limit_factor_o2_bubble),
        "lf_o2_d": float(bundle.limit_factor_o2_dense),
        "o2_supply_b": float(bundle.o2_supply_bubble),
        "o2_supply_d": float(bundle.o2_supply_dense),
        "o2_demand_d_lim": float(bundle.o2_demand_dense_limited),
        "o2_slack_d": float(bundle.o2_slack_dense),
        "extent_R1": float(bundle.extent_r1),
        "extent_R1_O2": float(bundle.extent_r1_o2),
        "extent_R5b": float(bundle.extent_r5_bubble),
        "extent_R5d": float(bundle.extent_r5_dense),
        "extent_R6b": float(bundle.extent_r6_bubble),
        "extent_R6d": float(bundle.extent_r6_dense),
        "O2_sink_R5b_over_R1": float(bundle.extent_r5_bubble / max(bundle.extent_r1_o2, 1e-30)),
        "top_class": top,
        "classes": class_rows,
    }


def run_solve(*, profile: str, T0_seed: float | None, max_iter: int):
    _orig = ips._install_temperature_fence_from_vorabrechnung

    def _install(reactor, T_est):
        _orig(reactor, T_est)
        if T0_seed is None:
            return
        beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if beds:
            b0 = beds[0]
            b0._nr_temperature_min_K = float(max(700.0, T0_seed - 200.0))
            b0._nr_temperature_max_K = float(min(2500.0, max(T0_seed + 400.0, 1400.0)))
            b0.T = float(T0_seed)
            b0._vorab_temperature_estimate_K = float(T0_seed)

    overrides: dict[str, object] = {
        "vorab_vm_devolatilization_zone_cells_thesis": 3,
        "vorab_vm_devolatilization_zone_profile_thesis": profile,
    }
    if T0_seed is not None:
        overrides["vorab_nr_temperature_seed_bed0_cold_anchor_thesis"] = False
        overrides["vorab_nr_temperature_seed_bottom_K_thesis"] = float(T0_seed)

    ips._install_temperature_fence_from_vorabrechnung = _install
    try:
        r = build_phase2_reactor(
            ladder_level="p2",
            builder_kwargs={"enable_vorab_march_staged": True},
            cfg_overrides=overrides,
        )
        out, wall = timed_phase2_solve(r, **phase2_solve_kwargs(max_global_iter=max_iter))
        return r, out, wall
    finally:
        ips._install_temperature_fence_from_vorabrechnung = _orig


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(p, default=8)
    add_json_output_arg(p, "phase2_r1_o2_phase_bottleneck.json")
    p.add_argument("--profile", default="front_loaded", choices=("front_loaded", "linear"))
    p.add_argument("--seed-T0", type=float, default=1250.0)
    args = p.parse_args()
    seed = None if float(args.seed_T0) < 0 else float(args.seed_T0)

    print(
        f"=== R1/O2 audit profile={args.profile} seed={seed} iter={args.max_global_iter} ===",
        flush=True,
    )
    reactor, result, wall = run_solve(
        profile=args.profile, T0_seed=seed, max_iter=int(args.max_global_iter)
    )
    beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    rows = []
    for i, cell in enumerate(beds[:3]):
        row = _audit_cell(cell, i)
        rows.append(row)
        top = row.get("top_class") or {}
        print(
            f"\nbed{i}: T={row['T_K']:.0f}K  yO2_b={row['y_O2_bubble']:.4f} "
            f"yO2_d={row['y_O2_dense']:.4f} (b/d={row['y_O2_ratio_b_over_d']:.2f})  "
            f"K_bd={row['K_bd_1_s']:.2f}  Xc={row['char_conversion_X']:.3f}  "
            f"A_char={row['A_char_total_m2']:.1f} m2",
            flush=True,
        )
        print(
            f"  O2 feed: zu_b={row['N_zu_b_O2']:.2f} zu_d={row['N_zu_d_O2']:.2f} "
            f"dense_frac={row['dense_frac_O2_feed']}  N_ex={row['N_ex_O2']:+.2f}  "
            f"supply_d={row['o2_supply_d']:.2f} slack_d={row['o2_slack_d']:.2f} lf_d={row['lf_o2_d']:.3f}",
            flush=True,
        )
        print(
            f"  R1={row['extent_R1']:.3f} molC/s  R1_O2={row['extent_R1_O2']:.3f}  "
            f"R5b={row['extent_R5b']:.3f}  R5b/R1_O2={row['O2_sink_R5b_over_R1']:.1f}",
            flush=True,
        )
        print(
            f"  假想 C_O2,b→R1×{row['R1_boost_if_bubble_C_O2']:.2f}；"
            f"mix→×{row['R1_extent_if_use_C_O2_mix']/max(row['R1_extent_manual'],1e-30):.2f}",
            flush=True,
        )
        if top:
            print(
                f"  主粒径: dp={top['d_p_m']*1e3:.2f}mm dc={top['d_core_m']*1e3:.2f}mm "
                f"regime={top['regime']} kin/diff={top['kin_over_diff']:.3f} "
                f"phi={top['phi_c']:.3f} res_g={top['res_gas']:.2e} res_ash={top['res_ash']:.2e}",
                flush=True,
            )

    payload = {
        "note": "R1 bottleneck vs bubble/dense O2; C_O2,b hypothetical is diagnostic only",
        "profile": args.profile,
        "T0_seed_K": seed,
        "wall_s": wall,
        "converged_fully": bool(result.get("converged_fully", False)),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
        "cells": rows,
        "hamel_notes": [
            "Hamel §5.2: heterogeneous char reactions ONLY in suspension (dense) phase",
            "Hamel p.89: O2 mainly consumed by gas-phase oxidation; char mainly via R2–R4",
            "Eq.5.17 Sh=2+0.6 Re^0.5 Sc^0.33 — code uses Sh=2 constant",
            "Eq.5.22 particle T: Hamel says effect small in gasifiers",
        ],
    }
    out = write_json(args.output, payload)
    print(f"\nWrote {out}", flush=True)


if __name__ == "__main__":
    main()
