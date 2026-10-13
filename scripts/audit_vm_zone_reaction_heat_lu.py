#!/usr/bin/env python3
"""下部 VM 释放区（默认 bed0–2）反应 extent + 放热去向审计。

在吸引子态（默认 front_loaded + soft seed1250）下：
1. 限幅后 R1–R12 extent [mol/s]
2. 各反应化学放热 Q = −extent·ΔH(T) [MW]（放热为正）
3. Eq.2.7 焓流分解：冷进料 / 出流 / 壁损 / 残差

Usage
-----
    PYTHONPATH=.:scripts python3 scripts/audit_vm_zone_reaction_heat_lu.py
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

import src.workflow.steps.init_precalc_step as ips
from _lu_harness import (
    add_json_output_arg,
    add_max_global_iter_arg,
    bed_cell_indices,
    build_phase2_reactor,
    phase2_solve_kwargs,
    timed_phase2_solve,
    write_json,
)
from src.core.cell import S_CHAR
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.species import GAS_SPECIES_INDEX as idx, enthalpy_molar, formation_enthalpy_298
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry


def _h(sp: str, T: float) -> float:
    return float(enthalpy_molar(sp, float(T)))


def _delta_h_gas(stoich: dict[str, float], T: float) -> float:
    """Σ ν_i h_i(T) [J/mol]，ν>0 为产物。"""
    return float(sum(float(nu) * _h(sp, T) for sp, nu in stoich.items()))


def _char_molar_enthalpy(cell: Any, T: float) -> float:
    """1 mol C（以 char 固体计）在 T 下的焓 [J/mol]。"""
    m = np.zeros_like(cell.m_solid, dtype=np.float64)
    if m.ndim == 1:
        m = m.reshape(1, -1)
    m = np.zeros((1, m.shape[1]), dtype=np.float64)
    m[0, S_CHAR] = float(ATOMIC_MASS_KG_PER_MOL["C"])
    return float(cell._calc_solid_enthalpy_flow(m, float(T)))


def _reaction_heat_table(cell: Any, bundle: Any) -> dict[str, dict[str, float]]:
    """各反应 extent 与放热（放热为正）[mol/s, MW]。"""
    T = float(cell.T)
    h_c = _char_molar_enthalpy(cell, T)
    alpha = (
        float(bundle.extent_r1_o2) / max(float(bundle.extent_r1), 1e-30)
        if float(bundle.extent_r1) > 1e-30
        else 1.0
    )
    # R1: C + α O2 → 2(1−α) CO + (2α−1) CO2
    dH_r1 = (
        2.0 * (1.0 - alpha) * _h("CO", T)
        + (2.0 * alpha - 1.0) * _h("CO2", T)
        - alpha * _h("O2", T)
        - h_c
    )
    dH_r2 = _h("CO", T) + _h("H2", T) - _h("H2O", T) - h_c  # C+H2O→CO+H2
    dH_r3 = _h("CH4", T) - 2.0 * _h("H2", T) - h_c  # C+2H2→CH4
    dH_r4 = 2.0 * _h("CO", T) - _h("CO2", T) - h_c  # C+CO2→2CO
    dH_r5 = _delta_h_gas({"CO": -2.0, "O2": -1.0, "CO2": 2.0}, T)
    dH_r6 = _delta_h_gas({"CH4": -1.0, "O2": -1.5, "CO": 1.0, "H2O": 2.0}, T)
    dH_r7 = _delta_h_gas({"CH4": -1.0, "H2O": -1.0, "CO": 1.0, "H2": 3.0}, T)
    dH_r8 = _delta_h_gas({"CO": -1.0, "H2O": -1.0, "CO2": 1.0, "H2": 1.0}, T)
    dH_r12 = _delta_h_gas({"H2": -1.0, "O2": -0.5, "H2O": 1.0}, T)

    fuel = str(getattr(cell, "fuel_type", "coal"))
    stoich10 = get_lumped_tar_stoichiometry("R10", fuel)
    stoich11 = get_lumped_tar_stoichiometry("R11", fuel)
    dH_r10 = _delta_h_gas({k: float(v) for k, v in stoich10.items() if k in idx}, T)
    dH_r11 = _delta_h_gas({k: float(v) for k, v in stoich11.items() if k in idx}, T)

    rows: dict[str, dict[str, float]] = {}

    def _add(name: str, extent: float, dH: float) -> None:
        ext = float(extent)
        q_w = -ext * float(dH)  # 放热为正
        rows[name] = {
            "extent_mol_s": ext,
            "dH_J_per_mol": float(dH),
            "Q_release_MW": q_w / 1e6,
            "dH_298_ref_J_per_mol": float("nan"),
        }

    _add("R1", bundle.extent_r1, dH_r1)
    _add("R2", bundle.extent_r2, dH_r2)
    _add("R3", bundle.extent_r3, dH_r3)
    _add("R4", bundle.extent_r4, dH_r4)
    _add("R5b", bundle.extent_r5_bubble, dH_r5)
    _add("R5d", bundle.extent_r5_dense, dH_r5)
    _add("R6b", bundle.extent_r6_bubble, dH_r6)
    _add("R6d", bundle.extent_r6_dense, dH_r6)
    _add("R7", bundle.extent_r7, dH_r7)
    _add("R8", bundle.extent_r8, dH_r8)
    _add("R10b", bundle.extent_r10_bubble, dH_r10)
    _add("R10d", bundle.extent_r10_dense, dH_r10)
    _add("R11b", bundle.extent_r11_bubble, dH_r11)
    _add("R11d", bundle.extent_r11_dense, dH_r11)
    _add("R12b", bundle.extent_r12_bubble, dH_r12)
    _add("R12d", bundle.extent_r12_dense, dH_r12)

    # 298 K 参考（仅气体反应，便于对照）
    for name, stoich in (
        ("R5b", {"CO": -2.0, "O2": -1.0, "CO2": 2.0}),
        ("R5d", {"CO": -2.0, "O2": -1.0, "CO2": 2.0}),
        ("R6b", {"CH4": -1.0, "O2": -1.5, "CO": 1.0, "H2O": 2.0}),
        ("R6d", {"CH4": -1.0, "O2": -1.5, "CO": 1.0, "H2O": 2.0}),
    ):
        if name in rows:
            rows[name]["dH_298_ref_J_per_mol"] = float(
                sum(float(nu) * formation_enthalpy_298(sp) for sp, nu in stoich.items())
            )
    return rows


def _enthalpy_parts(cell: Any) -> dict[str, float]:
    """Eq.2.7 焓流分量 [MW]。"""
    parts = {
        "gas_axial_in": cell._calc_gas_enthalpy_flow(cell.N_b_in + cell.N_d_in, cell.T_in_gas),
        "gas_zu": cell._calc_gas_enthalpy_flow(cell.N_zu_b + cell.N_zu_d, cell.T_zu_gas),
        "gas_rez": cell._calc_gas_enthalpy_flow(cell.N_rez_b + cell.N_rez_d, cell.T_rez_gas),
        "solid_in": cell._calc_solid_enthalpy_flow(cell.m_solid_in, cell.T_in_solid),
        "solid_zu": cell._calc_solid_enthalpy_flow(cell.m_solid_zu, cell.T_zu_solid),
        "solid_rez": cell._calc_solid_enthalpy_flow(cell.m_solid_rez, cell.T_rez_solid),
        "solid_auf_in": cell._calc_solid_enthalpy_flow(cell.m_solid_auf_in, cell.T_solid_auf_in),
        "solid_ab_in": cell._calc_solid_enthalpy_flow(cell.m_solid_ab_in, cell.T_solid_ab_in),
    }
    hin = sum(parts.values())
    hout_g = cell._calc_gas_enthalpy_flow(cell.N_b + cell.N_d, cell.T)
    hout_s = cell._calc_solid_enthalpy_flow(cell._solid_energy_outlet_rates(), cell.T)
    q = float(cell.Q_wall_W) + float(cell.Q_heat_exchanger_W)
    e = float(cell.calc_energy_balance())
    return {
        **{f"{k}_MW": v / 1e6 for k, v in parts.items()},
        "Hin_MW": hin / 1e6,
        "Hout_g_MW": hout_g / 1e6,
        "Hout_s_MW": hout_s / 1e6,
        "Q_wall_MW": q / 1e6,
        "E_residual_MW": e / 1e6,
        "Hin_minus_Hout_minus_Q_MW": (hin - hout_g - hout_s - q) / 1e6,
        "Tin_gas_K": float(cell.T_in_gas),
        "T_zu_gas_K": float(cell.T_zu_gas),
        "T_zu_solid_K": float(cell.T_zu_solid),
        "T_K": float(cell.T),
    }


def _vm_fuels(cell: Any) -> dict[str, float]:
    vm = np.asarray(cell._vm_gas_source_cache, dtype=np.float64)
    return {
        "sum_mol_s": float(np.sum(vm)),
        "CO": float(vm[idx["CO"]]),
        "H2": float(vm[idx["H2"]]),
        "CH4": float(vm[idx["CH4"]]),
        "TAR1": float(vm[idx["TAR1"]]),
        "TAR2": float(vm[idx["TAR2"]]),
        "H2O": float(vm[idx["H2O"]]),
        "combustible": float(
            vm[idx["CO"]] + vm[idx["H2"]] + vm[idx["CH4"]] + vm[idx["TAR1"]] + vm[idx["TAR2"]]
        ),
    }


def _get_bundle(cell: Any) -> Any:
    """与 cell.calc_reactions 同口径取 bundle（含 extent）。"""
    from src.core.cell_kinetics import build_reaction_sources
    from src.core.species import fuel_tar_target_hc_ratio

    cell.calc_hydrodynamics()
    thermo = cell._get_local_thermo_bundle()
    C_b = np.asarray(thermo["C_b"], dtype=np.float64)
    C_d = np.asarray(thermo["C_d"], dtype=np.float64)
    y_b = np.asarray(thermo["y_b"], dtype=np.float64)
    y_d = np.asarray(thermo["y_d"], dtype=np.float64)
    D_g = float(thermo["D_g"])
    gas_src_vm = cell._vm_gas_source_cache.copy()
    solid_sink_vm = cell._vm_solid_sink_cache.copy()
    areas = cell._calc_char_surface_area_per_class()
    target_tar_hc = None
    if cell.solid.C_dry > 1e-12 and cell.solid.H_dry >= 0.0:
        target_tar_hc = fuel_tar_target_hc_ratio(cell.solid.C_dry, cell.solid.H_dry)
    gibbs_minor = None
    if cell.use_gibbs_minor:
        use_nitrogen = not bool(getattr(cell, "_nr_solves_nh3_transport", False))
        gibbs_minor = cell.calc_minor_species_gibbs(use_nitrogen=use_nitrogen)
    return build_reaction_sources(
        T=float(cell.T),
        P=float(cell.P),
        fuel_type=cell.fuel_type,
        target_tar_hc_ratio=target_tar_hc,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=C_b,
        C_d=C_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=gas_src_vm,
        solid_sink_vm=solid_sink_vm,
        areas=areas,
        solid_d_p=cell.solid.d_p,
        solid_d_p_classes=cell.solid.d_p_classes,
        m_char_classes=np.maximum(cell.m_solid[:, S_CHAR], 0.0),
        D_g=D_g,
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        rho_ash_kg_m3=cell._ash_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=gibbs_minor,
        r4_scale=cell.r4_scale,
        r2_scale=cell.r2_scale,
        r5_scale=cell.r5_scale,
        r5_bubble_scale=cell.r5_bubble_scale,
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
        char_index=S_CHAR,
    )


def analyze_cell(cell: Any, bed_rank: int) -> dict[str, Any]:
    bundle = _get_bundle(cell)
    heats = _reaction_heat_table(cell, bundle)
    ranked = sorted(heats.items(), key=lambda kv: abs(kv[1]["Q_release_MW"]), reverse=True)
    q_ox = sum(
        heats[k]["Q_release_MW"]
        for k in ("R1", "R5b", "R5d", "R6b", "R6d", "R10b", "R10d", "R12b", "R12d")
    )
    q_endotherm = sum(
        heats[k]["Q_release_MW"] for k in ("R2", "R3", "R4", "R7", "R11b", "R11d") if heats[k]["Q_release_MW"] < 0
    )
    q_total = sum(v["Q_release_MW"] for v in heats.values())
    y = cell._mole_fractions("combined")
    nin_o2 = float(
        (cell.N_b_in + cell.N_d_in + cell.N_zu_b + cell.N_zu_d + cell.N_rez_b + cell.N_rez_d)[idx["O2"]]
    )
    nout_o2 = float((cell.N_b + cell.N_d)[idx["O2"]])
    row = {
        "bed_rank": bed_rank,
        "T_K": float(cell.T),
        "y_O2": float(y[idx["O2"]]),
        "y_CO": float(y[idx["CO"]]),
        "y_CO2": float(y[idx["CO2"]]),
        "y_CH4": float(y[idx["CH4"]]),
        "y_H2": float(y[idx["H2"]]),
        "dN_O2_mol_s": nin_o2 - nout_o2,
        "lf_o2_b": float(bundle.limit_factor_o2_bubble),
        "lf_o2_d": float(bundle.limit_factor_o2_dense),
        "vm": _vm_fuels(cell),
        "reactions": heats,
        "Q_total_release_MW": q_total,
        "Q_oxidation_MW": q_ox,
        "Q_endothermic_MW": q_endotherm,
        "top_reactions": [
            {"name": n, "extent": v["extent_mol_s"], "Q_MW": v["Q_release_MW"]}
            for n, v in ranked[:8]
            if abs(v["Q_release_MW"]) > 1e-6 or abs(v["extent_mol_s"]) > 1e-6
        ],
        "enthalpy": _enthalpy_parts(cell),
    }
    return row


def _print_cell(row: dict[str, Any]) -> None:
    e = row["enthalpy"]
    print(
        f"\nbed{row['bed_rank']}: T={row['T_K']:.0f} K  yO2={row['y_O2']:.3f} "
        f"dO2={row['dN_O2_mol_s']:.2f} mol/s  lf_d={row['lf_o2_d']:.3f}  "
        f"vm_comb={row['vm']['combustible']:.3f}",
        flush=True,
    )
    print(
        f"  Q_rxn total={row['Q_total_release_MW']:+.3f} MW  "
        f"(ox={row['Q_oxidation_MW']:+.3f}, endo={row['Q_endothermic_MW']:+.3f})",
        flush=True,
    )
    for t in row["top_reactions"][:6]:
        print(
            f"    {t['name']:5s}  ext={t['extent']:+8.3f} mol/s  Q={t['Q_MW']:+7.3f} MW",
            flush=True,
        )
    print(
        f"  enthalpy: Tin={e['Tin_gas_K']:.0f} Tzu_s={e['T_zu_solid_K']:.0f}  "
        f"Hin={e['Hin_MW']:+.2f} Hout_g={e['Hout_g_MW']:+.2f} Hout_s={e['Hout_s_MW']:+.2f} "
        f"Qwall={e['Q_wall_MW']:.3f} E={e['E_residual_MW']:+.3f} MW",
        flush=True,
    )
    sinks = sorted(
        (
            ("solid_zu", e["solid_zu_MW"]),
            ("gas_zu", e["gas_zu_MW"]),
            ("gas_axial_in", e["gas_axial_in_MW"]),
            ("solid_auf_in", e["solid_auf_in_MW"]),
            ("solid_ab_in", e["solid_ab_in_MW"]),
            ("solid_in", e["solid_in_MW"]),
        ),
        key=lambda kv: abs(kv[1]),
        reverse=True,
    )
    print(
        "  H_in parts (MW): "
        + ", ".join(f"{k}={v:+.2f}" for k, v in sinks if abs(v) > 0.01),
        flush=True,
    )
    # 去向粗账：放热 vs 把冷流加热到 T 的需求代理
    print(
        f"  对照: Q_rxn={row['Q_total_release_MW']:+.3f} vs "
        f"|solid_zu|={abs(e['solid_zu_MW']):.2f} |gas_zu|={abs(e['gas_zu_MW']):.2f} "
        f"|Hout_s−相关| 见上",
        flush=True,
    )


def run_solve(
    *,
    profile: str,
    T0_seed: float | None,
    max_iter: int,
) -> tuple[Any, dict[str, Any], float]:
    _orig = ips._install_temperature_fence_from_vorabrechnung

    def _install(reactor, T_est):
        _orig(reactor, T_est)
        if T0_seed is None:
            return
        beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if not beds:
            return
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
    add_json_output_arg(p, "phase2_vm_zone_reaction_heat_audit.json")
    p.add_argument("--profile", default="front_loaded", choices=("front_loaded", "linear"))
    p.add_argument(
        "--seed-T0",
        type=float,
        default=1250.0,
        help="bed0 软种子；负值表示不强制种子",
    )
    p.add_argument("--n-zone-cells", type=int, default=3)
    args = p.parse_args()
    seed = None if float(args.seed_T0) < 0 else float(args.seed_T0)

    print(
        f"=== solve profile={args.profile} seed={seed} iter={args.max_global_iter} ===",
        flush=True,
    )
    reactor, result, wall = run_solve(
        profile=str(args.profile), T0_seed=seed, max_iter=int(args.max_global_iter)
    )
    beds = bed_cell_indices(reactor)
    n_zone = int(args.n_zone_cells)
    cells_out = []
    for rank, bi in enumerate(beds[:n_zone]):
        row = analyze_cell(reactor.cells[bi], rank)
        cells_out.append(row)
        _print_cell(row)

    q_zone = sum(c["Q_total_release_MW"] for c in cells_out)
    q_ox_zone = sum(c["Q_oxidation_MW"] for c in cells_out)
    solid_zu0 = cells_out[0]["enthalpy"]["solid_zu_MW"] if cells_out else 0.0
    gas_zu0 = cells_out[0]["enthalpy"]["gas_zu_MW"] if cells_out else 0.0
    print(
        f"\nZONE SUMMARY: Q_rxn={q_zone:+.3f} MW (ox={q_ox_zone:+.3f})  "
        f"bed0 cold sinks solid_zu={solid_zu0:+.2f} gas_zu={gas_zu0:+.2f} MW  "
        f"ratio ox/|cold|={q_ox_zone / max(abs(solid_zu0) + abs(gas_zu0), 1e-9):.3f}",
        flush=True,
    )

    payload = {
        "note": (
            "VM-zone reaction extents + heat release vs Eq.2.7 enthalpy destinations; "
            "Q_release>0 is exothermic"
        ),
        "profile": args.profile,
        "T0_seed_K": seed,
        "wall_s": wall,
        "converged_fully": bool(result.get("converged_fully")),
        "rms_scaled_final": float(result.get("rms_scaled_final") or float("nan")),
        "exit_dry": {k: float(v) for k, v in (result.get("exit_gas_dry") or {}).items()},
        "zone_Q_rxn_MW": q_zone,
        "zone_Q_oxidation_MW": q_ox_zone,
        "bed0_cold_sink_MW": abs(solid_zu0) + abs(gas_zu0),
        "cells": cells_out,
    }
    write_json(args.output, payload)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
