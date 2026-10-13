#!/usr/bin/env python3
"""LU 工况下 bed-exit -> reactor-exit 反应预算审计。"""

from __future__ import annotations

import numpy as np

from _lu_harness import (
    build_phase2_config,
    load_case_lu,
    phase1_global_nr_solve_kwargs,
)
from _audit_common import compute_char_reaction_rates, prepare_cell_for_analysis
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.kinetics.gas_reactions import (
    rate_R12,
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R6,
    rate_R7,
    rate_R8,
    rate_R9,
)
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor

TRACK_RXNS_THESIS = ("R5", "R6", "R7", "R8", "R9", "R10", "R11")
TRACK_SPECIES = ("CO", "CO2", "H2", "H2O", "CH4", "O2")


def _term_dict() -> dict[str, float]:
    return {
        "pyrolysis": 0.0,
        "R2": 0.0,
        "R3": 0.0,
        "R6b": 0.0,
        "R6d": 0.0,
        "R7": 0.0,
        "R8": 0.0,
        "R11b": 0.0,
        "R11d": 0.0,
        "R12b": 0.0,
        "R12d": 0.0,
    }


def _analyze_cell(cell) -> dict[str, object]:
    """逐 cell 重算 CH4/H2O 来源去向（原 audit_ch4_h2o_paths_lu.py 内联）。"""
    prepare_cell_for_analysis(cell)

    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    gas_src_vm = cell._vm_gas_source_cache.copy()
    rho_cat = cell._catalyst_bulk_density()
    D_g = __import__("src.core.species", fromlist=["gas_diffusivity_correlation"]).gas_diffusivity_correlation(T, P)

    r5b = cell.r5_scale * rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b
    r6b = cell.r6_scale * rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]]) * cell.V_b
    r12b = rate_R12(T, C_b[idx["H2"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b if cell.enable_r12 else 0.0

    r5d = cell.r5_scale * rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * cell.V_d
    r6d = cell.r6_scale * rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]]) * cell.V_d
    r12d = rate_R12(T, C_d[idx["H2"]], C_d[idx["O2"]], P, y_d, idx) * cell.V_d if cell.enable_r12 else 0.0
    ext7 = cell.r7_scale * rate_R7(
        T, C_d[idx["CH4"]], C_d[idx["H2O"]], C_d[idx["CO"]], C_d[idx["H2"]],
        rho_ash_kg_m3=cell._ash_bulk_density(),
    ) * cell.V_d
    ext8 = rate_R8(T, P, float(y_d[idx["CO"]]), float(y_d[idx["H2O"]]), float(y_d[idx["CO2"]]), float(y_d[idx["H2"]])) * cell.V_d
    ext9 = 0.0 if cell.use_gibbs_minor else rate_R9(T, C_d[idx["H2S"]], C_d[idx["O2"]]) * cell.V_d

    C_tar_b = max(float(C_b[idx["TAR1"]] + C_b[idx["TAR2"]]), 0.0)
    C_tar_d = max(float(C_d[idx["TAR1"]] + C_d[idx["TAR2"]]), 0.0)
    ext10b = rate_R10(T, C_tar_b, C_b[idx["O2"]], P, cell.fuel_type) * cell.V_b
    ext10d = rate_R10(T, C_tar_d, C_d[idx["O2"]], P, cell.fuel_type) * cell.V_d
    ext11b = rate_R11_bubble(T, C_tar_b) * cell.V_b
    ext11d = rate_R11_suspension(T, C_tar_d, rho_cat) * cell.V_d

    stoich_r10 = get_lumped_tar_stoichiometry("R10", cell.fuel_type)
    stoich_r11 = get_lumped_tar_stoichiometry("R11", cell.fuel_type)
    nu_o2_r10 = max(0.0, -float(stoich_r10.get("O2", 0.0)))
    nu_h2o_r11 = max(0.0, -float(stoich_r11.get("H2O", 0.0)))

    char_rates = compute_char_reaction_rates(cell, D_g)
    r1 = char_rates["r1"]
    r2 = char_rates["r2"]
    r3 = char_rates["r3"]
    r4 = char_rates["r4"]
    alpha = char_rates["alpha"]
    total_area = char_rates["total_area"]

    n_ex = cell.N_ex
    o2_supply_bubble = max(
        float(cell.N_zu_b[idx["O2"]]) + float(cell.N_b_in[idx["O2"]]) + max(-float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_b[idx["O2"]]),
        1e-6,
    )
    o2_supply_dense = max(
        float(cell.N_zu_d[idx["O2"]]) + float(cell.N_d_in[idx["O2"]]) + max(float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_d[idx["O2"]]),
        1e-6,
    )
    # R5: CO + 1/2 O2 → CO2；extent = mol CO/s → O₂ demand = 0.5·r5
    bubble_o2_demand = 0.5 * r5b + 1.5 * r6b + 0.5 * r12b + nu_o2_r10 * ext10b
    dense_o2_demand = 0.5 * r5d + 1.5 * r6d + 0.5 * r12d + alpha * r1 * total_area + 1.5 * ext9 + nu_o2_r10 * ext10d
    limit_factor_o2_bubble = min(1.0, (0.995 * o2_supply_bubble) / max(bubble_o2_demand, 1e-9))
    limit_factor_o2_dense = min(1.0, (0.995 * o2_supply_dense) / max(dense_o2_demand, 1e-9))

    r6b_lim = r6b * limit_factor_o2_bubble
    r12b_lim = r12b * limit_factor_o2_bubble
    r6d_lim = r6d * limit_factor_o2_dense
    r12d_lim = r12d * limit_factor_o2_dense

    h2o_supply = max(
        float(cell.N_zu_d[idx["H2O"]]) + float(cell.N_d_in[idx["H2O"]]) + float(cell.N_zu_b[idx["H2O"]]) + float(cell.N_b_in[idx["H2O"]])
        + float(cell.N_rez_d[idx["H2O"]]) + float(cell.N_rez_b[idx["H2O"]]) + float(gas_src_vm[idx["H2O"]]),
        1e-6,
    )
    h2o_consume = r2 * total_area + ext7 + max(0.0, ext8) + nu_h2o_r11 * (ext11b + ext11d)
    limit_factor_h2o = min(1.0, (0.995 * h2o_supply) / max(h2o_consume, 1e-9))

    r2_lim = r2 * limit_factor_h2o
    r3_lim = r3
    ext7_lim = ext7 * limit_factor_h2o
    ext8_lim = ext8 * limit_factor_h2o if ext8 > 0.0 else ext8
    ext11b_lim = ext11b * limit_factor_h2o
    ext11d_lim = ext11d * limit_factor_h2o

    ch4 = _term_dict()
    h2o = _term_dict()
    ch4["pyrolysis"] += float(gas_src_vm[idx["CH4"]])
    h2o["pyrolysis"] += float(gas_src_vm[idx["H2O"]])

    ch4["R3"] += float(r3_lim * total_area)
    h2o["R2"] -= float(r2_lim * total_area)
    ch4["R6b"] -= float(r6b_lim)
    h2o["R6b"] += float(2.0 * r6b_lim)
    ch4["R6d"] -= float(r6d_lim)
    h2o["R6d"] += float(2.0 * r6d_lim)
    ch4["R7"] -= float(ext7_lim)
    h2o["R7"] -= float(ext7_lim)
    h2o["R8"] -= float(ext8_lim)
    ch4["R11b"] += float(ext11b_lim * float(stoich_r11.get("CH4", 0.0)))
    h2o["R11b"] += float(ext11b_lim * float(stoich_r11.get("H2O", 0.0)))
    ch4["R11d"] += float(ext11d_lim * float(stoich_r11.get("CH4", 0.0)))
    h2o["R11d"] += float(ext11d_lim * float(stoich_r11.get("H2O", 0.0)))
    h2o["R12b"] += float(r12b_lim)
    h2o["R12d"] += float(r12d_lim)

    return {
        "xi": float(cell.geo.h_center),
        "T": float(T),
        "wet": {
            "CH4": float(cell._mole_fractions("combined")[idx["CH4"]]),
            "H2O": float(cell._mole_fractions("combined")[idx["H2O"]]),
            "CO": float(cell._mole_fractions("combined")[idx["CO"]]),
            "CO2": float(cell._mole_fractions("combined")[idx["CO2"]]),
            "H2": float(cell._mole_fractions("combined")[idx["H2"]]),
        },
        "ch4_terms": ch4,
        "h2o_terms": h2o,
        "limit_factor_h2o": float(limit_factor_h2o),
        "limit_factor_o2_bubble": float(limit_factor_o2_bubble),
        "limit_factor_o2_dense": float(limit_factor_o2_dense),
        "raw": {
            "ext7": float(ext7),
            "ext8": float(ext8),
            "r2_area": float(r2 * total_area),
            "r3_area": float(r3 * total_area),
            "ext11_sum": float(ext11b + ext11d),
            "gas_src_ch4": float(gas_src_vm[idx["CH4"]]),
            "gas_src_h2o": float(gas_src_vm[idx["H2O"]]),
            "r4_area": float(r4 * total_area),
        },
    }


def main() -> int:
    case = load_case_lu()
    cfg = build_phase2_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        phase1_global_nr_solve_kwargs(),
        check_x0=True,
    )

    bed_exit_wet = result["bed_exit_gas"]
    bed_exit_dry = result["bed_exit_gas_dry"]
    exit_wet = result["exit_gas"]
    exit_dry = result["exit_gas_dry"]
    fb_diag_impl = result.get("freeboard_reaction_diag_impl", result["freeboard_reaction_diag"])
    fb_diag = result.get("freeboard_reaction_diag_thesis", fb_diag_impl)
    top_cells = []
    for i in range(max(len(reactor.cells) - 3, 0), len(reactor.cells)):
        row = _analyze_cell(reactor.cells[i])
        top_cells.append((i, row))

    print("=" * 176)
    print("LU bed-exit -> reactor-exit budget audit")
    print("=" * 176)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"secondary_agent={case['secondary_injection_agent']} "
        f"secondary_Nm3_h={case['secondary_injection_Nm3_h']:.1f} "
        f"xi_inj={case['secondary_injection_xi']:.2f} "
        f"Texit={result['reactor_exit_T']:.1f}K bed_exit_T={result['bed_T_profile'][-1]:.1f}K "
        f"n_iter={result['n_iter']} rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    numbering = result.get("reaction_numbering", {})
    if numbering:
        print(
            f"reaction_label_authority={numbering.get('authority', 'unknown')} | "
            "table=freeboard_reaction_diag_thesis (impl view in *_impl)"
        )
    print("-" * 176)
    print(f"{'species':>8} {'bed_exit_wet':>12} {'reactor_exit_wet':>18} {'delta_wet':>12} {'bed_dry':>12} {'exit_dry':>12}")
    print("-" * 176)
    for sp in TRACK_SPECIES:
        print(
            f"{sp:>8s} {float(bed_exit_wet.get(sp, 0.0)):>12.5f} "
            f"{float(exit_wet.get(sp, 0.0)):>18.5f} "
            f"{float(exit_wet.get(sp, 0.0) - bed_exit_wet.get(sp, 0.0)):>12.5f} "
            f"{float(bed_exit_dry.get(sp, 0.0)):>12.5f} {float(exit_dry.get(sp, 0.0)):>12.5f}"
        )

    print("-" * 176)
    print(f"{'freeboard rxn':>14} {'sum_extent':>14} {'max_slice':>14}")
    print("-" * 176)
    for rxn in TRACK_RXNS_THESIS:
        vals = [float(v) for v in fb_diag.get(rxn, [])]
        print(f"{rxn:>14s} {sum(vals):>14.6f} {(max(vals, key=abs) if vals else 0.0):>14.6f}")

    print("-" * 176)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'H2O':>8} {'CH4':>8} "
        f"{'R4_area':>10} {'R7':>10} {'R8':>10} {'CH4_net':>10} {'H2O_net':>10}"
    )
    print("-" * 176)
    for i, row in top_cells:
        ch4_net = sum(float(v) for v in row["ch4_terms"].values())
        h2o_net = sum(float(v) for v in row["h2o_terms"].values())
        print(
            f"{i:>4d} {float(reactor.cells[i].geo.h_center/cfg.H_bed):>6.2f} {row['T']:>8.1f} "
            f"{row['wet']['CO']:>8.4f} {row['wet']['CO2']:>8.4f} {row['wet']['H2']:>8.4f} "
            f"{row['wet']['H2O']:>8.4f} {row['wet']['CH4']:>8.4f} "
            f"{row['raw']['r4_area']:>10.3f} {row['raw']['ext7']:>10.3f} {row['raw']['ext8']:>10.3f} "
            f"{ch4_net:>10.3f} {h2o_net:>10.3f}"
        )

    print("-" * 176)
    print("readout:")
    print("  - 顶床层 3 个 cell 用来判断进入 freeboard 前，CO/H2/H2O/CH4 已经被床层 chemistry 推到什么方向。")
    print("  - freeboard reaction sum 用来判断 reactor-exit 主要是延续床顶组成，还是被 freeboard 显著再分配。")
    print("  - 若 freeboard Δ 远小于 bed-top composition bias，则后续优先级应回到 bed chemistry。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
