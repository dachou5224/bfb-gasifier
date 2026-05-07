#!/usr/bin/env python3
"""LU 工况下 CH4 / H2O 来源与去向审计。

固定 shared ``global_nr`` 口径，逐 cell 拆分以下对 `CH4` / `H2O` 的净贡献。
本文档与输出默认使用 **Hamel (1999) thesis 编号**；括号中的 `impl rN`
表示当前源码里的 legacy implementation label：

- `pyrolysis` / `drying`
- `R2` char-steam gasification
- `R4` hydrogasification（`impl r3`）
- `R7` CH4 oxidation（`impl r6`）
- `R9` methane reforming（`impl r7`）
- `R8` WGSR
- `R11` tar reforming / cracking
- `R6` H2 oxidation（`impl r12`）

所有反应量均按当前 `build_reaction_sources()` 的 limiter 口径重算，确保与主模型一致。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.kinetics.char_reactions import (
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)
from src.kinetics.gas_reactions import (
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R6,
    rate_R7,
    rate_R8,
    rate_R9,
    rate_R12,
)
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


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
    cell.calc_hydrodynamics()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau_val)
    cell.calc_exchange()

    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    gas_src_vm = cell._vm_gas_source_cache.copy()
    areas = cell._calc_char_surface_area_per_class()
    total_area = float(np.sum(areas))
    rho_cat = cell._catalyst_bulk_density()
    char_conversion = cell._compute_char_conversion()
    D_g = __import__("src.core.species", fromlist=["gas_diffusivity_correlation"]).gas_diffusivity_correlation(T, P)

    r5b = cell.r5_scale * rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b
    r6b = cell.r6_scale * rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]]) * cell.V_b
    r12b = rate_R12(T, C_b[idx["H2"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b if cell.enable_r12 else 0.0

    r5d = cell.r5_scale * rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * cell.V_d
    r6d = cell.r6_scale * rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]]) * cell.V_d
    r12d = rate_R12(T, C_d[idx["H2"]], C_d[idx["O2"]], P, y_d, idx) * cell.V_d if cell.enable_r12 else 0.0
    ext7 = cell.r7_scale * rate_R7(T, C_d[idx["CH4"]], C_d[idx["H2O"]], C_d[idx["CO"]], C_d[idx["H2"]]) * cell.V_d
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

    r1 = r2 = r3 = r4 = 0.0
    alpha = 1.0
    if total_area > 0.0:
        dc = d_core_from_spm_char_conversion(float(char_conversion), float(cell.solid.d_p))
        r1, alpha = rate_R1(T, C_d[idx["O2"]], float(cell.solid.d_p), D_g, dc, cell.fuel_type)
        r2 = rate_R2(T, C_d[idx["H2O"]], float(cell.solid.d_p), D_g, dc)
        r3 = rate_R3(T, C_d[idx["H2"]], float(cell.solid.d_p), D_g, dc)
        p_co2 = max(float(C_d[idx["CO2"]] * Rg * T), 0.0)
        p_co = max(float(C_d[idx["CO"]] * Rg * T), 0.0)
        r4 = cell.r4_scale * rate_R4_effective(T, p_co2, p_co, float(cell.solid.d_p), D_g, dc)

    n_ex = cell.N_ex
    o2_supply_bubble = max(
        float(cell.N_zu_b[idx["O2"]]) + float(cell.N_b_in[idx["O2"]]) + max(-float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_b[idx["O2"]]),
        1e-6,
    )
    o2_supply_dense = max(
        float(cell.N_zu_d[idx["O2"]]) + float(cell.N_d_in[idx["O2"]]) + max(float(n_ex[idx["O2"]]), 0.0) + float(cell.N_rez_d[idx["O2"]]),
        1e-6,
    )
    bubble_o2_demand = r5b + 1.5 * r6b + 0.5 * r12b + nu_o2_r10 * ext10b
    dense_o2_demand = r5d + 1.5 * r6d + 0.5 * r12d + alpha * r1 * total_area + 1.5 * ext9 + nu_o2_r10 * ext10d
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
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    rows = [_analyze_cell(cell) for cell in reactor.cells]
    print("=" * 180)
    print("LU CH4/H2O source-path audit (shared global NR)")
    print("=" * 180)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"init={result.get('nr_init_strategy')} jacobian={result.get('nr_jacobian_strategy')} "
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("labels: thesis 编号优先；R4i3=加氢气化，R7i6=CH4 oxidation，R9i7=methane reforming，R6i12=H2 oxidation")
    print("-" * 180)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'yCH4':>8} {'yH2O':>8} "
        f"{'CH4_py':>9} {'CH4_R4i3':>10} {'CH4_R7i6':>10} {'CH4_R9i7':>10} {'CH4_R11':>10} {'CH4_net':>10} "
        f"{'H2O_py':>9} {'H2O_R2':>9} {'H2O_R7i6':>10} {'H2O_R9i7':>10} {'H2O_R8':>9} {'H2O_R11':>10} {'H2O_R6i12':>11} {'H2O_net':>10}"
    )
    print("-" * 180)
    for i, row in enumerate(rows):
        ch4 = row["ch4_terms"]
        h2o = row["h2o_terms"]
        ch4_r6 = ch4["R6b"] + ch4["R6d"]
        ch4_r11 = ch4["R11b"] + ch4["R11d"]
        h2o_r6 = h2o["R6b"] + h2o["R6d"]
        h2o_r11 = h2o["R11b"] + h2o["R11d"]
        h2o_r12 = h2o["R12b"] + h2o["R12d"]
        ch4_net = sum(ch4.values())
        h2o_net = sum(h2o.values())
        print(
            f"{i:>4d} {float(reactor.cells[i].geo.h_center/cfg.H_bed):>6.2f} {row['T']:>8.1f} "
            f"{row['wet']['CH4']:>8.4f} {row['wet']['H2O']:>8.4f} "
            f"{ch4['pyrolysis']:>9.3f} {ch4['R3']:>10.3f} {ch4_r6:>10.3f} {ch4['R7']:>10.3f} {ch4_r11:>10.3f} {ch4_net:>10.3f} "
            f"{h2o['pyrolysis']:>9.3f} {h2o['R2']:>9.3f} {h2o_r6:>10.3f} {h2o['R7']:>10.3f} {h2o['R8']:>9.3f} {h2o_r11:>10.3f} {h2o_r12:>11.3f} {h2o_net:>10.3f}"
        )
        print(
            f"     raw: R9raw_i7={row['raw']['ext7']:.3f} R8raw={row['raw']['ext8']:.3f} "
            f"R2_area={row['raw']['r2_area']:.3f} R4_area_i3={row['raw']['r3_area']:.3f} "
            f"ext11_sum={row['raw']['ext11_sum']:.3f} gas_src(CH4,H2O)=({row['raw']['gas_src_ch4']:.3f},{row['raw']['gas_src_h2o']:.3f}) "
            f"limiters(H2O={row['limit_factor_h2o']:.3f}, O2b={row['limit_factor_o2_bubble']:.3f}, O2d={row['limit_factor_o2_dense']:.3f})"
        )
    print("-" * 180)
    print("Notes:")
    print("  - CH4_net / H2O_net 为按 build_reaction_sources 同口径重算的局部净源项近似。")
    print("  - 若 CH4 几乎全靠 R11 生成而 pyrolysis 很弱，应优先审 VM->CH4 分配。")
    print("  - 若 H2O 赤字主要由 thesis R2 / R8 / R9(impl r7) 驱动，则优先审蒸汽相关反应链而不是热峰。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
