"""单 cell 氧反应与收敛诊断。

目标：
1. 始终先在单 cell 级别检查是否可收敛；
2. 若单 cell 收敛较好，再看该 cell 的 O2 相关反应主导项；
3. 判断该 cell 是否适合作为整炉传播的稳定基元。

默认诊断：
- 稳定候选分支：cells 0-5
- 参考温峰分支：cell 0

用法：
    .venv/bin/python scripts/audit_single_cell_oxygen_convergence.py
"""

from __future__ import annotations

import copy
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import src.kinetics.char_reactions as cr
import src.kinetics.gas_reactions as gr
from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as IDX
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension, rate_R6, rate_R9, rate_R12
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10
from src.solvers.cell_solver import solve_cell
from tests.validation_case_utils import build_phase1_htw_lu_refined_config, load_case_LU


@dataclass
class Scenario:
    name: str
    dense: float
    heat_loss: float
    r5s: float = 1.0
    r6s: float = 1.0
    r4s: float = 1.0


@contextmanager
def kinetics_scaling(r5_scale: float, r6_scale: float, r4_scale: float):
    b_r5 = gr.R5_suspension_k_T05
    b_r6 = gr.R6_k0
    b_r4 = cr.R4_kf_k0
    try:
        gr.R5_suspension_k_T05 = b_r5 * r5_scale
        gr.R6_k0 = b_r6 * r6_scale
        cr.R4_kf_k0 = b_r4 * r4_scale
        yield
    finally:
        gr.R5_suspension_k_T05 = b_r5
        gr.R6_k0 = b_r6
        cr.R4_kf_k0 = b_r4


def _build_reactor(sc: Scenario) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)
    cfg.gas_inlet_dense_frac = sc.dense
    cfg.heat_loss_frac = sc.heat_loss
    return Reactor(cfg)


def _prepare_cell(sc: Scenario, cell_index: int):
    """按单 cell 顺序把入口传播到目标 cell，并求解其下方 cells。"""
    r = _build_reactor(sc)
    r._set_bottom_cell_feeds()

    for i in range(cell_index + 1):
        r._propagate_upstream(i)
        c = r.cells[i]
        if np.sum(np.maximum(c.N_d + c.N_b, 0.0)) < 1e-9:
            c.N_d[:] = np.maximum(c.N_d_in + c.N_zu_d + c.N_rez_d, 0.0)
            c.N_b[:] = np.maximum(c.N_b_in + c.N_zu_b + c.N_rez_b, 0.0)
            c.m_solid[:] = np.maximum(c.m_solid_in + c.m_solid_zu + c.m_solid_rez, 0.0)

        c.calc_hydrodynamics()
        c.compute_vorabrechnung(c.geo.dh / max(c.u_mf, 1e-3))

        for _ in range(2):
            c.calc_exchange()
            c.calc_reactions()
            mig = c._calc_size_migration()
            c.m_solid[:] = np.maximum(c.m_solid_in + c.m_solid_zu + c.m_solid_rez + c.R_solid + mig, 0.0)
            c.calc_hydrodynamics()

        if i < cell_index:
            solve_cell(c, stiff_stabilization=True, verbose=False)
            c.calc_exchange()
            c.calc_reactions()
            mig = c._calc_size_migration()
            c.m_solid[:] = np.maximum(c.m_solid_in + c.m_solid_zu + c.m_solid_rez + c.R_solid + mig, 0.0)

    return r.cells[cell_index]


def _oxygen_breakdown(cell) -> dict:
    idx = IDX
    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    from src.core.species import gas_diffusivity_correlation

    D_g = gas_diffusivity_correlation(T, P)

    r5b = rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b
    r6b = rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]]) * cell.V_b
    r12b = rate_R12(T, C_b[idx["H2"]], C_b[idx["O2"]], P, y_b, idx) * cell.V_b if cell.enable_r12 else 0.0
    r5d = rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * cell.V_d
    r6d = rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]]) * cell.V_d
    r12d = rate_R12(T, C_d[idx["H2"]], C_d[idx["O2"]], P, y_d, idx) * cell.V_d if cell.enable_r12 else 0.0
    ext9 = 0.0 if cell.use_gibbs_minor else rate_R9(T, C_d[idx["H2S"]], C_d[idx["O2"]]) * cell.V_d

    C_tar_b = max(C_b[idx["TAR1"]] + C_b[idx["TAR2"]], 0.0)
    C_tar_d = max(C_d[idx["TAR1"]] + C_d[idx["TAR2"]], 0.0)
    ext10b = rate_R10(T, C_tar_b, C_b[idx["O2"]], P, cell.fuel_type) * cell.V_b
    ext10d = rate_R10(T, C_tar_d, C_d[idx["O2"]], P, cell.fuel_type) * cell.V_d
    stoich_r10 = get_lumped_tar_stoichiometry("R10", cell.fuel_type)
    nu_o2_r10 = max(0.0, -float(stoich_r10.get("O2", 0.0)))

    areas = cell._calc_char_surface_area_per_class()
    total_area = float(np.sum(areas))
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    alpha = 1.0
    if total_area > 0:
        X = cell._compute_char_conversion()
        dc = d_core_from_spm_char_conversion(X, cell.solid.d_p)
        r1, alpha = rate_R1(T, C_d[idx["O2"]], cell.solid.d_p, D_g, dc, cell.fuel_type)
        r2 = rate_R2(T, C_d[idx["H2O"]], cell.solid.d_p, D_g, dc)
        r3 = rate_R3(T, C_d[idx["H2"]], cell.solid.d_p, D_g, dc)
        p_co2 = max(C_d[idx["CO2"]] * Rg * T, 0.0)
        p_co = max(C_d[idx["CO"]] * Rg * T, 0.0)
        r4 = rate_R4_effective(T, p_co2, p_co, cell.solid.d_p, D_g, dc)

    gas_src_vm = np.maximum(cell._vm_gas_source_cache, 0.0)
    vm_co = float(gas_src_vm[idx["CO"]])
    vm_ch4 = float(gas_src_vm[idx["CH4"]])
    vm_h2 = float(gas_src_vm[idx["H2"]])
    vm_h2s = float(gas_src_vm[idx["H2S"]])
    vm_tar = float(gas_src_vm[idx["TAR1"]] + gas_src_vm[idx["TAR2"]])
    vm_o2eq = 0.5 * vm_co + 1.5 * vm_ch4 + 0.5 * vm_h2 + 1.5 * vm_h2s + nu_o2_r10 * vm_tar

    terms = {
        "R10_tar_oxid": float(nu_o2_r10 * (ext10b + ext10d)),
        "R7i6_dense_CH4_oxid": float(1.5 * r6d),
        "R5b_CO_oxid": float(r5b),
        "R7i6_bubble_CH4_oxid": float(1.5 * r6b),
        "R6i12_bubble_H2_oxid": float(0.5 * r12b),
        "R1_char_comb": float(alpha * r1 * total_area),
        "R5d_CO_oxid": float(r5d),
        "R6i12_H2_oxid": float(0.5 * r12d),
        "implR9_H2S_placeholder": float(1.5 * ext9),
    }
    top_terms = sorted(terms.items(), key=lambda kv: abs(kv[1]), reverse=True)

    o2_supply = max(
        cell.N_zu_d[idx["O2"]] + cell.N_d_in[idx["O2"]]
        + cell.N_zu_b[idx["O2"]] + cell.N_b_in[idx["O2"]]
        + cell.N_rez_d[idx["O2"]] + cell.N_rez_b[idx["O2"]],
        1e-12,
    )
    homogeneous = (
        terms["R10_tar_oxid"]
        + terms["R7i6_dense_CH4_oxid"]
        + terms["R5b_CO_oxid"]
        + terms["R7i6_bubble_CH4_oxid"]
        + terms["R6i12_bubble_H2_oxid"]
        + terms["R5d_CO_oxid"]
        + terms["R6i12_H2_oxid"]
        + terms["implR9_H2S_placeholder"]
    )

    return {
        "o2_supply": float(o2_supply),
        "vm_o2eq": float(vm_o2eq),
        "homogeneous_o2": float(homogeneous),
        "r2_proxy": float(r2 * total_area),
        "r4_proxy": float(r4 * total_area),
        "top_terms": top_terms,
    }


def _run_single_cell(cell, stiff: bool):
    c = copy.deepcopy(cell)
    res = solve_cell(c, stiff_stabilization=stiff, verbose=False)
    c.calc_hydrodynamics()
    c.compute_vorabrechnung(c.geo.dh / max(c.u_mf, 1e-3))
    c.calc_exchange()
    diag = _oxygen_breakdown(c)
    return res, c, diag


def print_case(sc: Scenario, cell_indices: list[int]) -> None:
    print("\n" + "=" * 120)
    print(
        f"Scenario: {sc.name} | dense={sc.dense:.2f} heat_loss={sc.heat_loss:.2f} "
        f"impl_scales[R5={sc.r5s:.2f}, R6(thesis R7)={sc.r6s:.2f}, R4(thesis R3)={sc.r4s:.2f}]"
    )
    print("=" * 120)
    print("labels: thesis 编号优先；R7i6=CH4 oxidation，R6i12=H2 oxidation，implR9=sulfur placeholder")

    with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
        for i in cell_indices:
            base_cell = _prepare_cell(sc, i)
            res0, c0, d0 = _run_single_cell(base_cell, stiff=False)
            res1, c1, d1 = _run_single_cell(base_cell, stiff=True)

            print(f"\n[cell {i}] dh={base_cell.geo.dh:.3f} m  h_center={base_cell.geo.h_center:.3f} m")
            print(
                f"  default : converged={res0.get('converged')} phys={res0.get('physically_converged')} rms={res0.get('rms_scaled'):.3e} res={res0.get('residual'):.3e} T={c0.T:.1f}"
            )
            print(
                f"  stiff   : converged={res1.get('converged')} phys={res1.get('physically_converged')} rms={res1.get('rms_scaled'):.3e} res={res1.get('residual'):.3e} T={c1.T:.1f} attempted={res1.get('stiff_attempted')} accepted={res1.get('stiff_accepted')}"
            )
            print(
                f"  O2 supply={d1['o2_supply']:.3f}  VM_O2eq={d1['vm_o2eq']:.3f}  homogeneous_O2={d1['homogeneous_o2']:.3f}  R2proxy={d1['r2_proxy']:.3f}  R3proxy={d1['r4_proxy']:.3f}"
            )
            print("  dominant O2-related terms:")
            for name, val in d1["top_terms"][:5]:
                print(f"    - {name}: {val:.3f}")

            if bool(res1.get("physically_converged")):
                print("  judgement: 单 cell 物理收敛尚可，可继续用于传播诊断。")
            else:
                print("  judgement: 单 cell 物理收敛仍不足，先不要把该 cell 的问题归咎于整炉传播。")


def main() -> int:
    stable = Scenario(name="stable_window_candidate", dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)
    peaky = Scenario(name="peaky_branch_reference", dense=0.35, heat_loss=0.12, r5s=1.00, r6s=1.00, r4s=1.00)

    print_case(stable, [0, 1, 2, 3, 4, 5])
    print_case(peaky, [0, 1, 2])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
