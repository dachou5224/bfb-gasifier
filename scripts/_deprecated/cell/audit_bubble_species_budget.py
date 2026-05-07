"""Audit bubble CO/O2 budgets on locally converged single-cell states.

Purpose
-------
Compare the full-reactor accepted state against a re-solved single-cell state
for the same propagated upstream conditions, then decompose bubble CO/O2
budgets on the locally converged branch.

This is intended to answer whether the observed bubble-phase anomaly is
primarily a local constitutive issue or a propagated/accepted-state branch gap.

Usage
-----
    python3 scripts/audit_bubble_species_budget.py
"""

from __future__ import annotations

import copy
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from audit_single_cell_oxygen_convergence import Scenario, _prepare_cell, kinetics_scaling
from src.core.cell import S_CHAR
from src.core.cell_kinetics import build_reaction_sources
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as IDX, gas_diffusivity_correlation
from src.kinetics.gas_reactions import rate_R12, rate_R5_bubble, rate_R6
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10, rate_R11_bubble
from src.solvers.cell_solver import evaluate_cell_state, solve_cell
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_refined_config, load_case_LU


def build_full_reactor(sc: Scenario) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)
    cfg.gas_inlet_dense_frac = sc.dense
    cfg.heat_loss_frac = sc.heat_loss
    reactor = Reactor(cfg)
    reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)
    return reactor


def solve_local_cell(sc: Scenario, cell_index: int):
    base_cell = _prepare_cell(sc, cell_index)
    cell = copy.deepcopy(base_cell)
    result = solve_cell(cell, stiff_stabilization=True, verbose=False)
    metrics = evaluate_cell_state(cell)
    return cell, result, metrics


def _reaction_decomposition(cell):
    idx = IDX
    cell.calc_hydrodynamics()
    cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
    cell.calc_exchange()

    c_b = cell._concentrations("b")
    c_d = cell._concentrations("d")
    y_b = cell._mole_fractions("b")
    y_d = cell._mole_fractions("d")

    bundle = build_reaction_sources(
        T=cell.T,
        P=cell.P,
        fuel_type=cell.fuel_type,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=c_b,
        C_d=c_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=cell._vm_gas_source_cache,
        solid_sink_vm=cell._vm_solid_sink_cache,
        areas=cell._calc_char_surface_area_per_class(),
        solid_d_p=cell.solid.d_p,
        D_g=gas_diffusivity_correlation(cell.T, cell.P),
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=None,
        r4_scale=cell.r4_scale,
        r5_scale=cell.r5_scale,
        r6_scale=cell.r6_scale,
        r7_scale=cell.r7_scale,
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

    r5b = cell.r5_scale * rate_R5_bubble(cell.T, c_b[idx["CO"]], c_b[idx["O2"]], cell.P, y_b, idx) * cell.V_b
    r6b = cell.r6_scale * rate_R6(cell.T, c_b[idx["CH4"]], c_b[idx["O2"]]) * cell.V_b
    r12b = rate_R12(cell.T, c_b[idx["H2"]], c_b[idx["O2"]], cell.P, y_b, idx) * cell.V_b if cell.enable_r12 else 0.0
    c_tar_b = max(float(c_b[idx["TAR1"]] + c_b[idx["TAR2"]]), 0.0)
    ext10b = rate_R10(cell.T, c_tar_b, c_b[idx["O2"]], cell.P, cell.fuel_type) * cell.V_b
    ext11b = rate_R11_bubble(cell.T, c_tar_b) * cell.V_b

    lim_o2 = bundle.limit_factor_o2_bubble
    lim_h2o = bundle.limit_factor_h2o
    r5b_lim = r5b * lim_o2
    r6b_lim = r6b * lim_o2
    r12b_lim = r12b * lim_o2
    ext10b_lim = ext10b * lim_o2
    ext11b_lim = ext11b * lim_h2o

    sto10 = get_lumped_tar_stoichiometry("R10", cell.fuel_type)
    sto11 = get_lumped_tar_stoichiometry("R11", cell.fuel_type)

    co_terms = {
        "R5b": -2.0 * r5b_lim,
        "R6b": r6b_lim,
        "R10b": float(sto10.get("CO", 0.0)) * ext10b_lim,
        "R11b": float(sto11.get("CO", 0.0)) * ext11b_lim,
    }
    o2_terms = {
        "R5b": -r5b_lim,
        "R6b": -1.5 * r6b_lim,
        "R12b": -0.5 * r12b_lim,
        "R10b": float(sto10.get("O2", 0.0)) * ext10b_lim,
    }
    return bundle, co_terms, o2_terms


def _print_species_budget(cell, species: str, reaction_terms: dict[str, float]) -> None:
    j = IDX[species]
    inflow = float(cell.N_zu_b[j] + cell.N_b_in[j] + cell.N_rez_b[j])
    reaction_total = float(sum(reaction_terms.values()))
    exchange_to_bubble = float(-cell.N_ex[j])
    outflow = float(cell.N_b[j])
    residual = inflow + reaction_total + exchange_to_bubble - outflow

    print(
        f"  {species}: inflow={inflow:.4f} reaction={reaction_total:.4f} "
        f"exchange_to_bubble={exchange_to_bubble:.4f} outflow={outflow:.4f} residual={residual:.3e}"
    )
    print(
        "    reaction terms: "
        + ", ".join(f"{name}={value:.4f}" for name, value in reaction_terms.items())
    )


def main() -> int:
    sc = Scenario(name="stable_window_candidate", dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)
    cell_indices = [1, 2]

    with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
        reactor = build_full_reactor(sc)

        print("=" * 140)
        print(
            f"Scenario: {sc.name} | dense={sc.dense:.2f} heat_loss={sc.heat_loss:.2f} "
            f"R5={sc.r5s:.2f} R6={sc.r6s:.2f} R4={sc.r4s:.2f}"
        )
        print("=" * 140)

        for i in cell_indices:
            full_cell = reactor.cells[i]
            full_metrics = evaluate_cell_state(full_cell)
            local_cell, local_result, local_metrics = solve_local_cell(sc, i)
            bundle, co_terms, o2_terms = _reaction_decomposition(local_cell)

            print(f"\n[cell {i}]")
            print(
                f"  full-reactor accepted: T={full_cell.T:.1f}K "
                f"Nb(CO)={full_cell.N_b[IDX['CO']]:.4f} Nb(O2)={full_cell.N_b[IDX['O2']]:.4f} "
                f"rms={full_metrics['rms_scaled']:.3e}"
            )
            print(
                f"  single-cell resolved: T={local_cell.T:.1f}K "
                f"Nb(CO)={local_cell.N_b[IDX['CO']]:.4f} Nb(O2)={local_cell.N_b[IDX['O2']]:.4f} "
                f"phys={local_result['physically_converged']} rms={local_metrics['rms_scaled']:.3e}"
            )
            print(
                f"  limiter: bubble_o2={bundle.limit_factor_o2_bubble:.4f} "
                f"dense_o2={bundle.limit_factor_o2_dense:.4f} h2o={bundle.limit_factor_h2o:.4f}"
            )
            _print_species_budget(local_cell, "CO", co_terms)
            _print_species_budget(local_cell, "O2", o2_terms)

        print("\nInterpretation:")
        print("  If the single-cell resolved branch differs strongly from the full-reactor accepted state,")
        print("  the anomaly is at least partly a propagated/accepted-state branch issue, not just local kinetics.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
