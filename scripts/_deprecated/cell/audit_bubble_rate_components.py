"""Audit bubble-phase source terms behind the R5b/R10b dominance switch.

Purpose
-------
Decompose bubble-phase R5 and R10 into:
- concentration terms
- kinetic prefactors
- Gibbs driving force for R5
- per-volume rates and per-cell rates

This is meant to explain why stable-window cell1 is R10b-dominant while
cells 2-4 become R5b-dominant.

Usage
-----
    python3 scripts/audit_bubble_rate_components.py
"""

from __future__ import annotations

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
from src.core.species import GAS_SPECIES_INDEX as IDX, calc_tar_surrogate_fractions
from src.kinetics.arrhenius import k_hobbs, k_standard
from src.kinetics.gas_reactions import R5_bubble_E_Rg, R5_bubble_k0, rate_R5_bubble
from src.kinetics.tar_reactions import R10_E_Rg_AROM, R10_E_Rg_OLEF, R10_k0_AROM, R10_k0_OLEF, rate_R10
from src.thermodynamics.equilibrium import calc_gibbs_driving_force, calc_reaction_quotient, get_K_eq
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_refined_config, load_case_LU


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


def run_scenario(sc: Scenario):
    with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
        case = load_case_LU()
        cfg = build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)
        cfg.gas_inlet_dense_frac = sc.dense
        cfg.heat_loss_frac = sc.heat_loss
        reactor = Reactor(cfg)
        result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)
    return reactor, result


def weighted_r10_prefactor(T: float, P: float, fuel_type: str) -> float:
    fractions = calc_tar_surrogate_fractions(fuel_type)
    p_bar = P / 100_000.0
    total = 0.0
    for surrogate, weight in fractions.items():
        if weight <= 0.0:
            continue
        aromatic = surrogate in ("C6H6", "C10H8")
        k0 = R10_k0_AROM if aromatic else R10_k0_OLEF
        e_rg = R10_E_Rg_AROM if aromatic else R10_E_Rg_OLEF
        total += weight * k_hobbs(k0, e_rg * Rg, max(T, 300.0)) * (p_bar ** 0.3)
    return total


def main() -> int:
    sc = Scenario(name="stable_window_candidate", dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)
    reactor, result = run_scenario(sc)
    idx = IDX
    focus_cells = [0, 1, 2, 3, 4]
    rows: list[dict[str, float]] = []

    print("=" * 180)
    print(
        f"Scenario: {sc.name} | dense={sc.dense:.2f} heat_loss={sc.heat_loss:.2f} "
        f"R5={sc.r5s:.2f} R6={sc.r6s:.2f} R4={sc.r4s:.2f} | "
        f"Texit={result['T_profile'][-1]:.1f}K Tpeak={np.max(result['T_profile']):.1f}K"
    )
    print("-" * 180)
    print(
        f"{'cell':>4} {'T[K]':>8} {'C_CO':>9} {'C_O2':>9} {'C_tar':>9} {'Q/K':>9} {'drv5':>8} "
        f"{'k5':>10} {'k10':>10} {'phi5':>10} {'phi10':>10} {'r5v':>10} {'r10v':>10} {'R5/R10':>8}"
    )

    for i in focus_cells:
        cell = reactor.cells[i]
        cell.calc_hydrodynamics()
        cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
        cell.calc_exchange()

        c_b = cell._concentrations("b")
        y_b = cell._mole_fractions("b")

        c_co = float(c_b[idx["CO"]])
        c_o2 = float(c_b[idx["O2"]])
        c_tar = float(c_b[idx["TAR1"]] + c_b[idx["TAR2"]])

        q_p = calc_reaction_quotient("R5", y_b, cell.P, idx)
        k_eq = get_K_eq("R5", cell.T, cell.P)
        q_over_k = q_p / max(k_eq, 1e-300)
        drive_r5 = calc_gibbs_driving_force("R5", cell.T, cell.P, y_b, idx, clamp_irreversible=True)

        k5 = k_standard(R5_bubble_k0, R5_bubble_E_Rg * Rg, cell.T)
        k10 = weighted_r10_prefactor(cell.T, cell.P, cell.fuel_type)

        phi5 = c_co * np.sqrt(max(c_o2, 0.0))
        phi10 = np.sqrt(max(c_tar, 0.0)) * c_o2

        r5v = rate_R5_bubble(cell.T, c_co, c_o2, cell.P, y_b, idx)
        r10v = rate_R10(cell.T, c_tar, c_o2, cell.P, cell.fuel_type)

        rows.append(
            {
                "cell": float(i),
                "c_co": c_co,
                "c_o2": c_o2,
                "c_tar": c_tar,
                "q_over_k": q_over_k,
                "drive_r5": drive_r5,
                "k5": k5,
                "k10": k10,
                "phi5": phi5,
                "phi10": phi10,
                "r5v": r5v,
                "r10v": r10v,
            }
        )
        print(
            f"{i:4d} {cell.T:8.1f} {c_co:9.4f} {c_o2:9.4f} {c_tar:9.4f} {q_over_k:9.2e} {drive_r5:8.3f} "
            f"{k5:10.2f} {k10:10.2f} {phi5:10.4f} {phi10:10.4f} {r5v:10.2f} {r10v:10.2f} {r5v / max(r10v, 1e-30):8.2f}"
        )

    print("-" * 180)
    row1 = rows[1]
    row2 = rows[2]
    row4 = rows[4]
    k_ratio_24 = row4["k5"] / max(row4["k10"], 1e-30)

    print("Interpretation:")
    print(
        "  cell1: R5b is suppressed because bubble CO is effectively zero and "
        f"Qp/Keq={row1['q_over_k']:.2e}, so the irreversible Gibbs clamp forces drv5={row1['drive_r5']:.3f}."
    )
    print(
        "  cell1: bubble tar is still non-zero, so R10b remains finite even when R5b is shut off."
    )
    print(
        "  cells2-4: drv5 recovers to ~1 once bubble CO appears; from there phi5 and phi10 are same-order, "
        f"but k5/k10 is already about {k_ratio_24:.1f}x, so R5b overtakes R10b."
    )
    print(
        f"  switch point: cell2 already has C_CO={row2['c_co']:.3f}, C_tar={row2['c_tar']:.3f}, "
        f"r5v/r10v={row2['r5v'] / max(row2['r10v'], 1e-30):.2f}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
