"""Audit bubble-phase timescales and volumetric-throughput consistency.

Purpose
-------
Quantify whether bubble homogeneous oxidation is operating far faster than
cell transit / interphase exchange, and whether the tracked bubble molar flow
is consistent with the hydrodynamic bubble holdup.

Usage
-----
    python3 scripts/audit_bubble_timescales.py
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
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as IDX
from src.kinetics.gas_reactions import rate_R5_bubble
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


def main() -> int:
    sc = Scenario(name="stable_window_candidate", dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)
    reactor, result = run_scenario(sc)
    idx = IDX

    print("=" * 136)
    print(
        f"Scenario: {sc.name} | dense={sc.dense:.2f} heat_loss={sc.heat_loss:.2f} "
        f"R5={sc.r5s:.2f} R6={sc.r6s:.2f} R4={sc.r4s:.2f} | "
        f"Texit={result['T_profile'][-1]:.1f}K Tpeak={np.max(result['T_profile']):.1f}K"
    )
    print("-" * 136)
    print(
        f"{'cell':>4} {'T[K]':>8} {'tau_b[s]':>9} {'tau_ex[s]':>10} {'tau_r5[s]':>11} "
        f"{'F_hydro':>9} {'F_state':>9} {'ratio':>7} {'NexO2':>9} {'O2_b':>9} {'CO_b':>9}"
    )

    for i, cell in enumerate(reactor.cells):
        cell.calc_hydrodynamics()
        cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
        cell.calc_exchange()

        c_b = cell._concentrations("b")
        y_b = cell._mole_fractions("b")
        r5b_v = rate_R5_bubble(cell.T, c_b[idx["CO"]], c_b[idx["O2"]], cell.P, y_b, idx)

        tau_b = cell.geo.dh / max(cell.u_b, 1e-12)
        tau_ex = 1.0 / max(cell.K_bd, 1e-12)
        tau_r5 = c_b[idx["O2"]] / max(r5b_v, 1e-30)

        f_hydro = cell.V_b / max(tau_b, 1e-12)
        n_b_total = float(np.sum(np.maximum(cell.N_b, 0.0)))
        f_state = n_b_total * 8.314462618 * cell.T / cell.P
        ratio = f_state / max(f_hydro, 1e-12)

        print(
            f"{i:4d} {cell.T:8.1f} {tau_b:9.4f} {tau_ex:10.4f} {tau_r5:11.6f} "
            f"{f_hydro:9.4f} {f_state:9.4f} {ratio:7.3f} "
            f"{cell.N_ex[idx['O2']]:9.4f} {cell.N_b[idx['O2']]:9.4f} {cell.N_b[idx['CO']]:9.4f}"
        )

    print("-" * 136)
    print("Interpretation:")
    print("  tau_r5 << tau_b, tau_ex  => bubble R5 is effectively instantaneous relative to transport.")
    print("  ratio < 1                => tracked bubble molar flow is smaller than hydrodynamic bubble throughput.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
