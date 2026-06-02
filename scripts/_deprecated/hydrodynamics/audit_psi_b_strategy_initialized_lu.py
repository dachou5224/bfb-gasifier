#!/usr/bin/env python3
"""Audit psi_b strategy sensitivity on initialized LU Heinbockel chain."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.physics.bubble_dynamics import bubble_interaction_factor
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0
from tests.validation_case_utils import build_phase1_htw_lu_global_nr_reactor_config, load_case_LU


def _init_reactor(*, psi_b_strategy: str) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.hydrodynamics_u_d_closure = "wein_1992_eq312"
    cfg.hydrodynamics_bubble_diameter_model = "hilligardt_ode"
    cfg.hydrodynamics_psi_b_strategy = psi_b_strategy
    cfg.hydrodynamics_bubble_velocity_strategy = "heinbockel_eq343"
    cfg.hydrodynamics_bubble_ode_strategy = "heinbockel_eq341"
    reactor = Reactor(cfg)
    t_profile = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        moisture_wt=cfg.moisture_wt,
        P=cfg.P,
    )
    generate_initial_x0(
        cells=reactor.cells,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=cfg.C_dry,
        H_dry=cfg.H_dry,
        O_dry=cfg.O_dry,
        moisture_wt=cfg.moisture_wt,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
        T_profile=t_profile,
        fuel_type=cfg.fuel_type,
    )
    reactor._set_bottom_cell_feeds()
    for i, cell in enumerate(reactor.cells):
        reactor._propagate_upstream(i)
        cell.calc_hydrodynamics()
    return reactor


def _summarize(name: str, reactor: Reactor) -> None:
    top = reactor.cells[-1]
    bot = reactor.cells[0]
    psi_bot = bubble_interaction_factor(bot.u_mf, strategy=bot.psi_b_strategy)
    psi_top = bubble_interaction_factor(top.u_mf, strategy=top.psi_b_strategy)
    print(f"[{name}]")
    print(f"  psi_b_cell0 = {psi_bot:.4f}")
    print(f"  psi_b_top   = {psi_top:.4f}")
    print(f"  u_b_top     = {top.u_b:.4f} m/s")
    print(f"  d_b_top     = {top.d_b:.4f} m")
    print(f"  eps_b_top   = {top.eps_b:.4f}")
    print(f"  K_bd_top    = {top.K_bd:.4f} 1/s")


def main() -> int:
    print("LU initialized psi_b sensitivity on Heinbockel chain")
    print("=" * 72)
    for strategy in ("technical_distributor", "wein_1992"):
        reactor = _init_reactor(psi_b_strategy=strategy)
        _summarize(strategy, reactor)
    print()
    print("readout:")
    print("  - this isolates psi_b source-of-truth while keeping u_d=Wein Eq.3.12 and Heinbockel Eq.3.41/3.43/3.44.")
    print("  - Eq.3.15 is now locked to psi_b = 0.17 * u_mf^(-0.33); this script only compares it against the technical-distributor constant 0.76.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
