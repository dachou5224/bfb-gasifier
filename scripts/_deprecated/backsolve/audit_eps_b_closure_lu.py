#!/usr/bin/env python3
"""Audit epsilon_b closures on the initialized LU bed state."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.physics.phase_fractions import calc_epsilon_b, calc_visible_bubble_fraction
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0
from tests.validation_case_utils import build_phase1_htw_lu_global_nr_reactor_config, load_case_LU


def _init_reactor(*, bubble_model: str) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    cfg.hydrodynamics_bubble_diameter_model = bubble_model
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


def _print_branch(name: str, reactor: Reactor) -> None:
    cfg = reactor.config
    print(f"[{name}]")
    print(f"{'cell':>4} {'xi/H':>6} {'u0':>8} {'u_mf':>8} {'u_b':>8} {'u_d':>8} {'eps_ex':>8} {'eps_vis':>8} {'delta':>8}")
    for i, cell in enumerate(reactor.cells):
        eps_ex = float(calc_epsilon_b(cell.u0, cell.u_mf, cell.u_b))
        eps_vis = float(calc_visible_bubble_fraction(cell.u0, cell.u_b, cell.u_d))
        print(
            f"{i:>4d} {cell.geo.h_center/cfg.H_bed:>6.2f} {cell.u0:>8.4f} {cell.u_mf:>8.4f} "
            f"{cell.u_b:>8.4f} {cell.u_d:>8.4f} {eps_ex:>8.4f} {eps_vis:>8.4f} {eps_vis-eps_ex:>8.4e}"
        )
    print()


def main() -> int:
    print("LU epsilon_b closure audit on initialized bed state")
    print("=" * 96)
    _print_branch("mori_wen", _init_reactor(bubble_model="mori_wen"))
    _print_branch("hilligardt_ode", _init_reactor(bubble_model="hilligardt_ode"))
    print("readout:")
    print("  - if eps_excess and eps_visible collapse to the same value, the chosen u_d closure is enforcing that identity.")
    print("  - in that case the next source-of-truth question is not which epsilon formula is used, but whether u_d should be solved from it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
