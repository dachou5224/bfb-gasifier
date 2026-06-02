#!/usr/bin/env python3
"""Audit Heinbockel pressurized bubble-chain on initialized LU state."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0
from tests.validation_case_utils import build_phase1_htw_lu_global_nr_reactor_config, load_case_LU


def _init_reactor(*, velocity_strategy: str, ode_strategy: str) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.hydrodynamics_u_d_closure = "wein_1992_eq312"
    cfg.hydrodynamics_bubble_diameter_model = "hilligardt_ode"
    cfg.hydrodynamics_bubble_velocity_strategy = velocity_strategy
    cfg.hydrodynamics_bubble_ode_strategy = ode_strategy
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
    print(f"{'cell':>4} {'xi/H':>6} {'d_b':>8} {'u_b':>8} {'u_d':>8} {'eps_b':>8} {'eps_dv':>8} {'K_bd':>8}")
    for i, cell in enumerate(reactor.cells):
        print(
            f"{i:>4d} {cell.geo.h_center/cfg.H_bed:>6.2f} {cell.d_b:>8.4f} {cell.u_b:>8.4f} {cell.u_d:>8.4f} "
            f"{cell.eps_b:>8.4f} {cell.eps_d_voidage:>8.4f} {cell.K_bd:>8.4f}"
        )
    print()


def main() -> int:
    print("LU initialized Heinbockel audit (u_d = Wein Eq.3.12)")
    print("=" * 96)
    hill = _init_reactor(
        velocity_strategy="hilligardt_eq313",
        ode_strategy="hilligardt_eq333",
    )
    hamel_ocr = _init_reactor(
        velocity_strategy="heinbockel_eq343",
        ode_strategy="heinbockel_eq341",
    )
    hein = _init_reactor(
        velocity_strategy="heinbockel_eq343_legacy_2p14_p07",
        ode_strategy="heinbockel_eq341",
    )
    _print_branch("hilligardt_atm_chain", hill)
    _print_branch("heinbockel_pressurized_chain_hamel_ocr", hamel_ocr)
    _print_branch("heinbockel_pressurized_chain_legacy_2p14_p07", hein)
    print("readout:")
    print("  - this isolates the pressure-chain effect after fixing u_d to Wein Eq.3.12.")
    print("  - compare Hamel OCR Eq.3.43 against the previously over-aggressive 2.14*(P/P0)^0.7 legacy branch.")
    print("  - Eq.3.41 now uses n_B=5 (coalescing bubbles), while u_br still uses n_b=2.7.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
