#!/usr/bin/env python3
"""Compare bed hydrodynamics under Mori-Wen vs Hilligardt ODE on LU case.

This audit is intentionally hydrodynamics-first:
- keep chemistry/solver settings identical
- switch only the bed bubble-diameter model
- read out how d_b / eps_b / eps_d_voidage / K_bd shift along the bed
- operate on the Vorabrechnung-initialized LU state to avoid solver-side branch noise
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0
from tests.validation_case_utils import build_phase1_htw_lu_global_nr_reactor_config, load_case_LU


def _run(branch: str) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.hydrodynamics_u_d_closure = "backsolve_visible_epsb"
    cfg.hydrodynamics_bubble_diameter_model = branch
    cfg.hydrodynamics_psi_b_strategy = "technical_distributor"
    cfg.hydrodynamics_lambda_strategy = "hamel_280"
    cfg.hydrodynamics_xi_strategy = "hamel_regime"
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


def main() -> int:
    reactor_mw = _run("mori_wen")
    reactor_ode = _run("hilligardt_ode")

    print("=" * 178)
    print("LU initialized bed hydrodynamics audit: Mori-Wen vs Hilligardt ODE")
    print("=" * 178)
    print(
        f"{'cell':>4} {'xi/H':>6} "
        f"{'d_b MW':>10} {'d_b ODE':>10} "
        f"{'eps_b MW':>10} {'eps_b ODE':>10} "
        f"{'eps_dv MW':>10} {'eps_dv ODE':>10} "
        f"{'Kbd MW':>10} {'Kbd ODE':>10}"
    )
    print("-" * 178)
    for i, (cell_mw, cell_ode) in enumerate(zip(reactor_mw.cells, reactor_ode.cells)):
        xi = float(cell_mw.geo.h_center / reactor_mw.config.H_bed)
        print(
            f"{i:>4d} {xi:>6.2f} "
            f"{cell_mw.d_b:>10.4f} {cell_ode.d_b:>10.4f} "
            f"{cell_mw.eps_b:>10.4f} {cell_ode.eps_b:>10.4f} "
            f"{cell_mw.eps_d_voidage:>10.4f} {cell_ode.eps_d_voidage:>10.4f} "
            f"{cell_mw.K_bd:>10.4f} {cell_ode.K_bd:>10.4f}"
        )

    print("-" * 178)
    print("readout:")
    print("  - state = Vorabrechnung-initialized only; no global solve branch-selection noise is included.")
    print("  - default path still remains Mori-Wen; this script only exercises the new audit branch.")
    print("  - if ODE mainly changes d_b/K_bd while leaving eps_b collapse unresolved, u_d closure still dominates.")
    print("  - if ODE materially shifts top-bed eps_b/eps_d_voidage/K_bd together, the MW-vs-ODE structural gap is now active.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
