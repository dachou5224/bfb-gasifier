#!/usr/bin/env python3
"""Compare initialized LU bed hydrodynamics under multiple u_d closures."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0
from tests.validation_case_utils import build_phase1_htw_lu_global_nr_reactor_config, load_case_LU

_CLOSURES = (
    "current",
    "umf_over_epsmf",
    "hilligardt_eq311",
    "wein_1992_eq312",
    "backsolve_visible_epsb",
)


def _init_reactor(closure: str) -> Reactor:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.hydrodynamics_u_d_closure = closure
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
    print("=" * 132)
    print("LU initialized u_d-closure audit")
    print("=" * 132)
    print(f"{'closure':<24} {'eps_b top':>10} {'eps_dv top':>11} {'solid top':>11} {'u_d/u_mf top':>13}")
    print("-" * 132)
    for closure in _CLOSURES:
        reactor = _init_reactor(closure)
        cfg = reactor.config
        top_cells = [c for c in reactor.cells if float(c.geo.h_center / cfg.H_bed) >= 0.65]
        eps_b_top = float(sum(float(c.eps_b) for c in top_cells) / len(top_cells))
        eps_dv_top = float(sum(float(c.eps_d_voidage) for c in top_cells) / len(top_cells))
        solid_top = float(sum(1.0 - float(c.eps_d_voidage) for c in top_cells) / len(top_cells))
        ud_ratio_top = float(sum(float(c.u_d / max(c.u_mf, 1e-12)) for c in top_cells) / len(top_cells))
        print(
            f"{closure:<24} {eps_b_top:>10.4f} {eps_dv_top:>11.4f} {solid_top:>11.4f} {ud_ratio_top:>13.3f}"
        )
    print("-" * 132)
    print("readout:")
    print("  - Eq.3.11 / Eq.3.12 closures can be compared without solver-branch noise.")
    print("  - if a closure still yields top-bed eps_b far above ~0.6, it is likely outside bubbling-bed physics for HTW.")
    print("  - if a closure yields eps_d_voidage close to 0.99, it effectively wipes out dense-phase solids.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
