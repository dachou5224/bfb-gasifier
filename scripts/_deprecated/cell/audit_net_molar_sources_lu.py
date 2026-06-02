#!/usr/bin/env python3
"""Audit LU bed net molar gas source by reaction family on the shared NR path."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import S_CHAR
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from src.core.species import gas_diffusivity_correlation
from src.core.cell_kinetics import build_reaction_sources
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_case_LU,
)


def _bundle_for_cell(cell):
    cell.calc_hydrodynamics()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau_val)
    cell.calc_exchange()

    gas_src_vm = cell._vm_gas_source_cache.copy()
    solid_sink_vm = cell._vm_solid_sink_cache.copy()
    return build_reaction_sources(
        T=cell.T,
        P=cell.P,
        fuel_type=cell.fuel_type,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=cell._concentrations("b"),
        C_d=cell._concentrations("d"),
        y_b=cell._mole_fractions("b"),
        y_d=cell._mole_fractions("d"),
        gas_src_vm=gas_src_vm,
        solid_sink_vm=solid_sink_vm,
        areas=cell._calc_char_surface_area_per_class(),
        solid_d_p=cell.solid.d_p,
        D_g=gas_diffusivity_correlation(cell.T, cell.P),
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=cell.calc_minor_species_gibbs() if cell.use_gibbs_minor else None,
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


def main() -> int:
    case = load_case_LU()
    cfg = build_phase2_htw_lu_freeboard_reactor_config(case)
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    print("=" * 178)
    print("LU bed net molar source audit")
    print("=" * 178)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"converged={result.get('converged')} n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"bed_exit_T={float(result['bed_T_profile'][-1]):.1f}K"
    )
    print("-" * 178)
    print(
        f"{'cell':>4} {'xi/H':>6} {'u0':>8} {'Ntot':>10} {'net_total':>10} "
        f"{'vm':>8} {'R5':>8} {'R6':>8} {'R7':>8} {'R10':>8} {'R11':>8} {'R12':>8} {'char':>8}"
    )
    print("-" * 178)

    for i, cell in enumerate(reactor.cells):
        bundle = _bundle_for_cell(cell)
        n_tot = float(np.sum(cell.N_b + cell.N_d))
        xi = float(cell.geo.h_center / cfg.H_bed)
        print(
            f"{i:>4d} {xi:>6.2f} {cell.u0:>8.4f} {n_tot:>10.4f} {bundle.net_molar_gas_source_total:>10.4f} "
            f"{bundle.net_molar_gas_source_vm:>8.3f} {bundle.net_molar_gas_source_r5:>8.3f} "
            f"{bundle.net_molar_gas_source_r6:>8.3f} {bundle.net_molar_gas_source_r7:>8.3f} "
            f"{bundle.net_molar_gas_source_r10:>8.3f} {bundle.net_molar_gas_source_r11:>8.3f} "
            f"{bundle.net_molar_gas_source_r12:>8.3f} {bundle.net_molar_gas_source_char:>8.3f}"
        )

    print("-" * 178)
    print("readout:")
    print("  - `net_total = sum(R_gas_b + R_gas_d)` is the local chemistry-driven gas molar source seen by `u0`.")
    print("  - positive terms lift local superficial velocity; negative terms pull `u0` down.")
    print("  - on the Hamel-thesis path, `R7` should stay non-negative; upper-bed sinks should now be read mainly through char/WGSR/tar terms.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
