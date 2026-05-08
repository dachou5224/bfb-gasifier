#!/usr/bin/env python3
"""Audit CASE_HTW_WESSELING_1 JSON -> current ReactorConfig mapping.

This script intentionally uses ``tests.validation_case_utils`` builders instead
of a bare ``ReactorConfig`` so it reports the configuration actually exercised
by the LU validation scripts.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(".").resolve()))

from tests.validation_case_utils import (  # noqa: E402
    build_phase1_htw_lu_global_nr_reactor_config,
    build_phase1_htw_lu_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_case_LU,
)


SEP = "=" * 82


def _fmt(value: Any, digits: int = 5) -> str:
    if value is None:
        return "None"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def _pct(model: float, ref: float) -> float:
    return 100.0 * (float(model) - float(ref)) / (abs(float(ref)) + 1e-30)


def _check(label: str, model: float, ref: float | None, unit: str = "", tol_pct: float = 2.0) -> bool:
    if ref is None:
        print(f"  ?  {label:<38s} model={_fmt(model):>12s}{unit}  ref={'missing':>12s}")
        return True
    pct = _pct(model, ref)
    ok = abs(pct) <= tol_pct
    tag = "OK" if ok else "!!"
    print(
        f"  {tag} {label:<38s} model={_fmt(model):>12s}{unit}  "
        f"ref={_fmt(ref):>12s}{unit}  delta={pct:+.2f}%"
    )
    return ok


def _nm3h_air_to_mol_s(air_nm3_h: float) -> tuple[float, float]:
    mol_total = float(air_nm3_h) * 1000.0 / 22.414 / 3600.0
    return 0.21 * mol_total, 0.79 * mol_total


def main() -> int:
    case = load_case_LU()
    phase1 = build_phase1_htw_lu_reactor_config(case)
    thesis = build_phase1_htw_lu_global_nr_reactor_config(case)
    phase2 = build_phase2_htw_lu_freeboard_reactor_config(case)

    particle_range = tuple(case.get("particle_diameter_mm_range", (1.5, 3.0)))
    if len(particle_range) != 2:
        particle_range = (1.5, 3.0)
    d_p_ref_mm = 0.5 * (float(particle_range[0]) + float(particle_range[1]))
    n_size_ref = int(case.get("n_size_classes", 10))

    total_o2_ref, total_n2_ref = _nm3h_air_to_mol_s(float(case["total_air_Nm3_h"]))
    sec_o2_ref, sec_n2_ref = _nm3h_air_to_mol_s(float(case["secondary_agent_Nm3_h"]))

    print(SEP)
    print("CASE_HTW_WESSELING_1 input mapping audit")
    print(SEP)
    print("config builders:")
    print("  phase1 : build_phase1_htw_lu_reactor_config")
    print("  thesis : build_phase1_htw_lu_global_nr_reactor_config")
    print("  phase2 : build_phase2_htw_lu_freeboard_reactor_config")

    print("\n[1] Fuel / geometry / operating conditions")
    ok = True
    ok &= _check("fuel_feed", phase1.fuel_feed, float(case["fuel_feed"]) / 3600.0, " kg/s", 0.5)
    ok &= _check("moisture_wt", phase1.moisture_wt, float(case["moisture_wt"]), " wt%", 0.1)
    ok &= _check("ash_dry_wt", phase1.ash_dry_wt, float(case["ash_dry_wt"]), " wt%", 0.1)
    ok &= _check("VM_daf", phase1.VM_daf, float(case["VM_daf"]), " wt%", 0.1)
    ok &= _check("C_dry", phase1.C_dry, float(case["C_dry"]), " wt%", 0.1)
    ok &= _check("H_dry", phase1.H_dry, float(case["H_dry"]), " wt%", 0.1)
    ok &= _check("O_dry", phase1.O_dry, float(case["O_dry"]), " wt%", 0.1)
    ok &= _check("N_dry", phase1.nitrogen_fraction, float(case["N_dry"]), " wt%", 0.5)
    ok &= _check("S_dry", phase1.S_dry, float(case["S_dry"]), " wt%", 0.5)
    ok &= _check("H_bed", phase1.H_bed, float(case["H_bed"]), " m", 0.5)
    ok &= _check("H_total phase2", phase2.H_bed + phase2.H_freeboard, float(case["reactor_height_m"]), " m", 0.5)
    ok &= _check("D_bed", phase1.D_bed, float(case["bed_diameter_m"]), " m", 0.5)
    ok &= _check("P", phase1.P, float(case["P"]), " Pa", 0.1)
    ok &= _check("ER", phase1.ER, float(case["ER"]), "", 0.5)

    print("\n[2] Gasification-agent budget")
    ok &= _check("O2 total from ER", phase1.O2_feed, total_o2_ref, " mol/s", 1.0)
    ok &= _check("N2 total from air", phase1.N2_feed, total_n2_ref, " mol/s", 1.0)
    print(f"  ?  steam/H2O feed                         model={phase1.H2O_feed:12.5g} mol/s  ref={'missing':>12s}")
    print(f"     steam_to_o2_molar={phase1.steam_to_o2_molar:.3f}; CASE JSON has no explicit steam feed for Sim 1")
    ok &= _check("phase2 primary O2 after secondary split", phase2.O2_feed, total_o2_ref - sec_o2_ref, " mol/s", 1.0)
    ok &= _check("phase2 primary N2 after secondary split", phase2.N2_feed, total_n2_ref - sec_n2_ref, " mol/s", 1.0)
    ok &= _check("phase2 secondary O2", phase2.freeboard_secondary_O2_mol_s, sec_o2_ref, " mol/s", 1.0)
    ok &= _check("phase2 secondary N2", phase2.freeboard_secondary_N2_mol_s, sec_n2_ref, " mol/s", 1.0)
    ok &= _check("secondary injection xi", phase2.freeboard_secondary_injection_xi or 0.0, float(case["secondary_injection_xi"]), "", 0.5)

    print("\n[3] Particle / solid discretization")
    d_p_phase1_mm = 1000.0 * float(phase1.d_p)
    print(f"  ref particle range                         {particle_range[0]:.3g}-{particle_range[1]:.3g} mm; midpoint={d_p_ref_mm:.3g} mm")
    d_p_ok = float(particle_range[0]) <= d_p_phase1_mm <= float(particle_range[1])
    tag = "OK" if d_p_ok else "!!"
    print(f"  {tag} phase1 d_p                            model={d_p_phase1_mm:12.5g} mm  ref range={particle_range}")
    n_ok = int(phase1.n_age_classes) == n_size_ref
    tag = "OK" if n_ok else "!!"
    print(f"  {tag} phase1 n_age_classes                  model={phase1.n_age_classes:12d}      ref={n_size_ref:12d}")
    print(f"  info thesis d_p/n_age_classes               d_p={1000.0*float(thesis.d_p):.5g} mm, n={thesis.n_age_classes}")
    print("  note current LU builders now map the Sim 1 particle range to uniform size classes.")
    ok &= d_p_ok and n_ok

    print("\n[4] Local tuning / thesis-mode switches")
    print(f"  phase1 gas_inlet_dense_frac={phase1.gas_inlet_dense_frac:.3f}, r4={phase1.r4_scale:.3f}, r5={phase1.r5_scale:.3f}, r7={phase1.r7_scale:.3f}")
    print(f"  thesis gas_inlet_dense_frac={thesis.gas_inlet_dense_frac:.3f}, thesis_mode={thesis.thesis_mode}, use_gibbs_minor={thesis.use_gibbs_minor}")
    print(f"  recycle_gas={phase1.recycle_gas}, recirculation_frac={phase1.recirculation_frac:.3f}, top_solid_inlet_frac={phase1.top_solid_inlet_frac:.3f}")

    print("\n[5] Summary")
    if ok:
        print("  PASS: no high-risk mismatch detected by this input audit.")
        return 0
    print("  ATTENTION: current LU builder still differs from CASE JSON / Hamel extraction in particle discretization or other fields above.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
