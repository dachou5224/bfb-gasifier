#!/usr/bin/env python3
"""u_d 闭包 / 入口 split 敏感性：O2 气泡份额 vs 出口 CO/CO2（Phase2 LU 短 outer）。"""

from __future__ import annotations

from _lu_harness import (
    REF_CO,
    REF_CO2,
    add_json_output_arg,
    add_max_global_iter_arg,
    build_phase2_config,
    data_path,
    phase2_solve_kwargs,
    repo_root,
    run_phase2_init,
    write_json,
)

REPO_ROOT = repo_root()

import argparse
import copy
import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

from src.core.cell_kinetics import build_reaction_sources
from src.core.cell_thermo_helpers import gas_mole_fractions
from src.core.connectivity import _bottom_hydrodynamic_gas_inlet_dense_fraction
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
DRY = ("CO", "CO2", "H2", "CH4", "N2", "O2")


def _dry_y(N_b, N_d) -> dict[str, float]:
    y = gas_mole_fractions(N_b, N_d, "mixed")
    total = float(sum(y[idx[s]] for s in DRY))
    if total < 1e-15:
        return {s: 0.0 for s in DRY}
    return {s: float(y[idx[s]] / total) for s in DRY}


def _o2_zu_bubble_frac(cell) -> float:
    d = float(cell.N_zu_d[idx["O2"]])
    b = float(cell.N_zu_b[idx["O2"]])
    t = d + b
    return b / t if t > 1e-12 else float("nan")


def _r5_net_co2(cell) -> float:
    cell.calc_hydrodynamics()
    cell.calc_exchange()
    thermo = cell._get_local_thermo_bundle()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    if cell._vm_cache_valid:
        gas_src_vm = cell._vm_gas_source_cache.copy()
        solid_sink_vm = cell._vm_solid_sink_cache.copy()
    else:
        gas_src_vm = cell._calc_drying_pyrolysis_gas_source(tau_val)
        solid_sink_vm = cell.R_solid.copy()
    bundle = build_reaction_sources(
        T=cell.T,
        P=cell.P,
        fuel_type=cell.fuel_type,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=np.asarray(thermo["C_b"], dtype=np.float64),
        C_d=np.asarray(thermo["C_d"], dtype=np.float64),
        y_b=np.asarray(thermo["y_b"], dtype=np.float64),
        y_d=np.asarray(thermo["y_d"], dtype=np.float64),
        gas_src_vm=gas_src_vm,
        solid_sink_vm=solid_sink_vm,
        areas=cell._calc_char_surface_area_per_class(),
        solid_d_p=cell.solid.d_p,
        solid_d_p_classes=cell.solid.d_p_classes,
        D_g=float(thermo["D_g"]),
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        rho_ash_kg_m3=cell._ash_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=None,
        r4_scale=cell.r4_scale,
        r2_scale=cell.r2_scale,
        r5_scale=cell.r5_scale,
        r6_scale=cell.r6_scale,
        r7_scale=cell.r7_scale,
        r8_scale=cell.r8_scale,
        rate_multiplier=1.0,
        N_zu_d=cell.N_zu_d,
        N_d_in=cell.N_d_in,
        N_zu_b=cell.N_zu_b,
        N_b_in=cell.N_b_in,
        N_rez_d=cell.N_rez_d,
        N_rez_b=cell.N_rez_b,
        N_ex=cell.N_ex,
        solid_shape=cell.R_solid.shape,
        char_index=0,
    )
    j = idx["CO2"]
    return float(bundle.R_gas_b[j] + bundle.R_gas_d[j])


def _outlet_cell(reactor: Reactor):
    for cell in reversed(reactor.cells):
        if str(getattr(cell, "cell_type", "bed")) in {"freeboard", "bed"}:
            return cell
    return reactor.cells[-1]


def run_case(case_spec: dict, *, max_global_iter: int) -> dict:
    cfg = build_phase2_config(enable_vorab_march_staged=True)
    label = str(case_spec["label"])
    mode = str(case_spec.get("mode", "default"))

    if mode == "fixed_split":
        cfg.gas_inlet_split_strategy = "fixed"
        cfg.gas_inlet_dense_frac = float(case_spec["dense_frac"])
    elif mode == "u0_scale":
        area = np.pi / 4.0 * float(cfg.D_bed) ** 2
        primary = float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed)
        u0_ref = (primary * 8.314 * float(cfg.T_inlet) / float(cfg.P)) / max(area, 1e-12)
        cfg.u0_target = u0_ref * float(case_spec["u0_mult"])

    reactor = Reactor(cfg)
    u_d_override = case_spec.get("u_d_closure")
    if u_d_override:
        for cell in reactor.cells:
            cell.u_d_closure = str(u_d_override)

    init_strategy = str(phase2_solve_kwargs().get("nr_init_strategy", "vorabrechnung"))
    run_phase2_init(reactor, init_strategy=init_strategy, gs_warmup_steps=None)

    bed0 = reactor.cells[0]
    init_o2_b = _o2_zu_bubble_frac(bed0)
    init_dense = float(getattr(bed0, "_vorab_bottom_gas_inlet_dense_frac", np.nan))
    live_dense = _bottom_hydrodynamic_gas_inlet_dense_fraction(bed0)

    t0 = perf_counter()
    solve_kwargs = phase2_solve_kwargs()
    solve_kwargs["max_global_iter"] = int(max_global_iter)
    solve_kwargs["verbose"] = False
    solve_kwargs["skip_precalc"] = True
    result = reactor._solve_global_nr(
        max_iter=int(max_global_iter),
        tol=float(solve_kwargs.get("tol_global", 1.0)),
        verbose=False,
        init_strategy=init_strategy,
        jacobian_strategy=solve_kwargs.get("nr_jacobian_strategy"),
        skip_precalc=True,
    )
    wall_s = perf_counter() - t0

    outlet = _outlet_cell(reactor)
    y_out = _dry_y(outlet.N_b, outlet.N_d)
    y_b0 = _dry_y(bed0.N_b, bed0.N_d)
    solve_o2_b = _o2_zu_bubble_frac(bed0)

    return {
        "label": label,
        "mode": mode,
        "u_d_closure": u_d_override or "wein_locked",
        "max_global_iter": int(max_global_iter),
        "wall_s": float(wall_s),
        "converged": bool(result.get("converged")),
        "rms_scaled_final": float(result.get("rms_scaled_final", np.nan)),
        "init": {
            "u0": float(bed0.u0),
            "u_d": float(bed0.u_d),
            "eps_b": float(bed0.eps_b),
            "frozen_dense_frac": init_dense,
            "live_dense_frac": live_dense,
            "O2_zu_bubble_frac": init_o2_b,
            "O2_zu_d": float(bed0.N_zu_d[idx["O2"]]),
            "O2_zu_b": float(bed0.N_zu_b[idx["O2"]]),
        },
        "solve": {
            "O2_zu_bubble_frac": solve_o2_b,
            "bed0_y_CO": y_b0["CO"],
            "bed0_y_CO2": y_b0["CO2"],
            "bed0_y_O2": y_b0["O2"],
            "bed0_T_K": float(bed0.T),
            "bed0_R5_net_CO2_mol_s": _r5_net_co2(bed0),
            "outlet_y_CO": y_out["CO"],
            "outlet_y_CO2": y_out["CO2"],
            "outlet_T_K": float(outlet.T),
            "delta_CO_vs_ref": y_out["CO"] - REF_CO,
            "delta_CO2_vs_ref": y_out["CO2"] - REF_CO2,
        },
    }


CASES = [
    {"label": "wein_baseline", "mode": "default"},
    {"label": "hilligardt_ud", "mode": "default", "u_d_closure": "hilligardt_eq311"},
    {"label": "current_ud", "mode": "default", "u_d_closure": "current"},
    {"label": "fixed_dense30", "mode": "fixed_split", "dense_frac": 0.30},
    {"label": "fixed_dense50", "mode": "fixed_split", "dense_frac": 0.50},
    {"label": "fixed_dense70", "mode": "fixed_split", "dense_frac": 0.70},
    {"label": "u0_x0p9", "mode": "u0_scale", "u0_mult": 0.9},
    {"label": "u0_x1p1", "mode": "u0_scale", "u0_mult": 1.1},
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-global-iter", type=int, default=20)
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "data" / "phase2_ud_closure_o2_co_co2_sensitivity.json")
    parser.add_argument("--cases", nargs="*", default=None, help="subset of case labels")
    args = parser.parse_args()

    selected = CASES
    if args.cases:
        wanted = set(args.cases)
        selected = [c for c in CASES if c["label"] in wanted]

    results = []
    for spec in selected:
        print(f"Running {spec['label']} ...", flush=True)
        results.append(run_case(spec, max_global_iter=args.max_global_iter))

    payload = {
        "ref_CO": REF_CO,
        "ref_CO2": REF_CO2,
        "cases": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
