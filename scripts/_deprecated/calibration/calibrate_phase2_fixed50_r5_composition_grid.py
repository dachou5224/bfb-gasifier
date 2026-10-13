#!/usr/bin/env python3
"""Phase2 LU zone=0：fixed50 × r5_scale 组成 + 轴向 T 剖面标定。

目标：最小化 reactor exit dry CO/CO₂ 相对 Hamel CASE_HTW_WESSELING_1 偏差，
并记录 bed 顶温与 Hamel ξ 剖面的对照。

Usage
-----
    python3 scripts/calibrate_phase2_fixed50_r5_composition_grid.py --max-global-iter 40
    python3 scripts/calibrate_phase2_fixed50_r5_composition_grid.py --max-global-iter 150 --labels fixed50_r5_050
"""

from __future__ import annotations

from _lu_harness import (
    REF_CO,
    REF_CO2,
    LADDER_P2,
    add_max_global_iter_arg,
    apply_config_overrides,
    build_phase2_config,
    data_path,
    load_case_lu,
    phase2_solve_kwargs,
    repo_root,
    write_json,
)

from tests.validation_case_utils import load_validation_case_node

REPO_ROOT = repo_root()

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import numpy as np

from src.core.cell_thermo_helpers import gas_mole_fractions
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as GSI


@dataclass(frozen=True)
class GridCase:
    label: str
    fixed50: bool
    r5_scale: float


DEFAULT_GRID: tuple[GridCase, ...] = (
    GridCase("default70_r5_075", fixed50=False, r5_scale=0.75),
    GridCase("fixed50_r5_075", fixed50=True, r5_scale=0.75),
    GridCase("fixed50_r5_050", fixed50=True, r5_scale=0.50),
    GridCase("fixed50_r5_035", fixed50=True, r5_scale=0.35),
    GridCase("fixed50_r5_025", fixed50=True, r5_scale=0.25),
    GridCase("fixed50_r5_015", fixed50=True, r5_scale=0.15),
)


def _apply_fixed50(cfg) -> None:
    cfg.gas_inlet_split_strategy = "fixed"
    cfg.gas_inlet_dense_frac = 0.50
    cfg.vorab_staged_primary_gas_dense_frac_thesis = 0.50


def _score_exit(co: float, co2: float) -> float:
    return abs(co - REF_CO) / REF_CO + abs(co2 - REF_CO2) / REF_CO2


def _load_hamel_axial() -> dict[str, np.ndarray]:
    node = load_validation_case_node("CASE_HTW_WESSELING_1")
    prof = node["outputs"]["axial_profiles"]
    xi = np.asarray(prof["xi"], dtype=np.float64)
    t_k = np.asarray(prof["T_K"], dtype=np.float64)
    co_wet = np.asarray(prof["CO_mol_wet"], dtype=np.float64)
    co2_wet = np.asarray(prof["CO2_mol_wet"], dtype=np.float64)
    mask = np.isfinite(t_k)
    return {
        "xi": xi[mask],
        "T_K": t_k[mask],
        "CO_wet": co_wet[np.isfinite(co_wet)],
        "CO2_wet": co2_wet[np.isfinite(co2_wet)],
        "xi_co": xi[np.isfinite(co_wet)],
        "xi_co2": xi[np.isfinite(co2_wet)],
    }


def _hamel_interp(xi_query: float, xi_ref: np.ndarray, y_ref: np.ndarray) -> float:
    if xi_ref.size == 0:
        return float("nan")
    order = np.argsort(xi_ref)
    return float(np.interp(float(xi_query), xi_ref[order], y_ref[order]))


def _axial_profile(reactor: Reactor, case: dict) -> list[dict[str, float]]:
    h_total = float(case.get("reactor_height_m", case["bed_height_m"] + case.get("H_freeboard", 0.0)))
    if h_total <= 1.0e-12:
        h_total = float(reactor.config.H_bed + reactor.config.H_freeboard)
    rows: list[dict[str, float]] = []
    for i, cell in enumerate(reactor.cells):
        ctype = str(getattr(cell, "cell_type", "bed"))
        xi = float(cell.geo.h_center / h_total)
        y = gas_mole_fractions(cell.N_b, cell.N_d, "mixed")
        rows.append(
            {
                "cell": float(i),
                "cell_type": ctype,
                "xi_reactor": xi,
                "T_K": float(cell.T),
                "y_CO_wet": float(y[GSI["CO"]]),
                "y_CO2_wet": float(y[GSI["CO2"]]),
                "y_O2_wet": float(y[GSI["O2"]]),
            }
        )
    return rows


def _profile_temperature_rmse(profile: list[dict[str, float]], hamel: dict[str, np.ndarray]) -> float:
    bed_pts = [p for p in profile if p["cell_type"] == "bed" and p["xi_reactor"] <= 0.55]
    if not bed_pts:
        bed_pts = [p for p in profile if p["cell_type"] == "bed"]
    errs = []
    for p in bed_pts:
        t_ref = _hamel_interp(p["xi_reactor"], hamel["xi"], hamel["T_K"])
        if np.isfinite(t_ref):
            errs.append((float(p["T_K"]) - t_ref) ** 2)
    return float(np.sqrt(np.mean(errs))) if errs else float("nan")


def _run_grid_case(grid: GridCase, *, max_global_iter: int) -> dict[str, object]:
    case = load_case_lu()
    cfg = build_phase2_config(case)
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 0
    cfg.vorab_major_gibbs_x0 = False
    cfg.r5_scale = float(grid.r5_scale)
    if grid.fixed50:
        _apply_fixed50(cfg)
    for key, value in {
        "nr_line_search_gas_phase_split_merit_thesis": True,
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1180.0,
        "freeboard_trajectory_coeff_model": "exact_hamel",
    }.items():
        setattr(cfg, key, value)

    solve_kwargs = phase2_solve_kwargs()
    solve_kwargs["max_global_iter"] = int(max_global_iter)
    solve_kwargs["verbose"] = False

    t0 = perf_counter()
    reactor = Reactor(cfg)
    out = reactor.solve(
        max_global_iter=int(max_global_iter),
        tol_global=float(solve_kwargs.get("tol_global", 1.0)),
        verbose=False,
        nr_init_strategy=str(solve_kwargs.get("init_strategy", "vorabrechnung")),
        nr_jacobian_strategy=solve_kwargs.get("nr_jacobian_strategy"),
    )
    wall_s = perf_counter() - t0

    profile = _axial_profile(reactor, case)
    hamel = _load_hamel_axial()

    bed_cells = [p for p in profile if p["cell_type"] == "bed"]
    bed_top = bed_cells[-1] if bed_cells else profile[-1]
    bed0 = bed_cells[0] if bed_cells else profile[0]

    exit_dry = out.get("exit_gas_dry") or {}
    bed_dry = out.get("bed_exit_gas_dry") or {}
    co = float(exit_dry.get("CO", 0.0))
    co2 = float(exit_dry.get("CO2", 0.0))
    score_exit = _score_exit(co, co2)
    t_rmse_bed = _profile_temperature_rmse(profile, hamel)
    score_total = score_exit + 0.002 * t_rmse_bed

    return {
        "label": grid.label,
        "fixed50": bool(grid.fixed50),
        "r5_scale": float(grid.r5_scale),
        "gas_inlet_dense_frac": float(cfg.gas_inlet_dense_frac),
        "max_global_iter": int(max_global_iter),
        "wall_time_s": float(wall_s),
        "converged_fully": bool(out.get("converged_fully", False)),
        "converged": bool(out.get("converged", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", 0.0) or 0.0),
        "nr_outer_iters": int(out.get("nr_outer_iters", 0)),
        "reactor_exit_gas_dry": {"CO": co, "CO2": co2},
        "bed_exit_gas_dry": {
            "CO": float(bed_dry.get("CO", 0.0)),
            "CO2": float(bed_dry.get("CO2", 0.0)),
        },
        "score_exit_co_co2": float(score_exit),
        "score_total": float(score_total),
        "bed0_T_K": float(bed0["T_K"]),
        "bed_top_T_K": float(bed_top["T_K"]),
        "bed_top_xi": float(bed_top["xi_reactor"]),
        "T_rmse_bed_vs_hamel_K": float(t_rmse_bed),
        "ref_dry": {"CO": REF_CO, "CO2": REF_CO2},
        "axial_profile": profile,
        "hamel_xi": hamel["xi"].tolist(),
        "hamel_T_K": hamel["T_K"].tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=40)
    parser.add_argument("--labels", nargs="*", default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_fixed50_r5_composition_grid.json",
    )
    args = parser.parse_args()

    label_set = set(args.labels) if args.labels else None
    cases = [c for c in DEFAULT_GRID if label_set is None or c.label in label_set]
    if not cases:
        raise SystemExit("no matching grid cases")

    results = [_run_grid_case(c, max_global_iter=int(args.max_global_iter)) for c in cases]
    results.sort(key=lambda r: float(r["score_total"]))

    payload = {
        "max_global_iter": int(args.max_global_iter),
        "best_by_score_total": results[0]["label"] if results else None,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Wrote {args.output}")
    print(f"{'label':<22} {'conv':>4} {'rms':>7} {'CO':>6} {'CO2':>6} {'bedT':>6} {'T_rmse':>7} {'score':>6}")
    for row in results:
        conv = "Y" if row["converged_fully"] else "n"
        ex = row["reactor_exit_gas_dry"]
        print(
            f"{row['label']:<22} {conv:>4} {row['rms_scaled_final']:>7.4f} "
            f"{ex['CO']:>6.3f} {ex['CO2']:>6.3f} {row['bed_top_T_K']:>6.0f} "
            f"{row['T_rmse_bed_vs_hamel_K']:>7.1f} {row['score_total']:>6.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
