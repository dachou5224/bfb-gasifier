#!/usr/bin/env python3
"""Phase2 LU（fixed50 + r5=0.25 + vm_zone=0）：温度链路线扫描。

覆盖 Vorab 冷预算 march、NR 初温 seed 模式（adiabatic / budget / blend）、
床底 clip 窗、NR 温度围栏与全局 heat_loss 等杠杆，并对标 Hamel Fig 7.4 轴向 T。

Usage
-----
    python3 scripts/benchmark_phase2_temperature_chain_lu.py --max-global-iter 40
    python3 scripts/benchmark_phase2_temperature_chain_lu.py --labels seed_budget march_hl0_blend
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
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter

import numpy as np

from src.core.cell_thermo_helpers import gas_mole_fractions
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as GSI


@dataclass(frozen=True)
class TempChainCase:
    label: str
    overrides: dict[str, object] = field(default_factory=dict)


def _score_exit(co: float, co2: float) -> float:
    return abs(co - REF_CO) / REF_CO + abs(co2 - REF_CO2) / REF_CO2


def _load_hamel_axial() -> dict[str, np.ndarray]:
    node = load_validation_case_node("CASE_HTW_WESSELING_1")
    prof = node["outputs"]["axial_profiles"]
    xi = np.asarray(prof["xi"], dtype=np.float64)
    t_k = np.asarray(prof["T_K"], dtype=np.float64)
    mask = np.isfinite(t_k)
    return {"xi": xi[mask], "T_K": t_k[mask]}


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


def _profile_temperature_rmse(
    profile: list[dict[str, float]],
    hamel: dict[str, np.ndarray],
    *,
    scope: str,
) -> float:
    if scope == "bed":
        pts = [p for p in profile if p["cell_type"] == "bed" and p["xi_reactor"] <= 0.55]
        if not pts:
            pts = [p for p in profile if p["cell_type"] == "bed"]
    else:
        pts = [p for p in profile if p["xi_reactor"] <= 1.0]
    errs = []
    for p in pts:
        t_ref = _hamel_interp(p["xi_reactor"], hamel["xi"], hamel["T_K"])
        if np.isfinite(t_ref):
            errs.append((float(p["T_K"]) - t_ref) ** 2)
    return float(np.sqrt(np.mean(errs))) if errs else float("nan")


def _apply_baseline(cfg) -> None:
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 0
    cfg.vorab_major_gibbs_x0 = False
    cfg.r5_scale = 0.25
    cfg.gas_inlet_split_strategy = "fixed"
    cfg.gas_inlet_dense_frac = 0.50
    cfg.vorab_staged_primary_gas_dense_frac_thesis = 0.50
    for key, value in {
        "nr_line_search_gas_phase_split_merit_thesis": True,
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1180.0,
        "freeboard_trajectory_coeff_model": "exact_hamel",
    }.items():
        setattr(cfg, key, value)


DEFAULT_CASES: tuple[TempChainCase, ...] = (
    TempChainCase("baseline"),
    TempChainCase(
        "seed_budget",
        {
            "vorab_nr_temperature_seed_mode_thesis": "budget_progress",
        },
    ),
    TempChainCase(
        "seed_blend_50",
        {
            "vorab_nr_temperature_seed_mode_thesis": "blend",
            "vorab_nr_temperature_seed_adiabatic_blend_thesis": 0.50,
        },
    ),
    TempChainCase(
        "seed_hamel_1020_1100",
        {
            "vorab_nr_temperature_seed_mode_thesis": "budget_progress",
            "vorab_nr_temperature_seed_bottom_K_thesis": 1020.0,
            "vorab_nr_temperature_seed_top_K_thesis": 1100.0,
        },
    ),
    TempChainCase(
        "adiab_cap_1240",
        {
            "vorab_bottom_temperature_cap_K_thesis": 1240.0,
        },
    ),
    TempChainCase(
        "fence_lo_80",
        {
            "nr_temperature_fence_lower_margin_K_thesis": 80.0,
        },
    ),
    TempChainCase(
        "heat_loss_half",
        {
            "heat_loss_frac": 0.05,
        },
    ),
    TempChainCase(
        "march_hl0_blend",
        {
            "vorab_bed_temperature_march_thesis": True,
            "vorab_bed_temperature_march_heat_loss_frac_thesis": 0.0,
            "vorab_nr_temperature_seed_mode_thesis": "blend",
            "vorab_nr_temperature_seed_adiabatic_blend_thesis": 0.50,
        },
    ),
    TempChainCase(
        "march_budget_seed",
        {
            "vorab_bed_temperature_march_thesis": True,
            "vorab_bed_temperature_march_heat_loss_frac_thesis": 0.05,
            "vorab_nr_temperature_seed_mode_thesis": "budget_progress",
        },
    ),
    TempChainCase(
        "seed_peak1240",
        {
            "vorab_nr_temperature_seed_mode_thesis": "budget_progress",
            "vorab_nr_temperature_seed_bottom_K_thesis": 750.0,
            "vorab_nr_temperature_seed_top_K_thesis": 1240.0,
            "vorab_bottom_temperature_cap_K_thesis": 1325.0,
            "vorab_nr_temperature_fence_upper_anchor_K_thesis": 1500.0,
        },
    ),
)


def _run_case(spec: TempChainCase, *, max_global_iter: int) -> dict[str, object]:
    case = load_case_lu()
    cfg = build_phase2_config(
        case,
        enable_vorab_march_staged=True,
        enable_vorab_cold_budget_march=bool(
            spec.overrides.get("vorab_bed_temperature_march_thesis", False)
        ),
    )
    _apply_baseline(cfg)
    for key, value in spec.overrides.items():
        setattr(cfg, key, value)

    solve_kwargs = phase2_solve_kwargs()
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
    bed0 = bed_cells[0] if bed_cells else profile[0]
    bed_top = bed_cells[-1] if bed_cells else profile[-1]
    t_peak = max(profile, key=lambda p: p["T_K"])

    exit_dry = out.get("exit_gas_dry") or {}
    co = float(exit_dry.get("CO", 0.0))
    co2 = float(exit_dry.get("CO2", 0.0))
    score_exit = _score_exit(co, co2)
    t_rmse_bed = _profile_temperature_rmse(profile, hamel, scope="bed")
    t_rmse_full = _profile_temperature_rmse(profile, hamel, scope="full")
    score_total = score_exit + 0.002 * t_rmse_full

    return {
        "label": spec.label,
        "overrides": dict(spec.overrides),
        "max_global_iter": int(max_global_iter),
        "wall_time_s": float(wall_s),
        "converged_fully": bool(out.get("converged_fully", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", 0.0) or 0.0),
        "nr_outer_iters": int(out.get("nr_outer_iters", 0)),
        "reactor_exit_gas_dry": {"CO": co, "CO2": co2},
        "score_exit_co_co2": float(score_exit),
        "score_total": float(score_total),
        "bed0_T_K": float(bed0["T_K"]),
        "bed_top_T_K": float(bed_top["T_K"]),
        "T_peak_K": float(t_peak["T_K"]),
        "T_peak_xi": float(t_peak["xi_reactor"]),
        "T_peak_cell_type": str(t_peak["cell_type"]),
        "T_rmse_bed_vs_hamel_K": float(t_rmse_bed),
        "T_rmse_full_vs_hamel_K": float(t_rmse_full),
        "ref_dry": {"CO": REF_CO, "CO2": REF_CO2},
        "axial_T_summary": [
            {
                "xi": p["xi_reactor"],
                "T_K": p["T_K"],
                "cell_type": p["cell_type"],
            }
            for p in profile
            if p["cell_type"] == "bed" or abs(p["xi_reactor"] - 0.6) < 0.08 or p["xi_reactor"] > 0.95
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=40)
    parser.add_argument("--labels", nargs="*", default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_temperature_chain_lu.json",
    )
    args = parser.parse_args()

    cases = DEFAULT_CASES
    if args.labels:
        wanted = set(args.labels)
        cases = tuple(c for c in DEFAULT_CASES if c.label in wanted)
        if not cases:
            raise SystemExit(f"no matching labels in {wanted}")

    print(
        f"Phase2 temperature chain | fixed50 r5=0.25 vm_zone=0 | max_global_iter={args.max_global_iter}",
        flush=True,
    )
    print(f"ref CO={REF_CO} CO2={REF_CO2}", flush=True)
    print(f"{'label':<22} {'conv':>4} {'rms':>7} {'CO':>6} {'CO2':>6} {'bed0':>5} {'Tpk':>5} {'xi_pk':>5} {'T_bed':>6} {'T_all':>6} {'score':>6}", flush=True)
    print("=" * 105, flush=True)

    results: list[dict[str, object]] = []
    for spec in cases:
        row = _run_case(spec, max_global_iter=int(args.max_global_iter))
        results.append(row)
        print(
            f"{row['label']:<22} "
            f"{'Y' if row['converged_fully'] else 'N':>4} "
            f"{row['rms_scaled_final']:>7.4f} "
            f"{row['reactor_exit_gas_dry']['CO']:>6.3f} "
            f"{row['reactor_exit_gas_dry']['CO2']:>6.3f} "
            f"{row['bed0_T_K']:>5.0f} "
            f"{row['T_peak_K']:>5.0f} "
            f"{row['T_peak_xi']:>5.2f} "
            f"{row['T_rmse_bed_vs_hamel_K']:>6.1f} "
            f"{row['T_rmse_full_vs_hamel_K']:>6.1f} "
            f"{row['score_total']:>6.3f}",
            flush=True,
        )

    converged = [r for r in results if r["converged_fully"]]
    best_t = min(converged, key=lambda r: r["T_rmse_full_vs_hamel_K"]) if converged else None
    best_score = min(converged, key=lambda r: r["score_total"]) if converged else None

    payload = {
        "note": "temperature chain on fixed50+r5=0.25+vm_zone=0 baseline",
        "max_global_iter": int(args.max_global_iter),
        "best_converged_by_T_rmse_full": None if best_t is None else best_t["label"],
        "best_converged_by_score_total": None if best_score is None else best_score["label"],
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("=" * 105, flush=True)
    if best_t:
        print(f"best T_rmse (full): {best_t['label']} ({best_t['T_rmse_full_vs_hamel_K']:.1f} K)", flush=True)
    if best_score:
        print(f"best score: {best_score['label']} ({best_score['score_total']:.3f})", flush=True)
    print(f"Wrote {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
