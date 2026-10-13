#!/usr/bin/env python3
"""温度围栏下界 + freeboard secondary 影响审计（LU Phase2）。

验证：
- 相对围栏 ``t_ref - 150 K`` 是否把 bed0 钉死在 1030 K
- 绝对下界 700 K 是否允许床底继续下探
- secondary air（~0.54% 总 O2）对 ξ≈0.6 温度峰与出口组成的影响

Usage
-----
    python3 scripts/audit_temperature_fence_secondary_lu.py --max-global-iter 40
"""

from __future__ import annotations

from _lu_harness import (
    REF_CO,
    REF_CO2,
    LADDER_P2,
    add_json_output_arg,
    add_max_global_iter_arg,
    apply_config_overrides,
    build_phase2_config,
    load_case_lu,
    phase2_solve_kwargs,
    repo_root,
    run_phase2_init,
    write_json,
)

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
class AuditCase:
    label: str
    overrides: dict[str, object] = field(default_factory=dict)
    disable_secondary: bool = False


def _score(co: float, co2: float) -> float:
    return abs(co - REF_CO) / REF_CO + abs(co2 - REF_CO2) / REF_CO2


def _apply_baseline(cfg) -> None:
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 0
    cfg.vorab_major_gibbs_x0 = False
    cfg.r5_scale = 0.25
    cfg.gas_inlet_split_strategy = "fixed"
    cfg.gas_inlet_dense_frac = 0.50
    cfg.vorab_staged_primary_gas_dense_frac_thesis = 0.50
    cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
    cfg.vorab_nr_temperature_seed_bottom_K_thesis = 1020.0
    cfg.vorab_nr_temperature_seed_top_K_thesis = 1100.0
    for key, value in {
        "nr_line_search_gas_phase_split_merit_thesis": True,
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1180.0,
        "freeboard_trajectory_coeff_model": "exact_hamel",
    }.items():
        setattr(cfg, key, value)


def _fence_bounds(reactor: Reactor, bed_idx: int) -> tuple[float, float, float]:
    cell = reactor.cells[bed_idx]
    t_ref = float(getattr(cell, "_vorab_temperature_estimate_K", float("nan")))
    t_lo = float(getattr(cell, "_nr_temperature_min_K", float("nan")))
    t_hi = float(getattr(cell, "_nr_temperature_max_K", float("nan")))
    return t_ref, t_lo, t_hi


def _run(spec: AuditCase, *, max_global_iter: int) -> dict[str, object]:
    case = load_case_lu()
    cfg = build_phase2_config(case)
    _apply_baseline(cfg)
    for key, value in spec.overrides.items():
        setattr(cfg, key, value)
    if spec.disable_secondary:
        total_o2 = float(cfg.O2_feed + cfg.freeboard_secondary_O2_mol_s)
        total_n2 = float(cfg.N2_feed + cfg.freeboard_secondary_N2_mol_s)
        cfg.O2_feed = total_o2
        cfg.N2_feed = total_n2
        cfg.freeboard_secondary_O2_mol_s = 0.0
        cfg.freeboard_secondary_N2_mol_s = 0.0

    reactor = Reactor(cfg)
    init = run_phase2_init(reactor, init_strategy="vorabrechnung", gs_warmup_steps=None)
    bed_idx = next(i for i, c in enumerate(reactor.cells) if str(getattr(c, "cell_type", "bed")) == "bed")
    t_ref, t_lo, t_hi = _fence_bounds(reactor, bed_idx)
    bed0_init = float(init.T_est[bed_idx])

    kw = phase2_solve_kwargs()
    kw["max_global_iter"] = int(max_global_iter)
    kw["verbose"] = False
    t0 = perf_counter()
    out = reactor.solve(**kw)
    wall_s = perf_counter() - t0

    h_total = float(case["reactor_height_m"])
    profile: list[dict[str, object]] = []
    for i, cell in enumerate(reactor.cells):
        xi = float(cell.geo.h_center / h_total)
        y = gas_mole_fractions(cell.N_b, cell.N_d, "mixed")
        profile.append(
            {
                "cell": i,
                "cell_type": str(getattr(cell, "cell_type", "bed")),
                "xi": xi,
                "T_K": float(cell.T),
                "y_O2_wet": float(y[GSI["O2"]]),
            }
        )

    near_sec = min(profile, key=lambda p: abs(float(p["xi"]) - 0.6))
    t_peak = max(profile, key=lambda p: float(p["T_K"]))
    exit_dry = out.get("exit_gas_dry") or {}
    co = float(exit_dry.get("CO", 0.0))
    co2 = float(exit_dry.get("CO2", 0.0))

    return {
        "label": spec.label,
        "overrides": dict(spec.overrides),
        "disable_secondary": bool(spec.disable_secondary),
        "secondary_O2_mol_s": float(cfg.freeboard_secondary_O2_mol_s),
        "secondary_frac_of_total_O2_pct": 100.0
        * float(cfg.freeboard_secondary_O2_mol_s)
        / max(float(cfg.O2_feed + cfg.freeboard_secondary_O2_mol_s), 1e-12),
        "bed0_T_init_K": bed0_init,
        "bed0_fence_t_ref_K": t_ref,
        "bed0_fence_min_K": t_lo,
        "bed0_fence_max_K": t_hi,
        "bed0_T_solved_K": float(reactor.cells[bed_idx].T),
        "converged_fully": bool(out.get("converged_fully", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", 0.0) or 0.0),
        "reactor_exit_gas_dry": {"CO": co, "CO2": co2},
        "score_co_co2": _score(co, co2),
        "T_peak_K": float(t_peak["T_K"]),
        "T_peak_xi": float(t_peak["xi"]),
        "T_at_xi060_K": float(near_sec["T_K"]),
        "O2_at_xi060_wet": float(near_sec["y_O2_wet"]),
        "freeboard_secondary_injection_applied": bool(out.get("freeboard_secondary_injection_applied", False)),
        "wall_time_s": float(wall_s),
        "axial_T": [{"xi": p["xi"], "T_K": p["T_K"], "cell_type": p["cell_type"]} for p in profile],
    }


DEFAULT_CASES: tuple[AuditCase, ...] = (
    AuditCase("rel_fence_default"),
    AuditCase(
        "abs_fence_700",
        {
            "nr_temperature_fence_lower_mode_thesis": "absolute",
            "nr_temperature_fence_absolute_min_K_thesis": 700.0,
        },
    ),
    AuditCase(
        "rel_margin_480",
        {"nr_temperature_fence_lower_margin_K_thesis": 480.0},
    ),
    AuditCase("no_secondary", disable_secondary=True),
    AuditCase(
        "secondary_refine4",
        {
            "freeboard_secondary_local_refine": 4,
            "freeboard_secondary_injection_mode": "distributed_uniform",
        },
    ),
    AuditCase(
        "abs700_no_secondary",
        {
            "nr_temperature_fence_lower_mode_thesis": "absolute",
            "nr_temperature_fence_absolute_min_K_thesis": 700.0,
        },
        disable_secondary=True,
    ),
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=40)
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "data/phase2_fence_secondary_audit.json")
    args = parser.parse_args()

    case = load_case_lu()
    print("=== LU oxidant budget ===")
    print(f"secondary air: {case['secondary_air_Nm3_h']} Nm3/h @ xi={case['secondary_injection_xi']}")
    print(f"total air:     {case.get('total_air_Nm3_h')} Nm3/h")
    print(f"secondary O2:  {case['secondary_agent_O2_mol_s']:.4f} mol/s")
    print()
    print(f"{'label':<22} {'conv':>4} {'bed0':>5} {'f_min':>5} {'f_ref':>5} {'T@.6':>5} {'Tpk':>5} {'CO':>6} {'CO2':>6}")
    print("=" * 95)

    rows = [_run(spec, max_global_iter=int(args.max_global_iter)) for spec in DEFAULT_CASES]
    for row in rows:
        ex = row["reactor_exit_gas_dry"]
        print(
            f"{row['label']:<22} "
            f"{'Y' if row['converged_fully'] else 'N':>4} "
            f"{row['bed0_T_solved_K']:>5.0f} "
            f"{row['bed0_fence_min_K']:>5.0f} "
            f"{row['bed0_fence_t_ref_K']:>5.0f} "
            f"{row['T_at_xi060_K']:>5.0f} "
            f"{row['T_peak_K']:>5.0f} "
            f"{ex['CO']:>6.3f} "
            f"{ex['CO2']:>6.3f}"
        )

    payload = {
        "note": "temperature fence lower bound vs secondary air audit; baseline seed_hamel_1020_1100",
        "hamel_T_xi060_K": 1240.0,
        "max_global_iter": int(args.max_global_iter),
        "results": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("=" * 95)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
