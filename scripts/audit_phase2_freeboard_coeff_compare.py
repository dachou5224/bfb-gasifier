from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_validation_case_node,
)


def _run_mode(
    mode: str,
    *,
    n_bed: int,
    n_freeboard: int,
    max_global_iter: int,
    trajectory_model: str,
    closure_char_hetero: bool,
    explicit_char_hetero: bool,
) -> dict:
    cfg = build_phase2_htw_lu_freeboard_reactor_config(n_freeboard_cells=n_freeboard)
    cfg.n_cells = int(n_bed)
    cfg.n_freeboard_cells = int(n_freeboard)
    cfg.vorab_major_gibbs_x0 = False
    cfg.freeboard_trajectory_model = str(trajectory_model)
    cfg.freeboard_trajectory_coeff_model = str(mode)
    cfg.freeboard_closure_enable_char_hetero = bool(closure_char_hetero)
    cfg.freeboard_explicit_enable_char_hetero = bool(explicit_char_hetero)

    t0 = time.perf_counter()
    reactor = Reactor(cfg)
    result = reactor.solve(
        solver="global_nr",
        max_global_iter=int(max_global_iter),
        tol_global=1.0,
        nr_init_strategy="vorabrechnung",
        nr_jacobian_strategy=None,
    )
    wall_s = time.perf_counter() - t0
    dry = result["exit_gas_dry"]
    return {
        "mode": mode,
        "trajectory_model": str(trajectory_model),
        "closure_char_hetero": bool(closure_char_hetero),
        "explicit_char_hetero": bool(explicit_char_hetero),
        "wall_s": float(wall_s),
        "converged": bool(result.get("converged")),
        "converged_outer": bool(result.get("converged_outer")),
        "converged_fully": bool(result.get("converged_fully")),
        "rms_scaled_final": float(result.get("rms_scaled_final") or 0.0),
        "n_iter": int(result.get("n_iter") or 0),
        "nr_outer_iters": int(result.get("nr_outer_iters") or 0),
        "nr_linear_solver_backend": result.get("nr_linear_solver_backend"),
        "nr_jacobian_strategy": result.get("nr_jacobian_strategy"),
        "nr_timing": dict(result.get("nr_timing") or {}),
        "nr_counts": dict(result.get("nr_counts") or {}),
        "reactor_exit_T": float(result["reactor_exit_T"]),
        "carbon_conv": float(result["carbon_conv"]),
        "dry_gas": {sp: float(dry.get(sp, 0.0)) for sp in ("CO", "CO2", "H2", "CH4", "O2", "N2")},
        "freeboard_char_holdup_kg": [float(v) for v in result["freeboard_solid_holdup_char_profile_kg"]],
        "freeboard_ash_holdup_kg": [float(v) for v in result["freeboard_solid_holdup_ash_profile_kg"]],
        "freeboard_char_up_kg_s": [float(v) for v in result["freeboard_entrained_char_kg_s"]],
        "freeboard_char_down_kg_s": [float(v) for v in result["freeboard_entrained_return_char_profile_kg_s"]],
        "freeboard_char_rxn_mol_s": [float(v) for v in result["freeboard_char_reaction_source_profile_mol_s"]],
        "freeboard_closure_r1": [float(v) for v in result["freeboard_reaction_diag_impl"]["R1"]],
        "freeboard_closure_r2": [float(v) for v in result["freeboard_reaction_diag_impl"]["R2"]],
        "freeboard_explicit_r1": [float(v) for v in result["freeboard_explicit_reaction_diag_impl"]["R1"]],
        "freeboard_explicit_r2": [float(v) for v in result["freeboard_explicit_reaction_diag_impl"]["R2"]],
        "top_freeboard_char_up_kg_s": float(result["freeboard_entrained_char_kg_s"][-1]) if result["freeboard_entrained_char_kg_s"] else 0.0,
        "cyclone_in_char_kg_s": float(np.sum(np.maximum(reactor.cyclone_cell.m_solid_auf_in[:, 0], 0.0))) if reactor.cyclone_cell is not None else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare freeboard coefficient models on a lightweight phase-2 reactor case.")
    parser.add_argument("--n-bed", type=int, default=3)
    parser.add_argument("--n-freeboard", type=int, default=2)
    parser.add_argument("--max-global-iter", type=int, default=2)
    parser.add_argument("--json", action="store_true")
    parser.add_argument(
        "--trajectory-models",
        type=str,
        default="analytical_wirsum",
        help="Comma-separated freeboard trajectory models to compare (e.g. analytical_wirsum,global_projected_hamel)",
    )
    parser.add_argument(
        "--coeff-modes",
        type=str,
        default="stable_mixed_drag_split,exact_hamel",
        help="Comma-separated freeboard trajectory coefficient modes to compare.",
    )
    parser.add_argument(
        "--closure-char-hetero",
        type=str,
        default="off",
        help="Comma-separated closure char-hetero toggles: on,off",
    )
    parser.add_argument(
        "--explicit-char-hetero",
        type=str,
        default="on",
        help="Comma-separated explicit freeboard-cell char-hetero toggles: on,off",
    )
    args = parser.parse_args()

    trajectory_models = [item.strip() for item in str(args.trajectory_models).split(",") if item.strip()]
    coeff_modes = [item.strip() for item in str(args.coeff_modes).split(",") if item.strip()]
    closure_modes = [item.strip() for item in str(args.closure_char_hetero).split(",") if item.strip()]
    explicit_modes = [item.strip() for item in str(args.explicit_char_hetero).split(",") if item.strip()]

    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    report = {
        "config": {
            "n_bed": int(args.n_bed),
            "n_freeboard": int(args.n_freeboard),
            "max_global_iter": int(args.max_global_iter),
            "trajectory_models": trajectory_models,
            "coeff_modes": coeff_modes,
            "closure_char_hetero": closure_modes,
            "explicit_char_hetero": explicit_modes,
        },
        "reference": {
            "exit_temperature_K": ref.get("exit_temperature_K"),
            "carbon_conversion_pct": ref.get("carbon_conversion_pct"),
            "dry_gas_composition": ref.get("dry_gas_composition", {}),
        },
        "modes": [
            _run_mode(
                coeff_mode,
                n_bed=args.n_bed,
                n_freeboard=args.n_freeboard,
                max_global_iter=args.max_global_iter,
                trajectory_model=trajectory_model,
                closure_char_hetero=(closure_mode == "on"),
                explicit_char_hetero=(explicit_mode == "on"),
            )
            for trajectory_model in trajectory_models
            for coeff_mode in coeff_modes
            for closure_mode in closure_modes
            for explicit_mode in explicit_modes
        ],
    }

    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return

    print("Reference:")
    print(f"  T_exit={report['reference']['exit_temperature_K']} K | carbon_conversion={report['reference']['carbon_conversion_pct']}%")
    print(f"  dry_gas={report['reference']['dry_gas_composition']}")
    print(f"  config={report['config']}")
    print()
    for row in report["modes"]:
        jac_s = float(row["nr_timing"].get("jacobian_build_s", 0.0))
        lin_s = float(row["nr_timing"].get("linear_solve_s", 0.0))
        jac_frac = jac_s / max(float(row["wall_s"]), 1e-12)
        bottleneck = "jacobian_build" if jac_frac >= 0.6 else ("linear_solve" if lin_s > jac_s else "mixed")
        print(
            f"Mode: trajectory={row['trajectory_model']} coeff={row['mode']} "
            f"closure_char_hetero={row['closure_char_hetero']} "
            f"explicit_char_hetero={row['explicit_char_hetero']}"
        )
        print(
            f"  converged={row['converged']} outer={row['converged_outer']} fully={row['converged_fully']} "
            f"rms={row['rms_scaled_final']:.6f} wall_s={row['wall_s']:.2f} "
            f"n_iter={row['n_iter']} outer_iters={row['nr_outer_iters']}"
        )
        print(
            f"  nr_jacobian={row['nr_jacobian_strategy']} "
            f"linear={row['nr_linear_solver_backend']} "
            f"jac_s={jac_s:.2f} "
            f"lin_s={lin_s:.2f} "
            f"bottleneck={bottleneck}"
        )
        print(f"  T_exit={row['reactor_exit_T']:.3f} K | carbon_conv={row['carbon_conv']:.6f}")
        print(f"  dry_gas={row['dry_gas']}")
        print(f"  top_freeboard_char_up_kg_s={row['top_freeboard_char_up_kg_s']:.6f} | cyclone_in_char_kg_s={row['cyclone_in_char_kg_s']:.6f}")
        print(f"  freeboard_char_holdup_kg={row['freeboard_char_holdup_kg']}")
        print(f"  freeboard_char_up_kg_s={row['freeboard_char_up_kg_s']}")
        print(f"  freeboard_char_down_kg_s={row['freeboard_char_down_kg_s']}")
        print(f"  freeboard_char_rxn_mol_s={row['freeboard_char_rxn_mol_s']}")
        print(f"  freeboard_closure_R1={row['freeboard_closure_r1']}")
        print(f"  freeboard_closure_R2={row['freeboard_closure_r2']}")
        print(f"  freeboard_explicit_R1={row['freeboard_explicit_r1']}")
        print(f"  freeboard_explicit_R2={row['freeboard_explicit_r2']}")
        print()


if __name__ == "__main__":
    main()
