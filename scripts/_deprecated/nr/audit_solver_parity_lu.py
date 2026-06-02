"""LU 工况 GS vs global_nr 同参对照审计。

目标：在完全一致的基线参数下，对比两条求解路径的收敛与 KPI。

用法：
    /Users/liuzhen/AI-projects/bfb-gasifier/.venv/bin/python scripts/audit_solver_parity_lu.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_reactor_config,
    build_phase1_htw_lu_global_nr_reactor_config,
    load_validation_case_node,
)


def _rel_err(sim: float, ref_val: float) -> float:
    return abs(sim - ref_val) / max(abs(ref_val), 1e-12)


def run_solver(solver: str, enable_r12: bool = False) -> dict:
    cfg = (
        build_phase1_htw_lu_global_nr_reactor_config()
        if solver == "global_nr"
        else build_phase1_htw_lu_reactor_config()
    )
    if solver == "gauss_seidel":
        cfg.allow_legacy_gs = True
    cfg.enable_r12 = enable_r12

    reactor = Reactor(cfg)
    solve_kwargs = (
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
        if solver == "global_nr"
        else dict(PHASE1_HTW_LU_SOLVE_KWARGS)
    )
    solve_kwargs["solver"] = solver
    result, monitor = solve_with_nr_monitor(
        reactor,
        {**solve_kwargs, "verbose": False},
        check_x0=(solver == "global_nr"),
    )

    y_dry = result["exit_gas_dry"]
    return {
        "solver": solver,
        "enable_r12": enable_r12,
        "converged": bool(result.get("converged", False)),
        "converged_fully": result.get("converged_fully"),
        "nr_init_strategy": result.get("nr_init_strategy"),
        "nr_init_s_total": result.get("nr_init_s_total"),
        "nr_vorabrechnung_s": result.get("nr_vorabrechnung_s"),
        "nr_gs_warmup_steps": result.get("nr_gs_warmup_steps"),
        "nr_jacobian_strategy": result.get("nr_jacobian_strategy"),
        "n_iter": int(result.get("n_iter", -1)),
        "residual": float(result.get("residual", 0.0)),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "nr_total_s": result.get("nr_total_s"),
        "nr_gs_warmup_s": result.get("nr_gs_warmup_s"),
        "nr_inner_solve_s_total": result.get("nr_inner_solve_s_total"),
        "nr_jacobian_build_s": (result.get("nr_timing") or {}).get("jacobian_build_s"),
        "nr_global_residual_calls": (result.get("nr_counts") or {}).get("global_residual_calls"),
        "nr_jacobian_cell_residual_calls": (result.get("nr_counts") or {}).get("jacobian_cell_residual_calls"),
        "T_exit": float(result["T_profile"][-1]),
        "carbon_conv": float(result["carbon_conv"]),
        "CO_dry": float(y_dry.get("CO", 0.0)),
        "CO2_dry": float(y_dry.get("CO2", 0.0)),
        "H2_dry": float(y_dry.get("H2", 0.0)),
        "CH4_dry": float(y_dry.get("CH4", 0.0)),
        "wall_s": float(monitor.get("wall_time_s", float("nan"))),
        "x0_ok": (monitor.get("x0_sanity") or {}).get("ok"),
    }


def main() -> int:
    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    ref_dry = ref.get("exit_gas_dry_mol_frac", {})
    t_ref = float(ref.get("exit_temperature_K", 1100.0))
    xc_ref = float(ref.get("carbon_conversion_pct", 95.0)) / 100.0

    print("=" * 80)
    print("LU solver parity audit (GS baseline vs shared global NR dev config)")
    print("=" * 80)
    print("Shared assumptions: n_cells=10, recirculation=0.1, heat_loss=0.1, max_iter=20, tol=1.0")
    print("Mechanism baseline: enable_r12=False")
    print()

    gs = run_solver("gauss_seidel", enable_r12=False)
    nr = run_solver("global_nr", enable_r12=False)

    for row in (gs, nr):
        name = row["solver"]
        print(f"[{name}]")
        print(
            f"  converged={row['converged']} converged_fully={row['converged_fully']} "
            f"n_iter={row['n_iter']} residual={row['residual']:.3e}"
        )
        print(f"  wall_s={row['wall_s']:.2f}s x0_ok={row['x0_ok']}")
        if row["nr_init_strategy"] is not None:
            print(f"  nr_init_strategy={row['nr_init_strategy']}")
        if row["nr_gs_warmup_steps"] is not None:
            print(f"  nr_gs_warmup_steps={row['nr_gs_warmup_steps']}")
        if row["nr_jacobian_strategy"] is not None:
            print(f"  nr_jacobian_strategy={row['nr_jacobian_strategy']}")
        if row["rms_scaled_final"] is not None:
            print(f"  rms_scaled_final={row['rms_scaled_final']:.3e}")
        if row["nr_total_s"] is not None:
            print(
                f"  nr_init_s_total={row['nr_init_s_total']:.2f}s "
                f"vorabrechnung_s={row['nr_vorabrechnung_s']:.2f}s "
                f"  nr_total_s={row['nr_total_s']:.2f}s "
                f"gs_warmup_s={row['nr_gs_warmup_s']:.2f}s "
                f"inner_solve_s_total={row['nr_inner_solve_s_total']:.2f}s "
                f"jacobian_build_s={row['nr_jacobian_build_s']:.2f}s "
                f"global_residual_calls={row['nr_global_residual_calls']} "
                f"jacobian_cell_residual_calls={row['nr_jacobian_cell_residual_calls']}"
            )
        print(
            f"  T_exit={row['T_exit']:.1f} K | carbon_conv={row['carbon_conv']:.4f} | "
            f"CO={row['CO_dry']:.4f} CO2={row['CO2_dry']:.4f} H2={row['H2_dry']:.4f} CH4={row['CH4_dry']:.4f}"
        )
        print(
            "  rel_err: "
            f"T={_rel_err(row['T_exit'], t_ref):.2%}, "
            f"Xc={_rel_err(row['carbon_conv'], xc_ref):.2%}, "
            f"CO={_rel_err(row['CO_dry'], float(ref_dry.get('CO', 0.13))):.2%}, "
            f"CO2={_rel_err(row['CO2_dry'], float(ref_dry.get('CO2', 0.11))):.2%}, "
            f"H2={_rel_err(row['H2_dry'], float(ref_dry.get('H2', 0.12))):.2%}, "
            f"CH4={_rel_err(row['CH4_dry'], float(ref_dry.get('CH4', 0.028))):.2%}"
        )
        print()

    print("Decision gate:")
    print("- 若 NR 的 rms_scaled_final 明显下降且 KPI 优于 GS：优先深挖 NR 路径。")
    print("- 若两者同样失真：优先检查 reactor 耦合/边界传播与 GS 预叠加策略。")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
