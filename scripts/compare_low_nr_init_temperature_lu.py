#!/usr/bin/env python3
"""对比 NR 低温初温（800–900 K）与温度围栏是否剪掉低温解。"""

from __future__ import annotations

from _lu_harness import (
    apply_baseline_kinetics,
    bed_cell_indices,
    build_phase2_config,
    co_co2_score,
    data_path,
    ensure_repo_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from src.core.reactor import Reactor

ensure_repo_path()


def _apply_low_t(cfg, *, label: str) -> None:
    if label == "baseline_1180":
        return
    if label == "adiabatic_cap900":
        cfg.vorab_bottom_temperature_cap_K_thesis = 900.0
        cfg.vorab_bottom_temperature_floor_K_thesis = 720.0
        cfg.nr_temperature_fence_lower_margin_K_thesis = 80.0
        cfg.vorab_nr_temperature_seed_mode_thesis = "adiabatic_anchor"
        return
    if label == "seed850_progress":
        cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
        cfg.vorab_nr_temperature_seed_bottom_K_thesis = 850.0
        cfg.vorab_nr_temperature_seed_top_K_thesis = 1000.0
        cfg.vorab_bottom_temperature_cap_K_thesis = 1000.0
        cfg.vorab_bottom_temperature_floor_K_thesis = 720.0
        cfg.nr_temperature_fence_lower_margin_K_thesis = 100.0
        return
    if label == "seed750_progress":
        cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
        cfg.vorab_nr_temperature_seed_bottom_K_thesis = 750.0
        cfg.vorab_nr_temperature_seed_top_K_thesis = 980.0
        cfg.vorab_bottom_temperature_cap_K_thesis = 980.0
        cfg.vorab_bottom_temperature_floor_K_thesis = 650.0
        cfg.nr_temperature_fence_lower_margin_K_thesis = 90.0
        cfg.vorab_major_gibbs_min_temperature_K_thesis = 700.0
        return
    if label == "seed700_gibbs850":
        cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
        cfg.vorab_nr_temperature_seed_bottom_K_thesis = 700.0
        cfg.vorab_nr_temperature_seed_top_K_thesis = 920.0
        cfg.vorab_bottom_temperature_cap_K_thesis = 920.0
        cfg.vorab_bottom_temperature_floor_K_thesis = 600.0
        cfg.nr_temperature_fence_lower_margin_K_thesis = 80.0
        cfg.vorab_major_gibbs_min_temperature_K_thesis = 850.0
        return
    raise ValueError(label)


def run_case(label: str, *, max_global_iter: int) -> dict:
    cfg = build_phase2_config(enable_vorab_march_staged=True)
    apply_baseline_kinetics(cfg)
    _apply_low_t(cfg, label=label)

    reactor = Reactor(cfg)
    init = run_phase2_init(reactor)
    bed_idx = bed_cell_indices(reactor)
    bed0 = reactor.cells[bed_idx[0]]
    t_min = float(getattr(bed0, "_nr_temperature_min_K", float("nan")))
    t_max = float(getattr(bed0, "_nr_temperature_max_K", float("nan")))

    out, wall_s = timed_phase2_solve(
        reactor,
        **phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False),
    )

    exit_y = out.get("exit_gas_dry") or {}
    co = float(exit_y.get("CO", 0.0))
    co2 = float(exit_y.get("CO2", 0.0))
    bed_T_init = [float(init.T_est[i]) for i in bed_idx[:3]]
    bed_T_sol = [float(reactor.cells[i].T) for i in bed_idx[:3]]

    return {
        "label": label,
        "seed_mode": str(cfg.vorab_nr_temperature_seed_mode_thesis),
        "cap_K": float(cfg.vorab_bottom_temperature_cap_K_thesis),
        "floor_K": float(cfg.vorab_bottom_temperature_floor_K_thesis),
        "fence_lower_margin_K": float(cfg.nr_temperature_fence_lower_margin_K_thesis),
        "bed0_T_init_K": bed_T_init[0] if bed_T_init else float("nan"),
        "bed0_T_fence_min_K": t_min,
        "bed0_T_fence_max_K": t_max,
        "bed_T_init_K": bed_T_init,
        "bed_T_solved_K": bed_T_sol,
        "converged_fully": bool(out.get("converged_fully", False)),
        "rms_scaled_final": float(out.get("rms_scaled_final", float("nan"))),
        "co_exit": co,
        "co2_exit": co2,
        "score_co_co2": co_co2_score(co, co2),
        "wall_time_s": wall_s,
    }


def main() -> None:
    max_iter = 25
    rows = []
    for label in (
        "seed850_progress",
        "seed750_progress",
        "seed700_progress",
        "seed700_blend",
    ):
        print(f"=== {label} ===", flush=True)
        row = run_case(label, max_global_iter=max_iter)
        rows.append(row)
        print(
            f"init bed0/1/2={[round(x) for x in row['bed_T_init_K']]} "
            f"fence=[{row['bed0_T_fence_min_K']:.0f},{row['bed0_T_fence_max_K']:.0f}]K "
            f"solved={[round(x) for x in row['bed_T_solved_K']]} "
            f"conv={row['converged_fully']} rms={row['rms_scaled_final']:.4g} "
            f"CO={row['co_exit']:.4f} CO2={row['co2_exit']:.4f} score={row['score_co_co2']:.3f}",
            flush=True,
        )

    path = data_path("phase2_low_nr_init_seed700_compare.json")
    write_json(path, rows, ensure_ascii=False)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
