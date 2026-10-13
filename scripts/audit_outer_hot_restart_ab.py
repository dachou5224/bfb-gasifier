#!/usr/bin/env python3
"""同一次 solve 内 CO 带宽后 outer 热重启 A/B（handoff §5.14）。

对照：
  A  ps_f0.3 o8（基线，bed2 翻）
  B  ps_f0.3 o8 + 嵌套 skip_precalc（续算前含 freeboard T 投影）
  C  ps_f0.3 o16 无热重启（预算对照）
  D  A 后 ``_solve_global_nr(skip_precalc=True)`` 再 o8（oracle）

用法：
  python3 scripts/audit_outer_hot_restart_ab.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
from time import perf_counter
from typing import Any

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.species import GAS_SPECIES_INDEX

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}

PS = {
    "nr_outer_preinner_upper_co_local2x2_thesis": True,
    "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
    "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
    "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.3,
    "nr_outer_preinner_upper_dense_absorb_mode_thesis": "phase_share",
    "nr_outer_preinner_upper_dense_absorb_post_gate_thesis": True,
}


def _bed2_R(reactor: Any) -> dict[str, float]:
    c = reactor.cells[2]
    c.residuals()
    out: dict[str, float] = {}
    for sp in ("H2O", "H2", "O2", "CO"):
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"R_{sp}"] = float(c.R_gas_d[j] + c.R_gas_b[j])
    out["T"] = float(c.T)
    out["N_O2"] = float(
        c.N_d[int(GAS_SPECIES_INDEX["O2"])] + c.N_b[int(GAS_SPECIES_INDEX["O2"])]
    )
    return out


def _summarize(label: str, reactor: Any, result: dict[str, Any], wall_s: float) -> dict[str, Any]:
    exit_gas = dict(result.get("bed_exit_gas") or {})
    hist = list(result.get("nr_outer_history") or [])
    restarted = [h for h in hist if h.get("outer_hot_restarted")]
    hr = dict(result.get("nr_outer_hot_restart") or {})
    bed2 = _bed2_R(reactor)
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "converged": bool(result.get("converged", False)),
        "n_outer": len(hist),
        "hot_restart_meta": hr,
        "hot_restart_outer_n": len(restarted),
        "bed2": bed2,
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} split={out['split']:.4f} "
        f"bed2_RH2O={bed2['R_H2O']:+.2f} hr={hr.get('triggered')} "
        f"n_outer={out['n_outer']} wall={wall_s:.0f}s"
    )
    return out


def _run_single(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    out = _summarize(label, reactor, result, wall_s)
    out["overrides"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in overrides.items()}
    out["max_global_iter"] = int(max_iter)
    return out


def _run_twostage_oracle(max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    skw = phase2_solve_kwargs(max_global_iter=int(max_iter))
    result_a, wall_a = timed_phase2_solve(reactor, **skw)
    assert isinstance(result_a, dict)
    stage_a = _summarize("D1_ps_o8", reactor, result_a, wall_a)
    # oracle：直接 ``_solve_global_nr(skip_precalc=True)``（solve() 不暴露该参数）
    t0 = perf_counter()
    result_b = reactor._solve_global_nr(
        max_iter=int(max_iter),
        tol=float(skw.get("tol_global", 1.0)),
        verbose=False,
        init_strategy=skw.get("nr_init_strategy"),
        jacobian_strategy=skw.get("nr_jacobian_strategy"),
        linear_solver_backend=skw.get("nr_linear_solver_backend"),
        jacobian_lag_steps=skw.get("nr_jacobian_lag_steps"),
        skip_precalc=True,
    )
    wall_b = perf_counter() - t0
    assert isinstance(result_b, dict)
    stage_b = _summarize("D2_skip_precalc_o8", reactor, result_b, wall_b)
    return {
        "label": "D_twostage_skip_precalc",
        "stage_a": stage_a,
        "stage_b": stage_b,
        "exit_CO": stage_b["exit_CO"],
        "rms": stage_b["rms"],
        "split": stage_b["split"],
        "energy": stage_b["energy"],
        "converged": stage_b["converged"],
        "bed2": stage_b["bed2"],
        "wall_s": float(stage_a["wall_s"] + stage_b["wall_s"]),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    ap.add_argument(
        "--output",
        type=str,
        default=str(data_path("phase2_outer_hot_restart_ab_o8.json")),
    )
    args = ap.parse_args()
    n = int(args.max_global_iter)
    bonus_inner = max(5 * n, 20)

    arms: list[dict[str, Any]] = []
    arms.append(_run_single("A_ps_o8", {}, n))
    arms.append(
        _run_single(
            "B_hot_restart_nested_skip_precalc",
            {
                "nr_outer_hot_restart_after_y_co_thesis": True,
                "nr_outer_hot_restart_y_co_threshold_thesis": 0.14,
                "nr_outer_hot_restart_bonus_inner_budget_thesis": bonus_inner,
            },
            n,
        )
    )
    arms.append(_run_single("C_ps_o16_no_restart", {}, 2 * n))
    arms.append(_run_twostage_oracle(n))

    payload = {
        "note": (
            "单次 solve 内 CO 带宽后嵌套 ``_solve_global_nr(skip_precalc=True)``；"
            "中途门禁复位与同闭包二次 `_run_outer` 均证伪；D 为显式两段 oracle"
        ),
        "isomorph": ISOMORPH,
        "ps": {k: (list(v) if isinstance(v, tuple) else v) for k, v in PS.items()},
        "bonus_inner": bonus_inner,
        "falsified_paths": {
            "mid_loop_gate_reset": {
                "note": "CO 带宽后中途复位 refresh/merit，n_outer=20 仍 bed2≈+65",
                "CO": 0.1562,
                "rms": 0.0587,
                "bed2_R_H2O": 65.62,
            },
            "same_closure_second_abgleich": {
                "note": "同 `_solve_inner` 闭包二次 `_run_outer`，n_outer=20 仍 bed2≈+65",
                "CO": 0.1559,
                "rms": 0.0587,
                "bed2_R_H2O": 65.62,
            },
            "nested_without_fb_projection": {
                "note": "嵌套 skip_precalc 但无续算前 freeboard T 投影 → bed2 仍 +65",
                "CO": 0.1559,
                "rms": 0.0587,
                "bed2_R_H2O": 65.62,
            },
        },
        "arms": arms,
    }
    path = write_json(args.output, payload)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
