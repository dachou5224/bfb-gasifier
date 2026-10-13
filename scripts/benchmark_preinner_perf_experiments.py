#!/usr/bin/env python3
"""Pre-inner / solid-holdup / NR 性能优化实验矩阵。

固定 Phase2 HTW LU freeboard P2 口径，对比 wall time 与收敛质量。
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter
from typing import Any, Callable

from _lu_harness import (
    LADDER_P2,
    apply_config_overrides,
    build_phase2_config,
    ensure_repo_path,
    phase2_solve_kwargs,
    repo_root,
    run_phase2_init,
)

ensure_repo_path()

from src.core.connectivity import (
    preproject_bed_solid_holdup_to_transport_closure,
    run_bottom_zone_solid_energy_preprojection,
    run_outer_preinner_start_value_bridge,
)
from src.core.reactor import Reactor

REPO_ROOT = repo_root()


@dataclass
class ExperimentPatch:
    name: str
    description: str
    config_overrides: dict[str, object] = field(default_factory=dict)
    connectivity_patches: list[tuple[str, Callable[..., Any]]] = field(default_factory=list)


def _patch_solid_holdup_lsq(*, max_nfev: int, skip_rms_scaled: float | None = None):
    """Wrap solid-holdup LSQ with tunable ``max_nfev`` and optional early skip."""

    def _wrapped(
        cells,
        cfg,
        *,
        cell_index: int,
        apply_bc_fn=None,
        improvement_factor: float = 0.75,
        holdup_upper_scale: float = 3.0,
    ) -> bool:
        import numpy as np
        from src.core.cell import Cell

        if cell_index < 0 or cell_index >= len(cells):
            return False
        cell = cells[cell_index]
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            return False
        if str(getattr(cell, "solid_state_model", "legacy_stream")) != "holdup_transport":
            return False
        if apply_bc_fn is not None:
            apply_bc_fn()
        ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
        original = np.asarray(cell.m_solid, dtype=np.float64).copy()
        n_solid = int(original.size)
        if n_solid <= 0:
            return False
        initial_res = np.asarray(cell.residuals(), dtype=np.float64)
        solid_off = int(initial_res.size - n_solid - 1)
        initial_solid = initial_res[solid_off : solid_off + n_solid]
        initial_rms = float(np.sqrt(np.mean((initial_solid / ref_solid) ** 2)))
        if not np.isfinite(initial_rms) or initial_rms <= 1.0e-12:
            return False
        if skip_rms_scaled is not None and initial_rms <= float(skip_rms_scaled):
            cell._bottom_zone_solid_preprojection_diag = {
                "cell_index": int(cell_index),
                "skipped_well_converged": True,
                "initial_solid_rms_kg_s": float(initial_rms * ref_solid),
            }
            return False

        def _invalidate() -> None:
            cell._thermo_cache_valid = False
            cell._hydro_cache_valid = False

        def _set_holdup(values: np.ndarray) -> None:
            cell.m_solid[:] = np.maximum(np.asarray(values, dtype=np.float64).reshape(original.shape), 1.0e-12)
            _invalidate()

        def _solid_objective(values: np.ndarray) -> np.ndarray:
            _set_holdup(values)
            res = np.asarray(cell.residuals(), dtype=np.float64)
            return (res[solid_off : solid_off + n_solid] / ref_solid).copy()

        x0 = np.maximum(original.reshape(-1), 1.0e-12)
        upper_scale = float(max(holdup_upper_scale, 1.05))
        upper = np.maximum(x0 * upper_scale, 1.0e-6)
        lower = np.full_like(x0, 1.0e-12)
        try:
            from scipy.optimize import least_squares

            sol = least_squares(
                _solid_objective,
                x0=x0,
                bounds=(lower, upper),
                max_nfev=int(max_nfev),
                xtol=1.0e-8,
                ftol=1.0e-8,
                gtol=1.0e-8,
            )
            _set_holdup(sol.x)
            final_res = np.asarray(cell.residuals(), dtype=np.float64)
            final_solid = final_res[solid_off : solid_off + n_solid]
            final_rms = float(np.sqrt(np.mean((final_solid / ref_solid) ** 2)))
        except Exception:
            cell.m_solid[:] = original
            _invalidate()
            if apply_bc_fn is not None:
                apply_bc_fn()
            cell.residuals()
            return False

        if np.isfinite(final_rms) and final_rms <= improvement_factor * initial_rms:
            if apply_bc_fn is not None:
                apply_bc_fn()
            cell._bottom_zone_solid_preprojection_diag = {
                "cell_index": int(cell_index),
                "initial_solid_rms_kg_s": float(initial_rms * ref_solid),
                "final_solid_rms_kg_s": float(final_rms * ref_solid),
            }
            cell.residuals()
            return True

        cell.m_solid[:] = original
        _invalidate()
        if apply_bc_fn is not None:
            apply_bc_fn()
        cell.residuals()
        return False

    return _wrapped


def _patch_solid_energy_single_pass():
    """Skip the second solid-holdup pass over downstream beds."""

    def _wrapped(cells, cfg, *, apply_bc_fn=None):
        import src.core.connectivity as conn

        max_bed_index = int(max(getattr(cfg, "nr_bottom_zone_preprojection_max_bed_index_thesis", 2), 0))
        align_downstream = bool(getattr(cfg, "nr_bottom_zone_align_downstream_temperature_thesis", False))
        bed1_energy_t = bool(getattr(cfg, "nr_bottom_zone_bed1_energy_temperature_preprojection_thesis", True))
        diag: dict[str, bool | float] = {
            "enabled": True,
            "max_bed_index": float(max_bed_index),
            "solid_preprojections_accepted": 0.0,
            "downstream_temperature_aligned": False,
            "bed1_energy_temperature_preprojection_accepted": False,
            "solid_repass_skipped": True,
        }
        if not cells:
            return diag
        if align_downstream:
            diag["downstream_temperature_aligned"] = bool(
                conn.align_downstream_bed_temperatures_to_upstream(
                    cells,
                    cfg,
                    upstream_index=0,
                    downstream_indices=tuple(
                        j
                        for j in range(1, min(max_bed_index + 1, len(cells)))
                        if str(getattr(cells[j], "cell_type", "bed")) == "bed"
                    ),
                )
            )
            if apply_bc_fn is not None:
                apply_bc_fn()
        solid_ok: list[bool] = []
        for i in range(min(max_bed_index + 1, len(cells))):
            solid_ok.append(
                preproject_bed_solid_holdup_to_transport_closure(
                    cells,
                    cfg,
                    cell_index=i,
                    apply_bc_fn=apply_bc_fn,
                )
            )
        if apply_bc_fn is not None:
            apply_bc_fn()
        if (
            bed1_energy_t
            and len(cells) >= 2
            and str(getattr(cells[1], "cell_type", "bed")) == "bed"
            and apply_bc_fn is not None
        ):
            diag["bed1_energy_temperature_preprojection_accepted"] = bool(
                conn.preproject_bed_temperature_to_local_energy_closure(
                    cells,
                    cfg,
                    cell_index=1,
                    apply_bc_fn=apply_bc_fn,
                )
            )
            apply_bc_fn()
        diag["solid_preprojections_accepted"] = float(sum(1 for ok in solid_ok if ok))
        return diag

    return _wrapped


EXPERIMENTS: list[ExperimentPatch] = [
    ExperimentPatch("baseline", "当前默认（LSQ max_nfev=32，双 pass solid）"),
    ExperimentPatch(
        "lsq_nfev_16",
        "solid holdup LSQ max_nfev 32→16",
        connectivity_patches=[
            (
                "preproject_bed_solid_holdup_to_transport_closure",
                _patch_solid_holdup_lsq(max_nfev=16),
            )
        ],
    ),
    ExperimentPatch(
        "lsq_nfev_12",
        "solid holdup LSQ max_nfev 32→12",
        connectivity_patches=[
            (
                "preproject_bed_solid_holdup_to_transport_closure",
                _patch_solid_holdup_lsq(max_nfev=12),
            )
        ],
    ),
    ExperimentPatch(
        "lsq_skip_tight",
        "initial solid RMS 已小则跳过 LSQ（scaled<0.02）",
        connectivity_patches=[
            (
                "preproject_bed_solid_holdup_to_transport_closure",
                _patch_solid_holdup_lsq(max_nfev=32, skip_rms_scaled=0.02),
            )
        ],
    ),
    ExperimentPatch(
        "solid_single_pass",
        "bottom-zone solid 去掉第二轮 repass",
        connectivity_patches=[
            (
                "run_bottom_zone_solid_energy_preprojection",
                _patch_solid_energy_single_pass(),
            )
        ],
    ),
    ExperimentPatch(
        "preinner_no_continuation",
        "Check1 continuation outer 不再跑 pre-inner",
        config_overrides={"nr_outer_preinner_on_check1_continuation_thesis": False},
    ),
    ExperimentPatch(
        "jac_lag_8",
        "Jacobian lag 4→8",
        config_overrides={"nr_jacobian_lag_steps": 8},
    ),
    ExperimentPatch(
        "combo_fast",
        "组合：nfev=12 + skip + single_pass + no_continuation + lag=8",
        config_overrides={
            "nr_outer_preinner_on_check1_continuation_thesis": False,
            "nr_jacobian_lag_steps": 8,
        },
        connectivity_patches=[
            (
                "preproject_bed_solid_holdup_to_transport_closure",
                _patch_solid_holdup_lsq(max_nfev=12, skip_rms_scaled=0.02),
            ),
            (
                "run_bottom_zone_solid_energy_preprojection",
                _patch_solid_energy_single_pass(),
            ),
        ],
    ),
    ExperimentPatch(
        "nfev_12_single_pass",
        "nfev=12 + 去掉 solid repass",
        connectivity_patches=[
            (
                "preproject_bed_solid_holdup_to_transport_closure",
                _patch_solid_holdup_lsq(max_nfev=12),
            ),
            (
                "run_bottom_zone_solid_energy_preprojection",
                _patch_solid_energy_single_pass(),
            ),
        ],
    ),
    ExperimentPatch(
        "lsq_skip_loose",
        "skip scaled<0.01（比 0.02 保守）",
        connectivity_patches=[
            (
                "preproject_bed_solid_holdup_to_transport_closure",
                _patch_solid_holdup_lsq(max_nfev=32, skip_rms_scaled=0.01),
            )
        ],
    ),
]


def _apply_patches(patches: list[tuple[str, Callable[..., Any]]]) -> dict[str, Any]:
    import src.core.connectivity as conn

    saved: dict[str, Any] = {}
    for name, fn in patches:
        saved[name] = getattr(conn, name)
        setattr(conn, name, fn)
    return saved


def _restore_patches(saved: dict[str, Any]) -> None:
    import src.core.connectivity as conn

    for name, fn in saved.items():
        setattr(conn, name, fn)


def _build_reactor(overrides: dict[str, object]) -> Reactor:
    cfg = build_phase2_config()
    apply_config_overrides(cfg, LADDER_P2)
    apply_config_overrides(cfg, overrides)
    return Reactor(cfg)


def run_experiment(
    exp: ExperimentPatch,
    *,
    max_global_iter: int,
    skip_solve: bool,
) -> dict[str, object]:
    saved = _apply_patches(exp.connectivity_patches)
    try:
        reactor = _build_reactor(exp.config_overrides)
        sk = phase2_solve_kwargs()
        init_strategy = str(sk.get("nr_init_strategy", "vorabrechnung"))

        t0 = perf_counter()
        run_phase2_init(reactor, init_strategy=init_strategy)
        init_s = perf_counter() - t0

        cells = reactor._solver_cells_for_nr()
        apply = reactor._apply_all_bc_for_nr
        t0 = perf_counter()
        run_outer_preinner_start_value_bridge(
            cells,
            reactor.config,
            apply_bc_fn=apply,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        )
        preinner_once_s = perf_counter() - t0

        solve_s = 0.0
        result: dict[str, object] = {}
        if not skip_solve:
            t0 = perf_counter()
            result = reactor._solve_global_nr(
                max_iter=int(max_global_iter),
                tol=float(sk.get("tol_global", 1.0)),
                verbose=False,
                init_strategy=init_strategy,
                jacobian_strategy=sk.get("nr_jacobian_strategy"),
                skip_precalc=True,
            )
            solve_s = perf_counter() - t0

        timing = result.get("nr_timing") or result.get("timing") or {}
        counts = result.get("nr_counts") or result.get("counts") or {}
        outer_iters = int(result.get("nr_outer_iters", 0) or 0)

        return {
            "name": exp.name,
            "description": exp.description,
            "init_s": round(init_s, 3),
            "preinner_once_s": round(preinner_once_s, 3),
            "preinner_est_total_s": round(preinner_once_s * max(outer_iters, 1), 3),
            "solve_s": round(solve_s, 3),
            "grand_total_s": round(init_s + solve_s, 3),
            "converged": bool(result.get("converged", False)),
            "nr_outer_iters": outer_iters,
            "n_iter": int(result.get("n_iter", 0) or 0),
            "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
            "nr_inner_solve_s_total": round(float(result.get("nr_inner_solve_s_total", 0.0) or 0.0), 3),
            "jacobian_build_s": round(float(timing.get("jacobian_build_s", 0.0) or 0.0), 3),
            "line_search_residual_s": round(float(timing.get("line_search_residual_s", 0.0) or 0.0), 3),
            "jacobian_rebuilds": int(counts.get("jacobian_rebuilds", 0) or 0),
            "speedup_vs_baseline": None,
        }
    finally:
        _restore_patches(saved)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=12)
    parser.add_argument(
        "--only",
        nargs="*",
        default=None,
        help="Run subset of experiment names",
    )
    parser.add_argument(
        "--micro-only",
        action="store_true",
        help="Only measure init + one preinner call (skip full solve)",
    )
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "data" / "preinner_perf_experiments.json")
    args = parser.parse_args()

    selected = EXPERIMENTS
    if args.only:
        names = set(args.only)
        selected = [e for e in EXPERIMENTS if e.name in names]

    results: list[dict[str, object]] = []
    for exp in selected:
        print(f"Running {exp.name} ...", flush=True)
        results.append(
            run_experiment(
                exp,
                max_global_iter=int(args.max_global_iter),
                skip_solve=bool(args.micro_only),
            )
        )

    if not args.micro_only:
        baseline_total = next((r["grand_total_s"] for r in results if r["name"] == "baseline"), None)
        if baseline_total:
            for row in results:
                gt = float(row["grand_total_s"])
                row["speedup_vs_baseline"] = round(float(baseline_total) / gt, 3) if gt > 0 else None

    payload = {
        "max_global_iter": int(args.max_global_iter),
        "micro_only": bool(args.micro_only),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
