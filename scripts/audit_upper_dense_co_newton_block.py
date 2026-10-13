#!/usr/bin/env python3
"""平台态：上段 dense-CO 的 Newton dx / clip / fastox 阻塞诊断（handoff §5.6）。

同构 o8 末态后取一步线性解，报告 bed6–9 的 CO DOF：
1. F_hat（dense/bubble）与物理 res_sum
2. 原始 dx、clip 后 dx、全局 α
3. fastox（local_participating）对 N_d,CO / N_b,CO 的改写
4. J 列范数与对角近似 ∂F_d/∂N_d

用法：
  python3 scripts/audit_upper_dense_co_newton_block.py
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

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

from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX
from src.solvers.global_nr_solver import (
    _clip_dx,
    _clip_dx_diagnostics,
    _line_search_merit,
    _residual_gas_phase_metrics,
    _residual_group_metrics,
    _rms_norm,
    _solve_linear_step,
    build_equation_scales,
    build_jacobian_fd,
    global_residual,
    pack_reactor,
    unpack_reactor,
)
from src.solvers.nr_indexing import (
    gas_mode_for_nr,
    install_main_nr_gas_indexing_on_cells,
    main_nr_gas_idx_for_cell,
    n_gas_var,
    n_var,
)
from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure
from types import SimpleNamespace

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}


def _dof_map(cells: list[Any]) -> list[dict[str, Any]]:
    """Map (cell, phase, species) -> global DOF index."""
    rows: list[dict[str, Any]] = []
    offset = 0
    for i, cell in enumerate(cells):
        nv = int(n_var(cell))
        ng = int(n_gas_var(cell))
        if gas_mode_for_nr(cell) == "two_phase" and ng > 0:
            gas_idx = main_nr_gas_idx_for_cell(cell)
            half = ng // 2
            for p, gi in enumerate(gas_idx):
                sp = str(GAS_SPECIES[int(gi)])
                rows.append(
                    {
                        "cell": i,
                        "phase": "dense",
                        "sp": sp,
                        "dof": offset + p,
                        "res_dof": offset + p,
                    }
                )
                rows.append(
                    {
                        "cell": i,
                        "phase": "bubble",
                        "sp": sp,
                        "dof": offset + half + p,
                        "res_dof": offset + half + p,
                    }
                )
        offset += nv
    return rows


def _phys_co_res(cell: Any) -> dict[str, float]:
    j = int(GAS_SPECIES_INDEX["CO"])
    in_d = float(cell.N_zu_d[j] + cell.N_rez_d[j] + cell.N_d_in[j])
    in_b = float(cell.N_zu_b[j] + cell.N_rez_b[j] + cell.N_b_in[j])
    R_d = float(cell.R_gas_d[j])
    R_b = float(cell.R_gas_b[j])
    N_d = float(cell.N_d[j])
    N_b = float(cell.N_b[j])
    N_ex = float(cell.N_ex[j])
    res_d = in_d + R_d - N_d + N_ex
    res_b = in_b + R_b - N_b - N_ex
    return {
        "N_d": N_d,
        "N_b": N_b,
        "R_d": R_d,
        "R_b": R_b,
        "N_ex": N_ex,
        "res_d": res_d,
        "res_b": res_b,
        "res_sum": res_d + res_b,
        "res_split": res_d - res_b,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()

    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=ISOMORPH)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter))
    )
    assert isinstance(result, dict)

    cells = list(reactor._solver_cells_for_nr())
    install_main_nr_gas_indexing_on_cells(cells, reactor.config)
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)

    x0 = pack_reactor(cells)
    F0 = global_residual(x0, cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F0, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    metrics0 = {
        "rms": float(_rms_norm(F_hat)),
        "split": float(_residual_gas_phase_metrics(F_hat, cells)["gas_phase_split_rms"]),
        "energy": float(_residual_group_metrics(F_hat, cells)["energy_rms"]),
        "merit": float(
            _line_search_merit(F_hat, cells, use_gas_phase_split=True, use_energy=True)
        ),
    }

    J, jac_meta = build_jacobian_fd(
        x0,
        F0,
        cells,
        reactor._apply_all_bc_for_nr,
        scale,
        strategy="hybrid_band_plus_side_elements_structured",
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        affected_residual_cells_fn=reactor._affected_nr_residual_cells,
        connectivity_graph=reactor.connectivity_graph,
    )
    dx, method, _reg, _diag = _solve_linear_step(
        J,
        -F_hat,
        structured_jacobian=None,
        linear_solver_backend="sparse_direct_fallback",
        zero_rows=int(jac_meta.get("zero_rows", 0)),
        zero_cols=int(jac_meta.get("zero_cols", 0)),
    )
    dx = np.asarray(dx, dtype=np.float64)
    dx_clip = _clip_dx(dx, cells, ref_gas, ref_solid, t_step_limit_K=150.0)
    clip_diag = _clip_dx_diagnostics(dx, dx_clip, cells)

    # alpha estimate from max |dx_gas|/limit
    rg_limit = 0.2 * max(float(ref_gas), 1.0)
    gas_dof_mask = np.zeros(dx.size, dtype=bool)
    off = 0
    for cell in cells:
        nv = int(n_var(cell))
        ng = int(n_gas_var(cell))
        if ng > 0:
            gas_dof_mask[off : off + ng] = True
        off += nv
    gas_dx_max = float(np.max(np.abs(dx[gas_dof_mask]))) if np.any(gas_dof_mask) else 0.0
    alpha_gas = float(min(1.0, rg_limit / gas_dx_max)) if gas_dx_max > 0 else 1.0

    dof_rows = _dof_map(cells)
    Jcsc = J.tocsc()

    upper_rows = []
    for i in range(6, 10):
        cell = reactor.cells[i]
        cell.residuals()
        phys = _phys_co_res(cell)
        dense = next(r for r in dof_rows if r["cell"] == i and r["phase"] == "dense" and r["sp"] == "CO")
        bubble = next(r for r in dof_rows if r["cell"] == i and r["phase"] == "bubble" and r["sp"] == "CO")
        col_d = int(dense["dof"])
        col_b = int(bubble["dof"])
        row_d = int(dense["res_dof"])
        row_b = int(bubble["res_dof"])
        Jd = Jcsc.getcol(col_d)
        Jb = Jcsc.getcol(col_b)
        upper_rows.append(
            {
                "cell": i,
                **phys,
                "Fhat_d": float(F_hat[row_d]),
                "Fhat_b": float(F_hat[row_b]),
                "scale_d": float(scale[row_d]),
                "scale_b": float(scale[row_b]),
                "dx_raw_d": float(dx[col_d]),
                "dx_raw_b": float(dx[col_b]),
                "dx_clip_d": float(dx_clip[col_d]),
                "dx_clip_b": float(dx_clip[col_b]),
                "Jcol_norm_d": float(np.linalg.norm(Jd.toarray())),
                "Jcol_norm_b": float(np.linalg.norm(Jb.toarray())),
                "J_dd": float(Jcsc[row_d, col_d]),
                "J_db": float(Jcsc[row_d, col_b]),
                "J_bd": float(Jcsc[row_b, col_d]),
                "J_bb": float(Jcsc[row_b, col_b]),
                # ideal absorb step for res_sum (dense only)
                "ideal_absorb_d": float(max(phys["res_sum"], 0.0)),
            }
        )

    # Apply clip step + fastox; measure CO change on upper beds
    unpack_reactor(x0 + dx_clip, cells)
    reactor._apply_all_bc_for_nr()
    co_after_clip = []
    for i in range(6, 10):
        j = GAS_SPECIES_INDEX["CO"]
        co_after_clip.append(
            {
                "cell": i,
                "N_d": float(reactor.cells[i].N_d[j]),
                "N_b": float(reactor.cells[i].N_b[j]),
            }
        )
    project_bed_holdup_fast_oxidation_closure(
        [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"],
        SimpleNamespace(
            vorab_transport_x0_fast_oxidation_closure_thesis=True,
            fuel_type=str(getattr(reactor.cells[0], "fuel_type", "coal")),
        ),
        phase_update="local_participating",
    )
    co_after_fastox = []
    for i in range(6, 10):
        j = GAS_SPECIES_INDEX["CO"]
        co_after_fastox.append(
            {
                "cell": i,
                "N_d": float(reactor.cells[i].N_d[j]),
                "N_b": float(reactor.cells[i].N_b[j]),
            }
        )

    # Restore and counterfactual: force dense CO absorb = res_sum on dx, then clip+fastox
    unpack_reactor(x0, cells)
    reactor._apply_all_bc_for_nr()
    dx_forced = dx_clip.copy()
    for row in upper_rows:
        dense = next(
            r
            for r in dof_rows
            if r["cell"] == row["cell"] and r["phase"] == "dense" and r["sp"] == "CO"
        )
        dx_forced[int(dense["dof"])] = float(row["ideal_absorb_d"])
    unpack_reactor(x0 + dx_forced, cells)
    reactor._apply_all_bc_for_nr()
    F_force = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    Fh_force = np.asarray(F_force, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    metrics_force = {
        "rms": float(_rms_norm(Fh_force)),
        "split": float(_residual_gas_phase_metrics(Fh_force, cells)["gas_phase_split_rms"]),
        "energy": float(_residual_group_metrics(Fh_force, cells)["energy_rms"]),
        "merit": float(
            _line_search_merit(Fh_force, cells, use_gas_phase_split=True, use_energy=True)
        ),
    }

    # Summary ratios
    for row in upper_rows:
        row["dx_clip_over_ideal"] = (
            float(row["dx_clip_d"] / row["ideal_absorb_d"])
            if row["ideal_absorb_d"] > 1e-12
            else float("nan")
        )
        row["dx_raw_over_ideal"] = (
            float(row["dx_raw_d"] / row["ideal_absorb_d"])
            if row["ideal_absorb_d"] > 1e-12
            else float("nan")
        )

    payload = {
        "note": "upper dense-CO Newton/clip/fastox block at isomorphic plateau",
        "wall_s": float(wall_s),
        "exit_CO": float((result.get("bed_exit_gas") or {}).get("CO", float("nan"))),
        "metrics0": metrics0,
        "linear_method": str(method),
        "jac_zero_cols": int(jac_meta.get("zero_cols", 0)),
        "jac_zero_rows": int(jac_meta.get("zero_rows", 0)),
        "rg_limit": float(rg_limit),
        "gas_dx_max_raw": float(gas_dx_max),
        "alpha_gas_from_rg": float(alpha_gas),
        "clip_diag": {
            "clipped_fraction": float(clip_diag.get("clipped_fraction", float("nan"))),
            "clipped_components": int(clip_diag.get("clipped_components", -1)),
            "max_abs_dT_step": float(clip_diag.get("max_abs_dT_step", float("nan"))),
        },
        "upper_CO": upper_rows,
        "co_after_newton_clip": co_after_clip,
        "co_after_fastox": co_after_fastox,
        "metrics_force_full_absorb_on_clip_step": metrics_force,
        "mean_dx_clip_over_ideal": float(
            np.nanmean([r["dx_clip_over_ideal"] for r in upper_rows])
        ),
        "mean_dx_raw_over_ideal": float(
            np.nanmean([r["dx_raw_over_ideal"] for r in upper_rows])
        ),
    }
    out = data_path("phase2_upper_dense_co_newton_block_o8.json")
    write_json(out, payload)

    print(f"wrote {out}")
    print(f"metrics0={metrics0}")
    print(
        f"gas_dx_max={gas_dx_max:.3e} rg_limit={rg_limit:.3e} alpha_gas≈{alpha_gas:.4f} "
        f"method={method} zero_cols={payload['jac_zero_cols']}"
    )
    print(
        f"mean dx_raw/ideal={payload['mean_dx_raw_over_ideal']:.4f} "
        f"dx_clip/ideal={payload['mean_dx_clip_over_ideal']:.4f}"
    )
    print(
        f"{'cell':>4} {'Nd':>8} {'res_sum':>8} {'dx_raw':>9} {'dx_clip':>9} "
        f"{'ideal':>8} {'raw/id':>7} {'Jdd':>9} {'Fhat_d':>9}"
    )
    for r in upper_rows:
        print(
            f"{r['cell']:4d} {r['N_d']:8.3f} {r['res_sum']:8.3f} {r['dx_raw_d']:9.3f} "
            f"{r['dx_clip_d']:9.3f} {r['ideal_absorb_d']:8.3f} "
            f"{r['dx_raw_over_ideal']:7.3f} {r['J_dd']:9.3e} {r['Fhat_d']:9.3f}"
        )
    print("CO after clip vs fastox:")
    for a, b in zip(co_after_clip, co_after_fastox):
        print(
            f"  bed{a['cell']}: clip Nd={a['N_d']:.4f} Nb={a['N_b']:.4f} -> "
            f"fastox Nd={b['N_d']:.4f} Nb={b['N_b']:.4f} "
            f"ΔNd={b['N_d']-a['N_d']:+.4f}"
        )
    print("force full absorb on clipped-step DOFs ->", metrics_force)


if __name__ == "__main__":
    main()
