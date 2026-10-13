#!/usr/bin/env python3
"""local2+phase_share 栈：rms/energy 主导项 + 压残差窄 A/B（handoff §5.12）。

用法：
  python3 scripts/audit_local2_ps_abs_rms_stack.py --max-global-iter 20
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
    _residual_group_metrics,
    build_equation_scales,
    global_residual,
    pack_reactor,
)
from src.solvers.nr_indexing import gas_mode_for_nr, main_nr_gas_idx_for_cell, n_gas_var, n_var

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}

PS_STACK = {
    "nr_outer_preinner_upper_co_local2x2_thesis": True,
    "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
    "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
    "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 1.0,
    "nr_outer_preinner_upper_dense_absorb_mode_thesis": "phase_share",
    "nr_outer_preinner_upper_dense_absorb_post_gate_thesis": True,
}


def _topk_residuals(reactor: Any, k: int = 12) -> list[dict[str, Any]]:
    cfg = reactor.config
    cells = reactor.cells
    apply = reactor._apply_all_bc_for_nr
    apply()
    ref_gas = max(
        float(getattr(cfg, "O2_feed", 0.0) + getattr(cfg, "H2O_feed", 0.0) + getattr(cfg, "N2_feed", 0.0)),
        1.0,
    )
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    F = np.asarray(global_residual(pack_reactor(cells), cells, apply), dtype=np.float64)
    scale = np.asarray(build_equation_scales(cells, ref_gas, ref_solid, ref_energy), dtype=np.float64)
    F_hat = F / np.maximum(scale, 1e-30)
    groups = _residual_group_metrics(F_hat, cells)
    order = np.argsort(-np.abs(F_hat))[:k]
    rows: list[dict[str, Any]] = []
    for flat in order:
        offset = 0
        entry: dict[str, Any] = {
            "flat_k": int(flat),
            "F_hat": float(F_hat[flat]),
            "abs": float(abs(F_hat[flat])),
        }
        for i, cell in enumerate(cells):
            nv = int(n_var(cell))
            if offset <= int(flat) < offset + nv:
                local = int(flat) - offset
                ng = int(n_gas_var(cell))
                entry["cell"] = i
                entry["cell_type"] = str(getattr(cell, "cell_type", "?"))
                if local < ng and gas_mode_for_nr(cell) == "two_phase":
                    half = ng // 2
                    phase = "dense" if local < half else "bubble"
                    gi = local if local < half else local - half
                    gidx = main_nr_gas_idx_for_cell(cell)
                    sp = str(GAS_SPECIES[int(gidx[gi])]) if gi < len(gidx) else f"g{gi}"
                    entry["kind"] = "gas"
                    entry["phase"] = phase
                    entry["species"] = sp
                elif local == ng:
                    entry["kind"] = "energy"
                else:
                    entry["kind"] = "solid_or_other"
                    entry["local"] = local
                break
            offset += nv
        rows.append(entry)
    return rows, groups, float(np.sqrt(np.mean(F_hat**2))), float(np.max(np.abs(F_hat)))


def _history_summary(result: dict[str, Any]) -> dict[str, Any]:
    hist = list(result.get("nr_outer_history") or [])
    best = None
    traj = []
    for h in hist:
        row = {
            "outer": h.get("outer_iter"),
            "merit_best": h.get("outer_merit_best_so_far"),
            "split": h.get("gas_phase_split_rms"),
            "local2": h.get("upper_co_local2x2_accepted"),
            "absorb": h.get("upper_dense_absorb_accepted"),
            "addCO": h.get("upper_dense_absorb_added_CO_mol_s"),
        }
        traj.append(row)
        mb = h.get("outer_merit_best_so_far")
        if mb is not None and (best is None or float(mb) < float(best["merit_best"])):
            best = row
    return {"traj": traj, "best_so_far_row": best, "n_outer": len(hist)}


def _run(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    topk, groups, rms, max_abs = _topk_residuals(reactor)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    hist = _history_summary(result)
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "exit_H2": float(exit_gas.get("H2", float("nan"))),
        "exit_CO2": float(exit_gas.get("CO2", float("nan"))),
        "rms": float(result.get("rms_scaled_final", rms)),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", groups.get("energy_rms", float("nan")))),
        "max_abs": float(max_abs),
        "groups": {k: float(v) for k, v in groups.items()},
        "topk": topk,
        "history": hist,
        "overrides": dict(overrides),
        "n_acc": sum(
            1 for x in (result.get("nr_accepted_lambda_history") or []) if x is not None
        ),
        "ls_fail": int((result.get("nr_counts") or {}).get("line_search_failures", -1)),
        "converged": bool(result.get("converged", False)),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} split={out['split']:.4f} "
        f"E={out['energy']:.4f} max={out['max_abs']:.3f} "
        f"gas={groups.get('gas_rms', float('nan')):.4f} "
        f"sol={groups.get('solid_rms', float('nan')):.4f} "
        f"conv={out['converged']} n_acc={out['n_acc']} wall={wall_s:.0f}s"
    )
    for t in topk[:5]:
        print(
            f"  top |F|={t['abs']:.4f} cell{t.get('cell')} "
            f"{t.get('kind')} {t.get('phase','')}/{t.get('species','')} raw={t['F_hat']:+.4f}"
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=20)
    args = ap.parse_args()
    mi = int(args.max_global_iter)
    arms = [
        _run("baseline", {}, min(mi, 8)),
        _run("ps_stack_f1", dict(PS_STACK), mi),
        _run(
            "ps_stack_f0.5",
            {
                **PS_STACK,
                "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.5,
            },
            mi,
        ),
        _run(
            "ps_stack_f0.3",
            {
                **PS_STACK,
                "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.3,
            },
            mi,
        ),
        _run(
            "local2_only",
            {"nr_outer_preinner_upper_co_local2x2_thesis": True},
            mi,
        ),
    ]
    out = data_path(f"phase2_local2_ps_abs_rms_stack_o{mi}.json")
    write_json(
        out,
        {
            "note": "local2+phase_share rms/energy stack A/B; isomorphic",
            "max_global_iter": mi,
            "platform_y_CO": 0.154,
            "arms": arms,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
