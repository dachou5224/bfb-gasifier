#!/usr/bin/env python3
"""LSQ n=2 栈上 bed3 N2 split 审计 + targeted split A/B（handoff §5.24）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed3_n2_on_lsq_n2_ab.py
"""

from __future__ import annotations

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
from audit_bed9_fb_solid_enthalpy import _Fhat_T
from audit_bed9_holdup_closure_ab import _metrics
from audit_rate_gate_rms_stack_ab import GATE_AFTER, HOT, ISOMORPH, PS

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity.preprojection import (
    _bed_major_gas_phase_split_rms,
    preproject_bed_major_gas_phase_split_to_exchange_closure,
)
from src.core.species import GAS_SPECIES_INDEX

STACK = {**GATE_AFTER, **HOT}

LSQ2 = {
    "nr_outer_preinner_bed_top_solid_holdup_thesis": True,
    "nr_outer_preinner_bed_top_solid_holdup_after_y_co_thesis": True,
    "nr_outer_preinner_bed_top_solid_holdup_y_co_threshold_thesis": 0.14,
    "nr_outer_preinner_bed_top_solid_holdup_n_beds_thesis": 2,
}

ENERGY_PT = {
    "nr_bed_top_solid_energy_passthrough_thesis": True,
    "nr_bed_top_solid_energy_passthrough_after_y_co_thesis": True,
    "nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis": 0.14,
}


def _bed3_n2_decomp(reactor: Any) -> dict[str, float]:
    cell = reactor.cells[3]
    reactor._apply_all_bc_for_nr()
    cell.residuals()
    j = int(GAS_SPECIES_INDEX["N2"])
    n_d = float(cell.N_d[j])
    n_b = float(cell.N_b[j])
    n_tot = n_d + n_b
    gas_d = float(np.sum(cell.N_d))
    gas_b = float(np.sum(cell.N_b))
    in_d = float(cell.N_zu_d[j] + cell.N_rez_d[j] + cell.N_d_in[j])
    in_b = float(cell.N_zu_b[j] + cell.N_rez_b[j] + cell.N_b_in[j])
    r_d = float(cell.R_gas_d[j])
    r_b = float(cell.R_gas_b[j])
    n_ex = float(cell.N_ex[j])
    res_d = in_d + r_d - n_d - n_ex
    res_b = in_b + r_b - n_b + n_ex
    return {
        "T": float(cell.T),
        "N_d": n_d,
        "N_b": n_b,
        "y_d": n_d / max(gas_d, 1e-30),
        "y_b": n_b / max(gas_b, 1e-30),
        "dy": n_d / max(gas_d, 1e-30) - n_b / max(gas_b, 1e-30),
        "N_ex": n_ex,
        "R_d": r_d,
        "R_b": r_b,
        "in_d": in_d,
        "in_b": in_b,
        "res_d": res_d,
        "res_b": res_b,
        "res_sum": res_d + res_b,
        "res_split": res_d - res_b,
        "split_rms_major": float(_bed_major_gas_phase_split_rms(cell)),
        "Fhat_T": float(_Fhat_T(reactor, 3)),
        "N_tot": n_tot,
    }


def _n2_fhat(reactor: Any) -> dict[str, float]:
    """bed3 dense/bubble N2 的缩放残差。"""
    from src.core.species import N_GAS
    from src.solvers.global_nr_solver import build_equation_scales, global_residual, pack_reactor
    from src.solvers.nr_indexing import n_gas_var, n_var

    cells = reactor.cells
    cfg = reactor.config
    apply = reactor._apply_all_bc_for_nr
    apply()
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    F = np.asarray(global_residual(pack_reactor(cells), cells, apply), dtype=np.float64)
    scale = np.asarray(build_equation_scales(cells, ref_gas, ref_solid, ref_energy), dtype=np.float64)
    F_hat = F / np.maximum(scale, 1e-30)
    off = 0
    j = int(GAS_SPECIES_INDEX["N2"])
    for i, cell in enumerate(cells):
        nv = int(n_var(cell))
        ng = int(n_gas_var(cell))
        if i == 3 and ng >= 2 * N_GAS:
            half = ng // 2
            return {
                "Fhat_N2_dense": float(F_hat[off + j]),
                "Fhat_N2_bubble": float(F_hat[off + half + j]),
            }
        off += nv
    return {"Fhat_N2_dense": float("nan"), "Fhat_N2_bubble": float("nan")}


def _posthoc_bed3(reactor: Any) -> dict[str, Any]:
    # shallow state restore for gas only
    cell = reactor.cells[3]
    Nd0 = cell.N_d.copy()
    Nb0 = cell.N_b.copy()
    before = {**_bed3_n2_decomp(reactor), **_n2_fhat(reactor)}
    ok = preproject_bed_major_gas_phase_split_to_exchange_closure(
        reactor.cells,
        reactor.config,
        cell_index=3,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        improvement_factor=0.5,
    )
    after = {**_bed3_n2_decomp(reactor), **_n2_fhat(reactor), "accepted": bool(ok)}
    cell.N_d[:] = Nd0
    cell.N_b[:] = Nb0
    cell._thermo_cache_valid = False
    reactor._apply_all_bc_for_nr()
    return {"before": before, "after_split_lsq": after}


def _run(label: str, overrides: dict[str, object], max_iter: int = 8) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(STACK)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    met = _metrics(reactor, result)
    bed3 = {**_bed3_n2_decomp(reactor), **_n2_fhat(reactor)}
    print(
        f"[{label}] CO={met['exit_CO']:.4f} rms={met['rms']:.4f} "
        f"split={met['groups'].get('gas_rms', float('nan')):.4f} "
        f"E={met['groups']['energy_rms']:.4f} "
        f"F_N2d={bed3['Fhat_N2_dense']:+.3f} dy={bed3['dy']:+.4f} "
        f"res_split={bed3['res_split']:+.3e} wall={wall_s:.0f}s",
        flush=True,
    )
    print(
        f"  top3: "
        + ", ".join(
            f"c{t['cell']}:{t.get('species') or t.get('kind')}={t['F_hat']:+.3f}"
            for t in met["topk"][:3]
        ),
        flush=True,
    )
    return {
        "label": label,
        "wall_s": float(wall_s),
        "metrics": met,
        "bed3_n2": bed3,
        "overrides": dict(overrides),
        "_reactor": reactor,
    }


def main() -> None:
    arms: list[dict[str, Any]] = []
    posthoc: dict[str, Any] = {}
    for label, ov in (
        ("lsq_n2_base", LSQ2),
        (
            "lsq_n2_maxbed3",
            {**LSQ2, "nr_outer_preinner_phase_split_max_bed_index_thesis": 3},
        ),
        (
            "lsq_n2_target_bed3_after_y_co",
            {
                **LSQ2,
                "nr_outer_preinner_targeted_phase_split_beds_thesis": (3,),
                "nr_outer_preinner_targeted_phase_split_after_y_co_thesis": True,
                "nr_outer_preinner_targeted_phase_split_y_co_threshold_thesis": 0.14,
            },
        ),
        (
            "lsq_n2_target_bed3_5_after_y_co",
            {
                **LSQ2,
                "nr_outer_preinner_targeted_phase_split_beds_thesis": (3, 5),
                "nr_outer_preinner_targeted_phase_split_after_y_co_thesis": True,
            },
        ),
        ("energy_pt_after_y_co", ENERGY_PT),
        ("lsq_n2_plus_energy_pt", {**LSQ2, **ENERGY_PT}),
    ):
        arm = _run(label, ov)
        reactor = arm.pop("_reactor")
        if label == "lsq_n2_base":
            posthoc = _posthoc_bed3(reactor)
            b = posthoc["before"]
            a = posthoc["after_split_lsq"]
            print(
                f"  [posthoc bed3 split] accepted={a.get('accepted')} "
                f"F_N2d {b['Fhat_N2_dense']:+.3f}→{a['Fhat_N2_dense']:+.3f} "
                f"dy {b['dy']:+.4f}→{a['dy']:+.4f} "
                f"res_split {b['res_split']:+.3e}→{a['res_split']:+.3e}",
                flush=True,
            )
        arms.append(arm)

    out = data_path("phase2_bed3_n2_on_lsq_n2_ab_o8.json")
    write_json(
        out,
        {
            "note": "bed3 N2 on LSQ n=2 stack + targeted split A/B; §5.24",
            "arms": arms,
            "posthoc_on_lsq_n2": posthoc,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
