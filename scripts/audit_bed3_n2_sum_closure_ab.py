#!/usr/bin/env python3
"""bed3 N2 res_sum 分解 + 双向总量闭合 A/B（handoff §5.25）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed3_n2_sum_closure_ab.py
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
from audit_bed3_n2_on_lsq_n2_ab import LSQ2, STACK, _bed3_n2_decomp
from audit_bed9_holdup_closure_ab import _metrics
from audit_rate_gate_rms_stack_ab import ISOMORPH, PS

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity.preprojection import (
    oneshot_bed_dense_species_sum_absorb,
    preproject_bed_major_gas_phase_split_to_exchange_closure,
    preproject_bed_temperature_to_local_energy_closure,
)
from src.core.species import GAS_SPECIES_INDEX

SUM3 = {
    "nr_outer_preinner_targeted_species_sum_beds_thesis": (3,),
    "nr_outer_preinner_targeted_species_sum_species_thesis": ("N2",),
    "nr_outer_preinner_targeted_species_sum_fraction_thesis": 1.0,
    "nr_outer_preinner_targeted_species_sum_mode_thesis": "phase_share",
    "nr_outer_preinner_targeted_species_sum_bidirectional_thesis": True,
    "nr_outer_preinner_targeted_species_sum_after_y_co_thesis": True,
    "nr_outer_preinner_targeted_species_sum_y_co_threshold_thesis": 0.14,
}

SPLIT3 = {
    "nr_outer_preinner_targeted_phase_split_beds_thesis": (3,),
    "nr_outer_preinner_targeted_phase_split_after_y_co_thesis": True,
}


def _budget_row(reactor: Any) -> dict[str, float]:
    """N2 总量预算：入口 vs holdup（N_ex 在 sum 中对消）。"""
    d = _bed3_n2_decomp(reactor)
    n_hold = float(d["N_d"] + d["N_b"])
    n_in = float(d["in_d"] + d["in_b"])
    return {
        **d,
        "N_hold": n_hold,
        "N_in": n_in,
        "gap_hold_minus_in": n_hold - n_in,  # ≈ -res_sum when R=0
        "N_ex_abs": abs(float(d["N_ex"])),
    }


def _posthoc(reactor: Any) -> dict[str, Any]:
    cell = reactor.cells[3]
    Nd0, Nb0 = cell.N_d.copy(), cell.N_b.copy()
    T0 = float(cell.T)

    def _restore() -> None:
        cell.N_d[:] = Nd0
        cell.N_b[:] = Nb0
        cell.T = T0
        cell._thermo_cache_valid = False
        reactor._apply_all_bc_for_nr()

    out: dict[str, Any] = {"before": _budget_row(reactor)}

    # 单向 absorb（应无效，因 res_sum<0）
    _restore()
    cell.residuals()
    added = oneshot_bed_dense_species_sum_absorb(
        cell, species=("N2",), mode="phase_share", bidirectional=False
    )
    out["absorb_one_way"] = {**_budget_row(reactor), "added_N2": float(added.get("N2", 0.0))}

    # 双向放出
    _restore()
    cell.residuals()
    added = oneshot_bed_dense_species_sum_absorb(
        cell, species=("N2",), mode="phase_share", bidirectional=True, fraction=1.0
    )
    out["bleed_bidirectional"] = {
        **_budget_row(reactor),
        "added_N2": float(added.get("N2", 0.0)),
    }

    # 双向 + split
    _restore()
    cell.residuals()
    oneshot_bed_dense_species_sum_absorb(
        cell, species=("N2",), mode="phase_share", bidirectional=True, fraction=1.0
    )
    ok = preproject_bed_major_gas_phase_split_to_exchange_closure(
        reactor.cells,
        reactor.config,
        cell_index=3,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        improvement_factor=0.5,
    )
    out["bleed_then_split"] = {**_budget_row(reactor), "split_accepted": bool(ok)}

    # bed3 T 能量投影
    _restore()
    ok = preproject_bed_temperature_to_local_energy_closure(
        reactor.cells,
        reactor.config,
        cell_index=3,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
    )
    out["T_energy_bed3"] = {**_budget_row(reactor), "accepted": bool(ok)}

    _restore()
    return out


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
    bud = _budget_row(reactor)
    print(
        f"[{label}] CO={met['exit_CO']:.4f} rms={met['rms']:.4f} "
        f"E={met['groups']['energy_rms']:.4f} gas={met['groups']['gas_rms']:.4f} "
        f"res_sum={bud['res_sum']:+.3e} gap={bud['gap_hold_minus_in']:+.3e} "
        f"dy={bud['dy']:+.4f} wall={wall_s:.0f}s",
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
        "budget": bud,
        "overrides": dict(overrides),
        "_reactor": reactor,
    }


def main() -> None:
    arms: list[dict[str, Any]] = []
    posthoc: dict[str, Any] = {}
    for label, ov in (
        ("lsq_n2_base", LSQ2),
        ("lsq_n2_n2sum_bed3", {**LSQ2, **SUM3}),
        ("lsq_n2_n2sum_then_split", {**LSQ2, **SUM3, **SPLIT3}),
        (
            "lsq_n2_n2sum_f0p5",
            {
                **LSQ2,
                **SUM3,
                "nr_outer_preinner_targeted_species_sum_fraction_thesis": 0.5,
            },
        ),
    ):
        arm = _run(label, ov)
        reactor = arm.pop("_reactor")
        if label == "lsq_n2_base":
            posthoc = _posthoc(reactor)
            for name, ph in posthoc.items():
                if name == "before":
                    continue
                print(
                    f"  [posthoc {name}] res_sum={ph.get('res_sum', float('nan')):+.3e} "
                    f"gap={ph.get('gap_hold_minus_in', float('nan')):+.3e} "
                    f"dy={ph.get('dy', float('nan')):+.4f} "
                    f"added={ph.get('added_N2', ph.get('accepted', ''))}",
                    flush=True,
                )
        arms.append(arm)

    out = data_path("phase2_bed3_n2_sum_closure_ab_o8.json")
    write_json(
        out,
        {
            "note": "bed3 N2 res_sum budget + bidirectional sum closure; §5.25",
            "arms": arms,
            "posthoc_on_lsq_n2": posthoc,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
