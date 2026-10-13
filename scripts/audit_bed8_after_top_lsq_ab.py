#!/usr/bin/env python3
"""bed-top LSQ 后 bed8 能量/holdup 审计 + 双床顶 LSQ / split A/B（handoff §5.23）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed8_after_top_lsq_ab.py
"""

from __future__ import annotations

from typing import Any

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from audit_bed9_fb_solid_enthalpy import _enthalpy_gap, _Fhat_T
from audit_bed9_holdup_closure_ab import _bed_solid_row, _chain_table, _metrics
from audit_rate_gate_rms_stack_ab import GATE_AFTER, HOT, ISOMORPH, PS

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity.preprojection import (
    preproject_bed_solid_holdup_to_transport_closure,
    preproject_bed_temperature_to_local_energy_closure,
)

STACK = {**GATE_AFTER, **HOT}

LSQ1 = {
    "nr_outer_preinner_bed_top_solid_holdup_thesis": True,
    "nr_outer_preinner_bed_top_solid_holdup_after_y_co_thesis": True,
    "nr_outer_preinner_bed_top_solid_holdup_y_co_threshold_thesis": 0.14,
    "nr_outer_preinner_bed_top_solid_holdup_n_beds_thesis": 1,
}

LSQ2 = {
    **LSQ1,
    "nr_outer_preinner_bed_top_solid_holdup_n_beds_thesis": 2,
}

ENERGY_PT = {
    "nr_bed_top_solid_energy_passthrough_thesis": True,
    "nr_bed_top_solid_energy_passthrough_after_y_co_thesis": True,
    "nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis": 0.14,
}


def _upper_energy(reactor: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for i in (7, 8, 9):
        cell = reactor.cells[i]
        out[f"bed{i}"] = {
            **_bed_solid_row(cell, i),
            "enthalpy": _enthalpy_gap(cell),
            "Fhat_T": _Fhat_T(reactor, i),
            "T": float(cell.T),
        }
    return out


def _posthoc_on_lsq1(reactor: Any) -> dict[str, Any]:
    """在 LSQ1 末态上试 bed8 同态修复（不重跑 NR）。"""
    import numpy as np

    snap = [np.asarray(c.m_solid, dtype=np.float64).copy() for c in reactor.cells]
    T_snap = [float(c.T) for c in reactor.cells]

    def _restore() -> None:
        for c, m, T in zip(reactor.cells, snap, T_snap):
            c.m_solid[:] = m
            c.T = T
            c._thermo_cache_valid = False
            c._hydro_cache_valid = False
        reactor._apply_all_bc_for_nr()

    out: dict[str, Any] = {"before": _upper_energy(reactor)}

    _restore()
    ok = preproject_bed_solid_holdup_to_transport_closure(
        reactor.cells,
        reactor.config,
        cell_index=8,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        skip_rms_scaled=0.0,
        improvement_factor=0.5,
    )
    out["lsq_bed8"] = {**_upper_energy(reactor), "accepted": bool(ok)}

    _restore()
    ok = preproject_bed_temperature_to_local_energy_closure(
        reactor.cells,
        reactor.config,
        cell_index=8,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
    )
    out["T_energy_bed8"] = {**_upper_energy(reactor), "accepted": bool(ok)}

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
    upper = _upper_energy(reactor)
    print(
        f"[{label}] CO={met['exit_CO']:.4f} rms={met['rms']:.4f} "
        f"E={met['groups']['energy_rms']:.4f} "
        f"F8={upper['bed8']['Fhat_T']:+.3f} F9={upper['bed9']['Fhat_T']:+.3f} "
        f"res8={upper['bed8']['res_char']:+.3e} res9={upper['bed9']['res_char']:+.3e} "
        f"wall={wall_s:.0f}s",
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
        "upper": upper,
        "chain_tail": _chain_table(reactor)[7:],
        "overrides": dict(overrides),
        "_reactor": reactor,
    }


def main() -> None:
    arms: list[dict[str, Any]] = []
    posthoc: dict[str, Any] = {}

    for label, ov in (
        ("gate_hot_base", {}),
        ("lsq_n1_after_y_co", LSQ1),
        ("lsq_n2_after_y_co", LSQ2),
        ("lsq_n1_plus_energy_pt", {**LSQ1, **ENERGY_PT}),
        (
            "lsq_n1_plus_split_max3",
            {
                **LSQ1,
                "nr_outer_preinner_phase_split_max_bed_index_thesis": 3,
            },
        ),
    ):
        arm = _run(label, ov)
        reactor = arm.pop("_reactor")
        if label == "lsq_n1_after_y_co":
            posthoc = _posthoc_on_lsq1(reactor)
            for name, ph in posthoc.items():
                if name == "before":
                    continue
                b8 = ph.get("bed8") or {}
                print(
                    f"  [posthoc {name}] F8={b8.get('Fhat_T', float('nan')):+.3f} "
                    f"res8={b8.get('res_char', float('nan')):+.3e} "
                    f"resE8={((b8.get('enthalpy') or {}).get('res_Eout', float('nan'))):+.3e} "
                    f"acc={ph.get('accepted')}",
                    flush=True,
                )
        arms.append(arm)

    out = data_path("phase2_bed8_after_top_lsq_ab_o8.json")
    write_json(
        out,
        {
            "note": "bed8 after top-LSQ + n_beds/split A/B; §5.23",
            "arms": arms,
            "posthoc_on_lsq_n1": posthoc,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
