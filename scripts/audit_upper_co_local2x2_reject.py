#!/usr/bin/env python3
"""上段 CO local2x2：outer≥3 门禁拒因 + CO 预算缺口（handoff §5.8）。

区分：
1. local2 本地门禁（split/merit）拒
2. local2 本地过、全局 Pre-Inner 门禁回滚
3. 末态 CO 产/耗/库存残差 vs 出口 mole fraction

用法：
  python3 scripts/audit_upper_co_local2x2_reject.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
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

from src.core import connectivity as conn
from src.core.connectivity import preprojection as pp
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


def _co_budget(reactor: Any) -> dict[str, Any]:
    j = int(GAS_SPECIES_INDEX["CO"])
    rows: list[dict[str, float | int]] = []
    R_net = 0.0
    N_out = 0.0
    res_sum_abs = 0.0
    res_split_abs = 0.0
    n_upper = 0
    for i, cell in enumerate(reactor.cells):
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        cell.residuals()
        rd = float(cell.R_gas_d[j])
        rb = float(cell.R_gas_b[j])
        nd = float(cell.N_d[j])
        nb = float(cell.N_b[j])
        in_d = float(cell.N_zu_d[j] + cell.N_rez_d[j] + cell.N_d_in[j])
        in_b = float(cell.N_zu_b[j] + cell.N_rez_b[j] + cell.N_b_in[j])
        res_sum = in_d + in_b + rd + rb - nd - nb
        res_split = (in_d + rd - nd + float(cell.N_ex[j])) - (
            in_b + rb - nb - float(cell.N_ex[j])
        )
        R_net += rd + rb
        if i == len([c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]) - 1:
            N_out = nd + nb
        row = {
            "bed": i,
            "R_d": rd,
            "R_b": rb,
            "N_d": nd,
            "N_b": nb,
            "N_ex": float(cell.N_ex[j]),
            "res_sum": res_sum,
            "res_split": res_split,
            "T": float(getattr(cell, "T", getattr(cell, "T_g", float("nan")))),
        }
        rows.append(row)
        if i >= 6:
            res_sum_abs += abs(res_sum)
            res_split_abs += abs(res_split)
            n_upper += 1
    # bed-top exit from packed result path often uses top bed N; also y from totals
    top = max((r for r in rows if True), key=lambda r: int(r["bed"]))
    N_top = float(top["N_d"]) + float(top["N_b"])
    # total gas at top for mole fraction proxy
    top_cell = reactor.cells[int(top["bed"])]
    Ntot = float(sum(top_cell.N_d) + sum(top_cell.N_b))
    y_CO = N_top / max(Ntot, 1e-30)
    return {
        "R_net_CO_mol_s": float(R_net),
        "N_top_CO_mol_s": float(N_top),
        "y_CO_top": float(y_CO),
        "upper_mean_abs_res_sum": float(res_sum_abs / max(n_upper, 1)),
        "upper_mean_abs_res_split": float(res_split_abs / max(n_upper, 1)),
        "beds": rows,
    }


def _run(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(overrides)

    local_traces: list[dict[str, Any]] = []
    bridge_traces: list[dict[str, Any]] = []
    _orig_local2 = pp.run_outer_preinner_upper_co_local2x2
    _orig_bridge = pp.run_outer_preinner_start_value_bridge
    call_i = {"n": 0}

    def _wrap_local2(*args: Any, **kwargs: Any) -> dict[str, Any]:
        call_i["n"] += 1
        diag = dict(_orig_local2(*args, **kwargs))
        local_traces.append(
            {
                "call": int(call_i["n"]),
                "local_accepted": bool(diag.get("accepted", False)),
                "merit_before": float(diag.get("merit_before", float("nan"))),
                "merit_after": float(diag.get("merit_after", float("nan"))),
                "split_before": float(diag.get("split_rms_before", float("nan"))),
                "split_after": float(diag.get("split_rms_after", float("nan"))),
                "beds": list(diag.get("beds_touched", []) or []),
                "dx2_by_bed": {
                    k: [float(x) for x in v]
                    for k, v in dict(diag.get("dx2_by_bed") or {}).items()
                },
            }
        )
        return diag

    def _wrap_bridge(*args: Any, **kwargs: Any) -> dict[str, Any]:
        diag = dict(_orig_bridge(*args, **kwargs))
        # Infer reject class using last local2 trace if enabled this call
        loc = local_traces[-1] if local_traces else None
        local_ok = bool(loc and loc.get("local_accepted"))
        global_ok = bool(diag.get("global_gate_pass", True))
        reported = bool(diag.get("upper_co_local2x2_accepted", False))
        if not bool(diag.get("upper_co_local2x2_enabled", False)):
            reject_class = "disabled"
        elif reported:
            reject_class = "accepted"
        elif not local_ok:
            reject_class = "local_gate"
        elif local_ok and not global_ok:
            reject_class = "global_rollback"
        else:
            reject_class = "unknown"
        bridge_traces.append(
            {
                "call": len(bridge_traces) + 1,
                "global_gate_pass": global_ok,
                "split_only_fallback": bool(diag.get("split_only_fallback", False)),
                "merit_before": float(diag.get("merit_before", float("nan"))),
                "merit_after": float(diag.get("merit_after", float("nan"))),
                "split_before": float(diag.get("split_rms_before", float("nan"))),
                "split_after": float(diag.get("split_rms_after", float("nan"))),
                "solid_before": float(diag.get("solid_rms_before", float("nan"))),
                "solid_after": float(diag.get("solid_rms_after", float("nan"))),
                "energy_before": float(diag.get("energy_rms_before", float("nan"))),
                "energy_after": float(diag.get("energy_rms_after", float("nan"))),
                "local2_reported_accepted": reported,
                "local2_local_accepted": local_ok,
                "reject_class": reject_class,
                "d_merit": float(diag.get("merit_after", 0.0))
                - float(diag.get("merit_before", 0.0)),
                "d_split": float(diag.get("split_rms_after", 0.0))
                - float(diag.get("split_rms_before", 0.0)),
                "d_energy": float(diag.get("energy_rms_after", 0.0))
                - float(diag.get("energy_rms_before", 0.0)),
                "d_solid": float(diag.get("solid_rms_after", 0.0))
                - float(diag.get("solid_rms_before", 0.0)),
            }
        )
        return diag

    _orig_bridge_pkg = conn.run_outer_preinner_start_value_bridge
    pp.run_outer_preinner_upper_co_local2x2 = _wrap_local2  # type: ignore[assignment]
    # bridge 在 preprojection 内用局部名调 local2；reactor 从 connectivity 包导入 bridge
    pp.run_outer_preinner_start_value_bridge = _wrap_bridge  # type: ignore[assignment]
    conn.run_outer_preinner_start_value_bridge = _wrap_bridge  # type: ignore[assignment]
    try:
        reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
        run_phase2_init(reactor)
        result, wall_s = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
        )
    finally:
        pp.run_outer_preinner_upper_co_local2x2 = _orig_local2  # type: ignore[assignment]
        pp.run_outer_preinner_start_value_bridge = _orig_bridge  # type: ignore[assignment]
        conn.run_outer_preinner_start_value_bridge = _orig_bridge_pkg  # type: ignore[assignment]

    assert isinstance(result, dict)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    budget = _co_budget(reactor)
    # classify counts
    classes: dict[str, int] = {}
    for b in bridge_traces:
        c = str(b["reject_class"])
        classes[c] = classes.get(c, 0) + 1

    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "exit_H2": float(exit_gas.get("H2", float("nan"))),
        "exit_CO2": float(exit_gas.get("CO2", float("nan"))),
        "exit_CH4": float(exit_gas.get("CH4", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "reject_class_counts": classes,
        "local_traces": local_traces,
        "bridge_traces": bridge_traces,
        "co_budget": {
            k: v for k, v in budget.items() if k != "beds"
        },
        "co_budget_beds": budget["beds"],
        "overrides": dict(overrides),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} y_top={budget['y_CO_top']:.4f} "
        f"R_net={budget['R_net_CO_mol_s']:.3f} "
        f"|res_sum|={budget['upper_mean_abs_res_sum']:.3f} "
        f"|res_split|={budget['upper_mean_abs_res_split']:.3f} "
        f"classes={classes} wall={wall_s:.0f}s"
    )
    for b in bridge_traces:
        if b["reject_class"] != "disabled":
            print(
                f"  bridge#{b['call']}: {b['reject_class']} "
                f"dM={b['d_merit']:+.4f} dS={b['d_split']:+.4f} "
                f"dE={b['d_energy']:+.4f} dSol={b['d_solid']:+.4f} "
                f"local={b['local2_local_accepted']} global={b['global_gate_pass']}"
            )
    for loc in local_traces:
        if not loc["local_accepted"]:
            print(
                f"  local#{loc['call']}: REJECT merit "
                f"{loc['merit_before']:.4f}→{loc['merit_after']:.4f} "
                f"split {loc['split_before']:.4f}→{loc['split_after']:.4f}"
            )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()
    arms = [
        ("baseline", {}),
        (
            "local2x2",
            {"nr_outer_preinner_upper_co_local2x2_thesis": True},
        ),
        (
            "local2x2_no_global_gate",
            {
                "nr_outer_preinner_upper_co_local2x2_thesis": True,
                "nr_outer_preinner_global_merit_gate_thesis": False,
            },
        ),
    ]
    rows = [_run(label, ov, args.max_global_iter) for label, ov in arms]
    out = data_path("phase2_upper_co_local2x2_reject_o8.json")
    write_json(
        out,
        {
            "note": "local2x2 reject class + CO budget; isomorphic o8",
            "max_global_iter": int(args.max_global_iter),
            "platform_CO_ref": 0.154,
            "arms": rows,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
