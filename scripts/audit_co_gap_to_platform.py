#!/usr/bin/env python3
"""逐 bed CO 预算：local2 栈 vs baseline，相对平台对照 y_CO≈0.154 的缺口（handoff §5.11）。

对每个 bed 记录：
  R_d/R_b、N_d+N_b、res_sum、累计净产、出口 y 所需 N_top。

用法：
  python3 scripts/audit_co_gap_to_platform.py --max-global-iter 8
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

PLATFORM_Y_CO = 0.154


def _bed_co_rows(reactor: Any) -> list[dict[str, float | int]]:
    j = int(GAS_SPECIES_INDEX["CO"])
    rows: list[dict[str, float | int]] = []
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
        nex = float(cell.N_ex[j])
        res_sum = in_d + in_b + rd + rb - nd - nb
        res_split = (in_d + rd - nd + nex) - (in_b + rb - nb - nex)
        ntot = float(sum(cell.N_d) + sum(cell.N_b))
        rows.append(
            {
                "bed": i,
                "R_d": rd,
                "R_b": rb,
                "R_net": rd + rb,
                "N_d": nd,
                "N_b": nb,
                "N_tot": nd + nb,
                "in_tot": in_d + in_b,
                "N_ex": nex,
                "res_sum": res_sum,
                "res_split": res_split,
                "y_CO_cell": (nd + nb) / max(ntot, 1e-30),
                "N_gas_tot": ntot,
            }
        )
    return rows


def _summarize(label: str, result: dict[str, Any], reactor: Any, wall_s: float) -> dict[str, Any]:
    exit_gas = dict(result.get("bed_exit_gas") or {})
    rows = _bed_co_rows(reactor)
    top = rows[-1] if rows else {}
    n_top_co = float(top.get("N_tot", float("nan")))
    n_gas_top = float(top.get("N_gas_tot", float("nan")))
    y = float(exit_gas.get("CO", top.get("y_CO_cell", float("nan"))))
    n_need = PLATFORM_Y_CO * n_gas_top
    # remaining positive res_sum on upper beds = unabsorbed production
    upper = [r for r in rows if int(r["bed"]) >= 6]
    res_sum_pos = sum(max(float(r["res_sum"]), 0.0) for r in upper)
    res_sum_abs = sum(abs(float(r["res_sum"])) for r in upper)
    r_net_all = sum(float(r["R_net"]) for r in rows)
    r_net_upper = sum(float(r["R_net"]) for r in upper)
    r_net_lower = sum(float(r["R_net"]) for r in rows if int(r["bed"]) < 6)
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", y)),
        "exit_H2": float(exit_gas.get("H2", float("nan"))),
        "exit_CO2": float(exit_gas.get("CO2", float("nan"))),
        "rms": float(result.get("rms_scaled_final", float("nan"))),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "N_top_CO": n_top_co,
        "N_gas_top": n_gas_top,
        "N_top_CO_for_platform": float(n_need),
        "delta_N_top_vs_platform": float(n_top_co - n_need),
        "R_net_all": float(r_net_all),
        "R_net_lower0_5": float(r_net_lower),
        "R_net_upper6_9": float(r_net_upper),
        "upper_res_sum_pos": float(res_sum_pos),
        "upper_res_sum_abs": float(res_sum_abs),
        # 若把上段正 res_sum 全部吸入出口库存（粗估）
        "y_if_absorb_upper_res_sum": float((n_top_co + res_sum_pos) / max(n_gas_top, 1e-30)),
        "beds": rows,
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} N_top={n_top_co:.3f} need={n_need:.3f} "
        f"Δ={out['delta_N_top_vs_platform']:+.3f} "
        f"R_lo={r_net_lower:.2f} R_up={r_net_upper:.2f} "
        f"res+={res_sum_pos:.2f} y_if_abs={out['y_if_absorb_upper_res_sum']:.4f} "
        f"rms={out['rms']:.4f} split={out['split']:.4f}"
    )
    for r in rows:
        print(
            f"  bed{int(r['bed']):2d} R={float(r['R_net']):+7.3f} "
            f"N={float(r['N_tot']):7.3f} res_sum={float(r['res_sum']):+7.3f} "
            f"y={float(r['y_CO_cell']):.4f}"
        )
    return out


def _run(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    out = _summarize(label, result, reactor, wall_s)
    out["overrides"] = dict(overrides)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()
    arms = [
        _run("baseline", {}, args.max_global_iter),
        _run(
            "local2",
            {"nr_outer_preinner_upper_co_local2x2_thesis": True},
            args.max_global_iter,
        ),
    ]
    # delta table
    b0, b1 = arms[0], arms[1]
    delta_beds = []
    for r0, r1 in zip(b0["beds"], b1["beds"]):
        delta_beds.append(
            {
                "bed": r0["bed"],
                "dR": float(r1["R_net"]) - float(r0["R_net"]),
                "dN": float(r1["N_tot"]) - float(r0["N_tot"]),
                "d_res_sum": float(r1["res_sum"]) - float(r0["res_sum"]),
                "d_y": float(r1["y_CO_cell"]) - float(r0["y_CO_cell"]),
            }
        )
    print("local2 - baseline:")
    for d in delta_beds:
        print(
            f"  bed{int(d['bed']):2d} dR={d['dR']:+7.3f} dN={d['dN']:+7.3f} "
            f"d_res={d['d_res_sum']:+7.3f} d_y={d['d_y']:+.4f}"
        )

    out = data_path("phase2_co_gap_to_platform_o8.json")
    write_json(
        out,
        {
            "note": "per-bed CO budget vs platform y≈0.154; isomorphic o8",
            "platform_y_CO": PLATFORM_Y_CO,
            "max_global_iter": int(args.max_global_iter),
            "arms": arms,
            "local2_minus_baseline": delta_beds,
        },
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
