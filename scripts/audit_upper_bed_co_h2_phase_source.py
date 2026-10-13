#!/usr/bin/env python3
"""上段床 CO/H2（+H2O/N2）相源定位与 split 项分解（handoff §5.4）。

同构栈 o8 末态，对 bed6–9 报告：
1. 守恒项：in / out / R_d / R_b / N_ex / res_sum / res_split
2. 组成：y_d−y_b、C、holdup 流量
3. 反应源相份额：R_d/(R_d+R_b)（CO/H2）
4. 反事实：把密相净源按体积份额重分到两相后，代数重算 N_ex 与 res_split
   （诊断用，不改求解器）

用法：
  python3 scripts/audit_upper_bed_co_h2_phase_source.py
  python3 scripts/audit_upper_bed_co_h2_phase_source.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
import json
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

from src.core.constants import Rg
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX

ISOMORPH_CFG = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
}

FOCUS_SP = ("CO", "H2", "H2O", "N2", "CO2")
UPPER = (6, 7, 8, 9)


def _inflow(cell: Any, phase: str) -> np.ndarray:
    if phase == "d":
        return np.asarray(cell.N_zu_d, dtype=np.float64) + np.asarray(
            cell.N_rez_d, dtype=np.float64
        ) + np.asarray(cell.N_d_in, dtype=np.float64)
    return np.asarray(cell.N_zu_b, dtype=np.float64) + np.asarray(
        cell.N_rez_b, dtype=np.float64
    ) + np.asarray(cell.N_b_in, dtype=np.float64)


def _y_c(cell: Any, phase: str) -> np.ndarray:
    N = np.asarray(cell.N_d if phase == "d" else cell.N_b, dtype=np.float64)
    tot = float(np.sum(N))
    if tot <= 0.0:
        return np.zeros_like(N)
    return N / tot


def _row(cell: Any, i_cell: int, sp: str) -> dict[str, Any]:
    j = int(GAS_SPECIES_INDEX[sp])
    in_d = float(_inflow(cell, "d")[j])
    in_b = float(_inflow(cell, "b")[j])
    N_d = float(cell.N_d[j])
    N_b = float(cell.N_b[j])
    R_d = float(cell.R_gas_d[j])
    R_b = float(cell.R_gas_b[j])
    N_ex = float(cell.N_ex[j])
    res_d = in_d + R_d - N_d + N_ex
    res_b = in_b + R_b - N_b - N_ex
    y_d = _y_c(cell, "d")
    y_b = _y_c(cell, "b")
    T = float(cell.T)  # [K]
    P = float(cell.P)  # [Pa]
    C_d = float(y_d[j] * P / (Rg * T))  # [mol/m³]
    C_b = float(y_b[j] * P / (Rg * T))  # [mol/m³]
    K_bd = float(cell.K_bd)  # [1/s]
    V_b = float(cell.V_b)  # [m³]
    N_ex_chk = float(K_bd * V_b * (C_b - C_d))
    R_net = R_d + R_b
    dense_frac = float(R_d / R_net) if abs(R_net) > 1e-12 else float("nan")
    # split identity: res_d - res_b = (in_d-in_b) + (R_d-R_b) - (N_d-N_b) + 2 N_ex
    split_parts = {
        "in_diff": in_d - in_b,
        "R_diff": R_d - R_b,
        "out_diff": N_d - N_b,
        "two_N_ex": 2.0 * N_ex,
    }
    split_from_parts = (
        split_parts["in_diff"]
        + split_parts["R_diff"]
        - split_parts["out_diff"]
        + split_parts["two_N_ex"]
    )
    return {
        "cell": i_cell,
        "sp": sp,
        "T": T,
        "in_d": in_d,
        "in_b": in_b,
        "N_d": N_d,
        "N_b": N_b,
        "R_d": R_d,
        "R_b": R_b,
        "dense_source_frac": dense_frac,
        "N_ex": N_ex,
        "N_ex_chk": N_ex_chk,
        "y_d": float(y_d[j]),
        "y_b": float(y_b[j]),
        "dy": float(y_d[j] - y_b[j]),
        "C_d": C_d,
        "C_b": C_b,
        "K_bd": K_bd,
        "V_b": V_b,
        "eps_b": float(cell.eps_b),
        "res_d": res_d,
        "res_b": res_b,
        "res_sum": res_d + res_b,
        "res_split": res_d - res_b,
        "split_parts": split_parts,
        "split_from_parts": split_from_parts,
    }


def _vol_share_counterfactual(cell: Any, sp: str) -> dict[str, Any]:
    """把净源 R_net 按 V_d:V_b 重分，保持总量；用现 y 重算 N_ex（仅诊断）。"""
    j = int(GAS_SPECIES_INDEX[sp])
    R_net = float(cell.R_gas_d[j] + cell.R_gas_b[j])
    V_d = float(cell.V_d)  # [m³]
    V_b = float(cell.V_b)  # [m³]
    V = max(V_d + V_b, 1e-30)
    R_d_cf = R_net * (V_d / V)
    R_b_cf = R_net * (V_b / V)
    in_d = float(_inflow(cell, "d")[j])
    in_b = float(_inflow(cell, "b")[j])
    N_d = float(cell.N_d[j])
    N_b = float(cell.N_b[j])
    # 保持当前浓度 → N_ex 不变；只看源相重分对 split 的直接代数效果
    N_ex = float(cell.N_ex[j])
    res_d = in_d + R_d_cf - N_d + N_ex
    res_b = in_b + R_b_cf - N_b - N_ex
    base_split = float(
        (in_d + float(cell.R_gas_d[j]) - N_d + N_ex)
        - (in_b + float(cell.R_gas_b[j]) - N_b - N_ex)
    )
    return {
        "sp": sp,
        "R_d_cf": R_d_cf,
        "R_b_cf": R_b_cf,
        "res_split_cf": res_d - res_b,
        "res_split_base": base_split,
        "delta_split": (res_d - res_b) - base_split,
        "note": "same N_ex/y; only R phase share changed",
    }


def _extent_proxy_co_h2(cell: Any) -> dict[str, float]:
    """用 cell 已缓存的净源标量做粗分解（非完整 extent 重算）。"""
    # ReactionSourceBundle 标量若挂在 cell 上则读；否则只报 R_gas
    out: dict[str, float] = {
        "R_CO_d": float(cell.R_gas_d[GAS_SPECIES_INDEX["CO"]]),
        "R_CO_b": float(cell.R_gas_b[GAS_SPECIES_INDEX["CO"]]),
        "R_H2_d": float(cell.R_gas_d[GAS_SPECIES_INDEX["H2"]]),
        "R_H2_b": float(cell.R_gas_b[GAS_SPECIES_INDEX["H2"]]),
        "R_H2O_d": float(cell.R_gas_d[GAS_SPECIES_INDEX["H2O"]]),
        "R_H2O_b": float(cell.R_gas_b[GAS_SPECIES_INDEX["H2O"]]),
    }
    for name in (
        "net_molar_gas_source_r5",
        "net_molar_gas_source_r7",
        "net_molar_gas_source_r8",
        "net_molar_gas_source_char",
        "net_molar_gas_source_vm",
    ):
        if hasattr(cell, name):
            out[name] = float(getattr(cell, name))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    args = ap.parse_args()

    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=ISOMORPH_CFG)
    run_phase2_init(reactor)
    solve_kw = phase2_solve_kwargs(max_global_iter=int(args.max_global_iter))
    result, wall_s = timed_phase2_solve(reactor, **solve_kw)

    cells = list(reactor.cells)
    rows = []
    cf_rows = []
    extent_rows = []
    for i in UPPER:
        cell = cells[i]
        for sp in FOCUS_SP:
            rows.append(_row(cell, i, sp))
        for sp in ("CO", "H2"):
            cf = _vol_share_counterfactual(cell, sp)
            cf["cell"] = i
            cf_rows.append(cf)
        extent_rows.append({"cell": i, **_extent_proxy_co_h2(cell)})

    # 上段 split |rs|^2 份额
    ss = sum(r["res_split"] ** 2 for r in rows)
    share = {
        sp: sum(r["res_split"] ** 2 for r in rows if r["sp"] == sp) / ss if ss else 0.0
        for sp in FOCUS_SP
    }

    # CO 上段：split 主导部件均值
    co_u = [r for r in rows if r["sp"] == "CO"]
    parts_mean = {
        k: float(np.mean([r["split_parts"][k] for r in co_u]))
        for k in ("in_diff", "R_diff", "out_diff", "two_N_ex")
    }
    h2_u = [r for r in rows if r["sp"] == "H2"]
    h2_parts_mean = {
        k: float(np.mean([r["split_parts"][k] for r in h2_u]))
        for k in ("in_diff", "R_diff", "out_diff", "two_N_ex")
    }

    exit_gas = {}
    if isinstance(result, dict):
        exit_gas = dict(result.get("bed_exit_gas") or result.get("exit_gas") or {})
    exit_CO = float(exit_gas.get("CO", float("nan")))
    if not np.isfinite(exit_CO):
        top = cells[-1]
        Nt = np.asarray(top.N_d, dtype=np.float64) + np.asarray(top.N_b, dtype=np.float64)
        exit_CO = float(Nt[GAS_SPECIES_INDEX["CO"]] / max(float(np.sum(Nt)), 1e-30))

    metrics = {}
    if isinstance(result, dict):
        metrics = {
            "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
            "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
            "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
            "ls_fail": int((result.get("nr_counts") or {}).get("line_search_failures", -1)),
        }

    payload = {
        "note": "upper-bed CO/H2 phase-source + split-part decomp; isomorphic o8",
        "wall_s": float(wall_s),
        "exit_CO": exit_CO,
        "metrics": metrics,
        "upper_split_sq_share": share,
        "CO_split_parts_mean_bed6_9": parts_mean,
        "H2_split_parts_mean_bed6_9": h2_parts_mean,
        "dense_source_frac_CO_mean": float(np.nanmean([r["dense_source_frac"] for r in co_u])),
        "dense_source_frac_H2_mean": float(np.nanmean([r["dense_source_frac"] for r in h2_u])),
        "vol_share_cf": cf_rows,
        "extent_proxy": extent_rows,
        "rows": rows,
        "cfg": ISOMORPH_CFG,
        "max_global_iter": int(args.max_global_iter),
    }

    out = data_path("phase2_upper_bed_co_h2_phase_source_o8.json")
    write_json(out, payload)

    print(f"wrote {out}")
    print(f"exit_CO≈{payload['exit_CO']:.4f}  wall={wall_s:.1f}s")
    print("upper |rs|^2 share:", json.dumps(share, indent=2))
    print("CO dense_source_frac mean:", payload["dense_source_frac_CO_mean"])
    print("H2 dense_source_frac mean:", payload["dense_source_frac_H2_mean"])
    print("CO split_parts mean:", parts_mean)
    print("H2 split_parts mean:", h2_parts_mean)
    print(f"{'cell':>4} {'sp':>4} {'R_d':>9} {'R_b':>9} {'dens%':>6} {'N_ex':>9} "
          f"{'dy':>9} {'Rdiff':>9} {'2Nex':>9} {'rsplit':>9} {'rsum':>9}")
    for r in rows:
        if r["sp"] not in ("CO", "H2", "N2", "H2O"):
            continue
        p = r["split_parts"]
        df = r["dense_source_frac"]
        dfs = f"{100*df:5.1f}" if np.isfinite(df) else "  nan"
        print(
            f"{r['cell']:4d} {r['sp']:>4} {r['R_d']:9.3f} {r['R_b']:9.3f} {dfs}"
            f" {r['N_ex']:9.3f} {r['dy']:9.5f} {p['R_diff']:9.3f} {p['two_N_ex']:9.3f}"
            f" {r['res_split']:9.3f} {r['res_sum']:9.3f}"
        )
    print("vol-share CF delta_split (CO/H2):")
    for c in cf_rows:
        print(
            f"  bed{c['cell']} {c['sp']}: base={c['res_split_base']:+.3f} "
            f"cf={c['res_split_cf']:+.3f} Δ={c['delta_split']:+.3f}"
        )


if __name__ == "__main__":
    main()
