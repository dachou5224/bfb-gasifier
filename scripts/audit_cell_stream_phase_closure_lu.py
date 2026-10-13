#!/usr/bin/env python3
"""逐格进出物流 + 乳化相/气泡相摩尔闭合审计（Phase2 LU）。

口径（Hamel Eq.2-4 / 2-5；码内 ``N_ex = -ex_bd``）：
  悬浮相: 0 = N_zu_d + N_rez_d + N_d_in + R_d - N_d + N_ex
  气泡相: 0 = N_zu_b + N_rez_b + N_b_in + R_b - N_b - N_ex
  两相合计: 0 = Σ(zu+rez+in+R) - Σ(N_d+N_b)   （N_ex 相消）

轴向气链：bed 内 ``curr.N_*_in == prev.N_*``；FB 可退化为单相（总和进 N_d_in）。

用法：
  PYTHONPATH=scripts python3 scripts/audit_cell_stream_phase_closure_lu.py
  PYTHONPATH=scripts python3 scripts/audit_cell_stream_phase_closure_lu.py --max-global-iter 8
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

from _lu_harness import (
    build_phase2_reactor,
    phase2_solve_kwargs,
    timed_phase2_solve,
    write_json,
)
from src.core.cell import S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.species import GAS_SPECIES, N_GAS


MAJOR = ("CO", "CO2", "H2", "H2O", "O2", "CH4", "N2")


def _sp_idx(name: str) -> int:
    return GAS_SPECIES.index(name)


def _sum_arr(a: np.ndarray) -> float:
    return float(np.sum(np.asarray(a, dtype=np.float64)))


def _gas_terms(cell: Any) -> dict[str, np.ndarray]:
    return {
        "zu_d": np.asarray(cell.N_zu_d, dtype=np.float64),
        "rez_d": np.asarray(cell.N_rez_d, dtype=np.float64),
        "in_d": np.asarray(cell.N_d_in, dtype=np.float64),
        "R_d": np.asarray(cell.R_gas_d, dtype=np.float64),
        "out_d": np.asarray(cell.N_d, dtype=np.float64),
        "zu_b": np.asarray(cell.N_zu_b, dtype=np.float64),
        "rez_b": np.asarray(cell.N_rez_b, dtype=np.float64),
        "in_b": np.asarray(cell.N_b_in, dtype=np.float64),
        "R_b": np.asarray(cell.R_gas_b, dtype=np.float64),
        "out_b": np.asarray(cell.N_b, dtype=np.float64),
        "N_ex": np.asarray(cell.N_ex, dtype=np.float64),
    }


def _phase_residuals(t: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """返回 (res_d, res_b, res_combined) [mol/s]。"""
    res_d = t["zu_d"] + t["rez_d"] + t["in_d"] + t["R_d"] - t["out_d"] + t["N_ex"]
    res_b = t["zu_b"] + t["rez_b"] + t["in_b"] + t["R_b"] - t["out_b"] - t["N_ex"]
    res_c = (
        t["zu_d"]
        + t["zu_b"]
        + t["rez_d"]
        + t["rez_b"]
        + t["in_d"]
        + t["in_b"]
        + t["R_d"]
        + t["R_b"]
        - t["out_d"]
        - t["out_b"]
    )
    return res_d, res_b, res_c


def _major_dict(vec: np.ndarray) -> dict[str, float]:
    return {sp: float(vec[_sp_idx(sp)]) for sp in MAJOR if sp in GAS_SPECIES}


def audit_cell_gas(cell: Any, i: int) -> dict[str, Any]:
    t = _gas_terms(cell)
    res_d, res_b, res_c = _phase_residuals(t)
    # 与 cell.calc_gas_balance 交叉核对
    res_pkg = np.asarray(cell.calc_gas_balance(), dtype=np.float64)
    assert res_pkg.shape == (2 * N_GAS,)
    pkg_d, pkg_b = res_pkg[:N_GAS], res_pkg[N_GAS:]
    inflow = t["zu_d"] + t["zu_b"] + t["rez_d"] + t["rez_b"] + t["in_d"] + t["in_b"]
    outflow = t["out_d"] + t["out_b"]
    return {
        "i": i,
        "cell_type": str(cell.cell_type),
        "T_K": float(cell.T),
        "totals_mol_s": {
            "zu_d": _sum_arr(t["zu_d"]),
            "zu_b": _sum_arr(t["zu_b"]),
            "rez_d": _sum_arr(t["rez_d"]),
            "rez_b": _sum_arr(t["rez_b"]),
            "in_d": _sum_arr(t["in_d"]),
            "in_b": _sum_arr(t["in_b"]),
            "R_d": _sum_arr(t["R_d"]),
            "R_b": _sum_arr(t["R_b"]),
            "out_d": _sum_arr(t["out_d"]),
            "out_b": _sum_arr(t["out_b"]),
            "N_ex_net": _sum_arr(t["N_ex"]),
            "inflow": _sum_arr(inflow),
            "outflow": _sum_arr(outflow),
            "R_total": _sum_arr(t["R_d"] + t["R_b"]),
        },
        "residual_mol_s": {
            "dense_L2": float(np.linalg.norm(res_d)),
            "bubble_L2": float(np.linalg.norm(res_b)),
            "combined_L2": float(np.linalg.norm(res_c)),
            "dense_max_abs": float(np.max(np.abs(res_d))),
            "bubble_max_abs": float(np.max(np.abs(res_b))),
            "combined_max_abs": float(np.max(np.abs(res_c))),
            "dense_max_sp": GAS_SPECIES[int(np.argmax(np.abs(res_d)))],
            "bubble_max_sp": GAS_SPECIES[int(np.argmax(np.abs(res_b)))],
            "pkg_vs_manual_dense_max": float(np.max(np.abs(pkg_d - res_d))),
            "pkg_vs_manual_bubble_max": float(np.max(np.abs(pkg_b - res_b))),
        },
        "major_out_d": _major_dict(t["out_d"]),
        "major_out_b": _major_dict(t["out_b"]),
        "major_N_ex": _major_dict(t["N_ex"]),
        "major_res_d": _major_dict(res_d),
        "major_res_b": _major_dict(res_b),
        "major_res_combined": _major_dict(res_c),
    }


def audit_cell_solid(cell: Any, i: int) -> dict[str, Any]:
    model = str(getattr(cell, "solid_state_model", "legacy_stream"))
    zu = np.asarray(cell.m_solid_zu, dtype=np.float64)
    rez = np.asarray(cell.m_solid_rez, dtype=np.float64)
    inn = np.asarray(cell.m_solid_in, dtype=np.float64)
    auf = np.asarray(cell.m_solid_auf_in, dtype=np.float64)
    ab = np.asarray(cell.m_solid_ab_in, dtype=np.float64)
    R = np.asarray(cell.R_solid, dtype=np.float64)
    m = np.asarray(cell.m_solid, dtype=np.float64)
    mout = np.asarray(cell._solid_outflow_rates(), dtype=np.float64)
    mig = np.asarray(cell._calc_size_migration(), dtype=np.float64)
    res = np.asarray(cell.calc_solid_balance(), dtype=np.float64)

    def _comp_sum(arr: np.ndarray, col: int) -> float:
        return float(np.sum(arr[:, col]))

    return {
        "i": i,
        "cell_type": str(cell.cell_type),
        "solid_state_model": model,
        "char_kg_or_kg_s": {
            "zu": _comp_sum(zu, S_CHAR),
            "rez": _comp_sum(rez, S_CHAR),
            "in_vm_moist_channel": _comp_sum(inn, S_CHAR),
            "auf_in": _comp_sum(auf, S_CHAR),
            "ab_in": _comp_sum(ab, S_CHAR),
            "R": _comp_sum(R, S_CHAR),
            "migration": _comp_sum(mig, S_CHAR),
            "holdup_or_out_state": _comp_sum(m, S_CHAR),
            "outflow_Km": _comp_sum(mout, S_CHAR),
            "ash_holdup": _comp_sum(m, S_ASH),
            "vm_holdup": _comp_sum(m, S_VM),
            "moist_holdup": _comp_sum(m, S_MOISTURE),
        },
        "residual": {
            "L2": float(np.linalg.norm(res)),
            "max_abs": float(np.max(np.abs(res))),
            "char_sum": _comp_sum(res, S_CHAR),
            "ash_sum": _comp_sum(res, S_ASH),
            "vm_sum": _comp_sum(res, S_VM),
            "moist_sum": _comp_sum(res, S_MOISTURE),
        },
    }


def audit_axial_gas_chain(cells: list[Any]) -> list[dict[str, Any]]:
    """检查相邻单元气相入口是否等于上游出口（按路由语义）。"""
    rows: list[dict[str, Any]] = []
    for i in range(1, len(cells)):
        prev, curr = cells[i - 1], cells[i]
        ctype = str(curr.cell_type)
        # bed→bed: 保相；进入 FB/side：可能退化为单相总和进 dense
        gap_b = np.asarray(curr.N_b_in, dtype=np.float64) - np.asarray(prev.N_b, dtype=np.float64)
        gap_d = np.asarray(curr.N_d_in, dtype=np.float64) - np.asarray(prev.N_d, dtype=np.float64)
        prev_tot = np.asarray(prev.N_b, dtype=np.float64) + np.asarray(prev.N_d, dtype=np.float64)
        curr_tot_in = np.asarray(curr.N_b_in, dtype=np.float64) + np.asarray(
            curr.N_d_in, dtype=np.float64
        )
        gap_tot = curr_tot_in - prev_tot
        phase_preserve = bool(
            str(prev.cell_type) == "bed" and ctype == "bed"
        )
        rows.append(
            {
                "from_i": i - 1,
                "to_i": i,
                "from_type": str(prev.cell_type),
                "to_type": ctype,
                "phase_preserve_expected": phase_preserve,
                "gap_b_L2": float(np.linalg.norm(gap_b)),
                "gap_d_L2": float(np.linalg.norm(gap_d)),
                "gap_total_L2": float(np.linalg.norm(gap_tot)),
                "gap_total_max_abs": float(np.max(np.abs(gap_tot))),
                "gap_b_max_abs": float(np.max(np.abs(gap_b))),
                "gap_d_max_abs": float(np.max(np.abs(gap_d))),
                "prev_out_total": float(np.sum(prev_tot)),
                "curr_in_total": float(np.sum(curr_tot_in)),
            }
        )
    return rows


def audit_axial_solid_char_chain(cells: list[Any]) -> list[dict[str, Any]]:
    """char 轴向：上游 auf 出流 vs 下游 auf_in（holdup_transport）。"""
    rows: list[dict[str, Any]] = []
    for i in range(1, len(cells)):
        prev, curr = cells[i - 1], cells[i]
        if str(prev.cell_type) != "bed" or str(curr.cell_type) != "bed":
            continue
        prev_auf = float(np.sum(prev._solid_upflow_rates()[:, S_CHAR]))
        curr_auf_in = float(np.sum(np.maximum(curr.m_solid_auf_in[:, S_CHAR], 0.0)))
        prev_ab = float(np.sum(prev._solid_downflow_rates()[:, S_CHAR])) if hasattr(prev, "_solid_downflow_rates") else float(
            np.sum(np.asarray(prev.K_solid_ab, dtype=np.float64)[:, S_CHAR] * np.maximum(prev.m_solid[:, S_CHAR], 0.0))
        )
        # 下游 ab 来自更上层；此处只报本格 ab_in vs 上层 ab_out（若存在）
        above = cells[i + 1] if i + 1 < len(cells) and str(cells[i + 1].cell_type) == "bed" else None
        above_ab = (
            float(
                np.sum(
                    np.asarray(above.K_solid_ab, dtype=np.float64)[:, S_CHAR]
                    * np.maximum(above.m_solid[:, S_CHAR], 0.0)
                )
            )
            if above is not None
            else float("nan")
        )
        curr_ab_in = float(np.sum(np.maximum(curr.m_solid_ab_in[:, S_CHAR], 0.0)))
        rows.append(
            {
                "from_i": i - 1,
                "to_i": i,
                "prev_auf_out": prev_auf,
                "curr_auf_in": curr_auf_in,
                "gap_auf": curr_auf_in - prev_auf,
                "prev_ab_out": prev_ab,
                "curr_ab_in": curr_ab_in,
                "above_ab_out": above_ab,
                "gap_ab_vs_above": curr_ab_in - above_ab if above is not None else float("nan"),
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=8)
    ap.add_argument("--ladder", default="p2")
    ap.add_argument(
        "--out",
        default="data/phase2_cell_stream_phase_closure_o8.json",
    )
    ap.add_argument("--gas-tol", type=float, default=1e-4, help="mol/s max|res| 警戒阈值")
    ap.add_argument("--solid-tol", type=float, default=1e-3, help="kg/s max|res| 警戒阈值")
    ap.add_argument("--chain-tol", type=float, default=1e-6, help="轴向气链 gap mol/s")
    args = ap.parse_args()

    r = build_phase2_reactor(ladder_level=args.ladder)
    out, wall = timed_phase2_solve(
        r, **phase2_solve_kwargs(max_global_iter=args.max_global_iter, verbose=False)
    )
    ex = {
        k: float((out.get("exit_gas_dry") or {}).get(k, 0) or 0)
        for k in ("CO", "CO2", "H2", "O2")
    }
    print(f"exit {ex} wall={wall:.1f}s", flush=True)

    gas_rows = [audit_cell_gas(c, i) for i, c in enumerate(r.cells)]
    solid_rows = [audit_cell_solid(c, i) for i, c in enumerate(r.cells)]
    gas_chain = audit_axial_gas_chain(r.cells)
    solid_chain = audit_axial_solid_char_chain(r.cells)

    gas_flags = [
        {
            "i": row["i"],
            "type": row["cell_type"],
            "dense_max": row["residual_mol_s"]["dense_max_abs"],
            "bubble_max": row["residual_mol_s"]["bubble_max_abs"],
            "combined_max": row["residual_mol_s"]["combined_max_abs"],
        }
        for row in gas_rows
        if max(
            row["residual_mol_s"]["dense_max_abs"],
            row["residual_mol_s"]["bubble_max_abs"],
            row["residual_mol_s"]["combined_max_abs"],
        )
        > args.gas_tol
    ]
    solid_flags = [
        {
            "i": row["i"],
            "type": row["cell_type"],
            "max_abs": row["residual"]["max_abs"],
            "char_sum": row["residual"]["char_sum"],
        }
        for row in solid_rows
        if row["residual"]["max_abs"] > args.solid_tol
        and str(row["solid_state_model"]) != "freeboard_closure"
    ]
    chain_flags = [
        row
        for row in gas_chain
        if (
            row["phase_preserve_expected"]
            and max(row["gap_b_max_abs"], row["gap_d_max_abs"]) > args.chain_tol
        )
        or ((not row["phase_preserve_expected"]) and row["gap_total_max_abs"] > args.chain_tol)
    ]
    solid_chain_flags = [
        row
        for row in solid_chain
        if abs(row["gap_auf"]) > args.solid_tol
        or (np.isfinite(row["gap_ab_vs_above"]) and abs(row["gap_ab_vs_above"]) > args.solid_tol)
    ]

    # 打印摘要表
    print("\n=== GAS phase residuals (mol/s) ===", flush=True)
    print(
        f"{'i':>3} {'type':>10} {'in':>10} {'out':>10} {'R':>10} "
        f"{'|rd|':>10} {'|rb|':>10} {'|rc|':>10} {'Nex':>10}",
        flush=True,
    )
    for row in gas_rows:
        t = row["totals_mol_s"]
        rsd = row["residual_mol_s"]
        print(
            f"{row['i']:3d} {row['cell_type']:>10} {t['inflow']:10.4f} {t['outflow']:10.4f} "
            f"{t['R_total']:10.4f} {rsd['dense_max_abs']:10.2e} {rsd['bubble_max_abs']:10.2e} "
            f"{rsd['combined_max_abs']:10.2e} {t['N_ex_net']:10.4f}",
            flush=True,
        )

    print("\n=== SOLID char residuals ===", flush=True)
    print(
        f"{'i':>3} {'type':>10} {'auf_in':>10} {'ab_in':>10} {'out':>10} "
        f"{'R':>10} {'|res|':>10} {'char_res':>10}",
        flush=True,
    )
    for row in solid_rows:
        c = row["char_kg_or_kg_s"]
        print(
            f"{row['i']:3d} {row['cell_type']:>10} {c['auf_in']:10.4f} {c['ab_in']:10.4f} "
            f"{c['outflow_Km']:10.4f} {c['R']:10.4f} {row['residual']['max_abs']:10.2e} "
            f"{row['residual']['char_sum']:10.2e}",
            flush=True,
        )

    print("\n=== Axial GAS chain gaps ===", flush=True)
    for row in gas_chain:
        mark = "OK" if row not in chain_flags else "GAP"
        print(
            f"{mark} {row['from_i']}({row['from_type']})→{row['to_i']}({row['to_type']}): "
            f"gap_b={row['gap_b_L2']:.2e} gap_d={row['gap_d_L2']:.2e} "
            f"gap_tot={row['gap_total_L2']:.2e} phase_preserve={row['phase_preserve_expected']}",
            flush=True,
        )

    print("\n=== Axial SOLID char auf/ab gaps (bed) ===", flush=True)
    for row in solid_chain:
        mark = "OK" if row not in solid_chain_flags else "GAP"
        print(
            f"{mark} bed{row['from_i']}→{row['to_i']}: "
            f"auf {row['prev_auf_out']:.4f}→{row['curr_auf_in']:.4f} Δ={row['gap_auf']:.2e}  "
            f"ab_in={row['curr_ab_in']:.4f} above_ab={row['above_ab_out']}",
            flush=True,
        )

    payload = {
        "knife": "cell_stream_phase_closure",
        "ladder": args.ladder,
        "max_global_iter": args.max_global_iter,
        "wall_s": wall,
        "exit_dry": ex,
        "tol": {"gas_mol_s": args.gas_tol, "solid_kg_s": args.solid_tol, "chain_mol_s": args.chain_tol},
        "formula": {
            "dense": "zu_d+rez_d+in_d+R_d-N_d+N_ex = 0",
            "bubble": "zu_b+rez_b+in_b+R_b-N_b-N_ex = 0",
            "N_ex_sign": "code N_ex = -ex_bd (positive bubble→dense)",
            "axial_bed_gas": "curr.N_b_in=prev.N_b; curr.N_d_in=prev.N_d",
        },
        "gas_cells": gas_rows,
        "solid_cells": solid_rows,
        "axial_gas_chain": gas_chain,
        "axial_solid_char_chain": solid_chain,
        "flags": {
            "gas_residual_over_tol": gas_flags,
            "solid_residual_over_tol": solid_flags,
            "gas_chain_gaps": chain_flags,
            "solid_chain_gaps": solid_chain_flags,
        },
        "summary": {
            "n_gas_flags": len(gas_flags),
            "n_solid_flags": len(solid_flags),
            "n_gas_chain_flags": len(chain_flags),
            "n_solid_chain_flags": len(solid_chain_flags),
            "worst_gas_combined_max": max(
                (row["residual_mol_s"]["combined_max_abs"] for row in gas_rows), default=0.0
            ),
            "worst_gas_phase_split_max": max(
                (
                    max(
                        row["residual_mol_s"]["dense_max_abs"],
                        row["residual_mol_s"]["bubble_max_abs"],
                    )
                    for row in gas_rows
                ),
                default=0.0,
            ),
            "worst_solid_max": max(
                (
                    row["residual"]["max_abs"]
                    for row in solid_rows
                    if row["solid_state_model"] != "freeboard_closure"
                ),
                default=0.0,
            ),
        },
    }
    write_json(args.out, payload)
    print(f"\nflags: gas={len(gas_flags)} solid={len(solid_flags)} "
          f"gas_chain={len(chain_flags)} solid_chain={len(solid_chain_flags)}")
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
