#!/usr/bin/env python3
"""LU 工况专项审计：drying / pyrolysis 在轴向上持续多久、发生在哪些 cells。

输出两类信息：
1. 每个 cell 的局部 kinetics 视角：tau、x_dry、x_vm、局部释放量
2. 每个 cell 的 solved-state 视角：剩余 moisture / VM / char 组成

目的是回答：
- drying / pyrolysis 大约持续多少停留时间、覆盖多少 cells
- 何处开始进入“上部以 char 为主，几乎无进一步 drying/pyrolysis”的区域
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import S_CHAR, S_MOISTURE, S_VM
from src.core.cell_pyrolysis import calc_drying_pyrolysis_sources
from src.core.reactor import Reactor
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_reactor_config


REL_RELEASE_THRESHOLD = 0.001
REL_REMAIN_THRESHOLD = 0.01
CHAR_DOMINANT_THRESHOLD = 0.95


def _sum_positive(arr: np.ndarray) -> float:
    return float(np.sum(np.maximum(arr, 0.0)))


def _cell_row(cell, cfg, moist_feed: float, vm_feed: float) -> dict[str, Any]:
    cell.calc_hydrodynamics()
    tau = float(cell.geo.dh / max(cell.u_mf, 1e-3))

    m_vm_in = _sum_positive(cell.m_solid_zu[:, S_VM] + cell.m_solid_rez[:, S_VM] + cell.m_solid_in[:, S_VM])
    m_moist_in = _sum_positive(cell.m_solid_zu[:, S_MOISTURE] + cell.m_solid_rez[:, S_MOISTURE] + cell.m_solid_in[:, S_MOISTURE])

    bundle = calc_drying_pyrolysis_sources(
        tau=tau,
        T=float(cell.T),
        P=float(cell.P),
        T_init=float(cell.T_in_solid),
        d_p=cell.solid.d_p,
        moisture_wt=cell.solid.moisture_wt,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
        pyrolysis_tar_carbon_frac=cell.solid.pyrolysis_tar_carbon_frac,
        fuel_type=cell.fuel_type,
        m_vm_in=m_vm_in,
        m_moist_in=m_moist_in,
        solid_shape=cell.R_solid.shape,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )

    vm_release = float(max(-np.sum(bundle.solid_sink[:, S_VM]), 0.0))
    moist_release = float(max(-np.sum(bundle.solid_sink[:, S_MOISTURE]), 0.0))
    local_drying_frac = moist_release / max(m_moist_in, 1e-12)
    local_pyrolysis_frac = vm_release / max(m_vm_in, 1e-12)

    m_vm_out = _sum_positive(cell.m_solid[:, S_VM])
    m_moist_out = _sum_positive(cell.m_solid[:, S_MOISTURE])
    m_char_out = _sum_positive(cell.m_solid[:, S_CHAR])
    total_out = m_vm_out + m_moist_out + m_char_out
    char_share = m_char_out / max(total_out, 1e-12)

    drying_active = (
        local_drying_frac > REL_RELEASE_THRESHOLD
        or moist_release > REL_RELEASE_THRESHOLD * max(moist_feed, 1e-12)
    )
    pyrolysis_active = (
        local_pyrolysis_frac > REL_RELEASE_THRESHOLD
        or vm_release > REL_RELEASE_THRESHOLD * max(vm_feed, 1e-12)
    )
    char_dominant = (
        char_share >= CHAR_DOMINANT_THRESHOLD
        and moist_release <= REL_RELEASE_THRESHOLD * max(moist_feed, 1e-12)
        and vm_release <= REL_RELEASE_THRESHOLD * max(vm_feed, 1e-12)
    )

    return {
        "cell": int(round(cell.geo.h_center / cell.geo.dh - 0.5)),
        "xi": float(cell.geo.h_center / cfg.H_bed),
        "h_m": float(cell.geo.h_center),
        "T_K": float(cell.T),
        "tau_s": tau,
        "u_mf_m_s": float(cell.u_mf),
        "x_dry": float(bundle.x_dry),
        "x_vm": float(bundle.x_vm),
        "m_moist_in_kg_s": m_moist_in,
        "m_vm_in_kg_s": m_vm_in,
        "moist_release_kg_s": moist_release,
        "vm_release_kg_s": vm_release,
        "local_drying_frac": local_drying_frac,
        "local_pyrolysis_frac": local_pyrolysis_frac,
        "m_moist_out_kg_s": m_moist_out,
        "m_vm_out_kg_s": m_vm_out,
        "m_char_out_kg_s": m_char_out,
        "char_share_out": char_share,
        "drying_active": bool(drying_active),
        "pyrolysis_active": bool(pyrolysis_active),
        "char_dominant": bool(char_dominant),
    }


def _first_suffix(rows: list[dict[str, Any]], key: str) -> int | None:
    for i in range(len(rows)):
        if all(bool(row[key]) for row in rows[i:]):
            return i
    return None


def run_audit() -> dict[str, Any]:
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    bot = reactor.cells[0]
    moist_feed = _sum_positive(bot.m_solid_zu[:, S_MOISTURE])
    vm_feed = _sum_positive(bot.m_solid_zu[:, S_VM])

    rows = [_cell_row(cell, cfg, moist_feed, vm_feed) for cell in reactor.cells]

    drying_cells = [row["cell"] for row in rows if row["drying_active"]]
    pyrolysis_cells = [row["cell"] for row in rows if row["pyrolysis_active"]]
    char_suffix = _first_suffix(rows, "char_dominant")
    vm_done_cell = next((row["cell"] for row in rows if row["m_vm_out_kg_s"] <= REL_REMAIN_THRESHOLD * max(vm_feed, 1e-12)), None)
    moist_done_cell = next((row["cell"] for row in rows if row["m_moist_out_kg_s"] <= REL_REMAIN_THRESHOLD * max(moist_feed, 1e-12)), None)

    drying_tau_total = float(sum(row["tau_s"] for row in rows if row["drying_active"]))
    pyrolysis_tau_total = float(sum(row["tau_s"] for row in rows if row["pyrolysis_active"]))

    if drying_cells:
        drying_zone = {"start_cell": int(drying_cells[0]), "end_cell": int(drying_cells[-1]), "n_cells": int(len(drying_cells))}
    else:
        drying_zone = {"start_cell": None, "end_cell": None, "n_cells": 0}

    if pyrolysis_cells:
        pyrolysis_zone = {"start_cell": int(pyrolysis_cells[0]), "end_cell": int(pyrolysis_cells[-1]), "n_cells": int(len(pyrolysis_cells))}
    else:
        pyrolysis_zone = {"start_cell": None, "end_cell": None, "n_cells": 0}

    char_upper = {
        "from_cell": int(char_suffix) if char_suffix is not None else None,
        "n_cells": int(len(rows) - char_suffix) if char_suffix is not None else 0,
        "all_upper_cells_char_dominant": bool(char_suffix is not None),
    }

    return {
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "drying_zone": drying_zone,
        "pyrolysis_zone": pyrolysis_zone,
        "drying_tau_total_s": drying_tau_total,
        "pyrolysis_tau_total_s": pyrolysis_tau_total,
        "char_dominant_upper_zone": char_upper,
        "fresh_feed_completion": {
            "vm_done_cell": vm_done_cell,
            "moisture_done_cell": moist_done_cell,
        },
        "thresholds": {
            "relative_release_threshold": REL_RELEASE_THRESHOLD,
            "relative_remaining_threshold": REL_REMAIN_THRESHOLD,
            "char_dominant_threshold": CHAR_DOMINANT_THRESHOLD,
        },
        "rows": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="审计 LU 工况 drying / pyrolysis 的持续范围与上部 char 主导区")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    args = ap.parse_args()

    out = run_audit()

    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=2))
    else:
        print("=" * 132)
        print("LU drying / pyrolysis extent audit")
        print("=" * 132)
        print(
            f"converged={out['converged']} n_iter={out['n_iter']} "
            f"drying_cells={out['drying_zone']['n_cells']} tau={out['drying_tau_total_s']:.3f}s "
            f"pyro_cells={out['pyrolysis_zone']['n_cells']} tau={out['pyrolysis_tau_total_s']:.3f}s"
        )
        print(
            f"fresh-feed completion: vm_done_cell={out['fresh_feed_completion']['vm_done_cell']} "
            f"moisture_done_cell={out['fresh_feed_completion']['moisture_done_cell']}"
        )
        char_upper = out["char_dominant_upper_zone"]
        print(
            f"char-dominant upper zone: from_cell={char_upper['from_cell']} "
            f"n_cells={char_upper['n_cells']} all_upper_cells_char_dominant={char_upper['all_upper_cells_char_dominant']}"
        )
        print("-" * 132)
        print(
            f"{'cell':>4} {'xi':>6} {'T(K)':>8} {'tau(s)':>8} {'x_dry':>8} {'x_vm':>8} "
            f"{'dry_loc':>8} {'pyr_loc':>8} {'moist_rel':>10} {'vm_rel':>10} {'m_moist':>10} {'m_vm':>10} {'char%':>8} "
            f"{'dry?':>6} {'pyr?':>6} {'char?':>6}"
        )
        for row in out["rows"]:
            print(
                f"{row['cell']:>4d} {row['xi']:>6.3f} {row['T_K']:>8.1f} {row['tau_s']:>8.3f} "
                f"{row['x_dry']:>8.3f} {row['x_vm']:>8.3f} {row['local_drying_frac']:>8.3f} {row['local_pyrolysis_frac']:>8.3f} "
                f"{row['moist_release_kg_s']:>10.4e} "
                f"{row['vm_release_kg_s']:>10.4e} {row['m_moist_out_kg_s']:>10.4e} "
                f"{row['m_vm_out_kg_s']:>10.4e} {100.0 * row['char_share_out']:>7.2f}% "
                f"{str(row['drying_active']):>6} {str(row['pyrolysis_active']):>6} {str(row['char_dominant']):>6}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
