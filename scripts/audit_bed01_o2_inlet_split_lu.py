#!/usr/bin/env python3
"""Phase2 LU：bed0/bed1 O₂ 进料与相分配审计。

追踪底部两格 O₂ 从 ``N_zu``（主进料）、``N_rez``（循环）、``N_*_in``（轴向）
到 holdup ``N_d/N_b`` 的全链路，并对照 R5 氧化可用 O₂ 预算。

Usage
-----
    python3 scripts/audit_bed01_o2_inlet_split_lu.py
    python3 scripts/audit_bed01_o2_inlet_split_lu.py --output data/phase2_bed01_o2_split_audit.json
"""

from __future__ import annotations

from _lu_harness import (
    add_max_global_iter_arg,
    build_phase2_config,
    data_path,
    load_case_lu,
    phase2_solve_kwargs,
    repo_root,
    run_phase2_init,
    write_json,
)

REPO_ROOT = repo_root()

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from src.core.connectivity import resolve_bottom_gas_inlet_dense_fraction
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor

def _o2_row(cell, *, label: str, cell_index: int, cfg) -> dict[str, float | str | int]:
    """单 cell O₂ 链路快照。"""
    j = idx["O2"]
    n_ex_o2 = float(cell.N_ex[j])
    o2_zu_d = float(cell.N_zu_d[j])
    o2_zu_b = float(cell.N_zu_b[j])
    o2_d_in = float(cell.N_d_in[j])
    o2_b_in = float(cell.N_b_in[j])
    o2_rez_d = float(cell.N_rez_d[j])
    o2_rez_b = float(cell.N_rez_b[j])
    o2_d = float(cell.N_d[j])
    o2_b = float(cell.N_b[j])
    o2_total_hold = o2_d + o2_b

    # 与 cell_kinetics O₂ limiter 同结构
    o2_supply_b = max(o2_zu_b + o2_b_in + max(-n_ex_o2, 0.0) + o2_rez_b, 0.0)
    o2_supply_d = max(o2_zu_d + o2_d_in + max(n_ex_o2, 0.0) + o2_rez_d, 0.0)

    cell.calc_hydrodynamics()
    tau_val = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau_val)
    cell.calc_exchange()
    T, P = cell.T, cell.P
    C_b, C_d = cell._concentrations("b"), cell._concentrations("d")
    y_b, y_d = cell._mole_fractions("b"), cell._mole_fractions("d")
    r5b_raw = float(cell.r5_scale) * rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], C_b[idx["H2O"]], P, y_b, idx) * cell.V_b
    r5d_raw = float(cell.r5_scale) * rate_R5_suspension(
        T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx
    ) * cell.V_d
    # R5: CO + 1/2 O2 → CO2；extent = mol CO/s → O₂ demand = 0.5·r5
    bubble_o2_demand = 0.5 * r5b_raw
    dense_o2_demand = 0.5 * r5d_raw
    limit_o2_b = min(1.0, 0.995 * o2_supply_b / max(bubble_o2_demand, 1e-9))
    limit_o2_d = min(1.0, 0.995 * o2_supply_d / max(dense_o2_demand, 1e-9))

    y_combined = cell._mole_fractions("combined")
    dense_frac_hold = o2_d / max(o2_total_hold, 1e-12)

    out: dict[str, float | str | int] = {
        "label": label,
        "cell_index": cell_index,
        "xi_over_H": float(cell.geo.h_center / max(float(cfg.H_bed), 1e-12)),
        "T_K": float(T),
        "O2_feed_cfg": float(cfg.O2_feed) if cell_index == 0 else 0.0,
        "dense_frac_inlet_resolved": float(resolve_bottom_gas_inlet_dense_fraction(cell, cfg)) if cell_index == 0 else float("nan"),
        "vorab_frozen_dense_frac": float(getattr(cell, "_vorab_bottom_gas_inlet_dense_frac", float("nan"))),
        "u0_m_s": float(getattr(cell, "u0", 0.0)),
        "u_d_m_s": float(getattr(cell, "u_d", 0.0)),
        "eps_b": float(getattr(cell, "eps_b", 0.0)),
        "N_zu_d_O2": o2_zu_d,
        "N_zu_b_O2": o2_zu_b,
        "N_d_in_O2": o2_d_in,
        "N_b_in_O2": o2_b_in,
        "N_rez_d_O2": o2_rez_d,
        "N_rez_b_O2": o2_rez_b,
        "N_ex_O2": n_ex_o2,
        "N_d_O2": o2_d,
        "N_b_O2": o2_b,
        "y_O2_combined": float(y_combined[idx["O2"]]),
        "y_O2_bubble": float(y_b[idx["O2"]]),
        "y_O2_dense": float(y_d[idx["O2"]]),
        "y_CO_combined": float(y_combined[idx["CO"]]),
        "y_CO2_combined": float(y_combined[idx["CO2"]]),
        "holdup_dense_frac_O2": dense_frac_hold,
        "o2_supply_bubble_mol_s": o2_supply_b,
        "o2_supply_dense_mol_s": o2_supply_d,
        "r5b_raw_mol_s": r5b_raw,
        "r5d_raw_mol_s": r5d_raw,
        "limit_o2_bubble": limit_o2_b,
        "limit_o2_dense": limit_o2_d,
        "K_bd": float(getattr(cell, "K_bd", 0.0)),
    }
    return out


def _snapshot_bed01(reactor: Reactor, *, label: str) -> list[dict]:
    cfg = reactor.config
    rows = []
    for i in (0, 1):
        if i >= len(reactor.cells):
            break
        rows.append(_o2_row(reactor.cells[i], label=label, cell_index=i, cfg=cfg))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "phase2_bed01_o2_split_audit.json",
    )
    parser.add_argument("--max-global-iter", type=int, default=25)
    args = parser.parse_args()

    case = load_case_lu()
    cfg = build_phase2_config(case)
    cfg.vorab_major_gibbs_x0 = False
    reactor = Reactor(cfg)

    run_phase2_init(reactor, init_strategy="vorabrechnung")
    init_rows = _snapshot_bed01(reactor, label="after_init_precalc")

    kw = phase2_solve_kwargs()
    kw["max_global_iter"] = int(args.max_global_iter)
    kw["verbose"] = False
    result, monitor = solve_with_nr_monitor(reactor, kw, check_x0=False)
    final_rows = _snapshot_bed01(reactor, label="after_nr_converged")

    bed0 = reactor.cells[0]
    dense_resolved = resolve_bottom_gas_inlet_dense_fraction(bed0, cfg)
    payload = {
        "case_key": "CASE_HTW_WESSELING_1",
        "config": {
            "O2_feed_mol_s": float(cfg.O2_feed),
            "H2O_feed_mol_s": float(cfg.H2O_feed),
            "N2_feed_mol_s": float(cfg.N2_feed),
            "gas_inlet_dense_frac_cfg": float(cfg.gas_inlet_dense_frac),
            "gas_inlet_split_strategy": str(cfg.gas_inlet_split_strategy),
            "recirculation_frac": float(cfg.recirculation_frac),
            "recycle_gas": bool(cfg.recycle_gas),
            "r5_scale": float(cfg.r5_scale),
        },
        "converged_fully": bool(result.get("converged_fully")),
        "rms_scaled_final": float(result.get("rms_scaled_final", 0.0) or 0.0),
        "bed_exit_gas_dry": {
            sp: float((result.get("bed_exit_gas_dry") or {}).get(sp, 0.0))
            for sp in ("CO", "CO2", "O2")
        },
        "dense_frac_resolved_bed0": float(dense_resolved),
        "init": init_rows,
        "final": final_rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=" * 120)
    print("Phase2 LU bed0/bed1 O2 inlet split audit")
    print("=" * 120)
    print_nr_monitor(monitor, prefix="NR")
    print(
        f"O2_feed={cfg.O2_feed:.4f} mol/s  strategy={cfg.gas_inlet_split_strategy!r}  "
        f"cfg_dense_frac={cfg.gas_inlet_dense_frac:.3f}  resolved_dense_frac={dense_resolved:.3f}  "
        f"recirc={cfg.recirculation_frac:.2f}  conv={payload['converged_fully']}"
    )
    print("-" * 120)
    hdr = (
        f"{'stage':>12} {'cell':>4} {'xi/H':>5} {'N_zu_d':>8} {'N_zu_b':>8} "
        f"{'N_d_in':>8} {'N_b_in':>8} {'N_rez_d':>8} {'N_rez_b':>8} "
        f"{'N_d':>8} {'N_b':>8} {'yO2':>7} {'yCO':>7} {'lim_O2b':>7} {'lim_O2d':>7} {'R5b':>9} {'R5d':>9}"
    )
    print(hdr)
    print("-" * 120)
    for block, rows in (("init", init_rows), ("final", final_rows)):
        for row in rows:
            print(
                f"{block:>12} {row['cell_index']:>4d} {row['xi_over_H']:>5.2f} "
                f"{row['N_zu_d_O2']:>8.4f} {row['N_zu_b_O2']:>8.4f} "
                f"{row['N_d_in_O2']:>8.4f} {row['N_b_in_O2']:>8.4f} "
                f"{row['N_rez_d_O2']:>8.4f} {row['N_rez_b_O2']:>8.4f} "
                f"{row['N_d_O2']:>8.4f} {row['N_b_O2']:>8.4f} "
                f"{row['y_O2_combined']:>7.4f} {row['y_CO_combined']:>7.4f} "
                f"{row['limit_o2_bubble']:>7.3f} {row['limit_o2_dense']:>7.3f} "
                f"{row['r5b_raw_mol_s']:>9.2f} {row['r5d_raw_mol_s']:>9.2f}"
            )

    print("-" * 120)
    print("Readout:")
    if init_rows and final_rows:
        b0i, b0f = init_rows[0], final_rows[0]
        print(
            f"  bed0 primary O2 split (dense/bubble zu): init "
            f"{b0i['N_zu_d_O2']:.3f}/{b0i['N_zu_b_O2']:.3f} → final "
            f"{b0f['N_zu_d_O2']:.3f}/{b0f['N_zu_b_O2']:.3f}  "
            f"(holdup O2 dense/bubble: {b0f['N_d_O2']:.3f}/{b0f['N_b_O2']:.3f})"
        )
        if len(final_rows) > 1:
            b1f = final_rows[1]
            print(
                f"  bed1 axial O2 in (from bed0 holdup): N_d_in={b1f['N_d_in_O2']:.3f}  "
                f"N_b_in={b1f['N_b_in_O2']:.3f}  yO2={b1f['y_O2_combined']:.4f}"
            )
        print(
            f"  bed0 R5 demand vs O2 supply: bubble limit={b0f['limit_o2_bubble']:.3f} "
            f"(R5b_raw={b0f['r5b_raw_mol_s']:.1f}), dense limit={b0f['limit_o2_dense']:.3f} "
            f"(R5d_raw={b0f['r5d_raw_mol_s']:.1f})"
        )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
