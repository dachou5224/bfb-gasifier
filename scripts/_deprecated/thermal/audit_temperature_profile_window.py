"""Focused tuning-window audit around stable temperature branch.

Window:
- refined mesh (bottom-refined)
- gas_inlet_dense_frac around 0.35
- heat_loss fixed at 0.10
- small kinetic perturbations on R5 / R6 / R4

Goal:
- keep temperature profile smooth and bounded
- push O2 exhaustion earlier
- improve CO / CH4 direction without exploding temperature

Usage:
    .venv/bin/python scripts/audit_temperature_profile_window.py
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import src.kinetics.char_reactions as cr
import src.kinetics.gas_reactions as gr
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as IDX
from tests.validation_case_utils import (
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_refined_config,
    load_case_LU,
    load_validation_case_node,
)


@dataclass
class Row:
    dense: float
    r5s: float
    r6s: float
    r4s: float
    T_exit: float
    T_peak: float
    i_peak: int
    smooth_penalty: float
    o2_exhaust_cell: int | None
    Xc: float
    CO: float
    CO2: float
    H2: float
    CH4: float
    score: float


@contextmanager
def kinetics_scaling(r5_scale: float, r6_scale: float, r4_scale: float):
    b_r5 = gr.R5_suspension_k_T05
    b_r6 = gr.R6_k0
    b_r4 = cr.R4_kf_k0
    try:
        gr.R5_suspension_k_T05 = b_r5 * r5_scale
        gr.R6_k0 = b_r6 * r6_scale
        cr.R4_kf_k0 = b_r4 * r4_scale
        yield
    finally:
        gr.R5_suspension_k_T05 = b_r5
        gr.R6_k0 = b_r6
        cr.R4_kf_k0 = b_r4


def _o2_exhaust_cell(reactor: Reactor, cfg) -> int | None:
    o2_ref = max(float(cfg.O2_feed), 1e-12)
    for i, c in enumerate(reactor.cells):
        n = np.maximum(c.N_d + c.N_b, 0.0)
        if float(n[IDX["O2"]]) <= 0.01 * o2_ref:
            return i
    return None


def _smoothness_penalty(T: np.ndarray) -> float:
    # Penalize strong oscillations / spikes using second finite diff norm
    if len(T) < 3:
        return 0.0
    d2 = np.diff(T, n=2)
    return float(np.sqrt(np.mean(d2**2)))


def rel_err(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-12)


def run_one(dense: float, r5s: float, r6s: float, r4s: float) -> Row:
    with kinetics_scaling(r5s, r6s, r4s):
        case = load_case_LU()
        cfg = build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)
        cfg.gas_inlet_dense_frac = dense
        cfg.heat_loss_frac = 0.10

        r = Reactor(cfg)
        res = r.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    T = np.array(res["T_profile"], dtype=float)
    i_peak = int(np.argmax(T))
    T_peak = float(T[i_peak])
    T_exit = float(T[-1])
    smooth = _smoothness_penalty(T)

    y = res["exit_gas_dry"]
    xc = float(res["carbon_conv"])
    co = float(y.get("CO", 0.0))
    co2 = float(y.get("CO2", 0.0))
    h2 = float(y.get("H2", 0.0))
    ch4 = float(y.get("CH4", 0.0))
    o2_cell = _o2_exhaust_cell(r, cfg)

    # Reference targets (validation JSON)
    ref = load_validation_case_node()["outputs"]
    y_ref = ref.get("exit_gas_dry_mol_frac", {})
    T_ref = float(ref.get("exit_temperature_K", 1100.0))

    # Composite score: prioritize stable temperature profile first
    score = 0.0
    score += 2.0 * rel_err(T_exit, T_ref)
    score += 0.003 * max(T_peak - 1350.0, 0.0)  # strong penalty for hot spikes
    score += 0.01 * smooth

    # Secondary composition guidance
    score += 1.0 * rel_err(co, float(y_ref.get("CO", 0.13)))
    score += 0.7 * rel_err(co2, float(y_ref.get("CO2", 0.11)))
    score += 0.7 * rel_err(h2, float(y_ref.get("H2", 0.12)))
    score += 0.5 * rel_err(ch4, float(y_ref.get("CH4", 0.028)))

    # Reward earlier O2 exhaustion, but do not dominate score
    if o2_cell is None:
        score += 0.4
    else:
        score += 0.03 * max(o2_cell - 5, 0)

    # Reward conversion progress modestly
    score += 0.8 * rel_err(xc, 0.95)

    return Row(
        dense=dense,
        r5s=r5s,
        r6s=r6s,
        r4s=r4s,
        T_exit=T_exit,
        T_peak=T_peak,
        i_peak=i_peak,
        smooth_penalty=smooth,
        o2_exhaust_cell=o2_cell,
        Xc=xc,
        CO=co,
        CO2=co2,
        H2=h2,
        CH4=ch4,
        score=score,
    )


def main() -> int:
    dense_grid = [0.30, 0.35, 0.40]
    r5_grid = [0.5, 0.7, 1.0]
    r6_grid = [0.5, 0.7, 1.0]
    r4_grid = [1.0, 1.5, 2.0]

    rows: list[Row] = []
    n_total = len(dense_grid) * len(r5_grid) * len(r6_grid) * len(r4_grid)
    n = 0

    print("=" * 140)
    print("Temperature-profile focused tuning window")
    print("=" * 140)

    for dense in dense_grid:
        for r5s in r5_grid:
            for r6s in r6_grid:
                for r4s in r4_grid:
                    n += 1
                    row = run_one(dense, r5s, r6s, r4s)
                    rows.append(row)
                    o2 = "None" if row.o2_exhaust_cell is None else str(row.o2_exhaust_cell)
                    print(
                        f"[{n:03d}/{n_total}] dense={dense:.2f} r5={r5s:.2f} r6={r6s:.2f} r4={r4s:.2f} | "
                        f"Texit={row.T_exit:7.1f} Tpk={row.T_peak:7.1f}@{row.i_peak:02d} smooth={row.smooth_penalty:6.1f} "
                        f"Xc={row.Xc:5.3f} CO={row.CO:5.3f} CO2={row.CO2:5.3f} H2={row.H2:5.3f} CH4={row.CH4:5.3f} O2exh={o2:>4} "
                        f"score={row.score:6.3f}"
                    )

    rows.sort(key=lambda r: r.score)

    print("\n" + "=" * 140)
    print("Top 12 candidates (lower score is better)")
    print("=" * 140)
    print(
        f"{'rank':>4} {'dense':>6} {'r5':>5} {'r6':>5} {'r4':>5} {'Texit':>8} {'Tpeak':>8} {'i_pk':>5} {'smooth':>8} "
        f"{'Xc':>6} {'CO':>6} {'CO2':>6} {'H2':>6} {'CH4':>6} {'O2exh':>6} {'score':>8}"
    )
    for i, r in enumerate(rows[:12], start=1):
        o2 = "None" if r.o2_exhaust_cell is None else str(r.o2_exhaust_cell)
        print(
            f"{i:4d} {r.dense:6.2f} {r.r5s:5.2f} {r.r6s:5.2f} {r.r4s:5.2f} {r.T_exit:8.1f} {r.T_peak:8.1f} {r.i_peak:5d} {r.smooth_penalty:8.1f} "
            f"{r.Xc:6.3f} {r.CO:6.3f} {r.CO2:6.3f} {r.H2:6.3f} {r.CH4:6.3f} {o2:>6} {r.score:8.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
