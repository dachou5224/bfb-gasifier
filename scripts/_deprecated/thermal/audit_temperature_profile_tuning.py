"""Temperature-profile investigation under tuning scenarios (HTW LU, Phase-1).

Purpose
-------
Run a small set of physically motivated tuning scenarios and compare:
- axial temperature profile (per-cell T)
- profile shape metrics (peak location, bottom-zone gradient)
- key KPI side effects (carbon conversion, dry gas)
- O2 exhaustion location

Usage
-----
    .venv/bin/python scripts/audit_temperature_profile_tuning.py
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
    build_phase1_htw_lu_reactor_config,
    build_phase1_htw_lu_refined_config,
    load_case_LU,
)


@dataclass
class Scenario:
    name: str
    refined: bool = True
    n_fine: int = 3
    dh_fine: float = 0.15
    gas_inlet_dense_frac: float | None = None
    heat_loss_frac: float | None = None
    r5_scale: float = 1.0
    r6_scale: float = 1.0
    r4_scale: float = 1.0


@contextmanager
def kinetics_scaling(r5_scale: float = 1.0, r6_scale: float = 1.0, r4_scale: float = 1.0):
    """Temporary global scaling of selected kinetic constants."""
    base_r5 = gr.R5_suspension_k_T05
    base_r6 = gr.R6_k0
    base_r4 = cr.R4_kf_k0
    try:
        gr.R5_suspension_k_T05 = base_r5 * float(r5_scale)
        gr.R6_k0 = base_r6 * float(r6_scale)
        cr.R4_kf_k0 = base_r4 * float(r4_scale)
        yield
    finally:
        gr.R5_suspension_k_T05 = base_r5
        gr.R6_k0 = base_r6
        cr.R4_kf_k0 = base_r4


def _build_cfg(sc: Scenario):
    case = load_case_LU()
    if sc.refined:
        cfg = build_phase1_htw_lu_refined_config(case, n_fine=sc.n_fine, dh_fine=sc.dh_fine)
    else:
        cfg = build_phase1_htw_lu_reactor_config(case)

    if sc.gas_inlet_dense_frac is not None:
        cfg.gas_inlet_dense_frac = float(sc.gas_inlet_dense_frac)
    if sc.heat_loss_frac is not None:
        cfg.heat_loss_frac = float(sc.heat_loss_frac)
    return cfg


def _o2_exhaust_cell(reactor: Reactor, cfg) -> int | None:
    o2_ref = max(float(cfg.O2_feed), 1e-12)
    for i, c in enumerate(reactor.cells):
        N = np.maximum(c.N_d + c.N_b, 0.0)
        if float(N[IDX["O2"]]) <= 0.01 * o2_ref:
            return i
    return None


def run_scenario(sc: Scenario) -> dict:
    with kinetics_scaling(sc.r5_scale, sc.r6_scale, sc.r4_scale):
        cfg = _build_cfg(sc)
        reactor = Reactor(cfg)
        res = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    T = np.array(res["T_profile"], dtype=float)
    i_peak = int(np.argmax(T))
    o2_exhaust = _o2_exhaust_cell(reactor, cfg)
    y_dry = res["exit_gas_dry"]

    return {
        "name": sc.name,
        "n_cells": cfg.n_cells,
        "dh_head": [float(c.geo.dh) for c in reactor.cells[: min(4, len(reactor.cells))]],
        "T_profile": T.tolist(),
        "T_in": float(T[0]),
        "T_peak": float(T[i_peak]),
        "i_peak": i_peak,
        "T_exit": float(T[-1]),
        "dT_01": float(T[1] - T[0]) if len(T) > 1 else 0.0,
        "dT_12": float(T[2] - T[1]) if len(T) > 2 else 0.0,
        "carbon_conv": float(res["carbon_conv"]),
        "CO": float(y_dry.get("CO", 0.0)),
        "CO2": float(y_dry.get("CO2", 0.0)),
        "H2": float(y_dry.get("H2", 0.0)),
        "CH4": float(y_dry.get("CH4", 0.0)),
        "o2_exhaust_cell": o2_exhaust,
        "converged": bool(res.get("converged", False)),
        "residual": float(res.get("residual", np.nan)),
    }


def main() -> int:
    scenarios = [
        Scenario(name="A_baseline_uniform", refined=False),
        Scenario(name="B_refined_default", refined=True),
        Scenario(name="C_refined_hl0.12", refined=True, heat_loss_frac=0.12),
        Scenario(name="D_refined_dense0.35", refined=True, gas_inlet_dense_frac=0.35),
        Scenario(name="E_refined_dense0.35_hl0.12", refined=True, gas_inlet_dense_frac=0.35, heat_loss_frac=0.12),
        Scenario(name="F_refined_dense0.35_hl0.12_r5x0.3_r4x3", refined=True, gas_inlet_dense_frac=0.35, heat_loss_frac=0.12, r5_scale=0.3, r4_scale=3.0),
        Scenario(name="G_refined_dense0.35_hl0.12_r5x0.3_r6x0.5_r4x3", refined=True, gas_inlet_dense_frac=0.35, heat_loss_frac=0.12, r5_scale=0.3, r6_scale=0.5, r4_scale=3.0),
    ]

    print("=" * 132)
    print("Temperature profile tuning audit (HTW LU)")
    print("=" * 132)
    print(
        f"{'scenario':<46} {'n':>3} {'T0':>7} {'Tpk':>7} {'i_pk':>5} {'Texit':>7} "
        f"{'dT01':>7} {'dT12':>7} {'Xc':>6} {'CO':>6} {'CO2':>6} {'H2':>6} {'CH4':>6} {'O2exh':>6}"
    )
    print("-" * 132)

    rows: list[dict] = []
    for sc in scenarios:
        row = run_scenario(sc)
        rows.append(row)
        o2s = "None" if row["o2_exhaust_cell"] is None else str(row["o2_exhaust_cell"])
        print(
            f"{row['name']:<46} {row['n_cells']:>3d} {row['T_in']:>7.1f} {row['T_peak']:>7.1f} "
            f"{row['i_peak']:>5d} {row['T_exit']:>7.1f} {row['dT_01']:>7.1f} {row['dT_12']:>7.1f} "
            f"{row['carbon_conv']:>6.3f} {row['CO']:>6.3f} {row['CO2']:>6.3f} {row['H2']:>6.3f} {row['CH4']:>6.3f} {o2s:>6}"
        )

    print("-" * 132)
    print("Detailed T profiles:")
    for r in rows:
        prof = ", ".join(f"{t:.0f}" for t in r["T_profile"])
        print(f"  {r['name']}: [{prof}]")

    print("=" * 132)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
