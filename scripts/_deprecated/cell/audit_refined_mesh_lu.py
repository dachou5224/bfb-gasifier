"""Audit: compare baseline (uniform dh=0.6 m, n_cells=10) vs bottom-refined mesh
(3 fine cells dh=0.15 m + 10 upper cells dh=0.555 m) for the HTW LU case.

Reports per-cell rms_scaled, O2 profile, and top-level KPIs.

Usage::

    .venv/bin/python scripts/audit_refined_mesh_lu.py

"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.solvers.cell_solver import solve_cell
from tests.validation_case_utils import (
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_reactor_config,
    build_phase1_htw_lu_refined_config,
    load_case_LU,
    load_validation_case_node,
)

# ── helpers ─────────────────────────────────────────────────────────────────

def _cell_trace(reactor: Reactor) -> list[dict]:
    """Run ONE fresh GS sweep and collect per-cell diagnostics."""
    cfg = reactor.config
    reactor._set_bottom_cell_feeds()
    # Initialise recirculation stores
    from src.core.cell import N_SOLID_COMP, S_CHAR, S_ASH
    from src.core.species import N_GAS
    nk = max(1, cfg.n_age_classes)
    _m_solid_rez = np.zeros((nk, N_SOLID_COMP))
    _n_rez_d = np.zeros(N_GAS)
    _n_rez_b = np.zeros(N_GAS)

    trace = []
    for i, cell in enumerate(reactor.cells):
        reactor._propagate_upstream(i)
        total_gas = float(np.sum(np.maximum(cell.N_d + cell.N_b, 0.0)))
        if total_gas < 1e-9:
            cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + cell.N_rez_d, 0.0)
            cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + cell.N_rez_b, 0.0)
            cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez, 0.0)
        cell.calc_hydrodynamics()
        cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
        res = solve_cell(cell, verbose=False)
        y_d = cell._mole_fractions("d")
        y_b = cell._mole_fractions("b")
        trace.append({
            "cell":        i,
            "dh":          cell.geo.dh,
            "h_center":    cell.geo.h_center,
            "T":           cell.T,
            "rms_scaled":  res["rms_scaled"],
            "converged":   res["converged"],
            "O2_d":        float(y_d[idx["O2"]]),
            "O2_b":        float(y_b[idx["O2"]]),
            "CO_d":        float(y_d[idx["CO"]]),
            "CO2_d":       float(y_d[idx["CO2"]]),
            "H2_d":        float(y_d[idx["H2"]]),
        })
    return trace


def _print_trace(label: str, trace: list[dict]) -> None:
    print(f"\n{'─'*70}")
    print(f"  {label}")
    print(f"{'─'*70}")
    hdr = f"{'cell':>4}  {'dh[m]':>6}  {'h[m]':>6}  {'T[K]':>7}  {'rms_sc':>7}  "
    hdr += f"{'O2_d':>7}  {'O2_b':>7}  {'CO_d':>7}  {'CO2_d':>7}  {'H2_d':>7}"
    print(hdr)
    print("─" * 70)
    for r in trace:
        print(
            f"  {r['cell']:>2}  {r['dh']:>6.3f}  {r['h_center']:>6.3f}  {r['T']:>7.1f}  "
            f"{r['rms_scaled']:>7.3f}  "
            f"{r['O2_d']:>7.4f}  {r['O2_b']:>7.4f}  {r['CO_d']:>7.4f}  "
            f"{r['CO2_d']:>7.4f}  {r['H2_d']:>7.4f}"
        )


def _run_full_solve(reactor: Reactor) -> dict:
    kwargs = dict(PHASE1_HTW_LU_SOLVE_KWARGS)
    return reactor.solve(**kwargs)


def _print_kpis(label: str, result: dict, case_node: dict) -> None:
    outputs = case_node.get("outputs", {})
    T_exp = float(outputs.get("T_exit_K", 1120.0))
    cc_exp = float(outputs.get("carbon_conversion", 0.95))
    gas_exp = outputs.get("exit_gas_dry_mole_frac", {})

    T_out = result["T_profile"][-1]
    cc_out = result["carbon_conv"]
    dry = result["exit_gas_dry"]

    print(f"\n  KPIs [{label}]")
    print(f"    T_exit  : {T_out:7.1f} K   target={T_exp:.0f} K   err={abs(T_out-T_exp)/T_exp*100:.1f}%")
    print(f"    carbon_conv: {cc_out:.3f}   target={cc_exp:.3f}   err={abs(cc_out-cc_exp)/cc_exp*100:.1f}%")
    for sp in ["CO", "CO2", "H2", "CH4"]:
        y_out = dry.get(sp, 0.0)
        y_exp = float(gas_exp.get(sp, 0.0))
        if y_exp > 1e-6:
            err = abs(y_out - y_exp) / y_exp * 100
            print(f"    {sp:4s}   : {y_out:.4f}  target={y_exp:.4f}   err={err:.1f}%")
        else:
            print(f"    {sp:4s}   : {y_out:.4f}  (no target)")


# ── main ────────────────────────────────────────────────────────────────────

def main() -> None:
    case = load_case_LU()
    case_node = load_validation_case_node()

    # ── Baseline ──
    cfg_base = build_phase1_htw_lu_reactor_config(case)
    r_base = Reactor(cfg_base)
    print("Running baseline full solve …", flush=True)
    res_base = _run_full_solve(r_base)

    # ── Refined mesh ──
    cfg_ref = build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)
    r_ref = Reactor(cfg_ref)
    print("Running refined-mesh full solve …", flush=True)
    res_ref = _run_full_solve(r_ref)

    # Per-cell diagnostics (one fresh sweep from best-state)
    trace_base = _cell_trace(r_base)
    trace_ref  = _cell_trace(r_ref)

    _print_trace("BASELINE  (n=10, dh=0.60 m uniform)", trace_base)
    _print_trace(f"REFINED   (n={cfg_ref.n_cells}, 3×dh=0.15 m + {cfg_ref.n_cells-3}×dh≈0.555 m)", trace_ref)

    print(f"\n{'═'*70}")
    _print_kpis("BASELINE", res_base, case_node)
    _print_kpis("REFINED", res_ref, case_node)
    print()


if __name__ == "__main__":
    main()
