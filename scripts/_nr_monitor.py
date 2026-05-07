#!/usr/bin/env python3
"""Shared NR monitor helpers: wall-time / RMS / Vorabrechnung x0 sanity."""

from __future__ import annotations

from time import perf_counter
from typing import Any

import numpy as np

from src.core.reactor import Reactor
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0


def build_vorab_x0_sanity(reactor: Reactor) -> dict[str, Any]:
    """Build Hamel-style Vorabrechnung x0 and report physical sanity metrics."""
    cfg = reactor.config
    T_est = estimate_axial_T_profile(
        n_cells=int(cfg.n_cells),
        T_inlet=float(cfg.T_inlet),
        O2_feed=float(cfg.O2_feed),
        H2O_feed=float(cfg.H2O_feed),
        N2_feed=float(cfg.N2_feed),
        fuel_feed_kg_s=float(cfg.fuel_feed),
        C_dry=float(cfg.C_dry),
        H_dry=float(cfg.H_dry),
        moisture_wt=float(cfg.moisture_wt),
        P=float(cfg.P),
    )
    for i, cell in enumerate(reactor.cells):
        cell.T = float(T_est[i])
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    for cell in reactor.cells:
        cell.calc_hydrodynamics()
    generate_initial_x0(
        cells=reactor.cells,
        O2_feed=float(cfg.O2_feed),
        H2O_feed=float(cfg.H2O_feed),
        N2_feed=float(cfg.N2_feed),
        fuel_feed_kg_s=float(cfg.fuel_feed),
        C_dry=float(cfg.C_dry),
        H_dry=float(cfg.H_dry),
        O_dry=float(cfg.O_dry),
        moisture_wt=float(cfg.moisture_wt),
        ash_dry_wt=float(cfg.ash_dry_wt),
        VM_daf=float(cfg.VM_daf),
        T_profile=T_est,
        fuel_type=cfg.fuel_type,
        use_hamel_major_gibbs_x0=bool(cfg.vorab_major_gibbs_x0),
        strict_hamel_major_gibbs_x0=bool(cfg.thesis_mode),
        major_gibbs_solver_mode=str(cfg.major_gibbs_solver_mode),
    )

    # Re-evaluate hydrodynamics on generated x0 to catch physically non-credible seeds
    # that may still be finite/non-negative but destabilize inner NR.
    for cell in reactor.cells:
        cell.calc_hydrodynamics()

    all_Nd = np.concatenate([c.N_d for c in reactor.cells]) if reactor.cells else np.zeros(0)
    all_Nb = np.concatenate([c.N_b for c in reactor.cells]) if reactor.cells else np.zeros(0)
    all_ms = np.concatenate([c.m_solid.reshape(-1) for c in reactor.cells]) if reactor.cells else np.zeros(0)
    eps_b = np.array([float(c.eps_b) for c in reactor.cells], dtype=float) if reactor.cells else np.zeros(0)
    u_b = np.array([float(c.u_b) for c in reactor.cells], dtype=float) if reactor.cells else np.zeros(0)
    kbd = np.array([float(c.K_bd) for c in reactor.cells], dtype=float) if reactor.cells else np.zeros(0)
    cell_total_mol = np.array([float(np.sum(c.N_d + c.N_b)) for c in reactor.cells], dtype=float) if reactor.cells else np.zeros(0)
    inlet_total_mol = float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed)
    max_cell_total_mol = float(np.max(cell_total_mol)) if cell_total_mol.size else 0.0
    max_cell_total_to_inlet = max_cell_total_mol / max(inlet_total_mol, 1e-12)

    finite_ok = bool(
        np.all(np.isfinite(all_Nd))
        and np.all(np.isfinite(all_Nb))
        and np.all(np.isfinite(all_ms))
        and np.all(np.isfinite(T_est))
    )
    nonneg_ok = bool(np.min(all_Nd, initial=0.0) >= -1e-12 and np.min(all_Nb, initial=0.0) >= -1e-12 and np.min(all_ms, initial=0.0) >= -1e-12)
    T_ok = bool(np.min(T_est, initial=300.0) >= 250.0 and np.max(T_est, initial=1200.0) <= 2600.0)
    epsb_ok = bool(np.min(eps_b, initial=0.0) >= -1e-8 and np.max(eps_b, initial=0.0) <= 0.98)
    ub_ok = bool(np.min(u_b, initial=0.0) >= -1e-8 and np.max(u_b, initial=0.0) <= 20.0)
    kbd_ok = bool(np.min(kbd, initial=0.0) >= -1e-10)
    molar_scale_ok = bool(max_cell_total_to_inlet <= 100.0)

    t_min = float(np.min(T_est)) if T_est.size else float("nan")
    t_max = float(np.max(T_est)) if T_est.size else float("nan")
    eps_min = float(np.min(eps_b)) if eps_b.size else float("nan")
    eps_max = float(np.max(eps_b)) if eps_b.size else float("nan")
    ub_min = float(np.min(u_b)) if u_b.size else float("nan")
    ub_max = float(np.max(u_b)) if u_b.size else float("nan")
    kbd_min = float(np.min(kbd)) if kbd.size else float("nan")
    kbd_max = float(np.max(kbd)) if kbd.size else float("nan")

    return {
        "ok": bool(finite_ok and nonneg_ok and T_ok and epsb_ok and ub_ok and kbd_ok and molar_scale_ok),
        "finite_ok": finite_ok,
        "nonneg_ok": nonneg_ok,
        "T_ok": T_ok,
        "epsb_ok": epsb_ok,
        "ub_ok": ub_ok,
        "kbd_ok": kbd_ok,
        "molar_scale_ok": molar_scale_ok,
        "inlet_total_mol_s": inlet_total_mol,
        "max_cell_total_mol_s": max_cell_total_mol,
        "max_cell_total_to_inlet": max_cell_total_to_inlet,
        "T_min": t_min,
        "T_max": t_max,
        "eps_b_min": eps_min,
        "eps_b_max": eps_max,
        "u_b_min": ub_min,
        "u_b_max": ub_max,
        "K_bd_min": kbd_min,
        "K_bd_max": kbd_max,
    }


def solve_with_nr_monitor(
    reactor: Reactor,
    solve_kwargs: dict[str, Any],
    *,
    check_x0: bool = True,
) -> tuple[dict[str, Any], dict[str, Any]]:
    x0 = build_vorab_x0_sanity(reactor) if check_x0 else None
    t0 = perf_counter()
    result = reactor.solve(**solve_kwargs)
    wall_s = perf_counter() - t0
    monitor = {
        "wall_time_s": float(wall_s),
        "rms_scaled_final": float(result.get("rms_scaled_final", np.nan)),
        "converged": bool(result.get("converged", False)),
        "converged_outer": result.get("converged_outer"),
        "converged_inner_nr": result.get("converged_inner_nr"),
        "n_iter": int(result.get("n_iter", 0)),
        "nr_total_s": float(result.get("nr_total_s", np.nan)) if result.get("nr_total_s") is not None else np.nan,
        "x0_sanity": x0,
    }
    return result, monitor


def print_nr_monitor(m: dict[str, Any], *, prefix: str = "NR monitor") -> None:
    print(
        f"{prefix}: wall={m['wall_time_s']:.2f}s nr_total={m['nr_total_s']:.2f}s "
        f"rms={m['rms_scaled_final']:.3e} n_iter={m['n_iter']} "
        f"conv={m['converged']} outer={m['converged_outer']} inner={m['converged_inner_nr']}"
    )
    x0 = m.get("x0_sanity")
    if isinstance(x0, dict):
        print(
            "x0 sanity: "
            f"ok={x0.get('ok')} finite={x0.get('finite_ok')} nonneg={x0.get('nonneg_ok')} "
            f"molar_scale={x0.get('molar_scale_ok')} ratio={x0.get('max_cell_total_to_inlet'):.3e} "
            f"T=[{x0.get('T_min'):.1f},{x0.get('T_max'):.1f}] "
            f"eps_b=[{x0.get('eps_b_min'):.4f},{x0.get('eps_b_max'):.4f}] "
            f"u_b=[{x0.get('u_b_min'):.4f},{x0.get('u_b_max'):.4f}] "
            f"Kbd=[{x0.get('K_bd_min'):.4f},{x0.get('K_bd_max'):.4f}]"
        )
