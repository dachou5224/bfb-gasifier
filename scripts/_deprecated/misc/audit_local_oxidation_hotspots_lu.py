#!/usr/bin/env python3
"""LU 工况 primary / secondary inlet 局部氧化热点审计。

目的：
1. 对 `cell0`（primary inlet）输出当前单元级 O2 预算与主要快氧化项
2. 对 freeboard `secondary inlet` 段重放当前实现里的 12 个 reaction substeps
3. 判断 `O2 spike` 是否在同段内被瞬时吃掉
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from audit_oxygen_reaction_trace import _cell_o2_trace
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.reactor import Reactor, _resolve_axial_heat_loss_distribution
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, N_GAS
from src.core.freeboard_segment import _freeboard_reaction_step, _solve_T_from_enthalpy, simulate_freeboard
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    apply_secondary_air_repartition,
    build_phase2_htw_lu_freeboard_reactor_config,
)


SCENARIOS = (
    ("baseline", 0.00),
    ("sec_air_05", 0.05),
    ("sec_air_10", 0.10),
)


def _run_freeboard(cfg, reactor: Reactor):
    top = reactor.cells[-1]
    total_height = float(cfg.H_bed + max(cfg.H_freeboard, 0.0))
    _, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)
    return simulate_freeboard(
        N_in=np.maximum(top.N_d + top.N_b, 0.0),
        T_in=float(top.T),
        P=cfg.P,
        D_bed=cfg.D_bed,
        H_freeboard=float(cfg.H_freeboard),
        n_cells=int(cfg.n_freeboard_cells),
        u_b_bed_top=float(top.u_b),
        d_b_bed_top=float(top.d_b),
        eps_b_bed_top=float(top.eps_b),
        eps_d_void_bed_top=float(top.eps_d_voidage),
        rho_solid_bed_top=float(top.solid.rho_s),
        d_p_classes_bed_top=np.array(top.solid.d_p_classes, dtype=np.float64, copy=True),
        m_char_classes_bed_top=np.maximum(top.m_solid[:, 2], 0.0),
        m_ash_classes_bed_top=np.maximum(top.m_solid[:, 3], 0.0),
        fuel_type=cfg.fuel_type,
        heat_loss_frac=float(freeboard_loss),
        trajectory_model=str(cfg.freeboard_trajectory_model),
        u_gb_scale=float(max(cfg.freeboard_u_gb_scale, 1e-6)),
        beta_a_scale=float(max(cfg.freeboard_beta_a_scale, 0.0)),
        velocity_sigma=float(max(cfg.freeboard_velocity_sigma, 0.0)),
        velocity_bins=int(max(cfg.freeboard_velocity_bins, 1)),
        phi_s_bed_top=float(top.solid.phi_s),
        z_base_m=float(cfg.H_bed),
        total_height_m=total_height,
        secondary_injection_xi=cfg.freeboard_secondary_injection_xi,
        secondary_O2_mol_s=float(max(cfg.freeboard_secondary_O2_mol_s, 0.0)),
        secondary_N2_mol_s=float(max(cfg.freeboard_secondary_N2_mol_s, 0.0)),
        secondary_H2O_mol_s=float(max(cfg.freeboard_secondary_H2O_mol_s, 0.0)),
        secondary_T_K=float(max(cfg.freeboard_secondary_T_K, 1.0)),
        enabled_reactions=tuple(cfg.freeboard_enabled_reactions),
    )


def _replay_injection_segment(cfg, reactor: Reactor, fb: dict) -> list[dict[str, float]]:
    seg = fb["secondary_injection_segment"]
    if seg is None:
        return []

    states = fb["states"]
    idx = GAS_SPECIES_INDEX
    top = reactor.cells[-1]
    if seg == 0:
        N = np.maximum(top.N_d + top.N_b, 0.0).copy()
        T = float(top.T)
    else:
        N = states[seg - 1].N.copy()
        T = float(states[seg - 1].T)

    secondary_N = np.zeros(N_GAS, dtype=np.float64)
    secondary_N[idx["O2"]] = float(max(cfg.freeboard_secondary_O2_mol_s, 0.0))
    secondary_N[idx["N2"]] = float(max(cfg.freeboard_secondary_N2_mol_s, 0.0))
    secondary_N[idx["H2O"]] = float(max(cfg.freeboard_secondary_H2O_mol_s, 0.0))

    H_mix = calc_gas_enthalpy_flow(N, T, h_cache={}) + calc_gas_enthalpy_flow(
        secondary_N, float(cfg.freeboard_secondary_T_K), h_cache={}
    )
    N = N + secondary_N
    T = _solve_T_from_enthalpy(N, H_mix, max(T, float(cfg.freeboard_secondary_T_K)))

    dh = float(cfg.H_freeboard) / int(cfg.n_freeboard_cells)
    area = math.pi * (cfg.D_bed ** 2) / 4.0
    v_seg = area * dh
    tau = float(states[seg].tau)
    dt = tau / 12.0
    _, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)
    seg_loss = 1.0 - math.pow(max(1.0 - float(freeboard_loss), 0.0), 1.0 / max(int(cfg.n_freeboard_cells), 1))
    substep_loss = 1.0 - math.pow(max(1.0 - seg_loss, 0.0), 1.0 / 12.0)
    rho_cat = float(states[seg].rho_cat)

    rows: list[dict[str, float]] = []
    for k in range(12):
        h_before = calc_gas_enthalpy_flow(N, T, h_cache={})
        N, diag = _freeboard_reaction_step(
            N=N,
            T=T,
            P=cfg.P,
            V_seg=v_seg,
            dt=dt,
            fuel_type=cfg.fuel_type,
            rho_cat=rho_cat,
            enabled_reactions=tuple(cfg.freeboard_enabled_reactions),
        )
        T = _solve_T_from_enthalpy(N, h_before * (1.0 - substep_loss), T)
        n_tot = max(float(np.sum(np.maximum(N, 0.0))), 1e-12)
        rows.append(
            {
                "step": float(k + 1),
                "T": float(T),
                "O2": float(max(N[idx["O2"]], 0.0) / n_tot),
                "CO": float(max(N[idx["CO"]], 0.0) / n_tot),
                "CO2": float(max(N[idx["CO2"]], 0.0) / n_tot),
                "H2": float(max(N[idx["H2"]], 0.0) / n_tot),
                "CH4": float(max(N[idx["CH4"]], 0.0) / n_tot),
                "R5": float(diag.get("R5", 0.0)),
                "R6": float(diag.get("R6", 0.0)),
                "R12": float(diag.get("R12", 0.0)),
            }
        )
    return rows


def main() -> int:
    print("=" * 196)
    print("LU local oxidation hotspot audit")
    print("=" * 196)
    for name, sec_frac in SCENARIOS:
        cfg = build_phase2_htw_lu_freeboard_reactor_config()
        if sec_frac > 0.0:
            apply_secondary_air_repartition(
                cfg,
                secondary_air_frac=sec_frac,
                injection_xi=cfg.freeboard_secondary_injection_xi,
                secondary_T_K=cfg.T_inlet,
            )
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
            check_x0=True,
        )
        fb = _run_freeboard(cfg, reactor)

        print("-" * 196)
        print(
            f"{name}: sec_air_frac={sec_frac:.2f} Texit={result['reactor_exit_T']:.1f}K "
            f"inj_seg={fb['secondary_injection_segment']} "
            f"inj_xi={result['freeboard_axial_xi'][fb['secondary_injection_segment']] if fb['secondary_injection_segment'] is not None else float('nan'):.3f}"
        )
        print_nr_monitor(monitor, prefix=f"NR monitor [{name}]")
        print("-" * 196)

        cell0 = reactor.cells[0]
        cell0.calc_hydrodynamics()
        cell0.compute_vorabrechnung(cell0.geo.dh / max(cell0.u_mf, 1e-3))
        cell0.calc_exchange()
        row0 = _cell_o2_trace(cell0)
        print("cell0 / primary inlet:")
        print(
            f"  T={row0['T']:.1f}K O2sup_b={row0['O2_supply_bubble']:.3f} O2sup_d={row0['O2_supply_dense']:.3f} "
            f"O2dem_b={row0['bubble_o2_unlimited']:.3f} O2dem_d={row0['dense_o2_unlimited']:.3f} "
            f"lim_b={row0['limit_factor_o2_bubble']:.3f} lim_d={row0['limit_factor_o2_dense']:.3f} "
            f"short_b={row0['O2_shortfall_bubble']:.3f} short_d={row0['O2_shortfall_dense']:.3f}"
        )
        print(
            f"  fast oxidation: R5b={row0['r5b_o2']:.3f} R6b={row0['r6b_o2']:.3f} R12b={row0['r12b_o2']:.3f} "
            f"R5d={row0['r5d_o2']:.3f} R6d={row0['r6d_o2']:.3f} R12d={row0['r12d_o2']:.3f} "
            f"R10={row0['r10_o2']:.3f} R1={row0['r1_char_comb_o2']:.3f}"
        )

        steps = _replay_injection_segment(cfg, reactor, fb)
        if not steps:
            print("freeboard / secondary inlet: not active")
            continue

        print("freeboard secondary inlet / per-substep replay:")
        print(
            f"{'step':>4} {'T[K]':>8} {'O2':>8} {'CO':>8} {'CO2':>8} {'H2':>8} {'CH4':>8} {'R5':>9} {'R6':>9} {'R12':>9}"
        )
        for r in steps:
            print(
                f"{int(r['step']):>4d} {r['T']:>8.1f} {r['O2']:>8.4f} {r['CO']:>8.4f} {r['CO2']:>8.4f} "
                f"{r['H2']:>8.4f} {r['CH4']:>8.4f} {r['R5']:>9.4f} {r['R6']:>9.4f} {r['R12']:>9.4f}"
            )
        print(
            "  note: 若 step1-step2 内 O2 就掉到接近 0，则当前 coarse freeboard 口径会天然洗平局部 O2 spike。"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
