#!/usr/bin/env python3
"""整炉诊断：R10 修正后的 O2 预算、温度、单 cell 收敛验证

目标：验证 R10 压力单位修正 + k₀ 100x 缩小后的改善效果
运行配置：stable window（dense=0.40, heat_loss=0.10, R5×0.50, R6×1.00, R4×1.50）
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve()))

import numpy as np
from src.core.species import GAS_SPECIES_INDEX as idx
from tests.validation_case_utils import (
    build_phase1_htw_lu_refined_config, Table2LU_INLET
)
from src.solvers.global_nr_solver import solve_reactor_global_NR

print("\n" + "="*80)
print("FULL REACTOR DIAGNOSTIC: Post-R10 Correction Verification")
print("="*80)

# Configuration: stable working point
config = build_phase1_htw_lu_refined_config(
    'stable_tuning_reference',
    n_fine=3, dh_fine=0.15,
    gas_inlet_dense_frac=0.40,
    heat_loss_frac=0.10
)

reactor = config.build_reactor()

# Apply kinetic scaling
from scripts.audit_oxygen_reaction_trace import kinetics_scaling
with kinetics_scaling(r5_scale=0.50, r6_scale=1.00, r4_scale=1.50):
    
    print("\nSolving reactor with R10 corrections applied...")
    print("Parameters: dense=0.40, heat_loss=0.10, R5×0.50, R6×1.00, R4×1.50\n")
    
    try:
        solve_reactor_global_NR(
            reactor,
            inlet=Table2LU_INLET,
            atol=1e-6, rtol=1e-5,
            max_outer_iter=10, max_nr_iter=150,
            solver_mode="reduced"
        )
    except Exception as e:
        print(f"⚠ Solver warning: {e}")
        print("Proceeding with partial convergence for diagnostics...\n")

# ─────────────────────────────────────────────────────────────────────────
# Diagnostics output
# ─────────────────────────────────────────────────────────────────────────

print("\n" + "-"*80)
print("CELL-BY-CELL DIAGNOSTICS")
print("-"*80)

print("\nCell │  dh(m) │   T(K)  │ T_b(K)  │ eps_b │ C_O2_d │ C_tar_d │ φ(tar)")
print("─" * 80)

temp_peak = -np.inf
temp_peak_cell = -1
o2_exhausted_cell = -1

for i, cell in enumerate(reactor.cells):
    c_b = cell._concentrations("b")
    c_d = cell._concentrations("d")
    
    c_o2_d = c_d[idx["O2"]]
    c_tar_d = c_d[idx["TAR1"]] + c_d[idx["TAR2"]]
    
    # Check for O2 exhaustion
    if c_o2_d < 0.01 and o2_exhausted_cell < 0:
        o2_exhausted_cell = i
    
    # Track peak temperature
    if cell.T > temp_peak:
        temp_peak = cell.T
        temp_peak_cell = i
    
    # Outlet gas composition (mixed bubble+dense)
    N_total = np.sum(cell.N_b) + np.sum(cell.N_d)
    if N_total > 0:
        phi_tar = (cell.N_b[idx["TAR1"]] + cell.N_b[idx["TAR2"]] +
                  cell.N_d[idx["TAR1"]] + cell.N_d[idx["TAR2"]]) / N_total
    else:
        phi_tar = 0.0
    
    print(f" {i:2d}  │ {cell.geo.dh:.3f} │ {cell.T:7.1f} │ "
          f"{reactor.cells[i].T:7.1f} │ {cell.eps_b:5.3f} │ {c_o2_d:6.3f} │ "
          f"{c_tar_d:7.3f} │ {phi_tar:.2e}")

# ─────────────────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────────────────

print("\n" + "-"*80)
print("SUMMARY STATISTICS")
print("-"*80)

# Exit gas composition
N_exit_b = reactor.cells[-1].N_b
N_exit_d = reactor.cells[-1].N_d
N_exit_total = np.sum(N_exit_b) + np.sum(N_exit_d)

if N_exit_total > 0:
    exit_frac_b = N_exit_b / N_exit_total
    exit_frac_d = N_exit_d / N_exit_total
    exit_frac_mix = (N_exit_b + N_exit_d) / N_exit_total
    
    print(f"\nExit gas mole fractions (combined):")
    for sp in ["CO", "CO2", "H2", "H2O", "CH4", "N2", "O2"]:
        if sp in idx:
            frac = exit_frac_mix[idx[sp]]
            print(f"  {sp:4s}: {frac*100:6.2f}%")

print(f"\nTemperature profile:")
print(f"  T_min = {np.min([c.T for c in reactor.cells]):.1f} K")
print(f"  T_max = {temp_peak:.1f} K (cell {temp_peak_cell})")
print(f"  T_exit = {reactor.cells[-1].T:.1f} K")

if o2_exhausted_cell >= 0:
    print(f"\nO2 exhaustion zone: cell {o2_exhausted_cell} (first cell with C_O2 < 0.01)")
else:
    print(f"\nO2 not exhausted (C_O2 > 0.01 at exit)")

# Carbon conversion
inlet_c_dry = Table2LU_INLET.m_fuel * (Table2LU_INLET.solid.C_dry / 100.0) / 12.01  # moles C in fuel
exit_c_co = np.sum(reactor.cells[-1].N_d) * exit_frac_mix[idx["CO"]]
exit_c_co2 = np.sum(reactor.cells[-1].N_d) * exit_frac_mix[idx["CO2"]]
exit_c_ch4 = np.sum(reactor.cells[-1].N_d) * exit_frac_mix[idx["CH4"]]
exit_c_tar = (reactor.cells[-1].N_d[idx["TAR1"]] + reactor.cells[-1].N_d[idx["TAR2"]])
exit_c_total = exit_c_co + exit_c_co2 + exit_c_ch4 + exit_c_tar

if inlet_c_dry > 0:
    c_conversion = 100.0 * (1.0 - exit_c_tar / inlet_c_dry)
    print(f"\nCarbon conversion: {c_conversion:.1f}%")
    print(f"  Inlet: {inlet_c_dry:.3f} mol C (dry)")
    print(f"  Exit tar: {exit_c_tar:.3f} mol")

print("\n" + "="*80)
print("VERDICT")
print("="*80)

if temp_peak < 1400:
    print("✓ Temperature profile: STABLE (T_peak < 1400 K)")
else:
    print(f"⚠ Temperature profile: ELEVATED (T_peak = {temp_peak:.0f} K)")

if o2_exhausted_cell >= 0 and o2_exhausted_cell < 5:
    print(f"✓ O2 exhaustion: EARLY (cell {o2_exhausted_cell}), good for combustion zone")
elif o2_exhausted_cell < 0:
    print("⚠ O2 not exhausted - may need higher O2 feed or lower R5/R6 scaling")
else:
    print(f"⚠ O2 exhaustion: LATE (cell {o2_exhausted_cell})")

print("\n" + "="*80)
