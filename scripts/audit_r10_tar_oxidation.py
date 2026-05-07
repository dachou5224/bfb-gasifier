#!/usr/bin/env python3
"""R10 焦油氧化反应 —— 压力尺度、速率量级审计

目标：检查 R10 在单 cell 中是否存在物理不合理的量级问题。

关键问题：
1. P^0.3 在 25 bar (2.5e6 Pa) 的贡献有多大？
2. 芳香 vs 烷烃的 k0 相差大约多少倍？
3. 在实际 tar 浓度（~0.2 mol/m³）和 O2 浓度（~0.5–2 mol/m³）下，R10 速率应该是多少量级？
4. 单 cell 中 R10 的实际 O2 消耗量是否合理（与其他反应相比）？

Source: Hamel (1999) Eq.5.59, Table 5.4; Siminski (1972)
"""

from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve()))

import numpy as np
from src.kinetics.tar_reactions import (
    rate_R10, _r10_single_class, calc_tar_surrogate_fractions,
    R10_k0_AROM, R10_E_Rg_AROM, R10_k0_OLEF, R10_E_Rg_OLEF
)
from src.kinetics.arrhenius import k_hobbs
from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)
from contextlib import contextmanager

print("\n" + "="*70)
print("R10 TAR OXIDATION AUDIT: Pressure Scaling, Rate Magnitude, O2 Budget")
print("="*70)

# ─────────────────────────────────────────────────────────────────────────
# Part 1: Pressure scaling factor analysis
# ─────────────────────────────────────────────────────────────────────────
print("\n[Part 1] Pressure Scaling Factor P^0.3")
print("-" * 70)

P_Pa = 2_500_000.0  # 25 bar, nominal bed pressure
P_factor = P_Pa ** 0.3
print(f"P = {P_Pa:.2e} Pa (≈ {P_Pa/101325:.1f} bar)")
print(f"P^0.3 = {P_factor:.6f}")
print(f"  → If P were in bar: (25)^0.3 = {25**0.3:.6f}")
print(f"  → If P were in atm: (24.7)^0.3 = {24.7**0.3:.6f}")
print()
print("⚠ NOTE: P^0.3 对单位敏感（Pa 与 bar/atm 会产生数量级差异）")
print("   在 2.5e6 Pa 条件下，P^0.3 ≈ 83；若误用 25(bar) 仅为 ≈2.63。")
print("   For comparison:")
for P_test, label in [(1e5, "1 bar (Pa)"), (1e6, "10 bar (Pa)"), (2.5e6, "25 bar (Pa)"),
                       (1.0, "1 bar (dimensionless)"), (25.0, "25 bar (dimensionless)")]:
    print(f"     P = {P_test:.2e} → P^0.3 = {P_test**0.3:.6f}")

# ─────────────────────────────────────────────────────────────────────────
# Part 2: Kinetic constant analysis
# ─────────────────────────────────────────────────────────────────────────
print("\n[Part 2] Pre-exponential Factors and Activation Energies")
print("-" * 70)
print(f"Aromatic (C6H6, C10H8):")
print(f"  k0 = {R10_k0_AROM:.2e}  (Siminski 1972 via Hamel)")
print(f"  E/Rg = {R10_E_Rg_AROM:.1f} K")
print()
print(f"Olefin/Alkane (C16H34):")
print(f"  k0 = {R10_k0_OLEF:.2e}")
print(f"  E/Rg = {R10_E_Rg_OLEF:.1f} K")
print()
print(f"Ratio: k0_arom / k0_olef = {R10_k0_AROM / R10_k0_OLEF:.0f}x")

# ─────────────────────────────────────────────────────────────────────────
# Part 3: Rate constant as function of temperature
# ─────────────────────────────────────────────────────────────────────────
print("\n[Part 3] Kinetic Rate Constant k(T) for Aromatics")
print("-" * 70)
print("Evaluating k = k_hobbs(k0, E, T) * P^0.3 across typical bed temperatures:\n")

T_test_list = [900, 1100, 1200, 1300, 1400]
P_factor = P_Pa ** 0.3

for T_test in T_test_list:
    k_arom = k_hobbs(R10_k0_AROM, R10_E_Rg_AROM * Rg, T_test)  # Note: E in J/mol
    k_arom_with_P = k_arom * P_factor
    k_olef = k_hobbs(R10_k0_OLEF, R10_E_Rg_OLEF * Rg, T_test)
    k_olef_with_P = k_olef * P_factor
    print(f"T = {T_test:4d} K:")
    print(f"  k_arom(T) = {k_arom:.4e},  k_arom(T)·P^0.3 = {k_arom_with_P:.4e}")
    print(f"  k_olef(T) = {k_olef:.4e},  k_olef(T)·P^0.3 = {k_olef_with_P:.4e}")
    print()

# ─────────────────────────────────────────────────────────────────────────
# Part 4: Synthetic concentration scenarios
# ─────────────────────────────────────────────────────────────────────────
print("\n[Part 4] Synthetic Rate Calculations Across Concentration Scenarios")
print("-" * 70)

T = 1200.0  # K, typical lower bed
C_tar_scenarios = [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]  # mol/m³
C_O2_scenarios = [0.1, 0.5, 1.0, 2.0]  # mol/m³

print(f"\nTemperature: T = {T} K, P = {P_Pa:.2e} Pa, fuel_type = 'coal'\n")
print("Aromatic tar (default coal surrogate mixture):")
print("R10 rate [mol/(m³·s)] as function of C_tar and C_O2:\n")

header = "C_tar (mol/m³) │ " + " │ ".join([f"C_O2={c:.1f}" for c in C_O2_scenarios])
print(header)
print("-" * len(header))

for C_tar in C_tar_scenarios:
    rates = []
    for C_O2 in C_O2_scenarios:
        r = rate_R10(T, C_tar, C_O2, P_Pa, "coal")
        rates.append(f"{r:.3e}")
    line = f"{C_tar:>14.3f} │ " + " │ ".join([f"{r:>12s}" for r in rates])
    print(line)

print("\n" + "="*70)
print("INTERPRETATION:")
print("="*70)
print("""
If R10 rates in the above table are in the range 1e4 - 1e6 mol/(m³·s), then:
  → R10 is UNREASONABLY FAST (single cell would be O2-limited instantly)
  
If R10 rates are in the range 1e-2 - 1e2 mol/(m³·s), then:
  → R10 is REASONABLE (comparable to char combustion R1, CO oxidation R5)
  
If R10 rates are <1e-4 mol/(m³·s), then:
  → R10 is NEGLIGIBLE (tar oxidation doesn't drive the process)
""")

# ─────────────────────────────────────────────────────────────────────────
# Part 5: Single-cell context from actual gasifier run
# ─────────────────────────────────────────────────────────────────────────
print("\n[Part 5] Single-Cell Rates from Actual Gasifier Simulation")
print("-" * 70)

def _prepare_cell_minimal(scenario, i: int):
    """Load cell state from scenario using current global-NR reactor API."""
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    cfg.gas_inlet_dense_frac = float(scenario.dense)
    cfg.heat_loss_frac = float(scenario.heat_loss)
    cfg.r5_scale = float(scenario.r5s)
    cfg.r6_scale = float(scenario.r6s)
    cfg.r4_scale = float(scenario.r4s)
    reactor = Reactor(cfg)
    solve_kwargs = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    solve_kwargs["max_global_iter"] = 6
    solve_kwargs["tol_global"] = 1.0
    reactor.solve(**solve_kwargs)

    # Extract cell i
    return reactor.cells[i], reactor

from scripts.audit_single_cell_oxygen_convergence import Scenario as ScenarioAudit

# Stable working point from earlier audits
sc_stable = ScenarioAudit(
    name='stable_working_point',
    dense=0.40, heat_loss=0.10,
    r5s=0.50, r6s=1.00, r4s=1.50
)

print(f"Testing scenario: {sc_stable.name}")
print(f"  dense_frac={sc_stable.dense}, heat_loss={sc_stable.heat_loss}")
print(f"  R5×{sc_stable.r5s}, R6×{sc_stable.r6s}, R4×{sc_stable.r4s}\n")

print("Lower-cell R10 rate prediction in actual gasifier context:\n")
print("Cell │ dh (m) │   T (K)  │ C_tar_b │ C_O2_b │ R10_b rate │ C_tar_d │ C_O2_d │ R10_d rate")
print("─" * 95)

try:
    for i in range(3):
        cell, _ = _prepare_cell_minimal(sc_stable, i)
        
        # Concentrations
        C_b = cell._concentrations("b")
        C_d = cell._concentrations("d")
        C_tar_b = cell.N_b[idx["TAR1"]] + cell.N_b[idx["TAR2"]]
        C_tar_d = cell.N_d[idx["TAR1"]] + cell.N_d[idx["TAR2"]]
        C_O2_b = C_b[idx["O2"]]
        C_O2_d = C_d[idx["O2"]]
        
        # R10 rates
        r10_b = rate_R10(cell.T, C_tar_b, C_O2_b, cell.P, "coal")
        r10_d = rate_R10(cell.T, C_tar_d, C_O2_d, cell.P, "coal")
        
        print(f"  {i}  │ {cell.geo.dh:.3f} │ {cell.T:8.1f} │ {C_tar_b:7.3e} │ {C_O2_b:6.3e} │ {r10_b:10.3e} │ {C_tar_d:7.3e} │ {C_O2_d:6.3e} │ {r10_d:10.3e}")

except Exception as e:
    print(f"⚠ Could not load single-cell data: {e}")
    print("   (This is okay—will diagnose based on synthetic scenarios above)")

print("\n" + "="*70)
print("SUMMARY & DIAGNOSIS")
print("="*70)
print("""
The R10 audit will help identify:

1. Is P^0.3 calculation correct (Pa vs dimensionless)?
2. Is the aromatic tar surrogate weight reasonable for coal?
3. In the 1200–1400 K range typical of HTW gasifier, what is k(T)·P^0.3?
4. Given measured tar and O2 concentrations in single cells, is R10 rate
   physically plausible, or does it dominate unexpectedly?

If R10 rates at single-cell conditions (e.g., 0.2 mol/m³ tar, 1 mol/m³ O2)
are found to exceed 1e4 mol/(m³·s):
  → Pressure scaling likely uses wrong units (should be Pa^0.3 ≠ bar^0.3)
  → OR k0 values for aromatics are too large
  → OR the model itself has conceptual flaws (tar shouldn't react this fast)

Next step: Commit findings and propose corrected R10 if needed.
""")
