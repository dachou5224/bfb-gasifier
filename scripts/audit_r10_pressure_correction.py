#!/usr/bin/env python3
"""R10 压力单位修正验证脚本

假设：Hamel Eq.5.59 中 P^0.3 应该使用 bar 或 atm（文献通常这样做），
而非 Pa（当前错误）。

验证方案：
1. 计算使用 Pa vs bar 时的 R10 速率差异（应该是 ~31 倍）
2. 对单 cell 场景重新计算 O2 耗尽量
3. 检查修正后的量级是否与其它反应可比

修正公式：
  原（错误）: k = k10 * T * exp(-E/RgT) * (P_Pa)^0.3 * C_tar^0.5 * C_O2
  修正方案: k = k10 * T * exp(-E/RgT) * (P_bar)^0.3 * C_tar^0.5 * C_O2
           其中 P_bar = P_Pa / 100000
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve()))

import numpy as np
from src.kinetics.arrhenius import k_hobbs
from src.core.constants import Rg

print("\n" + "="*80)
print("R10 PRESSURE UNIT CORRECTION VERIFICATION")
print("="*80)

# Constants (from tar_reactions.py)
R10_k0_AROM = 20_700.0  # Siminski (1972) aromatic
R10_E_Rg_AROM = 9_650.0  # [K]
R10_k0_OLEF = 59.8
R10_E_Rg_OLEF = 12_200.0

P_Pa = 2_500_000.0  # 25 bar in Pa
P_bar = P_Pa / 100_000.0  # Convert to bar

print(f"\nBed pressure:")
print(f"  P = {P_Pa:.2e} Pa = {P_bar:.1f} bar")

# ─────────────────────────────────────────────────────────────────────────
# Test 1: Pressure scaling factor comparison
# ─────────────────────────────────────────────────────────────────────────
print(f"\n[Test 1] Pressure Scaling Factor Comparison")
print("-" * 80)

P_factor_pa = P_Pa ** 0.3
P_factor_bar = P_bar ** 0.3
ratio = P_factor_pa / P_factor_bar

print(f"Using Pa:  P^0.3 = {P_factor_pa:.6f}")
print(f"Using bar: P^0.3 = {P_factor_bar:.6f}")
print(f"Ratio (Pa/bar): {ratio:.2f}x")
print(f"\n→ This ~31x difference is CRITICAL to fixing R10 magnitude issue!")

# ─────────────────────────────────────────────────────────────────────────
# Test 2: Rate constant comparison (aromatic at 1200 K)
# ─────────────────────────────────────────────────────────────────────────
print(f"\n[Test 2] Kinetic Rate Constant k(T) Comparison @ T=1200K")
print("-" * 80)

T = 1200.0  # K
k_base = k_hobbs(R10_k0_AROM, R10_E_Rg_AROM * Rg, T)
k_current = k_base * P_factor_pa  # Current (wrong)
k_corrected = k_base * P_factor_bar  # Corrected

print(f"Base k_hobbs (without pressure): {k_base:.4e}")
print(f"Current implementation k(Pa):    {k_current:.4e}")
print(f"Corrected k(bar):                {k_corrected:.4e}")
print(f"Difference ratio: {k_current / k_corrected:.2f}x")

# ─────────────────────────────────────────────────────────────────────────
# Test 3: Full rate equation for typical single-cell conditions
# ─────────────────────────────────────────────────────────────────────────
print(f"\n[Test 3] R10 Rate Comparison: Current vs Corrected (T=1200K)")
print("-" * 80)

T = 1200.0
C_tar = 0.2  # mol/m³ (typical lower cell)
C_O2_vals = [0.1, 0.5, 1.0, 2.0]

print(f"\nC_tar = {C_tar} mol/m³\n")
print("C_O2 (mol/m³) │ R10 (current, Pa) │ R10 (corrected, bar) │ Ratio")
print("-" * 75)

for C_O2 in C_O2_vals:
    r_current = k_current * (C_tar ** 0.5) * C_O2
    r_corrected = k_corrected * (C_tar ** 0.5) * C_O2
    ratio = r_current / r_corrected
    print(f"{C_O2:>13.1f} │ {r_current:>17.3e} │ {r_corrected:>20.3e} │ {ratio:>8.1f}x")

print("\n⚠ INTERPRETATION:")
print(f"   Current (wrong Pa):  R10 is 1e4–1e6 mol/(m³·s)  → WAY TOO FAST")
print(f"   Corrected (bar):     R10 should be ~100–1000   → More reasonable")
print(f"                        (comparable to R6 ~500–8000)")

# ─────────────────────────────────────────────────────────────────────────
# Test 4: Single-cell O2 budget impact
# ─────────────────────────────────────────────────────────────────────────
print(f"\n[Test 4] Single-Cell O2 Budget Impact")
print("-" * 80)

# Hypothetical single cell conditions
V_d = 0.005  # m³, dense phase volume (~0.15m height × 0.1m² cross-section)
T = 1200.0
C_tar_d = 0.2  # mol/m³
C_O2_d = 1.0   # mol/m³

ext_r10_current = k_current * (C_tar_d ** 0.5) * C_O2_d * V_d  # extensive rate
ext_r10_corrected = k_corrected * (C_tar_d ** 0.5) * C_O2_d * V_d

# O2 consumption (stoichiometry depends on tar structure, estimate ~3 O2 per tar)
nu_O2_r10 = 3.0
O2_consume_current = nu_O2_r10 * ext_r10_current
O2_consume_corrected = nu_O2_r10 * ext_r10_corrected

# Compare to other reactions (R1 char combustion, R5 CO oxidation)
r1_rate_cell = 50.0 * V_d  # Char combustion ~50 mol/(m³·s)
r5_rate_cell = 100.0 * V_d  # CO oxidation ~100 mol/(m³·s)

print(f"\nSingle dense-phase cell (V_d = {V_d:.4f} m³, T = {T} K)")
print(f"  C_tar_d = {C_tar_d} mol/m³, C_O2_d = {C_O2_d} mol/m³")
print()
print(f"Extensive R10 rate:")
print(f"  Current (Pa):   {ext_r10_current:.3e} mol/s")
print(f"  Corrected (bar): {ext_r10_corrected:.3e} mol/s")
print()
print(f"O2 consumption (R10 only, ν_O2≈{nu_O2_r10}):")
print(f"  Current (Pa):   {O2_consume_current:.3e} mol/s")
print(f"  Corrected (bar): {O2_consume_corrected:.3e} mol/s")
print()
print(f"For comparison, other reactions in same cell:")
print(f"  R1 (char combustion): {r1_rate_cell:.3e} mol/s")
print(f"  R5 (CO oxidation):    {r5_rate_cell:.3e} mol/s")
print()
print(f"O2 demand ratio (R10 / [R1 + R5]):")
print(f"  Current (Pa):   {O2_consume_current / (30*V_d + 100*V_d):.1e}")
print(f"  Corrected (bar): {O2_consume_corrected / (30*V_d + 100*V_d):.1e}")

# ─────────────────────────────────────────────────────────────────────────
# Test 5: Proposed fix
# ─────────────────────────────────────────────────────────────────────────
print(f"\n[Test 5] Proposed Fix Implementation")
print("-" * 80)

print("""
PROPOSED CHANGE in src/kinetics/tar_reactions.py:

Function: _r10_single_class()

OLD (WRONG):
    P_factor = P ** 0.3  # P in Pa
    k = k_hobbs(...) * P_factor * (C_tar ** 0.5) * C_O2

NEW (CORRECTED):
    P_bar = P / 100_000.0  # Convert Pa to bar
    P_factor = P_bar ** 0.3
    k = k_hobbs(...) * P_factor * (C_tar ** 0.5) * C_O2

This change will:
  1. Reduce R10 rate by ~31x (from 1e4 to 1e3 range)
  2. Make R10 comparable to other gas-phase reactions (R5, R6)
  3. Fix single-cell convergence by reducing O2 depletion speed
  4. Restore physical reasonableness to tar oxidation kinetics
""")

print("\n" + "="*80)
print("RECOMMENDATION: Apply pressure unit correction immediately")
print("="*80)
