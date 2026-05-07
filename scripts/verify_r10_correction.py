#!/usr/bin/env python3
"""验证 R10 100x 缩小后的效果"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve()))

from src.kinetics.tar_reactions import R10_k0_AROM, R10_k0_OLEF
from src.kinetics.arrhenius import k_hobbs
from src.core.constants import Rg

print("\n" + "="*70)
print("R10 AFTER 100x REDUCTION: Verify Rate Magnitude")
print("="*70 + "\n")

# Test new values
P_Pa = 2_500_000.0
T = 1200.0
P_bar = P_Pa / 100_000.0
C_tar = 0.2
C_O2 = 1.0
E_Rg_arom = 9_650.0

k_base_arom = k_hobbs(R10_k0_AROM, E_Rg_arom * Rg, T)
P_factor = P_bar ** 0.3
rate_arom = k_base_arom * P_factor * (C_tar ** 0.5) * C_O2

print(f"After 100x reduction:")
print(f"  k10_AROM = {R10_k0_AROM:.2f} (was 20700, now {R10_k0_AROM})")
print(f"  k10_OLEF = {R10_k0_OLEF:.4f} (was 59.8, now {R10_k0_OLEF})")
print()
print(f"Test conditions (typical lower cell): T={T}K, P={P_bar:.0f}bar, C_tar={C_tar}, C_O2={C_O2}")
print()
print(f"  k_hobbs(T={T}) = {k_base_arom:.4e}")
print(f"  P_bar^0.3 = {P_factor:.4f}")
print(f"  R10_arom @ C_tar={C_tar}, C_O2={C_O2} = {rate_arom:.4e} mol/(m³·s)")
print()

if 50 < rate_arom < 1000:
    print(f"✓ PASS: R10 now in reasonable range ({rate_arom:.0f} mol/(m³·s))")
    print("        Comparable to R6 (typically 200-1000) and R5 (typically 50-200)")
    print("        This should allow single cells to converge better")
else:
    print(f"⚠ WARNING: R10 = {rate_arom:.1e} might still be off-target")

print("\n" + "="*70)
print("Next: Run full single-cell diagnostic")
print("="*70)
