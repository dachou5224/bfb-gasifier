#!/usr/bin/env python3
"""Analyze k10 units and magnitude"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve()))

from src.kinetics.tar_reactions import R10_k0_AROM, R10_E_Rg_AROM
from src.kinetics.arrhenius import k_hobbs
from src.core.constants import Rg

print("R10 Kinetic Analysis: Is k10 too large?")
print("=" * 70)

k10 = R10_k0_AROM
E_Rg = R10_E_Rg_AROM
T = 1200.0
P_bar = 25.0
C_tar = 0.2
C_O2 = 1.0

# Calculate
k_base = k_hobbs(k10, E_Rg * Rg, T)
p_factor = P_bar ** 0.3
rate = k_base * p_factor * (C_tar ** 0.5) * C_O2

print(f"\nNominal R10 calculation (aromatic, 1200K, 25 bar):")
print(f"  k_hobbs(k10={k10}, E/Rg={E_Rg}K, T={T}K) = {k_base:.4e}")
print(f"  P_bar^0.3 = {p_factor:.4f}")
print(f"  C_tar^0.5 = {C_tar**0.5:.4f}")
print(f"  C_O2 = {C_O2:.4f}")
print(f"  → R10 = {rate:.4e} mol/(m³·s)")

print(f"\nDimensional analysis:")
print(f"  Formula: R = k_hobbs(k10, E, T) · P^0.3 · C_tar^0.5 · C_O2")
print(f"  If R should be [mol/(m³·s)], then:")
print(f"  k_hobbs · P^0.3 · C_tar^0.5 · C_O2")
print(f"  = (k10·T·exp) · [bar]^0.3 · [mol/m³]^0.5 · [mol/m³]")
print(f"  = k10 · [K] · [bar]^0.3 · [mol^1.5/m^4.5]")
print(f"\n  For this to equal [mol/(m³·s)], we need:")
print(f"  [k10] = [mol/(m³·s)] / ([K] · [bar]^0.3 · [mol^1.5/m^4.5])")
print(f"  [k10] = 1 / ([K] · [bar]^0.3) · [m^1.5/mol^0.5]")
print(f"\n  But k10 = {k10} is a pure number!")
print(f"  → k10 must have hidden units OR the formula is missing something")

print(f"\n\nHypothesis 1: k10 should be much smaller")
target_rate = 100.0
required_k_hobbs = target_rate / (p_factor * (C_tar**0.5) * C_O2)
required_k10 = required_k_hobbs / T

print(f"  If target R10 = {target_rate} mol/(m³·s), then:")
print(f"    Required k_hobbs ~ {required_k_hobbs:.2e}")
print(f"    Required k10 ~ {required_k10:.2e}")
print(f"    Current k10 / Required k10 ~ {k10 / required_k10:.0e}x")

print(f"\nHypothesis 2: Hamel formula interpretation error")
print(f"  Perhaps Eq.5.59 uses different concentration/pressure units?")
print(f"  Or maybe tar model (TAR1+TAR2) has different meaning than Hamel's definition?")

print(f"\nRecommendation for NEXT STEP:")
print(f"  1. Verify Hamel dissertation §5.2.6, Table 5.4 for exact formula & units")
print(f"  2. If k10=20700 is correct, then maybe R10 SHOULD be ~ 1e4?")
print(f"     (meaning tar oxidation is TOO FAST and shouldn't be in model)")
print(f"  3. Or reduce k10 by 100-200x to make R10 reasonable")
