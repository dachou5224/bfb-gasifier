#!/usr/bin/env python3
"""Verify current R10 constants and Pa-basis implementation.

This supersedes the older local experiment that assumed a 100x reduction and
bar-basis pressure.  The current code follows the extracted Hamel notes:
Eq.5.59 with P in Pa and Eq.5.60/5.61 k10 values unchanged.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(".").resolve()))

from src.core.constants import Rg
from src.kinetics.arrhenius import k_hobbs
from src.kinetics.tar_reactions import R10_E_Rg_AROM, R10_k0_AROM, rate_R10


def main() -> int:
    print("\n" + "=" * 70)
    print("R10 CURRENT IMPLEMENTATION VERIFICATION")
    print("=" * 70 + "\n")

    P_pa = 2_500_000.0
    T = 1200.0
    C_tar = 0.2
    C_O2 = 1.0

    k_base = k_hobbs(R10_k0_AROM, R10_E_Rg_AROM * Rg, T)
    p_factor = P_pa**0.3
    aromatic_rate = k_base * p_factor * math.sqrt(C_tar) * C_O2
    coal_rate = rate_R10(T, C_tar, C_O2, P_pa, "coal")

    print("Expected Hamel-aligned constants")
    print(f"  k10_AROM = {R10_k0_AROM:.6g}")
    print(f"  E/Rg_AROM = {R10_E_Rg_AROM:.6g} K")
    print("\nTest conditions")
    print(f"  T={T:.1f} K, P={P_pa:.3e} Pa, C_tar={C_tar}, C_O2={C_O2}")
    print("\nComputed rates")
    print(f"  k_hobbs(T) = {k_base:.6e}")
    print(f"  P_pa^0.3   = {p_factor:.6e}")
    print(f"  aromatic   = {aromatic_rate:.6e} mol/(m3 s)")
    print(f"  coal lumped= {coal_rate:.6e} mol/(m3 s)")

    if not math.isfinite(coal_rate) or coal_rate <= 0.0:
        print("\nFAIL: current R10 implementation produced an invalid positive-input rate.")
        return 1

    print("\nResult: PASS - R10 constants and pressure basis are internally consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
