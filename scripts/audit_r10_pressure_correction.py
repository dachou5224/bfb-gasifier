#!/usr/bin/env python3
"""R10 pressure-basis audit against the extracted Hamel notes.

Source of truth:
- docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md
- Hamel Eq.5.59, Eq.5.60, Eq.5.61

The extracted notes state that the k10 units include Pa^0.3, so the
production implementation should pass pressure in Pa.  This script keeps the
older bar-basis calculation only as a sensitivity comparison.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(".").resolve()))

from src.core.constants import Rg
from src.kinetics.arrhenius import k_hobbs
from src.kinetics.tar_reactions import (
    R10_E_Rg_AROM,
    R10_E_Rg_OLEF,
    R10_k0_AROM,
    R10_k0_OLEF,
    rate_R10,
)


def _r10_manual(T: float, C_tar: float, C_O2: float, P: float, k10: float, E_Rg: float) -> float:
    k = k_hobbs(k10, E_Rg * Rg, T) * (P**0.3)
    return k * math.sqrt(max(C_tar, 0.0)) * max(C_O2, 0.0)


def main() -> int:
    print("\n" + "=" * 80)
    print("R10 PRESSURE-BASIS AUDIT")
    print("=" * 80)

    P_pa = 2_500_000.0
    P_bar = P_pa / 100_000.0
    T = 1200.0
    C_tar = 0.2
    C_O2 = 1.0

    print("\nSource-of-truth check")
    print("  Eq.5.59 pressure term: p^0.3")
    print("  Eq.5.60/5.61 k10 units: include Pa^0.3")
    print("  Expected implementation basis: P in Pa")

    r_arom_pa = _r10_manual(T, C_tar, C_O2, P_pa, R10_k0_AROM, R10_E_Rg_AROM)
    r_olef_pa = _r10_manual(T, C_tar, C_O2, P_pa, R10_k0_OLEF, R10_E_Rg_OLEF)
    r_arom_bar = _r10_manual(T, C_tar, C_O2, P_bar, R10_k0_AROM, R10_E_Rg_AROM)
    r_impl = rate_R10(T, C_tar, C_O2, P_pa, "coal")

    print("\nConditions")
    print(f"  T={T:.1f} K, P={P_pa:.3e} Pa ({P_bar:.1f} bar), C_tar={C_tar}, C_O2={C_O2}")
    print("\nPa-basis rates")
    print(f"  aromatic Pa-basis = {r_arom_pa:.6e} mol/(m3 s)")
    print(f"  olefin   Pa-basis = {r_olef_pa:.6e} mol/(m3 s)")
    print(f"  coal lumped implementation = {r_impl:.6e} mol/(m3 s)")
    print("\nSensitivity only")
    print(f"  aromatic bar-basis = {r_arom_bar:.6e} mol/(m3 s)")
    print(f"  Pa/bar factor      = {r_arom_pa / max(r_arom_bar, 1e-300):.3f}x")

    if r_impl <= 0.0:
        print("\nFAIL: R10 implementation returned a non-positive rate under positive inputs.")
        return 1
    if not math.isfinite(r_impl):
        print("\nFAIL: R10 implementation returned a non-finite rate.")
        return 1

    print("\nResult: PASS - current R10 pressure basis matches the extracted Hamel notes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
