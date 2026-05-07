"""Audit lambda_b pressure-correction over a representative pressure window.

目的：
1. 比较 current(-0.2) 与 hamel_280(-0.7) 在 0.1–2.5 MPa 范围内的差异；
2. 直接观察它们在 ODE 衰减项 d_b/(3*lambda_b*u_b) 上造成的量级变化；
3. 给后续锁定 source-of-truth 提供一个固定 pressure-window 证据面。
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.core.constants import P0_HAMEL, g
from src.physics.bubble_dynamics import bubble_lifetime


def main() -> None:
    pressures = [
        P0_HAMEL,
        5.0e5,
        1.0e6,
        1.5e6,
        2.0e6,
        2.5e6,
    ]
    d_b = 0.12
    u_b = 1.2
    u_mf = 0.05

    print("lambda_b Pressure Window Audit")
    print("=" * 84)
    print(
        f"{'P[MPa]':>8} {'P/P0':>8} {'lam_cur[s]':>12} {'lam_ham[s]':>12} "
        f"{'ham/cur':>10} {'dec_cur[1/m]':>14} {'dec_ham[1/m]':>14} {'ham/cur':>10}"
    )

    for p in pressures:
        lam_cur = float(bubble_lifetime(d_b, u_b, p, strategy="current", u_mf=u_mf))
        lam_ham = float(bubble_lifetime(d_b, u_b, p, strategy="hamel_280", u_mf=u_mf))
        dec_cur = d_b / (3.0 * lam_cur * u_b)
        dec_ham = d_b / (3.0 * lam_ham * u_b)
        print(
            f"{p/1e6:8.3f} {p/P0_HAMEL:8.3f} {lam_cur:12.5f} {lam_ham:12.5f} "
            f"{lam_ham/max(lam_cur,1e-12):10.4f} {dec_cur:14.5f} {dec_ham:14.5f} "
            f"{dec_ham/max(dec_cur,1e-12):10.4f}"
        )

    print()
    print("Reference formulas:")
    print("  current   : lambda_b = d_b/(0.5*u_b) * (P/P0)^(-0.2)")
    print("  hamel_280 : lambda_b = 280*u_mf/g * (P/P0)^(-0.7)")
    print()
    print("Interpretation:")
    print("  - 若 hamel_280 在高压下把 lambda_b 压到远低于 current，则 ODE 衰减项会同步被显著放大。")
    print("  - 这会直接限制 d_b(h) 的可达上限，因此 pressure correction 不只是二级细节，而是 bubble ODE 主项。")
    print(f"  - At u_mf={u_mf:.3f} m/s, Hamel base lifetime 280*u_mf/g = {280.0*u_mf/g:.5f} s at 1 atm.")


if __name__ == "__main__":
    main()
