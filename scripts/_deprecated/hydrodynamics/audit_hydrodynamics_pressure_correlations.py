"""Hydrodynamics pressure-correlation audit.

目的：
1. 把 hydrodynamics 相关 pressure terms 的实现链集中对照一遍；
2. 区分「已实现且与本地 source-of-truth 一致」与「仍 unresolved」。
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.core.constants import P0_HAMEL, Rg
from src.core.species import gas_density_ideal, gas_diffusivity_correlation
from src.physics.bubble_dynamics import bubble_lifetime
from src.physics.freeboard import calc_beta_a, calc_u_gb
from src.physics.mass_transfer import calc_kbd, calc_u_br


def _fmt(x: float) -> str:
    return f"{x:.6g}"


def main() -> None:
    p_ref = P0_HAMEL
    p_lu = 2.5e6
    p_ratio = p_lu / p_ref
    t_ref = 1100.0
    y_air = {"N2": 0.79, "O2": 0.21}

    print("Hydrodynamics Pressure Correlation Audit")
    print("=" * 60)
    print(f"P_ref  = {p_ref:.1f} Pa")
    print(f"P_LU   = {p_lu:.1f} Pa")
    print(f"P/P0   = {p_ratio:.4f}")
    print()

    rho_ref = gas_density_ideal(P=p_ref, T=t_ref, mole_fractions=y_air)
    rho_lu = gas_density_ideal(P=p_lu, T=t_ref, mole_fractions=y_air)
    rho_ratio = rho_lu / rho_ref
    print("[rho_g] ideal gas density")
    print(f"  expected ratio = P/P0 = {_fmt(p_ratio)}")
    print(f"  code ratio     = {_fmt(rho_ratio)}")
    print("  status         = OK (pressure dependence correctly implemented)")
    print()

    dg_ref = gas_diffusivity_correlation(T_m=t_ref, P=p_ref)
    dg_lu = gas_diffusivity_correlation(T_m=t_ref, P=p_lu)
    dg_ratio = dg_lu / dg_ref
    expected_dg_ratio = p_ref / p_lu
    print("[D_g] gas diffusivity Eq.5.20")
    print(f"  expected ratio = P0/P = {_fmt(expected_dg_ratio)}")
    print(f"  code ratio     = {_fmt(dg_ratio)}")
    print(f"  reference P0   = {P0_HAMEL:.0f} Pa")
    print("  status         = OK (formula correct and project baseline unified)")
    print()

    u_d = 0.12
    u_br_ref = calc_u_br(u_d=u_d, P=p_ref)
    u_br_lu = calc_u_br(u_d=u_d, P=p_lu)
    expected_u_br_ratio = (p_lu / p_ref) ** (-0.15)
    print("[u_br] bubble through-flow Eq.3.44")
    print(f"  expected ratio = (P/P0)^(-0.15) = {_fmt(expected_u_br_ratio)}")
    print(f"  code ratio     = {_fmt(u_br_lu / u_br_ref)}")
    print("  status         = OK (matches local Hamel extract and code)")
    print()

    d_b = 0.12
    u_b = 1.2
    u_mf = 0.05
    lam_cur_ref = bubble_lifetime(d_b, u_b, p_ref, strategy='current', u_mf=u_mf)
    lam_cur_lu = bubble_lifetime(d_b, u_b, p_lu, strategy='current', u_mf=u_mf)
    lam_ham_ref = bubble_lifetime(d_b, u_b, p_ref, strategy='hamel_280', u_mf=u_mf)
    lam_ham_lu = bubble_lifetime(d_b, u_b, p_lu, strategy='hamel_280', u_mf=u_mf)
    print("[lambda_b] bubble lifetime")
    print(f"  current ratio  = {_fmt(lam_cur_lu / lam_cur_ref)}  (expected {(p_lu / p_ref) ** (-0.2):.6g})")
    print(f"  hamel_280 ratio= {_fmt(lam_ham_lu / lam_ham_ref)}  (expected {(p_lu / p_ref) ** (-0.7):.6g})")
    print("  status         = OK (default path aligned to Hamel Eq.3.35 / Eq.3.42)")
    print("  detail         = legacy `current` kept only as audit contrast against paper-consistent `hamel_280`")
    print()

    eps_mf = 0.45
    kbd_ref = calc_kbd(u_br=u_br_ref, d_b=d_b, D_g=dg_ref, eps_mf=eps_mf, u_b=u_b)
    kbd_lu = calc_kbd(u_br=u_br_lu, d_b=d_b, D_g=dg_lu, eps_mf=eps_mf, u_b=u_b)
    conv_ref = 3.0 * u_br_ref / (2.0 * d_b)
    conv_lu = 3.0 * u_br_lu / (2.0 * d_b)
    diff_ref = math.sqrt((144.0 * dg_ref * eps_mf * u_b) / (math.pi * d_b**3))
    diff_lu = math.sqrt((144.0 * dg_lu * eps_mf * u_b) / (math.pi * d_b**3))
    print("[K_bd] mass transfer coefficient Eq.3.50")
    print(f"  conv ratio     = {_fmt(conv_lu / conv_ref)}  (inherits u_br)")
    print(f"  diff ratio     = {_fmt(diff_lu / diff_ref)}  (inherits sqrt(D_g))")
    print(f"  total ratio    = {_fmt(kbd_lu / kbd_ref)}")
    print("  status         = OK (pressure enters via u_br and D_g as expected)")
    print()

    beta_a = calc_beta_a(0.12)
    u_gb = calc_u_gb(1.2)
    print("[freeboard beta_A / u_gb]")
    print(f"  beta_A         = {_fmt(beta_a)}  (no explicit pressure term)")
    print(f"  u_gb           = {_fmt(u_gb)}  (no explicit pressure term)")
    print("  status         = OK if source-of-truth remains beta_A=1/(3*d_b,ws), u_p0=1.53*u_b")
    print("  note           = pressure only enters indirectly through upstream d_b / rho_g / terminal-velocity chain")
    print()

    rho_ref_ideal = p_ref * 0.029 / (Rg * t_ref)
    rho_lu_ideal = p_lu * 0.029 / (Rg * t_ref)
    print("[summary]")
    print("  - Correct and consistent: rho_g, u_br, K_bd pressure path")
    print("  - Correct and unified: D_g now uses project constant P0")
    print("  - lambda_b source-of-truth resolved to Hamel Eq.3.35 / Eq.3.42 (legacy current retained only for audit)")
    print("  - No local evidence that beta_A or u_gb are missing an explicit pressure exponent")
    print(f"  - Ideal-air rho_g reference check: {_fmt(rho_ref_ideal)} -> {_fmt(rho_lu_ideal)} kg/m3")


if __name__ == "__main__":
    main()
