"""数量级验证（Phase 5.4）。

每个模块完成后必须通过以下 6 项检验。
用 assert 形式编写，可直接 python tests/sanity_checks.py 运行。

Source: docs/CLAUDE.md §数量级验证
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np


def check_u_mf():
    """u_mf: 褐煤 d_p=0.5mm, rho_s=1000, T=900K, P=2.5MPa => 0.02-0.08 m/s。"""
    from src.physics.minimum_fluidization import compute_u_mf
    from src.core.species import gas_viscosity_power_law, gas_density_ideal

    P, T = 2.5e6, 1173.15
    rho_g = P * 0.029 / (8.314 * T)  # 空气近似
    mu_g = gas_viscosity_power_law(T, 1.8e-5)
    u_mf = compute_u_mf(rho_g, 1000.0, 0.5e-3, mu_g)
    print(f"  u_mf = {u_mf:.4f} m/s  [0.02-0.08]", end=" ")
    assert 0.02 <= u_mf <= 0.08, f"u_mf={u_mf} out of range"
    print("PASS")


def check_d_b():
    """d_b(H=5m): u0=0.3, P=2.5MPa => 0.05-0.3 m。"""
    from src.physics.bubble_dynamics import integrate_bubble_diameter

    _, db = integrate_bubble_diameter(0.3, 0.05, 2.5e6, 5.0, D_bed=0.6)
    db_top = db[-1]
    print(f"  d_b(5m) = {db_top:.4f} m  [0.05-0.3]", end=" ")
    assert 0.05 <= db_top <= 0.3, f"d_b={db_top} out of range"
    print("PASS")


def check_K_bd():
    """K_bd: d_b=0.1m, P=2.5MPa, T=1000K => 1-15 s⁻¹。"""
    from src.physics.mass_transfer import calc_kbd, calc_u_br
    from src.core.species import gas_diffusivity_correlation

    P, T = 2.5e6, 1000.0
    D_g = gas_diffusivity_correlation(T, P)
    u_d = 0.05 / 0.45
    u_br = calc_u_br(u_d, P)

    from src.physics.bubble_dynamics import bubble_rise_velocity
    u_b = bubble_rise_velocity(0.3, 0.05, 0.1)
    K = calc_kbd(u_br, 0.1, D_g, 0.45, u_b)
    print(f"  K_bd = {K:.2f} 1/s  [1-15]", end=" ")
    assert 1.0 <= K <= 15.0, f"K_bd={K} out of range"
    print("PASS")


def check_R_boudouard():
    """R4 Boudouard: T=1073K, P_CO2=5e4 Pa => ~1e-5 mol/(m²·s)。"""
    from src.kinetics.char_reactions import rate_R4

    r = rate_R4(1073.0, 5e4, 1e4)
    print(f"  R_boud = {r:.2e} mol/(m²·s)  [1e-7, 1e-3]", end=" ")
    assert 1e-7 < r < 1e-3, f"R4={r} out of range"
    print("PASS")


def check_R8_direction():
    """R8 WGSR: T=1073K, P=2.5MPa, y_CO=0.25 => 正向（朝平衡方向）。"""
    from src.kinetics.gas_reactions import rate_R8

    r = rate_R8(1073.0, 2.5e6, 0.25, 0.10, 0.05, 0.15)
    print(f"  R8 = {r:.4f} mol/(m³·s)  [>0 = forward]", end=" ")
    assert r > 0, f"R8={r} should be positive (forward)"
    print("PASS")


def check_DAEM():
    """DAEM: 升温至 900°C，100s => 挥发分释放率 > 80%。"""
    from src.thermal.devolatilization import daem_conversion, _linear_heating_profile

    T_h, t_h = _linear_heating_profile(300.0, 1173.15, 100.0, 500)
    X = daem_conversion(T_h, t_h)
    print(f"  DAEM X_VM = {X*100:.1f}%  [>80%]", end=" ")
    assert X > 0.80, f"X_VM={X} < 0.80"
    print("PASS")


def check_hv_corrected():
    """Eq.4.4: 修正蒸发焓 h_v' 应大于 h_v。"""
    from src.thermal.drying import corrected_evaporation_enthalpy, H_EVAP

    hvp = corrected_evaporation_enthalpy(w0_tr=0.169, T0=300.0, T_e=373.15)
    print(f"  h_v' = {hvp:.3e} J/kg  [> h_v={H_EVAP:.3e}]", end=" ")
    assert hvp > H_EVAP, "corrected evaporation enthalpy must be larger than latent heat"
    print("PASS")


def check_nusselt_htc():
    """Eq.4.9: Nu>2 且 h_conv 量级合理。"""
    from src.thermal.drying import nusselt_particle, convective_htc_from_nusselt

    Re, Pr = 50.0, 0.7
    Nu = nusselt_particle(Re=Re, Pr=Pr)
    h = convective_htc_from_nusselt(
        d_p=1.0e-3,
        u_rel=1.0,
        rho_g=0.35,
        mu_g=4.0e-5,
        cp_g=1200.0,
        lambda_g=0.08,
    )
    print(f"  Nu = {Nu:.2f}  [>2], h_conv = {h:.1f} W/(m²·K)  [10-5000]", end=" ")
    assert Nu > 2.0, f"Nu={Nu} should be >2"
    assert 10.0 <= h <= 5000.0, f"h_conv={h} out of range"
    print("PASS")


def check_daem_radial():
    """Eq.4.12: 径向体积分 DAEM 返回 [0,1]。"""
    from src.thermal.devolatilization import daem_conversion_radial

    r = np.linspace(0.0, 1.0e-3, 8)
    t = np.linspace(0.0, 50.0, 120)
    T_rt = np.vstack([700.0 + 300.0 * (ri / r[-1]) + 120.0 * (t / t[-1]) for ri in r])
    x = daem_conversion_radial(T_rt, t, r)
    print(f"  X_VM,radial = {x*100:.1f}%  [0-100%]", end=" ")
    assert 0.0 <= x <= 1.0, f"X_VM_radial={x} out of range"
    print("PASS")


def check_pyrolysis_elemental_allocator():
    """热解元素守恒分配：C/H/O 不超支。"""
    from src.core.cell import Cell, SolidProps
    from src.core.species import get_tar_component_mapping, TAR_SURROGATE_FORMULA

    c = Cell(solid=SolidProps())
    c.fuel_type = "coal"
    nC_in, nH_in, nO_in = 10.0, 12.0, 4.0
    prod = c._allocate_pyrolysis_products_elemental(nC=nC_in, nH=nH_in, nO=nO_in)

    tar_map = get_tar_component_mapping("coal")
    c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map["TAR1"]]
    nC_out = prod["CO"] + prod["CO2"] + prod["CH4"] + c_tar * prod["TAR1"]
    nH_out = 2.0 * prod["H2"] + 2.0 * prod["H2O"] + 4.0 * prod["CH4"] + h_tar * prod["TAR1"]
    nO_out = prod["CO"] + 2.0 * prod["CO2"] + prod["H2O"]

    print(
        f"  C/H/O out = ({nC_out:.3f},{nH_out:.3f},{nO_out:.3f}) "
        f"<= in ({nC_in:.3f},{nH_in:.3f},{nO_in:.3f})",
        end=" "
    )
    assert nC_out <= nC_in + 1e-9
    assert nH_out <= nH_in + 1e-9
    assert nO_out <= nO_in + 1e-9
    assert all(v >= 0.0 for v in prod.values())
    print("PASS")


def check_drying_pyrolysis_source_hook():
    """Cell 源项链：干燥/热解应产生气相源项并扣减固相。"""
    from src.core.cell import Cell, SolidProps

    c = Cell(solid=SolidProps())
    c.T = 1173.15
    c.P = 2.5e6
    c.m_solid_zu[0] = 0.1  # 给定固体进料以触发源项
    src = c._calc_drying_pyrolysis_gas_source(tau_cell=5.0)

    src_sum = float(np.sum(np.maximum(src, 0.0)))
    rs_sum = float(np.sum(c.R_solid))
    print(f"  gas_source_sum = {src_sum:.3e} mol/s, R_solid_sum = {rs_sum:.3e} kg/s", end=" ")
    assert src_sum > 0.0, "drying/pyrolysis gas source should be positive"
    assert rs_sum < 0.0, "solid source should be negative (mass removed)"
    print("PASS")


def main():
    checks = [
        ("1. u_mf", check_u_mf),
        ("2. d_b(H_bed)", check_d_b),
        ("3. K_bd", check_K_bd),
        ("4. R_boudouard", check_R_boudouard),
        ("5. R8 WGSR direction", check_R8_direction),
        ("6. DAEM", check_DAEM),
        ("7. h_v corrected (Eq.4.4)", check_hv_corrected),
        ("8. Nu & h_conv (Eq.4.9)", check_nusselt_htc),
        ("9. DAEM radial (Eq.4.12)", check_daem_radial),
        ("10. Pyrolysis elemental allocator", check_pyrolysis_elemental_allocator),
        ("11. Drying-pyrolysis source hook", check_drying_pyrolysis_source_hook),
    ]

    all_pass = True
    for name, fn in checks:
        print(f"[{name}]")
        try:
            fn()
        except AssertionError as e:
            print(f"  FAIL: {e}")
            all_pass = False
        except Exception as e:
            print(f"  ERROR: {e}")
            all_pass = False

    print("\n" + "=" * 40)
    if all_pass:
        print(f"ALL {len(checks)} SANITY CHECKS PASSED")
    else:
        print("SOME CHECKS FAILED")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
