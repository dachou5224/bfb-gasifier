"""数量级验证（Sanity Gate）。

定位：
- 模块级 smoke gate（快速发现量级/接线/默认口径回退）
- 不替代 cell 守恒审计、LU 端到端验证、profile 拟合
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def check_u_mf() -> None:
    """u_mf 量级检查（HTW 2.5 MPa，热态近似）"""
    from src.core.constants import Rg
    from src.core.species import gas_viscosity_power_law
    from src.physics.minimum_fluidization import compute_u_mf

    p, t = 2.5e6, 1200.0
    rho_g = p * 0.029 / (Rg * t)
    mu_g = gas_viscosity_power_law(t, 1.8e-5)
    u_mf = compute_u_mf(rho_g, 1200.0, 1.0e-3, mu_g, eps_mf=0.45, phi_s=0.75)
    print(f"  u_mf = {u_mf:.4f} m/s  [0.10-0.25]", end=" ")
    assert 0.10 <= u_mf <= 0.25
    print("PASS")


def check_d_b() -> None:
    """d_b(H=5m) 量级检查"""
    from src.physics.bubble_dynamics import integrate_bubble_diameter

    _, db = integrate_bubble_diameter(0.3, 0.05, 2.5e6, 5.0, D_bed=0.6)
    db_top = float(db[-1])
    print(f"  d_b(5m) = {db_top:.4f} m  [0.05-0.3]", end=" ")
    assert 0.05 <= db_top <= 0.3
    print("PASS")


def check_kbd_hamel_chain() -> None:
    """K_bd 使用 Hamel 链路（u_d closure + Eq.3.44 + Eq.3.50）"""
    from src.core.cell_hydrodynamics import calc_cell_hydrodynamics
    from src.core.species import N_GAS

    n_d = np.zeros(N_GAS, dtype=np.float64)
    n_b = np.zeros(N_GAS, dtype=np.float64)
    # 给少量主要组分，保证气体物性可解
    n_d[0], n_d[1], n_d[2], n_d[3], n_d[6] = 1.0, 0.6, 1.2, 1.4, 3.0

    bundle = calc_cell_hydrodynamics(
        T=1073.15,
        P=2.5e6,
        N_b=n_b,
        N_d=n_d,
        D_bed=0.6,
        dh=0.5,
        h_center=2.5,
        N_or=120,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=0.30,
        u_d_closure="wein_1992_eq312",
        bubble_diameter_model="hilligardt_ode",
        psi_b_strategy="wein_1992",
        lambda_strategy="hamel_280",
        xi_strategy="hamel_regime",
        bubble_velocity_strategy="heinbockel_eq343",
        bubble_ode_strategy="heinbockel_eq341",
    )
    print(f"  u_d = {bundle.u_d:.4f} m/s, K_bd = {bundle.K_bd:.2f} 1/s  [1-15]", end=" ")
    assert bundle.u_d > 0.0
    assert 1.0 <= bundle.K_bd <= 15.0
    print("PASS")


def check_r_boudouard() -> None:
    """R3 Boudouard（实现入口 rate_R4_effective）量级"""
    from src.kinetics.char_reactions import rate_R4

    r = rate_R4(1073.0, 5e4, 1e4)
    print(f"  R3_boud (impl r4) = {r:.2e} mol/(m²·s)  [1e-7, 1e-3]", end=" ")
    assert 1e-7 < r < 1e-3
    print("PASS")


def check_r8_direction() -> None:
    """R8 WGSR 在给定组分下应正向"""
    from src.kinetics.gas_reactions import rate_R8

    r = rate_R8(1073.0, 2.5e6, 0.25, 0.10, 0.05, 0.15)
    print(f"  R8 = {r:.4f} mol/(m³·s)  [>0 = forward]", end=" ")
    assert r > 0.0
    print("PASS")


def check_daem() -> None:
    """DAEM 挥发分释放率量级"""
    from src.thermal.devolatilization import _linear_heating_profile, daem_conversion

    t_h, tau_h = _linear_heating_profile(300.0, 1173.15, 100.0, 500)
    x_vm = daem_conversion(t_h, tau_h)
    print(f"  DAEM X_VM = {x_vm*100:.1f}%  [>80%]", end=" ")
    assert x_vm > 0.80
    print("PASS")


def check_hv_corrected() -> None:
    """Eq.4.4 修正蒸发焓应高于潜热"""
    from src.thermal.drying import H_EVAP, corrected_evaporation_enthalpy

    hvp = corrected_evaporation_enthalpy(w0_tr=0.169, T0=300.0, T_e=373.15)
    print(f"  h_v' = {hvp:.3e} J/kg  [> h_v={H_EVAP:.3e}]", end=" ")
    assert hvp > H_EVAP
    print("PASS")


def check_nusselt_htc() -> None:
    """Eq.4.9 Nu/h_conv 量级"""
    from src.thermal.drying import convective_htc_from_nusselt, nusselt_particle

    nu = nusselt_particle(Re=50.0, Pr=0.7)
    h = convective_htc_from_nusselt(
        d_p=1.0e-3,
        u_rel=1.0,
        rho_g=0.35,
        mu_g=4.0e-5,
        cp_g=1200.0,
        lambda_g=0.08,
    )
    print(f"  Nu = {nu:.2f}  [>2], h_conv = {h:.1f} W/(m²·K)  [10-5000]", end=" ")
    assert nu > 2.0
    assert 10.0 <= h <= 5000.0
    print("PASS")


def check_daem_radial() -> None:
    """Eq.4.12 径向积分输出应在 [0,1]"""
    from src.thermal.devolatilization import daem_conversion_radial

    r = np.linspace(0.0, 1.0e-3, 8)
    t = np.linspace(0.0, 50.0, 120)
    t_rt = np.vstack([700.0 + 300.0 * (ri / r[-1]) + 120.0 * (t / t[-1]) for ri in r])
    x = daem_conversion_radial(t_rt, t, r)
    print(f"  X_VM,radial = {x*100:.1f}%  [0-100%]", end=" ")
    assert 0.0 <= x <= 1.0
    print("PASS")


def check_pyrolysis_elemental_allocator() -> None:
    """热解元素分配不超支"""
    from src.core.cell import Cell, SolidProps
    from src.core.species import TAR_SURROGATE_FORMULA, get_tar_component_mapping

    c = Cell(solid=SolidProps())
    c.fuel_type = "coal"
    n_c_in, n_h_in, n_o_in = 10.0, 12.0, 4.0
    prod = c._allocate_pyrolysis_products_elemental(nC=n_c_in, nH=n_h_in, nO=n_o_in)

    tar_map = get_tar_component_mapping("coal")
    n_c_out = prod["CO"] + prod["CO2"] + prod["CH4"]
    n_h_out = 2.0 * prod["H2"] + 2.0 * prod["H2O"] + 4.0 * prod["CH4"]
    for tar_label in ("TAR1", "TAR2"):
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map[tar_label]]
        n_c_out += c_tar * prod[tar_label]
        n_h_out += h_tar * prod[tar_label]
    n_o_out = prod["CO"] + 2.0 * prod["CO2"] + prod["H2O"]
    print(f"  C/H/O out = ({n_c_out:.3f},{n_h_out:.3f},{n_o_out:.3f}) <= in ({n_c_in:.3f},{n_h_in:.3f},{n_o_in:.3f})", end=" ")
    assert n_c_out <= n_c_in + 1e-9
    assert n_h_out <= n_h_in + 1e-9
    assert n_o_out <= n_o_in + 1e-9
    print("PASS")


def check_drying_pyrolysis_source_hook() -> None:
    """Cell 干燥/热解源项接线"""
    from src.core.cell import Cell, S_MOISTURE, S_VM, SolidProps

    c = Cell(solid=SolidProps())
    c.T = 1173.15
    c.P = 2.5e6
    c.m_solid_zu[0, S_VM] = 0.05
    c.m_solid_zu[0, S_MOISTURE] = 0.02
    src = c._calc_drying_pyrolysis_gas_source(tau=5.0)

    src_sum = float(np.sum(np.maximum(src, 0.0)))
    rs_sum = float(np.sum(c.R_solid))
    print(f"  gas_source_sum = {src_sum:.3e} mol/s, R_solid_sum = {rs_sum:.3e} kg/s", end=" ")
    assert src_sum > 0.0
    assert rs_sum < 0.0
    print("PASS")


def check_hamel_psi_b_eq315() -> None:
    """Eq.3.15: psi_b = 0.17 * u_mf^-0.33"""
    from src.physics.bubble_dynamics import bubble_interaction_factor

    u_mf = 0.05
    psi = bubble_interaction_factor(u_mf, strategy="wein_1992")
    psi_ref = 0.17 * u_mf ** (-0.33)
    print(f"  psi_b = {psi:.4f}, ref = {psi_ref:.4f}", end=" ")
    assert abs(psi - psi_ref) <= 1e-12
    print("PASS")


def check_hamel_lambda_eq342_trend() -> None:
    """Eq.3.42: lambda_b ~ (P/P0)^-0.7，随压升高应降低"""
    from src.physics.bubble_dynamics import bubble_lifetime

    u_mf = 0.05
    lam_low = bubble_lifetime(d_b=0.1, u_b=1.0, P=101300.0, strategy="hamel_280", u_mf=u_mf)
    lam_high = bubble_lifetime(d_b=0.1, u_b=1.0, P=2.5e6, strategy="hamel_280", u_mf=u_mf)
    ratio = lam_high / max(lam_low, 1e-12)
    print(f"  lambda_low={lam_low:.4e}, lambda_high={lam_high:.4e}, ratio={ratio:.3f} [<1]", end=" ")
    assert lam_high < lam_low
    assert ratio < 0.2
    print("PASS")


def check_visible_epsb_physical_range() -> None:
    """visible epsilon_b（Eq.3.24）在当前主线 clip 口径内"""
    from src.physics.phase_fractions import calc_visible_bubble_fraction

    eps_b = calc_visible_bubble_fraction(u0=0.30, u_b=1.20, u_d=0.12)
    print(f"  eps_b_visible = {eps_b:.4f} [0.05-0.70]", end=" ")
    assert 0.05 <= eps_b <= 0.70
    print("PASS")


def check_freeboard_analytical_default_and_metadata() -> None:
    """默认 freeboard 轨迹应为 analytical_wirsum，且诊断字段存在"""
    from src.core.freeboard_segment import simulate_freeboard
    from src.core.reactor import GAS_SPECIES, GAS_SPECIES_INDEX, ReactorConfig
    from src.core.species import configure_tar_components_by_fuel

    cfg = ReactorConfig()
    print(f"  default freeboard_trajectory_model = {cfg.freeboard_trajectory_model}", end=" ")
    assert cfg.freeboard_trajectory_model == "analytical_wirsum"
    print("PASS")

    configure_tar_components_by_fuel("coal")
    n = np.zeros(len(GAS_SPECIES), dtype=np.float64)
    idx = GAS_SPECIES_INDEX
    n[idx["CO"]] = 1.0
    n[idx["H2"]] = 1.0
    n[idx["H2O"]] = 1.0
    n[idx["N2"]] = 3.0
    out = simulate_freeboard(
        N_in=n,
        T_in=1100.0,
        P=2.5e6,
        D_bed=0.6,
        H_freeboard=2.0,
        n_cells=2,
        u_b_bed_top=1.2,
        d_b_bed_top=0.12,
        eps_b_bed_top=0.3,
        eps_d_void_bed_top=0.5,
        rho_solid_bed_top=1700.0,
        d_p_classes_bed_top=np.array([7e-4]),
        m_char_classes_bed_top=np.array([0.02]),
        m_ash_classes_bed_top=np.array([0.01]),
        fuel_type="coal",
        trajectory_model="analytical_wirsum",
        enabled_reactions=("R5", "R6", "R7", "R8", "R12"),
    )
    print(f"  solver = {out['trajectory_solver']}, diag = {out['trajectory_diag']}", end=" ")
    assert out["trajectory_solver"] == "wirsum_analytical"
    diag = out["trajectory_diag"]
    assert {
        "delta_pos",
        "delta_zero",
        "delta_neg",
        "disc_pos",
        "disc_zero",
        "disc_neg",
        "sign_switch",
        "fallbacks",
        "returns",
    } <= set(diag.keys())
    assert (diag["delta_pos"] + diag["delta_zero"] + diag["delta_neg"]) > 0
    print("PASS")


def check_reaction_numbering_dual_labels() -> None:
    """诊断层必须能同时给出 impl/thesis 编号映射。"""
    from src.core.reaction_numbering import map_impl_diag_to_thesis, reaction_numbering_metadata

    meta = reaction_numbering_metadata()
    assert meta["authority"] == "hamel_1999_thesis"
    assert meta["impl_to_thesis"]["R12"] == "R6"
    assert meta["impl_to_thesis"]["R7"] == "R9"

    mapped = map_impl_diag_to_thesis({"R11b": 1.0, "R11d": 2.0, "R12": 3.0, "R9": 4.0})
    print(f"  mapped keys = {sorted(mapped.keys())}", end=" ")
    assert mapped["R11"] == 3.0
    assert mapped["R6"] == 3.0
    assert mapped["_impl_only:R9"] == 4.0
    print("PASS")


def check_bulk_solid_ratio_definition() -> None:
    """床层总体固含率定义应为 (1-eps_b)*(1-eps_d_void)。"""
    from src.physics.phase_fractions import calc_bulk_solid_holdup

    eps_b = 0.28
    eps_d_void = 0.52
    bulk = calc_bulk_solid_holdup(eps_b, eps_d_void)
    ref = (1.0 - eps_b) * (1.0 - eps_d_void)
    print(f"  bulk_solid = {bulk:.4f}, ref = {ref:.4f}", end=" ")
    assert abs(bulk - ref) < 1e-12
    print("PASS")


def main() -> int:
    checks = [
        ("1. u_mf", check_u_mf),
        ("2. d_b(H_bed)", check_d_b),
        ("3. K_bd (Hamel chain)", check_kbd_hamel_chain),
        ("4. R3 Boudouard (impl r4)", check_r_boudouard),
        ("5. R8 WGSR direction", check_r8_direction),
        ("6. DAEM", check_daem),
        ("7. h_v corrected (Eq.4.4)", check_hv_corrected),
        ("8. Nu & h_conv (Eq.4.9)", check_nusselt_htc),
        ("9. DAEM radial (Eq.4.12)", check_daem_radial),
        ("10. Pyrolysis elemental allocator", check_pyrolysis_elemental_allocator),
        ("11. Drying-pyrolysis source hook", check_drying_pyrolysis_source_hook),
        ("12. psi_b Eq.3.15", check_hamel_psi_b_eq315),
        ("13. lambda_b Eq.3.42 trend", check_hamel_lambda_eq342_trend),
        ("14. visible epsilon_b range", check_visible_epsb_physical_range),
        ("15. freeboard analytical default/meta", check_freeboard_analytical_default_and_metadata),
        ("16. reaction numbering dual labels", check_reaction_numbering_dual_labels),
        ("17. bulk solid ratio definition", check_bulk_solid_ratio_definition),
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
