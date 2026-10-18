# -*- coding: utf-8 -*-
"""§3.1 与床内固体输运。期望值按印刷式手算，不用被测函数回代。"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from hamel.hydrodynamics.bed import (
    F_W,
    G,
    N_B,
    P_0,
    BedConditions,
    HydrodynamicsError,
    alpha_b,
    archimedes_number,
    bed_profile,
    bed_upflow,
    bubble_drag_coefficient,
    bubble_growth_rate,
    bubble_lifetime,
    bubble_rise_velocity,
    coalescence_parameter,
    ergun_archimedes,
    exchange_coefficient,
    growth_start_height,
    initial_bubble_diameter,
    isolated_bubble_velocity,
    minimum_fluidization_velocity,
    orifice_volumetric_flow,
    particle_density_from_bed_mass,
    psi_b,
    richardson_zaki_exponent,
    settling_reynolds,
    sphericity_diameter_product,
    state_at_height,
    suspension_dynamic_viscosity,
    suspension_kinematic_viscosity,
    suspension_superficial_velocity,
    suspension_voidage,
    terminal_velocity,
    throughflow_velocity,
    upsilon_b,
    visible_bubble_fraction,
    wake_solids_upflow,
)


def test_minimum_fluidization_hand_calculation():
    # ε=0.4，ν=1.5e-5，φd=7e-4，ρ_P−ρ_g=1999，ρ_g=1，g=9.81
    # inside = 3.1e-4 * 0.4**3 * 9.81 * (7e-4)**3 * 1999 / (0.6**2 * (1.5e-5)**2 * 1)
    inside = 3.1e-4 * 0.4**3 * G * (7e-4) ** 3 * 1999.0 / (0.6**2 * (1.5e-5) ** 2)
    expected = 42.9 * 0.6 * 1.5e-5 / 7e-4 * (math.sqrt(1.0 + inside) - 1.0)
    got = minimum_fluidization_velocity(0.4, 1.5e-5, 0.7, 1e-3, 2000.0, 1.0)
    assert got == pytest.approx(expected)
    assert got == pytest.approx(0.3459027656268236)


def test_archimedes_and_ergun():
    # Ar = d³ g (ρ_P−ρ_g) / (ν² ρ_g)
    expected = (1e-3) ** 3 * G * 1999.0 / ((1.5e-5) ** 2 * 1.0)
    assert archimedes_number(1e-3, 1.5e-5, 2000.0, 1.0) == pytest.approx(expected)
    # 式 (3.1) 右端：150(1−ε)/ε³ Re + 1.75/ε³ Re²，Re=2.5，ε=0.4
    reynolds = 2.5
    voidage = 0.4
    ergun = 150.0 * (1.0 - voidage) / voidage**3 * reynolds + 1.75 / voidage**3 * reynolds**2
    assert ergun_archimedes(0.4, 1.0, 2.5) == pytest.approx(ergun)


def test_sphericity_diameter_product():
    # 式 (3.3) 的正根。u=0.2，μ=1.8e-5，ε=0.4，ρ_g=1，ρ_P=2000
    # buoyancy = ε³ (ρ_P−ρ_g) g
    # linear = 1.75 ρ_g u² / (2 buoyancy)
    # under = 150 (1−ε) μ u / buoyancy
    # φd = linear + sqrt(under + linear²)
    buoyancy = 0.4**3 * 1999.0 * G
    linear = 1.75 * 1.0 * 0.2**2 / (2.0 * buoyancy)
    under = 150.0 * 0.6 * 1.8e-5 * 0.2 / buoyancy
    expected = linear + math.sqrt(under + linear**2)
    got = sphericity_diameter_product(0.2, 0.4, 1.8e-5, 1.0, 2000.0)
    assert got == pytest.approx(expected)
    assert got == pytest.approx(0.0005367431626620126)


def test_suspension_velocity_allows_d_and_rejects_other_beds():
    # 式 (3.11)：u_d = u_mf + (u_0−u_mf)/3
    assert suspension_superficial_velocity(0.2, 0.5, "technical", "D", True) == pytest.approx(
        0.2 + 0.3 / 3.0
    )
    with pytest.raises(HydrodynamicsError):
        suspension_superficial_velocity(0.2, 0.5, "porous", "B", True)
    with pytest.raises(HydrodynamicsError):
        suspension_superficial_velocity(0.2, 0.5, "technical", "B", False)
    with pytest.raises(HydrodynamicsError):
        suspension_superficial_velocity(0.2, 0.5, "technical", "A", True)
    with pytest.raises(HydrodynamicsError):
        suspension_superficial_velocity(0.2, 0.2, "technical", "B", True)


def test_psi_and_upsilon_only_for_b():
    assert psi_b("technical", "B") == pytest.approx(0.76)
    assert psi_b("porous", "B") == pytest.approx(0.67)
    assert upsilon_b(0.05, "B") == pytest.approx(0.64)
    assert upsilon_b(0.15, "B") == pytest.approx(2.0 * math.sqrt(0.15))
    assert upsilon_b(1.2, "B") == pytest.approx(2.0)
    with pytest.raises(HydrodynamicsError):
        psi_b("technical", "D")
    with pytest.raises(HydrodynamicsError):
        upsilon_b(0.15, "D")


def test_throughflow_and_visible_fraction():
    # P=P_0 时 u_b,r = 2.7 u_d
    assert throughflow_velocity(0.4, P_0) == pytest.approx(N_B * 0.4)
    # 0.5 MPa：2.7 * 0.4 * (5e5/101300)^(−0.15)
    expected = 2.7 * 0.4 * (5.0e5 / P_0) ** (-0.15)
    assert throughflow_velocity(0.4, 5.0e5) == pytest.approx(expected)
    # 式 (3.7)：u_0 = ε_b u_b + ε_b u_b,r + (1−ε_b) u_d
    u_0, u_d, u_b, u_br = 0.9, 0.4, 1.2, 0.5
    fraction = visible_bubble_fraction(u_0, u_d, u_b, u_br)
    assert u_0 == pytest.approx(fraction * u_b + fraction * u_br + (1.0 - fraction) * u_d)
    # 常压分母用 n_b u_d，与式 (3.24) 的印刷形式相同
    printed = (0.8 - 0.3) / (1.1 + N_B * 0.3 - 0.3)
    assert visible_bubble_fraction(0.8, 0.3, 1.1, N_B * 0.3) == pytest.approx(printed)


def test_alpha_is_the_adopted_reading_only():
    # 采用：u_b ε_d / u_d = 0.4*0.5/0.3
    assert alpha_b(0.4, 0.3, 0.5) == pytest.approx(0.4 * 0.5 / 0.3)
    # 印刷式 (u_b/u_d)ε_b 与定义句 u_b/u_d 都不是这个数
    assert alpha_b(0.4, 0.3, 0.5) != pytest.approx((0.4 / 0.3) * 0.2)
    assert alpha_b(0.4, 0.3, 0.5) != pytest.approx(0.4 / 0.3)
    assert coalescence_parameter(0.5) == pytest.approx(1.0 - 0.5**3)
    assert coalescence_parameter(1.2) == pytest.approx(0.0)
    with pytest.raises(HydrodynamicsError):
        coalescence_parameter(-0.1)


def test_suspension_voidage_and_settling():
    # Re_S=0.2 → 4.65；0.5 → 4.4*0.5^(−0.03)；10 → 4.4*10^(−0.1)；800 → 2.4
    assert richardson_zaki_exponent(0.2) == pytest.approx(4.65)
    assert richardson_zaki_exponent(0.5) == pytest.approx(4.4 * 0.5 ** (-0.03))
    assert richardson_zaki_exponent(10.0) == pytest.approx(4.4 * 10.0 ** (-0.1))
    assert richardson_zaki_exponent(800.0) == pytest.approx(2.4)
    # ε_d = 0.46 * (0.3/0.2)^(1/4)
    assert suspension_voidage(0.46, 0.3, 0.2, 4.0) == pytest.approx(0.46 * (0.3 / 0.2) ** 0.25)
    # 式 (3.31)：d=1e-3，ν=1.5e-5，ρ_P=2000，ρ_g=1
    archimedes = (1e-3) ** 3 * G * 1999.0 / ((1.5e-5) ** 2)
    factor = math.sqrt(1.0 + archimedes**0.5 / 9.0) - 1.0
    expected = 18.0 * (1.5e-5 / 1e-3) * factor**2
    terminal = terminal_velocity(1e-3, 1.5e-5, 2000.0, 1.0)
    assert terminal == pytest.approx(expected)
    # 式 (3.28) 的 Re 就是式 (3.26) 的 Re_S
    reynolds = settling_reynolds(terminal, 1e-3, 1.0, 1.5e-5)
    assert reynolds == pytest.approx(terminal * 1e-3 * 1.0 / 1.5e-5)


def test_suspension_viscosity():
    fine = 1.5e-4
    coarse = 5e-4
    assert suspension_dynamic_viscosity(fine) == pytest.approx(232.0 * fine ** (2.0 / 3.0))
    assert suspension_dynamic_viscosity(coarse) == pytest.approx(2.5e5 * coarse**1.5)
    assert suspension_kinematic_viscosity(0.02, 1.2, 0.5) == pytest.approx(0.02 / (1.2 * 0.5))
    with pytest.raises(HydrodynamicsError):
        suspension_kinematic_viscosity(0.02, 1.2, 1.0)


def test_bubble_rise_fixed_point():
    # 先取 u_b=0.8、d_b=0.05、ν_d=1e-4，反解 excess，再让求解器回到 0.8
    velocity = 0.8
    diameter = 0.05
    kinematic = 1e-4
    upsilon = 0.64
    reynolds = velocity * diameter / kinematic
    drag = 16.0 / reynolds + 2.64
    isolated = math.sqrt((1.0 / drag) * (4.0 / 3.0) * diameter * G)
    excess = velocity - upsilon * isolated
    assert bubble_drag_coefficient(reynolds) == pytest.approx(drag)
    assert isolated_bubble_velocity(diameter, drag) == pytest.approx(isolated)
    assert bubble_rise_velocity(excess, upsilon, diameter, kinematic) == pytest.approx(velocity)


def test_lifetime_growth_and_exchange():
    # λ_b(P_0) = 280 * 0.2 / 9.81；加压再乘 (P/P_0)^(−0.7)
    assert bubble_lifetime(0.2, P_0) == pytest.approx(280.0 * 0.2 / G)
    assert bubble_lifetime(0.2, 5.0e5) == pytest.approx(280.0 * 0.2 / G * (5.0e5 / P_0) ** (-0.7))
    # 式 (3.50)，S. 35：1.5 u_b,r/d_b + sqrt(144 D_g ε_mf u_b / (π d_b³))
    convective = 1.5 * 0.8 / 0.2
    diffusive = math.sqrt(144.0 * 2e-5 * 0.46 * 0.9 / (math.pi * 0.2**3))
    assert exchange_coefficient(0.8, 0.2, 2e-5, 0.46, 0.9) == pytest.approx(convective + diffusive)
    # 式 (3.41)，P=P_0：聚并项减分裂项。分母 1 − ξ (6/π)^(1/3) ε^(1/3)
    denominator = 1.0 - 0.5 * (6.0 / math.pi) ** (1.0 / 3.0) * 0.2 ** (1.0 / 3.0)
    coalescence = (2.0 / (9.0 * math.pi)) ** (1.0 / 3.0) * (0.2 ** (1.0 / 3.0) / denominator)
    splitting = 0.1 / (3.0 * 1.0 * 0.8)
    assert bubble_growth_rate(0.1, 0.2, 0.5, 1.0, 0.8, P_0) == pytest.approx(coalescence - splitting)


def test_distributor_diameters():
    # V̇_0 = 0.5 * 0.02 / 10 = 0.001；d_b0 = 1.3 (V̇²/g)^0.2
    flow = 0.5 * 0.02 / 10.0
    expected = 1.3 * (flow**2 / G) ** 0.2
    assert orifice_volumetric_flow(0.5, 0.02, 10) == pytest.approx(flow)
    assert initial_bubble_diameter("technical", 0.5, 0.02, 10, 0.2) == pytest.approx(expected)
    # 多孔板：2.78/g (u_0−u_mf)²
    assert initial_bubble_diameter("porous", 0.5, 0.02, 10, 0.2) == pytest.approx(
        2.78 / G * (0.3) ** 2
    )
    assert growth_start_height("technical", 0.04, 0.01) == pytest.approx(0.01 + 0.02)
    assert growth_start_height("porous", 0.04, 0.01) == pytest.approx(0.02)
    with pytest.raises(HydrodynamicsError):
        orifice_volumetric_flow(0.5, 0.02, 10.5)


def test_atmospheric_state_matches_printed_bubble_fraction():
    conditions = _conditions(pressure=P_0, rho_p=2000.0)
    state = state_at_height(0.4, 0.55, 0.2, 0.02, 40, 0.0, 40, conditions)
    printed = (0.55 - state.u_d) / (state.u_b + N_B * state.u_d - state.u_d)
    assert state.u_br == pytest.approx(N_B * state.u_d)
    assert state.epsilon_b == pytest.approx(printed)
    # 采用的 α_b 不是 (u_b/u_d)ε_b
    assert state.alpha == pytest.approx(state.u_b * state.epsilon_d / state.u_d)


def test_pressurized_fraction_uses_corrected_throughflow():
    conditions = _conditions(pressure=5.0e5, rho_p=2000.0)
    state = state_at_height(0.4, 0.8, 0.2, 0.02, 40, 0.0, 60, conditions)
    balance = (0.8 - state.u_d) / (state.u_b + state.u_br - state.u_d)
    printed = (0.8 - state.u_d) / (state.u_b + N_B * state.u_d - state.u_d)
    assert state.u_br != pytest.approx(N_B * state.u_d)
    assert state.epsilon_b == pytest.approx(balance)
    assert state.epsilon_b != pytest.approx(printed)


def test_visible_fraction_vanishes_as_u0_approaches_umf():
    conditions = _conditions(pressure=P_0, rho_p=2000.0)
    # 取 h=0，直径停在 d_b,0，避免长大项在极小过剩气速下把直径积成非正数。
    state = state_at_height(0.0, 0.2001, 0.2, 0.02, 40, 0.0, 20, conditions)
    assert 0.0 < state.epsilon_b < 1.0e-3


def test_pressure_rise_shrinks_bubbles_in_the_stated_window():
    # 物性随压力成比例，过剩气速固定为 0.6 m/s。0.5 MPa 到 2.0 MPa 是叙述成立的区间。
    low = _pressurized_state(0.5e6)
    high = _pressurized_state(2.0e6)
    assert high.d_b < low.d_b
    assert high.u_b > low.u_b


def test_inventory_converges_when_the_mesh_is_refined():
    edges = np.linspace(0.0, 1.2, 9)
    velocities = np.full(8, 0.6)
    coarse = bed_profile(edges, velocities, 0.2, 0.02, 80, 0.0, 8, _conditions(P_0, 900.0))
    fine = bed_profile(edges, velocities, 0.2, 0.02, 80, 0.0, 16, _conditions(P_0, 900.0))
    assert coarse.solids_mass.sum() > 0.0
    assert fine.solids_mass.sum() == pytest.approx(coarse.solids_mass.sum(), rel=1e-4)


def test_bed_surface_is_supplied_by_the_caller():
    internal = wake_solids_upflow([0.2, 0.25], [0.4, 0.5], 800.0, [0.5, 0.5], 0.02)
    surface = 0.01
    full = bed_upflow(surface, internal)
    assert full[0] == pytest.approx(surface)
    assert full[1:] == pytest.approx(internal)
    # 床面处若误用尾涡公式，数值与调用方给定值不同
    mistaken = float(wake_solids_upflow([0.2], [0.4], 800.0, [0.5], 0.02)[0])
    assert full[0] != pytest.approx(mistaken)
    assert F_W == pytest.approx(0.25)


def test_density_inversion_recovers_the_mass():
    edges = np.linspace(0.0, 0.6, 5)
    conditions = _conditions(P_0, 900.0)
    u_mf = minimum_fluidization_velocity(0.45, conditions.nu_g, 0.7, 1.0e-3, 900.0, 1.2)
    velocities = np.full(4, u_mf + 0.4)
    profile = bed_profile(edges, velocities, u_mf, 0.02, 40, 0.0, 4, conditions)
    recovered = particle_density_from_bed_mass(
        float(profile.solids_mass.sum()),
        edges,
        velocities,
        0.02,
        40,
        0.0,
        4,
        0.7,
        conditions,
    )
    assert recovered == pytest.approx(900.0, rel=1e-3)


def test_hamel_has_no_reading_switch():
    text = Path("hamel/hydrodynamics/bed.py").read_text(encoding="utf-8")
    assert "printed_minus_one" not in text
    assert "alpha_reading" not in text


def test_vtt_case_1_smoke_with_temporary_inputs():
    # X5 未裁定，仅作冒烟。n_0、h_A、粒径、空隙率、球形度、温度和气体组成都是临时值。
    # Geldart 取 B：原文没有印出分类边界，这不是裁定。
    result = _vtt_case_1_smoke()
    assert result["u_mf"] < result["u_0"]
    assert result["rho_p"] > result["rho_g"]
    assert np.all(result["epsilon_b"] > 0.0)
    assert np.all(result["epsilon_b"] < 1.0)
    assert np.all(result["solids_mass"] > 0.0)
    assert result["solids_mass"].sum() == pytest.approx(result["bed_mass"], rel=1e-6)
    assert result["internal_wake_upflow"].size == 7
    assert result["bed_surface"] == "待 T23 接入"


def _conditions(pressure, rho_p):
    return BedConditions(
        1.0e-3,
        0.45,
        rho_p,
        1.2,
        2.8e-5,
        2.8e-5 / 1.2,
        pressure,
        2.0e-5,
        0.15,
        "technical",
        "B",
        True,
    )


def _pressurized_state(pressure):
    rho_g = 1.2 * pressure / P_0
    mu_g = 2.8e-5
    conditions = BedConditions(
        0.7e-3,
        0.46,
        2600.0,
        rho_g,
        mu_g,
        mu_g / rho_g,
        pressure,
        2.0e-5 * P_0 / pressure,
        0.15,
        "technical",
        "B",
        True,
    )
    u_mf = minimum_fluidization_velocity(0.46, mu_g / rho_g, 0.8, 0.7e-3, 2600.0, rho_g)
    return state_at_height(0.55, u_mf + 0.6, u_mf, 0.02, 100, 0.0, 40, conditions)


def _vtt_case_1_smoke():
    """工况 1 的临时输入。返回值供测试和 PR 记录，不是工况文件里的裁定。"""
    temperature = 1123.15
    gas_constant = 8.314462618
    air_mol_s = 0.6419598861822472
    steam_kg_s = 0.001888888888888889
    steam_mol_s = steam_kg_s / 0.01801528
    molar_flow = air_mol_s + steam_mol_s
    molar_mass = (air_mol_s * 0.0289647 + steam_mol_s * 0.01801528) / molar_flow
    pressure = 500_000.0
    area = math.pi * 0.075**2
    rho_g = pressure * molar_mass / (gas_constant * temperature)
    mu_g = 1.716e-5 * (temperature / 273.15) ** 1.5 * (273.15 + 110.4) / (temperature + 110.4)
    nu_g = mu_g / rho_g
    diffusivity = 2.0e-5 * (temperature / 273.15) ** 1.5 * (101_325.0 / pressure)
    u_0 = molar_flow * gas_constant * temperature / pressure / area
    diameter = 5.0e-4
    epsilon_mf = 0.45
    phi_s = 0.75
    n_0 = 100
    h_a = 0.01
    bed_mass = 2200.0 * area / G
    edges = np.linspace(0.0, 1.2, 9)
    velocities = np.full(8, u_0)
    seed = BedConditions(
        diameter, epsilon_mf, 100.0, rho_g, mu_g, nu_g, pressure, diffusivity, 0.15, "technical", "B", True
    )
    rho_p = particle_density_from_bed_mass(
        bed_mass, edges, velocities, area, n_0, h_a, 6, phi_s, seed
    )
    u_mf = minimum_fluidization_velocity(epsilon_mf, nu_g, phi_s, diameter, rho_p, rho_g)
    conditions = BedConditions(
        diameter, epsilon_mf, rho_p, rho_g, mu_g, nu_g, pressure, diffusivity, 0.15, "technical", "B", True
    )
    profile = bed_profile(edges, velocities, u_mf, area, n_0, h_a, 6, conditions)
    return {
        "u_0": u_0,
        "u_mf": u_mf,
        "rho_g": rho_g,
        "rho_p": rho_p,
        "bed_mass": bed_mass,
        "epsilon_b": profile.epsilon_b,
        "u_b": profile.u_b,
        "d_b": profile.d_b,
        "epsilon_d": profile.epsilon_d,
        "k_bd": profile.k_bd,
        "solids_mass": profile.solids_mass,
        "internal_wake_upflow": profile.internal_wake_upflow,
        "bed_surface": "待 T23 接入",
    }
