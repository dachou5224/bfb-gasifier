"""P3：Wirsum 1997 exact 自由板核（抛射 / Roof 3.59 / 解析轨迹）。

不依赖 Reactor 或生产 ``simulate_freeboard``。
"""

from __future__ import annotations

import numpy as np
import pytest

from src.physics.wirsum_1997_freeboard import (
    FREEBOARD_MODEL_ENGINEERING,
    FREEBOARD_MODEL_HAMEL_EXTENDED,
    FREEBOARD_MODEL_LEGACY_ALIAS,
    FREEBOARD_MODEL_WIRSUM_1997_EXACT,
    MARTENS_U_P0_OVER_U_B,
    WIRSUM_1997_U_P0_OVER_U_B,
    WIRSUM_1997_VELOCITY_SIGMA,
    WIRSUM_1997_Z_DRAG_DEFAULT,
    WIRSUM_1997_Z_DRAG_MAX,
    WIRSUM_1997_Z_DRAG_MIN,
    build_launch_velocity_samples,
    calc_roof_mass_flux_eq359,
    calc_u_p0_martens,
    calc_u_p0_wirsum_1997,
    calc_wake_mass_flux_engineering,
    classify_particle_fate,
    classify_particle_fate_branching,
    eject_size_class_mass_flows,
    hamel_eq386_asymptotic_z,
    integrate_quadratic_analytical,
    integrate_quadratic_rk4,
    particle_acceleration_cw,
    resolve_freeboard_model,
    wirsum_1997_abc,
    wirsum_1997_delta,
    wirsum_1997_drag_z,
)


def test_resolve_freeboard_model_splits_exact_from_engineering() -> None:
    exact = resolve_freeboard_model(FREEBOARD_MODEL_WIRSUM_1997_EXACT)
    eng = resolve_freeboard_model(FREEBOARD_MODEL_ENGINEERING)
    hamel = resolve_freeboard_model(FREEBOARD_MODEL_HAMEL_EXTENDED)
    alias = resolve_freeboard_model(FREEBOARD_MODEL_LEGACY_ALIAS)

    assert exact.u_p0_over_u_b == WIRSUM_1997_U_P0_OVER_U_B
    assert exact.ejection == "roof_eq359"
    assert exact.drag == "wirsum_1997_cw"
    assert eng.u_p0_over_u_b == MARTENS_U_P0_OVER_U_B
    assert eng.ejection == "roof_plus_wake_renormalize"
    assert hamel.ejection == "roof_plus_wake"
    assert alias.name == FREEBOARD_MODEL_ENGINEERING
    with pytest.raises(ValueError, match="freeboard_model"):
        resolve_freeboard_model("analytical_wirsum_exact")


def test_launch_mean_is_1p2_ub_and_weights_sum_to_one() -> None:
    u_b = 1.0  # [m/s]
    u_p0 = calc_u_p0_wirsum_1997(u_b)
    assert u_p0 == pytest.approx(1.2)
    assert calc_u_p0_martens(u_b) == pytest.approx(1.53)
    assert u_p0 != calc_u_p0_martens(u_b)

    velocities, weights = build_launch_velocity_samples(
        u_p0, WIRSUM_1997_VELOCITY_SIGMA, n_bins=21
    )
    assert velocities.shape == weights.shape == (21,)
    assert velocities.min() == pytest.approx(0.0)
    assert float(np.sum(weights)) == pytest.approx(1.0)
    # 1.2 − 2·0.6 = 0，高斯对称，质量加权均值就是 u_p0
    assert float(np.dot(weights, velocities)) == pytest.approx(u_p0, abs=1e-12)


def test_launch_distribution_is_symmetric_gaussian_in_xi() -> None:
    velocities, weights = build_launch_velocity_samples(2.0, 0.5, n_bins=11)
    assert float(np.sum(weights)) == pytest.approx(1.0)
    assert weights[0] == pytest.approx(weights[-1])
    assert velocities[5] == pytest.approx(2.0)


def test_roof_eq359_hand_calculated_flux() -> None:
    # F = 1e-3 · 1400 · 0.5 · (3 · 0.2 · 0.8 / 0.1) = 3.36 kg/(m² s)
    flux = calc_roof_mass_flux_eq359(
        d_p=1.0e-3,
        rho_p=1400.0,
        eps_d=0.50,
        eps_b=0.20,
        u_b=0.80,
        d_b=0.10,
    )
    assert flux == pytest.approx(3.36, rel=0.0, abs=1e-12)


def test_exact_ejection_is_roof_only_and_not_renormalized() -> None:
    d_p = np.array([5.0e-4, 1.0e-3], dtype=np.float64)  # [m]
    feed = np.array([0.2, 0.8], dtype=np.float64)  # [kg/s]
    area = 0.25  # [m²]
    kwargs = dict(
        d_p_classes=d_p,
        feed_classes=feed,
        rho_p=1400.0,
        eps_d=0.50,
        eps_b=0.20,
        u_b=0.80,
        d_b=0.10,
        area=area,
    )
    exact = eject_size_class_mass_flows(mode=FREEBOARD_MODEL_WIRSUM_1997_EXACT, **kwargs)
    hamel = eject_size_class_mass_flows(mode=FREEBOARD_MODEL_HAMEL_EXTENDED, **kwargs)
    eng = eject_size_class_mass_flows(mode="engineering_conservative", **kwargs)

    d_p_mean = 0.2 * 5.0e-4 + 0.8 * 1.0e-3
    roof = calc_roof_mass_flux_eq359(
        d_p=d_p_mean, rho_p=1400.0, eps_d=0.50, eps_b=0.20, u_b=0.80, d_b=0.10
    )
    wake = calc_wake_mass_flux_engineering(
        rho_p=1400.0, eps_d=0.50, eps_b=0.20, u_b=0.80
    )
    assert float(np.sum(exact)) == pytest.approx(roof * area)
    assert wake > 0.0
    assert float(np.sum(hamel)) == pytest.approx((roof + wake) * area)
    assert float(np.sum(eng)) == pytest.approx(float(np.sum(feed)))
    assert float(np.sum(exact)) != pytest.approx(float(np.sum(feed)))
    np.testing.assert_allclose(exact / np.sum(exact), feed / np.sum(feed))


def test_abc_matches_force_law_on_both_velocity_branches() -> None:
    common = dict(u_g=1.2, rho_g=0.35, rho_p=1700.0, d_p=7.0e-4, mu_g=2.0e-5, z=0.44)
    for u_p in (2.4, 0.6):
        a, b, c, branch = wirsum_1997_abc(u_p=u_p, **common)
        if u_p > 1.2:
            assert branch == "faster"
        else:
            assert branch == "slower"
        poly = a * u_p * u_p + b * u_p + c
        force = particle_acceleration_cw(u_p, **common)
        assert poly == pytest.approx(force, rel=0.0, abs=1e-12)


def test_delta_sign_branches_are_both_reachable() -> None:
    # Δ>0：构造 A=1, B=0, C=1 → 4AC−B²=4
    assert wirsum_1997_delta(1.0, 0.0, 1.0) > 0.0
    # Δ<0：A=0.5, B=0, C=-0.5 → 4AC−B²=-1
    assert wirsum_1997_delta(0.5, 0.0, -0.5) < 0.0

    faster = wirsum_1997_abc(
        u_g=1.0,
        u_p=3.0,
        rho_g=0.35,
        rho_p=1700.0,
        d_p=2.0e-3,
        mu_g=2.0e-5,
        z=0.8,
    )
    slower = wirsum_1997_abc(
        u_g=1.5,
        u_p=0.4,
        rho_g=0.80,
        rho_p=1400.0,
        d_p=2.0e-4,
        mu_g=3.0e-5,
        z=0.44,
    )
    assert faster[3] == "faster"
    assert slower[3] == "slower"
    # 物理点上至少覆盖一个非零判别式；符号由 4AC−B² 决定，不由 u_p≷u_g 绑定
    assert wirsum_1997_delta(*faster[:3]) != 0.0
    assert wirsum_1997_delta(*slower[:3]) != 0.0


def test_analytical_trajectory_matches_rk4_on_frozen_abc() -> None:
    cases = (
        (1.2, 1.0, 0.0, 1.0, 0.02),   # Δ>0
        (1.2, 0.5, 0.0, -0.5, 0.02),  # Δ<0
        (2.4, *wirsum_1997_abc(
            u_g=1.2, u_p=2.4, rho_g=0.35, rho_p=1700.0, d_p=7e-4, mu_g=2e-5, z=0.44
        )[:3], 0.04),
        (0.6, *wirsum_1997_abc(
            u_g=1.2, u_p=0.6, rho_g=0.35, rho_p=1700.0, d_p=7e-4, mu_g=2e-5, z=0.44
        )[:3], 0.04),
    )
    errors: list[dict[str, float]] = []
    for u0, a, b, c, t_end in cases:
        t_an, u_an, z_an = integrate_quadratic_analytical(u0, a, b, c, t_end, n_steps=40)
        t_rk, u_rk, z_rk = integrate_quadratic_rk4(u0, a, b, c, t_end, n_steps=4000)
        u_rk_i = np.interp(t_an, t_rk, u_rk)
        z_rk_i = np.interp(t_an, t_rk, z_rk)
        err_u = float(np.max(np.abs(u_an - u_rk_i)))
        err_z = float(np.max(np.abs(z_an - z_rk_i)))
        errors.append(
            {
                "delta": wirsum_1997_delta(a, b, c),
                "err_u_m_per_s": err_u,
                "err_z_m": err_z,
            }
        )
        assert err_u < 2.0e-6
        assert err_z < 2.0e-7
    # 量化误差写入断言消息，便于 A/B 对照
    assert errors[0]["delta"] > 0.0
    assert errors[1]["delta"] < 0.0


def test_particle_reaches_next_cell_or_returns() -> None:
    # 几乎无阻力：高初速、矮格子 → 穿出；低初速、高格子 → 回落
    a, b, c = 0.0, 0.0, -9.81
    assert classify_particle_fate(5.0, a, b, c, cell_height=0.20) == "reaches_next_cell"
    assert classify_particle_fate(0.40, a, b, c, cell_height=2.00) == "returns"

    a_f, b_f, c_f, _ = wirsum_1997_abc(
        u_g=0.8, u_p=4.0, rho_g=0.35, rho_p=1700.0, d_p=5e-4, mu_g=2e-5, z=0.44
    )
    a_s, b_s, c_s, _ = wirsum_1997_abc(
        u_g=1.2, u_p=0.5, rho_g=0.35, rho_p=2500.0, d_p=2e-3, mu_g=2e-5, z=0.80
    )
    assert classify_particle_fate(4.0, a_f, b_f, c_f, cell_height=0.05) == "reaches_next_cell"
    assert classify_particle_fate(0.5, a_s, b_s, c_s, cell_height=1.50) == "returns"


def test_default_z_is_inside_wirsum_range() -> None:
    z = wirsum_1997_drag_z()
    assert WIRSUM_1997_Z_DRAG_MIN <= z <= WIRSUM_1997_Z_DRAG_MAX
    assert z == WIRSUM_1997_Z_DRAG_DEFAULT
    with pytest.raises(ValueError, match="drag z"):
        wirsum_1997_drag_z(0.0)


def test_hamel_eq386_z_depends_on_phi() -> None:
    z08 = hamel_eq386_asymptotic_z(0.8)
    z10 = hamel_eq386_asymptotic_z(1.0)
    z03 = hamel_eq386_asymptotic_z(0.3)
    assert z08 == pytest.approx(73.69 * np.exp(-5.0784 * 0.8), rel=0.0, abs=1e-12)
    assert z08 > 1.2
    assert z10 < z08 < z03
    hamel = resolve_freeboard_model(FREEBOARD_MODEL_HAMEL_EXTENDED)
    assert hamel.drag == "hamel_z_exp"
    with pytest.raises(ValueError, match="phi_s"):
        hamel_eq386_asymptotic_z(0.2)


def test_branch_acceleration_matches_force_across_ug() -> None:
    common = dict(u_g=1.2, rho_g=0.35, rho_p=1700.0, d_p=7.0e-4, mu_g=2.0e-5, z=1.27)
    for u_p in (1.2 + 1e-6, 1.2 - 1e-6):
        a, b, c, _ = wirsum_1997_abc(u_p=u_p, **common)
        poly = a * u_p * u_p + b * u_p + c
        force = particle_acceleration_cw(u_p, **common)
        assert poly == pytest.approx(force, rel=0.0, abs=1e-10)


def test_return_event_is_bed_surface_not_apex() -> None:
    # 无阻力竖直上抛：顶点 u=0 时 z>0，必须继续落到 z=0
    fate = classify_particle_fate_branching(
        0.40,
        u_g=0.0,
        rho_g=1.0,
        rho_p=2000.0,
        d_p=2.0e-3,
        mu_g=2.0e-5,
        z=0.44,
        cell_height=2.0,
        dt=5.0e-4,
        t_max=2.0,
    )
    assert fate == "returns"
    high = classify_particle_fate_branching(
        8.0,
        u_g=2.0,
        rho_g=0.3,
        rho_p=800.0,
        d_p=5.0e-5,
        mu_g=2.0e-5,
        z=0.44,
        cell_height=0.05,
        dt=5.0e-4,
        t_max=1.0,
    )
    assert high == "reaches_next_cell"


def test_gaussian_weights_are_named_truncation() -> None:
    from src.physics.wirsum_1997_freeboard import GAUSSIAN_SAMPLE_MODE

    assert GAUSSIAN_SAMPLE_MODE == "truncated_pm_2sigma_uniform_xi"
    u5, w5 = build_launch_velocity_samples(1.2, 0.6, 5)
    u21, w21 = build_launch_velocity_samples(1.2, 0.6, 21)
    assert w5.sum() == pytest.approx(1.0)
    assert w21.sum() == pytest.approx(1.0)
    mean5 = float(np.dot(u5, w5))
    mean21 = float(np.dot(u21, w21))
    assert abs(mean5 - mean21) < 0.05
    assert WIRSUM_1997_VELOCITY_SIGMA == pytest.approx(0.6)
