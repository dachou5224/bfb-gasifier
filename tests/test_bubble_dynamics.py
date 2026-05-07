from __future__ import annotations

import numpy as np

from src.core.constants import P0_HAMEL, g
from src.physics.bubble_dynamics import (
    bubble_interaction_factor,
    bubble_diameter_ode,
    bubble_lifetime,
    bubble_rise_velocity,
    classify_bubble_regime,
    integrate_bubble_diameter,
    resolve_xi_b,
)


def test_bubble_lifetime_default_matches_hamel_pressure_form() -> None:
    u_mf = 0.045
    P = 2.5e6
    expected = 280.0 * u_mf / g * (P / P0_HAMEL) ** (-0.7)
    got = bubble_lifetime(0.2, 2.0, P, u_mf=u_mf)
    assert np.isclose(got, expected)


def test_bubble_lifetime_current_matches_legacy_formula() -> None:
    d_b = 0.2
    u_b = 2.0
    P = 2.5e6
    expected = (d_b / (0.5 * u_b)) * (P / P0_HAMEL) ** (-0.2)
    got = bubble_lifetime(d_b, u_b, P, strategy="current")
    assert np.isclose(got, expected)


def test_bubble_lifetime_hamel_280_uses_umf_pressure_form() -> None:
    u_mf = 0.045
    P = 2.5e6
    expected = 280.0 * u_mf / g * (P / P0_HAMEL) ** (-0.7)
    got = bubble_lifetime(0.2, 2.0, P, strategy="hamel_280", u_mf=u_mf)
    assert np.isclose(got, expected)


def test_integrate_bubble_diameter_hamel_lambda_changes_ode_path() -> None:
    kwargs = {
        "u0": 0.3,
        "u_mf": 0.05,
        "P": 2.5e6,
        "H_bed": 5.0,
        "D_bed": 0.6,
        "n_points": 12,
        "method": "hilligardt_ode",
    }
    _, db_current = integrate_bubble_diameter(**kwargs, lambda_strategy="current")
    _, db_hamel = integrate_bubble_diameter(**kwargs, lambda_strategy="hamel_280")

    assert db_current.shape == db_hamel.shape
    assert np.all(db_current > 0.0)
    assert np.all(db_hamel > 0.0)
    assert not np.allclose(db_current, db_hamel)


def test_resolve_xi_b_supports_fixed_and_hamel_regime() -> None:
    assert np.isclose(resolve_xi_b(0.5, strategy="fixed_035"), 0.35)
    assert np.isclose(resolve_xi_b(0.5, strategy="hamel_regime"), 1.0 - 0.5**3)
    assert np.isclose(resolve_xi_b(2.0, strategy="hamel_regime"), 0.0)


def test_classify_bubble_regime_uses_ud_ratio() -> None:
    assert classify_bubble_regime(0.8, 1.0) == "slow"
    assert classify_bubble_regime(1.2, 1.0) == "fast"


def test_bubble_interaction_factor_supports_hamel_sources() -> None:
    assert np.isclose(bubble_interaction_factor(strategy="technical_distributor"), 0.76)
    assert np.isclose(bubble_interaction_factor(strategy="porous_plate"), 0.67)
    assert np.isclose(
        bubble_interaction_factor(0.05, strategy="wein_1992"),
        0.17 * 0.05 ** (-0.33),
    )


def test_bubble_rise_velocity_supports_heinbockel_pressure_form() -> None:
    u0 = 0.30
    u_mf = 0.05
    d_b = 0.10
    P = 2.5e6
    u_bs = 0.711 * np.sqrt(g * d_b)
    hill = bubble_rise_velocity(u0, u_mf, d_b, strategy="hilligardt_eq313")
    hein = bubble_rise_velocity(u0, u_mf, d_b, P=P, strategy="heinbockel_eq343")
    legacy_plus = bubble_rise_velocity(u0, u_mf, d_b, P=P, strategy="heinbockel_eq343_legacy_plus")
    legacy = bubble_rise_velocity(u0, u_mf, d_b, P=P, strategy="heinbockel_eq343_legacy_2p14_p07")
    expected_hamel = 0.76 * (u0 - u_mf) * ((P / P0_HAMEL) ** 0.2 - 1.0) + u_bs
    assert np.isclose(hein, expected_hamel)
    assert legacy_plus > hein
    assert legacy > hein
    assert hein >= u_bs


def test_integrate_bubble_diameter_supports_heinbockel_ode_branch() -> None:
    kwargs = {
        "u0": 0.3,
        "u_mf": 0.05,
        "P": 2.5e6,
        "H_bed": 5.0,
        "D_bed": 0.6,
        "n_points": 12,
        "method": "hilligardt_ode",
        "u_d": 1.45 * 0.05,
    }
    _, db_hill = integrate_bubble_diameter(
        **kwargs,
        velocity_strategy="hilligardt_eq313",
        ode_strategy="hilligardt_eq333",
    )
    _, db_hein = integrate_bubble_diameter(
        **kwargs,
        velocity_strategy="heinbockel_eq343",
        ode_strategy="heinbockel_eq341",
    )
    assert db_hill.shape == db_hein.shape
    assert np.all(db_hill > 0.0)
    assert np.all(db_hein > 0.0)
    assert not np.allclose(db_hill, db_hein)


def test_integrate_bubble_diameter_heinbockel_ode_matches_transcribed_eq341() -> None:
    """heinbockel_eq341 follows docs/hamel_submodels Eq.3.41 transcription."""
    d_b = 0.2
    u0 = 0.3
    u_mf = 0.05
    u_d = 1.45 * u_mf
    P = 2.5e6
    psi_b = 0.76
    u_b = bubble_rise_velocity(u0, u_mf, d_b, psi_b=psi_b, P=P, strategy="heinbockel_eq343")
    eps_b = np.clip((u0 - u_mf) / u_b, 1e-6, 0.95)
    lam_b = bubble_lifetime(d_b, u_b, P, strategy="hamel_280", u_mf=u_mf)
    eps_b_13 = eps_b ** (1.0 / 3.0)
    pressure_ratio = P / P0_HAMEL
    denom = 1.0 - eps_b * (pressure_ratio ** (1.0 / 3.0)) * eps_b_13
    expected = (
        ((2.0 / (9.0 * np.pi)) ** (1.0 / 3.0))
        * eps_b_13
        * (pressure_ratio ** (P0_HAMEL / P))
        / denom
        * d_b
        / (3.0 * lam_b * u_b)
    )
    got = bubble_diameter_ode(
        0.0,
        np.array([d_b]),
        u0=u0,
        u_mf=u_mf,
        P=P,
        u_d=u_d,
        psi_b=psi_b,
        velocity_strategy="heinbockel_eq343",
        ode_strategy="heinbockel_eq341",
    )[0]
    assert np.isclose(got, expected)
    assert got > 0.0
