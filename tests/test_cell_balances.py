from __future__ import annotations

import numpy as np
import pytest

from src.core import species as species_mod
from src.core.cell_balances import (
    assemble_cell_residual_vector,
    calc_energy_balance_residual,
    calc_gas_balance_residual,
    calc_size_migration,
    calc_solid_balance_residual,
)
from src.core.cell_hydrodynamics import calc_cell_hydrodynamics, calc_phase_exchange
from src.core.species import configure_tar_components_by_fuel
from src.physics.freeboard import calc_beta_a, calc_u_gb


def test_gas_balance_uses_shared_exchange_with_opposite_sign() -> None:
    n_ex = np.array([1.0, -2.0, 0.5] + [0.0] * 8, dtype=np.float64)
    zero = np.zeros_like(n_ex)

    res = calc_gas_balance_residual(
        N_zu_d=zero,
        N_rez_d=zero,
        N_d_in=zero,
        R_gas_d=zero,
        N_d=zero,
        N_ex=n_ex,
        N_zu_b=zero,
        N_rez_b=zero,
        N_b_in=zero,
        R_gas_b=zero,
        N_b=zero,
    )

    np.testing.assert_allclose(res[:11], n_ex)
    np.testing.assert_allclose(res[11:], -n_ex)
    np.testing.assert_allclose(res[:11] + res[11:], 0.0)


def test_size_migration_transfers_mass_from_larger_to_smaller_class() -> None:
    m_solid = np.array(
        [
            [1.0, 0.5, 0.0, 0.1],
            [2.0, 1.0, 0.0, 0.2],
            [3.0, 1.5, 0.0, 0.3],
        ],
        dtype=np.float64,
    )
    r_solid = np.zeros_like(m_solid)
    r_solid[:, 0] = np.array([-0.03, -0.06, -0.09], dtype=np.float64)
    areas = np.array([0.4, 0.8, 1.2], dtype=np.float64)
    d_p = np.array([1.0e-3, 2.0e-3, 3.0e-3], dtype=np.float64)

    res = calc_size_migration(
        m_solid=m_solid,
        R_solid=r_solid,
        areas=areas,
        d_p_classes=d_p,
        rho_s=1400.0,
        eps_mf=0.45,
        V_d=0.08,
        V_cell=0.1,
        char_index=0,
    )

    assert res.shape == m_solid.shape
    assert np.sum(res[0]) > 0.0
    assert np.sum(res[-1]) < 0.0
    assert res[0, 0] > 0.0
    assert res[-1, 0] < 0.0


def test_solid_balance_adds_migration_and_sources() -> None:
    shape = (2, 4)
    m_solid = np.full(shape, 1.0, dtype=np.float64)
    res = calc_solid_balance_residual(
        m_solid_zu=np.full(shape, 0.3, dtype=np.float64),
        m_solid_rez=np.full(shape, 0.2, dtype=np.float64),
        m_solid_in=np.full(shape, 0.1, dtype=np.float64),
        R_solid=np.full(shape, -0.05, dtype=np.float64),
        size_migration=np.full(shape, 0.02, dtype=np.float64),
        m_solid=m_solid,
    )
    np.testing.assert_allclose(res, -0.43)


def test_solid_balance_holdup_transport_uses_frozen_neighbor_inflows_and_k_out() -> None:
    shape = (2, 4)
    m_solid = np.full(shape, 2.0, dtype=np.float64)
    res = calc_solid_balance_residual(
        m_solid_zu=np.full(shape, 0.3, dtype=np.float64),
        m_solid_rez=np.full(shape, 0.2, dtype=np.float64),
        m_solid_in=np.zeros(shape, dtype=np.float64),
        m_solid_auf_in=np.full(shape, 0.4, dtype=np.float64),
        m_solid_ab_in=np.full(shape, 0.1, dtype=np.float64),
        R_solid=np.full(shape, -0.05, dtype=np.float64),
        size_migration=np.full(shape, 0.02, dtype=np.float64),
        m_solid=m_solid,
        K_solid_auf=np.full(shape, 0.1, dtype=np.float64),
        K_solid_ab=np.full(shape, 0.2, dtype=np.float64),
        solid_state_model="holdup_transport",
    )
    expected = np.full(shape, 0.47, dtype=np.float64)
    expected[:, 0] = 0.37  # char
    expected[:, 3] = 0.37  # ash
    np.testing.assert_allclose(res, expected)


def test_solid_balance_holdup_transport_anchors_disconnected_holdup_buckets() -> None:
    shape = (2, 4)
    m_solid = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    res = calc_solid_balance_residual(
        m_solid_zu=np.zeros(shape, dtype=np.float64),
        m_solid_rez=np.zeros(shape, dtype=np.float64),
        m_solid_in=np.zeros(shape, dtype=np.float64),
        m_solid_auf_in=np.zeros(shape, dtype=np.float64),
        m_solid_ab_in=np.zeros(shape, dtype=np.float64),
        R_solid=np.zeros(shape, dtype=np.float64),
        size_migration=np.zeros(shape, dtype=np.float64),
        m_solid=m_solid,
        K_solid_auf=np.zeros(shape, dtype=np.float64),
        K_solid_ab=np.zeros(shape, dtype=np.float64),
        solid_state_model="holdup_transport",
    )
    np.testing.assert_allclose(res[:, 0], -m_solid[:, 0])
    np.testing.assert_allclose(res[:, 3], 0.0)


def test_solid_balance_holdup_transport_anchors_source_without_transport_channel() -> None:
    shape = (1, 4)
    m_solid = np.zeros(shape, dtype=np.float64)
    m_solid[0, 0] = 0.25
    m_auf_in = np.zeros(shape, dtype=np.float64)
    m_auf_in[0, 0] = 1.0
    res = calc_solid_balance_residual(
        m_solid_zu=np.zeros(shape, dtype=np.float64),
        m_solid_rez=np.zeros(shape, dtype=np.float64),
        m_solid_in=np.zeros(shape, dtype=np.float64),
        m_solid_auf_in=m_auf_in,
        m_solid_ab_in=np.zeros(shape, dtype=np.float64),
        R_solid=np.zeros(shape, dtype=np.float64),
        size_migration=np.zeros(shape, dtype=np.float64),
        m_solid=m_solid,
        K_solid_auf=np.zeros(shape, dtype=np.float64),
        K_solid_ab=np.zeros(shape, dtype=np.float64),
        solid_state_model="holdup_transport",
    )
    np.testing.assert_allclose(res[0, 0], 0.75)
    np.testing.assert_allclose(res[0, 3], 0.0)


def test_energy_balance_zero_for_matched_inlet_and_outlet_streams() -> None:
    saved_mw = {
        "TAR1": species_mod.MOLECULAR_WEIGHT["TAR1"],
        "TAR2": species_mod.MOLECULAR_WEIGHT["TAR2"],
    }
    saved_nasa = {
        key: species_mod._NASA_DATA.get(key)  # type: ignore[attr-defined]
        for key in ("TAR1", "TAR2")
    }
    configure_tar_components_by_fuel("coal")
    gas = np.zeros(11, dtype=np.float64)
    gas[0] = 1.5
    gas[4] = 2.0
    solid = np.array([[0.4, 0.2, 0.1, 0.3]], dtype=np.float64)

    try:
        residual = calc_energy_balance_residual(
            N_b_in=gas * 0.3,
            N_d_in=gas * 0.7,
            T_in_gas=1100.0,
            N_zu_b=np.zeros_like(gas),
            N_zu_d=np.zeros_like(gas),
            T_zu_gas=1100.0,
            N_rez_b=np.zeros_like(gas),
            N_rez_d=np.zeros_like(gas),
            T_rez_gas=1100.0,
            m_solid_rez=np.zeros_like(solid),
            T_rez_solid=1100.0,
            m_solid_in=np.zeros_like(solid),
            T_in_solid=1100.0,
            m_solid_zu=solid,
            T_zu_solid=1100.0,
            N_b=gas * 0.3,
            N_d=gas * 0.7,
            m_solid=solid,
            T=1100.0,
            heat_loss_frac=0.0,
            ash_dry_wt=11.41,
            VM_daf=53.42,
            h_f_dry=-5.0e5,
            h_cache={},
        )
        assert abs(residual) < 1e-9
    finally:
        species_mod.MOLECULAR_WEIGHT["TAR1"] = saved_mw["TAR1"]
        species_mod.MOLECULAR_WEIGHT["TAR2"] = saved_mw["TAR2"]
        for key in ("TAR1", "TAR2"):
            if saved_nasa[key] is None:
                species_mod._NASA_DATA.pop(key, None)  # type: ignore[attr-defined]
            else:
                species_mod._NASA_DATA[key] = saved_nasa[key]  # type: ignore[attr-defined]


def test_assemble_cell_residual_vector_packs_sections() -> None:
    out = np.zeros(8, dtype=np.float64)
    packed = assemble_cell_residual_vector(
        res_gas=np.array([1.0, 2.0, 3.0], dtype=np.float64),
        res_solid=np.array([[4.0, 5.0], [6.0, 7.0]], dtype=np.float64),
        res_energy=8.0,
        out=out,
    )
    np.testing.assert_allclose(packed, np.arange(1.0, 9.0))


def test_hydrodynamics_bundle_stays_in_physical_bounds() -> None:
    N_d = np.zeros(11, dtype=np.float64)
    N_b = np.zeros(11, dtype=np.float64)
    N_d[5] = 1.5
    N_d[6] = 10.0
    N_b[5] = 0.3
    N_b[6] = 3.0
    bundle = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
    )
    assert 0.0 < bundle.u_mf < bundle.u_b
    assert 0.01 <= bundle.eps_b <= 0.7
    assert np.isclose(bundle.eps_b + bundle.eps_d, 1.0)
    assert bundle.u_d > 0.0
    assert bundle.n_rz > 0.0
    assert bundle.eps_d_voidage >= 0.45
    assert bundle.V_b > 0.0 and bundle.V_d > 0.0
    assert bundle.K_bd > 0.0


def test_ud_closure_alternatives_prevent_visible_bubble_collapse() -> None:
    N_d = np.zeros(11, dtype=np.float64)
    N_b = np.zeros(11, dtype=np.float64)
    N_d[5] = 1.5
    N_d[6] = 10.0
    N_b[5] = 0.3
    N_b[6] = 3.0

    bundle_current = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="current",
    )
    bundle_umf = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="umf_over_epsmf",
    )
    bundle_backsolve = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="backsolve_visible_epsb",
    )
    bundle_hilligardt = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="hilligardt_eq311",
    )
    bundle_wein = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="wein_1992_eq312",
    )

    assert np.isclose(bundle_current.eps_b, bundle_umf.eps_b)
    assert np.isclose(bundle_current.eps_b, bundle_hilligardt.eps_b)
    assert np.isclose(bundle_current.eps_b, bundle_wein.eps_b)
    assert np.isclose(bundle_backsolve.eps_b, bundle_current.eps_b)
    assert bundle_umf.eps_d_voidage < bundle_current.eps_d_voidage
    assert bundle_hilligardt.u_d > bundle_wein.u_d
    assert np.isclose(bundle_wein.u_d, 1.45 * bundle_wein.u_mf)
    assert np.isclose(
        bundle_hilligardt.u_d,
        bundle_hilligardt.u_mf + (bundle_hilligardt.u0 - bundle_hilligardt.u_mf) / 3.0,
    )
    assert bundle_backsolve.eps_d_voidage <= bundle_umf.eps_d_voidage


def test_wein_ud_closure_uses_micro_bubbling_guard_when_ud_would_exceed_u0() -> None:
    N_d = np.zeros(11, dtype=np.float64)
    N_b = np.zeros(11, dtype=np.float64)
    N_d[5] = 1.5
    N_d[6] = 10.0
    N_b[5] = 0.3
    N_b[6] = 3.0

    bundle_wein = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=0.25,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=0.18,
        u_d_closure="wein_1992_eq312",
    )

    u_d_potential = 1.45 * bundle_wein.u_mf
    u_d_guard = bundle_wein.u_mf + (bundle_wein.u0 - bundle_wein.u_mf) * 0.999

    assert u_d_potential > bundle_wein.u0
    assert bundle_wein.u_d < bundle_wein.u0
    assert np.isclose(bundle_wein.u_d, u_d_guard)
    assert bundle_wein.u_d < u_d_potential


def test_freeboard_cell_hydrodynamics_uses_ghost_bubble_decay_and_zero_exchange() -> None:
    N_d = np.zeros(11, dtype=np.float64)
    N_d[5] = 1.5
    N_d[6] = 10.0

    bundle = calc_cell_hydrodynamics(
        T=1125.0,
        P=2.5e6,
        N_b=np.zeros(11, dtype=np.float64),
        N_d=N_d,
        D_bed=0.6,
        dh=0.4,
        h_center=0.2,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        cell_type="freeboard",
        freeboard_u_bed_top=1.8,
        freeboard_d_b_bed_top=0.08,
        freeboard_height_from_bed=0.4,
        freeboard_eps_d_voidage=0.995,
    )

    u_gb0 = calc_u_gb(1.8)
    beta_a = calc_beta_a(0.08)
    expected_u_g = bundle.u0 + (u_gb0 - bundle.u0) * np.exp(-beta_a * 0.4)

    assert bundle.u_mf == 0.0
    assert bundle.eps_b == 0.0
    assert bundle.K_bd == 0.0
    assert bundle.V_b == 0.0
    assert np.isclose(bundle.V_d, np.pi / 4.0 * 0.6**2 * 0.4)
    assert np.isclose(bundle.eps_d_voidage, 0.995)
    assert np.isclose(bundle.u_d, expected_u_g)
    assert np.isclose(bundle.u_b, expected_u_g)


def test_freeboard_velocity_samples_keep_particle_launch_speed_distinct_from_ghost_bubble_speed() -> None:
    from src.core.freeboard_segment import _build_velocity_samples

    u_b_ws = 1.8
    u_gb0 = calc_u_gb(u_b_ws)
    samples, _ = _build_velocity_samples(1.53 * u_b_ws, sigma=0.6, n_bins=5)

    assert u_gb0 == pytest.approx(u_b_ws)
    assert float(np.mean(samples)) > u_gb0
    assert float(np.mean(samples)) == pytest.approx(1.53 * u_b_ws, rel=0.0, abs=1e-12)


def test_cyclone_and_return_leg_hydrodynamics_degenerate_from_bed_chain() -> None:
    N_d = np.zeros(11, dtype=np.float64)
    N_d[5] = 1.5
    N_d[6] = 10.0

    cyclone = calc_cell_hydrodynamics(
        T=1125.0,
        P=2.5e6,
        N_b=np.zeros(11, dtype=np.float64),
        N_d=N_d,
        D_bed=0.6,
        dh=0.4,
        h_center=0.2,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        cell_type="cyclone",
    )
    return_leg = calc_cell_hydrodynamics(
        T=1125.0,
        P=2.5e6,
        N_b=np.zeros(11, dtype=np.float64),
        N_d=N_d,
        D_bed=0.6,
        dh=0.4,
        h_center=0.2,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        cell_type="return_leg",
    )

    assert cyclone.eps_b == 0.0
    assert cyclone.K_bd == 0.0
    assert cyclone.u_d == cyclone.u0
    assert cyclone.eps_d_voidage == 1.0

    assert return_leg.eps_b == 0.0
    assert return_leg.K_bd == 0.0
    assert return_leg.u0 == 0.0
    assert return_leg.u_d == 0.0
    assert np.isclose(return_leg.eps_d_voidage, 0.45)


def test_hilligardt_ode_bubble_branch_is_available_for_hydrodynamics_audit() -> None:
    N_d = np.zeros(11, dtype=np.float64)
    N_b = np.zeros(11, dtype=np.float64)
    N_d[5] = 1.5
    N_d[6] = 10.0
    N_b[5] = 0.3
    N_b[6] = 3.0

    bundle_mw = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=2.75,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="backsolve_visible_epsb",
        bubble_diameter_model="mori_wen",
    )
    bundle_ode = calc_cell_hydrodynamics(
        T=1150.0,
        P=2.5e6,
        N_b=N_b,
        N_d=N_d,
        D_bed=0.6,
        dh=0.5,
        h_center=2.75,
        N_or=100,
        rho_s=1400.0,
        d_p=1.0e-3,
        eps_mf=0.45,
        phi_s=0.86,
        u0_target=1.2,
        u_d_closure="backsolve_visible_epsb",
        bubble_diameter_model="hilligardt_ode",
        psi_b_strategy="technical_distributor",
        lambda_strategy="hamel_280",
        xi_strategy="hamel_regime",
    )

    assert bundle_ode.d_b > 0.0
    assert bundle_ode.u_b > 0.0
    assert bundle_ode.u_d > 0.0
    assert 0.01 <= bundle_ode.eps_b <= 0.7
    assert bundle_ode.K_bd > 0.0
    assert not np.isclose(bundle_ode.d_b, bundle_mw.d_b)


def test_phase_exchange_is_zero_when_phase_concentrations_match() -> None:
    C = np.array([1.0, 2.0, 3.0] + [0.0] * 8, dtype=np.float64)
    n_ex = calc_phase_exchange(K_bd=3.0, V_b=0.1, C_b=C, C_d=C)
    np.testing.assert_allclose(n_ex, 0.0)
