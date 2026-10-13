from __future__ import annotations

import numpy as np

from src.core.freeboard_segment import (
    _advance_particle_samples_analytical_wirsum,
    _freeboard_reaction_step,
    _quadratic_ode_step,
    _wirsum_abc_exact_audit,
    simulate_freeboard,
)
from src.core.reactor import GAS_SPECIES, GAS_SPECIES_INDEX
from src.core.species import configure_tar_components_by_fuel


def test_quadratic_ode_step_handles_negative_thesis_delta() -> None:
    u1, delta = _quadratic_ode_step(u0=1.2, a=0.5, b=0.0, c=-0.5, dt=0.02)
    assert delta < 0.0
    assert np.isfinite(u1)


def test_quadratic_ode_step_handles_zero_thesis_delta() -> None:
    u1, delta = _quadratic_ode_step(u0=0.8, a=1.0, b=2.0, c=1.0, dt=0.02)
    assert abs(delta) <= 1e-12
    assert np.isfinite(u1)


def test_quadratic_ode_step_handles_positive_thesis_delta() -> None:
    u1, delta = _quadratic_ode_step(u0=0.6, a=1.0, b=0.0, c=1.0, dt=0.01)
    assert delta > 0.0
    assert np.isfinite(u1)


def test_hamel_freeboard_initial_velocity_and_eq397() -> None:
    from src.core.freeboard.trajectory_common import (
        calc_critical_diameter_eq397,
        hamel_delta_pos_height_from_velocity,
        _wirsum_abc_exact_audit,
    )
    from src.physics.freeboard import calc_u_p0_demich_bohnet, calc_u_p0_martens, calc_u_p0_son

    u_b = 0.8
    d_b = 0.25
    d_p = 1.0e-3
    rho_g = 0.35
    rho_p = 2600.0
    u70 = calc_u_p0_demich_bohnet(u_b, d_b, d_p, rho_g, rho_p)
    u71 = calc_u_p0_son(u_b, d_b)
    u_m = calc_u_p0_martens(u_b)
    assert u70 > 0.0 and u71 > 0.0 and u_m == 1.53 * u_b

    d_star = calc_critical_diameter_eq397(mu_g=4e-5, rho_g=rho_g, rho_p=rho_p, c_w=0.44)
    assert d_star is not None and d_star > 0.0

    # Eq. (3.96)：使用 Δ>0 且 q,q0>0 的构造系数（与上方 u_p0 公式独立）
    h = hamel_delta_pos_height_from_velocity(
        u=2.0,
        u0=1.0,
        h0=0.0,
        a=1.0,
        b=0.0,
        c=1.0,
    )
    assert h is not None and np.isfinite(h)


def test_exact_hamel_coefficients_do_not_tie_solution_branch_to_velocity_branch() -> None:
    a, b, c, re_p, c_w, z_drag = _wirsum_abc_exact_audit(
        u_g=1.595,
        u_p=3.1302,
        rho_g=7.511309531141224,
        rho_p=1700.0,
        d_p=7e-4,
        mu_g=1.902999373563184e-05,
        phi_s=0.86,
        branch="faster",
    )
    delta = 4.0 * a * c - b * b
    assert re_p > 0.0
    assert c_w > 0.0
    assert z_drag > 0.0
    assert np.isfinite(delta)


def test_simulate_freeboard_reports_analytical_solver_metadata() -> None:
    configure_tar_components_by_fuel("coal")
    idx = GAS_SPECIES_INDEX
    n = len(GAS_SPECIES)
    n_in = np.zeros(n, dtype=np.float64)
    n_in[idx["CO"]] = 1.0
    n_in[idx["H2"]] = 1.0
    n_in[idx["H2O"]] = 1.0
    n_in[idx["N2"]] = 3.0

    out = simulate_freeboard(
        N_in=n_in,
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

    assert out["trajectory_model"] == "analytical_wirsum"
    assert out["trajectory_solver"] == "wirsum_analytical_exact_hamel_coeffs"
    assert out["trajectory_coeff_model"] == "exact_hamel"
    assert isinstance(out["trajectory_diag"], dict)
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
    } <= set(out["trajectory_diag"].keys())
    assert isinstance(out["trajectory_coeff_diag"], dict)
    assert {
        "re_p_min",
        "re_p_max",
        "cd_min",
        "cd_max",
        "re_p_active_min",
        "re_p_active_max",
        "cd_active_min",
        "cd_active_max",
        "z_drag",
        "exact_delta_pos",
        "exact_delta_zero",
        "exact_delta_neg",
        "exact_a_rel_max",
        "exact_b_rel_max",
        "exact_c_rel_max",
        "exact_dp_grav_min",
        "exact_dp_grav_max",
        "exact_dp_ratio_min",
        "exact_dp_ratio_max",
    } <= set(out["trajectory_coeff_diag"].keys())
    assert out["trajectory_coeff_diag"]["re_p_max"] >= out["trajectory_coeff_diag"]["re_p_min"] >= 0.0
    assert out["trajectory_coeff_diag"]["cd_max"] >= out["trajectory_coeff_diag"]["cd_min"] >= 0.0
    assert out["trajectory_coeff_diag"]["re_p_active_max"] >= out["trajectory_coeff_diag"]["re_p_active_min"] >= 0.0
    assert out["trajectory_coeff_diag"]["cd_active_max"] >= out["trajectory_coeff_diag"]["cd_active_min"] >= 0.0
    assert out["trajectory_coeff_diag"]["z_drag"] > 0.0
    assert out["trajectory_coeff_diag"]["exact_delta_pos"] + out["trajectory_coeff_diag"]["exact_delta_zero"] + out["trajectory_coeff_diag"]["exact_delta_neg"] > 0
    assert out["trajectory_coeff_diag"]["exact_a_rel_max"] >= 0.0
    assert out["trajectory_coeff_diag"]["exact_b_rel_max"] >= 0.0
    assert out["trajectory_coeff_diag"]["exact_c_rel_max"] >= 0.0
    assert out["trajectory_coeff_diag"]["exact_dp_grav_max"] >= out["trajectory_coeff_diag"]["exact_dp_grav_min"] >= 0.0
    assert out["trajectory_coeff_diag"]["exact_dp_ratio_max"] >= out["trajectory_coeff_diag"]["exact_dp_ratio_min"] >= 0.0
    assert len(out["profiles"]["solid_K_auf_classes_1_s"]) == 2
    assert len(out["profiles"]["solid_K_ab_classes_1_s"]) == 2
    assert len(out["profiles"]["solid_K_auf_classes_1_s"][0]) == 1
    assert len(out["profiles"]["solid_holdup_char_classes_kg"][0]) == 1
    assert out["profiles"]["solid_K_auf_classes_1_s"][0][0] >= 0.0
    assert out["profiles"]["solid_K_ab_classes_1_s"][0][0] >= 0.0


def test_simulate_freeboard_accepts_exact_hamel_coefficient_mode() -> None:
    configure_tar_components_by_fuel("coal")
    idx = GAS_SPECIES_INDEX
    n = len(GAS_SPECIES)
    n_in = np.zeros(n, dtype=np.float64)
    n_in[idx["CO"]] = 1.0
    n_in[idx["H2"]] = 1.0
    n_in[idx["H2O"]] = 1.0
    n_in[idx["N2"]] = 3.0

    out = simulate_freeboard(
        N_in=n_in,
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
        trajectory_coeff_model="exact_hamel",
        enabled_reactions=("R5", "R6", "R7", "R8", "R12"),
    )

    assert out["trajectory_model"] == "analytical_wirsum"
    assert out["trajectory_coeff_model"] == "exact_hamel"
    assert out["trajectory_solver"] == "wirsum_analytical_exact_hamel_coeffs"
    assert isinstance(out["trajectory_coeff_diag"], dict)


def test_delta_neg_slow_upward_particle_can_cross_cell_boundary() -> None:
    u_p, m_dot, hold, ret, slower_u, slower_z, diag = _advance_particle_samples_analytical_wirsum(
        dh=8.5 / 8.0,
        u_g_prev=1.2524,
        u_g_next=1.2016,
        rho_g=7.6,
        rho_p=1200.0,
        d_p=np.array([1.0e-3], dtype=np.float64),
        phi_s=0.75,
        mu_g=1.9e-5,
        u_p_samples_prev=np.array([0.2359813463609112], dtype=np.float64),
        m_dot_samples_prev=np.array([1.0], dtype=np.float64),
        char_frac_samples=np.array([0.75], dtype=np.float64),
        coeff_model="exact_hamel",
    )

    assert diag["fallbacks"] == 0
    assert diag["fallback_no_event"] == 0
    assert float(np.sum(m_dot)) > 0.0
    assert float(np.sum(ret)) == 0.0
    assert float(np.sum(slower_u)) > 0.0
    assert np.all(np.isfinite(slower_z))


def test_freeboard_reaction_step_activates_char_heterogeneous_reactions() -> None:
    idx = GAS_SPECIES_INDEX
    n = np.zeros(len(GAS_SPECIES), dtype=np.float64)
    n[idx["O2"]] = 0.2
    n[idx["H2O"]] = 0.8
    n[idx["CO2"]] = 0.5
    n[idx["CO"]] = 0.2
    n[idx["H2"]] = 0.3
    n[idx["CH4"]] = 0.1
    n[idx["N2"]] = 1.0

    out, diag = _freeboard_reaction_step(
        N=n,
        T=1150.0,
        P=2.5e6,
        V_seg=0.2,
        dt=0.5,
        fuel_type="coal",
        rho_cat=20.0,
        char_area_total=3.0,
        solid_d_p=7e-4,
        D_g=2e-4,
        char_conversion=0.0,
        enabled_reactions=("R1", "R2", "R3", "R4"),
    )

    assert diag["R1"] > 0.0
    assert diag["R2"] > 0.0
    assert diag["R3"] >= 0.0
    assert diag["R4"] >= 0.0
    assert out[idx["O2"]] < n[idx["O2"]]
    assert out[idx["H2O"]] < n[idx["H2O"]]
    assert out[idx["CO"]] > n[idx["CO"]]
