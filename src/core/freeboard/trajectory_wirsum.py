"""Analytical Wirsum trajectory advance.

Ref: Hamel (1999) Ch.3 freeboard particle trajectory; OPT-003 split from freeboard_segment.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.constants import Rg, g
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    N_GAS,
    gas_density_ideal,
    gas_diffusivity_correlation,
    gas_viscosity_power_law,
)
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R6, rate_R7, rate_R8, rate_R12
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from src.physics.freeboard import (
    calc_beta_a,
    calc_cd_hamel_eq384,
    calc_cd_haider,
    calc_terminal_velocity_haider,
    calc_u_gb,
)
from src.core.freeboard.trajectory_common import (
    _estimate_exact_hamel_critical_diameter,
    _find_first_event_time,
    _find_first_event_time_with_growth,
    _haider_high_re_drag_constant,
    _hamel_delta_neg_boundary_velocity,
    _hamel_delta_neg_max_height,
    _hamel_delta_neg_real_roots,
    _hamel_delta_neg_time_between_velocities,
    _hamel_delta_neg_zero_reachable,
    _hamel_exact_high_re_drag_constant,
    _quadratic_ode_state,
    _wirsum_abc,
    _wirsum_abc_exact_audit,
)

def _advance_particle_samples_analytical_wirsum(
    *,
    dh: float,
    z_lower_abs: float = 0.0,
    u_g_prev: float,
    u_g_next: float,
    rho_g: float,
    rho_p: float,
    d_p: np.ndarray,
    phi_s: float,
    mu_g: float,
    u_p_samples_prev: np.ndarray,
    m_dot_samples_prev: np.ndarray,
    char_frac_samples: np.ndarray,
    slower_ref_u_prev: np.ndarray | None = None,
    slower_ref_z_prev: np.ndarray | None = None,
    coeff_model: str = "stable_mixed_drag_split",
    n_subsegments: int = 6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Advance solids with boundary-event analytical trajectory updates.

    ``n_subsegments`` is intentionally ignored now. Hamel/Wirsum uses
    boundary-to-boundary analytical residence times rather than internal
    pseudo-spatial stepping within a freeboard cell.
    """
    u_p_samples = np.array(u_p_samples_prev, dtype=np.float64, copy=True)
    m_dot_samples = np.array(m_dot_samples_prev, dtype=np.float64, copy=True)
    hold_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    return_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    slower_ref_u_samples = np.array(
        u_p_samples_prev if slower_ref_u_prev is None else slower_ref_u_prev,
        dtype=np.float64,
        copy=True,
    )
    slower_ref_z_samples = np.array(
        np.zeros_like(u_p_samples_prev) if slower_ref_z_prev is None else slower_ref_z_prev,
        dtype=np.float64,
        copy=True,
    )
    diag = {
        "delta_pos": 0,
        "delta_zero": 0,
        "delta_neg": 0,
        "disc_pos": 0,
        "disc_zero": 0,
        "disc_neg": 0,
        "sign_switch": 0,
        "fallbacks": 0,
        "returns": 0,
        "fallback_no_event": 0,
        "fallback_boundary_velocity": 0,
        "fallback_delta_neg_time": 0,
        "re_p_min": float("inf"),
        "re_p_max": 0.0,
        "cd_min": float("inf"),
        "cd_max": 0.0,
        "re_p_active_min": float("inf"),
        "re_p_active_max": 0.0,
        "cd_active_min": float("inf"),
        "cd_active_max": 0.0,
        "z_drag": 0.0,
        "exact_delta_pos": 0,
        "exact_delta_zero": 0,
        "exact_delta_neg": 0,
        "exact_a_rel_max": 0.0,
        "exact_b_rel_max": 0.0,
        "exact_c_rel_max": 0.0,
        "exact_dp_grav_min": float("inf"),
        "exact_dp_grav_max": 0.0,
        "exact_dp_ratio_min": float("inf"),
        "exact_dp_ratio_max": 0.0,
    }
    u_g = float(0.5 * (u_g_prev + u_g_next))
    z_drag = (
        _hamel_exact_high_re_drag_constant(float(phi_s))
        if coeff_model == "exact_hamel"
        else _haider_high_re_drag_constant(float(phi_s))
    )
    diag["z_drag"] = float(z_drag)

    for j in range(len(u_p_samples)):
        if m_dot_samples_prev[j] <= 0.0:
            u_p_samples[j] = max(u_p_samples_prev[j], 1e-6)
            m_dot_samples[j] = 0.0
            continue

        z = float(z_lower_abs)
        t_seg = 0.0
        u_p = max(float(u_p_samples_prev[j]), 1e-6)
        survived = False
        branch_mode = "slower" if u_p < u_g else "faster"
        slower_ref_u = float(slower_ref_u_samples[j]) if branch_mode == "slower" else None
        slower_ref_z = float(slower_ref_z_samples[j]) if branch_mode == "slower" else None
        for _ in range(6):
            if z >= float(z_lower_abs + dh) - 1e-10:
                survived = True
                break

            rel = abs(float(u_g) - float(u_p))
            re_p = max(float(rho_g) * rel * max(float(d_p[j]), 1e-9) / max(float(mu_g), 1e-12), 1e-9)
            cd = calc_cd_haider(re_p, float(np.clip(phi_s, 1e-3, 1.0)))
            diag["re_p_min"] = min(float(diag["re_p_min"]), float(re_p))
            diag["re_p_max"] = max(float(diag["re_p_max"]), float(re_p))
            diag["cd_min"] = min(float(diag["cd_min"]), float(cd))
            diag["cd_max"] = max(float(diag["cd_max"]), float(cd))
            if rel >= 1e-3:
                diag["re_p_active_min"] = min(float(diag["re_p_active_min"]), float(re_p))
                diag["re_p_active_max"] = max(float(diag["re_p_active_max"]), float(re_p))
                diag["cd_active_min"] = min(float(diag["cd_active_min"]), float(cd))
                diag["cd_active_max"] = max(float(diag["cd_active_max"]), float(cd))

            a_ex, b_ex, c_ex, _, _, _ = _wirsum_abc_exact_audit(
                u_g=u_g,
                u_p=u_p,
                rho_g=rho_g,
                rho_p=rho_p,
                d_p=float(d_p[j]),
                mu_g=mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            if coeff_model == "exact_hamel":
                a, b, c = a_ex, b_ex, c_ex
            elif coeff_model == "stable_mixed_drag_split":
                a, b, c = _wirsum_abc(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=mu_g,
                    phi_s=phi_s,
                    branch=branch_mode,
                )
            else:
                raise ValueError(f"Unknown analytical_wirsum coeff_model: {coeff_model}")
            delta_exact = 4.0 * a_ex * c_ex - b_ex * b_ex
            if delta_exact > 1e-14:
                diag["exact_delta_pos"] += 1
            elif delta_exact < -1e-14:
                diag["exact_delta_neg"] += 1
            else:
                diag["exact_delta_zero"] += 1
            d_p_grav = _estimate_exact_hamel_critical_diameter(
                u_g=u_g,
                u_p=u_p,
                rho_g=rho_g,
                rho_p=rho_p,
                mu_g=mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            if d_p_grav is not None:
                diag["exact_dp_grav_min"] = min(float(diag["exact_dp_grav_min"]), float(d_p_grav))
                diag["exact_dp_grav_max"] = max(float(diag["exact_dp_grav_max"]), float(d_p_grav))
                ratio = float(d_p[j]) / max(float(d_p_grav), 1e-12)
                diag["exact_dp_ratio_min"] = min(float(diag["exact_dp_ratio_min"]), ratio)
                diag["exact_dp_ratio_max"] = max(float(diag["exact_dp_ratio_max"]), ratio)
            diag["exact_a_rel_max"] = max(
                float(diag["exact_a_rel_max"]),
                float(abs(a_ex - a) / max(abs(a), abs(a_ex), 1e-12)),
            )
            diag["exact_b_rel_max"] = max(
                float(diag["exact_b_rel_max"]),
                float(abs(b_ex - b) / max(abs(b), abs(b_ex), 1e-12)),
            )
            diag["exact_c_rel_max"] = max(
                float(diag["exact_c_rel_max"]),
                float(abs(c_ex - c) / max(abs(c), abs(c_ex), 1e-12)),
            )
            _, _, delta = _quadratic_ode_state(u_p, a, b, c, 0.0)
            if delta > 1e-14:
                diag["delta_pos"] += 1
                diag["disc_neg"] += 1
            elif delta < -1e-14:
                diag["delta_neg"] += 1
                diag["disc_pos"] += 1
            else:
                diag["delta_zero"] += 1
                diag["disc_zero"] += 1

            z_upper_abs = float(z_lower_abs + dh)

            t_scale = dh / max(max(abs(u_p), abs(u_g), 1e-3), 1e-6)
            t_max = max(20.0 * t_scale, 1.0)

            def vel_at(dt: float) -> float:
                u_dt, _, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(u_dt)

            def z_at(dt: float) -> float:
                _, dz_dt, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(z + dz_dt)

            t_switch = None
            if branch_mode == "faster":
                t_switch = _find_first_event_time(vel_at, target=u_g, t_max=t_max)
            t_turn = None
            if u_p > 0.0:
                t_turn = _find_first_event_time(vel_at, target=0.0, t_max=t_max)
            t_top = _find_first_event_time(z_at, target=z_upper_abs, t_max=t_max)

            candidates = [t for t in (t_switch, t_turn, t_top) if t is not None and t > 1e-10]
            if not candidates:
                if branch_mode == "slower" and u_p > 0.0:
                    ref_u = float(slower_ref_u if slower_ref_u is not None else u_p)
                    ref_z = float(slower_ref_z if slower_ref_z is not None else z)
                    t_to_current = _hamel_delta_neg_time_between_velocities(
                        u_start=ref_u,
                        u_end=float(u_p),
                        a=a,
                        b=b,
                        c=c,
                    )
                    zero_reachable = _hamel_delta_neg_zero_reachable(a, b, c)
                    roots = _hamel_delta_neg_real_roots(a, b, c)

                    if not zero_reachable:
                        u_lo = 1e-9
                        if roots is not None:
                            u_lo = max(float(roots[0]) + 1e-9, 1e-9)
                        # Hamel Fall 2 with two positive roots corresponds to
                        # a particle that never reaches u_p = 0. In that case
                        # Eq. 3.103 / h_max is not the right discriminant;
                        # the particle remains elutriated and we must recover
                        # the top-boundary crossing from the time domain.
                        if coeff_model == "exact_hamel":
                            t_cap_long = max(120.0, 8.0 * t_max)
                            if roots is not None:
                                u_lim = max(float(roots[0]), 1e-6)
                                t_cap_long = max(t_cap_long, 6.0 * dh / u_lim)
                            t_top_long = _find_first_event_time_with_growth(
                                z_at,
                                target=z_upper_abs,
                                t_init=t_max,
                                t_cap=t_cap_long,
                            )
                            if t_top_long is not None and t_top_long > 1e-10:
                                u_evt, dz_evt, _ = _quadratic_ode_state(u_p, a, b, c, t_top_long)
                                t_seg += float(t_top_long)
                                u_p = max(float(u_evt), 1e-9)
                                z = float(z_upper_abs)
                                survived = True
                                break
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=ref_u,
                            h0=ref_z,
                            a=a,
                            b=b,
                            c=c,
                            u_lo=u_lo,
                            u_hi=max(float(u_p), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_total = _hamel_delta_neg_time_between_velocities(
                                u_start=ref_u,
                                u_end=float(u_bound_up),
                                a=a,
                                b=b,
                                c=c,
                            )
                            if (
                                t_total is not None
                                and t_to_current is not None
                                and math.isfinite(t_total)
                                and math.isfinite(t_to_current)
                            ):
                                t_cross = max(float(t_total - t_to_current), 0.0)
                                t_seg += float(t_cross)
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                survived = True
                                break
                    h_max = _hamel_delta_neg_max_height(
                        u0=ref_u,
                        h0=ref_z,
                        a=a,
                        b=b,
                        c=c,
                    )
                    if h_max is not None and h_max < float(z_upper_abs) - 1e-10:
                        t_to_apex_total = _hamel_delta_neg_time_between_velocities(
                            u_start=ref_u,
                            u_end=0.0,
                            a=a,
                            b=b,
                            c=c,
                        )
                        u_bound_down = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_lower_abs),
                            u0=0.0,
                            h0=float(h_max),
                            a=a,
                            b=b,
                            c=c,
                        )
                        t_down = _hamel_delta_neg_time_between_velocities(
                            u_start=0.0,
                            u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                            a=a,
                            b=b,
                            c=c,
                        )
                        if (
                            u_bound_down is not None
                            and t_to_apex_total is not None
                            and t_down is not None
                            and t_to_current is not None
                            and math.isfinite(t_to_apex_total)
                            and math.isfinite(t_down)
                            and math.isfinite(t_to_current)
                        ):
                            t_back = max(float(t_to_apex_total - t_to_current), 0.0) + max(float(t_down), 0.0)
                            t_seg += float(t_back)
                            diag["returns"] += 1
                            return_samples[j] = float(m_dot_samples_prev[j])
                            u_p_samples[j] = 0.0
                            m_dot_samples[j] = 0.0
                            break
                    elif h_max is not None and h_max >= float(z_upper_abs) - 1e-10:
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=ref_u,
                            h0=ref_z,
                            a=a,
                            b=b,
                            c=c,
                            u_lo=1e-9,
                            u_hi=max(float(u_p), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_total = _hamel_delta_neg_time_between_velocities(
                                u_start=ref_u,
                                u_end=float(u_bound_up),
                                a=a,
                                b=b,
                                c=c,
                            )
                            if (
                                t_total is not None
                                and t_to_current is not None
                                and math.isfinite(t_total)
                                and math.isfinite(t_to_current)
                            ):
                                t_cross = max(float(t_total - t_to_current), 0.0)
                                t_seg += float(t_cross)
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                survived = True
                                break
                diag["fallbacks"] += 1
                diag["fallback_no_event"] += 1
                return_samples[j] = float(m_dot_samples_prev[j])
                u_p_samples[j] = 0.0
                m_dot_samples[j] = 0.0
                break

            t_evt = min(candidates)
            u_evt, dz_evt, _ = _quadratic_ode_state(u_p, a, b, c, t_evt)
            z = float(z + dz_evt)
            t_seg += float(t_evt)

            if t_top is not None and abs(t_evt - t_top) <= 1e-8 * max(1.0, t_top):
                u_p = max(float(u_evt), 1e-6)
                z = float(z_upper_abs)
                survived = True
                break

            if t_switch is not None and abs(t_evt - t_switch) <= 1e-8 * max(1.0, t_switch):
                diag["sign_switch"] += 1
                # Preserve the analytical event state continuously at the
                # branch switch instead of hard-resetting the particle speed.
                u_p = max(float(u_g), 1e-9)
                slower_ref_u = float(u_g)
                slower_ref_z = float(z)
                slower_ref_u_samples[j] = float(u_g)
                slower_ref_z_samples[j] = float(z)
                branch_mode = "slower"
                continue

            if t_turn is not None and abs(t_evt - t_turn) <= 1e-8 * max(1.0, t_turn):
                # Hamel Fall 2 uses the signed particle velocity directly.
                # After the apex, solve the lower-boundary return using
                # Eq. 3.103 (height as a function of signed velocity) and
                # the exact Delta<0 time primitive instead of restarting
                # a pseudo-time march from u=0 and relying on fallback.
                u_p = min(float(u_evt), 0.0)
                a2, b2, c2 = _wirsum_abc(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=mu_g,
                    phi_s=phi_s,
                    branch="slower",
                ) if coeff_model == "stable_mixed_drag_split" else _wirsum_abc_exact_audit(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=mu_g,
                    phi_s=phi_s,
                    branch="slower",
                )[:3]
                apex_z = float(z)
                u_bound = _hamel_delta_neg_boundary_velocity(
                    h_target=float(z_lower_abs),
                    u0=0.0,
                    h0=apex_z,
                    a=a2,
                    b=b2,
                    c=c2,
                )
                if u_bound is None:
                    diag["fallbacks"] += 1
                    diag["fallback_boundary_velocity"] += 1
                    return_samples[j] = float(m_dot_samples_prev[j])
                    u_p_samples[j] = 0.0
                    m_dot_samples[j] = 0.0
                    break
                t_back = _hamel_delta_neg_time_between_velocities(
                    u_start=0.0,
                    u_end=float(u_bound),
                    a=a2,
                    b=b2,
                    c=c2,
                )
                if (
                    t_back is None
                    or not math.isfinite(t_back)
                ):
                    diag["fallbacks"] += 1
                    diag["fallback_delta_neg_time"] += 1
                    return_samples[j] = float(m_dot_samples_prev[j])
                    u_p_samples[j] = 0.0
                    m_dot_samples[j] = 0.0
                    break
                t_seg += max(float(t_back), 0.0)
                diag["returns"] += 1
                return_samples[j] = float(m_dot_samples_prev[j])
                u_p_samples[j] = 0.0
                m_dot_samples[j] = 0.0
                break

        hold_samples[j] = float(m_dot_samples_prev[j]) * max(t_seg, 0.0)
        if survived:
            u_p_samples[j] = max(u_p, 1e-6)
            if branch_mode == "slower" and slower_ref_u is not None and slower_ref_z is not None:
                # Crossing into the next cell changes the local background-gas
                # state and therefore the frozen Fall-2 coefficients. Hamel's
                # boundary-to-boundary construction carries the *boundary*
                # state forward as the next cell's initial condition.
                slower_ref_u_samples[j] = max(float(u_p), 1e-6)
                slower_ref_z_samples[j] = float(z)
            else:
                slower_ref_u_samples[j] = max(float(u_p), 1e-6)
                slower_ref_z_samples[j] = 0.0
        else:
            return_samples[j] = float(m_dot_samples_prev[j])
            u_p_samples[j] = 0.0
            m_dot_samples[j] = 0.0
            slower_ref_u_samples[j] = max(float(u_p_samples_prev[j]), 1e-6)
            slower_ref_z_samples[j] = 0.0

    return (
        u_p_samples,
        m_dot_samples,
        hold_samples,
        return_samples,
        slower_ref_u_samples,
        slower_ref_z_samples,
        diag,
    )

