"""Global projected Hamel trajectory + segment driver.

Ref: Hamel (1999) Ch.3 freeboard particle trajectory; OPT-003 split from freeboard_segment.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.energy_balance_defaults import (
    DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    apply_freeboard_substep_gas_enthalpy_heat_loss,
)
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
from src.core.freeboard.reaction import (
    _build_freeboard_background_profile,
    _estimate_rho_cat_profile_from_projected_transport,
    _freeboard_reaction_step,
)
from src.core.freeboard.trajectory_common import (
    _estimate_exact_hamel_critical_diameter,
    _find_first_event_time,
    _find_first_event_time_with_growth,
    _hamel_delta_neg_boundary_velocity,
    _hamel_delta_neg_max_height,
    _hamel_delta_neg_real_roots,
    _hamel_delta_neg_time_between_velocities,
    _hamel_delta_neg_zero_reachable,
    _quadratic_ode_state,
    _wirsum_abc,
    _wirsum_abc_exact_audit,
)
from src.core.freeboard.transport import (
    _aggregate_class_transport,
    _apply_char_sink_to_class_transport,
    bed_face_class_flows_kg_s,
    _build_size_class_bundles,
    _build_velocity_samples,
    _char_area_and_representative_diameter,
    _empty_class_transport,
    _freeboard_inventory_char_sink_diag,
    _representative_char_conversion,
    _resolve_transport_class_mapping,
    _solve_T_from_enthalpy,
    _wet_gas_profile,
)
from src.core.freeboard.types import FreeboardBackground, FreeboardState

def _advance_particle_samples_global_projected_hamel(
    *,
    backgrounds: list[FreeboardBackground],
    rho_p: float,
    d_p: np.ndarray,
    phi_s: float,
    u_p_samples_prev: np.ndarray,
    m_dot_samples_prev: np.ndarray,
    char_frac_samples: np.ndarray,
    char_conversion_samples: np.ndarray,
    class_idx_samples: np.ndarray,
    coeff_model: str,
) -> dict[str, object]:
    """Track each particle sample continuously across the whole freeboard and project back to cells."""
    n_cells = len(backgrounds)
    n_classes = int(np.max(class_idx_samples)) + 1 if len(class_idx_samples) else 0
    hold_samples_by_seg = [np.zeros_like(m_dot_samples_prev, dtype=np.float64) for _ in range(n_cells)]
    upflow_samples_by_seg = [np.zeros_like(m_dot_samples_prev, dtype=np.float64) for _ in range(n_cells)]
    downflow_samples_by_seg = [np.zeros_like(m_dot_samples_prev, dtype=np.float64) for _ in range(n_cells)]
    exit_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    exit_u_samples = np.zeros_like(u_p_samples_prev, dtype=np.float64)
    total_return_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
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
        "projection_crossings": 0,
    }
    coeff_diag = {
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

    for j in range(len(u_p_samples_prev)):
        if m_dot_samples_prev[j] <= 0.0:
            continue
        cell_idx = 0
        z = 0.0
        u_p = max(float(u_p_samples_prev[j]), 1e-6)
        ref_u = float(u_p)
        ref_z = 0.0
        steps = 0
        while 0 <= cell_idx < n_cells and steps < max(12 * n_cells, 32):
            steps += 1
            bg = backgrounds[cell_idx]
            z = min(max(float(z), float(bg.z_lower_m)), float(bg.z_upper_m))
            u_g = float(bg.u_gb)
            branch_mode = "slower" if u_p < u_g else "faster"
            rel = abs(float(u_g) - float(u_p))
            re_p = max(float(bg.rho_g) * rel * max(float(d_p[j]), 1e-9) / max(float(bg.mu_g), 1e-12), 1e-9)
            cd = calc_cd_haider(re_p, float(np.clip(phi_s, 1e-3, 1.0)))
            coeff_diag["re_p_min"] = min(float(coeff_diag["re_p_min"]), float(re_p))
            coeff_diag["re_p_max"] = max(float(coeff_diag["re_p_max"]), float(re_p))
            coeff_diag["cd_min"] = min(float(coeff_diag["cd_min"]), float(cd))
            coeff_diag["cd_max"] = max(float(coeff_diag["cd_max"]), float(cd))
            if rel >= 1e-3:
                coeff_diag["re_p_active_min"] = min(float(coeff_diag["re_p_active_min"]), float(re_p))
                coeff_diag["re_p_active_max"] = max(float(coeff_diag["re_p_active_max"]), float(re_p))
                coeff_diag["cd_active_min"] = min(float(coeff_diag["cd_active_min"]), float(cd))
                coeff_diag["cd_active_max"] = max(float(coeff_diag["cd_active_max"]), float(cd))

            a_ex, b_ex, c_ex, _, _, z_drag = _wirsum_abc_exact_audit(
                u_g=u_g,
                u_p=u_p,
                rho_g=bg.rho_g,
                rho_p=rho_p,
                d_p=float(d_p[j]),
                mu_g=bg.mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            coeff_diag["z_drag"] = float(max(float(coeff_diag["z_drag"]), float(z_drag)))
            if coeff_model == "exact_hamel":
                a, b, c = a_ex, b_ex, c_ex
            else:
                a, b, c = _wirsum_abc(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=bg.rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=bg.mu_g,
                    phi_s=phi_s,
                    branch=branch_mode,
                )
            delta_exact = 4.0 * a_ex * c_ex - b_ex * b_ex
            if delta_exact > 1e-14:
                coeff_diag["exact_delta_pos"] += 1
            elif delta_exact < -1e-14:
                coeff_diag["exact_delta_neg"] += 1
            else:
                coeff_diag["exact_delta_zero"] += 1
            d_p_grav = _estimate_exact_hamel_critical_diameter(
                u_g=u_g,
                u_p=u_p,
                rho_g=bg.rho_g,
                rho_p=rho_p,
                mu_g=bg.mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            if d_p_grav is not None:
                coeff_diag["exact_dp_grav_min"] = min(float(coeff_diag["exact_dp_grav_min"]), float(d_p_grav))
                coeff_diag["exact_dp_grav_max"] = max(float(coeff_diag["exact_dp_grav_max"]), float(d_p_grav))
                ratio = float(d_p[j]) / max(float(d_p_grav), 1e-12)
                coeff_diag["exact_dp_ratio_min"] = min(float(coeff_diag["exact_dp_ratio_min"]), ratio)
                coeff_diag["exact_dp_ratio_max"] = max(float(coeff_diag["exact_dp_ratio_max"]), ratio)

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

            z_lower_abs = float(bg.z_lower_m)
            z_upper_abs = float(bg.z_upper_m)
            t_scale = (z_upper_abs - z) / max(max(abs(u_p), abs(u_g), 1e-3), 1e-6)
            t_max = max(20.0 * max(float(t_scale), 1e-6), 1.0)

            def vel_at(dt: float) -> float:
                u_dt, _, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(u_dt)

            def z_at(dt: float) -> float:
                _, dz_dt, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(z + dz_dt)

            zero_reachable = _hamel_delta_neg_zero_reachable(a, b, c) if branch_mode == "slower" and u_p > 0.0 else True
            t_switch = _find_first_event_time(vel_at, target=u_g, t_max=t_max) if branch_mode == "faster" else None
            t_turn = _find_first_event_time(vel_at, target=0.0, t_max=t_max) if u_p > 0.0 and zero_reachable else None
            t_top = _find_first_event_time(z_at, target=z_upper_abs, t_max=t_max)
            candidates = [t for t in (t_switch, t_turn, t_top) if t is not None and t > 1e-10]

            if not candidates and branch_mode == "slower" and u_p > 0.0:
                roots = _hamel_delta_neg_real_roots(a, b, c)
                if not zero_reachable:
                    u_lim = max(float(roots[0]), 1e-6) if roots is not None else 1e-6
                    t_cap = max(120.0, 6.0 * max(z_upper_abs - z, 1e-9) / u_lim)
                    t_top = _find_first_event_time_with_growth(
                        z_at,
                        target=z_upper_abs,
                        t_init=t_max,
                        t_cap=t_cap,
                    )
                    if t_top is not None and t_top > 1e-10:
                        candidates = [float(t_top)]

            if not candidates:
                if branch_mode == "slower" and u_p > 0.0:
                    h_max = _hamel_delta_neg_max_height(
                        u0=float(ref_u),
                        h0=float(ref_z),
                        a=a,
                        b=b,
                        c=c,
                    )
                    if h_max is not None and h_max < z_upper_abs - 1e-10:
                        t_to_current = _hamel_delta_neg_time_between_velocities(
                            u_start=float(ref_u),
                            u_end=float(u_p),
                            a=a,
                            b=b,
                            c=c,
                        )
                        t_to_apex_total = _hamel_delta_neg_time_between_velocities(
                            u_start=float(ref_u),
                            u_end=0.0,
                            a=a,
                            b=b,
                            c=c,
                        )
                        a2, b2, c2, _, _, _ = _wirsum_abc_exact_audit(
                            u_g=u_g,
                            u_p=0.0,
                            rho_g=bg.rho_g,
                            rho_p=rho_p,
                            d_p=float(d_p[j]),
                            mu_g=bg.mu_g,
                            phi_s=phi_s,
                            branch="slower",
                        )
                        u_bound_down = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_lower_abs),
                            u0=0.0,
                            h0=float(h_max),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        t_down = _hamel_delta_neg_time_between_velocities(
                            u_start=0.0,
                            u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        if (
                            u_bound_down is not None
                            and t_to_current is not None
                            and t_to_apex_total is not None
                            and t_down is not None
                            and math.isfinite(t_to_current)
                            and math.isfinite(t_to_apex_total)
                            and math.isfinite(t_down)
                        ):
                            t_cell = max(float(t_to_apex_total - t_to_current), 0.0) + max(float(t_down), 0.0)
                            hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(t_cell, 0.0)
                            downflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                            total_return_samples[j] += float(m_dot_samples_prev[j])
                            diag["returns"] += 1
                            break
                diag["fallbacks"] += 1
                diag["fallback_no_event"] += 1
                break

            t_evt = min(candidates)
            u_evt, dz_evt, _ = _quadratic_ode_state(u_p, a, b, c, t_evt)
            hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_evt), 0.0)
            z = float(z + dz_evt)

            if t_top is not None and abs(t_evt - t_top) <= 1e-8 * max(1.0, t_top):
                upflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                u_p = max(float(u_evt), 1e-9)
                z = float(z_upper_abs)
                if branch_mode == "slower":
                    ref_u = max(float(u_p), 1e-9)
                    ref_z = float(z)
                diag["projection_crossings"] += 1
                cell_idx += 1
                if cell_idx >= n_cells:
                    exit_samples[j] = float(m_dot_samples_prev[j])
                    exit_u_samples[j] = float(u_p)
                    break
                continue

            if t_switch is not None and abs(t_evt - t_switch) <= 1e-8 * max(1.0, t_switch):
                diag["sign_switch"] += 1
                ref_u = float(u_g)
                ref_z = float(z)
                # In the thesis-aligned exact path, ``u_p = u_g`` is the
                # analytical branch handoff point. Resolve the remaining part
                # of the current cell with Fall 2 boundary formulas instead of
                # re-entering the time-domain marcher from a numerically
                # nudged velocity. This preserves per-sample boundary states
                # and avoids collapsing multiple trajectories onto a single
                # lower-root attractor inside the same cell.
                if coeff_model == "exact_hamel":
                    a2, b2, c2, _, _, _ = _wirsum_abc_exact_audit(
                        u_g=u_g,
                        u_p=float(u_g),
                        rho_g=bg.rho_g,
                        rho_p=rho_p,
                        d_p=float(d_p[j]),
                        mu_g=bg.mu_g,
                        phi_s=phi_s,
                        branch="slower",
                    )
                    zero_reachable_2 = _hamel_delta_neg_zero_reachable(a2, b2, c2)
                    roots_2 = _hamel_delta_neg_real_roots(a2, b2, c2)
                    if not zero_reachable_2:
                        u_lo = 1e-9
                        if roots_2 is not None:
                            u_lo = max(float(roots_2[0]) + 1e-9, 1e-9)
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=float(ref_u),
                            h0=float(ref_z),
                            a=a2,
                            b=b2,
                            c=c2,
                            u_lo=u_lo,
                            u_hi=max(float(ref_u), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_cross = _hamel_delta_neg_time_between_velocities(
                                u_start=float(ref_u),
                                u_end=float(u_bound_up),
                                a=a2,
                                b=b2,
                                c=c2,
                            )
                            if t_cross is not None and math.isfinite(t_cross):
                                hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_cross), 0.0)
                                upflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                ref_u = float(u_p)
                                ref_z = float(z)
                                diag["projection_crossings"] += 1
                                cell_idx += 1
                                if cell_idx >= n_cells:
                                    exit_samples[j] = float(m_dot_samples_prev[j])
                                    exit_u_samples[j] = float(u_p)
                                    break
                                continue
                    h_max = _hamel_delta_neg_max_height(
                        u0=float(ref_u),
                        h0=float(ref_z),
                        a=a2,
                        b=b2,
                        c=c2,
                    )
                    if h_max is not None and h_max < float(z_upper_abs) - 1e-10:
                        t_to_apex_total = _hamel_delta_neg_time_between_velocities(
                            u_start=float(ref_u),
                            u_end=0.0,
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        u_bound_down = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_lower_abs),
                            u0=0.0,
                            h0=float(h_max),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        t_down = _hamel_delta_neg_time_between_velocities(
                            u_start=0.0,
                            u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        if (
                            u_bound_down is not None
                            and t_to_apex_total is not None
                            and t_down is not None
                            and math.isfinite(t_to_apex_total)
                            and math.isfinite(t_down)
                        ):
                            t_cell = max(float(t_to_apex_total), 0.0) + max(float(t_down), 0.0)
                            hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_cell), 0.0)
                            downflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                            total_return_samples[j] += float(m_dot_samples_prev[j])
                            diag["returns"] += 1
                            break
                    elif h_max is not None and h_max >= float(z_upper_abs) - 1e-10:
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=float(ref_u),
                            h0=float(ref_z),
                            a=a2,
                            b=b2,
                            c=c2,
                            u_lo=1e-9,
                            u_hi=max(float(ref_u), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_cross = _hamel_delta_neg_time_between_velocities(
                                u_start=float(ref_u),
                                u_end=float(u_bound_up),
                                a=a2,
                                b=b2,
                                c=c2,
                            )
                            if t_cross is not None and math.isfinite(t_cross):
                                hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_cross), 0.0)
                                upflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                ref_u = float(u_p)
                                ref_z = float(z)
                                diag["projection_crossings"] += 1
                                cell_idx += 1
                                if cell_idx >= n_cells:
                                    exit_samples[j] = float(m_dot_samples_prev[j])
                                    exit_u_samples[j] = float(u_p)
                                    break
                                continue
                if float(u_g) > 1e-9:
                    u_p = max(float(np.nextafter(float(u_g), -math.inf)), 1e-9)
                else:
                    u_p = 1e-9
                continue

            if t_turn is not None and abs(t_evt - t_turn) <= 1e-8 * max(1.0, t_turn):
                a2, b2, c2, _, _, _ = _wirsum_abc_exact_audit(
                    u_g=u_g,
                    u_p=0.0,
                    rho_g=bg.rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=bg.mu_g,
                    phi_s=phi_s,
                    branch="slower",
                )
                u_bound_down = _hamel_delta_neg_boundary_velocity(
                    h_target=float(z_lower_abs),
                    u0=0.0,
                    h0=float(z),
                    a=a2,
                    b=b2,
                    c=c2,
                )
                t_back = _hamel_delta_neg_time_between_velocities(
                    u_start=0.0,
                    u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                    a=a2,
                    b=b2,
                    c=c2,
                )
                if u_bound_down is None or t_back is None or not math.isfinite(t_back):
                    diag["fallbacks"] += 1
                    diag["fallback_delta_neg_time"] += 1
                    break
                hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_back), 0.0)
                downflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                total_return_samples[j] += float(m_dot_samples_prev[j])
                diag["returns"] += 1
                break

    class_transport_by_seg: list[dict[str, np.ndarray]] = []
    hold_time_mean_by_seg: list[float] = []
    hold_time_max_by_seg: list[float] = []
    for i in range(n_cells):
        if n_classes <= 0:
            class_transport_by_seg.append(_empty_class_transport(0))
            hold_time_mean_by_seg.append(0.0)
            hold_time_max_by_seg.append(0.0)
            continue
        positive_mask = np.maximum(m_dot_samples_prev, 0.0) > 1e-12
        hold_time_samples = np.divide(
            hold_samples_by_seg[i],
            np.maximum(m_dot_samples_prev, 1e-12),
            out=np.zeros_like(hold_samples_by_seg[i]),
            where=positive_mask,
        )
        active_times = hold_time_samples[(hold_time_samples > 1e-12) & positive_mask]
        hold_time_mean_by_seg.append(float(np.mean(active_times)) if active_times.size else 0.0)
        hold_time_max_by_seg.append(float(np.max(active_times)) if active_times.size else 0.0)
        class_transport_by_seg.append(
            _aggregate_class_transport(
                n_classes=n_classes,
                class_idx_samples=class_idx_samples,
                hold_samples=hold_samples_by_seg[i],
                upflow_samples=upflow_samples_by_seg[i],
                downflow_samples=downflow_samples_by_seg[i],
                char_frac_samples=char_frac_samples,
                char_conversion_samples=char_conversion_samples,
            )
        )
    return {
        "class_transport_by_seg": class_transport_by_seg,
        "diag": diag,
        "coeff_diag": coeff_diag,
        "hold_time_mean_by_seg": hold_time_mean_by_seg,
        "hold_time_max_by_seg": hold_time_max_by_seg,
        "exit_samples": exit_samples,
        "exit_u_samples": exit_u_samples,
        "return_samples": total_return_samples,
    }


def _simulate_freeboard_global_projected_hamel(
    *,
    N_in: np.ndarray,
    T_in: float,
    P: float,
    D_bed: float,
    H_freeboard: float,
    n_cells: int,
    u_b_bed_top: float,
    d_b_bed_top: float,
    eps_b_bed_top: float,
    eps_d_void_bed_top: float,
    rho_solid_bed_top: float,
    d_p_classes_bed_top: np.ndarray,
    m_char_classes_bed_top: np.ndarray,  # [kg/s]，已是床顶上行流，不再乘 K
    m_ash_classes_bed_top: np.ndarray,  # [kg/s]
    char_conversion_bed_top: float,
    class_aggregate_indices_bed_top: np.ndarray | None,
    n_output_size_classes: int | None,
    fuel_type: str,
    heat_loss_frac: float,
    trajectory_coeff_model: str,
    u_gb_scale: float,
    beta_a_scale: float,
    velocity_sigma: float,
    velocity_bins: int,
    phi_s_bed_top: float,
    z_base_m: float,
    total_height_m: float | None,
    secondary_injection_xi: float | None,
    secondary_O2_mol_s: float,
    secondary_N2_mol_s: float,
    secondary_H2O_mol_s: float,
    secondary_T_K: float,
    secondary_injection_mode: str,
    enabled_reactions: tuple[str, ...] | None,
    inventory_char_sink_enabled: bool,
    heat_loss_mode: str = DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    Q_wall_W: float = 0.0,
    Q_heat_exchanger_W: float = 0.0,
    extent_limiters_enabled: bool = False,
) -> dict:
    A = math.pi * (D_bed**2) / 4.0
    dh = H_freeboard / n_cells
    beta_a = calc_beta_a(max(d_b_bed_top, 1e-6)) * float(max(beta_a_scale, 0.0))
    u_gb0 = calc_u_gb(max(u_b_bed_top, 1e-6), scale=float(max(u_gb_scale, 1e-6)))
    u_p0_mean = 1.53 * max(u_b_bed_top, 1e-6)
    u_p_seed, velocity_weights = _build_velocity_samples(
        u_p0_mean,
        sigma=float(max(velocity_sigma, 0.0)),
        n_bins=int(max(velocity_bins, 1)),
    )
    d_p_class, char_frac_class, m_dot_class0, m_dot_eject0 = _build_size_class_bundles(
        d_p_classes=np.asarray(d_p_classes_bed_top, dtype=np.float64),
        m_char_classes=np.asarray(m_char_classes_bed_top, dtype=np.float64),
        m_ash_classes=np.asarray(m_ash_classes_bed_top, dtype=np.float64),
        rho_p=float(rho_solid_bed_top),
        eps_b=float(eps_b_bed_top),
        eps_d_void=float(eps_d_void_bed_top),
        u_b=float(u_b_bed_top),
        d_b=float(d_b_bed_top),
        area=A,
    )
    aggregate_idx, n_transport_classes, d_p_transport_class = _resolve_transport_class_mapping(
        d_p_class=d_p_class,
        m_dot_class=m_dot_class0,
        class_aggregate_indices=class_aggregate_indices_bed_top,
        n_output_size_classes=n_output_size_classes,
    )
    m_dot_samples0 = np.repeat(m_dot_class0, len(velocity_weights)) * np.tile(velocity_weights, len(m_dot_class0))
    char_frac_samples = np.repeat(char_frac_class, len(velocity_weights))
    char_conversion_samples = np.full_like(m_dot_samples0, float(np.clip(char_conversion_bed_top, 0.0, 1.0 - 1e-9)), dtype=np.float64)
    class_idx_samples = np.repeat(aggregate_idx, len(velocity_weights))
    d_p_samples = np.repeat(d_p_class, len(velocity_weights))
    u_p_samples0 = np.tile(u_p_seed, len(m_dot_class0))
    secondary_N = np.zeros(N_GAS, dtype=np.float64)
    secondary_N[GAS_SPECIES_INDEX["O2"]] = float(max(secondary_O2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["N2"]] = float(max(secondary_N2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["H2O"]] = float(max(secondary_H2O_mol_s, 0.0))

    backgrounds = _build_freeboard_background_profile(
        N_in=N_in,
        T_in=T_in,
        P=P,
        A=A,
        H_freeboard=H_freeboard,
        n_cells=n_cells,
        u_gb0=u_gb0,
        beta_a=beta_a,
        heat_loss_frac=heat_loss_frac,
        fuel_type=fuel_type,
        z_base_m=z_base_m,
        total_height_m=total_height_m,
        secondary_injection_xi=secondary_injection_xi,
        secondary_N=secondary_N,
        secondary_T_K=secondary_T_K,
        secondary_injection_mode=secondary_injection_mode,
        enabled_reactions=enabled_reactions,
        heat_loss_mode=heat_loss_mode,
        Q_wall_W=Q_wall_W,
        Q_heat_exchanger_W=Q_heat_exchanger_W,
        extent_limiters_enabled=extent_limiters_enabled,
    )
    traj = _advance_particle_samples_global_projected_hamel(
        backgrounds=backgrounds,
        rho_p=max(float(rho_solid_bed_top), 1e-6),
        d_p=d_p_samples,
        phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
        u_p_samples_prev=u_p_samples0,
        m_dot_samples_prev=m_dot_samples0,
        char_frac_samples=char_frac_samples,
        char_conversion_samples=char_conversion_samples,
        class_idx_samples=class_idx_samples,
        coeff_model=str(trajectory_coeff_model),
    )
    rho_cat_profile = _estimate_rho_cat_profile_from_projected_transport(
        class_transport_by_seg=traj["class_transport_by_seg"],
        backgrounds=backgrounds,
    )
    if np.any(rho_cat_profile > 0.0):
        backgrounds = _build_freeboard_background_profile(
            N_in=N_in,
            T_in=T_in,
            P=P,
            A=A,
            H_freeboard=H_freeboard,
            n_cells=n_cells,
            u_gb0=u_gb0,
            beta_a=beta_a,
            heat_loss_frac=heat_loss_frac,
            fuel_type=fuel_type,
            z_base_m=z_base_m,
            total_height_m=total_height_m,
            secondary_injection_xi=secondary_injection_xi,
            secondary_N=secondary_N,
            secondary_T_K=secondary_T_K,
            secondary_injection_mode=secondary_injection_mode,
            enabled_reactions=enabled_reactions,
            rho_cat_profile=rho_cat_profile,
            heat_loss_mode=heat_loss_mode,
            Q_wall_W=Q_wall_W,
            Q_heat_exchanger_W=Q_heat_exchanger_W,
        )
        traj = _advance_particle_samples_global_projected_hamel(
            backgrounds=backgrounds,
            rho_p=max(float(rho_solid_bed_top), 1e-6),
            d_p=d_p_samples,
            phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
            u_p_samples_prev=u_p_samples0,
            m_dot_samples_prev=m_dot_samples0,
            char_frac_samples=char_frac_samples,
            char_conversion_samples=char_conversion_samples,
            class_idx_samples=class_idx_samples,
            coeff_model=str(trajectory_coeff_model),
        )

    N = np.array(N_in, dtype=np.float64, copy=True)
    T = float(T_in)
    states: list[FreeboardState] = []
    return_char_by_cell: list[float] = []
    return_ash_by_cell: list[float] = []
    secondary_injection_applied = False
    secondary_injection_segment: int | None = None
    h_total = None if total_height_m is None else float(max(total_height_m, 1e-12))
    xi_inj = None if secondary_injection_xi is None else float(np.clip(secondary_injection_xi, 0.0, 1.0))
    secondary_observation = {
        "z_m": [],
        "xi_reactor": [],
        "slice_global_index": [],
        "segment_index": [],
        "slice_index_within_segment": [],
        "pre_mix_T": [],
        "post_mix_pre_rxn_T": [],
        "post_rxn_T": [],
        "pre_mix_wet_gas": {sp: [] for sp in GAS_SPECIES},
        "post_mix_pre_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
        "post_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
    }

    for i, bg in enumerate(backgrounds):
        class_transport = traj["class_transport_by_seg"][i]
        hold_char = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
        hold_ash = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64)
        x_before = np.asarray(class_transport["char_conversion_classes"], dtype=np.float64).copy()
        up_char = np.asarray(class_transport["m_dot_auf_char_classes"], dtype=np.float64)
        up_ash = np.asarray(class_transport["m_dot_auf_ash_classes"], dtype=np.float64)
        down_char = np.asarray(class_transport["m_dot_ab_char_classes"], dtype=np.float64)
        down_ash = np.asarray(class_transport["m_dot_ab_ash_classes"], dtype=np.float64)
        char_area_total, solid_d_p_eff = _char_area_and_representative_diameter(
            hold_char,
            d_p_transport_class,
            max(float(rho_solid_bed_top), 1e-6),
        )
        char_conversion_seg = _representative_char_conversion(class_transport)
        m_hold_seg = float(np.sum(hold_char + hold_ash))
        return_char_seg = float(np.sum(np.maximum(down_char, 0.0)))
        return_ash_seg = float(np.sum(np.maximum(down_ash, 0.0)))
        return_char_by_cell.append(return_char_seg)
        return_ash_by_cell.append(return_ash_seg)
        m_dot_char = float(np.sum(np.maximum(up_char, 0.0)))
        m_dot_ash = float(np.sum(np.maximum(up_ash, 0.0)))
        m_dot_s_avg = max(m_dot_char + m_dot_ash, 1e-12)
        ash_hold_frac = m_dot_ash / m_dot_s_avg
        rho_cat = (m_hold_seg * ash_hold_frac) / max(bg.V_seg, 1e-12)
        u_t_samples = np.array(
            [
                calc_terminal_velocity_haider(
                    rho_g=bg.rho_g,
                    rho_p=max(float(rho_solid_bed_top), 1e-6),
                    d_p=max(float(dp), 1e-9),
                    phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                    mu_g=bg.mu_g,
                    g_acc=g,
                )
                for dp in d_p_samples
            ],
            dtype=np.float64,
        )
        flow_avg = 0.5 * (np.maximum(up_char + up_ash, 0.0) + np.maximum(down_char + down_ash, 0.0))
        u_p_mean = float(bg.u_gb)
        u_t_mean = float(np.mean(u_t_samples)) if len(u_t_samples) else 0.0
        carry_ratio = float(bg.u_gb / max(u_t_mean, 1e-12))
        xi_seg = (
            float((float(z_base_m) + bg.z_center_m) / h_total)
            if h_total is not None
            else float(bg.z_center_m / max(H_freeboard, 1e-12))
        )
        if (
            not secondary_injection_applied
            and xi_inj is not None
            and float(np.sum(secondary_N)) > 0.0
            and xi_seg >= xi_inj
        ):
            H_mix = (
                calc_gas_enthalpy_flow(N, T, h_cache={})
                + calc_gas_enthalpy_flow(secondary_N, float(secondary_T_K), h_cache={})
            )
            N = N + secondary_N
            T = _solve_T_from_enthalpy(N, H_mix, max(T, float(secondary_T_K)))
            secondary_injection_applied = True
            secondary_injection_segment = i

        pre_mix_T = float(T)
        pre_mix_wet = _wet_gas_profile(N)
        post_mix_pre_rxn_T = float(T)
        post_mix_pre_rxn_wet = _wet_gas_profile(N)
        dt_slice = bg.tau / 12.0
        diag_acc: dict[str, float] = {}
        for k in range(12):
            H_before = calc_gas_enthalpy_flow(N, T, h_cache={})
            N, diag_step = _freeboard_reaction_step(
                N=N,
                T=T,
                P=P,
                V_seg=bg.V_seg,
                dt=dt_slice,
                fuel_type=fuel_type,
                rho_cat=rho_cat,
                char_area_total=char_area_total,
                m_char_total=float(np.sum(np.maximum(hold_char, 0.0))),
                solid_d_p=solid_d_p_eff,
                D_g=gas_diffusivity_correlation(T, P),
                char_conversion=char_conversion_seg,
                enabled_reactions=enabled_reactions,
                extent_limiters_enabled=extent_limiters_enabled,
            )
            H_after = apply_freeboard_substep_gas_enthalpy_heat_loss(
                H_before,
                heat_loss_mode=heat_loss_mode,
                heat_loss_frac_total=heat_loss_frac,
                n_cells=n_cells,
                n_substeps=12,
                n_slices=1,
                Q_wall_W=Q_wall_W,
                Q_heat_exchanger_W=Q_heat_exchanger_W,
            )
            T = _solve_T_from_enthalpy(N, H_after, T)
            if k == 0:
                diag_acc = {name: 0.0 for name in diag_step}
            for name, value in diag_step.items():
                diag_acc[name] += float(value)

        hold_char_before = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64).copy()
        hold_ash_before = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64).copy()
        inventory_sink_diag = (
            _freeboard_inventory_char_sink_diag(
                N=N,
                T=T,
                P=P,
                V_seg=bg.V_seg,
                dt=bg.tau,
                fuel_type=fuel_type,
                char_area_total=char_area_total,
                solid_d_p=solid_d_p_eff,
                D_g=gas_diffusivity_correlation(T, P),
                char_conversion=char_conversion_seg,
                m_char_total=float(np.sum(np.maximum(hold_char_before, 0.0))),
            )
            if inventory_char_sink_enabled
            else {"R1": 0.0, "R2": 0.0, "R3": 0.0, "R4": 0.0}
        )
        class_transport = _apply_char_sink_to_class_transport(class_transport, diag_acc=inventory_sink_diag)
        hold_char = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
        hold_ash = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64)

        secondary_observation["z_m"].append(float(bg.z_center_m))
        secondary_observation["xi_reactor"].append(float(xi_seg))
        secondary_observation["slice_global_index"].append(len(states))
        secondary_observation["segment_index"].append(int(i))
        secondary_observation["slice_index_within_segment"].append(0)
        secondary_observation["pre_mix_T"].append(float(pre_mix_T))
        secondary_observation["post_mix_pre_rxn_T"].append(float(post_mix_pre_rxn_T))
        secondary_observation["post_rxn_T"].append(float(T))
        post_rxn_wet = _wet_gas_profile(N)
        for sp in GAS_SPECIES:
            secondary_observation["pre_mix_wet_gas"][sp].append(float(pre_mix_wet[sp]))
            secondary_observation["post_mix_pre_rxn_wet_gas"][sp].append(float(post_mix_pre_rxn_wet[sp]))
            secondary_observation["post_rxn_wet_gas"][sp].append(float(post_rxn_wet[sp]))

        states.append(
            FreeboardState(
                segment_index=int(i),
                z_center_m=float(bg.z_center_m),
                T=float(T),
                N=N.copy(),
                u0=float(bg.u0),
                u_gb=float(bg.u_gb),
                u_p_mean=float(u_p_mean),
                u_t_mean=float(u_t_mean),
                carry_ratio=float(carry_ratio),
                tau=float(bg.tau),
                m_dot_char=float(m_dot_char),
                m_dot_ash=float(m_dot_ash),
                m_dot_return_char=float(return_char_seg),
                m_dot_return_ash=float(return_ash_seg),
                m_hold_char_before_classes=np.array(hold_char_before, copy=True),
                m_hold_ash_before_classes=np.array(hold_ash_before, copy=True),
                char_conversion_before_classes=np.array(x_before, copy=True),
                char_conversion_classes=np.array(class_transport["char_conversion_classes"], copy=True),
                m_hold_char_classes=np.array(hold_char, copy=True),
                m_hold_ash_classes=np.array(hold_ash, copy=True),
                m_char_sink_applied_classes=np.array(class_transport["m_char_sink_applied_classes"], copy=True),
                m_dot_auf_char_classes=np.array(up_char, copy=True),
                m_dot_auf_ash_classes=np.array(up_ash, copy=True),
                m_dot_ab_char_classes=np.array(down_char, copy=True),
                m_dot_ab_ash_classes=np.array(down_ash, copy=True),
                K_auf_classes=np.array(class_transport["K_auf_classes"], copy=True),
                K_ab_classes=np.array(class_transport["K_ab_classes"], copy=True),
                rho_solid=float(m_hold_seg / max(bg.V_seg, 1e-12)),
                rho_cat=float(rho_cat),
                reaction_diag=diag_acc,
            )
        )

    exit_samples = np.asarray(traj["exit_samples"], dtype=np.float64)
    entrained_exit_char = float(np.sum(exit_samples * char_frac_samples))
    entrained_exit_ash = float(np.sum(exit_samples * (1.0 - char_frac_samples)))
    return {
        "states": states,
        "exit_T": float(T),
        "exit_N": N,
        "exit_gas": _wet_gas_profile(N),
        "beta_A": float(beta_a),
        "trajectory_model": "global_projected_hamel",
        "trajectory_solver": "global_projected_hamel",
        "trajectory_coeff_model": str(trajectory_coeff_model),
        "trajectory_diag": {k: int(v) for k, v in traj["diag"].items()},
        "trajectory_coeff_diag": {
            "re_p_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["re_p_min"]) else traj["coeff_diag"]["re_p_min"]),
            "re_p_max": float(traj["coeff_diag"]["re_p_max"]),
            "cd_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["cd_min"]) else traj["coeff_diag"]["cd_min"]),
            "cd_max": float(traj["coeff_diag"]["cd_max"]),
            "re_p_active_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["re_p_active_min"]) else traj["coeff_diag"]["re_p_active_min"]),
            "re_p_active_max": float(traj["coeff_diag"]["re_p_active_max"]),
            "cd_active_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["cd_active_min"]) else traj["coeff_diag"]["cd_active_min"]),
            "cd_active_max": float(traj["coeff_diag"]["cd_active_max"]),
            "z_drag": float(traj["coeff_diag"]["z_drag"]),
            "exact_delta_pos": int(traj["coeff_diag"]["exact_delta_pos"]),
            "exact_delta_zero": int(traj["coeff_diag"]["exact_delta_zero"]),
            "exact_delta_neg": int(traj["coeff_diag"]["exact_delta_neg"]),
            "exact_a_rel_max": float(traj["coeff_diag"]["exact_a_rel_max"]),
            "exact_b_rel_max": float(traj["coeff_diag"]["exact_b_rel_max"]),
            "exact_c_rel_max": float(traj["coeff_diag"]["exact_c_rel_max"]),
            "exact_dp_grav_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["exact_dp_grav_min"]) else traj["coeff_diag"]["exact_dp_grav_min"]),
            "exact_dp_grav_max": float(traj["coeff_diag"]["exact_dp_grav_max"]),
            "exact_dp_ratio_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["exact_dp_ratio_min"]) else traj["coeff_diag"]["exact_dp_ratio_min"]),
            "exact_dp_ratio_max": float(traj["coeff_diag"]["exact_dp_ratio_max"]),
        },
        "secondary_injection_applied": bool(secondary_injection_applied),
        "secondary_injection_segment": secondary_injection_segment,
        "secondary_local_refine": 1,
        "secondary_injection_mode": str(secondary_injection_mode),
        "entrained_eject_flux_kg_m2_s": float(m_dot_eject0 / max(A, 1e-12)),
        "entrained_eject_char_ash_kg_s": float(m_dot_eject0),
        "entrained_exit_char_kg_s": entrained_exit_char,
        "entrained_exit_ash_kg_s": entrained_exit_ash,
        "entrained_return_char_kg_s": float(bed_face_class_flows_kg_s(return_char_by_cell)[0]) if return_char_by_cell else 0.0,
        "entrained_return_ash_kg_s": float(bed_face_class_flows_kg_s(return_ash_by_cell)[0]) if return_ash_by_cell else 0.0,
        "secondary_observation": secondary_observation,
        "profiles": {
            "T": [float(st.T) for st in states],
            "z_m": [float(st.z_center_m) for st in states],
            "tau_s": [float(st.tau) for st in states],
            "u0_m_s": [float(st.u0) for st in states],
            "u_gb_m_s": [float(st.u_gb) for st in states],
            "u_p_mean_m_s": [float(st.u_p_mean) for st in states],
            "u_t_mean_m_s": [float(st.u_t_mean) for st in states],
            "carry_ratio": [float(st.carry_ratio) for st in states],
            "hold_time_mean_s": [float(v) for v in traj["hold_time_mean_by_seg"]],
            "hold_time_max_s": [float(v) for v in traj["hold_time_max_by_seg"]],
            "reaction_diag": {
                name: [float(st.reaction_diag.get(name, 0.0)) for st in states]
                for name in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R11b", "R11d", "R12")
            },
            "entrained_char_kg_s": [float(st.m_dot_char) for st in states],
            "entrained_ash_kg_s": [float(st.m_dot_ash) for st in states],
            "entrained_return_char_kg_s": [float(st.m_dot_return_char) for st in states],
            "entrained_return_ash_kg_s": [float(st.m_dot_return_ash) for st in states],
            "entrained_solid_density_kg_m3": [float(st.rho_solid) for st in states],
            "entrained_catalyst_density_kg_m3": [float(st.rho_cat) for st in states],
            "solid_holdup_char_before_classes_kg": [np.asarray(st.m_hold_char_before_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_ash_before_classes_kg": [np.asarray(st.m_hold_ash_before_classes, dtype=np.float64).tolist() for st in states],
            "char_conversion_before_classes": [np.asarray(st.char_conversion_before_classes, dtype=np.float64).tolist() for st in states],
            "char_conversion_classes": [np.asarray(st.char_conversion_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_char_classes_kg": [np.asarray(st.m_hold_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_ash_classes_kg": [np.asarray(st.m_hold_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_char_sink_applied_classes_kg": [np.asarray(st.m_char_sink_applied_classes, dtype=np.float64).tolist() for st in states],
            "solid_upflow_char_classes_kg_s": [np.asarray(st.m_dot_auf_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_upflow_ash_classes_kg_s": [np.asarray(st.m_dot_auf_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_downflow_char_classes_kg_s": [np.asarray(st.m_dot_ab_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_downflow_ash_classes_kg_s": [np.asarray(st.m_dot_ab_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_K_auf_classes_1_s": [np.asarray(st.K_auf_classes, dtype=np.float64).tolist() for st in states],
            "solid_K_ab_classes_1_s": [np.asarray(st.K_ab_classes, dtype=np.float64).tolist() for st in states],
            "wet_gas": {
                sp: [
                    float(max(st.N[GAS_SPECIES_INDEX[sp]], 0.0) / max(np.sum(np.maximum(st.N, 0.0)), 1e-12))
                    for st in states
                ]
                for sp in GAS_SPECIES
            },
        },
    }

