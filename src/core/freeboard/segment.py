"""simulate_freeboard orchestration.

Ref: Hamel (1999) Ch.3 freeboard particle trajectory; OPT-003 split from freeboard_segment.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.constants import Rg, g
from src.core.energy_balance_defaults import (
    DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    apply_freeboard_substep_gas_enthalpy_heat_loss,
)
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
from src.core.freeboard.reaction import _freeboard_reaction_step
from src.core.freeboard.trajectory_hamel import _simulate_freeboard_global_projected_hamel
from src.core.freeboard.trajectory_surrogate import (
    _advance_particle_samples_force_balance,
    _advance_particle_samples_surrogate,
)
from src.core.freeboard.trajectory_wirsum import _advance_particle_samples_analytical_wirsum
from src.core.freeboard.transport import (
    _aggregate_class_transport,
    _apply_char_sink_to_class_transport,
    bed_face_class_flows_kg_s,
    _build_size_class_bundles,
    _build_velocity_samples,
    _char_area_and_representative_diameter,
    _freeboard_inventory_char_sink_diag,
    _representative_char_conversion,
    _resolve_transport_class_mapping,
    _solve_T_from_enthalpy,
    _wet_gas_profile,
)
from src.core.freeboard.types import FreeboardState

def simulate_freeboard(
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
    char_conversion_bed_top: float = 0.0,
    class_aggregate_indices_bed_top: np.ndarray | None = None,
    n_output_size_classes: int | None = None,
    fuel_type: str,
    heat_loss_frac: float = 0.0,
    heat_loss_mode: str = DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    Q_wall_W: float = 0.0,
    Q_heat_exchanger_W: float = 0.0,
    trajectory_model: str = "analytical_wirsum",
    trajectory_coeff_model: str = "exact_hamel",
    u_gb_scale: float = 1.0,
    beta_a_scale: float = 1.0,
    velocity_sigma: float = 0.60,
    velocity_bins: int = 5,
    phi_s_bed_top: float = 0.86,
    z_base_m: float = 0.0,
    total_height_m: float | None = None,
    secondary_injection_xi: float | None = None,
    secondary_O2_mol_s: float = 0.0,
    secondary_N2_mol_s: float = 0.0,
    secondary_H2O_mol_s: float = 0.0,
    secondary_T_K: float = 293.15,
    secondary_local_refine: int = 1,
    secondary_injection_mode: str = "lumped",
    enabled_reactions: tuple[str, ...] | None = None,
    inventory_char_sink_enabled: bool = False,
    extent_limiters_enabled: bool = False,
) -> dict:
    if H_freeboard <= 0.0 or n_cells <= 0:
        return {
            "states": [],
            "exit_T": float(T_in),
            "exit_N": np.array(N_in, dtype=np.float64, copy=True),
            "exit_gas": {sp: float(v) for sp, v in zip(GAS_SPECIES, np.maximum(N_in, 0.0) / max(np.sum(np.maximum(N_in, 0.0)), 1e-12))},
        }

    if trajectory_model not in {"surrogate_exp", "force_balance", "analytical_wirsum", "global_projected_hamel"}:
        raise ValueError(f"未知的 freeboard trajectory_model: {trajectory_model}")
    if trajectory_coeff_model not in {"stable_mixed_drag_split", "exact_hamel"}:
        raise ValueError(f"未知的 analytical_wirsum coeff_model: {trajectory_coeff_model}")
    if secondary_injection_mode not in {"lumped", "distributed_uniform"}:
        raise ValueError(f"未知的 secondary_injection_mode: {secondary_injection_mode}")

    if trajectory_model == "global_projected_hamel":
        return _simulate_freeboard_global_projected_hamel(
            N_in=N_in,
            T_in=T_in,
            P=P,
            D_bed=D_bed,
            H_freeboard=H_freeboard,
            n_cells=n_cells,
            u_b_bed_top=u_b_bed_top,
            d_b_bed_top=d_b_bed_top,
            eps_b_bed_top=eps_b_bed_top,
            eps_d_void_bed_top=eps_d_void_bed_top,
            rho_solid_bed_top=rho_solid_bed_top,
            d_p_classes_bed_top=d_p_classes_bed_top,
            m_char_classes_bed_top=m_char_classes_bed_top,
            m_ash_classes_bed_top=m_ash_classes_bed_top,
            char_conversion_bed_top=char_conversion_bed_top,
            class_aggregate_indices_bed_top=class_aggregate_indices_bed_top,
            n_output_size_classes=n_output_size_classes,
            fuel_type=fuel_type,
            heat_loss_frac=heat_loss_frac,
            heat_loss_mode=heat_loss_mode,
            Q_wall_W=Q_wall_W,
            Q_heat_exchanger_W=Q_heat_exchanger_W,
            trajectory_coeff_model=trajectory_coeff_model,
            u_gb_scale=u_gb_scale,
            beta_a_scale=beta_a_scale,
            velocity_sigma=velocity_sigma,
            velocity_bins=velocity_bins,
            phi_s_bed_top=phi_s_bed_top,
            z_base_m=z_base_m,
            total_height_m=total_height_m,
            secondary_injection_xi=secondary_injection_xi,
            secondary_O2_mol_s=secondary_O2_mol_s,
            secondary_N2_mol_s=secondary_N2_mol_s,
            secondary_H2O_mol_s=secondary_H2O_mol_s,
            secondary_T_K=secondary_T_K,
            secondary_injection_mode=secondary_injection_mode,
            enabled_reactions=enabled_reactions,
            inventory_char_sink_enabled=inventory_char_sink_enabled,
            extent_limiters_enabled=extent_limiters_enabled,
        )

    A = math.pi * (D_bed**2) / 4.0
    dh = H_freeboard / n_cells
    beta_a = calc_beta_a(max(d_b_bed_top, 1e-6)) * float(max(beta_a_scale, 0.0))
    u_gb0 = calc_u_gb(max(u_b_bed_top, 1e-6), scale=float(max(u_gb_scale, 1e-6)))
    u_gb_prev = float(u_gb0)
    # Hamel Eq. 3.73/3.76: ghost-bubble gas speed and particle launch speed
    # are distinct. u_gb,0 = u_b,ws while mean particle launch speed is
    # u_p,0 ~= 1.53 * u_b,ws.
    u_p0_mean = 1.53 * max(u_b_bed_top, 1e-6)
    u_p_samples_prev, velocity_weights = _build_velocity_samples(
        u_p0_mean,
        sigma=float(max(velocity_sigma, 0.0)),
        n_bins=int(max(velocity_bins, 1)),
    )
    N = np.array(N_in, dtype=np.float64, copy=True)
    T = float(T_in)
    states: list[FreeboardState] = []
    d_p_classes = np.asarray(d_p_classes_bed_top, dtype=np.float64)
    m_char_classes = np.asarray(m_char_classes_bed_top, dtype=np.float64)
    m_ash_classes = np.asarray(m_ash_classes_bed_top, dtype=np.float64)
    d_p_class, char_frac_class, m_dot_class0, m_dot_eject0 = _build_size_class_bundles(
        d_p_classes=d_p_classes,
        m_char_classes=m_char_classes,
        m_ash_classes=m_ash_classes,
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
    u_p_samples_prev = np.tile(u_p_samples_prev, len(m_dot_class0))
    slower_ref_u_prev = np.array(u_p_samples_prev, dtype=np.float64, copy=True)
    slower_ref_z_prev = np.zeros_like(u_p_samples_prev, dtype=np.float64)
    f0 = m_dot_eject0 / max(A, 1e-12)
    return_char_by_cell: list[float] = []
    return_ash_by_cell: list[float] = []
    secondary_injection_applied = False
    secondary_injection_segment: int | None = None
    secondary_N = np.zeros(N_GAS, dtype=np.float64)
    secondary_N[GAS_SPECIES_INDEX["O2"]] = float(max(secondary_O2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["N2"]] = float(max(secondary_N2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["H2O"]] = float(max(secondary_H2O_mol_s, 0.0))
    xi_inj = None if secondary_injection_xi is None else float(np.clip(secondary_injection_xi, 0.0, 1.0))
    h_total = None if total_height_m is None else float(max(total_height_m, 1e-12))
    refine_local = int(max(secondary_local_refine, 1))
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
    trajectory_diag = {
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
    }
    trajectory_coeff_diag = {
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

    for i in range(n_cells):
        n_tot = max(float(np.sum(np.maximum(N, 0.0))), 1e-12)
        u0 = (n_tot * Rg * T / P) / max(A, 1e-12)
        u_gb_in = float(u_gb_prev if i > 0 else u_gb0)
        exp_arg_out = float(np.clip(-beta_a * dh, -200.0, 200.0))
        exp_arg_center = float(np.clip(-beta_a * 0.5 * dh, -200.0, 200.0))
        u_gb_out = u0 + (u_gb_in - u0) * math.exp(exp_arg_out)
        u_gb = u0 + (u_gb_in - u0) * math.exp(exp_arg_center)
        z_center = (i + 0.5) * dh
        xi_local = z_center / max(H_freeboard, 1e-12)
        xi_seg = xi_local if h_total is None else float((float(z_base_m) + z_center) / h_total)
        tau = dh / max(0.5 * (u_gb_in + u_gb_out), 1e-6)
        V_seg = A * dh
        dt = tau / 12.0

        injection_applied_here = False
        if (
            not secondary_injection_applied
            and xi_inj is not None
            and float(np.sum(secondary_N)) > 0.0
            and xi_seg >= xi_inj
        ):
            secondary_injection_applied = True
            injection_applied_here = True
            secondary_injection_segment = i

        y_map = {sp: float(max(N[GAS_SPECIES_INDEX[sp]], 0.0) / n_tot) for sp in GAS_SPECIES}
        rho_g = gas_density_ideal(P, T, y_map)
        mu_g = gas_viscosity_power_law(T, 1.8e-5)
        u_t_samples = np.array(
            [
                calc_terminal_velocity_haider(
                    rho_g=rho_g,
                    rho_p=max(float(rho_solid_bed_top), 1e-6),
                    d_p=max(float(dp), 1e-9),
                    phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                    mu_g=mu_g,
                    g_acc=g,
                )
                for dp in d_p_samples
            ],
            dtype=np.float64,
        )
        if trajectory_model == "analytical_wirsum":
            u_p_samples, m_dot_samples, hold_samples, return_samples, slower_ref_u, slower_ref_z, diag_seg = _advance_particle_samples_analytical_wirsum(
                dh=dh,
                z_lower_abs=float(i * dh),
                u_g_prev=u_gb_in,
                u_g_next=u_gb_out,
                rho_g=rho_g,
                rho_p=max(float(rho_solid_bed_top), 1e-6),
                phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                mu_g=mu_g,
                u_p_samples_prev=u_p_samples_prev,
                m_dot_samples_prev=m_dot_samples0,
                d_p=d_p_samples,
                char_frac_samples=char_frac_samples,
                slower_ref_u_prev=slower_ref_u_prev,
                slower_ref_z_prev=slower_ref_z_prev,
                coeff_model=str(trajectory_coeff_model),
            )
            for key, val in diag_seg.items():
                if key in trajectory_diag:
                    trajectory_diag[key] += int(val)
                elif key in trajectory_coeff_diag:
                    if key.endswith("_min"):
                        trajectory_coeff_diag[key] = min(float(trajectory_coeff_diag[key]), float(val))
                    elif key.endswith("_max"):
                        trajectory_coeff_diag[key] = max(float(trajectory_coeff_diag[key]), float(val))
                    elif key.startswith("exact_delta_"):
                        trajectory_coeff_diag[key] += int(val)
                    else:
                        trajectory_coeff_diag[key] = float(val)
        elif trajectory_model == "force_balance":
            u_p_samples, m_dot_samples, hold_samples, return_samples = _advance_particle_samples_force_balance(
                dh=dh,
                u_g_prev=u_gb_in,
                u_g_next=u_gb_out,
                rho_g=rho_g,
                rho_p=max(float(rho_solid_bed_top), 1e-6),
                phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                mu_g=mu_g,
                u_p_samples_prev=u_p_samples_prev,
                m_dot_samples_prev=m_dot_samples0,
                d_p=d_p_samples,
                char_frac_samples=char_frac_samples,
            )
            trajectory_diag["fallbacks"] += 0
        else:
            u_p_samples, m_dot_samples, hold_samples, return_samples = _advance_particle_samples_surrogate(
                dh=dh,
                beta_a=beta_a,
                u_g_asym=u0,
                u_g_prev=u_gb_in,
                u_p_samples_prev=u_p_samples_prev,
                m_dot_samples_prev=m_dot_samples0,
            )
        m_hold_seg = float(np.sum(hold_samples))
        return_char_seg = float(np.sum(np.maximum(return_samples, 0.0) * char_frac_samples))
        return_ash_seg = float(np.sum(np.maximum(return_samples, 0.0) * (1.0 - char_frac_samples)))
        class_transport = _aggregate_class_transport(
            n_classes=n_transport_classes,
            class_idx_samples=class_idx_samples,
            hold_samples=hold_samples,
            upflow_samples=m_dot_samples,
            downflow_samples=return_samples,
            char_frac_samples=char_frac_samples,
            char_conversion_samples=char_conversion_samples,
        )
        char_area_total, solid_d_p_eff = _char_area_and_representative_diameter(
            np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64),
            d_p_transport_class,
            max(float(rho_solid_bed_top), 1e-6),
        )
        char_conversion_seg = _representative_char_conversion(class_transport)
        flow_avg = 0.5 * (m_dot_samples0 + m_dot_samples)
        m_dot_s_avg = float(np.sum(flow_avg))
        rho_solid_seg = m_hold_seg / max(V_seg, 1e-12)
        m_dot_char = float(np.sum(flow_avg * char_frac_samples))
        m_dot_ash = float(np.sum(flow_avg * (1.0 - char_frac_samples)))
        ash_hold_frac = m_dot_ash / max(m_dot_s_avg, 1e-12)
        rho_cat = (m_hold_seg * ash_hold_frac) / max(V_seg, 1e-12)
        u_p_mean = float(np.sum(flow_avg * np.maximum(u_p_samples, 0.0)) / max(m_dot_s_avg, 1e-12))
        u_t_mean = float(np.sum(flow_avg * u_t_samples) / max(m_dot_s_avg, 1e-12))
        carry_ratio = float(u_gb / max(u_t_mean, 1e-12))
        return_char_by_cell.append(float(return_char_seg))
        return_ash_by_cell.append(float(return_ash_seg))

        n_slices = refine_local if injection_applied_here else 1
        dh_slice = dh / float(n_slices)
        tau_slice = tau / float(n_slices)
        v_slice = A * dh_slice
        secondary_slice_N = secondary_N / float(n_slices) if (injection_applied_here and secondary_injection_mode == "distributed_uniform") else secondary_N

        for j in range(n_slices):
            slice_z_center = i * dh + (j + 0.5) * dh_slice
            slice_xi = (
                float((float(z_base_m) + slice_z_center) / h_total)
                if h_total is not None
                else float(slice_z_center / max(H_freeboard, 1e-12))
            )
            pre_mix_T = float(T)
            pre_mix_wet = _wet_gas_profile(N)
            mixed_here = False
            if injection_applied_here:
                if secondary_injection_mode == "lumped":
                    if j == 0:
                        H_mix = (
                            calc_gas_enthalpy_flow(N, T, h_cache={})
                            + calc_gas_enthalpy_flow(secondary_N, float(secondary_T_K), h_cache={})
                        )
                        N = N + secondary_N
                        T = _solve_T_from_enthalpy(N, H_mix, max(T, float(secondary_T_K)))
                        mixed_here = True
                else:
                    H_mix = (
                        calc_gas_enthalpy_flow(N, T, h_cache={})
                        + calc_gas_enthalpy_flow(secondary_slice_N, float(secondary_T_K), h_cache={})
                    )
                    N = N + secondary_slice_N
                    T = _solve_T_from_enthalpy(N, H_mix, max(T, float(secondary_T_K)))
                    mixed_here = True
            post_mix_pre_rxn_T = float(T)
            post_mix_pre_rxn_wet = _wet_gas_profile(N)
            dt_slice = tau_slice / 12.0
            diag_acc = {}
            for k in range(12):
                H_before = calc_gas_enthalpy_flow(N, T, h_cache={})
                N, diag_step = _freeboard_reaction_step(
                    N=N,
                    T=T,
                    P=P,
                    V_seg=v_slice,
                    dt=dt_slice,
                    fuel_type=fuel_type,
                    rho_cat=rho_cat,
                    char_area_total=char_area_total,
                    m_char_total=float(
                        np.sum(np.maximum(class_transport["m_hold_char_classes"], 0.0))
                    ),
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
                    n_slices=n_slices,
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
            x_before = np.asarray(class_transport["char_conversion_classes"], dtype=np.float64).copy()
            inventory_sink_diag = (
                _freeboard_inventory_char_sink_diag(
                    N=N,
                    T=T,
                    P=P,
                    V_seg=v_slice,
                    dt=tau_slice,
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

            if mixed_here:
                secondary_observation["z_m"].append(float(slice_z_center))
                secondary_observation["xi_reactor"].append(float(slice_xi))
                secondary_observation["slice_global_index"].append(len(states))
                secondary_observation["segment_index"].append(int(i))
                secondary_observation["slice_index_within_segment"].append(int(j))
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
                    z_center_m=slice_z_center,
                    T=float(T),
                    N=N.copy(),
                    u0=float(u0),
                    u_gb=float(u_gb),
                    u_p_mean=float(u_p_mean),
                    u_t_mean=float(u_t_mean),
                    carry_ratio=float(carry_ratio),
                    tau=float(tau_slice),
                    m_dot_char=float(m_dot_char),
                    m_dot_ash=float(m_dot_ash),
                    m_dot_return_char=float(return_char_seg / float(n_slices)),
                    m_dot_return_ash=float(return_ash_seg / float(n_slices)),
                    m_hold_char_before_classes=np.array(hold_char_before, copy=True),
                    m_hold_ash_before_classes=np.array(hold_ash_before, copy=True),
                    char_conversion_before_classes=np.array(x_before, copy=True),
                    char_conversion_classes=np.array(class_transport["char_conversion_classes"], copy=True),
                    m_hold_char_classes=np.array(class_transport["m_hold_char_classes"], copy=True),
                    m_hold_ash_classes=np.array(class_transport["m_hold_ash_classes"], copy=True),
                    m_char_sink_applied_classes=np.array(class_transport["m_char_sink_applied_classes"], copy=True),
                    m_dot_auf_char_classes=np.array(class_transport["m_dot_auf_char_classes"], copy=True),
                    m_dot_auf_ash_classes=np.array(class_transport["m_dot_auf_ash_classes"], copy=True),
                    m_dot_ab_char_classes=np.array(class_transport["m_dot_ab_char_classes"], copy=True),
                    m_dot_ab_ash_classes=np.array(class_transport["m_dot_ab_ash_classes"], copy=True),
                    K_auf_classes=np.array(class_transport["K_auf_classes"], copy=True),
                    K_ab_classes=np.array(class_transport["K_ab_classes"], copy=True),
                    rho_solid=float(rho_solid_seg),
                    rho_cat=float(rho_cat),
                    reaction_diag=diag_acc,
                )
            )
        u_p_samples_prev = u_p_samples
        m_dot_samples0 = m_dot_samples
        if trajectory_model == "analytical_wirsum":
            slower_ref_u_prev = slower_ref_u
            slower_ref_z_prev = slower_ref_z
        u_gb_prev = float(u_gb_out)

    exit_gas = _wet_gas_profile(N)
    entrained_exit_total = float(np.sum(m_dot_samples0))
    entrained_exit_char = float(np.sum(m_dot_samples0 * char_frac_samples))
    entrained_exit_ash = float(np.sum(m_dot_samples0 * (1.0 - char_frac_samples)))
    return {
        "states": states,
        "exit_T": float(T),
        "exit_N": N,
        "exit_gas": exit_gas,
        "beta_A": float(beta_a),
        "trajectory_model": str(trajectory_model),
        "trajectory_solver": (
            ("wirsum_analytical_exact_hamel_coeffs" if trajectory_coeff_model == "exact_hamel" else "wirsum_analytical")
            if trajectory_model == "analytical_wirsum"
            else ("force_balance_numeric" if trajectory_model == "force_balance" else "surrogate_exp")
        ),
        "trajectory_coeff_model": (
            str(trajectory_coeff_model)
            if trajectory_model == "analytical_wirsum"
            else ("full_force_balance" if trajectory_model == "force_balance" else "surrogate_decay")
        ),
        "trajectory_diag": {k: int(v) for k, v in trajectory_diag.items()},
        "trajectory_coeff_diag": {
            "re_p_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["re_p_min"]) else trajectory_coeff_diag["re_p_min"]),
            "re_p_max": float(trajectory_coeff_diag["re_p_max"]),
            "cd_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["cd_min"]) else trajectory_coeff_diag["cd_min"]),
            "cd_max": float(trajectory_coeff_diag["cd_max"]),
            "re_p_active_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["re_p_active_min"]) else trajectory_coeff_diag["re_p_active_min"]),
            "re_p_active_max": float(trajectory_coeff_diag["re_p_active_max"]),
            "cd_active_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["cd_active_min"]) else trajectory_coeff_diag["cd_active_min"]),
            "cd_active_max": float(trajectory_coeff_diag["cd_active_max"]),
            "z_drag": float(trajectory_coeff_diag["z_drag"]),
            "exact_delta_pos": int(trajectory_coeff_diag["exact_delta_pos"]),
            "exact_delta_zero": int(trajectory_coeff_diag["exact_delta_zero"]),
            "exact_delta_neg": int(trajectory_coeff_diag["exact_delta_neg"]),
            "exact_a_rel_max": float(trajectory_coeff_diag["exact_a_rel_max"]),
            "exact_b_rel_max": float(trajectory_coeff_diag["exact_b_rel_max"]),
            "exact_c_rel_max": float(trajectory_coeff_diag["exact_c_rel_max"]),
            "exact_dp_grav_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["exact_dp_grav_min"]) else trajectory_coeff_diag["exact_dp_grav_min"]),
            "exact_dp_grav_max": float(trajectory_coeff_diag["exact_dp_grav_max"]),
            "exact_dp_ratio_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["exact_dp_ratio_min"]) else trajectory_coeff_diag["exact_dp_ratio_min"]),
            "exact_dp_ratio_max": float(trajectory_coeff_diag["exact_dp_ratio_max"]),
        },
        "secondary_injection_applied": bool(secondary_injection_applied),
        "secondary_injection_segment": secondary_injection_segment,
        "secondary_local_refine": int(refine_local),
        "secondary_injection_mode": str(secondary_injection_mode),
        "entrained_eject_flux_kg_m2_s": float(f0),
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
