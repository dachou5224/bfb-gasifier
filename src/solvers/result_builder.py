"""Result assembly helpers for Reactor NR orchestration.

Keep `src/core/reactor.py` focused on flow orchestration (Config -> Vorabrechnung -> Cell Model -> NR),
while this module builds reporting payloads from solved states.
"""

from __future__ import annotations

import numpy as np

from src.core.composition import wet_to_dry_mole_fractions
from src.core.connectivity import build_thesis_connectivity_topology
from src.core.cell_kinetics import build_reaction_sources
from src.core.freeboard_bridge import bed_top_entrained_size_proxy, effective_bed_top_entrained_d_p_classes
from src.core.reaction_numbering import map_impl_profile_diag_to_thesis, reaction_numbering_metadata
from src.core.species import GAS_SPECIES
from src.core.cell import S_ASH, S_CHAR
from src.core.freeboard_segment import simulate_freeboard


def compute_carbon_conversion(reactor, cell_solid_outflow_component_fn) -> float:
    bot = reactor.cells[0]
    top = reactor.cells[-1]
    m_char_in = float(
        np.sum(
            np.maximum(
                bot.m_solid_zu[:, S_CHAR] + bot.m_solid_in[:, S_CHAR],
                0.0,
            )
        )
    )
    if reactor._use_explicit_side_block_cells() and reactor.cyclone_cell is not None:
        m_char_out = float(cell_solid_outflow_component_fn(reactor.cyclone_cell, S_CHAR, direction="up"))
    else:
        top_up_char = float(cell_solid_outflow_component_fn(top, S_CHAR, direction="up"))
        recycle_char = float(np.sum(np.maximum(bot.m_solid_rez[:, S_CHAR], 0.0)))
        m_char_out = max(top_up_char - recycle_char, 0.0)
    carbon_conv = 1.0 - (m_char_out / max(m_char_in, 1e-12))
    return float(np.clip(carbon_conv, 0.0, 1.0))


def refresh_explicit_freeboard_closure_for_result(reactor) -> dict:
    """Report explicit freeboard closure freshness without mutating NR state."""
    if not reactor._use_explicit_freeboard_solver_graph():
        return {
            "freeboard_result_closure_refreshed": False,
            "freeboard_result_closure_refresh_reason": "not_explicit_freeboard_graph",
        }
    if not getattr(reactor, "freeboard_cells", None):
        return {
            "freeboard_result_closure_refreshed": False,
            "freeboard_result_closure_refresh_reason": "no_freeboard_cells",
        }

    current_bed_top_initial = float(np.sum(np.maximum(reactor.cells[-1]._solid_upflow_rates()[:, S_CHAR], 0.0)))
    previous_fb = getattr(reactor, "_last_explicit_freeboard_closure", None)
    previous_bed_top = (
        float(previous_fb.get("bed_top_up_char_kg_s", current_bed_top_initial))
        if isinstance(previous_fb, dict)
        else current_bed_top_initial
    )
    # Do not refresh/sync the closure here. The NR solver has just produced a
    # self-consistent state vector; mutating explicit freeboard/side solid states
    # during result assembly invalidates the final residual metrics and can turn
    # side-element solid DOFs into apparent post-solve residual spikes.
    final_bed_top = current_bed_top_initial
    closure_gap = current_bed_top_initial - previous_bed_top
    return {
        "freeboard_result_closure_refreshed": False,
        "freeboard_result_boundary_refreshed": False,
        "freeboard_result_closure_refresh_reason": "diagnostic_only_no_state_mutation",
        "freeboard_pre_result_refresh_bed_top_up_char_kg_s": previous_bed_top,
        "freeboard_pre_result_refresh_current_bed_top_up_char_kg_s": current_bed_top_initial,
        "freeboard_current_bed_top_up_char_kg_s": final_bed_top,
        "freeboard_pre_result_refresh_gap_current_bed_top_kg_s": closure_gap,
        "freeboard_post_result_refresh_bed_top_up_char_kg_s": previous_bed_top,
        "freeboard_post_result_refresh_gap_current_bed_top_kg_s": closure_gap,
        "freeboard_result_closure_refresh_iters": 0,
    }


def build_exit_summary(
    reactor,
    resolve_axial_heat_loss_distribution_fn,
    cell_solid_outflow_component_fn,
) -> dict:
    cfg = reactor.config
    top = reactor.cells[-1]
    freeboard_bed_top_proxy = bed_top_entrained_size_proxy(reactor)
    freeboard_d_p_eff_bed_top = np.array(freeboard_bed_top_proxy["d_p_eff"], dtype=np.float64, copy=True)
    freeboard_bed_top_char_conversion = float(freeboard_bed_top_proxy["bed_spm_shrink_proxy"])
    _, freeboard_loss = resolve_axial_heat_loss_distribution_fn(cfg)
    y_bed = top._mole_fractions("combined")
    bed_exit_gas = {sp: float(y_bed[j]) for j, sp in enumerate(GAS_SPECIES)}
    bed_exit_N = np.maximum(top.N_d + top.N_b, 0.0)
    bed_z = [float(c.geo.h_center) for c in reactor.cells]
    bed_T = [float(c.T) for c in reactor.cells]
    bed_u0 = [float(c.u0) for c in reactor.cells]
    bed_eps_b = [float(c.eps_b) for c in reactor.cells]
    bed_eps_d_void = [float(c.eps_d_voidage) for c in reactor.cells]
    bed_bulk_solid = [float((1.0 - c.eps_b) * (1.0 - c.eps_d_voidage)) for c in reactor.cells]
    total_height = float(cfg.H_bed + max(cfg.H_freeboard, 0.0))
    reactor_bed_char_inventory = float(
        sum(float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0))) for c in reactor.cells)
    )
    reactor_visible_char_inventory = float(
        sum(float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0))) for c in reactor.cells)
        + sum(float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0))) for c in reactor.freeboard_cells)
        + (
            float(np.sum(np.maximum(reactor.cyclone_cell.m_solid[:, S_CHAR], 0.0)))
            if reactor.cyclone_cell is not None
            else 0.0
        )
        + (
            float(np.sum(np.maximum(reactor.return_leg_cell.m_solid[:, S_CHAR], 0.0)))
            if reactor.return_leg_cell is not None
            else 0.0
        )
    )

    summary = {
        "bed_T_profile": bed_T,
        "bed_u0_profile_m_s": bed_u0,
        "bed_eps_b_profile": bed_eps_b,
        "bed_eps_d_void_profile": bed_eps_d_void,
        "bed_bulk_solid_fraction_profile": bed_bulk_solid,
        "bed_exit_gas": bed_exit_gas,
        "bed_exit_gas_dry": wet_to_dry_mole_fractions(bed_exit_gas),
        "bed_axial_z_m": bed_z,
        "bed_axial_xi": [float(z / max(cfg.H_bed, 1e-12)) for z in bed_z],
        "reaction_numbering": reaction_numbering_metadata(),
        "thesis_connectivity_topology": None,
        "side_block_active": bool(reactor._use_explicit_side_block_cells()),
    }
    bot = reactor.cells[0]
    summary.update(
        {
            "bed0_energy_residual_W": float(bot.calc_energy_balance()),
            "bed0_recycle_solid_char_kg_s": float(np.sum(np.maximum(bot.m_solid_rez[:, S_CHAR], 0.0))),
            "bed0_recycle_solid_ash_kg_s": float(np.sum(np.maximum(bot.m_solid_rez[:, S_ASH], 0.0))),
            "bed0_recycle_solid_T_K": float(bot.T_rez_solid),
            "bed0_recycle_solid_enthalpy_W": float(bot._calc_solid_enthalpy_flow(bot.m_solid_rez, bot.T_rez_solid)),
            "bed0_fresh_solid_enthalpy_W": float(bot._calc_solid_enthalpy_flow(bot.m_solid_zu, bot.T_zu_solid)),
            "bed_top_downflow_solid_enthalpy_W": float(top._calc_solid_enthalpy_flow(top.m_solid_ab_in, top.T_solid_ab_in)),
            "bed_top_downflow_solid_T_K": float(top.T_solid_ab_in),
        }
    )

    def _explicit_freeboard_reaction_bundle(cell):
        if getattr(cell, "cell_type", "") != "freeboard":
            return None
        if float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) <= 0.0:
            return None
        thermo = cell._get_local_thermo_bundle()
        areas = cell._calc_char_surface_area_per_class()
        if not bool(getattr(cfg, "freeboard_explicit_enable_char_hetero", True)):
            areas = np.zeros_like(areas)
        return build_reaction_sources(
            T=cell.T,
            P=cell.P,
            fuel_type=cell.fuel_type,
            V_b=cell.V_b,
            V_d=cell.V_d,
            C_b=thermo["C_b"],
            C_d=thermo["C_d"],
            y_b=thermo["y_b"],
            y_d=thermo["y_d"],
            gas_src_vm=cell._vm_gas_source_cache,
            solid_sink_vm=cell._vm_solid_sink_cache,
            areas=areas,
            solid_d_p=cell.solid.d_p,
            D_g=float(thermo["D_g"]),
            char_conversion=cell._compute_char_conversion(),
            rho_cat=cell._catalyst_bulk_density(),
            enable_r12=cell.enable_r12,
            use_gibbs_minor=cell.use_gibbs_minor,
            gibbs_minor_sources=None,
            r4_scale=cell.r4_scale,
            r5_scale=cell.r5_scale,
            r6_scale=cell.r6_scale,
            r7_scale=cell.r7_scale,
            rate_multiplier=1.0,
            N_zu_d=cell.N_zu_d,
            N_d_in=cell.N_d_in,
            N_zu_b=cell.N_zu_b,
            N_b_in=cell.N_b_in,
            N_rez_d=cell.N_rez_d,
            N_rez_b=cell.N_rez_b,
            N_ex=cell.N_ex,
            solid_shape=cell.R_solid.shape,
            char_index=S_CHAR,
        )

    if cfg.H_freeboard > 0.0 and int(cfg.n_freeboard_cells) > 0:
        fb = reactor._last_explicit_freeboard_closure if reactor._use_explicit_freeboard_solver_graph() else None
        enabled_reactions = tuple(cfg.freeboard_enabled_reactions)
        if not bool(getattr(cfg, "freeboard_closure_enable_char_hetero", True)):
            enabled_reactions = tuple(name for name in enabled_reactions if name not in {"R1", "R2", "R3", "R4"})
        if fb is None:
            fb = simulate_freeboard(
                N_in=bed_exit_N,
                T_in=float(top.T),
                P=cfg.P,
                D_bed=cfg.D_bed,
                H_freeboard=float(cfg.H_freeboard),
                n_cells=int(cfg.n_freeboard_cells),
                u_b_bed_top=float(top.u_b),
                d_b_bed_top=float(top.d_b),
                eps_b_bed_top=float(top.eps_b),
                eps_d_void_bed_top=float(top.eps_d_voidage),
                rho_solid_bed_top=float(top.solid.rho_s),
                d_p_classes_bed_top=np.array(freeboard_d_p_eff_bed_top, dtype=np.float64, copy=True),
                m_char_classes_bed_top=np.maximum(top._solid_upflow_rates()[:, S_CHAR], 0.0),
                m_ash_classes_bed_top=np.maximum(top._solid_upflow_rates()[:, S_ASH], 0.0),
                char_conversion_bed_top=float(freeboard_bed_top_proxy["bed_total_char_conversion_proxy"]),
                fuel_type=cfg.fuel_type,
                heat_loss_frac=float(freeboard_loss),
                trajectory_model=str(cfg.freeboard_trajectory_model),
                trajectory_coeff_model=str(cfg.freeboard_trajectory_coeff_model),
                u_gb_scale=float(max(cfg.freeboard_u_gb_scale, 1e-6)),
                beta_a_scale=float(max(cfg.freeboard_beta_a_scale, 0.0)),
                velocity_sigma=float(max(cfg.freeboard_velocity_sigma, 0.0)),
                velocity_bins=int(max(cfg.freeboard_velocity_bins, 1)),
                phi_s_bed_top=float(top.solid.phi_s),
                z_base_m=float(cfg.H_bed),
                total_height_m=total_height,
                secondary_injection_xi=cfg.freeboard_secondary_injection_xi,
                secondary_O2_mol_s=float(max(cfg.freeboard_secondary_O2_mol_s, 0.0)),
                secondary_N2_mol_s=float(max(cfg.freeboard_secondary_N2_mol_s, 0.0)),
                secondary_H2O_mol_s=float(max(cfg.freeboard_secondary_H2O_mol_s, 0.0)),
                secondary_T_K=float(max(cfg.freeboard_secondary_T_K, 1.0)),
                secondary_local_refine=int(max(cfg.freeboard_secondary_local_refine, 1)),
                secondary_injection_mode=str(cfg.freeboard_secondary_injection_mode),
                enabled_reactions=enabled_reactions,
                inventory_char_sink_enabled=bool(getattr(cfg, "freeboard_explicit_enable_char_hetero", True)),
            )
            if not reactor._use_explicit_freeboard_solver_graph():
                reactor._sync_freeboard_cells_from_closure(fb)
                reactor._initialize_explicit_side_block_states()
        if reactor._use_explicit_freeboard_solver_graph() and reactor.freeboard_cells:
            entrained_char_profile = [
                cell_solid_outflow_component_fn(cell, S_CHAR, direction="up") for cell in reactor.freeboard_cells
            ]
            entrained_ash_profile = [
                cell_solid_outflow_component_fn(cell, S_ASH, direction="up") for cell in reactor.freeboard_cells
            ]
            entrained_return_char_profile = [
                cell_solid_outflow_component_fn(cell, S_CHAR, direction="down") for cell in reactor.freeboard_cells
            ]
            entrained_return_ash_profile = [
                cell_solid_outflow_component_fn(cell, S_ASH, direction="down") for cell in reactor.freeboard_cells
            ]
            entrained_eject_char_ash = float(fb["entrained_eject_char_ash_kg_s"])
            entrained_exit_char = float(entrained_char_profile[-1])
            entrained_exit_ash = float(entrained_ash_profile[-1])
            entrained_return_char = float(entrained_return_char_profile[0])
            entrained_return_ash = float(entrained_return_ash_profile[0])
            if reactor.cyclone_cell is not None:
                cyclone_capture_char = cell_solid_outflow_component_fn(reactor.cyclone_cell, S_CHAR, direction="down")
                cyclone_capture_ash = cell_solid_outflow_component_fn(reactor.cyclone_cell, S_ASH, direction="down")
            else:
                cyclone_capture_char = 0.0
                cyclone_capture_ash = 0.0
        else:
            entrained_char_profile = fb["profiles"]["entrained_char_kg_s"]
            entrained_ash_profile = fb["profiles"]["entrained_ash_kg_s"]
            entrained_return_char_profile = fb["profiles"]["entrained_return_char_kg_s"]
            entrained_return_ash_profile = fb["profiles"]["entrained_return_ash_kg_s"]
            entrained_eject_char_ash = float(fb["entrained_eject_char_ash_kg_s"])
            entrained_exit_char = float(fb["entrained_exit_char_kg_s"])
            entrained_exit_ash = float(fb["entrained_exit_ash_kg_s"])
            entrained_return_char = float(fb["entrained_return_char_kg_s"])
            entrained_return_ash = float(fb["entrained_return_ash_kg_s"])
            cyclone_capture_char = entrained_exit_char * float(np.clip(cfg.freeboard_cyclone_capture_char_frac, 0.0, 1.0))
            cyclone_capture_ash = entrained_exit_ash * float(np.clip(cfg.freeboard_cyclone_capture_ash_frac, 0.0, 1.0))
        use_explicit_freeboard_profiles = reactor._use_explicit_freeboard_cells() and len(reactor.freeboard_cells) == len(
            fb["states"]
        )
        if use_explicit_freeboard_profiles:
            freeboard_z = [float(c.geo.h_center) for c in reactor.freeboard_cells]
            freeboard_T = [float(c.T) for c in reactor.freeboard_cells]
            freeboard_u0 = [float(c.u0) for c in reactor.freeboard_cells]
            freeboard_eps_b = [float(c.eps_b) for c in reactor.freeboard_cells]
            freeboard_eps_d_void = [float(c.eps_d_voidage) for c in reactor.freeboard_cells]
            freeboard_hold_char_before = [
                float(np.sum(np.maximum(np.asarray(st.m_hold_char_before_classes, dtype=np.float64), 0.0))) for st in fb["states"]
            ]
            freeboard_hold_ash_before = [
                float(np.sum(np.maximum(np.asarray(st.m_hold_ash_before_classes, dtype=np.float64), 0.0))) for st in fb["states"]
            ]
            freeboard_hold_char = [
                float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0))) for c in reactor.freeboard_cells
            ]
            freeboard_hold_ash = [
                float(np.sum(np.maximum(c.m_solid[:, S_ASH], 0.0))) for c in reactor.freeboard_cells
            ]
            freeboard_explicit_char_sink_applied = [
                float(np.sum(np.maximum(np.asarray(getattr(c, "freeboard_explicit_char_sink_applied_kg", 0.0), dtype=np.float64), 0.0)))
                for c in reactor.freeboard_cells
            ]
            explicit_bundles = [_explicit_freeboard_reaction_bundle(c) for c in reactor.freeboard_cells]
            freeboard_char_rxn = [
                float(bundle.net_molar_gas_source_char) if bundle is not None else 0.0
                for bundle in explicit_bundles
            ]
            freeboard_explicit_reaction_diag_impl = {
                "R1": [float(bundle.extent_r1) if bundle is not None else 0.0 for bundle in explicit_bundles],
                "R2": [float(bundle.extent_r2) if bundle is not None else 0.0 for bundle in explicit_bundles],
                "R3": [float(bundle.extent_r3) if bundle is not None else 0.0 for bundle in explicit_bundles],
                "R4": [float(bundle.extent_r4) if bundle is not None else 0.0 for bundle in explicit_bundles],
            }
            freeboard_transport_char_delta = [after - before for before, after in zip(freeboard_hold_char_before, freeboard_hold_char)]
        else:
            freeboard_z = [float(cfg.H_bed + st.z_center_m) for st in fb["states"]]
            freeboard_T = [float(st.T) for st in fb["states"]]
            freeboard_u0 = fb["profiles"]["u0_m_s"]
            freeboard_eps_b = [0.0] * len(fb["states"])
            freeboard_eps_d_void = [1.0] * len(fb["states"])
            freeboard_hold_char = [
                float(np.sum(np.maximum(np.asarray(st.m_hold_char_classes, dtype=np.float64), 0.0))) for st in fb["states"]
            ]
            freeboard_hold_ash = [
                float(np.sum(np.maximum(np.asarray(st.m_hold_ash_classes, dtype=np.float64), 0.0))) for st in fb["states"]
            ]
            freeboard_hold_char_before = [
                float(np.sum(np.maximum(np.asarray(st.m_hold_char_before_classes, dtype=np.float64), 0.0))) for st in fb["states"]
            ]
            freeboard_hold_ash_before = [
                float(np.sum(np.maximum(np.asarray(st.m_hold_ash_before_classes, dtype=np.float64), 0.0))) for st in fb["states"]
            ]
            freeboard_explicit_char_sink_applied = [0.0] * len(fb["states"])
            freeboard_char_rxn = [0.0] * len(fb["states"])
            freeboard_explicit_reaction_diag_impl = {"R1": [0.0] * len(fb["states"]), "R2": [0.0] * len(fb["states"]), "R3": [0.0] * len(fb["states"]), "R4": [0.0] * len(fb["states"])}
            freeboard_transport_char_delta = [after - before for before, after in zip(freeboard_hold_char_before, freeboard_hold_char)]
        freeboard_xi = [float(z / max(total_height, 1e-12)) for z in freeboard_z]
        freeboard_equiv_char_res_time = [
            float(after / max(up, 1e-12))
            for after, up in zip(freeboard_hold_char, entrained_char_profile)
        ]
        freeboard_equiv_char_res_time_before = [
            float(before / max(up, 1e-12))
            for before, up in zip(freeboard_hold_char_before, entrained_char_profile)
        ]
        freeboard_hold_time_mean = fb["profiles"].get("hold_time_mean_s", [0.0] * len(freeboard_T))
        freeboard_hold_time_max = fb["profiles"].get("hold_time_max_s", [0.0] * len(freeboard_T))
        freeboard_hold_time_to_tau_mean = [
            float(hmean / max(tau, 1e-12))
            for hmean, tau in zip(freeboard_hold_time_mean, fb["profiles"]["tau_s"])
        ]
        freeboard_hold_time_to_tau_max = [
            float(hmax / max(tau, 1e-12))
            for hmax, tau in zip(freeboard_hold_time_max, fb["profiles"]["tau_s"])
        ]
        freeboard_char_holdup_to_visible_inventory = [
            float(after / max(reactor_visible_char_inventory, 1e-12))
            for after in freeboard_hold_char
        ]
        freeboard_char_holdup_before_to_visible_inventory = [
            float(before / max(reactor_visible_char_inventory, 1e-12))
            for before in freeboard_hold_char_before
        ]
        freeboard_char_holdup_to_bed_inventory = [
            float(after / max(reactor_bed_char_inventory, 1e-12))
            for after in freeboard_hold_char
        ]
        freeboard_char_holdup_before_to_bed_inventory = [
            float(before / max(reactor_bed_char_inventory, 1e-12))
            for before in freeboard_hold_char_before
        ]
        trajectory_coeff_diag = dict(fb.get("trajectory_coeff_diag", {}))
        d_p_input_arr = np.asarray(fb.get("bed_top_d_p_input_m", freeboard_bed_top_proxy["d_p_input"]), dtype=np.float64)
        d_p_eff_arr = np.asarray(fb.get("bed_top_d_p_eff_m", freeboard_d_p_eff_bed_top), dtype=np.float64)
        d_p_eff_min = float(np.min(d_p_eff_arr)) if d_p_eff_arr.size else 0.0
        d_p_eff_max = float(np.max(d_p_eff_arr)) if d_p_eff_arr.size else 0.0
        d_p_input_min = float(np.min(d_p_input_arr)) if d_p_input_arr.size else 0.0
        exact_dp_grav_max = float(trajectory_coeff_diag.get("exact_dp_grav_max", 0.0) or 0.0)
        exact_dp_grav_min = float(trajectory_coeff_diag.get("exact_dp_grav_min", 0.0) or 0.0)
        required_dp_scale = (
            float(exact_dp_grav_max / d_p_eff_min)
            if d_p_eff_min > 0.0 and exact_dp_grav_max > 0.0
            else 0.0
        )
        required_spm_conversion = (
            float(np.clip(1.0 - (exact_dp_grav_max / d_p_input_min) ** 3, 0.0, 1.0))
            if d_p_input_min > 0.0 and exact_dp_grav_max > 0.0
            else 0.0
        )
        summary.update(
            {
                "freeboard_active": True,
                "freeboard_closure_mode": (
                    "explicit_cells_with_external_closure"
                    if reactor._use_explicit_freeboard_solver_graph()
                    else (
                        "explicit_cells_with_external_closure"
                        if use_explicit_freeboard_profiles
                        else (
                            "external_closure_with_explicit_cells"
                            if reactor._use_explicit_freeboard_cells()
                            else "external_closure_only"
                        )
                    )
                ),
                "freeboard_beta_A": float(fb["beta_A"]),
                "freeboard_trajectory_model": str(fb["trajectory_model"]),
                "freeboard_trajectory_solver": str(fb.get("trajectory_solver", "unknown")),
                "freeboard_trajectory_coeff_model": str(fb.get("trajectory_coeff_model", "unknown")),
                "freeboard_trajectory_diag": dict(fb.get("trajectory_diag", {})),
                "freeboard_trajectory_coeff_diag": trajectory_coeff_diag,
                "freeboard_exact_dp_grav_min_m": exact_dp_grav_min,
                "freeboard_exact_dp_grav_max_m": exact_dp_grav_max,
                "freeboard_bed_top_d_p_eff_min_m": d_p_eff_min,
                "freeboard_bed_top_d_p_eff_max_m": d_p_eff_max,
                "freeboard_required_dp_scale_for_smallest_eff_to_grav": required_dp_scale,
                "freeboard_required_spm_conversion_for_smallest_input_to_grav": required_spm_conversion,
                "freeboard_secondary_injection_applied": bool(fb["secondary_injection_applied"]),
                "freeboard_secondary_injection_segment": fb["secondary_injection_segment"],
                "freeboard_secondary_local_refine": int(fb["secondary_local_refine"]),
                "freeboard_secondary_injection_mode": str(fb["secondary_injection_mode"]),
                "freeboard_secondary_observation_z_m": [float(z) for z in fb["secondary_observation"]["z_m"]],
                "freeboard_secondary_observation_xi": [float(x) for x in fb["secondary_observation"]["xi_reactor"]],
                "freeboard_secondary_observation_slice_global_index": [
                    int(j) for j in fb["secondary_observation"]["slice_global_index"]
                ],
                "freeboard_secondary_observation_segment_index": [
                    int(j) for j in fb["secondary_observation"]["segment_index"]
                ],
                "freeboard_secondary_observation_slice_index_within_segment": [
                    int(j) for j in fb["secondary_observation"]["slice_index_within_segment"]
                ],
                "freeboard_secondary_observation_pre_mix_T": [float(v) for v in fb["secondary_observation"]["pre_mix_T"]],
                "freeboard_secondary_observation_post_mix_pre_rxn_T": [
                    float(v) for v in fb["secondary_observation"]["post_mix_pre_rxn_T"]
                ],
                "freeboard_secondary_observation_post_rxn_T": [float(v) for v in fb["secondary_observation"]["post_rxn_T"]],
                "freeboard_secondary_observation_pre_mix_wet_gas": {
                    sp: [float(v) for v in fb["secondary_observation"]["pre_mix_wet_gas"][sp]] for sp in GAS_SPECIES
                },
                "freeboard_secondary_observation_post_mix_pre_rxn_wet_gas": {
                    sp: [float(v) for v in fb["secondary_observation"]["post_mix_pre_rxn_wet_gas"][sp]] for sp in GAS_SPECIES
                },
                "freeboard_secondary_observation_post_rxn_wet_gas": {
                    sp: [float(v) for v in fb["secondary_observation"]["post_rxn_wet_gas"][sp]] for sp in GAS_SPECIES
                },
                "freeboard_entrained_eject_flux_kg_m2_s": float(fb["entrained_eject_flux_kg_m2_s"]),
                "freeboard_entrained_eject_char_ash_kg_s": entrained_eject_char_ash,
                "freeboard_bed_top_char_conversion": float(fb.get("bed_top_char_conversion", freeboard_bed_top_char_conversion)),
                "freeboard_bed_top_char_conversion_used": float(fb.get("bed_top_char_conversion_used", freeboard_bed_top_proxy["bed_total_char_conversion_proxy"])),
                "freeboard_bed_top_char_conversion_local": float(fb.get("bed_top_char_conversion_local", freeboard_bed_top_proxy["top_local_char_conversion"])),
                "freeboard_bed_top_char_conversion_proxy_total": float(fb.get("bed_top_char_conversion_proxy_total", freeboard_bed_top_proxy["bed_total_char_conversion_proxy"])),
                "freeboard_bed_top_combustion_share_proxy": float(fb.get("bed_top_combustion_share_proxy", freeboard_bed_top_proxy["bed_combustion_share_proxy"])),
                "freeboard_bed_top_combustion_conversion_proxy": float(fb.get("bed_top_combustion_conversion_proxy", freeboard_bed_top_proxy["bed_combustion_conversion_proxy"])),
                "freeboard_bed_top_gasification_age_proxy": float(fb.get("bed_top_gasification_age_proxy", freeboard_bed_top_proxy["bed_gasification_age_proxy"])),
                "freeboard_bed_top_age_quadrature_diag": dict(fb.get("bed_top_age_quadrature_diag", {})),
                "freeboard_bed_top_r1_char_consumption_kg_s": float(fb.get("bed_top_r1_char_consumption_kg_s", freeboard_bed_top_proxy["bed_r1_char_consumption_kg_s"])),
                "freeboard_bed_top_hetero_char_consumption_kg_s": float(fb.get("bed_top_hetero_char_consumption_kg_s", freeboard_bed_top_proxy["bed_hetero_char_consumption_kg_s"])),
                "freeboard_bed_top_up_char_kg_s": float(fb.get("bed_top_up_char_kg_s", freeboard_bed_top_proxy["bed_top_up_char_kg_s"])),
                "freeboard_bed_top_d_p_input_m": [
                    float(v) for v in d_p_input_arr
                ],
                "freeboard_bed_top_d_p_eff_m": [
                    float(v) for v in d_p_eff_arr
                ],
                "freeboard_entrained_exit_char_kg_s": entrained_exit_char,
                "freeboard_entrained_exit_ash_kg_s": entrained_exit_ash,
                "freeboard_entrained_return_char_kg_s": entrained_return_char,
                "freeboard_entrained_return_ash_kg_s": entrained_return_ash,
                "freeboard_cyclone_capture_char_kg_s": cyclone_capture_char,
                "freeboard_cyclone_capture_ash_kg_s": cyclone_capture_ash,
                "freeboard_cyclone_recycle_candidate_char_ash_kg_s": cyclone_capture_char + cyclone_capture_ash,
                "freeboard_T_profile": freeboard_T,
                "freeboard_axial_z_m": freeboard_z,
                "freeboard_axial_xi": freeboard_xi,
                "freeboard_tau_profile_s": [float(st.tau) for st in fb["states"]],
                "freeboard_u0_profile_m_s": freeboard_u0,
                "freeboard_eps_b_profile": freeboard_eps_b,
                "freeboard_eps_d_void_profile": freeboard_eps_d_void,
                "freeboard_solid_holdup_char_profile_kg": freeboard_hold_char,
                "freeboard_solid_holdup_ash_profile_kg": freeboard_hold_ash,
                "freeboard_solid_holdup_char_before_profile_kg": freeboard_hold_char_before,
                "freeboard_solid_holdup_ash_before_profile_kg": freeboard_hold_ash_before,
                "freeboard_transport_char_delta_profile_kg": freeboard_transport_char_delta,
                "reactor_bed_char_inventory_kg": reactor_bed_char_inventory,
                "reactor_visible_char_inventory_kg": reactor_visible_char_inventory,
                "freeboard_char_holdup_to_bed_inventory_ratio_profile": freeboard_char_holdup_to_bed_inventory,
                "freeboard_char_holdup_before_to_bed_inventory_ratio_profile": freeboard_char_holdup_before_to_bed_inventory,
                "freeboard_char_holdup_to_visible_inventory_ratio_profile": freeboard_char_holdup_to_visible_inventory,
                "freeboard_char_holdup_before_to_visible_inventory_ratio_profile": freeboard_char_holdup_before_to_visible_inventory,
                "freeboard_equiv_char_residence_time_before_profile_s": freeboard_equiv_char_res_time_before,
                "freeboard_equiv_char_residence_time_profile_s": freeboard_equiv_char_res_time,
                "freeboard_explicit_char_sink_applied_profile_kg": freeboard_explicit_char_sink_applied,
                "freeboard_char_reaction_source_profile_mol_s": freeboard_char_rxn,
                "freeboard_u_gb_profile_m_s": fb["profiles"]["u_gb_m_s"],
                "freeboard_u_p_mean_profile_m_s": fb["profiles"]["u_p_mean_m_s"],
                "freeboard_u_t_mean_profile_m_s": fb["profiles"]["u_t_mean_m_s"],
                "freeboard_carry_ratio_profile": fb["profiles"]["carry_ratio"],
                "freeboard_hold_time_mean_profile_s": freeboard_hold_time_mean,
                "freeboard_hold_time_max_profile_s": freeboard_hold_time_max,
                "freeboard_hold_time_to_tau_mean_ratio_profile": freeboard_hold_time_to_tau_mean,
                "freeboard_hold_time_to_tau_max_ratio_profile": freeboard_hold_time_to_tau_max,
                "freeboard_gas_profiles_wet": fb["profiles"]["wet_gas"],
                "freeboard_reaction_diag": fb["profiles"]["reaction_diag"],
                "freeboard_reaction_diag_impl": fb["profiles"]["reaction_diag"],
                "freeboard_reaction_diag_thesis": map_impl_profile_diag_to_thesis(fb["profiles"]["reaction_diag"]),
                "freeboard_reaction_diag_labeling": "impl_ids",
                "freeboard_explicit_reaction_diag_impl": freeboard_explicit_reaction_diag_impl,
                "freeboard_entrained_char_kg_s": entrained_char_profile,
                "freeboard_entrained_ash_kg_s": entrained_ash_profile,
                "freeboard_entrained_return_char_profile_kg_s": entrained_return_char_profile,
                "freeboard_entrained_return_ash_profile_kg_s": entrained_return_ash_profile,
                "freeboard_entrained_solid_density_kg_m3": fb["profiles"]["entrained_solid_density_kg_m3"],
                "freeboard_entrained_catalyst_density_kg_m3": fb["profiles"]["entrained_catalyst_density_kg_m3"],
                "T_profile": bed_T + freeboard_T,
                "axial_z_m": bed_z + freeboard_z,
                "axial_xi_reactor": [float(z / max(total_height, 1e-12)) for z in (bed_z + freeboard_z)],
                "exit_gas": (
                    {sp: float(reactor.cyclone_cell._mole_fractions("combined")[j]) for j, sp in enumerate(GAS_SPECIES)}
                    if reactor._use_explicit_freeboard_solver_graph() and reactor.cyclone_cell is not None
                    else fb["exit_gas"]
                ),
                "exit_gas_dry": wet_to_dry_mole_fractions(
                    {sp: float(reactor.cyclone_cell._mole_fractions("combined")[j]) for j, sp in enumerate(GAS_SPECIES)}
                    if reactor._use_explicit_freeboard_solver_graph() and reactor.cyclone_cell is not None
                    else fb["exit_gas"]
                ),
                "reactor_exit_T": float(reactor.cyclone_cell.T)
                if reactor._use_explicit_freeboard_solver_graph() and reactor.cyclone_cell is not None
                else float(fb["exit_T"]),
            }
        )
    else:
        exit_gas = bed_exit_gas
        reactor_exit_T = float(top.T)
        if reactor._use_explicit_side_block_cells() and reactor.cyclone_cell is not None and reactor.return_leg_cell is not None:
            y_cyc = reactor.cyclone_cell._mole_fractions("combined")
            exit_gas = {sp: float(y_cyc[j]) for j, sp in enumerate(GAS_SPECIES)}
            reactor_exit_T = float(reactor.cyclone_cell.T)
        summary.update(
            {
                "freeboard_active": False,
                "freeboard_closure_mode": None,
                "freeboard_beta_A": None,
                "freeboard_trajectory_model": None,
                "freeboard_trajectory_solver": None,
                "freeboard_trajectory_coeff_model": None,
                "freeboard_trajectory_diag": {},
                "freeboard_trajectory_coeff_diag": {},
                "freeboard_exact_dp_grav_min_m": 0.0,
                "freeboard_exact_dp_grav_max_m": 0.0,
                "freeboard_bed_top_d_p_eff_min_m": 0.0,
                "freeboard_bed_top_d_p_eff_max_m": 0.0,
                "freeboard_required_dp_scale_for_smallest_eff_to_grav": 0.0,
                "freeboard_required_spm_conversion_for_smallest_input_to_grav": 0.0,
                "freeboard_secondary_injection_applied": False,
                "freeboard_secondary_injection_segment": None,
                "freeboard_secondary_local_refine": 1,
                "freeboard_secondary_injection_mode": "lumped",
                "freeboard_secondary_observation_z_m": [],
                "freeboard_secondary_observation_xi": [],
                "freeboard_secondary_observation_slice_global_index": [],
                "freeboard_secondary_observation_segment_index": [],
                "freeboard_secondary_observation_slice_index_within_segment": [],
                "freeboard_secondary_observation_pre_mix_T": [],
                "freeboard_secondary_observation_post_mix_pre_rxn_T": [],
                "freeboard_secondary_observation_post_rxn_T": [],
                "freeboard_secondary_observation_pre_mix_wet_gas": {sp: [] for sp in GAS_SPECIES},
                "freeboard_secondary_observation_post_mix_pre_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
                "freeboard_secondary_observation_post_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
                "freeboard_entrained_eject_flux_kg_m2_s": 0.0,
                "freeboard_entrained_eject_char_ash_kg_s": 0.0,
                "freeboard_entrained_exit_char_kg_s": 0.0,
                "freeboard_entrained_exit_ash_kg_s": 0.0,
                "freeboard_entrained_return_char_kg_s": 0.0,
                "freeboard_entrained_return_ash_kg_s": 0.0,
                "freeboard_cyclone_capture_char_kg_s": 0.0,
                "freeboard_cyclone_capture_ash_kg_s": 0.0,
                "freeboard_cyclone_recycle_candidate_char_ash_kg_s": 0.0,
                "freeboard_T_profile": [],
                "freeboard_axial_z_m": [],
                "freeboard_axial_xi": [],
                "freeboard_tau_profile_s": [],
                "freeboard_u0_profile_m_s": [],
                "freeboard_eps_b_profile": [],
                "freeboard_eps_d_void_profile": [],
                "freeboard_solid_holdup_char_profile_kg": [],
                "freeboard_solid_holdup_ash_profile_kg": [],
                "freeboard_char_reaction_source_profile_mol_s": [],
                "freeboard_u_gb_profile_m_s": [],
                "freeboard_u_p_mean_profile_m_s": [],
                "freeboard_u_t_mean_profile_m_s": [],
                "freeboard_carry_ratio_profile": [],
                "freeboard_gas_profiles_wet": {sp: [] for sp in GAS_SPECIES},
                "freeboard_reaction_diag": {name: [] for name in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R11b", "R11d", "R12")},
                "freeboard_reaction_diag_impl": {
                    name: [] for name in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R11b", "R11d", "R12")
                },
                "freeboard_reaction_diag_thesis": {name: [] for name in ("R5", "R6", "R7", "R8", "R9", "R10", "R11")},
                "freeboard_reaction_diag_labeling": "impl_ids",
                "freeboard_explicit_reaction_diag_impl": {name: [] for name in ("R1", "R2", "R3", "R4")},
                "freeboard_entrained_char_kg_s": [],
                "freeboard_entrained_ash_kg_s": [],
                "freeboard_entrained_return_char_profile_kg_s": [],
                "freeboard_entrained_return_ash_profile_kg_s": [],
                "freeboard_entrained_solid_density_kg_m3": [],
                "freeboard_entrained_catalyst_density_kg_m3": [],
                "T_profile": bed_T,
                "axial_z_m": bed_z,
                "axial_xi_reactor": [float(z / max(total_height, 1e-12)) for z in bed_z],
                "exit_gas": exit_gas,
                "exit_gas_dry": wet_to_dry_mole_fractions(exit_gas),
                "reactor_exit_T": reactor_exit_T,
            }
        )
    if reactor._use_explicit_side_block_cells() and reactor.cyclone_cell is not None and reactor.return_leg_cell is not None:
        cyclone_char_flow_up = cell_solid_outflow_component_fn(reactor.cyclone_cell, S_CHAR, direction="up")
        return_leg_char_flow_down = cell_solid_outflow_component_fn(reactor.return_leg_cell, S_CHAR, direction="down")
        summary["cyclone_T"] = float(reactor.cyclone_cell.T)
        summary["return_leg_T"] = float(reactor.return_leg_cell.T)
        # Prefer explicit *_flow_kg_s names; keep legacy inventory keys for compatibility.
        summary["cyclone_char_flow_up_kg_s"] = cyclone_char_flow_up
        summary["return_leg_char_flow_down_kg_s"] = return_leg_char_flow_down
        summary["cyclone_char_inventory_kg_s"] = cyclone_char_flow_up
        summary["return_leg_char_inventory_kg_s"] = return_leg_char_flow_down
    if bool(cfg.thesis_mode):
        summary["thesis_connectivity_topology"] = build_thesis_connectivity_topology(
            n_bed_cells=int(cfg.n_cells),
            n_freeboard_cells=int(cfg.n_freeboard_cells) if summary["freeboard_active"] else 0,
            H_bed=float(cfg.H_bed),
            H_freeboard=float(cfg.H_freeboard) if summary["freeboard_active"] else 0.0,
            recycle_gas=bool(cfg.recycle_gas),
            recirculation_frac=float(cfg.recirculation_frac),
            cyclone_capture_char_frac=float(cfg.freeboard_cyclone_capture_char_frac),
            cyclone_capture_ash_frac=float(cfg.freeboard_cyclone_capture_ash_frac),
            secondary_injection_applied=bool(summary["freeboard_secondary_injection_applied"]),
            secondary_injection_segment=summary["freeboard_secondary_injection_segment"],
            explicit_side_blocks=bool(reactor._use_explicit_side_block_cells()),
        )
    return summary


def finalize_global_nr_result(
    reactor,
    nr_result: dict,
    resolve_axial_heat_loss_distribution_fn,
    cell_solid_outflow_component_fn,
) -> dict:
    freeboard_result_refresh = refresh_explicit_freeboard_closure_for_result(reactor)
    carbon_conv = compute_carbon_conversion(reactor, cell_solid_outflow_component_fn)
    exit_summary = build_exit_summary(
        reactor,
        resolve_axial_heat_loss_distribution_fn=resolve_axial_heat_loss_distribution_fn,
        cell_solid_outflow_component_fn=cell_solid_outflow_component_fn,
    )
    return {
        "converged": nr_result["converged"],
        "converged_outer": nr_result.get("converged_outer"),
        "converged_inner_nr": nr_result.get("converged_inner_nr"),
        "converged_fully": nr_result.get("converged_outer") and nr_result.get("converged_inner_nr"),
        "rms_scaled_final": nr_result.get("rms_scaled_final"),
        "rms_scaled_gas_final": nr_result.get("rms_scaled_gas_final"),
        "rms_scaled_solid_final": nr_result.get("rms_scaled_solid_final"),
        "rms_scaled_energy_final": nr_result.get("rms_scaled_energy_final"),
        "rms_scaled_gas_combined_final": nr_result.get("rms_scaled_gas_combined_final"),
        "rms_scaled_gas_phase_split_final": nr_result.get("rms_scaled_gas_phase_split_final"),
        "rms_scaled_component_max_final": nr_result.get("rms_scaled_component_max_final"),
        "max_abs_scaled_gas_final": nr_result.get("max_abs_scaled_gas_final"),
        "max_abs_scaled_solid_final": nr_result.get("max_abs_scaled_solid_final"),
        "max_abs_scaled_energy_final": nr_result.get("max_abs_scaled_energy_final"),
        "max_abs_scaled_gas_combined_final": nr_result.get("max_abs_scaled_gas_combined_final"),
        "max_abs_scaled_gas_phase_split_final": nr_result.get("max_abs_scaled_gas_phase_split_final"),
        "max_abs_scaled_final": nr_result.get("max_abs_scaled_final"),
        "n_iter": nr_result["n_iter"],
        "residual": nr_result.get("residual", 0.0),
        "norm_history": nr_result.get("norm_history", []),
        "nr_init_strategy": nr_result.get("nr_init_strategy"),
        "nr_init_s_total": nr_result.get("nr_init_s_total"),
        "nr_vorabrechnung_s": nr_result.get("nr_vorabrechnung_s"),
        "nr_vorabrechnung_policy": nr_result.get("nr_vorabrechnung_policy"),
        "nr_gs_warmup_steps": nr_result.get("nr_gs_warmup_steps", 0),
        "nr_outer_max": nr_result.get("nr_outer_max"),
        "nr_outer_iters": nr_result.get("nr_outer_iters"),
        "nr_outer_history": nr_result.get("nr_outer_history"),
        "nr_vorabrechnung_signatures": nr_result.get("nr_vorabrechnung_signatures", []),
        "nr_inner_tol_rms": nr_result.get("nr_inner_tol_rms"),
        "nr_jacobian_strategy": nr_result.get("nr_jacobian_strategy", nr_result.get("jacobian_strategy")),
        "nr_layout_audit": nr_result.get("nr_layout_audit", []),
        "nr_jacobian_structure": nr_result.get("jacobian_structure"),
        "nr_linear_solver_backend": nr_result.get("nr_linear_solver_backend", nr_result.get("linear_solver_backend")),
        "nr_linear_solver_backend_last": nr_result.get("linear_solver_backend_last"),
        "nr_side_element_count": (nr_result.get("jacobian_structure") or {}).get("side_element_count"),
        "nr_band_block_count": (nr_result.get("jacobian_structure") or {}).get("band_block_count"),
        "nr_structure_validation_ok": (nr_result.get("jacobian_structure") or {}).get("structure_validation_ok"),
        "nr_linear_fallback_used": bool((nr_result.get("counts") or {}).get("structured_fallbacks", 0)),
        "nr_schur_size": (nr_result.get("counts") or {}).get("nr_schur_size_last"),
        "nr_jacobian_lag_steps": nr_result.get("nr_jacobian_lag_steps", 1),
        "nr_reaction_continuation_active": nr_result.get("nr_reaction_continuation_active", False),
        "nr_reaction_continuation_stages": nr_result.get("nr_reaction_continuation_stages", [1.0]),
        "nr_reaction_continuation_target_cells": nr_result.get("nr_reaction_continuation_target_cells", 0),
        "nr_gs_warmup_s": nr_result.get("nr_gs_warmup_s", 0.0),
        "nr_outer_refresh_s_total": nr_result.get("nr_outer_refresh_s_total"),
        "nr_inner_solve_s_total": nr_result.get("nr_inner_solve_s_total"),
        "nr_inner_budget_total": nr_result.get("nr_inner_budget_total"),
        "nr_inner_budget_used": nr_result.get("nr_inner_budget_used"),
        "nr_total_s": nr_result.get("nr_total_s"),
        "nr_timing": nr_result.get("timing"),
        "nr_counts": nr_result.get("counts"),
        "nr_accepted_lambda_history": nr_result.get("accepted_lambda_history", []),
        "nr_line_search_trial_counts": nr_result.get("line_search_trial_counts", []),
        "nr_clip_history": nr_result.get("clip_history", []),
        "carbon_conv": carbon_conv,
        **exit_summary,
        **freeboard_result_refresh,
    }
