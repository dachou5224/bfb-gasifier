"""``scripts/_lu_harness`` 契约：与 ``validation_case_utils`` 直接构建等价。"""

from __future__ import annotations

import sys
from pathlib import Path

_scripts = Path(__file__).resolve().parent.parent / "scripts"
if str(_scripts) not in sys.path:
    sys.path.insert(0, str(_scripts))

from _lu_harness import (
    HAMEL_WIRSUM_CORE,
    LADDER_LEVELS,
    LADDER_P2,
    apply_baseline_kinetics,
    bed_cell_indices,
    build_phase2_config,
    build_phase2_reactor,
    co_co2_score,
    load_case_lu,
    phase1_global_nr_solve_kwargs,
    phase2_solve_kwargs,
    repo_root,
)
from tests.validation_case_utils import (
    PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
)


def test_ladder_p2_matches_convergence_ladder_definition() -> None:
    assert LADDER_P2 == LADDER_LEVELS["p2"]


def test_build_phase2_config_matches_validation_utils() -> None:
    direct = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    via_harness = build_phase2_config(enable_vorab_march_staged=True)
    assert direct.H_bed == via_harness.H_bed
    assert direct.O2_feed == via_harness.O2_feed
    assert bool(direct.vorab_bed_temperature_march_thesis) == bool(
        via_harness.vorab_bed_temperature_march_thesis
    )


def test_phase2_solve_kwargs_matches_validation_utils() -> None:
    assert phase2_solve_kwargs(max_global_iter=7)["max_global_iter"] == 7
    assert phase2_solve_kwargs()["solver"] == PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS["solver"]


def test_build_phase2_reactor_ladder_baseline() -> None:
    reactor = build_phase2_reactor(ladder_level="baseline")
    cfg = reactor.config
    assert cfg.nr_line_search_gas_phase_split_merit_thesis is False
    assert cfg.nr_line_search_energy_merit_thesis is False
    assert cfg.freeboard_trajectory_coeff_model == "stable_mixed_drag_split"
    assert len(bed_cell_indices(reactor)) >= 1


def test_ladder_p2_uses_wirsum_rms_line_search_merit() -> None:
    """Wirsum Eq.2.23：p0/p1/p2 默认关 split/energy max-merit。"""
    for level in ("p0", "p1", "p2"):
        assert LADDER_LEVELS[level]["nr_line_search_gas_phase_split_merit_thesis"] is False
        assert LADDER_LEVELS[level]["nr_line_search_energy_merit_thesis"] is False
    reactor = build_phase2_reactor(ladder_level="p2")
    cfg = reactor.config
    assert cfg.nr_line_search_gas_phase_split_merit_thesis is False
    assert cfg.nr_line_search_energy_merit_thesis is False
    assert cfg.nr_refresh_hydrodynamics_on_accepted_step_thesis is False


def test_hamel_wirsum_core_strips_platform_bridges() -> None:
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=HAMEL_WIRSUM_CORE)
    cfg = reactor.config
    assert cfg.extent_limiters_enabled_thesis is False
    assert cfg.vorab_transport_x0_fast_oxidation_closure_thesis is True
    assert cfg.vorab_transport_x0_char_oxidation_o2_closure_thesis is False
    assert cfg.vorab_bed0_eq24_two_phase_startwert_thesis is True
    assert cfg.vorab_init_solve_bottom_cell_thesis is False
    assert cfg.nr_line_search_gas_phase_split_merit_thesis is False
    assert cfg.nr_line_search_energy_merit_thesis is False
    assert cfg.nr_line_search_fast_oxidation_projection_thesis is False
    assert cfg.nr_line_search_per_phase_fastox_thesis is True
    assert cfg.nr_clip_same_phase_oxidizer_fuel_step_thesis is True
    assert cfg.nr_clip_allow_oxidizer_into_fuel_thesis is False
    assert cfg.vorab_bed0_r1_oxidizer_seed_mol_s_thesis == 2.0
    assert cfg.nr_inner_t_step_default_cap_K_thesis == 40.0
    assert cfg.nr_temperature_fence_bed0_lower_margin_K_thesis == 350.0
    assert cfg.nr_temperature_fence_upper_bed_lower_margin_K_thesis is None
    assert cfg.vorab_bed0_nr_startwert_T_K_thesis == 820.0
    assert cfg.vorab_vm_devolatilization_zone_cells_thesis == 3
    assert cfg.vorab_a_tier_post_staged_h2o_passthrough_thesis is False
    assert cfg.vorab_a_tier_post_staged_syngas_passthrough_thesis is False
    assert cfg.nr_line_search_syngas_collapse_guard_thesis is False
    assert cfg.nr_outer_preinner_upper_co_local2x2_thesis is False
    assert cfg.nr_outer_hot_restart_after_y_co_thesis is False
    assert cfg.nr_outer_preinner_bed_top_solid_holdup_thesis is False
    assert cfg.freeboard_trajectory_coeff_model == "exact_hamel"


def test_co_co2_score_at_reference_is_zero() -> None:
    assert co_co2_score(0.13, 0.11) == 0.0


def test_repo_root_points_at_project() -> None:
    assert (repo_root() / "src" / "core" / "reactor.py").is_file()


def test_apply_baseline_kinetics_sets_unity_scales() -> None:
    cfg = build_phase2_config()
    cfg.r2_scale = 0.5
    apply_baseline_kinetics(cfg)
    assert cfg.r2_scale == 1.0
    assert cfg.r8_scale == 1.0


def test_load_case_lu_has_validation_key_fields() -> None:
    case = load_case_lu()
    assert "secondary_injection_agent" in case
    assert "bed_height_m" in case or "H_bed" in str(case)


def test_phase1_global_nr_solve_kwargs_has_solver() -> None:
    kw = phase1_global_nr_solve_kwargs()
    assert kw["solver"] == "global_nr"
