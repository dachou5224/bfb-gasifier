"""阶段 A：两条基线必须分路径标注，不能把 toy 核写成生产已修。"""

from __future__ import annotations

from src.solvers.wirsum_1997_baseline_contract import (
    build_phase2_init_contract,
    build_reduced_3cell_contract,
)


def test_reduced_fixture_is_x0_only_and_not_thesis() -> None:
    snap = build_reduced_3cell_contract()
    assert snap["path"] == "reduced_3cell_x0"
    assert snap["init_kind"] == "generate_initial_x0_only"
    assert snap["not_a_correctness_claim"] is True
    assert snap["knobs"]["thesis_mode"] is False
    assert snap["n_bed_cells"] == 3
    assert snap["n_freeboard_cells"] == 0
    assert snap["units"]["gas_holdup"] == "mol/s"
    assert snap["units"]["solid_holdup"] == "kg"
    assert snap["nr_solid_bed"][0] == ["char", "ash"]
    assert snap["vm_cache_valid"] == [False, False, False]
    assert snap["thesis_vorab_sources_single_shot"] is False
    assert snap["vorab_outer_source_policy"] == "hamel_outer_refresh"
    assert snap["vorab_outer_source_policy_configured"] == "auto"
    assert snap["outer_may_refresh_drying_pyro"] is True
    assert snap["vorab_tau_model"] == "engineering_H_over_umf"
    assert snap["freeboard_trajectory"]["outer_updates_freeboard_trajectory"] is False


def test_phase2_init_is_not_the_reduced_fixture() -> None:
    snap = build_phase2_init_contract()
    assert snap["path"] == "phase2_3bed_2fb_init"
    assert snap["init_kind"] == "run_init_and_precalc_for_global_nr"
    assert snap["knobs"]["thesis_mode"] is True
    assert snap["n_bed_cells"] == 3
    assert snap["n_freeboard_cells"] == 2
    assert snap["thesis_vorab_sources_single_shot"] is True
    assert snap["vorab_outer_source_policy"] == "engineering_single_shot"
    assert snap["outer_may_refresh_drying_pyro"] is False
    assert snap["nr_solid_freeboard"] == [[], []]
    assert "H2S" in snap["nr_gas_excluded"]


def test_q8_does_not_open_production_defaults() -> None:
    """Q8：本轮不切生产默认；hamel_outer_refresh 不因 Q9 入口而开放。"""
    from src.solvers.global_nr_solver import global_residual
    from src.solvers.nr_indexing import unpack_reactor

    phase2 = build_phase2_init_contract()
    assert phase2["thesis_vorab_sources_single_shot"] is True
    assert phase2["vorab_outer_source_policy"] == "engineering_single_shot"
    assert phase2["outer_may_refresh_drying_pyro"] is False
    reduced = build_reduced_3cell_contract()
    assert reduced["knobs"]["thesis_mode"] is False
    assert global_residual.__kwdefaults__["unpack_fn"] is unpack_reactor
