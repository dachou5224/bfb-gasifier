"""阶段 D：冻结预计算上下文。未冻结非零库存可历史依赖；冻结后 F 纯。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.solvers.vorabrechnung.vorab_budgets import vorabrechnung_tau_for_cell
from src.solvers.vorabrechnung.vorab_policy import (
    VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT,
    VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH,
    VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF,
    resolve_vorab_outer_source_policy,
    vorab_outer_may_refresh_drying_pyro_sources,
    vorab_tau_model_label,
)
from src.solvers.wirsum_1997_vorab_context import (
    VORAB_CONTEXT_VERSION,
    build_stage_cd_coupled_payload,
    build_stage_d_mapping_payload,
    build_stage_d_purity_payload,
    build_stage_energy_outer_payload,
    eval_F_with_frozen_context,
    freeboard_trajectory_contract,
    snapshot_vorab_frozen_context,
)
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    MAPPING_EXPR,
    example_lu_like_feed,
    map_one_state,
    occupancy_plug_flow,
    map_release_by_occupancy,
)
from src.solvers.wirsum_1997_strict import apply_wirsum_1997_strict_config


def test_auto_policy_matches_existing_single_shot_flags() -> None:
    cfg = ReactorConfig()
    cfg.thesis_mode = True
    cfg.thesis_vorab_sources_single_shot = True
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT
    assert not vorab_outer_may_refresh_drying_pyro_sources(cfg)

    cfg.thesis_vorab_sources_single_shot = False
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH
    assert vorab_outer_may_refresh_drying_pyro_sources(cfg)


def test_explicit_hamel_outer_refresh_overrides_thesis_single_shot() -> None:
    cfg = ReactorConfig(
        thesis_mode=True,
        thesis_vorab_sources_single_shot=True,
        vorab_outer_source_policy="hamel_outer_refresh",
    )
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH
    assert vorab_outer_may_refresh_drying_pyro_sources(cfg)


def test_explicit_engineering_single_shot_overrides_nonthesis() -> None:
    cfg = ReactorConfig(
        thesis_mode=False,
        vorab_outer_source_policy="engineering_single_shot",
    )
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT
    assert not vorab_outer_may_refresh_drying_pyro_sources(cfg)


def test_unknown_vorab_outer_source_policy_rejected() -> None:
    with pytest.raises(ValueError, match="vorab_outer_source_policy"):
        ReactorConfig(vorab_outer_source_policy="clip")


def test_strict_config_does_not_open_hamel_outer_refresh() -> None:
    cfg = ReactorConfig(thesis_mode=True, thesis_vorab_sources_single_shot=True)
    apply_wirsum_1997_strict_config(cfg)
    assert cfg.vorab_outer_source_policy == "auto"
    assert cfg.thesis_vorab_sources_single_shot is True
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT


def test_tau_is_labeled_engineering_h_over_umf() -> None:
    src = inspect.getsource(vorabrechnung_tau_for_cell)
    assert "engineering" in src.lower()
    assert "Hamel 单次预算语义" not in src
    assert vorab_tau_model_label() == VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF


def test_freeboard_trajectory_not_owned_by_outer() -> None:
    fb = freeboard_trajectory_contract()
    assert fb["outer_updates_freeboard_trajectory"] is False
    assert fb["owner"] == "init_precalc"


def test_unfrozen_nonzero_holdup_F_is_history_dependent() -> None:
    payload = build_stage_d_purity_payload()
    assert payload["unfrozen_max_abs_F_xx"] > 1.0


def test_frozen_context_makes_F_repeatable_after_holdup_startwert() -> None:
    payload = build_stage_d_purity_payload()
    assert payload["tau_model"] == VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF
    assert payload["outer_updates_freeboard_trajectory"] is False
    assert payload["frozen_max_abs_F_xx"] == pytest.approx(0.0, abs=0.0)
    assert payload["frozen_max_abs_F_xyx"] == pytest.approx(0.0, abs=0.0)
    assert int(payload["n_unknowns"]) > 0
    assert VORAB_CONTEXT_VERSION == 1


def test_two_outer_maps_conserve_atoms_and_name_enthalpy() -> None:
    payload = build_stage_d_mapping_payload()
    assert payload["rewrites_unknown_N"] is False
    assert payload["default_policy_unchanged"] is True
    assert payload["mapping_expr"] == MAPPING_EXPR
    assert payload["is_hamel_closed_form"] is False
    assert payload["atom_C_gap_a"] < 1.0e-9
    assert payload["atom_C_gap_b"] < 1.0e-9
    assert payload["gas_hf_counted_once"] is True
    assert payload["vm_pool_basis"] == "fuel_waf_c"
    assert np.isfinite(payload["H_map_gap_a_W"])
    assert np.isfinite(payload["H_map_gap_b_W"])


def test_cd_coupled_payload_does_not_claim_furnace() -> None:
    payload = build_stage_cd_coupled_payload()
    assert payload["solved"] is False
    assert payload["n_unknowns"] == 39
    assert payload["physical_convergence"] is True
    assert payload["q"]["rewrites_unknown_N"] is False
    assert VORAB_CONTEXT_VERSION == 1


def test_energy_outer_payload_does_not_open_production_refresh() -> None:
    payload = build_stage_energy_outer_payload()
    assert payload["solved"] is False
    assert payload["outer_rewrites_N"] is False
    assert "production_hamel_outer_refresh" in payload["omitted"]
    assert VORAB_CONTEXT_VERSION == 1


def test_occupancy_integral_totals_match_conversion() -> None:
    feed = example_lu_like_feed()
    tau = np.array([0.5, 1.5], dtype=np.float64)
    occ = occupancy_plug_flow(tau, source="named_test_age_bins")
    X = np.array([0.0, 0.25, 1.0], dtype=np.float64)
    n_i = map_release_by_occupancy(feed.nC_vm, X, occ)
    assert float(np.sum(n_i)) == pytest.approx(feed.nC_vm)
    hid = occupancy_plug_flow(tau, source=VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF)
    assert hid["hides_H_over_umf"] is True
    qa = map_one_state(
        feed, P_Pa=2.5e6, T_K=1150.0, X_vm=X, tau_s=tau, occupancy_source="named_test_age_bins"
    )
    qb = map_one_state(
        feed, P_Pa=2.5e6, T_K=1150.0, X_vm=X, tau_s=tau, occupancy_source="named_test_age_bins"
    )
    assert qa["atom_C_gap"] == pytest.approx(qb["atom_C_gap"])
    assert qa["H_in_W"] == pytest.approx(qb["H_in_W"])


def test_frozen_context_version_mismatch_rejected() -> None:
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    ctx = snapshot_vorab_frozen_context(list(reactor.cells), reactor=reactor)
    ctx.version = VORAB_CONTEXT_VERSION + 1
    with pytest.raises(ValueError, match="version"):
        eval_F_with_frozen_context(
            np.zeros(1),
            list(reactor.cells),
            reactor._apply_all_bc_for_nr,
            ctx,
        )
