"""Merit-driven outer pre-inner tier scheduling."""

from __future__ import annotations

import numpy as np

from src.core.connectivity import resolve_preinner_tier_plan
from src.core.reactor import ReactorConfig


def _cfg(**overrides: object) -> ReactorConfig:
    cfg = ReactorConfig(thesis_mode=True)
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def test_merit_driven_skips_tiers_when_residuals_are_small():
    cfg = _cfg(
        nr_outer_preinner_merit_driven_thesis=True,
        nr_outer_preinner_split_tier_min_rms_thesis=0.05,
        nr_outer_preinner_energy_tier_min_rms_thesis=0.08,
        nr_outer_preinner_solid_tier_min_rms_thesis=0.02,
    )
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.01,
            "merit_scaled": 0.012,
            "gas_phase_split_rms": 0.02,
            "energy_rms": 0.01,
            "solid_rms": 0.005,
        },
        cfg,
    )
    assert plan["merit_driven"] is True
    assert plan["run_split_tier"] is False
    assert plan["run_solid_energy_tier"] is False
    assert plan["run_holdup_reconcile"] is False


def test_merit_driven_runs_split_when_split_rms_elevated():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.04,
            "merit_scaled": 0.06,
            "gas_phase_split_rms": 0.06,
            "energy_rms": 0.02,
            "solid_rms": 0.005,
        },
        cfg,
    )
    assert plan["run_split_tier"] is True
    assert plan["run_solid_energy_tier"] is False


def test_merit_driven_runs_split_when_max_abs_hot_despite_low_rms():
    """Global S can sit just below tau_split while bed0 CO₂ max|split| stays large."""
    cfg = _cfg(
        nr_outer_preinner_merit_driven_thesis=True,
        nr_outer_preinner_split_tier_min_rms_thesis=0.03,
        nr_outer_preinner_split_tier_min_max_abs_thesis=0.05,
    )
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.02,
            "merit_scaled": 0.028,
            "gas_phase_split_rms": 0.028,
            "gas_phase_split_max_abs": 0.098,
            "energy_rms": 0.015,
            "solid_rms": 0.005,
        },
        cfg,
    )
    assert plan["run_split_tier"] is True


def test_merit_driven_runs_energy_tiers_when_energy_elevated():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.04,
            "merit_scaled": 0.10,
            "gas_phase_split_rms": 0.02,
            "energy_rms": 0.10,
            "solid_rms": 0.005,
        },
        cfg,
    )
    assert plan["run_holdup_reconcile"] is True
    assert plan["run_solid_energy_tier"] is True
    assert plan["run_split_tier"] is False


def test_merit_driven_runs_holdup_when_flux_mismatch_elevated():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.02,
            "merit_scaled": 0.025,
            "gas_phase_split_rms": 0.02,
            "energy_rms": 0.005,
            "solid_rms": 0.005,
            "max_bed_energy_abs_W": 1.0e5,
        },
        cfg,
        holdup_flux_mismatch=0.25,
        holdup_stream_imbalance=0.25,
    )
    assert plan["run_holdup_reconcile"] is True
    assert plan["run_solid_energy_tier"] is False


def test_merit_driven_skips_holdup_when_stream_balanced_despite_auf_gap():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.02,
            "merit_scaled": 0.025,
            "gas_phase_split_rms": 0.02,
            "energy_rms": 0.005,
            "solid_rms": 0.005,
            "max_bed_energy_abs_W": 1.0e5,
        },
        cfg,
        holdup_flux_mismatch=0.72,
        holdup_stream_imbalance=0.02,
    )
    assert plan["run_holdup_reconcile"] is False


def test_merit_driven_skips_zone_on_split_shock():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.20,
            "merit_scaled": 0.22,
            "gas_phase_split_rms": 0.22,
            "energy_rms": 0.10,
            "solid_rms": 0.005,
            "max_bed_energy_abs_W": 5.0e6,
        },
        cfg,
    )
    assert plan["run_split_tier"] is True
    assert plan["split_shock"] is True
    assert plan["run_solid_energy_tier"] is False
    assert plan["run_holdup_reconcile"] is False


def test_merit_driven_runs_holdup_on_split_shock_only_for_transport_mismatch():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.22,
            "merit_scaled": 0.22,
            "gas_phase_split_rms": 0.22,
            "energy_rms": 0.10,
            "solid_rms": 0.005,
            "max_bed_energy_abs_W": 5.0e6,
        },
        cfg,
        holdup_flux_mismatch=0.30,
        holdup_stream_imbalance=0.20,
    )
    assert plan["run_holdup_reconcile"] is True


def test_merit_driven_skips_holdup_on_split_shock_when_only_global_energy_elevated():
    cfg = _cfg(nr_outer_preinner_merit_driven_thesis=True)
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.22,
            "merit_scaled": 0.22,
            "gas_phase_split_rms": 0.22,
            "energy_rms": 0.10,
            "solid_rms": 0.005,
            "max_bed_energy_abs_W": 1.0e6,
        },
        cfg,
        holdup_flux_mismatch=0.72,
        holdup_stream_imbalance=0.08,
    )
    assert plan["split_shock"] is True
    assert plan["run_holdup_reconcile"] is False
    assert plan["run_split_tier"] is True


def test_outer_preinner_phase_split_bed_indices_skips_bed1_by_default():
    from src.core.connectivity.preprojection import _outer_preinner_phase_split_bed_indices

    class _Cell:
        cell_type = "bed"

    cells = [_Cell(), _Cell(), _Cell()]
    assert _outer_preinner_phase_split_bed_indices(
        cells, max_bed_index=2, include_bed1=False
    ) == [0, 2]
    assert _outer_preinner_phase_split_bed_indices(
        cells, max_bed_index=2, include_bed1=True
    ) == [0, 1, 2]
    # Flag alone must enable bed1 when max_bed_index=0 (Phase2 / VTT contract).
    assert _outer_preinner_phase_split_bed_indices(
        cells, max_bed_index=0, include_bed1=True
    ) == [0, 1]
    assert _outer_preinner_phase_split_bed_indices(
        cells, max_bed_index=0, include_bed1=False
    ) == [0]


def test_bed1_hotspot_followup_should_run_after_bed0_cleared():
    from src.core.connectivity.preprojection import _bed1_hotspot_followup_should_run

    cfg = _cfg(
        nr_outer_preinner_bed1_hotspot_followup_thesis=True,
        nr_outer_preinner_bed1_hotspot_min_rms_thesis=0.05,
        nr_outer_preinner_bed1_hotspot_max_bed0_rms_thesis=0.01,
        nr_outer_preinner_bed1_hotspot_max_global_rms_thesis=0.04,
    )
    # Late floor: bed0 cleared, bed1 H₂ still hot, global S already low.
    assert _bed1_hotspot_followup_should_run(
        {0: 1e-9, 1: 0.059},
        cfg,
        already_selected=False,
        global_split_rms=0.021,
    )
    # Mid-solve: global S still high — do not follow up.
    assert not _bed1_hotspot_followup_should_run(
        {0: 1e-9, 1: 0.059},
        cfg,
        already_selected=False,
        global_split_rms=0.10,
    )
    # bed0 still hot — do not follow up.
    assert not _bed1_hotspot_followup_should_run(
        {0: 0.061, 1: 0.058},
        cfg,
        already_selected=False,
        global_split_rms=0.021,
    )
    assert not _bed1_hotspot_followup_should_run(
        {0: 1e-9, 1: 0.059},
        cfg,
        already_selected=True,
        global_split_rms=0.021,
    )


def test_select_outer_preinner_phase_split_bed_indices_prefers_hot_upper_beds():
    from src.core.connectivity.preprojection import _select_outer_preinner_phase_split_bed_indices

    class _Cfg:
        nr_outer_preinner_phase_split_max_bed_index_thesis = 2
        nr_outer_preinner_bed1_phase_split_thesis = False
        nr_outer_preinner_phase_split_upper_bed_min_rms_thesis = 0.06
        O2_feed = 1.0
        H2O_feed = 1.0
        N2_feed = 1.0
        fuel_feed = 1.0

    class _Cell:
        cell_type = "bed"

        def __init__(self) -> None:
            self.N_d = np.zeros(7)
            self.N_b = np.zeros(7)

    cells = [_Cell(), _Cell(), _Cell()]
    split_map = {0: 0.04, 1: 0.02, 2: 0.08}

    def _fake_map(_cells, _cfg, *, apply_bc_fn):
        apply_bc_fn()
        return split_map

    import src.core.connectivity.preprojection as pp

    original = pp._per_bed_nr_phase_split_rms_map
    pp._per_bed_nr_phase_split_rms_map = _fake_map
    try:
        assert _select_outer_preinner_phase_split_bed_indices(
            cells, _Cfg(), apply_bc_fn=lambda: None
        ) == [0, 2]
        # 平台期（max local≥0.03）τ_eff=0.015：bed2=0.05 仍入选
        split_map[2] = 0.05
        assert _select_outer_preinner_phase_split_bed_indices(
            cells, _Cfg(), apply_bc_fn=lambda: None
        ) == [0, 2]
        # Near-tie with bed0
        split_map[0] = 0.064
        split_map[2] = 0.067
        assert _select_outer_preinner_phase_split_bed_indices(
            cells, _Cfg(), apply_bc_fn=lambda: None
        ) == [0, 2]
        # Explicit bed1 flag on plateau
        cfg_bed1 = _Cfg()
        cfg_bed1.nr_outer_preinner_bed1_phase_split_thesis = True
        cfg_bed1.nr_outer_preinner_phase_split_max_bed_index_thesis = 0
        split_map[0] = 0.04
        split_map[1] = 0.025
        assert _select_outer_preinner_phase_split_bed_indices(
            cells, cfg_bed1, apply_bc_fn=lambda: None
        ) == [0, 1]
        split_map[0] = 0.01
        split_map[1] = 0.01
        split_map[2] = 0.01
        # 全局冷：非平台，bed1=0.01 < τ=0.06
        assert _select_outer_preinner_phase_split_bed_indices(
            cells, cfg_bed1, apply_bc_fn=lambda: None
        ) == [0]
    finally:
        pp._per_bed_nr_phase_split_rms_map = original


def test_legacy_mode_keeps_config_flags():
    cfg = _cfg(
        nr_outer_preinner_merit_driven_thesis=False,
        nr_outer_preinner_phase_split_preprojection_thesis=True,
        nr_outer_preinner_solid_energy_preprojection_thesis=False,
    )
    plan = resolve_preinner_tier_plan(
        {
            "rms_scaled": 0.5,
            "merit_scaled": 0.5,
            "gas_phase_split_rms": 0.5,
            "energy_rms": 0.5,
            "solid_rms": 0.5,
        },
        cfg,
    )
    assert plan["merit_driven"] is False
    assert plan["run_split_tier"] is True
    assert plan["run_solid_energy_tier"] is False


def test_check1_continuation_skips_holdup_reconcile_by_default():
    """Bild 2.2: Check1 continuation must not algebraically remap holdup (F4)."""
    from contextlib import ExitStack
    from unittest.mock import patch

    import src.core.connectivity.preprojection as pp

    cfg = _cfg(
        nr_outer_preinner_phase_split_preprojection_thesis=True,
        nr_outer_preinner_solid_energy_preprojection_thesis=False,
        nr_outer_preinner_merit_driven_thesis=True,
        nr_outer_preinner_holdup_on_check1_continuation_thesis=False,
    )

    class _Cell:
        cell_type = "bed"

        def __init__(self) -> None:
            self.N_d = np.ones(7)
            self.N_b = np.ones(7)
            self.T = 1000.0
            self.m_solid = np.ones((4, 5))

    cells = [_Cell(), _Cell(), _Cell()]
    metrics = {
        "rms_scaled": 0.12,
        "merit_scaled": 0.12,
        "gas_phase_split_rms": 0.12,
        "energy_rms": 0.05,
        "solid_rms": 0.002,
        "max_abs_scaled": 0.3,
        "max_bed_energy_abs_W": 0.0,
    }
    holdup_calls = {"n": 0}

    def _fake_holdup(*_a, **_k):
        holdup_calls["n"] += 1
        return True

    with ExitStack() as stack:
        stack.enter_context(
            patch.object(pp, "_measure_nr_init_scaled_metrics", return_value=dict(metrics))
        )
        stack.enter_context(
            patch.object(pp, "max_upper_bed_solid_auf_flux_ratio_error", return_value=0.7)
        )
        stack.enter_context(
            patch.object(pp, "max_upper_bed_solid_stream_imbalance", return_value=0.2)
        )
        stack.enter_context(
            patch.object(pp, "reconcile_upper_bed_solid_holdup_for_init", side_effect=_fake_holdup)
        )
        stack.enter_context(
            patch.object(
                pp,
                "run_outer_preinner_gas_phase_split_preprojection",
                return_value={"enabled": True, "bed0_accepted": False},
            )
        )
        diag = pp.run_outer_preinner_start_value_bridge(
            cells,
            cfg,
            apply_bc_fn=lambda: None,
            check1_continuation=True,
        )
    assert diag["holdup_skipped_check1_continuation"] is True
    assert diag["tier_plan"]["run_holdup_reconcile"] is False
    assert holdup_calls["n"] == 0


def test_preinner_gate_rollback_restores_upper_bed_holdup():
    """F1: holdup edits on bed1+ must roll back with the global merit gate."""
    from contextlib import ExitStack
    from unittest.mock import patch

    import src.core.connectivity.preprojection as pp

    cfg = _cfg(
        nr_outer_preinner_phase_split_preprojection_thesis=True,
        nr_outer_preinner_solid_energy_preprojection_thesis=False,
        nr_outer_preinner_merit_driven_thesis=True,
        nr_outer_preinner_global_merit_gate_thesis=True,
        nr_line_search_gas_phase_split_merit_thesis=True,
        nr_bottom_joint_energy_aware_gate_thesis=True,
    )

    class _Cell:
        cell_type = "bed"

        def __init__(self, m0: float) -> None:
            self.N_d = np.ones(7)
            self.N_b = np.ones(7)
            self.T = 1000.0
            self.m_solid = np.full((4, 5), m0)

    cells = [_Cell(1.0), _Cell(2.0), _Cell(3.0)]
    m_before = {
        "rms_scaled": 0.10,
        "merit_scaled": 0.10,
        "gas_phase_split_rms": 0.10,
        "energy_rms": 0.04,
        "solid_rms": 0.001,
        "max_abs_scaled": 0.2,
        "max_bed_energy_abs_W": 0.0,
    }
    m_worse = {
        "rms_scaled": 0.20,
        "merit_scaled": 0.20,
        "gas_phase_split_rms": 0.20,
        "energy_rms": 0.18,
        "solid_rms": 0.01,
        "max_abs_scaled": 0.4,
        "max_bed_energy_abs_W": 0.0,
    }
    measure_calls = {"n": 0}

    def _measure(*_a, **_k):
        measure_calls["n"] += 1
        if measure_calls["n"] == 1:
            return dict(m_before)
        if measure_calls["n"] == 2:
            return dict(m_worse)
        return dict(m_before)

    def _fake_holdup(bed_cells, *_a, **_k):
        bed_cells[1].m_solid[:] = 99.0
        bed_cells[2].m_solid[:] = 88.0
        return True

    with ExitStack() as stack:
        stack.enter_context(
            patch.object(pp, "_measure_nr_init_scaled_metrics", side_effect=_measure)
        )
        stack.enter_context(
            patch.object(pp, "max_upper_bed_solid_auf_flux_ratio_error", return_value=0.7)
        )
        stack.enter_context(
            patch.object(pp, "max_upper_bed_solid_stream_imbalance", return_value=0.2)
        )
        stack.enter_context(
            patch.object(pp, "reconcile_upper_bed_solid_holdup_for_init", side_effect=_fake_holdup)
        )
        stack.enter_context(
            patch.object(
                pp,
                "run_outer_preinner_gas_phase_split_preprojection",
                return_value={"enabled": True, "bed0_accepted": True},
            )
        )
        stack.enter_context(patch.object(pp, "_bottom_joint_global_gate_pass", return_value=False))
        stack.enter_context(
            patch.object(pp, "_preinner_split_only_fallback_gate_pass", return_value=False)
        )
        diag = pp.run_outer_preinner_start_value_bridge(
            cells,
            cfg,
            apply_bc_fn=lambda: None,
            check1_continuation=False,
        )

    assert diag["global_gate_pass"] is False
    assert diag["bed0_accepted"] is False
    assert diag["holdup_reconciled"] is False
    assert float(cells[1].m_solid[0, 0]) == 2.0
    assert float(cells[2].m_solid[0, 0]) == 3.0
    assert measure_calls["n"] >= 3
