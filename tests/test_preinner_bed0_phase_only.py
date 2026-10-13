"""bed0 Pre-Inner split：phase-only majors（禁止产物凭空造量）。"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import src.core.connectivity.preprojection as pp


def test_bed0_phase_only_passes_preserve_primary_false():
    seen: dict[str, bool] = {}

    class _Cell:
        cell_type = "bed"

        def __init__(self) -> None:
            self.N_d = np.ones(3, dtype=np.float64)
            self.N_b = np.ones(3, dtype=np.float64)
            self._thermo_cache_valid = True
            self._hydro_cache_valid = True

    cfg = SimpleNamespace(
        nr_outer_preinner_phase_split_preprojection_thesis=True,
        nr_outer_preinner_phase_split_max_bed_index_thesis=0,
        nr_outer_preinner_bed1_phase_split_thesis=False,
        nr_outer_preinner_bed0_phase_only_majors_thesis=True,
        nr_outer_preinner_split_global_residual_guard_thesis=True,
        nr_outer_preinner_split_global_rms_factor_thesis=1.25,
        nr_outer_preinner_split_global_max_abs_factor_thesis=1.25,
        nr_outer_preinner_split_global_max_abs_delta_thesis=0.5,
        nr_outer_preinner_split_plateau_bed0_improvement_thesis=0.85,
        nr_outer_preinner_bed1_hotspot_followup_thesis=False,
        O2_feed=1.0,
        H2O_feed=1.0,
        N2_feed=1.0,
        fuel_feed=1.0,
    )

    def _fake_lsq(*_a, **kwargs):
        seen["preserve"] = bool(kwargs.get("preserve_primary_inlet_totals"))
        return False

    metrics = {
        "gas_phase_split_rms": 0.05,
        "rms_scaled": 0.06,
        "max_abs_scaled": 0.9,
        "merit_scaled": 0.06,
        "solid_rms": 0.01,
        "energy_rms": 0.04,
    }

    with patch.object(pp, "_measure_nr_init_scaled_metrics", return_value=dict(metrics)), patch.object(
        pp, "_select_outer_preinner_phase_split_bed_indices", return_value=[0]
    ), patch.object(
        pp, "preproject_bed_major_gas_phase_split_to_exchange_closure", side_effect=_fake_lsq
    ), patch.object(pp, "_per_bed_nr_phase_split_rms_map", return_value={0: 0.05}):
        pp.run_outer_preinner_gas_phase_split_preprojection(
            [_Cell()], cfg, apply_bc_fn=lambda: None
        )

    assert seen["preserve"] is False


def test_bed0_legacy_preserve_primary_when_phase_only_off():
    seen: dict[str, bool] = {}

    class _Cell:
        cell_type = "bed"

        def __init__(self) -> None:
            self.N_d = np.ones(3, dtype=np.float64)
            self.N_b = np.ones(3, dtype=np.float64)
            self._thermo_cache_valid = True
            self._hydro_cache_valid = True

    cfg = SimpleNamespace(
        nr_outer_preinner_phase_split_preprojection_thesis=True,
        nr_outer_preinner_phase_split_max_bed_index_thesis=0,
        nr_outer_preinner_bed1_phase_split_thesis=False,
        nr_outer_preinner_bed0_phase_only_majors_thesis=False,
        nr_outer_preinner_split_global_residual_guard_thesis=True,
        nr_outer_preinner_split_global_rms_factor_thesis=1.25,
        nr_outer_preinner_split_global_max_abs_factor_thesis=1.25,
        nr_outer_preinner_split_global_max_abs_delta_thesis=0.5,
        nr_outer_preinner_split_plateau_bed0_improvement_thesis=0.85,
        nr_outer_preinner_bed1_hotspot_followup_thesis=False,
        O2_feed=1.0,
        H2O_feed=1.0,
        N2_feed=1.0,
        fuel_feed=1.0,
    )

    def _fake_lsq(*_a, **kwargs):
        seen["preserve"] = bool(kwargs.get("preserve_primary_inlet_totals"))
        return False

    metrics = {
        "gas_phase_split_rms": 0.05,
        "rms_scaled": 0.06,
        "max_abs_scaled": 0.9,
        "merit_scaled": 0.06,
        "solid_rms": 0.01,
        "energy_rms": 0.04,
    }

    with patch.object(pp, "_measure_nr_init_scaled_metrics", return_value=dict(metrics)), patch.object(
        pp, "_select_outer_preinner_phase_split_bed_indices", return_value=[0]
    ), patch.object(
        pp, "preproject_bed_major_gas_phase_split_to_exchange_closure", side_effect=_fake_lsq
    ), patch.object(pp, "_per_bed_nr_phase_split_rms_map", return_value={0: 0.05}):
        pp.run_outer_preinner_gas_phase_split_preprojection(
            [_Cell()], cfg, apply_bc_fn=lambda: None
        )

    assert seen["preserve"] is True
