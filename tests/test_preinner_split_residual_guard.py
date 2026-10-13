"""Pre-Inner split LSQ：全局残差守卫回滚（handoff §5.9）。"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import src.core.connectivity.preprojection as pp


def test_split_lsq_rolls_back_when_global_max_abs_detonates():
    """局部 LSQ 返回 True 且 split 略降，但 max_abs 爆炸 → 必须回滚 N_d/N_b。"""

    class _Cell:
        cell_type = "bed"

        def __init__(self) -> None:
            self.N_d = np.array([1.0, 2.0, 3.0], dtype=np.float64)
            self.N_b = np.array([1.0, 2.0, 3.0], dtype=np.float64)
            self._thermo_cache_valid = True
            self._hydro_cache_valid = True

    cell0 = _Cell()
    cells = [cell0]
    snap_d = cell0.N_d.copy()
    snap_b = cell0.N_b.copy()

    cfg = SimpleNamespace(
        nr_outer_preinner_phase_split_preprojection_thesis=True,
        nr_outer_preinner_phase_split_max_bed_index_thesis=0,
        nr_outer_preinner_bed1_phase_split_thesis=False,
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

    metrics_seq = [
        # before
        {
            "gas_phase_split_rms": 0.05,
            "rms_scaled": 0.06,
            "max_abs_scaled": 0.9,
            "merit_scaled": 0.06,
            "solid_rms": 0.01,
            "energy_rms": 0.04,
        },
        # after bed0 LSQ (detonation)
        {
            "gas_phase_split_rms": 0.04,
            "rms_scaled": 13.5,
            "max_abs_scaled": 148.0,
            "merit_scaled": 13.5,
            "solid_rms": 0.01,
            "energy_rms": 0.05,
        },
        # after rollback remeasure (final)
        {
            "gas_phase_split_rms": 0.05,
            "rms_scaled": 0.06,
            "max_abs_scaled": 0.9,
            "merit_scaled": 0.06,
            "solid_rms": 0.01,
            "energy_rms": 0.04,
        },
    ]
    call_i = {"n": 0}

    def _fake_metrics(*_a, **_k):
        i = min(call_i["n"], len(metrics_seq) - 1)
        call_i["n"] += 1
        return dict(metrics_seq[i])

    def _fake_lsq(*_a, **_k):
        cell0.N_d[:] = 99.0
        cell0.N_b[:] = 0.01
        return True

    with patch.object(pp, "_measure_nr_init_scaled_metrics", side_effect=_fake_metrics), patch.object(
        pp, "_select_outer_preinner_phase_split_bed_indices", return_value=[0]
    ), patch.object(
        pp, "preproject_bed_major_gas_phase_split_to_exchange_closure", side_effect=_fake_lsq
    ), patch.object(pp, "_per_bed_nr_phase_split_rms_map", return_value={0: 0.05}):
        diag = pp.run_outer_preinner_gas_phase_split_preprojection(
            cells, cfg, apply_bc_fn=lambda: None
        )

    assert diag["enabled"] is True
    assert diag["bed0_accepted"] is False
    assert float(diag["global_residual_rollbacks"]) >= 1.0
    assert np.allclose(cell0.N_d, snap_d)
    assert np.allclose(cell0.N_b, snap_b)
