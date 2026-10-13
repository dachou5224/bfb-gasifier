"""targeted phase-split post-gate 接线（§5.24）。"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import GAS_SPECIES, N_GAS
import numpy as np


def _bed() -> Cell:
    c = Cell.__new__(Cell)
    c.cell_type = "bed"
    c.solid_state_model = "holdup_transport"
    c.m_solid = np.ones((1, N_SOLID_COMP))
    c.N_d = np.zeros(N_GAS)
    c.N_b = np.zeros(N_GAS)
    c.T = 1100.0
    return c


def test_targeted_phase_split_calls_bed3_when_armed():
    from src.core.connectivity import preprojection as pp

    cells = [_bed() for _ in range(5)]
    j_co = GAS_SPECIES.index("CO")
    j_n2 = GAS_SPECIES.index("N2")
    cells[-1].N_d[j_co] = 0.2
    cells[-1].N_d[j_n2] = 0.8
    cfg = SimpleNamespace(
        nr_outer_preinner_phase_split_preprojection_thesis=False,
        nr_outer_preinner_solid_energy_preprojection_thesis=False,
        nr_outer_preinner_upper_dense_syngas_absorb_thesis=False,
        nr_outer_preinner_upper_co_local2x2_thesis=False,
        nr_outer_preinner_bed_top_solid_holdup_thesis=False,
        nr_outer_preinner_targeted_phase_split_beds_thesis=(3,),
        nr_outer_preinner_targeted_phase_split_after_y_co_thesis=True,
        nr_outer_preinner_targeted_phase_split_y_co_threshold_thesis=0.14,
        nr_outer_freeze_upper_co_bridge_after_y_co_thesis=False,
        nr_outer_freeze_upper_co_bridge_after_local2_accepts_thesis=0,
        fuel_feed=1.0,
        nr_outer_preinner_merit_driven_thesis=True,
        nr_outer_preinner_global_merit_gate_thesis=False,
        nr_outer_preinner_reconcile_holdup_before_zone_thesis=False,
    )

    def _fake_metrics(*_a, **_k):
        return {
            "merit_scaled": 0.1,
            "gas_phase_split_rms": 0.0,
            "gas_phase_split_max_abs": 0.0,
            "energy_rms": 0.0,
            "solid_rms": 0.0,
            "max_bed_energy_abs_W": 0.0,
        }

    with (
        patch.object(pp, "_measure_nr_init_scaled_metrics", side_effect=_fake_metrics),
        patch.object(pp, "max_upper_bed_solid_auf_flux_ratio_error", return_value=0.0),
        patch.object(pp, "max_upper_bed_solid_stream_imbalance", return_value=0.0),
        patch.object(pp, "resolve_preinner_tier_plan") as tier,
        patch.object(
            pp, "preproject_bed_major_gas_phase_split_to_exchange_closure", return_value=True
        ) as lsq,
        patch.object(pp, "_snapshot_nr_bottom_bed_states", return_value=[]),
    ):
        tier.return_value = {
            "merit_driven": True,
            "run_split_tier": False,
            "run_solid_energy_tier": False,
            "run_holdup_reconcile": False,
            "run_bed0_resplit": False,
        }
        diag = pp.run_outer_preinner_start_value_bridge(cells, cfg, apply_bc_fn=lambda: None)
    assert diag["targeted_phase_split_enabled"] is True
    assert diag["targeted_phase_split_armed"] is True
    assert diag["targeted_phase_split_accepted"] is True
    assert lsq.call_count == 1
    assert int(lsq.call_args.kwargs["cell_index"]) == 3
