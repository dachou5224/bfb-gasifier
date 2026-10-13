"""Vorab single-shot 隔离策略单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import Cell
from src.core.reactor import ReactorConfig
from src.solvers.vorab_cell_abgleich import run_outer_vorab_cell_abgleich_bridge
from src.solvers.vorabrechnung.vorab_budgets import VorabrechnungCellBudget
from src.solvers.vorabrechnung.vorab_policy import (
    VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT,
    VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH,
    clear_vorab_drying_pyro_source_freeze_for_cells,
    freeze_vorab_drying_pyro_sources_for_cells,
    resolve_vorab_outer_source_policy,
    thesis_vorab_sources_single_shot_enabled,
    vorab_outer_may_refresh_drying_pyro_sources,
)
from src.solvers.vorabrechnung.vorab_refresh import refresh_cell_vorabrechnung


def test_thesis_single_shot_policy_flags() -> None:
    cfg = ReactorConfig()
    cfg.thesis_mode = True
    cfg.thesis_vorab_sources_single_shot = True
    assert thesis_vorab_sources_single_shot_enabled(cfg)
    assert not vorab_outer_may_refresh_drying_pyro_sources(cfg)

    cfg.thesis_vorab_sources_single_shot = False
    assert not thesis_vorab_sources_single_shot_enabled(cfg)
    assert vorab_outer_may_refresh_drying_pyro_sources(cfg)
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH


def test_configured_policy_overrides_auto_single_shot() -> None:
    cfg = ReactorConfig()
    cfg.thesis_mode = True
    cfg.thesis_vorab_sources_single_shot = True
    cfg.vorab_outer_source_policy = "hamel_outer_refresh"
    assert thesis_vorab_sources_single_shot_enabled(cfg)
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH
    assert vorab_outer_may_refresh_drying_pyro_sources(cfg)

    cfg.vorab_outer_source_policy = "engineering_single_shot"
    cfg.thesis_mode = False
    assert resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT
    assert not vorab_outer_may_refresh_drying_pyro_sources(cfg)


def test_frozen_vorab_sources_skip_invalidate_and_full_refresh(monkeypatch) -> None:
    cell = Cell()
    cell._vm_cache_valid = True
    cell._vm_gas_source_cache = np.ones(11)
    cell._vm_solid_sink_cache = np.zeros((1, 4))
    freeze_vorab_drying_pyro_sources_for_cells([cell])

    compute_calls = {"n": 0}

    def _compute(tau: float) -> None:
        compute_calls["n"] += 1

    monkeypatch.setattr(cell, "compute_vorabrechnung", _compute)
    monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: None)
    monkeypatch.setattr(cell, "_snapshot_vorabrechnung_hydrodynamics", lambda: None)

    cell.invalidate_vorabrechnung_cache()
    assert cell._vm_cache_valid is True
    assert compute_calls["n"] == 0

    refresh_cell_vorabrechnung(cell, force=True)
    assert compute_calls["n"] == 0

    clear_vorab_drying_pyro_source_freeze_for_cells([cell])
    cell.invalidate_vorabrechnung_cache()
    assert cell._vm_cache_valid is False


def test_outer_abgleich_bridge_skips_vm_zone_source_refresh_when_single_shot(monkeypatch) -> None:
    from src.core.reactor import Reactor

    cfg = ReactorConfig()
    cfg.thesis_mode = True
    cfg.thesis_vorab_sources_single_shot = True
    cfg.vorab_a_tier_outer_abgleich_bridge_thesis = True
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 2

    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    cell._vm_cache_valid = True
    freeze_vorab_drying_pyro_sources_for_cells([cell])

    refresh_calls = {"n": 0}

    monkeypatch.setattr(
        "src.solvers.vorabrechnung.vorab_budgets.refresh_vm_zone_fixed_pyrolysis_sources",
        lambda *args, **kwargs: refresh_calls.__setitem__("n", refresh_calls["n"] + 1)
        or {"vm_zone_sources_refreshed": 1.0},
    )
    monkeypatch.setattr(
        "src.solvers.vorab_cell_abgleich.assess_vorabrechnung_cell_alignment",
        lambda *args, **kwargs: {"aligned": False},
    )
    monkeypatch.setattr(
        "src.solvers.vorab_cell_abgleich.apply_vorabrechnung_bed_stream_mapping",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        "src.core.connectivity.preproject_bottom_primary_gas_state_to_exchange_closure",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        "src.core.connectivity.preproject_bottom_major_gas_state_to_local_balance",
        lambda *args, **kwargs: None,
    )

    budget = VorabrechnungCellBudget(
        cell_index=0,
        axial_xi=0.0,
        tau_cumulative_s=1.0,
        t_budget_K=1000.0,
        o2_remaining_mol_s=1.0,
        x_dry_cumulative=0.1,
        x_vm_cumulative=0.2,
        nC_vm_mol_s=0.01,
        nH_vm_mol_s=0.02,
        nO_vm_mol_s=0.005,
    )
    run_outer_vorab_cell_abgleich_bridge(
        reactor.cells,
        [budget],
        cfg,
        apply_bc_fn=lambda: None,
    )
    assert refresh_calls["n"] == 0
