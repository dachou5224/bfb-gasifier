from __future__ import annotations

import math

import pytest

from src.core.reactor import (
    Reactor,
    ReactorConfig,
    _resolve_axial_explicit_heat_loss_watts,
    _resolve_axial_heat_loss_distribution,
)
from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    SOLID_ENTHALPY_MODE_HAMEL_EQ28,
)


def _compound_total(losses: list[float], freeboard_loss: float) -> float:
    remain = 1.0
    for frac in losses:
        remain *= 1.0 - float(frac)
    remain *= 1.0 - float(freeboard_loss)
    return 1.0 - remain


def test_bed_only_total_heat_loss_is_evenly_compounded_across_cells():
    cfg = ReactorConfig(n_cells=10, H_bed=6.0, H_freeboard=0.0, heat_loss_frac=0.1)
    bed_losses, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)

    assert len(bed_losses) == 10
    assert all(math.isclose(loss, bed_losses[0], rel_tol=0.0, abs_tol=1e-12) for loss in bed_losses)
    assert math.isclose(freeboard_loss, 0.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(_compound_total(bed_losses, freeboard_loss), 0.1, rel_tol=0.0, abs_tol=1e-12)


def test_total_heat_loss_is_split_across_bed_and_freeboard():
    cfg = ReactorConfig(n_cells=10, H_bed=6.0, H_freeboard=9.0, heat_loss_frac=0.1)
    bed_losses, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)

    assert len(bed_losses) == 10
    assert freeboard_loss > 0.0
    assert math.isclose(_compound_total(bed_losses, freeboard_loss), 0.1, rel_tol=0.0, abs_tol=1e-12)


def test_reactor_applies_distributed_bed_heat_loss_to_cells():
    cfg = ReactorConfig(n_cells=5, H_bed=6.0, H_freeboard=4.0, heat_loss_frac=0.1)
    expected_bed_losses, _ = _resolve_axial_heat_loss_distribution(cfg)

    reactor = Reactor(cfg)

    assert [cell.heat_loss_frac for cell in reactor.cells] == expected_bed_losses


def test_reactor_applies_freeboard_heat_loss_frac_to_explicit_cells():
    """NR freeboard cells must receive the freeboard inlet_fraction share (not 0)."""
    cfg = ReactorConfig(
        n_cells=4,
        n_freeboard_cells=4,
        H_bed=6.0,
        H_freeboard=8.0,
        heat_loss_frac=0.1,
        thesis_mode=True,
        explicit_freeboard_solver_graph_enabled=True,
    )
    _, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)
    reactor = Reactor(cfg)
    assert reactor.freeboard_cells
    fb_fracs = [float(c.heat_loss_frac) for c in reactor.freeboard_cells]
    assert all(f > 0.0 for f in fb_fracs)
    # Per-cell compound weights reconstruct the freeboard height share.
    remain = 1.0
    for frac in fb_fracs:
        remain *= 1.0 - frac
    assert math.isclose(1.0 - remain, freeboard_loss, rel_tol=0.0, abs_tol=1e-9)


def test_explicit_heat_loss_watts_split_by_axial_height():
    cfg = ReactorConfig(
        n_cells=4,
        H_bed=6.0,
        H_freeboard=3.0,
        Q_wall_W_total=900.0,
        Q_heat_exchanger_W_total=300.0,
    )
    bed_q_wall, bed_q_hex, fb_q_wall, fb_q_hex = _resolve_axial_explicit_heat_loss_watts(cfg)
    assert len(bed_q_wall) == 4
    assert sum(bed_q_wall) == pytest.approx(600.0)
    assert sum(bed_q_hex) == pytest.approx(200.0)
    assert fb_q_wall == pytest.approx(300.0)
    assert fb_q_hex == pytest.approx(100.0)


def test_reactor_applies_energy_balance_modes_to_cells():
    cfg = ReactorConfig(
        n_cells=3,
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        solid_enthalpy_mode=SOLID_ENTHALPY_MODE_HAMEL_EQ28,
        Q_wall_W_total=120.0,
        Q_heat_exchanger_W_total=30.0,
        H_bed=3.0,
        H_freeboard=0.0,
    )
    bed_q_wall, bed_q_hex, _, _ = _resolve_axial_explicit_heat_loss_watts(cfg)
    reactor = Reactor(cfg)
    for cell, q_wall, q_hex in zip(reactor.cells, bed_q_wall, bed_q_hex):
        assert cell.heat_loss_mode == HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT
        assert cell.solid_enthalpy_mode == SOLID_ENTHALPY_MODE_HAMEL_EQ28
        assert cell.Q_wall_W == pytest.approx(q_wall)
        assert cell.Q_heat_exchanger_W == pytest.approx(q_hex)
