from __future__ import annotations

import math

from src.core.reactor import Reactor, ReactorConfig, _resolve_axial_heat_loss_distribution


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
