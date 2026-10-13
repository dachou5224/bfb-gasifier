"""Eq. (2.7) 显式热损默认与 Q_W 解析。"""

from __future__ import annotations

import pytest

from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    HEAT_LOSS_MODE_INLET_FRACTION,
    resolve_reactor_Q_wall_W_total,
)
from src.core.reactor import ReactorConfig


def test_resolve_reactor_Q_wall_uses_explicit_total_when_set() -> None:
    cfg = ReactorConfig(
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        Q_wall_W_total=500_000.0,
    )
    assert resolve_reactor_Q_wall_W_total(cfg) == pytest.approx(500_000.0)


def test_resolve_reactor_Q_wall_bootstraps_from_fuel_hhv_when_unset() -> None:
    cfg = ReactorConfig(
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        heat_loss_frac=0.1,
        fuel_feed=1.0,
        HHV_dry=20.0,
        Q_wall_W_total=0.0,
    )
    assert resolve_reactor_Q_wall_W_total(cfg) == pytest.approx(2.0e6)


def test_resolve_reactor_Q_wall_zero_for_inlet_fraction_mode() -> None:
    cfg = ReactorConfig(
        energy_balance_heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION,
        heat_loss_frac=0.1,
        fuel_feed=1.0,
        HHV_dry=20.0,
    )
    assert resolve_reactor_Q_wall_W_total(cfg) == 0.0
