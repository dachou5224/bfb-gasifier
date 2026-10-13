from __future__ import annotations

from src.core.cell import Cell
from src.core.hydrodynamics_defaults import (
    HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL,
    HAMEL_HYDRO_BUBBLE_ODE_STRATEGY,
    HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY,
    HAMEL_HYDRO_U_D_CLOSURE,
)
from src.core.reactor import Reactor
from tests.validation.config_phase1 import build_phase1_htw_lu_reactor_config


def test_standalone_cell_defaults_to_hamel_hydrodynamics_chain() -> None:
    cell = Cell()
    assert cell.u_d_closure == HAMEL_HYDRO_U_D_CLOSURE
    assert cell.bubble_diameter_model == HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL
    assert cell.bubble_velocity_strategy == HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY
    assert cell.bubble_ode_strategy == HAMEL_HYDRO_BUBBLE_ODE_STRATEGY


def test_non_thesis_reactor_propagates_hamel_hydrodynamics_chain() -> None:
    cfg = build_phase1_htw_lu_reactor_config()
    assert cfg.thesis_mode is False
    reactor = Reactor(cfg)
    assert {cell.u_d_closure for cell in reactor.cells} == {HAMEL_HYDRO_U_D_CLOSURE}
    assert {cell.bubble_diameter_model for cell in reactor.cells} == {HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL}
    assert {cell.bubble_velocity_strategy for cell in reactor.cells} == {HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY}
    assert {cell.bubble_ode_strategy for cell in reactor.cells} == {HAMEL_HYDRO_BUBBLE_ODE_STRATEGY}
