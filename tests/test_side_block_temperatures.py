from __future__ import annotations

from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config
from src.core.reactor import Reactor
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr


def test_side_block_boundary_refresh_inherits_cyclone_and_return_leg_temperatures() -> None:
    cfg = build_phase2_htw_lu_freeboard_reactor_config(n_freeboard_cells=2)
    cfg.n_cells = 3
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )

    assert reactor.cyclone_cell is not None
    assert reactor.return_leg_cell is not None

    top = reactor.freeboard_cells[-1]
    top.T = 912.0
    reactor.cyclone_cell.T = 777.0
    reactor.return_leg_cell.T = 666.0

    reactor._apply_all_bc_for_nr()

    assert reactor.cyclone_cell.T == top.T
    assert reactor.return_leg_cell.T == reactor.cyclone_cell.T
    assert reactor.return_leg_cell.T_solid_ab_in == reactor.cyclone_cell.T
