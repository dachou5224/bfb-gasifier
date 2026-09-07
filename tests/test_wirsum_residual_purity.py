"""P0.3：给定 packed x，global_residual 必须是纯函数。"""

from __future__ import annotations

import numpy as np

from src.core.cell import S_CHAR
from src.solvers.global_nr_solver import global_residual, pack_reactor
from tests.global_nr_fixtures import call_generate_initial_x0_phase1, phase1_axial_t_profile
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config
from src.core.reactor import Reactor


def _two_cell_reactor() -> Reactor:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    call_generate_initial_x0_phase1(reactor, cfg, phase1_axial_t_profile(cfg))
    return reactor


def _cell_physics_snapshot(cells) -> dict[str, np.ndarray]:
    return {
        "N_d": np.stack([np.array(c.N_d, dtype=np.float64, copy=True) for c in cells]),
        "N_b": np.stack([np.array(c.N_b, dtype=np.float64, copy=True) for c in cells]),
        "T": np.array([float(c.T) for c in cells], dtype=np.float64),
        "m_solid_char": np.stack(
            [np.array(c.m_solid[:, S_CHAR], dtype=np.float64, copy=True) for c in cells]
        ),
        "eps_b": np.array([float(c.eps_b) for c in cells], dtype=np.float64),
        "K_bd": np.array([float(c.K_bd) for c in cells], dtype=np.float64),
    }


def test_global_residual_does_not_mutate_packed_x() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x0 = pack_reactor(cells)
    x_in = x0.copy()
    _ = global_residual(x_in, cells, reactor._apply_all_bc_for_nr)
    np.testing.assert_array_equal(x_in, x0)


def test_global_residual_is_repeatable_on_same_x() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x0 = pack_reactor(cells)
    f1 = global_residual(x0.copy(), cells, reactor._apply_all_bc_for_nr)
    snap1 = _cell_physics_snapshot(cells)
    f2 = global_residual(x0.copy(), cells, reactor._apply_all_bc_for_nr)
    snap2 = _cell_physics_snapshot(cells)
    np.testing.assert_allclose(f1, f2, rtol=0.0, atol=0.0)
    for key in snap1:
        np.testing.assert_allclose(snap1[key], snap2[key], rtol=0.0, atol=0.0)
