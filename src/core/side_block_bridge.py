"""侧支路（cyclone/return-leg）状态桥接服务。

将主链（bed/freeboard）顶端状态同步到 side-block cells，供 NR 边界路径使用。
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from src.core.connectivity import route_auxiliary_side_blocks as _route_auxiliary_side_blocks_exec

if TYPE_CHECKING:
    from src.core.reactor import Reactor


def initialize_explicit_side_block_states(reactor: "Reactor") -> None:
    """Initialize side-block cells from current top source chain."""
    if not reactor._use_explicit_side_block_cells():
        return
    assert reactor.cyclone_cell is not None and reactor.return_leg_cell is not None
    source_cells = (
        reactor.freeboard_cells
        if reactor.freeboard_cells and any(float(np.sum(c.N_d + c.N_b)) > 0.0 for c in reactor.freeboard_cells)
        else reactor.cells
    )
    top = source_cells[-1]
    _route_auxiliary_side_blocks_exec(source_cells, reactor.cyclone_cell, reactor.return_leg_cell, reactor.config)
    reactor.cyclone_cell.N_d[:] = np.maximum(top.N_d, 0.0)
    reactor.cyclone_cell.N_b[:] = np.maximum(top.N_b, 0.0)
    reactor.cyclone_cell.m_solid.fill(0.0)
    reactor.cyclone_cell.T = float(top.T)
    reactor.return_leg_cell.N_d.fill(0.0)
    reactor.return_leg_cell.N_b.fill(0.0)
    reactor.return_leg_cell.m_solid.fill(0.0)
    reactor.return_leg_cell.T = float(reactor.cyclone_cell.T)

