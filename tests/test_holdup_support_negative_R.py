"""holdup support 计入负 R_solid（§5.22）。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import Cell, N_SOLID_COMP, S_CHAR
from src.core.connectivity.solid_transport import (
    _transport_inflow_support_holdup,
    align_bed_solid_holdup_to_transport_inflow_support,
)


def _cell() -> Cell:
    c = Cell.__new__(Cell)
    c.solid_state_model = "holdup_transport"
    nk = 1
    c.m_solid = np.zeros((nk, N_SOLID_COMP))
    c.m_solid_zu = np.zeros((nk, N_SOLID_COMP))
    c.m_solid_rez = np.zeros((nk, N_SOLID_COMP))
    c.m_solid_in = np.zeros((nk, N_SOLID_COMP))
    c.m_solid_auf_in = np.zeros((nk, N_SOLID_COMP))
    c.m_solid_ab_in = np.zeros((nk, N_SOLID_COMP))
    c.K_solid_auf = np.zeros((nk, N_SOLID_COMP))
    c.K_solid_ab = np.zeros((nk, N_SOLID_COMP))
    c.R_solid = np.zeros((nk, N_SOLID_COMP))
    return c


def test_default_support_ignores_negative_R():
    c = _cell()
    c.m_solid_auf_in[0, S_CHAR] = 1.0
    c.K_solid_auf[0, S_CHAR] = 0.5
    c.R_solid[0, S_CHAR] = -0.2
    s_default = _transport_inflow_support_holdup(c)
    s_neg = _transport_inflow_support_holdup(c, include_negative_R=True)
    assert float(s_default[0, S_CHAR]) == pytest.approx(1.0 / 0.5)
    assert float(s_neg[0, S_CHAR]) == pytest.approx(0.8 / 0.5)
    assert float(s_neg[0, S_CHAR]) < float(s_default[0, S_CHAR])


def test_align_reduce_with_negative_R_targets_lower_holdup():
    c = _cell()
    c.m_solid_auf_in[0, S_CHAR] = 1.0
    c.K_solid_auf[0, S_CHAR] = 0.5
    c.K_solid_ab[0, S_CHAR] = 0.0
    c.R_solid[0, S_CHAR] = -0.2
    c.m_solid[0, S_CHAR] = 3.0
    changed = align_bed_solid_holdup_to_transport_inflow_support(
        c, allow_reduce=True, include_negative_R=True
    )
    assert changed is True
    assert float(c.m_solid[0, S_CHAR]) == pytest.approx(0.8 / 0.5)
