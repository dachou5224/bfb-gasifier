"""holdup_transport 下炭转化 X 须计入 auf/ab 轴向入口。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import Cell, CellGeometry, S_CHAR, SolidProps


def _minimal_holdup_cell(*, n_size: int = 2) -> Cell:
    geo = CellGeometry(D_bed=1.0, dh=0.1, h_center=0.05)
    d_p_classes = np.linspace(1.5e-3, 3.0e-3, n_size)
    mass_fractions = np.zeros(n_size)
    mass_fractions[-1] = 1.0
    solid = SolidProps(
        rho_s=1200.0,
        d_p=2.0e-3,
        phi_s=0.8,
        eps_mf=0.45,
        n_size_classes=n_size,
        d_p_classes=d_p_classes,
        mass_fractions=mass_fractions,
    )
    cell = Cell(geo=geo, solid=solid, fuel_type="coal")
    cell.solid_state_model = "holdup_transport"
    cell.cell_type = "bed"
    return cell


def test_holdup_char_conversion_includes_auf_ab_inflows() -> None:
    """上层床炭入口在 auf/ab；忽略时 X 会被错误 clip 为 0。"""
    cell = _minimal_holdup_cell()
    cell.m_solid_zu.fill(0.0)
    cell.m_solid_in.fill(0.0)
    cell.m_solid_rez.fill(0.0)
    cell.m_solid_auf_in.fill(0.0)
    cell.m_solid_ab_in.fill(0.0)
    cell.m_solid.fill(0.0)
    cell.K_solid_auf.fill(0.0)
    cell.K_solid_ab.fill(0.0)

    # 轴向入口 1.0 kg/s，出流 0.9 kg/s ⇒ X=0.1
    cell.m_solid_auf_in[-1, S_CHAR] = 0.4
    cell.m_solid_ab_in[-1, S_CHAR] = 0.6
    cell.m_solid[-1, S_CHAR] = 1.0
    cell.K_solid_auf[-1, S_CHAR] = 0.5
    cell.K_solid_ab[-1, S_CHAR] = 0.4
    # out = (0.5+0.4)*1.0 = 0.9

    assert cell._compute_char_conversion() == pytest.approx(0.1)


def test_legacy_stream_char_conversion_ignores_auf_ab() -> None:
    """legacy_stream 仍只用 zu/in/rez；auf/ab 不参与（与旧口径兼容）。"""
    cell = _minimal_holdup_cell()
    cell.solid_state_model = "legacy_stream"
    cell.m_solid_zu[-1, S_CHAR] = 1.0
    cell.m_solid_auf_in[-1, S_CHAR] = 5.0  # 应被忽略
    cell.m_solid[-1, S_CHAR] = 0.8  # legacy out = m_solid
    assert cell._compute_char_conversion() == pytest.approx(0.2)


def test_holdup_without_auf_ab_still_uses_feed() -> None:
    """bed0 仅有 zu 时与旧式一致。"""
    cell = _minimal_holdup_cell()
    cell.m_solid_zu[-1, S_CHAR] = 1.0
    cell.m_solid[-1, S_CHAR] = 1.0
    cell.K_solid_auf[-1, S_CHAR] = 0.9
    cell.K_solid_ab[-1, S_CHAR] = 0.0
    # out = 0.9
    assert cell._compute_char_conversion() == pytest.approx(0.1)
