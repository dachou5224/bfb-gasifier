from __future__ import annotations

import pytest

from src.physics.mass_transfer import calc_kbd, calc_u_br


def test_calc_u_br_decreases_with_pressure():
    u_d = 0.12
    u_low = calc_u_br(u_d=u_d, P=101_325.0)
    u_high = calc_u_br(u_d=u_d, P=2_500_000.0)
    assert u_low > 0.0
    assert u_high > 0.0
    assert u_high < u_low


def test_calc_kbd_positive_and_increases_with_u_br():
    base = calc_kbd(u_br=0.03, d_b=0.08, D_g=2.0e-5, eps_mf=0.45, u_b=0.9)
    higher = calc_kbd(u_br=0.06, d_b=0.08, D_g=2.0e-5, eps_mf=0.45, u_b=0.9)
    assert base > 0.0
    assert higher > base


def test_calc_kbd_invalid_epsmf_raises():
    with pytest.raises(AssertionError):
        calc_kbd(u_br=0.03, d_b=0.08, D_g=2.0e-5, eps_mf=1.1, u_b=0.9)
