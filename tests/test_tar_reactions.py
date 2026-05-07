from __future__ import annotations

import math

from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    get_tar_component_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)


def test_lumped_stoichiometry_contains_tar_consumption_terms():
    s = get_lumped_tar_stoichiometry("R10", fuel_type="coal")
    assert s["TAR1"] < 0.0
    assert s["TAR2"] < 0.0
    assert s["CO"] > 0.0
    assert s["H2"] > 0.0


def test_component_stoichiometry_has_unit_tar_consumption():
    s1 = get_tar_component_stoichiometry("R11", fuel_type="coal", component="TAR1")
    s2 = get_tar_component_stoichiometry("R11", fuel_type="coal", component="TAR2")
    assert s1["TAR1"] == -1.0
    assert s2["TAR2"] == -1.0
    assert "CO" in s1 and "CO" in s2


def test_rate_r10_increases_with_oxygen():
    r_low = rate_R10(T=1200.0, C_tar=0.2, C_O2=0.01, P=2.5e6, fuel_type="coal")
    r_high = rate_R10(T=1200.0, C_tar=0.2, C_O2=0.08, P=2.5e6, fuel_type="coal")
    assert r_low >= 0.0
    assert r_high > r_low


def test_rate_r11_paths_finite_and_nonnegative():
    rb = rate_R11_bubble(T=1150.0, C_tar=0.15)
    rs = rate_R11_suspension(T=1150.0, C_tar=0.15, rho_cat=800.0)
    assert rb >= 0.0 and math.isfinite(rb)
    assert rs >= 0.0 and math.isfinite(rs)
