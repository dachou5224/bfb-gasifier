"""feed_inlet：ER + 元素分析 → O2/H2O/N2 摩尔流率。"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.core.feed_inlet import compute_gas_feeds_mol_s, stoichiometric_o2_mol_s


def test_stoich_o2_decreases_with_fuel_o():
    """燃料氧降低外供化学计量氧。"""
    base = stoichiometric_o2_mol_s(1.0, 10.0, 60.0, 5.0, 5.0, 0.0)
    with_o = stoichiometric_o2_mol_s(1.0, 10.0, 60.0, 5.0, 25.0, 0.0)
    assert with_o < base


def test_air_steam_n2_scales_like_air():
    """air_steam：N2/O2 ≈ 79/21。"""
    o2, h2o, n2 = compute_gas_feeds_mol_s(
        fuel_feed_kg_s=1.0,
        moisture_wt=10.0,
        C_dry=60.0,
        H_dry=5.0,
        O_dry=20.0,
        ER=0.3,
        primary_agent="air_steam",
        S_dry=0.0,
        steam_to_o2_molar=0.8,
    )
    assert abs(n2 / o2 - (0.79 / 0.21)) < 1e-9
    assert abs(h2o / o2 - 0.8) < 1e-9


def test_o2_steam_small_n2():
    o2, h2o, n2 = compute_gas_feeds_mol_s(
        fuel_feed_kg_s=1.0,
        moisture_wt=10.0,
        C_dry=60.0,
        H_dry=5.0,
        O_dry=20.0,
        ER=0.3,
        primary_agent="o2_steam",
        S_dry=0.0,
        steam_to_o2_molar=1.0,
        o2_steam_n2_frac_of_o2=0.01,
    )
    assert abs(n2 / o2 - 0.01) < 1e-9


if __name__ == "__main__":
    test_stoich_o2_decreases_with_fuel_o()
    test_air_steam_n2_scales_like_air()
    test_o2_steam_small_n2()
    print("test_feed_inlet: OK")
