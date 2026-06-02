"""气化剂入口摩尔流率：由 ER、燃料元素分析与气化剂类型计算 O2/H2O/N2。

与 Hamel Table 2 / CASE_LU 一致：化学计量氧按干基燃料完全氧化至 CO2/H2O/SO2，
再按当量比 ER 缩放；**空气** 工况下 N2 随 O2 按体积比 79/21 配平。

Source: ``docs/hamel_submodels/00_readme_and_citation_rules.md``;
``data/validation_cases.json`` CASE_HTW_WESSELING_1（Table 7.1 / LU）
"""

from __future__ import annotations

from typing import Literal

# 干空气中 O2、N2 摩尔分数（用于 air / air_steam）
X_O2_AIR: float = 0.21
X_N2_AIR: float = 0.79

PrimaryAgent = Literal["air_steam", "o2_steam", "air"]


def stoichiometric_o2_mol_s(
    fuel_feed_kg_s: float,
    moisture_wt: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    S_dry: float = 0.0,
) -> float:
    """完全氧化至 CO2、H2O(l)、SO2 所需外供 O2 [mol/s]。

    干基：C、H、O、S 为质量分数 [%]；燃料中氧减少外供 O2：
    n_O2,ext = max(0, n_C + n_H2/2 + n_S - n_O/2)，其中 n_O 为燃料中氧原子摩尔流量。

    Parameters
    ----------
    fuel_feed_kg_s
        湿基燃料质量流率 [kg/s]
    moisture_wt
        湿基水分 [wt%]
    C_dry, H_dry, O_dry, S_dry
        干基元素质量分数 [wt%]
    """
    moisture_frac = max(0.0, min(1.0, moisture_wt / 100.0))
    dry_fuel = max(fuel_feed_kg_s * (1.0 - moisture_frac), 0.0)

    M_C = 12.011e-3
    M_H2 = 2.016e-3
    M_O = 15.999e-3
    M_S = 32.065e-3

    C_mass = dry_fuel * C_dry / 100.0
    H_mass = dry_fuel * H_dry / 100.0
    O_mass = dry_fuel * O_dry / 100.0
    S_mass = dry_fuel * S_dry / 100.0

    n_C = C_mass / M_C
    n_H2 = H_mass / M_H2
    n_O = O_mass / M_O
    n_S = S_mass / M_S

    # C + O2 -> CO2 ; H2 + 1/2 O2 -> H2O ; S + O2 -> SO2 ；燃料 O 为原子氧
    n_o2 = n_C + 0.5 * n_H2 + n_S - 0.5 * n_O
    return max(n_o2, 1e-30)


def compute_gas_feeds_mol_s(
    fuel_feed_kg_s: float,
    moisture_wt: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    ER: float,
    primary_agent: PrimaryAgent | str = "air_steam",
    S_dry: float = 0.0,
    steam_to_o2_molar: float = 0.8,
    o2_steam_n2_frac_of_o2: float = 0.01,
) -> tuple[float, float, float]:
    """由 ER 与燃料分析计算入口 O2、H2O、N2 摩尔流率 [mol/s]。

    Parameters
    ----------
    ER
        当量比（与化学计量氧之比，Hamel Table 2 / CASE_LU 常用记号）
    primary_agent
        - **air_steam**：O2 来自干空气；N2 = (X_N2/X_O2)*O2_feed；蒸汽单独 H2O = steam_to_o2_molar * O2_feed
        - **o2_steam**：纯氧 + 蒸汽；N2 取 o2_steam_n2_frac_of_o2 * O2（载气/微量）
        - **air**：仅干空气，H2O_feed = steam_to_o2_molar * O2_feed 仍可表示雾化/漏入蒸汽（默认可调 0）
    steam_to_o2_molar
        蒸汽与氧的摩尔比 H2O/O2（HTW 等需与工况表对齐时可改）
    o2_steam_n2_frac_of_o2
        **o2_steam** 时 N2_feed = 该因子 * O2_feed（旧占位 0.01 量级）

    Returns
    -------
    (O2_feed, H2O_feed, N2_feed) 均为 [mol/s]
    """
    n_o2_stoich = stoichiometric_o2_mol_s(
        fuel_feed_kg_s, moisture_wt, C_dry, H_dry, O_dry, S_dry
    )
    o2 = ER * n_o2_stoich

    agent = str(primary_agent).lower().replace("-", "_")
    if agent in ("air_steam", "airsteam", "luft_dampf"):
        n2 = o2 * (X_N2_AIR / X_O2_AIR)
        h2o = steam_to_o2_molar * o2
    elif agent in ("o2_steam", "oxygen_steam", "o2steam"):
        n2 = o2 * o2_steam_n2_frac_of_o2
        h2o = steam_to_o2_molar * o2
    elif agent == "air":
        n2 = o2 * (X_N2_AIR / X_O2_AIR)
        h2o = steam_to_o2_molar * o2
    else:
        raise ValueError(
            f"未知气化剂类型 primary_agent={primary_agent!r}，"
            f"请使用 air_steam | o2_steam | air"
        )

    return (float(o2), float(h2o), float(n2))
