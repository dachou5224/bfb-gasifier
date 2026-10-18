#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""按 ADR 0013 比较 Bild 3.5 / Bild 3.3 上的 6 个读法。

6 个组合只在本工具里：式 (3.6) 的三种 α_b，乘式 (3.43) 的两种压力系数。
hamel.hydrodynamics 里没有这些开关。判定规则：

- E = sqrt(mean(((d_calc − d_read)/σ_d)^2))，E′ 对 u_b,r 同样定义；
- E ≤ 3 为合格；合格组合里 E ≤ E_min + 1 的全部进入 C；
- C 只有一个就选定；否则对 C 中每一个用同样规则看 E′，得到 C′；
- C′ 只有一个就选定。没有合格组合、C′ 仍不唯一、或某个组合无解，都不选取。

图注没有给出 n_0、h_A、容器直径、ε_mf、φ_S。这些量由调用方传入，本工具不设默认值。
读数在 data/digitized/fig_3_5.json 与 fig_3_3.json。Bild 3.3 里分不开的点是 null，不进入 E′。
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from hamel.hydrodynamics.bed import (
    BedConditions,
    HydrodynamicsError,
    alpha_b,
    bubble_growth_rate,
    bubble_lifetime,
    bubble_rise_velocity,
    coalescence_parameter,
    growth_start_height,
    initial_bubble_diameter,
    integrate_growth,
    minimum_fluidization_velocity,
    psi_b,
    richardson_zaki_exponent,
    settling_reynolds,
    suspension_dynamic_viscosity,
    suspension_kinematic_viscosity,
    suspension_superficial_velocity,
    suspension_voidage,
    terminal_velocity,
    throughflow_velocity,
    upsilon_b,
    visible_bubble_fraction,
    P_0,
)

ALPHA_READINGS = ("printed", "ratio", "interstitial")
PRESSURE_READINGS = ("printed_minus_one", "power")
QUALIFIED = 3.0
BAND = 1.0


def combinations():
    """式 (3.6) 三种读法 × 式 (3.43) 两种压力系数。"""
    return [(alpha, pressure) for alpha in ALPHA_READINGS for pressure in PRESSURE_READINGS]


def _alpha(reading, u_b, u_d, epsilon_b, epsilon_d):
    if reading == "printed":
        return (u_b / u_d) * epsilon_b
    if reading == "ratio":
        return u_b / u_d
    if reading == "interstitial":
        return alpha_b(u_b, u_d, epsilon_d)
    raise ValueError(f"未知的 α_b 读法：{reading}")


def _pressure_factor(reading, pressure):
    power = (pressure / P_0) ** 0.2
    if reading == "printed_minus_one":
        return power - 1.0
    if reading == "power":
        return power
    raise ValueError(f"未知的压力系数读法：{reading}")


def _variant_growth(u_0, d_b, u_mf, conditions, alpha_reading, pressure_reading):
    suspension = suspension_superficial_velocity(
        u_mf, u_0, conditions.distributor, conditions.geldart, conditions.three_dimensional
    )
    through = throughflow_velocity(suspension, conditions.pressure)
    terminal = terminal_velocity(conditions.d_p, conditions.nu_g, conditions.rho_p, conditions.rho_g)
    n_rz = richardson_zaki(terminal, conditions)
    voidage = suspension_voidage(conditions.epsilon_mf, suspension, u_mf, n_rz)
    kinematic = suspension_kinematic_viscosity(
        suspension_dynamic_viscosity(conditions.d_p), conditions.rho_g, voidage
    )
    factor = _pressure_factor(pressure_reading, conditions.pressure)
    excess = psi_b(conditions.distributor, conditions.geldart) * (u_0 - u_mf) * factor
    velocity = bubble_rise_velocity(
        excess, upsilon_b(conditions.vessel_diameter, conditions.geldart), d_b, kinematic
    )
    fraction = visible_bubble_fraction(u_0, suspension, velocity, through)
    if fraction <= 0.0 or fraction >= 1.0:
        raise HydrodynamicsError("ε_b 不在 (0, 1) 内")
    characteristic = _alpha(alpha_reading, velocity, suspension, fraction, voidage)
    coalescence = coalescence_parameter(characteristic)
    growth = bubble_growth_rate(
        d_b,
        fraction,
        coalescence,
        bubble_lifetime(u_mf, conditions.pressure),
        velocity,
        conditions.pressure,
    )
    return growth, through


def richardson_zaki(terminal, conditions):
    return richardson_zaki_exponent(
        settling_reynolds(terminal, conditions.d_p, conditions.rho_g, conditions.mu_g)
    )


def bubble_diameter(height_m, excess_velocity, pressure, alpha_reading, pressure_reading, model):
    """在图注高度上计算一个组合的 d_b。无解时抛出 HydrodynamicsError。"""
    rho_g, mu_g, nu_g = model["properties_at"](pressure)
    u_mf = minimum_fluidization_velocity(
        model["epsilon_mf"], nu_g, model["phi_s"], model["d_p"], model["rho_p"], rho_g
    )
    u_0 = u_mf + excess_velocity
    conditions = BedConditions(
        model["d_p"],
        model["epsilon_mf"],
        model["rho_p"],
        rho_g,
        mu_g,
        nu_g,
        pressure,
        model["diffusivity"],
        model["vessel_diameter"],
        "technical",
        "B",
        True,
    )
    d_b0 = initial_bubble_diameter("technical", u_0, model["area_m2"], model["n_0"], u_mf)
    h_0 = growth_start_height("technical", d_b0, model["h_a"])
    if height_m <= h_0:
        return d_b0

    def derivative(diameter):
        growth, _through = _variant_growth(
            u_0, diameter, u_mf, conditions, alpha_reading, pressure_reading
        )
        return growth

    return integrate_growth(d_b0, h_0, height_m, model["n_steps"], derivative)


def throughflow(excess_velocity, pressure, model):
    """式 (3.44) 的 u_b,r。它不依赖 α_b，也不依赖式 (3.43) 的减 1。"""
    rho_g, mu_g, nu_g = model["properties_at"](pressure)
    u_mf = minimum_fluidization_velocity(
        model["epsilon_mf"], nu_g, model["phi_s"], model["d_p"], model["rho_p"], rho_g
    )
    u_d = suspension_superficial_velocity(u_mf, u_mf + excess_velocity, "technical", "B", True)
    return throughflow_velocity(u_d, pressure)


def normalized_rms(calculated, read, sigma):
    if sigma <= 0.0:
        raise ValueError("读数误差必须大于 0")
    if len(calculated) != len(read) or not calculated:
        raise ValueError("计算值与读数必须非空且一一对应")
    squares = [((calc - observed) / sigma) ** 2 for calc, observed in zip(calculated, read)]
    return math.sqrt(sum(squares) / len(squares))


def decide(rows):
    """ADR 0013 第 1、2 项的判定。rows 含 name、E、E_prime、solvable。"""
    names = [row["name"] for row in rows]
    if len(names) != len(set(names)):
        raise ValueError("组合名称重复")
    failed = [row["name"] for row in rows if not row["solvable"]]
    if failed:
        return _undecided(rows, "某组合在图注条件下无解", failed)
    qualified = [row for row in rows if row["E"] <= QUALIFIED]
    if not qualified:
        return _undecided(rows, "没有合格组合")
    minimum = min(row["E"] for row in qualified)
    primary = [row for row in qualified if row["E"] <= minimum + BAND]
    if len(primary) == 1:
        return _selected(rows, primary, primary[0]["name"])
    secondary_pool = [row for row in primary if row["E_prime"] <= QUALIFIED]
    if not secondary_pool:
        return _undecided(rows, "辅图没有合格组合", [row["name"] for row in primary])
    prime_minimum = min(row["E_prime"] for row in secondary_pool)
    secondary = [row for row in secondary_pool if row["E_prime"] <= prime_minimum + BAND]
    if len(secondary) == 1:
        return _selected(rows, primary, secondary[0]["name"])
    return _undecided(rows, "辅图候选仍不唯一", [row["name"] for row in secondary])


def _selected(rows, primary, name):
    return {
        "status": "selected",
        "selected": name,
        "reason": "候选集合只剩一个组合",
        "primary": [row["name"] for row in primary],
        "rows": rows,
    }


def _undecided(rows, reason, pending=None):
    return {
        "status": "unresolved",
        "selected": None,
        "reason": reason,
        "primary": list(pending or []),
        "rows": rows,
    }


def compare(diameter_figure, throughflow_figure, model):
    """对图中各点计算 6 个组合的 E 与 E′，并给出候选集合。"""
    sigma_d = float(diameter_figure["est_error"]["bubble_diameter_m"])
    sigma_u = float(throughflow_figure["est_error"]["throughflow_velocity_m_s"])
    diameter_points = list(diameter_figure["points"])
    through_points = list(throughflow_figure["points"])
    rows = []
    for alpha_reading, pressure_reading in combinations():
        name = f"{alpha_reading}/{pressure_reading}"
        calculated_d = []
        calculated_u = []
        solvable = True
        try:
            for point in diameter_points:
                calculated_d.append(
                    bubble_diameter(
                        float(diameter_figure["height_m"]),
                        float(point["excess_velocity_m_s"]),
                        float(point["pressure_Pa"]),
                        alpha_reading,
                        pressure_reading,
                        model,
                    )
                )
            for point in through_points:
                calculated_u.append(
                    throughflow(
                        float(point["excess_velocity_m_s"]),
                        float(point["pressure_Pa"]),
                        model,
                    )
                )
        except (HydrodynamicsError, ValueError, ZeroDivisionError):
            solvable = False
        row = {
            "name": name,
            "alpha_reading": alpha_reading,
            "pressure_reading": pressure_reading,
            "solvable": solvable,
            "bubble_diameter_m": calculated_d,
            "throughflow_velocity_m_s": calculated_u,
        }
        if solvable:
            row["E"] = normalized_rms(
                calculated_d,
                [float(point["bubble_diameter_m"]) for point in diameter_points],
                sigma_d,
            )
            row["E_prime"] = normalized_rms(
                calculated_u,
                [float(point["throughflow_velocity_m_s"]) for point in through_points],
                sigma_u,
            )
        else:
            row["E"] = None
            row["E_prime"] = None
        rows.append(row)
    return decide(rows)


def load_figure(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def normalize_diameter_figure(document):
    """把 Bild 3.5 的 series 收成 compare() 用的点。已是点列时原样返回。"""
    if "series" not in document:
        return document
    points = []
    for series in document["series"]:
        excess = float(series["u0_minus_umf_m_s"])
        for pressure_mpa, diameter in series["points"]:
            if diameter is None:
                continue
            points.append(
                {
                    "excess_velocity_m_s": excess,
                    "pressure_Pa": float(pressure_mpa) * 1.0e6,
                    "bubble_diameter_m": float(diameter),
                }
            )
    return {
        "height_m": float(document["conditions"]["h_m"]),
        "est_error": {"bubble_diameter_m": float(document["est_error"]["d_b"])},
        "points": points,
    }


def normalize_throughflow_figure(document):
    """把 Bild 3.3 的 series 收成 compare() 用的点。null 表示线分不开，不进入 E′。"""
    if "series" not in document:
        return document
    points = []
    for series in document["series"]:
        excess = float(series["u0_minus_umf_m_s"])
        for pressure_mpa, velocity in series["points"]:
            if velocity is None:
                continue
            points.append(
                {
                    "excess_velocity_m_s": excess,
                    "pressure_Pa": float(pressure_mpa) * 1.0e6,
                    "throughflow_velocity_m_s": float(velocity),
                }
            )
    return {
        "est_error": {"throughflow_velocity_m_s": float(document["est_error"]["u_b_r"])},
        "points": points,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="比较 Bild 3.5 与 Bild 3.3 的 6 个读法。")
    parser.add_argument("--fig-3-5", type=Path, required=True)
    parser.add_argument("--fig-3-3", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True, help="模型条件 JSON，不含默认物性")
    args = parser.parse_args(argv)
    model = load_figure(args.model)
    if "properties_at" not in model and "pressures" not in model:
        raise SystemExit("模型 JSON 必须给出各压力下的 rho_g、mu_g、nu_g")
    model = _bind_properties(model)
    result = compare(
        normalize_diameter_figure(load_figure(args.fig_3_5)),
        normalize_throughflow_figure(load_figure(args.fig_3_3)),
        model,
    )
    printable = {
        "status": result["status"],
        "selected": result["selected"],
        "reason": result["reason"],
        "primary": result["primary"],
        "rows": [
            {
                "name": row["name"],
                "solvable": row["solvable"],
                "E": row["E"],
                "E_prime": row["E_prime"],
            }
            for row in result["rows"]
        ],
    }
    print(json.dumps(printable, ensure_ascii=False, indent=2))
    return 0


def _bind_properties(model):
    table = {float(item["pressure_Pa"]): item for item in model["pressures"]}

    def properties_at(pressure):
        item = table[float(pressure)]
        return float(item["rho_g"]), float(item["mu_g"]), float(item["nu_g"])

    bound = dict(model)
    bound["properties_at"] = properties_at
    return bound


if __name__ == "__main__":
    raise SystemExit(main())
