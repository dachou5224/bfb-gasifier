"""首次运行前的报告记账、闭式比对和数值加密。

不进入残差，不改冻结值，不设通过门。依据 ADR 0040 第 5、6 项。
"""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError, account_freeboard, cell_state, discriminant, ejection_velocity_samples,
    motion_coefficients, terminal_settling_velocity,
)
from hamel.properties.data import GAS_SPECIES, atomic_weight, atoms, molar_mass_kg_mol
from hamel.properties.thermo import density_kg_m3, volume_flow_m3_s
from hamel.properties.transport import mixture_viscosity_pa_s
from hamel.pyrolysis import lignite_particle_properties, particle_release_times
from hamel.solver.gates import block_ratios
from hamel.solver.residual import _bed_ash_by_class, evaluate_residual, particle_temperatures
from hamel.solver.vector import pack_state

# T23 单元测试里 Δ > 0 的常气速构造。气速为 0，所以速度尺度不用 u_g。
_BENCHMARK = {
    "diameter_m": 8.0e-3, "rho_g": 1.1, "mu_g": 1.8e-5, "rho_p": 900.0, "phi_s": 0.8,
    "gas_m_s": 0.0, "initial_m_s": 2.0, "height_m": 0.3, "duration_s": 0.02,
}
_ACCEPTANCE = ("T_FRB_K", "X_C", "CO", "CO2", "H2", "CH4", "H2O")
_T_FRB = "Σ T_i·Δh_i / Σ Δh_i，自由空间全部 cell，含锥段。解释性选择：ADR 0040 第 1 项。"
_X_C = "出口气体（含焦油蒸气）中的碳 / 燃料碳。分子不计带出焦炭与排灰。解释性选择：ADR 0040 第 1 项。"


def fuel_element_fractions(physical) -> dict:
    """1 kg 收到基燃料里的 C、H、O 质量。与出口碳转化率的分母同一套元素。"""
    moisture = float(physical.moisture_mass_fraction)
    dry = 1.0 - moisture
    water = molar_mass_kg_mol("H2O")
    hydrogen = dry * float(physical.hydrogen_mass_fraction_wf)
    hydrogen += 2.0 * moisture * atomic_weight("H") * 1.0e-3 / water
    oxygen = dry * float(physical.oxygen_mass_fraction_wf)
    oxygen += moisture * atomic_weight("O") * 1.0e-3 / water
    return {
        "C": dry * float(physical.carbon_mass_fraction_wf),
        "H": hydrogen,
        "O": oxygen,
    }


def dry_element_fractions(physical) -> dict:
    """1 kg 干料（已释放水分的 P）里的 C、H、O。不含水分带入的氢和氧。"""
    return {
        "C": float(physical.carbon_mass_fraction_wf),
        "H": float(physical.hydrogen_mass_fraction_wf),
        "O": float(physical.oxygen_mass_fraction_wf),
    }


def _blank() -> dict:
    return {"mass": 0.0, "C": 0.0, "H": 0.0, "O": 0.0}


def _add(total, part) -> None:
    for key in total:
        total[key] += part[key]


def _gas(flow) -> dict:
    packet = _blank()
    for index, name in enumerate(GAS_SPECIES):
        moles = float(flow[index])
        if moles == 0.0:
            continue
        packet["mass"] += moles * molar_mass_kg_mol(name)
        for element, count in atoms(name).items():
            if element in packet:
                packet[element] += moles * count * atomic_weight(element) * 1.0e-3
    return packet


def _of(kind: str, kilograms) -> dict:
    mass = float(np.sum(kilograms))
    packet = _blank()
    packet["mass"] = mass
    if kind == "char":
        packet["C"] = mass
    return packet


def _fuel_packet(mass_w, mass_p, received, dried) -> dict:
    """W 按收到基，P 的干料质量按干基。两者不再乘同一套分数。"""
    wet = float(np.sum(mass_w))
    dry = float(np.sum(mass_p))
    packet = _blank()
    packet["mass"] = wet + dry
    for element in ("C", "H", "O"):
        packet[element] = wet * float(received[element]) + dry * float(dried[element])
    return packet


def _split_mass(coeff, undried, dried_stock, moisture: float):
    undried = np.asarray(undried, dtype=float)
    dried_stock = np.asarray(dried_stock, dtype=float)
    return coeff * undried, coeff * (1.0 - moisture) * dried_stock


def fuel_parts(frozen, reported, physical) -> dict:
    """边界上的 W 与 P。P 的质量是 (1−水分) 乘当量库存，已不含水分。"""
    moisture = float(reported["moisture_mass_fraction"])
    received, dried = fuel_element_fractions(physical), dry_element_fractions(physical)
    up = frozen.a_free_per_s[0]
    down = frozen.b_free_per_s[-1]
    up_w, up_p = _split_mass(up, reported["undried_free_kg"][0], reported["dried_free_kg"][0], moisture)
    down_w, down_p = _split_mass(down, reported["undried_free_kg"][-1], reported["dried_free_kg"][-1], moisture)
    bed_w, bed_p = _split_mass(frozen.k_auf_per_s[0], reported["undried_bed_kg"][0], reported["dried_bed_kg"][0], moisture)
    redirect_w = np.zeros(frozen.n_class)
    redirect_p = np.zeros(frozen.n_class)
    if frozen.hover_closure in {"B", "C"}:
        # 进入悬停状态的流：c↑、c↓ 由自由状态转入。ADR 0048 第 2 项。
        for klass in range(frozen.n_class):
            for cell in range(frozen.n_freeboard):
                for coeff in (frozen.c_up_per_s[cell, klass], frozen.c_down_per_s[cell, klass]):
                    piece_w, piece_p = _split_mass(
                        coeff, reported["undried_free_kg"][cell, klass],
                        reported["dried_free_kg"][cell, klass], moisture,
                    )
                    redirect_w[klass] += piece_w
                    redirect_p[klass] += piece_p
    return {
        "received": received, "dried": dried,
        "bed_w": bed_w, "bed_p": bed_p, "up_w": up_w, "up_p": up_p,
        "down_w": down_w, "down_p": down_p, "redirect_w": redirect_w, "redirect_p": redirect_p,
    }


def _mobile_ash(state, frozen, cell) -> np.ndarray:
    """自由空间的输运灰。ADR 0044 起灰不因 m_C = 0 停止输运。"""
    del frozen
    return np.array(state.freeboard_ash_free_kg[cell], dtype=float, copy=True)


def _combine(*parts) -> dict:
    total = _blank()
    for part in parts:
        _add(total, part)
    return total


def freeboard_closure(state, frozen, assembly, fuel) -> dict:
    """自由空间全部相的边界流。反应不另计。B、C 的改送是出流。

    F^ej,t 在抛出入流里计一次，在改送出流里再计一次。ADR 0040 第 5 项。
    燃料的 W 与 P 分开记账。P 已释放水分，不用收到基组成。
    """
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    share = _bed_ash_by_class(state.bed_char_kg, inert)
    thrown_char = frozen.k_auf_per_s[0] * state.bed_char_kg[0]
    thrown_ash = frozen.k_auf_per_s[0] * share[0]
    received, dried = fuel["received"], fuel["dried"]
    gas_in = _gas(state.bed_gas_mol_s[0, 0] + state.bed_gas_mol_s[0, 1])
    side_in = _gas(frozen.secondary_feed_mol_s)
    solids_in = _combine(
        _of("char", thrown_char), _of("ash", thrown_ash),
        _fuel_packet(fuel["bed_w"], fuel["bed_p"], received, dried),
    )
    gas_out = _gas(state.freeboard_gas_mol_s[0])
    up = frozen.a_free_per_s[0]
    down = frozen.b_free_per_s[-1]
    carried = _combine(
        _of("char", up * state.freeboard_char_free_kg[0]),
        _of("ash", up * _mobile_ash(state, frozen, 0)),
        _fuel_packet(fuel["up_w"], fuel["up_p"], received, dried),
    )
    returned = _combine(
        _of("char", down * state.freeboard_char_free_kg[-1]),
        _of("ash", down * _mobile_ash(state, frozen, -1)),
        _fuel_packet(fuel["down_w"], fuel["down_p"], received, dried),
    )
    redirected = _combine(
        _of("char", assembly.redirect_char_kg_s),
        _of("ash", assembly.redirect_ash_kg_s),
        _fuel_packet(fuel["redirect_w"], fuel["redirect_p"], received, dried),
    )
    names = (
        ("床区最上一格气体", "入", gas_in),
        ("二次风及其它进入自由空间的气体", "入", side_in),
        ("床面抛出固体，含 F^ej,f 与 F^ej,t", "入", solids_in),
        ("自由空间出口气体", "出", gas_out),
        ("普通带出固体", "出", carried),
        ("普通回落固体", "出", returned),
        (f"{frozen.hover_closure} 的改送流", "出", redirected),
    )
    incoming = _blank()
    outgoing = _blank()
    streams = []
    for name, direction, packet in names:
        _add(incoming if direction == "入" else outgoing, packet)
        streams.append({"name": name, "direction": direction, **packet})
    totals = {}
    for key in incoming:
        gap = incoming[key] - outgoing[key]
        totals[key] = {
            "in": incoming[key], "out": outgoing[key], "difference": gap,
            "relative": None if incoming[key] == 0.0 else gap / incoming[key],
        }
    hover_char = float(
        np.sum((frozen.c_up_per_s + frozen.c_down_per_s) * state.freeboard_char_free_kg)
    )
    return {
        "streams": streams,
        "totals": totals,
        "reaction_char_kg_s": float(np.sum(assembly.char_consumption_free_kg_s) + np.sum(assembly.char_consumption_hover_kg_s)),
        "hover_inflow_char_kg_s": hover_char,
        "redirect_char_kg_s": float(assembly.redirect_char_kg_s),
    }


def _hold(reported, name, cell):
    values = (reported or {}).get(name)
    if values is None:
        return None
    return np.asarray(values[cell], dtype=float)


def _fuel_hold(reported, cell, kind: str) -> np.ndarray:
    """未干燥 W 加上已干燥 P 的实际干料质量。没有库存时这一项为 0。"""
    undried = _hold(reported, f"undried_{kind}_kg", cell)
    dried = _hold(reported, f"dried_{kind}_kg", cell)
    width = 0 if undried is None else undried.size
    if undried is None or dried is None:
        return np.zeros(width)
    moisture = float((reported or {}).get("moisture_mass_fraction") or 0.0)
    return undried + (1.0 - moisture) * dried


def _redirect_class(assembly, name: str, klass: int) -> float:
    values = getattr(assembly, f"redirect_{name}_by_class_kg_s", None)
    if values is None:
        return 0.0
    return float(values[klass])


def stock_report(state, frozen, assembly, reported) -> dict:
    """收敛存量，以及存量乘输运系数得到的回落和带出。

    B 的改送计入回落，C 的改送计入带出。顶格 a 与底格 b 已在普通边界流里，不再加一次。
    """
    classes = []
    for klass in range(frozen.n_class):
        upward = float(frozen.a_free_per_s[0, klass] + frozen.a_hover_per_s[0, klass])
        downward = float(frozen.b_free_per_s[-1, klass] + frozen.b_hover_per_s[-1, klass])
        top = float(state.freeboard_char_free_kg[0, klass] + _mobile_ash(state, frozen, 0)[klass])
        bottom = float(state.freeboard_char_free_kg[-1, klass] + _mobile_ash(state, frozen, -1)[klass])
        redirected = sum(_redirect_class(assembly, name, klass) for name in ("char", "ash", "fuel"))
        returned = downward * bottom + float(frozen.fuel_down_free_kg_s[-1, klass])
        carried = upward * top + float(frozen.fuel_up_free_kg_s[0, klass])
        if frozen.hover_closure == "B":
            returned += redirected
        elif frozen.hover_closure == "C":
            carried += redirected
        classes.append({
            "class": klass + 1,
            "reference_return_share": _listed(reported, "reference_return_share", klass),
            "reference_elutriation_share": _listed(reported, "reference_elutriation_share", klass),
            "return_kg_s": returned,
            "elutriation_kg_s": carried,
        })
    cells = []
    hover = []
    for cell in range(frozen.n_freeboard):
        free_fuel = _fuel_hold(reported, cell, "free")
        hover_fuel = _fuel_hold(reported, cell, "hover")
        holds = []
        for klass in range(frozen.n_class):
            free = float(state.freeboard_char_free_kg[cell, klass] + state.freeboard_ash_free_kg[cell, klass])
            held = float(state.freeboard_char_hover_kg[cell, klass] + state.freeboard_ash_hover_kg[cell, klass])
            if free_fuel.size:
                free += float(free_fuel[klass])
            if hover_fuel.size:
                held += float(hover_fuel[klass])
            holds.append((free + held, held))
        solid = float(sum(item[0] for item in holds))
        consumed = float(assembly.char_consumption_free_kg_s[cell].sum() + assembly.char_consumption_hover_kg_s[cell].sum())
        cells.append({"cell": cell, "solid_kg": solid, "char_consumption_kg_s": consumed})
        for klass, (_total, held) in enumerate(holds):
            hover.append({
                "cell": cell, "class": klass + 1, "hover_kg": held, "cell_solid_kg": solid,
                "hover_fraction": None if solid == 0.0 else held / solid,
            })
    return {"classes": classes, "cells": cells, "hover": hover}


def _listed(reported, name, index):
    values = (reported or {}).get(name)
    if values is None or index >= len(values):
        return None
    return float(values[index])


def used_values(frozen, reported, measurements) -> dict:
    """这一次运行实际用到的冻结值。口径两项每次都回显。"""
    measured = measurements or {}
    reported = reported or {}
    return {
        "epsilon_mf": float(frozen.epsilon_mf),
        "phi_s": reported.get("phi_s"),
        "n_0": reported.get("n_0"),
        "epsilon_p": float(frozen.particle_emissivity),
        "serio_frequency_per_s": float(frozen.serio_frequency_per_s),
        "serio_activation_j_mol": float(frozen.serio_activation_j_mol),
        "hover_closure": frozen.hover_closure,
        "cell_division": reported.get("cell_division"),
        "k_w_w_m2_k": float(frozen.k_w_w_m2_k),
        "rho_p_kg_m3": float(frozen.rho_p_kg_m3),
        "T_FRB_definition": _T_FRB,
        "X_C_definition": _X_C,
        "T_FRB_K": measured.get("T_FRB_K"),
        "X_C": measured.get("X_C"),
    }


def attach_converged(outcome, state, frozen, assembly, physical, reported, config, refs) -> None:
    """收敛之后把记账、闭式比对和加密放进诊断。求解已经结束，这里不改选。"""
    reported = reported or {}
    fuel = fuel_parts(frozen, reported, physical) if "undried_free_kg" in reported else None
    outcome.diagnostics["stocks"] = stock_report(state, frozen, assembly, reported)
    if fuel is not None:
        outcome.diagnostics["freeboard_closure"] = freeboard_closure(state, frozen, assembly, fuel)
    outcome.diagnostics["used"] = used_values(frozen, reported, outcome.measurements)
    outcome.diagnostics["bed_gas_temperature_deviation_k"] = reported.get("bed_gas_temperature_deviation_k")
    if "trajectory" in reported and "pyrolysis" in reported:
        outcome.diagnostics["closed_form"] = closed_form_from_state(state, frozen, reported, config)
        outcome.diagnostics["refinement"] = refinement_from_case(config, reported, state, frozen, refs)


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def _rk4(a, b, c, u0, h0, duration_s, steps=20000):
    """常系数的参考积分。只用于与闭式比对，不进入求解。求解手段：ADR 0047 的测试。"""
    step = float(duration_s) / steps
    velocity, height = float(u0), float(h0)
    f = lambda value: a * value * value + b * value + c
    for _ in range(steps):
        k1 = f(velocity)
        k2 = f(velocity + 0.5 * step * k1)
        k3 = f(velocity + 0.5 * step * k2)
        k4 = f(velocity + step * k3)
        l1 = velocity
        l2 = velocity + 0.5 * step * k1
        l3 = velocity + 0.5 * step * k2
        l4 = velocity + step * k3
        velocity += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        height += step * (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
    return velocity, height


def fixed_closed_form(samples: int) -> dict:
    """常气速、Δ > 0 的固定构造，与高精度 RK4 比。样本数不改变这一条轨迹。"""
    item = _BENCHMARK
    coefficient = motion_coefficients(item["gas_m_s"], True, item["rho_g"], item["mu_g"], item["diameter_m"], item["rho_p"], item["phi_s"])
    delta = discriminant(*coefficient)
    analytic_u, analytic_h = cell_state(item["duration_s"], *coefficient, item["initial_m_s"], item["height_m"])
    numeric_u, numeric_h = _rk4(*coefficient, item["initial_m_s"], item["height_m"], item["duration_s"], 200000)
    velocity = _relative(numeric_u, analytic_u)
    height = _relative(numeric_h, analytic_h)
    return {
        "delta": delta, "samples": int(samples),
        "velocity_relative": velocity, "height_relative": height, "max_relative": max(velocity, height),
    }


def applicability_row(diameter_m, rho_g, mu_g, rho_p, phi_s, gas_m_s, bubble_m_s, height_m, time_limit_s) -> dict:
    """一个粒级。Δ ≤ 0 记 N/A。Δ > 0 时在常系数下与高精度 RK4 比到时间上限。"""
    sample_speeds, _shares, direct = ejection_velocity_samples(float(bubble_m_s), 1)
    initial = float(sample_speeds[0])
    branches = {}
    for name, faster in (("faster", True), ("slower", False)):
        coefficient = motion_coefficients(gas_m_s, faster, rho_g, mu_g, diameter_m, rho_p, phi_s)
        branches[name] = (discriminant(*coefficient), coefficient)
    branch = "faster" if initial > gas_m_s else "slower"
    delta, coefficient = branches[branch]
    row = {
        "diameter_m": float(diameter_m), "branch": branch, "delta": delta,
        "delta_faster": branches["faster"][0], "delta_slower": branches["slower"][0],
        "initial_velocity_m_s": initial, "direct_return_share": float(direct),
        "status": "N/A" if delta <= 0.0 else "compared",
        "stopped": None, "max_velocity_relative": None, "max_height_relative": None,
    }
    if delta <= 0.0 or gas_m_s == 0.0 or height_m <= 0.0:
        if delta > 0.0 and (gas_m_s == 0.0 or height_m <= 0.0):
            row["stopped"] = "scale"
        return row
    duration = float(time_limit_s)
    try:
        closed_u, closed_h = cell_state(duration, *coefficient, initial, 0.0)
    except FreeboardError:
        row["stopped"] = "real_domain"
        return row
    numeric_u, numeric_h = _rk4(*coefficient, initial, 0.0, duration)
    row["stopped"] = "branch_change" if (initial - gas_m_s) * (numeric_u - gas_m_s) <= 0.0 else "time_limit"
    row["max_velocity_relative"] = abs(numeric_u - closed_u) / abs(gas_m_s)
    row["max_height_relative"] = abs(numeric_h - closed_h) / height_m
    return row


def toy_applicability_rows() -> list:
    """玩具用的两个粒级。一个预期 Δ > 0，一个用更小的颗粒去碰 Δ ≤ 0。"""
    return [
        (8.0e-3, 1.1, 1.8e-5, 900.0, 0.8, 1.5, 1.0, 0.4),
        (1.5e-4, 1.2, 2.0e-5, 2000.0, 1.0, 0.5, 1.0, 0.5),
    ]


def closed_form_report(config: SolverConfig, rows: list) -> dict:
    limit = float(config.freeboard_time_limit_s)
    table = [
        applicability_row(*item, limit) if not isinstance(item, dict) else item
        for item in rows
    ]
    return {"benchmark": fixed_closed_form(int(config.freeboard_ejection_samples)), "classes": table}


def _velocity(column):
    gas = column.get("gas")
    if callable(gas):
        return gas
    if column.get("gas_m_s") is not None:
        speed = float(column["gas_m_s"])
        return lambda _cell, _height: speed
    from hamel.freeboard.transport import ghost_bubble_velocity

    edges, diameters = column["edges_m"], column["edge_diameter_m"]
    flows = column["volume_flow_m3_s"]
    bubble, width = float(column["bubble_m_s"]), float(column["bubble_diameter_m"])
    return lambda cell, height: ghost_bubble_velocity(cell, height, edges, diameters, flows, bubble, width)


def _books(samples: int, column, config: SolverConfig):
    tuned = replace(config, freeboard_ejection_samples=int(samples))
    speeds, shares, direct = ejection_velocity_samples(float(column["bubble_m_s"]), int(samples))
    rho = np.asarray(column["rho_g"], dtype=float)
    mu = np.asarray(column["mu_g"], dtype=float)
    return account_freeboard(
        np.asarray(column["edges_m"], dtype=float), _velocity(column), rho, mu,
        float(column["diameter_m"]), float(column["rho_p"]), float(column["phi_s"]),
        speeds, shares, direct, tuned,
    )


def toy_trajectory() -> dict:
    """玩具柱。气速取终端速度加 3 m/s，使系数存在。换直径或气速会改变 a、b。"""
    diameter, rho_g, mu, rho_p, phi = 2.0e-4, 1.2, 2.0e-5, 1200.0, 0.8
    settling = terminal_settling_velocity(rho_g, mu, diameter, rho_p, phi)
    return {
        "edges_m": [0.0, 0.25, 0.5], "gas_m_s": settling + 3.0, "bubble_m_s": 1.0,
        "rho_g": [rho_g, rho_g], "mu_g": [mu, mu], "diameter_m": diameter,
        "rho_p": rho_p, "phi_s": phi,
    }


def toy_particle() -> dict:
    return {
        "diameter_m": 2.0e-4, "temperature_k": 1100.0, "pressure_pa": 1.0e5,
        "moisture_mass_fraction": 0.02, "feed_temperature_k": 298.15,
        "gas": {"N2": 0.79, "O2": 0.21}, "relative_velocity_m_s": 0.2,
        "ultimate_volatile_mass_fraction": 0.4,
    }


def _book_change(base, fine) -> tuple[float, bool]:
    """两边都没有的系数是同一种自由状态缺失，不比较。只有一边缺失才记存在性变化。"""
    worst = 0.0
    changed = False
    for left, right in zip(base.cells, fine.cells):
        for name in ("a_f", "b_f", "c_up", "c_down"):
            first, second = getattr(left, name), getattr(right, name)
            if first is None and second is None:
                continue
            if first is None or second is None:
                changed = True
                continue
            worst = max(worst, _relative(float(first), float(second)))
    return worst, changed


def _coefficients(book) -> list:
    return [None if cell.a_f is None else float(cell.a_f) for cell in book.cells]


def refine_trajectory(config: SolverConfig, column: dict) -> dict:
    """样本数加倍。门是 a、b 的相对变化 1e-4。ADR 0016、ADR 0047。

    column 是这一次运行的几何、气体和粒径。换这些输入会改变 a、b。
    """
    samples = int(config.freeboard_ejection_samples)
    base = _books(samples, column, config)
    fine = _books(samples * 2, column, config)
    change, existence_changed = _book_change(base, fine)
    return {
        "samples": samples, "diameter_m": float(column["diameter_m"]),
        "relative_change": change, "a_f": _coefficients(base), "existence_changed": existence_changed,
        "gate": 1.0e-4, "met": (not existence_changed) and change < 1.0e-4,
    }


def refine_pyrolysis(config: SolverConfig, particle: dict) -> dict:
    """径向加倍、时间容差收到十分之一、活化能节点增到 24。门 1e-3。ADR 0027。"""
    material = particle.get("material") or lignite_particle_properties()
    args = (
        float(particle["diameter_m"]), float(particle["temperature_k"]), float(particle["pressure_pa"]),
        float(particle["moisture_mass_fraction"]), float(particle["feed_temperature_k"]), dict(particle["gas"]),
        float(particle["relative_velocity_m_s"]), float(particle["ultimate_volatile_mass_fraction"]), material,
    )
    base = particle_release_times(*args, config)
    variants = {
        "radial": replace(config, pyrolysis_radial_intervals=int(config.pyrolysis_radial_intervals) * 2),
        "time": replace(config, pyrolysis_time_rtol=float(config.pyrolysis_time_rtol) / 10.0),
        "activation": replace(config, pyrolysis_activation_nodes=24),
    }
    rows = []
    for name, other in variants.items():
        refined = particle_release_times(*args, other)
        change = max(
            _relative(base.t_w_s, refined.t_w_s),
            _relative(base.t_v_s, refined.t_v_s),
            _relative(base.tau_pyr_s, refined.tau_pyr_s),
        )
        rows.append({"name": name, "relative_change": change, "gate": 1.0e-3, "met": change < 1.0e-3})
    return {"diameter_m": float(particle["diameter_m"]), "temperature_k": float(particle["temperature_k"]), "tau_pyr_s": float(base.tau_pyr_s), "rows": rows}


def _temperature_change(base: dict, fine: dict) -> float:
    worst = 0.0
    for key in base:
        gap = np.abs(base[key] - fine[key])
        finite = gap[np.isfinite(gap)]
        if finite.size:
            worst = max(worst, float(np.max(finite)))
    return worst


def refine_particle(state, frozen, n_ref: float, m_ref: float, q_ref: float) -> dict:
    """容差收到十分之一，迭代上限加倍，重算残差。ADR 0040 第 6 项。"""
    from hamel.solver.driver import _blocks

    tight = replace(
        frozen,
        particle_tolerance_w=float(frozen.particle_tolerance_w) / 10.0,
        coupling_tolerance_k=float(frozen.coupling_tolerance_k) / 10.0,
        particle_max_iterations=int(frozen.particle_max_iterations) * 2,
        coupling_max_iterations=int(frozen.coupling_max_iterations) * 2,
    )
    change = _temperature_change(particle_temperatures(frozen, state), particle_temperatures(tight, state))
    residual = evaluate_residual(pack_state(state, tight), tight)
    ratios = block_ratios(residual.values, _blocks(tight, n_ref, m_ref, q_ref))
    norm = max(float(value) for value in ratios.values())
    return {
        "temperature_change_k": change,
        "residual_relative_norm": norm,
        "ratios": ratios,
        "temperature_gate_k": 0.1,
        "residual_gate": 1.0e-6,
        "met": change < 0.1 and norm <= 1.0e-6,
    }


def constructed_books() -> dict:
    """构造例：自由空间有气固反应，闭合取 B。数字不是工况 1。"""
    from tests.test_assembly import make_frozen, make_state

    hover_rate = np.zeros((2, 2))
    hover_rate[1, 0] = 0.4
    frozen = make_frozen(
        hover_closure="B", a_hover_per_s=hover_rate,
        ejection_to_free=np.array([0.75, 0.75]), ejection_to_hover=np.array([0.25, 0.25]),
        fuel_up_bed_kg_s=np.full((2, 2), 1.0e-4), fuel_up_free_kg_s=np.full((2, 2), 2.0e-5),
        fuel_down_free_kg_s=np.full((2, 2), 1.0e-5),
    )
    state = make_state(frozen)
    assembly = evaluate_residual(pack_state(state, frozen), frozen)
    present = np.zeros((2, 2), dtype=bool)
    present[0, 1] = True
    theta = np.zeros((2, 2))
    theta[0, 1] = 1.0
    quiet = np.zeros((2, 2))
    held = make_frozen(
        hover_closure="A", hover=present, free_absent=present, theta=theta,
        a_free_per_s=quiet, a_hover_per_s=quiet, b_free_per_s=quiet, b_hover_per_s=quiet,
    )
    held_state = make_state(held)
    held_assembly = evaluate_residual(pack_state(held_state, held), held)
    received = {"C": 0.5, "H": 0.06, "O": 0.2}
    undried_hover = np.zeros((2, 2))
    undried_hover[0, 1] = 1.0e-3
    reported = {
        "phi_s": 0.8, "n_0": 139.5, "cell_division": "2×0.15 m + 2×0.2875 m",
        "reference_return_share": [0.2, 0.8], "reference_elutriation_share": [0.1, 0.05],
        "bed_gas_temperature_deviation_k": 0.0, "moisture_mass_fraction": 0.0,
        "undried_hover_kg": undried_hover, "dried_hover_kg": np.zeros((2, 2)),
        "undried_free_kg": np.zeros((2, 2)), "dried_free_kg": np.zeros((2, 2)),
    }
    stocks = stock_report(held_state, held, held_assembly, reported)
    stocks["note"] = "悬停占比来自闭合 A 的构造。B、C 不设悬停状态。悬停质量含未转化燃料。"
    fuel = {
        "received": received, "dried": dict(received),
        "bed_w": frozen.fuel_up_bed_kg_s[0], "bed_p": np.zeros(2),
        "up_w": frozen.fuel_up_free_kg_s[0], "up_p": np.zeros(2),
        "down_w": frozen.fuel_down_free_kg_s[-1], "down_p": np.zeros(2),
        "redirect_w": np.zeros(2) if assembly.redirect_fuel_by_class_kg_s is None else assembly.redirect_fuel_by_class_kg_s,
        "redirect_p": np.zeros(2),
    }
    measurements = {
        "T_FRB_K": 1100.0, "X_C": 50.0, "CO": 0.1, "CO2": 0.1, "H2": 0.1, "CH4": 0.01, "H2O": 0.1,
    }
    return {
        "state": state, "frozen": frozen,
        "diagnostics": {
            "stocks": stocks,
            "freeboard_closure": freeboard_closure(state, frozen, assembly, fuel),
            "used": used_values(frozen, reported, measurements),
            "bed_gas_temperature_deviation_k": 0.0,
        },
        "measurements": measurements,
    }


def closed_form_from_state(state, frozen, reported, config: SolverConfig) -> dict:
    """逐粒级适用性。气体取自由空间最下一格，初速取床面抛出分布的中位数。"""
    flow = state.freeboard_gas_mol_s[-1]
    total = float(np.sum(flow))
    temperature = float(state.freeboard_temperature_k[-1])
    pressure = float(frozen.freeboard_pressure_pa[-1])
    fractions = {name: float(flow[index]) / total for index, name in enumerate(GAS_SPECIES)}
    rho = density_kg_m3(temperature, pressure, fractions)
    mu = mixture_viscosity_pa_s(fractions, temperature, pressure)
    area = math.pi * (0.5 * float(frozen.bed_diameter_m)) ** 2
    speed = volume_flow_m3_s(total, temperature, pressure) / area
    height = float(frozen.n_freeboard) * float(frozen.freeboard_height_m)
    bubble = float(reported["trajectory"]["bubble_m_s"])
    phi = float(reported.get("phi_s") or reported["trajectory"]["phi_s"])
    limit = float(config.freeboard_time_limit_s)
    rows = [
        applicability_row(
            float(diameter), rho, mu, float(frozen.rho_p_kg_m3), phi, speed, bubble, height, limit,
        )
        for diameter in frozen.diameter_m
    ]
    for row in rows:
        row["rho_g"] = rho
        row["gas_m_s"] = speed
    return {"benchmark": fixed_closed_form(int(config.freeboard_ejection_samples)), "classes": rows}


def refinement_from_case(config: SolverConfig, reported, state, frozen, refs) -> dict:
    """用这一次运行保存的轨迹和热解输入做加密。粒级逐个算，取最大相对变化。"""
    column = dict(reported["trajectory"])
    diameters = [float(item) for item in np.asarray(column["diameter_m"], dtype=float)]
    worst = 0.0
    existence = False
    rows = []
    for diameter in diameters:
        one = dict(column)
        one["diameter_m"] = diameter
        result = refine_trajectory(config, one)
        worst = max(worst, result["relative_change"])
        existence = existence or result["existence_changed"]
        rows.append({
            "diameter_m": diameter, "relative_change": result["relative_change"], "a_f": result["a_f"],
            "existence_changed": result["existence_changed"], "met": result["met"],
        })
    pyrolysis = dict(reported["pyrolysis"])
    pyro_rows = []
    worst_pyro = 0.0
    for diameter in [float(item) for item in np.asarray(pyrolysis.pop("diameter_m"), dtype=float)]:
        particle = dict(pyrolysis)
        particle["diameter_m"] = diameter
        result = refine_pyrolysis(config, particle)
        change = max(item["relative_change"] for item in result["rows"])
        worst_pyro = max(worst_pyro, change)
        pyro_rows.append({"diameter_m": diameter, "temperature_k": result["temperature_k"], "tau_pyr_s": result["tau_pyr_s"], "rows": result["rows"]})
    report = {
        "trajectory": {
            "rows": rows, "relative_change": worst, "existence_changed": existence,
            "gate": 1.0e-4, "met": (not existence) and worst < 1.0e-4,
        },
        "pyrolysis": {"rows": pyro_rows, "relative_change": worst_pyro, "gate": 1.0e-3, "met": worst_pyro < 1.0e-3},
    }
    if state is not None and frozen is not None and refs is not None:
        report["particle"] = refine_particle(state, frozen, *refs)
    return report


def refinement_report(config: SolverConfig | None = None, trajectory=None, pyrolysis=None, particle=None) -> dict:
    """加密入口。不传柱和颗粒时用玩具输入；传入的是这次运行自己的输入。"""
    config = config or SolverConfig()
    report = {
        "trajectory": refine_trajectory(config, toy_trajectory() if trajectory is None else trajectory),
        "pyrolysis": refine_pyrolysis(config, toy_particle() if pyrolysis is None else pyrolysis),
    }
    if particle is not None:
        report["particle"] = refine_particle(*particle)
    return report
