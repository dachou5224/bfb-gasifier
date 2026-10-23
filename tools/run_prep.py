"""首次运行前的报告记账、闭式比对和数值加密。

不进入残差，不改冻结值，不设通过门。依据 ADR 0040 第 5、6 项。
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError, account_freeboard, closed_form_state, discriminant, ejection_velocity_samples,
    integrate_motion, motion_coefficients, terminal_settling_velocity,
)
from hamel.properties.data import GAS_SPECIES, atomic_weight, atoms, molar_mass_kg_mol
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


def _of(kind: str, kilograms, fractions=None) -> dict:
    mass = float(np.sum(kilograms))
    packet = _blank()
    packet["mass"] = mass
    if kind == "char":
        packet["C"] = mass
    elif kind == "fuel":
        for element, share in fractions.items():
            packet[element] = mass * float(share)
    return packet


def _mobile_ash(state, frozen, cell) -> np.ndarray:
    ash = np.array(state.freeboard_ash_free_kg[cell], dtype=float, copy=True)
    char = state.freeboard_char_free_kg[cell]
    for klass in range(1, frozen.n_class):
        if char[klass] <= 0.0:
            ash[klass] = 0.0
    return ash


def _combine(*parts) -> dict:
    total = _blank()
    for part in parts:
        _add(total, part)
    return total


def freeboard_closure(state, frozen, assembly, fractions) -> dict:
    """自由空间全部相的边界流。反应不另计。B、C 的改送是出流。

    F^ej,t 在抛出入流里计一次，在改送出流里再计一次。ADR 0040 第 5 项。
    """
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    share = _bed_ash_by_class(state.bed_char_kg, inert)
    thrown_char = frozen.k_auf_per_s[0] * state.bed_char_kg[0]
    thrown_ash = frozen.k_auf_per_s[0] * share[0]
    gas_in = _gas(state.bed_gas_mol_s[0, 0] + state.bed_gas_mol_s[0, 1])
    side_in = _gas(frozen.secondary_feed_mol_s)
    solids_in = _combine(_of("char", thrown_char), _of("ash", thrown_ash), _of("fuel", frozen.fuel_up_bed_kg_s[0], fractions))
    gas_out = _gas(state.freeboard_gas_mol_s[0])
    up = frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]
    down = frozen.b_free_per_s[-1] + frozen.b_hover_per_s[-1]
    carried = _combine(
        _of("char", up * state.freeboard_char_free_kg[0]),
        _of("ash", up * _mobile_ash(state, frozen, 0)),
        _of("fuel", frozen.fuel_up_free_kg_s[0], fractions),
    )
    returned = _combine(
        _of("char", down * state.freeboard_char_free_kg[-1]),
        _of("ash", down * _mobile_ash(state, frozen, -1)),
        _of("fuel", frozen.fuel_down_free_kg_s[-1], fractions),
    )
    redirected = _combine(
        _of("char", assembly.redirect_char_kg_s),
        _of("ash", assembly.redirect_ash_kg_s),
        _of("fuel", assembly.redirect_fuel_kg_s, fractions),
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
    hover_char = float(np.sum(thrown_char * frozen.ejection_to_hover))
    return {
        "streams": streams,
        "totals": totals,
        "reaction_char_kg_s": float(np.sum(assembly.char_consumption_free_kg_s) + np.sum(assembly.char_consumption_hover_kg_s)),
        "ejection_hover_char_kg_s": hover_char,
        "redirect_char_kg_s": float(assembly.redirect_char_kg_s),
    }


def stock_report(state, frozen, assembly, reported) -> dict:
    """收敛存量，以及存量乘输运系数得到的回落和带出。"""
    classes = []
    for klass in range(frozen.n_class):
        upward = float(frozen.a_free_per_s[0, klass] + frozen.a_hover_per_s[0, klass])
        downward = float(frozen.b_free_per_s[-1, klass] + frozen.b_hover_per_s[-1, klass])
        top = float(state.freeboard_char_free_kg[0, klass] + _mobile_ash(state, frozen, 0)[klass])
        bottom = float(state.freeboard_char_free_kg[-1, klass] + _mobile_ash(state, frozen, -1)[klass])
        classes.append({
            "class": klass + 1,
            "reference_return_share": _listed(reported, "reference_return_share", klass),
            "reference_elutriation_share": _listed(reported, "reference_elutriation_share", klass),
            "return_kg_s": downward * bottom + float(frozen.fuel_down_free_kg_s[-1, klass]),
            "elutriation_kg_s": upward * top + float(frozen.fuel_up_free_kg_s[0, klass]),
        })
    cells = []
    hover = []
    for cell in range(frozen.n_freeboard):
        holds = []
        for klass in range(frozen.n_class):
            free = float(state.freeboard_char_free_kg[cell, klass] + state.freeboard_ash_free_kg[cell, klass])
            held = float(state.freeboard_char_hover_kg[cell, klass] + state.freeboard_ash_hover_kg[cell, klass])
            holds.append((free + held, held))
        solid = float(sum(item[0] for item in holds))
        consumed = float(assembly.char_consumption_free_kg_s[cell].sum() + assembly.char_consumption_hover_kg_s[cell].sum())
        cells.append({"cell": cell, "solid_kg": solid, "char_consumption_kg_s": consumed})
        for klass, (total, held) in enumerate(holds):
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


def attach_converged(outcome, state, frozen, assembly, physical, reported) -> None:
    """收敛之后把记账放进诊断。求解已经结束。"""
    fractions = fuel_element_fractions(physical)
    outcome.diagnostics["stocks"] = stock_report(state, frozen, assembly, reported)
    outcome.diagnostics["freeboard_closure"] = freeboard_closure(state, frozen, assembly, fractions)
    outcome.diagnostics["used"] = used_values(frozen, reported, outcome.measurements)
    outcome.diagnostics["bed_gas_temperature_deviation_k"] = (reported or {}).get("bed_gas_temperature_deviation_k")


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def fixed_closed_form(step_s: float, samples: int) -> dict:
    """常气速、Δ > 0 的固定构造。步长用这次运行的值。样本数不改变这一条轨迹。"""
    item = _BENCHMARK
    gas = item["gas_m_s"]
    coefficient = motion_coefficients(gas, True, item["rho_g"], item["mu_g"], item["diameter_m"], item["rho_p"], item["phi_s"])
    delta = discriminant(*coefficient)
    analytic_u, analytic_h = closed_form_state(item["duration_s"], *coefficient, item["initial_m_s"], item["height_m"])
    numeric_h, numeric_u = integrate_motion(
        item["height_m"], item["initial_m_s"], item["duration_s"], gas, item["rho_g"], item["mu_g"],
        item["diameter_m"], item["rho_p"], item["phi_s"], step_s,
    )
    velocity = _relative(numeric_u, analytic_u)
    height = _relative(numeric_h, analytic_h)
    return {
        "delta": delta, "step_s": float(step_s), "samples": int(samples),
        "velocity_relative": velocity, "height_relative": height, "max_relative": max(velocity, height),
    }


def applicability_row(diameter_m, rho_g, mu_g, rho_p, phi_s, gas_m_s, bubble_m_s, height_m, step_s, time_limit_s) -> dict:
    """一个粒级。Δ ≤ 0 记 N/A。Δ > 0 时比到分支、实数域、出口、床面或时间上限。"""
    speeds, _shares, direct = ejection_velocity_samples(float(bubble_m_s), 1)
    initial = float(speeds[0])
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
            row["status"] = "N/A"
            row["stopped"] = "scale"
        return row
    height = 0.0
    velocity = initial
    elapsed = 0.0
    worst_u = 0.0
    worst_h = 0.0
    stopped = "time_limit"
    while elapsed < time_limit_s:
        piece = min(float(step_s), time_limit_s - elapsed)
        new_h, new_u = integrate_motion(height, velocity, piece, gas_m_s, rho_g, mu_g, diameter_m, rho_p, phi_s, piece)
        elapsed += piece
        crossed = (initial - gas_m_s) * (new_u - gas_m_s) <= 0.0
        try:
            closed_u, closed_h = closed_form_state(elapsed, *coefficient, initial, 0.0)
        except FreeboardError:
            stopped = "real_domain"
            break
        if not crossed:
            worst_u = max(worst_u, abs(new_u - closed_u) / abs(gas_m_s))
            worst_h = max(worst_h, abs(new_h - closed_h) / height_m)
            height, velocity = new_h, new_u
        if crossed:
            stopped = "branch_change"
            break
        if height >= height_m:
            stopped = "exit"
            break
        if height <= 0.0:
            stopped = "bed"
            break
    row["stopped"] = stopped
    row["max_velocity_relative"] = worst_u
    row["max_height_relative"] = worst_h
    return row


def toy_applicability_rows() -> list:
    """玩具用的两个粒级。一个预期 Δ > 0，一个用更小的颗粒去碰 Δ ≤ 0。"""
    return [
        (8.0e-3, 1.1, 1.8e-5, 900.0, 0.8, 1.5, 1.0, 0.4),
        (1.5e-4, 1.2, 2.0e-5, 2000.0, 1.0, 0.5, 1.0, 0.5),
    ]


def closed_form_report(config: SolverConfig, rows: list) -> dict:
    step = float(config.freeboard_max_step_s)
    limit = float(config.freeboard_time_limit_s)
    table = [
        applicability_row(*item, step, limit) if not isinstance(item, dict) else item
        for item in rows
    ]
    return {"benchmark": fixed_closed_form(step, int(config.freeboard_ejection_samples)), "classes": table}


def _books(samples: int, step_s: float):
    diameter = 2.0e-4
    rho_g, mu, rho_p, phi = 1.2, 2.0e-5, 1200.0, 0.8
    settling = terminal_settling_velocity(rho_g, mu, diameter, rho_p, phi)
    edges = np.array([0.0, 0.25, 0.5])
    def gas(_cell, _height, speed=settling + 3.0):
        return speed
    speeds, shares, direct = ejection_velocity_samples(1.0, samples)
    config = replace(SolverConfig(), freeboard_ejection_samples=samples, freeboard_max_step_s=step_s)
    return account_freeboard(edges, gas, np.full(2, rho_g), np.full(2, mu), diameter, rho_p, phi, speeds, shares, direct, config)


def _book_change(base, fine) -> float:
    worst = 0.0
    for left, right in zip(base.cells, fine.cells):
        for name in ("a_f", "a_t", "b_f", "b_t"):
            first, second = getattr(left, name), getattr(right, name)
            if first is None or second is None:
                raise RuntimeError("加密改变了系数是否存在")
            worst = max(worst, _relative(float(first), float(second)))
    return worst


def refine_trajectory(config: SolverConfig) -> dict:
    """样本数加倍且步长减半。门是 a、b 的相对变化 1e-4。ADR 0016、ADR 0028。"""
    step = float(config.freeboard_max_step_s)
    samples = int(config.freeboard_ejection_samples)
    change = _book_change(_books(samples, step), _books(samples * 2, step / 2.0))
    return {"samples": samples, "step_s": step, "relative_change": change, "gate": 1.0e-4, "met": change < 1.0e-4}


def refine_pyrolysis(config: SolverConfig) -> dict:
    """径向加倍、时间容差收到十分之一、活化能节点 16 增到 24。门 1e-3。ADR 0027。"""
    material = lignite_particle_properties()
    gas = {"N2": 0.79, "O2": 0.21}
    args = (2.0e-4, 1100.0, 1.0e5, 0.02, 298.15, gas, 0.2, 0.4, material)
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
    return {"diameter_m": 2.0e-4, "rows": rows}


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
    reported = {
        "phi_s": 0.8, "n_0": 139.5, "cell_division": "2×0.15 m + 2×0.2875 m",
        "reference_return_share": [0.2, 0.8], "reference_elutriation_share": [0.1, 0.05],
        "bed_gas_temperature_deviation_k": 0.0,
    }
    stocks = stock_report(held_state, held, held_assembly, reported)
    stocks["note"] = "悬停占比来自闭合 A 的构造。B、C 不设悬停状态。"
    measurements = {
        "T_FRB_K": 1100.0, "X_C": 50.0, "CO": 0.1, "CO2": 0.1, "H2": 0.1, "CH4": 0.01, "H2O": 0.1,
    }
    return {
        "state": state, "frozen": frozen,
        "diagnostics": {
            "stocks": stocks,
            "freeboard_closure": freeboard_closure(state, frozen, assembly, {"C": 0.5, "H": 0.06, "O": 0.2}),
            "used": used_values(frozen, reported, measurements),
            "bed_gas_temperature_deviation_k": 0.0,
        },
        "measurements": measurements,
    }


def refinement_report(config: SolverConfig | None = None, particle=None) -> dict:
    """玩具入口。逐根 δ 和床层子步由该次运行带入，这里不重算工况 1。"""
    config = config or SolverConfig()
    report = {"trajectory": refine_trajectory(config), "pyrolysis": refine_pyrolysis(config)}
    if particle is not None:
        report["particle"] = refine_particle(*particle)
    return report
