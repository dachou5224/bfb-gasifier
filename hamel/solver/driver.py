"""外层对齐与工况求解。每轮重算预计算，内层从上一轮单元解出发。
求解手段：ADR 0012；ADR 0017 第 10、14 项；ADR 0035 第 2 项；ADR 0045。
"""
from __future__ import annotations

import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.cases.balance import ElementStream, external_balance
from hamel.properties import (
    ash_specific_heat_j_per_kg_k, char_specific_enthalpy_j_per_kg,
    enthalpy_j_mol, gas_constant_j_mol_k, molar_mass_kg_mol,
)
from hamel.properties.data import GAS_SPECIES, atomic_weight, atoms
from hamel.solver.frozen import FreeboardDomainError, TrialNotEvaluable
from hamel.solver.sensitivity import bubble_balance_ratio, redirect_feedback
from hamel.solver.gates import (
    OuterMonitor, ResidualBlocks, arrays_allowed, balance_met, block_ratios, inner_met, outer_met,
)
from hamel.solver.initial import DomainRejected, InitialStateError, construct
from hamel.solver.newton import SolveTimeout, damped_newton, normalized_merit
from hamel.solver.precalc import PrecalcError, precalculate
from hamel.solver.residual import evaluate_residual, wall_area_m2
from hamel.solver.vector import jacobian_layout, pack_state, unknown_slots, unpack_state

# ADR 0029。不是新的阈值。
_AMBIENT_K = 298.15
_WALL_RELATIVE = 1.0e-4
_DRY_GAS = ("CO", "CO2", "CH4", "C2H4", "C6H6", "C10H8")


@dataclass
class SolveResult:
    status: str
    reason: str
    rounds: int = 0
    iterations: int = 0
    ratios: dict = field(default_factory=dict)
    balance: dict | None = None
    nonnegative: bool | None = None
    measurements: dict | None = None
    diagnostics: dict = field(default_factory=dict)
    x: np.ndarray | None = None


def solve_aligned(refresh, evaluate, prepare, monitor_from, consistent, config, deadline=None, journal=None) -> SolveResult:
    """refresh(x, round, consistency) 给出本轮起点。x 为 None 时是第一次。
    (a) 比较相邻两轮结束时的监测量。(b) 用最终单元解重算预计算，再查内层门。
    """
    diagnostics = {"rounds": []}
    x = None
    previous = None
    refs = None
    iterations = 0
    for round_index in range(1, int(config.outer_round_limit) + 1):
        _clock(deadline)
        try:
            opened = refresh(x, round_index, False)
        except PrecalcError as exc:
            return SolveResult("DOMAIN_STOP", str(exc), round_index, iterations, diagnostics=diagnostics, x=x)
        blocks, scales, feasible, refs = prepare(opened)
        mark = 0 if journal is None else len(journal.setdefault("iterations", []))
        # 核对绑定首次实际构造的雅可比，不绑定外层轮号：0 次内层迭代不消耗这次核对。ADR 0046 第 5 项。
        layout = None if journal is None else journal.get("jacobian_layout")
        verify = layout is not None and not journal.get("jacobian_verified")
        try:
            inner = damped_newton(
                opened.x, lambda trial, opened=opened: evaluate(trial, opened), scales,
                lambda residual, blocks=blocks: inner_met(residual, blocks, config.residual_tolerance),
                feasible, lambda residual, blocks=blocks: normalized_merit(residual, blocks),
                config, deadline, journal, blocks, layout, verify,
            )
        except SolveTimeout:
            raise
        except PrecalcError as exc:
            return SolveResult("DOMAIN_STOP", str(exc), round_index, iterations, diagnostics=diagnostics, x=x)
        iterations += inner.iterations
        row = {"round": round_index, "inner": inner.status, "reason": inner.reason, "info": opened.info}
        if journal is not None:
            row["newton"] = journal["iterations"][mark:]
            if journal.get("failure_npz"):
                diagnostics["failure_npz"] = journal["failure_npz"]
        diagnostics["rounds"].append(row)
        if inner.status != "CONVERGED":
            return SolveResult(inner.status, inner.reason, round_index, iterations, diagnostics=diagnostics, x=inner.x)
        try:
            monitor = monitor_from(opened, inner.x)
        except TrialNotEvaluable as exc:
            return SolveResult("NOT_CONVERGED", str(exc), round_index, iterations, diagnostics=diagnostics, x=inner.x)
        aligned = previous is not None and outer_met(
            monitor, previous, m_dot_ref=refs[1], n_dot_ref=refs[0], m_tot=refs[3],
            temperature_tolerance=config.outer_temperature_tolerance_k, relative_tolerance=config.outer_relative_tolerance,
        )
        previous = monitor
        x = inner.x
        if not aligned:
            continue
        try:
            checked = consistent(inner.x, round_index)
        except PrecalcError as exc:
            return SolveResult("DOMAIN_STOP", str(exc), round_index, iterations, diagnostics=diagnostics, x=inner.x)
        ratios = block_ratios(checked.residual, checked.blocks)
        if inner_met(checked.residual, checked.blocks, config.residual_tolerance):
            diagnostics["consistent"] = checked.info
            return SolveResult("CONVERGED", "", round_index, iterations, ratios, diagnostics=diagnostics, x=inner.x)
    return SolveResult("NOT_CONVERGED", "外层超过 outer_round_limit", int(config.outer_round_limit), iterations, diagnostics=diagnostics, x=x)


def solve_case(path, case_id: str = "1", *, adiabatic: bool = False, deadline=None, config: SolverConfig | None = None, artifact_dir=None) -> SolveResult:
    """第 1 关入口。绝热对照传 adiabatic=True：k_W = 0，初值不扣散热。

    config 是 SolverConfig 副本。缺省时用工况装配的默认配置。
    """
    assembled = assemble_case(path, case_id)
    config = assembled.solver if config is None else config
    physical = assembled.physical
    n_ref, m_ref, q_ref = _references(physical)
    holder = {"state": None, "reaction": None, "drain": 0.0, "fuel_return": 0.0, "records": [], "redirect": (0.0, 0.0, 0.0, 0.0)}
    journal = {"iterations": []}
    if artifact_dir is not None:
        folder = Path(artifact_dir)
        folder.mkdir(parents=True, exist_ok=True)
        journal["failure_path"] = str(folder / "failure.npz")

    def _feedback(state):
        """从已解固相与它对应的冻结量更新床区反应汇与改送反馈。
        求解手段：ADR 0052 第 1 项第 5 步；PR #130 审查 P2。
        """
        if not holder["records"]:
            return
        reaction, redirect = solid_feedback(state, holder["records"][-1].frozen)
        holder["reaction"] = reaction
        holder["redirect"] = redirect

    def _precalc(iteration, consistency, state=None):
        used = holder["state"] if state is None else state
        if state is not None:
            _feedback(state)
        bed_return, ash_out, fuel_out = redirect_feedback(config.hover_closure, holder["redirect"])
        result = precalculate(
            path, case_id, used, reaction_solid_kg_s=holder["reaction"],
            previous_drain_kg_s=holder["drain"], previous_fuel_return_kg_s=holder["fuel_return"],
            adiabatic=adiabatic, iteration=iteration, consistency=consistency, solver=config,
            redirected_return_kg_s=bed_return, redirected_ash_out_kg_s=ash_out, redirected_fuel_out_kg_s=fuel_out,
        )
        holder["drain"] = float(result.frozen.drain_ash_kg_s + result.frozen.drain_char_kg_s.sum())
        holder["fuel_return"] = float(result.frozen.fuel_down_free_kg_s[-1].sum())
        holder["records"].append(result)
        return result

    def _initial(iteration):
        """初值：焦油另加，固相、排灰与床区向下流按同一状态联立牛顿，再与预计算交替。求解手段：ADR 0052"""
        first = _precalc(iteration, False)
        loops: list = []
        # 先挂上空列表：交替中途停止时，已完成的记录仍写进报告。
        holder["initial"] = {"loops": loops, "note": "停止"}
        try:
            result, loops, note = construct(
                first, lambda state: _precalc(iteration, False, state), m_ref,
                gate=float(config.residual_tolerance), config=config, records=loops,
            )
        except DomainRejected as exc:
            holder["initial"] = {"loops": loops, "note": f"停止：{exc.phase}的有效域拒绝"}
            raise PrecalcError(f"初值构造停止：{exc}") from exc
        except FreeboardDomainError as exc:
            holder["initial"] = {"loops": loops, "note": "停止：有效域拒绝"}
            raise PrecalcError(f"初值构造停止：有效域拒绝：{exc}") from exc
        except InitialStateError as exc:
            holder["initial"] = {"loops": loops, "note": "停止"}
            raise PrecalcError(f"初值构造停止：{exc}") from exc
        holder["initial"] = {"loops": loops, "note": note}
        return result

    def refresh(x, iteration, consistency):
        del x
        if holder["state"] is None and not consistency:
            result = _initial(iteration)
        else:
            result = _precalc(iteration, consistency)
        return _Opened(pack_state(result.state, result.frozen), _info(result))

    def evaluate(trial, opened):
        return evaluate_residual(trial, holder["records"][-1].frozen).values

    def prepare(opened):
        frozen = holder["records"][-1].frozen
        blocks = _blocks(frozen, n_ref, m_ref, q_ref)
        m_tot = float(frozen.geometric_mass_kg.sum())
        slots = unknown_slots(frozen)
        journal["blocks"] = blocks
        journal["jacobian_layout"] = jacobian_layout(slots, frozen)
        journal["names"] = [slot.name for slot in slots]
        journal["explain"] = lambda trial, frozen=frozen: _infeasible_clause(trial, frozen)
        bottom = opened.info.get("bottom_residual_kg_s")
        journal["bottom_norm"] = None if bottom is None else abs(float(bottom)) / m_ref
        return blocks, _scales(frozen, n_ref, m_tot), lambda trial, frozen=frozen: _feasible(trial, frozen), (n_ref, m_ref, q_ref, m_tot)

    def monitor_from(opened, values):
        frozen = holder["records"][-1].frozen
        state = unpack_state(values, frozen)
        holder["state"] = state
        assembly = evaluate_residual(values, frozen)
        holder["reaction"] = -assembly.char_consumption_bed_kg_s.sum(axis=1)
        holder["redirect"] = _redirect_tuple(assembly)
        return _monitor(state, frozen)

    def consistent(values, iteration):
        state = unpack_state(values, holder["records"][-1].frozen)
        holder["state"] = state
        result = _precalc(iteration, True)
        residual = evaluate_residual(pack_state(state, result.frozen), result.frozen)
        holder["redirect"] = _redirect_tuple(residual)
        holder["assembly"] = residual
        holder["final"] = result
        holder["state"] = state
        return _Checked(residual.values, _blocks(result.frozen, n_ref, m_ref, q_ref), _info(result))

    outcome = solve_aligned(refresh, evaluate, prepare, monitor_from, consistent, config, deadline, journal)
    outcome.diagnostics["initial"] = holder.get("initial")
    if artifact_dir is not None and outcome.status == "DAMPING_FAILED" and holder["records"]:
        frozen_path = Path(artifact_dir) / "frozen.pkl"
        with frozen_path.open("wb") as handle:
            pickle.dump(holder["records"][-1].frozen, handle, protocol=pickle.HIGHEST_PROTOCOL)
        outcome.diagnostics["frozen_path"] = str(frozen_path)
    outcome.diagnostics["geldart"] = [item.geldart_record for item in holder["records"]]
    outcome.diagnostics["adiabatic"] = adiabatic
    if holder["records"]:
        # ADR 0050 第 5 项：本次预计算的最大段数、归宿时刻最大值与 99 % 分位、最长的 (粒级, d, φ_S, 起点)。
        outcome.diagnostics["trajectory"] = holder["records"][-1].reported.get("trajectory")
        # ADR 0051 第 2 项：第 1 次预计算的 δ、μ 与逐 φ_p 的去向，只报告。
        outcome.diagnostics["hover_drift"] = holder["records"][0].reported.get("hover_drift")
    if outcome.status != "CONVERGED":
        return outcome
    result = holder["final"]
    state = holder["state"]
    heat = _wall_heat(state, result.frozen)
    gap_in, gap_out, inlet_mol = _streams(state, result.frozen, physical, holder.get("redirect"))
    gap = external_balance(gap_in, gap_out)
    outcome.balance = balance_met(gap, inlet_mol, heat, q_ref, config.element_balance_tolerance, config.energy_balance_tolerance)
    outcome.balance["heat_w"] = heat
    outcome.balance["k_w"] = float(result.frozen.k_w_w_m2_k)
    target = 0.0 if adiabatic else float(_case_heat(path, case_id))
    outcome.balance["heat_target_w"] = target
    outcome.balance["heat_relative"] = None if target <= 0.0 else abs(heat - target) / target
    bottom = abs(float(result.bottom_residual_kg_s)) / m_ref
    outcome.nonnegative = _nonnegative(state, result.frozen) and bottom <= config.residual_tolerance and np.all(result.frozen.m_dot_ab_kg_s >= 0.0)
    outcome.diagnostics["bottom_relative"] = bottom
    outcome.diagnostics["final"] = _with_sensitivity(_info(result), state, result.frozen, holder["assembly"])
    heights = np.full(state.freeboard_temperature_k.shape, float(result.frozen.freeboard_height_m))
    outcome.measurements = _measurements(state, physical, heights)
    heat_ok = adiabatic or (outcome.balance["heat_relative"] is not None and outcome.balance["heat_relative"] <= _WALL_RELATIVE)
    if not outcome.balance["ok"] or not outcome.nonnegative or not heat_ok:
        outcome.status = "NOT_CONVERGED"
        outcome.reason = "外部平衡、非负或壁面散热未通过"
        return outcome
    from tools.run_prep import attach_converged
    attach_converged(
        outcome, state, result.frozen, holder["assembly"], physical, result.reported,
        config, (n_ref, m_ref, q_ref),
    )
    return outcome


@dataclass
class _Opened:
    x: np.ndarray
    info: dict


@dataclass
class _Checked:
    residual: np.ndarray
    blocks: ResidualBlocks
    info: dict


def solid_feedback(state, frozen):
    """从已解固相与冻结量取床区反应汇与改送流，供下一次预计算回馈。

    求解手段：ADR 0052 第 1 项第 5 步；PR #130 审查 P2。
    """
    assembly = evaluate_residual(pack_state(state, frozen), frozen)
    return -assembly.char_consumption_bed_kg_s.sum(axis=1), _redirect_tuple(assembly)


def _references(physical):
    dry = float(physical.fuel_mass_flow) * (1.0 - float(physical.moisture_mass_fraction))
    steam = float(physical.steam_mass_flow) / molar_mass_kg_mol("H2O")
    gas = float(physical.primary_air_amount_flow) + float(physical.secondary_air_amount_flow)
    gas += float(physical.nitrogen_purge_volume_flow) + steam
    return gas, dry, float(physical.heating_value) * float(physical.fuel_mass_flow)


def _blocks(frozen, n_ref, m_ref, q_ref) -> ResidualBlocks:
    gas, solid, energy = [], [], []
    for slot in unknown_slots(frozen):
        (gas if slot.kind == "gas_outflow" else energy if slot.kind == "temperature" else solid).append(slot.index)
    return ResidualBlocks(np.array(gas), np.array(solid), np.array(energy), n_ref, m_ref, q_ref)


def _scales(frozen, n_ref, m_tot) -> np.ndarray:
    """固体尺度取 M_tot。解释性选择：ADR 0037。"""
    scale = np.empty(len(unknown_slots(frozen)))
    for slot in unknown_slots(frozen):
        if slot.kind == "gas_outflow":
            scale[slot.index] = n_ref
        elif slot.kind == "temperature":
            scale[slot.index] = 1.0
        else:
            scale[slot.index] = m_tot
    return scale


def _feasible(values, frozen) -> bool:
    state = unpack_state(values, frozen)
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    return arrays_allowed(
        (state.bed_gas_mol_s, state.freeboard_gas_mol_s),
        (state.bed_char_kg, state.freeboard_char_free_kg, state.freeboard_char_hover_kg, state.freeboard_ash_free_kg, state.freeboard_ash_hover_kg),
        np.concatenate([state.bed_temperature_k, state.freeboard_temperature_k]),
        inert,
    )


def _infeasible_clause(values, frozen) -> str:
    """与 arrays_allowed 同一组条件，按 ADR 0017 第 7 项的顺序报告第一条。"""
    state = unpack_state(values, frozen)
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    groups = (
        ("Ṅ < 0", state.bed_gas_mol_s), ("Ṅ < 0", state.freeboard_gas_mol_s),
        ("m < 0", state.bed_char_kg), ("m < 0", state.freeboard_char_free_kg),
        ("m < 0", state.freeboard_char_hover_kg), ("m < 0", state.freeboard_ash_free_kg),
        ("m < 0", state.freeboard_ash_hover_kg), ("m_inert < 0", inert),
    )
    for label, array in groups:
        hit = _bound_hit(array, label)
        if hit:
            return hit
    temperature = np.concatenate([state.bed_temperature_k, state.freeboard_temperature_k])
    if not np.isfinite(temperature).all():
        return "T 非有限"
    if np.any(temperature <= 0.0):
        return f"T ≤ 0（下标 {int(np.flatnonzero(temperature <= 0.0)[0])}）"
    return "feasible 不满足"


def _bound_hit(array, label):
    flat = np.asarray(array, dtype=float).ravel()
    if flat.size == 0 or (np.isfinite(flat).all() and np.all(flat >= 0.0)):
        return None
    if not np.isfinite(flat).all():
        return f"{label} 的量非有限"
    return f"{label}（下标 {int(np.flatnonzero(flat < 0.0)[0])}）"


def _monitor(state, frozen) -> OuterMonitor:
    temperature = np.concatenate([state.freeboard_temperature_k, state.bed_temperature_k])
    elutriation = (frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * (
        state.freeboard_char_free_kg[0] + state.freeboard_ash_free_kg[0]
    )
    solid = np.concatenate([
        frozen.k_auf_per_s * frozen.geometric_mass_kg, frozen.m_dot_ab_kg_s, elutriation,
        np.zeros(1), frozen.drain_char_kg_s, np.array([frozen.drain_ash_kg_s]),
    ])
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    gas = np.concatenate([state.bed_gas_mol_s.ravel(), state.freeboard_gas_mol_s.ravel(), _exchange(state, frozen).ravel()])
    return OuterMonitor(temperature, solid, np.concatenate([frozen.geometric_mass_kg, inert]), gas)


def _exchange(state, frozen):
    """与残差中的交换流同一式。ADR 0011 第 3 条。"""
    gas_r = gas_constant_j_mol_k()
    exchange = np.zeros_like(state.bed_gas_mol_s[:, 0])
    for cell in range(frozen.n_bed):
        exchange[cell] = frozen.k_bd_per_s[cell] * frozen.v_bubble_exchange_m3[cell] * (
            _concentration(state.bed_gas_mol_s[cell, 0], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell], gas_r)
            - _concentration(state.bed_gas_mol_s[cell, 1], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell], gas_r)
        )
    return exchange


def _concentration(flow, temperature, pressure, gas_r):
    total = float(np.sum(flow))
    if total <= 0.0:
        raise TrialNotEvaluable("气相总摩尔流不是正数，交换流没有定义")
    return flow / total * pressure / (gas_r * temperature)


def _wall_heat(state, frozen) -> float:
    total = 0.0
    pairs = (("bed", state.bed_temperature_k, frozen.bed_height_m), ("freeboard", state.freeboard_temperature_k, frozen.freeboard_height_m))
    for zone, temperatures, height in pairs:
        for index, temperature in enumerate(temperatures):
            area = wall_area_m2(zone, index, len(temperatures), height, frozen.bed_diameter_m, frozen.freeboard_diameter_m)
            total += frozen.k_w_w_m2_k * area * (float(temperature) - _AMBIENT_K)
    return float(total)


def _nonnegative(state, frozen) -> bool:
    return _feasible(pack_state(state, frozen), frozen) and bool(np.all(np.isfinite(pack_state(state, frozen))))


def _measurements(state, physical, freeboard_dh_m) -> dict:
    """验收量。解释性选择：ADR 0040。

    T̄_FRB = Σ T_i·Δh_i / Σ Δh_i，对自由空间全部 cell（含锥段）按高度加权。
    等高时公因子约去，用算术平均，与现行求值逐位相同。
    X_C = 自由空间出口气体（含焦油蒸气）中的碳 / 加入燃料中的碳。
    分母只计燃料碳；分子不计带出焦炭与排灰中的碳。
    """
    exit_flow = state.freeboard_gas_mol_s[0]
    total = float(np.sum(exit_flow))
    named = {name: float(exit_flow[GAS_SPECIES.index(name)]) / total for name in ("CO", "CO2", "H2", "CH4", "H2O")}
    carbon_out = sum(float(exit_flow[index]) * atoms(name).get("C", 0) for index, name in enumerate(GAS_SPECIES) if name in _DRY_GAS or atoms(name).get("C", 0))
    carbon_in = _fuel_moles(physical, float(physical.fuel_mass_flow))[0]
    temperatures = np.asarray(state.freeboard_temperature_k, dtype=float)
    heights = np.asarray(freeboard_dh_m, dtype=float)
    if heights.shape != temperatures.shape or heights.size < 1 or np.any(heights <= 0.0):
        raise TrialNotEvaluable("自由空间格子高度必须与气体温度一一对应，且为正")
    if np.all(heights == heights[0]):
        named["T_FRB_K"] = float(np.mean(temperatures))
    else:
        named["T_FRB_K"] = float(np.dot(temperatures, heights) / np.sum(heights))
    named["X_C"] = 100.0 * carbon_out / carbon_in
    return named


def _redirect_tuple(assembly) -> tuple:
    return (
        float(assembly.redirect_char_kg_s), float(assembly.redirect_ash_kg_s),
        float(assembly.redirect_fuel_kg_s), float(assembly.redirect_enthalpy_w),
    )
def _streams(state, frozen, physical, redirect=None):
    inlet = {name: 0.0 for name in ("C", "H", "O")}
    fuel_mass = float(frozen.fuel_feed_kg_s.sum())
    per_kg = _fuel_moles(physical, 1.0)
    for index, name in enumerate(("C", "H", "O")):
        inlet[name] += per_kg[index] * fuel_mass
    _add_gas(inlet, frozen.primary_feed_mol_s)
    _add_gas(inlet, frozen.secondary_feed_mol_s)
    outlet = {name: 0.0 for name in ("C", "H", "O")}
    _add_gas(outlet, state.freeboard_gas_mol_s[0])
    char_top = float(((frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * state.freeboard_char_free_kg[0]).sum())
    char_drain = float(frozen.drain_char_kg_s.sum())
    outlet["C"] += (char_top + char_drain) / (atomic_weight("C") * 1.0e-3)
    fuel_out = float(frozen.fuel_up_free_kg_s[0].sum() + frozen.fuel_exit_bed_kg_s.sum())
    for index, name in enumerate(("C", "H", "O")):
        outlet[name] += per_kg[index] * fuel_out
    inlet_h = _gas_enthalpy(frozen.primary_feed_mol_s, frozen.primary_temperature_k)
    inlet_h += _gas_enthalpy(frozen.secondary_feed_mol_s, frozen.secondary_temperature_k)
    inlet_h += _fuel_enthalpy(fuel_mass, frozen.fuel_feed_temperature_k, frozen)
    top = float(state.freeboard_temperature_k[0])
    outlet_h = _gas_enthalpy(state.freeboard_gas_mol_s[0], top)
    outlet_h += char_top * char_specific_enthalpy_j_per_kg(top)
    ash_top = float(((frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * state.freeboard_ash_free_kg[0]).sum())
    outlet_h += ash_top * ash_specific_heat_j_per_kg_k() * (top - _AMBIENT_K)
    outlet_h += _fuel_enthalpy(float(frozen.fuel_up_free_kg_s[0].sum()), top, frozen)
    bottom = float(state.bed_temperature_k[-1])
    outlet_h += char_drain * char_specific_enthalpy_j_per_kg(bottom)
    outlet_h += float(frozen.drain_ash_kg_s) * ash_specific_heat_j_per_kg_k() * (bottom - _AMBIENT_K)
    outlet_h += _fuel_enthalpy(float(frozen.fuel_exit_bed_kg_s.sum()), bottom, frozen)
    if frozen.hover_closure == "C" and redirect is not None:
        char_mass, _ash_mass, fuel_mass, enthalpy = redirect
        outlet["C"] += char_mass / (atomic_weight("C") * 1.0e-3)
        for index, name in enumerate(("C", "H", "O")):
            outlet[name] += per_kg[index] * fuel_mass
        outlet_h += enthalpy
    return (
        [ElementStream(inlet, inlet_h)], [ElementStream(outlet, outlet_h)],
        (inlet["C"], inlet["H"], inlet["O"]),
    )


def _fuel_moles(physical, mass):
    moisture = float(physical.moisture_mass_fraction)
    dry = mass * (1.0 - moisture)
    water = mass * moisture
    carbon = dry * float(physical.carbon_mass_fraction_wf) / (atomic_weight("C") * 1.0e-3)
    hydrogen = dry * float(physical.hydrogen_mass_fraction_wf) / (atomic_weight("H") * 1.0e-3)
    hydrogen += 2.0 * water / molar_mass_kg_mol("H2O")
    oxygen = dry * float(physical.oxygen_mass_fraction_wf) / (atomic_weight("O") * 1.0e-3)
    oxygen += water / molar_mass_kg_mol("H2O")
    return carbon, hydrogen, oxygen


def _add_gas(total, flow) -> None:
    for index, name in enumerate(GAS_SPECIES):
        for element, count in atoms(name).items():
            if element in total:
                total[element] += float(flow[index]) * count


def _gas_enthalpy(flow, temperature) -> float:
    return float(sum(float(flow[index]) * enthalpy_j_mol(name, temperature) for index, name in enumerate(GAS_SPECIES) if float(flow[index])))


def _fuel_enthalpy(mass, temperature, frozen) -> float:
    sensible = frozen.fuel_formation_j_per_kg + frozen.fuel_heat_capacity_j_per_kg_k * (float(temperature) - _AMBIENT_K)
    return float(mass) * sensible

def _info(result) -> dict:
    frozen = result.frozen
    names = (
        "t_ad_k", "t_init_k", "y_ad", "geldart", "d_bar_m", "d_sv_m", "bottom_residual_kg_s", "f0_a_kg_s",
        "hover_balance", "release_to_gas", "hydro_group", "geldart_record",
        "yield_notes", "slip_ratios", "bed_intervals", "substep_history",
    )
    copied = (
        "rho_p_kg_m3", "k_w_w_m2_k", "drain_ash_kg_s", "drain_char_kg_s", "ejection_to_free",
        "ejection_to_hover", "a_free_per_s", "a_hover_per_s", "b_free_per_s", "b_hover_per_s",
        "c_up_per_s", "c_down_per_s", "theta_up", "theta_down",
        "hover", "theta", "u_r_bed_m_s", "u_r_freeboard_m_s", "u_r_hover_m_s", "m_dot_ab_kg_s",
    )
    info = {name: getattr(result, name) for name in names}
    info.update({name: getattr(frozen, name) for name in copied})
    info["sensitivity"] = result.reported
    return info


def _with_sensitivity(info, state, frozen, assembly) -> dict:
    reported = dict(info.get("sensitivity") or {})
    visible = reported.get("visible_bubble_m3_s")
    if visible is not None:
        reported["bubble_balance_over_visible"] = [
            bubble_balance_ratio(
                float(state.bed_gas_mol_s[cell, 1].sum()), float(state.bed_temperature_k[cell]),
                float(frozen.bed_pressure_pa[cell]), float(visible[cell]),
            )
            for cell in range(frozen.n_bed)
        ]
    consumed = float(
        assembly.char_consumption_bed_kg_s.sum()
        + assembly.char_consumption_free_kg_s.sum()
        + assembly.char_consumption_hover_kg_s.sum()
    )
    hover = float(assembly.char_consumption_hover_kg_s.sum())
    reported["hover_char_consumption_fraction"] = 0.0 if consumed == 0.0 else hover / consumed
    info["sensitivity"] = reported
    return info


def _case_heat(path, case_id) -> float:
    from hamel.cases.assemble import load_case_file, resolve_case_key
    doc = load_case_file(path)
    key = resolve_case_key(doc, case_id)
    return float(doc["cases"][key]["wall_heat_loss"]["value"])


def _clock(deadline) -> None:
    if deadline is not None and deadline():
        raise SolveTimeout("运行超过墙钟上限")
