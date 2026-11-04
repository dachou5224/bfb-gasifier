# -*- coding: utf-8 -*-
"""ADR 0052：初值构造的固相、排灰与床区向下流联立牛顿，以及交替与停止条件。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.properties.data import GAS_SPECIES, molar_mass_kg_mol
from hamel.pyrolysis.yields import PYROLYSIS_GAS_SPECIES
from hamel.solver import initial
from hamel.solver.initial import (
    DomainRejected, InitialStateError, SolidSolve, construct, recompute_flows, solve_solids,
)
from hamel.solver.newton import IDENTITY, build_rule
from hamel.solver.precalc import PrecalcError, SolidContext, _drain, _initial_tar, bed_downflow, bed_release_terms
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots, unpack_state

from test_assembly import make_frozen, make_state

# 工况 1 的燃料干基质量流，作 ADR 0012 的 ṁ_ref。
_M_REF = 0.013027777777777777 * (1.0 - 0.118)
_ASH_WF = 0.043
_MOISTURE = 0.118


def _context(n_bed=2, n_free=2, n_class=1, fuel_kg_s=0.013027777777777777, a_f=0.0, ash_out=0.0, fuel_out=0.0, loss_kg_s=0.0, fuel_exit_kg_s=0.0):
    """小例子的排灰与向下流输入。字段与预计算调用 _drain 时相同。

    loss_kg_s 以最底床格的水分流给出（ADR 0053 第 1 项），fuel_exit_kg_s 是最底格的燃料出口。
    """
    case = SimpleNamespace(
        physical=SimpleNamespace(moisture_mass_fraction=_MOISTURE, ash_mass_fraction_wf=_ASH_WF),
        fuel_kg_s=fuel_kg_s, n_free=n_free,
    )
    zeros = np.zeros((n_free + n_bed, n_class))
    released = SimpleNamespace(
        fuel_up_kg_s=np.zeros((1, n_class)), undried_exit_kg_s=zeros, dried_exit_kg_s=zeros.copy(),
        water_kg_s=np.zeros(n_free + n_bed),
        gas_mol_s={name: np.zeros(n_free + n_bed) for name in PYROLYSIS_GAS_SPECIES},
        benzene_mol_s=np.zeros(n_free + n_bed), naphthalene_mol_s=np.zeros(n_free + n_bed),
    )
    released.water_kg_s[-1] = float(loss_kg_s)
    released.undried_exit_kg_s[-1] = float(fuel_exit_kg_s)
    books = tuple(
        SimpleNamespace(cells=[SimpleNamespace(a_f=a_f) for _ in range(n_free)]) for _ in range(n_class)
    )
    return SolidContext(case, released, books, 0.0, ash_out, fuel_out)


def _system(generation=1.2e-3):
    """两格床区、两格自由空间、一个粒级的小例子。冻结的 Ṁ_ab 取 0.01 kg/s。

    温度 900 K、O₂ 取 1·10⁻⁴ 使单位焦炭消耗速率与生成量同量级，床区焦炭在几何存量附近。
    """
    frozen = make_frozen(
        n_bed=2, n_free=2, n_class=1,
        k_auf_per_s=np.array([0.025, 0.125]),
        m_dot_ab_kg_s=np.array([0.01, 0.0]),
        geometric_mass_kg=np.array([0.4, 0.4]),
        b_free_per_s=np.array([[0.0], [1.0]]),
        a_free_per_s=np.full((2, 1), 0.1),
        char_gen_bed_kg_s=np.full((2, 1), generation),
    )
    state = make_state(frozen, temperature=900.0, char=2.0e-3, ash=0.0)
    oxygen = GAS_SPECIES.index("O2")
    state.bed_gas_mol_s[:, :, oxygen] = 1.0e-4
    state.freeboard_gas_mol_s[:, oxygen] = 1.0e-4
    return frozen, state, _context()


def _inert(state, frozen):
    return frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)


def _char_rows(frozen):
    return np.array([slot.index for slot in unknown_slots(frozen) if slot.kind == "char"], dtype=int)


def _frozen_group_solution(frozen, state, context):
    """只按状态重算排灰、冻结 Ṁ_ab 的焦炭块解。

    焦炭行对焦炭未知量线性，故一趟逐列差分即精确系数，一次线性求解给出该组下的稳态。
    """
    rows = _char_rows(frozen)
    base = pack_state(state, frozen)
    scale = float(np.sum(frozen.geometric_mass_kg))

    def values(y):
        x = base.copy()
        x[rows] = y
        current = unpack_state(x, frozen)
        reaction = -evaluate_residual(x, frozen).char_consumption_bed_kg_s.sum(axis=1)
        _, drain_ash, drain_char, _, _, _ = recompute_flows(current, frozen, context, reaction)
        updated = replace(frozen, drain_ash_kg_s=drain_ash, drain_char_kg_s=drain_char)
        return evaluate_residual(x, updated).values[rows]

    r0 = values(base[rows])
    matrix = np.zeros((rows.size, rows.size))
    for position, index in enumerate(rows):
        step = np.sqrt(np.finfo(float).eps) * max(abs(float(base[index])), scale)
        trial = base.copy()
        trial[index] += step
        matrix[:, position] = (values(trial[rows]) - r0) / step
    x = base.copy()
    x[rows] += np.linalg.solve(matrix, -r0)
    return unpack_state(x, frozen)


def test_recomputed_drain_and_downflow_match_the_precalculation_under_the_same_state():
    """ADR 0052 实现约束：按当前状态重算的排灰与递推与预计算在同一状态上逐位相同。"""
    frozen, state, context = _system()
    reaction = -evaluate_residual(pack_state(state, frozen), frozen).char_consumption_bed_kg_s.sum(axis=1)
    inert, drain_ash, drain_char, down, closure, returned = recompute_flows(state, frozen, context, reaction)
    # 先按 ADR 0014 A5 与 ADR 0009 第 2 项手算。
    ash_in = context.case.fuel_kg_s * _ASH_WF * (1.0 - _MOISTURE)
    fraction = float(inert[-1] / frozen.geometric_mass_kg[-1])
    total = ash_in / fraction
    assert drain_ash == pytest.approx(ash_in)
    assert drain_char[0] == pytest.approx(total * state.bed_char_kg[-1, 0] / frozen.geometric_mass_kg[-1])
    up = np.array([0.01, 0.05])
    hand_returned = 1.0 * (state.freeboard_char_free_kg[-1, 0] + state.freeboard_ash_free_kg[-1, 0])
    assert returned == pytest.approx(hand_returned)
    assert down[0] == pytest.approx(up[1] - up[0] + hand_returned + reaction[0])
    assert down[-1] == 0.0
    assert closure == pytest.approx(
        context.case.fuel_kg_s + down[0] - up[1] - (drain_ash + drain_char.sum()) + reaction[-1]
    )
    # 再与预计算的同一函数、同一输入逐位对照。
    expected_ash, expected_char = _drain(
        context.case, state, context.released, inert, frozen.geometric_mass_kg, context.books, 0.0, 0.0,
    )
    expected_down, expected_closure = bed_downflow(
        up, context.case.fuel_kg_s, reaction, hand_returned, expected_ash + expected_char.sum(),
    )
    assert drain_ash == expected_ash
    assert np.array_equal(drain_char, expected_char)
    assert np.array_equal(down, np.asarray(expected_down, dtype=float))
    assert closure == expected_closure


def test_recomputed_downflow_adds_release_loss_and_the_bottom_fuel_exit():
    """ADR 0053 第 1、2 项：两处共用同一函数，同一状态、同一释放分布下逐位相同。

    补入两项后闭合到舍入；去掉任一项，闭合差恰为该项。
    """
    frozen, state, _ignored = _system()
    context = _context(loss_kg_s=1.0e-3, fuel_exit_kg_s=2.0e-4)
    reaction = -evaluate_residual(pack_state(state, frozen), frozen).char_consumption_bed_kg_s.sum(axis=1)
    _inert, drain_ash, drain_char, down, closure, returned = recompute_flows(state, frozen, context, reaction)
    loss, fuel_exit = bed_release_terms(context.released, frozen.n_freeboard)
    assert loss[-1] == pytest.approx(1.0e-3)
    assert fuel_exit == pytest.approx(2.0e-4)
    up = frozen.k_auf_per_s * frozen.geometric_mass_kg
    drain = float(drain_ash + np.sum(drain_char))
    expected_down, expected_closure = bed_downflow(up, context.case.fuel_kg_s, reaction, returned, drain, loss, fuel_exit)
    assert np.array_equal(down, np.asarray(expected_down, dtype=float))
    assert closure == expected_closure
    _down, without_loss = bed_downflow(up, context.case.fuel_kg_s, reaction, returned, drain, None, fuel_exit)
    _down, without_fuel = bed_downflow(up, context.case.fuel_kg_s, reaction, returned, drain, loss, 0.0)
    assert without_loss - closure == pytest.approx(float(np.sum(loss)))
    assert without_fuel - closure == pytest.approx(fuel_exit)


def test_frozen_downflow_gives_negative_bed_inert_and_recomputation_positive():
    """ADR 0052 实现约束：冻结 Ṁ_ab 时床格惰性物为负，按状态重算后为正。"""
    frozen, state, context = _system()
    frozen_state = _frozen_group_solution(frozen, state, context)
    assert float(np.min(_inert(frozen_state, frozen))) < 0.0
    solved = solve_solids(state, frozen, context, _M_REF, config=SolverConfig())
    assert float(np.min(solved.bed_inert_kg)) > 0.0
    assert np.all(solved.state.bed_char_kg > 0.0)
    assert solved.solid_residual <= 1.0e-6
    assert solved.m_dot_ab_kg_s[0] > frozen.m_dot_ab_kg_s[0]
    assert np.all(np.asarray(solved.drain_char_kg_s) >= 0.0) and solved.drain_ash_kg_s >= 0.0
    assert all(0.0 <= value < 1.0 for value in initial._volume_ratio(solved.state, frozen))
    # 气相与温度固定，只解固相。
    assert np.array_equal(solved.state.bed_gas_mol_s, state.bed_gas_mol_s)
    assert np.array_equal(solved.state.freeboard_gas_mol_s, state.freeboard_gas_mol_s)
    assert np.array_equal(solved.state.bed_temperature_k, state.bed_temperature_k)
    assert len(solved.merit_history) == solved.iterations + 1


def test_trial_feasibility_only_requires_a_positive_bottom_ash_fraction():
    """ADR 0054 第 2 项第 2 层：试探点只查幂变量 s ≥ 0、温度 > 0 与（初值构造）最底床格灰分 > 0。

    存量、床格惰性物与不变换气体出流的非负不在试探点上检查，只在收敛解上按 −τ、−τ_n 检查。
    """
    frozen, state, context = _system()
    rows = np.array([slot.index for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")], dtype=int)
    scale = float(np.sum(frozen.geometric_mass_kg))
    rule = build_rule(
        unknown_slots(frozen), frozen, IDENTITY,
        tau=SolverConfig().stock_roundoff_relative * scale, unknowns=rows, bottom_ash=True,
    )
    assert rule.trial(pack_state(state, frozen))
    # 负存量与负的床格惰性物在试探点上都可以出现。
    negative = replace(state, freeboard_ash_free_kg=np.full_like(state.freeboard_ash_free_kg, -1.0e-9))
    assert rule.trial(pack_state(negative, frozen))
    heavy = replace(state, bed_char_kg=state.bed_char_kg.copy())
    heavy.bed_char_kg[0, 0] = 0.5
    assert rule.trial(pack_state(heavy, frozen))
    carried = make_state(frozen, temperature=900.0, char=2.0e-3, ash=0.0)
    carried.freeboard_ash_free_kg[0] = 1.0
    assert rule.trial(pack_state(carried, frozen))
    # 最底床格灰分为零时排灰总流率没有定义，仍是不可行。
    ash_only = replace(state, bed_char_kg=np.zeros_like(state.bed_char_kg))
    ash_only.bed_char_kg[-1] = frozen.geometric_mass_kg[-1]
    assert not rule.trial(pack_state(ash_only, frozen))
    # 初值构造的未知量只有固相：温度与气体是固定的，不在试探点约束里重复检查。
    cold = replace(state, bed_temperature_k=np.full_like(state.bed_temperature_k, 0.0))
    assert rule.trial(pack_state(cold, frozen))


def test_negative_drain_balance_is_evaluated_at_trial_points(monkeypatch):
    """ADR 0052 修订：试探点上的负排灰余额照原样参与原残差，只在收敛解上检查非负。

    预计算仍用带符号保护的 ``_drain``（A5 口径不变）；初值固相牛顿走不带保护的余额式。
    """
    frozen, state, _ignored = _system()
    # 排灰用的 a^f 与残差里的输运系数取同一个值，避免两处不一致。
    context = _context(a_f=float(frozen.a_free_per_s[0, 0]))
    start = replace(state, freeboard_ash_free_kg=state.freeboard_ash_free_kg.copy())
    # 带出灰 = a^f × 最上一格自由灰 = 0,1 /s × 6·10⁻³ kg = 6·10⁻⁴ kg/s，略大于燃料带入灰 4,94·10⁻⁴ kg/s。
    start.freeboard_ash_free_kg[0, 0] = 6.0e-3
    inert, drain_ash, drain_char, _down, _closure, _returned = recompute_flows(start, frozen, context)
    assert drain_ash < 0.0
    assert float(np.sum(drain_char)) < 0.0
    # 同一状态在预计算里仍按 ADR 0014 A5 停止，余额不被截断。
    with pytest.raises(PrecalcError, match="带出灰多于燃料带入灰"):
        _drain(context.case, start, context.released, inert, frozen.geometric_mass_kg, context.books)

    balances = []
    original = initial._drain_terms

    def probe(*args, **kwargs):
        terms = original(*args, **kwargs)
        balances.append(terms["ash"])
        return terms

    monkeypatch.setattr(initial, "_drain_terms", probe)
    solved = solve_solids(start, frozen, context, _M_REF, config=SolverConfig())
    # 迭代中的试探点确实出现过负余额，且没有被截断或拒评。
    assert min(balances) < 0.0
    assert solved.drain_ash_kg_s >= 0.0
    assert float(np.sum(solved.drain_char_kg_s)) >= 0.0


def test_converged_solution_violations_stop_and_are_reported():
    """ADR 0054 第 2 项第 3 层：收敛解违反任一条即停止并报告；存量带舍入容差 τ。"""
    frozen, state, context = _system()
    solved = solve_solids(state, frozen, context, _M_REF, config=SolverConfig())
    rows = np.array([slot.index for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")], dtype=int)
    scale = float(np.sum(frozen.geometric_mass_kg))
    tau = SolverConfig().stock_roundoff_relative * scale
    rule = build_rule(unknown_slots(frozen), frozen, IDENTITY, tau=tau, unknowns=rows, bottom_ash=True)
    assert rule.final(pack_state(solved.state, frozen), solved.drain_ash_kg_s)["ok"]
    # 舍入容差的边界：−τ 通过，−2τ 停止。
    at_tau = replace(solved.state, bed_char_kg=solved.state.bed_char_kg.copy())
    at_tau.bed_char_kg[0, 0] = -tau
    assert rule.final(pack_state(at_tau, frozen), solved.drain_ash_kg_s)["ok"]
    below_tau = replace(solved.state, bed_char_kg=solved.state.bed_char_kg.copy())
    below_tau.bed_char_kg[0, 0] = -2.0 * tau
    detail = rule.final(pack_state(below_tau, frozen), solved.drain_ash_kg_s)
    assert not detail["ok"] and "负存量" in detail["message"]
    # 惰性物为负：几何存量扣去焦炭与燃料后为负。
    heavy = replace(solved.state, bed_char_kg=solved.state.bed_char_kg.copy())
    heavy.bed_char_kg[0, 0] = 1.2 * float(frozen.geometric_mass_kg[0])
    detail = rule.final(pack_state(heavy, frozen), solved.drain_ash_kg_s)
    assert not detail["ok"] and "惰性物为负" in detail["message"]
    # 排灰灰分量为负、最底床格灰分为零。
    negative = replace(solved.state, bed_char_kg=np.zeros_like(solved.state.bed_char_kg))
    detail = rule.final(pack_state(negative, frozen), -1.0e-6)
    assert not detail["ok"] and "排灰灰分量" in detail["message"]
    ash_only = replace(solved.state, bed_char_kg=np.zeros_like(solved.state.bed_char_kg))
    ash_only.bed_char_kg[-1] = frozen.geometric_mass_kg[-1]
    detail = rule.final(pack_state(ash_only, frozen), solved.drain_ash_kg_s)
    assert not detail["ok"] and "灰分为零" in detail["message"]
    # 固相门仍由 solve_solids 在收敛解上复核（converged 用同一个 gate），不再重复断言。


def test_signed_drain_allows_a_negative_a5_balance():
    """ADR 0054 第 3 项第 5 步：交替中间的预计算取同式的负值；外层与 guarded 调用仍按符号停止。"""
    frozen, state, _ignored = _system()
    context = _context(a_f=float(frozen.a_free_per_s[0, 0]))
    start = replace(state, freeboard_ash_free_kg=state.freeboard_ash_free_kg.copy())
    start.freeboard_ash_free_kg[0, 0] = 6.0e-3
    inert = frozen.geometric_mass_kg - start.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    with pytest.raises(PrecalcError, match="带出灰多于燃料带入灰"):
        _drain(context.case, start, context.released, inert, frozen.geometric_mass_kg, context.books)
    ash, char = _drain(
        context.case, start, context.released, inert, frozen.geometric_mass_kg, context.books,
        0.0, 0.0, guard_sign=False,
    )
    assert ash < 0.0 and float(np.sum(char)) < 0.0


def test_construct_stops_when_the_last_precalculation_has_a_negative_drain(monkeypatch):
    """ADR 0054 第 3 项第 5 步：交替结束时的排灰灰分量必须非负，否则停止并报告。"""
    frozen = make_frozen()
    state = make_state(frozen)

    def precalc(previous, drain):
        return SimpleNamespace(frozen=replace(frozen, drain_ash_kg_s=-1.0e-5), state=previous, solid_context=None)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    with pytest.raises(InitialStateError, match="排灰灰分量"):
        construct(_opened(frozen, state), precalc, 1.0e-3)


def test_domain_rejection_at_the_start_keeps_the_state_for_the_report():
    """ADR 0052 第 3 项：起点落在 ADR 0026 第 5 项的有效域外时停止，并留下未截断的固相。"""
    frozen, state, context = _system()
    with pytest.raises(DomainRejected) as caught:
        solve_solids(make_state(frozen, temperature=900.0, char=50.0, ash=50.0), frozen, context, _M_REF)
    assert caught.value.phase == "起点"
    assert caught.value.state.freeboard_char_free_kg.shape == state.freeboard_char_free_kg.shape


def test_construct_records_the_alternation_when_the_precalc_stops(monkeypatch):
    """ADR 0052 报告项：重做预计算报错时，那次交替的固相解与 A5 各项仍写进交替记录。"""
    frozen = make_frozen()
    state = make_state(frozen)
    terms = {
        "ash_in": 4.9e-4, "carried": 8.7e-4, "fuel_ash": 5.7e-6, "redirected": 0.0,
        "ash_out": 8.7e-4, "fraction": 0.07, "ash": -3.8e-4, "char": np.array([-1.0e-3]),
    }

    def precalc(previous, drain):
        error = PrecalcError("带出灰多于燃料带入灰")
        error.drain_terms = terms
        raise error

    records: list = []
    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    with pytest.raises(PrecalcError, match="带出灰多于燃料带入灰"):
        construct(_opened(frozen, state), precalc, 1.0e-3, records=records)
    assert len(records) == 1
    record = records[0]
    # 固相牛顿的那一次已收敛，记录里保留它的迭代数与解出的存量。
    assert record["iterations"] == 1 and record["bed_char_kg"] == pytest.approx(1.0)
    assert record["precalc_error"] == "带出灰多于燃料带入灰"
    assert record["a5"]["ash"] == pytest.approx(-3.8e-4)
    assert record["a5"]["fuel_ash"] == pytest.approx(5.7e-6)
    assert record["a5"]["char"] == pytest.approx([-1.0e-3])


def _stub(result):
    """一次固相牛顿的最小返回。排灰取非零值，供交替回调核对回馈。"""
    return SolidSolve(
        result, 1, (1.0, 1.0e-9), 1.0e-9, 1.5e-4, np.array([2.5e-4]),
        np.zeros(2), np.zeros(2), 0.0, 1.0, np.zeros(2),
    )


def _opened(frozen, state):
    return SimpleNamespace(frozen=frozen, state=state, solid_context=None)


def test_construct_stops_after_max_loops(monkeypatch):
    """交替至多 5 次；未达判据时以最后一次的固相为初值并写明。ADR 0045 第 3 项第 5 步。"""
    frozen = make_frozen()
    state = make_state(frozen)
    calls = {"n": 0}

    def precalc(previous, drain):
        calls["n"] += 1
        return _opened(replace(frozen, k_auf_per_s=frozen.k_auf_per_s + 1.0e-4 * calls["n"]), previous)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(_opened(frozen, state), precalc, 1.0e-3)
    assert len(records) == 5 and calls["n"] == 5
    assert "交替 5 次" in note
    assert result.state is state


def test_construct_stops_when_flows_settle(monkeypatch):
    """相邻两次的界面流与回落相对变化小于 1·10⁻² 时停止；收敛固相的排灰传给下一次预计算。"""
    frozen = make_frozen()
    state = make_state(frozen)
    seen = []

    def precalc(previous, drain):
        seen.append(drain)
        return _opened(frozen, previous)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(_opened(frozen, state), precalc, 1.0e-3)
    assert note == "" and len(seen) == 1
    assert all(value == 0.0 for value in records[-1]["change"].values())
    assert records[-1]["iterations"] == 1
    # ADR 0052 第 1 项第 5 步：排灰随收敛固相回馈，不用上一轮预计算的冻结值。
    assert seen[0][0] == pytest.approx(1.5e-4)
    assert np.array_equal(seen[0][1], np.array([2.5e-4]))
    assert result.state is state


def test_construct_second_precalc_carries_the_bed_sink_without_resetting(monkeypatch):
    """集成：重做预计算的回调已带非零床区反应汇与改送，来自 driver.solid_feedback。"""
    from hamel.solver.driver import solid_feedback

    frozen = make_frozen(char_gen_bed_kg_s=np.full((2, 2), 1.0e-3))
    state = make_state(frozen, char=2.0e-3, ash=0.0)
    seen = []

    def precalc(previous, drain):
        reaction, redirect = solid_feedback(previous, frozen)
        seen.append({"reaction": reaction, "redirect": redirect})
        return _opened(frozen, previous)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(_opened(frozen, state), precalc, 1.0e-3)
    assert len(seen) == 1
    assert float(seen[0]["reaction"].sum()) < 0.0
    assert len(seen[0]["redirect"]) == 4
    assert note == "" and result.state is state


def test_solid_feedback_supplies_the_bed_sink_and_the_redirect():
    """求解手段：ADR 0052 第 1 项第 5 步；PR #130 审查 P2。数据流与驱动共用同一个函数。"""
    from hamel.solver.driver import solid_feedback

    frozen = make_frozen(char_gen_bed_kg_s=np.full((2, 2), 1.0e-3))
    state = make_state(frozen, char=2.0e-3, ash=0.0)
    reaction, redirect = solid_feedback(state, frozen)
    assert reaction.shape == (frozen.n_bed,)
    # 床区焦炭被消耗，递推的输入为负。
    assert float(reaction.sum()) < 0.0
    assert len(redirect) == 4
    assert all(isinstance(value, float) for value in redirect)


@dataclass
class _Case:
    n_bed: int
    n_free: int
    bed_cell_m: float
    free_cell_m: float
    fuel_kg_s: float
    physical: SimpleNamespace


def test_tar_is_added_at_and_above_the_feed_cell():
    """ADR 0045 第 2 项：燃料加入格及以上的苯与萘按 T_init 下的产率另加，不经反应。"""
    case = _Case(
        2, 2, 0.2, 0.2, 1.0,
        SimpleNamespace(fuel_feed_height=0.5, moisture_mass_fraction=0.1, ash_mass_fraction_wf=0.05),
    )
    yields = SimpleNamespace(benzene_mol_per_kg=0.3, naphthalene_mol_per_kg=0.1)
    bed, free = _initial_tar(case, yields)
    # 床面 0.4 m，0.5 m 落在自由空间最下格；该格及其上方各格都另加焦油。
    assert bed.sum() == 0.0
    scale = 1.0 * (1.0 - 0.1) * (1.0 - 0.05)
    assert free[-1].sum() == pytest.approx(0.4 * scale)
    assert free[0].sum() == pytest.approx(0.4 * scale)


def test_tar_is_not_added_below_the_feed_cell():
    case = _Case(
        2, 2, 0.2, 0.2, 1.0,
        SimpleNamespace(fuel_feed_height=0.3, moisture_mass_fraction=0.0, ash_mass_fraction_wf=0.0),
    )
    yields = SimpleNamespace(benzene_mol_per_kg=0.3, naphthalene_mol_per_kg=0.1)
    bed, free = _initial_tar(case, yields)
    # 0.3 m 在床区顶格（0.2–0.4 m）；其下的床格不加焦油。
    assert free.sum() == pytest.approx(2 * 0.4)
    assert bed[0].sum() == pytest.approx(0.4)
    assert bed[-1].sum() == 0.0
