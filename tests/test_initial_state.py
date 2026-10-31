# -*- coding: utf-8 -*-
"""ADR 0052：初值构造的固相、排灰与床区向下流联立牛顿，以及交替与停止条件。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.properties.data import GAS_SPECIES
from hamel.solver import initial
from hamel.solver.initial import (
    DomainRejected, InitialStateError, SolidSolve, construct, recompute_flows,
    require_converged_solution, solid_feasible, solve_solids,
)
from hamel.solver.precalc import SolidContext, _drain, _initial_tar, bed_downflow
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots, unpack_state

from test_assembly import make_frozen, make_state

# 工况 1 的燃料干基质量流，作 ADR 0012 的 ṁ_ref。
_M_REF = 0.013027777777777777 * (1.0 - 0.118)
_ASH_WF = 0.043
_MOISTURE = 0.118


def _context(n_bed=2, n_free=2, n_class=1, fuel_kg_s=0.013027777777777777, a_f=0.0, ash_out=0.0, fuel_out=0.0):
    """小例子的排灰与向下流输入。字段与预计算调用 _drain 时相同。"""
    case = SimpleNamespace(
        physical=SimpleNamespace(moisture_mass_fraction=_MOISTURE, ash_mass_fraction_wf=_ASH_WF),
        fuel_kg_s=fuel_kg_s, n_free=n_free,
    )
    zeros = np.zeros((n_free + n_bed, n_class))
    released = SimpleNamespace(
        fuel_up_kg_s=np.zeros((1, n_class)), undried_exit_kg_s=zeros, dried_exit_kg_s=zeros.copy(),
    )
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


def test_solid_feasibility_rejects_negative_inert_and_negative_drain_ash():
    """ADR 0052 第 1 项第 3 步：可行条件为 ADR 0017 第 7 项，另加排灰灰分量非负。"""
    frozen, state, context = _system()
    assert solid_feasible(state, frozen, context)
    heavy = replace(state, bed_char_kg=np.full_like(state.bed_char_kg, 0.5))
    assert not solid_feasible(heavy, frozen, context)
    negative = replace(state, freeboard_ash_free_kg=np.full_like(state.freeboard_ash_free_kg, -1.0e-9))
    assert not solid_feasible(negative, frozen, context)
    carried = make_state(frozen, temperature=900.0, char=2.0e-3, ash=0.0)
    carried.freeboard_ash_free_kg[0] = 1.0
    heavy_ash = _context(a_f=1.0)
    assert not solid_feasible(carried, frozen, heavy_ash)


def test_converged_solution_violations_stop_and_are_reported():
    """ADR 0052 第 3 项：收敛解违反任一条即停止并报告。"""
    frozen, state, context = _system()
    solved = solve_solids(state, frozen, context, _M_REF, config=SolverConfig())
    flows = (
        solved.bed_inert_kg, solved.drain_ash_kg_s, solved.drain_char_kg_s,
        solved.m_dot_ab_kg_s, solved.bottom_closure_kg_s, 0.0,
    )
    rows = [slot.index for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")]
    residual = np.zeros(len(rows))
    assert require_converged_solution(solved.state, frozen, residual, flows, _M_REF, 1.0e-6) == 0.0
    bad_stock = replace(solved.state, bed_char_kg=solved.state.bed_char_kg.copy())
    bad_stock.bed_char_kg[0, 0] = -1.0e-3
    with pytest.raises(InitialStateError, match="负存量"):
        require_converged_solution(bad_stock, frozen, residual, flows, _M_REF, 1.0e-6)
    bad_inert = (solved.bed_inert_kg.copy(),) + flows[1:]
    bad_inert[0][0] = -1.0e-3
    with pytest.raises(InitialStateError, match="惰性物为负"):
        require_converged_solution(solved.state, frozen, residual, bad_inert, _M_REF, 1.0e-6)
    bad_drain = flows[:1] + (-1.0e-6,) + flows[2:]
    with pytest.raises(InitialStateError, match="排灰灰分量"):
        require_converged_solution(solved.state, frozen, residual, bad_drain, _M_REF, 1.0e-6)
    zero_fraction = (np.zeros_like(solved.bed_inert_kg),) + flows[1:]
    with pytest.raises(InitialStateError, match="灰分为零"):
        require_converged_solution(solved.state, frozen, residual, zero_fraction, _M_REF, 1.0e-6)
    with pytest.raises(InitialStateError, match="固相门"):
        require_converged_solution(solved.state, frozen, np.full(len(rows), 2.0 * _M_REF), flows, _M_REF, 1.0e-6)


def test_domain_rejection_at_the_start_keeps_the_state_for_the_report():
    """ADR 0052 第 3 项：起点落在 ADR 0026 第 5 项的有效域外时停止，并留下未截断的固相。"""
    frozen, state, context = _system()
    with pytest.raises(DomainRejected) as caught:
        solve_solids(make_state(frozen, temperature=900.0, char=50.0, ash=50.0), frozen, context, _M_REF)
    assert caught.value.phase == "起点"
    assert caught.value.state.freeboard_char_free_kg.shape == state.freeboard_char_free_kg.shape


def _stub(result):
    return SolidSolve(result, 1, (1.0, 1.0e-9), 1.0e-9, 0.0, np.zeros(1), np.zeros(2), np.zeros(2), 0.0, 1.0, np.zeros(2))


def _opened(frozen, state):
    return SimpleNamespace(frozen=frozen, state=state, solid_context=None)


def test_construct_stops_after_max_loops(monkeypatch):
    """交替至多 5 次；未达判据时以最后一次的固相为初值并写明。ADR 0045 第 3 项第 5 步。"""
    frozen = make_frozen()
    state = make_state(frozen)
    calls = {"n": 0}

    def precalc(previous):
        calls["n"] += 1
        return _opened(replace(frozen, k_auf_per_s=frozen.k_auf_per_s + 1.0e-4 * calls["n"]), previous)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(_opened(frozen, state), precalc, 1.0e-3)
    assert len(records) == 5 and calls["n"] == 5
    assert "交替 5 次" in note
    assert result.state is state


def test_construct_stops_when_flows_settle(monkeypatch):
    """相邻两次的界面流与回落相对变化小于 1·10⁻² 时停止。"""
    frozen = make_frozen()
    state = make_state(frozen)
    calls = {"n": 0}

    def precalc(previous):
        calls["n"] += 1
        return _opened(frozen, previous)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(_opened(frozen, state), precalc, 1.0e-3)
    assert note == "" and calls["n"] == 1
    assert all(value == 0.0 for value in records[-1]["change"].values())
    assert records[-1]["iterations"] == 1
    assert result.state is state


def test_construct_second_precalc_carries_the_bed_sink_without_resetting(monkeypatch):
    """集成：重做预计算的回调已带非零床区反应汇与改送，来自 driver.solid_feedback。"""
    from hamel.solver.driver import solid_feedback

    frozen = make_frozen(char_gen_bed_kg_s=np.full((2, 2), 1.0e-3))
    state = make_state(frozen, char=2.0e-3, ash=0.0)
    seen = []

    def precalc(previous):
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
