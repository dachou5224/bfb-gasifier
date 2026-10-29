# -*- coding: utf-8 -*-
"""ADR 0045：焦炭块与灰块的顺序线性求解、块结构与停止条件。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.solver import initial
from hamel.solver.initial import (
    DomainRejected, InitialStateError, SolidSolve, construct, solve_solids,
)
from hamel.solver.precalc import _initial_tar
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots

from test_assembly import make_frozen, make_state


def test_char_rows_do_not_depend_on_ash_unknowns():
    """ADR 0045 第 3 项：焦炭行对任一灰未知量的导数逐位为零。"""
    frozen = make_frozen(char_gen_bed_kg_s=np.full((2, 2), 1.0e-3))
    state = make_state(frozen)
    x = pack_state(state, frozen)
    base = evaluate_residual(x, frozen).values
    rows = [slot.index for slot in unknown_slots(frozen) if slot.kind == "char"]
    for slot in unknown_slots(frozen):
        if slot.kind != "ash":
            continue
        trial = x.copy()
        trial[slot.index] += 1.0e-9
        assert np.array_equal(evaluate_residual(trial, frozen).values[rows], base[rows])
    solved = solve_solids(state, frozen, m_dot_ref=0.5)
    assert solved.char_ash_derivative == 0.0


def test_block_order_solution_satisfies_the_solid_rows():
    """小例子上块顺序解与回代残差在舍入误差内。"""
    frozen = make_frozen(char_gen_bed_kg_s=np.full((2, 2), 1.0e-3))
    solved = solve_solids(make_state(frozen), frozen, m_dot_ref=0.5)
    assert solved.solid_residual < 1.0e-9
    assert solved.char_residual < 1.0e-9
    assert np.all(solved.state.bed_char_kg >= 0.0)
    assert np.all(solved.state.freeboard_ash_free_kg >= 0.0)
    # 回代后完整固相向量在 evaluate_residual 下不过门之外没有别的差异。
    check = evaluate_residual(pack_state(solved.state, frozen), frozen).values
    solid = [slot.index for slot in unknown_slots(frozen) if slot.kind in ("char", "ash")]
    assert np.max(np.abs(check[solid])) < 1.0e-9


def test_bed_char_exceeding_the_geometric_mass_stops():
    """解出 m_inert < 0 属负存量，按 ADR 0045 第 4 项停止。"""
    frozen = make_frozen(
        char_gen_bed_kg_s=np.full((2, 2), 1.0e-3), geometric_mass_kg=np.full(2, 1.0e-4),
    )
    with pytest.raises(InitialStateError, match="惰性物为负"):
        solve_solids(make_state(frozen), frozen, m_dot_ref=0.5)


def test_zero_bed_char_stops():
    """床区焦炭总量为零时灰分配没有定义，停止上报。"""
    frozen = make_frozen()
    state = make_state(frozen, char=0.0, ash=0.0)
    with pytest.raises(InitialStateError):
        solve_solids(state, frozen, m_dot_ref=0.5)


def test_domain_rejection_keeps_the_state_for_the_report():
    """有效域拒绝带上未截断的固相与阶段，便于报告。"""
    frozen = make_frozen(
        char_gen_bed_kg_s=np.full((2, 2), 1.0e-3), ash_gen_free_kg_s=np.full((2, 2), 2.0),
    )
    with pytest.raises(DomainRejected) as caught:
        solve_solids(make_state(frozen), frozen, m_dot_ref=0.5)
    assert caught.value.phase in {"焦炭块", "自由空间灰块", "回代"}
    assert caught.value.state.freeboard_char_free_kg.shape == (2, 2)


def _stub(result):
    return SolidSolve(result, 0.0, 0.0, 0.0, 1.0e-3, 1.0e-2)


def test_construct_stops_after_max_loops(monkeypatch):
    """交替至多 5 次；未达判据时以最后一次的固相为初值并写明。"""
    frozen = make_frozen()
    state = make_state(frozen)
    calls = {"n": 0}

    def precalc(previous):
        calls["n"] += 1
        return SimpleNamespace(
            frozen=replace(frozen, k_auf_per_s=frozen.k_auf_per_s + 1.0e-4 * calls["n"]),
            state=previous,
        )

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(SimpleNamespace(frozen=frozen, state=state), precalc, 1.0e-3)
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
        return SimpleNamespace(frozen=frozen, state=previous)

    monkeypatch.setattr(initial, "solve_solids", lambda *args, **kwargs: _stub(state))
    result, records, note = construct(SimpleNamespace(frozen=frozen, state=state), precalc, 1.0e-3)
    assert note == "" and calls["n"] == 1
    assert all(value == 0.0 for value in records[0]["change"].values())
    assert result.state is state


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
