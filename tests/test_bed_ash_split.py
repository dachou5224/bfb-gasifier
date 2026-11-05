# -*- coding: utf-8 -*-
"""ADR 0055：床灰按飞灰份额分为细灰与粗灰、排灰与 E(w) 的惰性料衡算校正。"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.solver import campaign, driver
from hamel.solver.gates import OuterMonitor, ash_met
from hamel.solver.inert_balance import (
    ASH_RELATIVE_TOLERANCE, InertBalanceError, fine_share, free_ash_state, top_ash, update_fine_share,
)
from hamel.solver.precalc import PrecalcError, _drain_terms
from hamel.solver.residual import bed_ash_by_class, evaluate_residual
from hamel.solver.vector import pack_state, unpack_state, unknown_slots

from test_assembly import make_frozen, make_state

_CASE = Path("data/cases/vtt_pressurized.json")
_M_REF = 1.0e-3


def _ash_frozen(phi=0.7, share=0.0, target=7.0e-4, k_auf=5.0e-3):
    """一个能解出 E(w) 的冻结量：床灰抛出只进最小粒级的带出快、最大粒级的带出慢。"""
    return make_frozen(
        n_bed=1, n_free=2, n_class=2,
        k_auf_per_s=np.array([k_auf]),
        m_dot_ab_kg_s=np.array([0.0]),
        geometric_mass_kg=np.array([0.4]),
        a_free_per_s=np.array([[0.5, 0.05], [0.5, 0.05]]),
        b_free_per_s=np.array([[0.1, 0.02], [0.4, 0.4]]),
        fly_ash_fraction=phi, bed_ash_fine_share=share, ash_target_kg_s=target,
    )


def _ash_state(frozen):
    return make_state(frozen, char=1.0e-2, ash=1.0e-3)


def _response(frozen, state):
    """E(w) 的响应函数，记下子求解的迭代数。"""

    def call(w):
        trial = replace(frozen, bed_ash_fine_share=float(w))
        x, _count = free_ash_state(trial, state, SolverConfig(), _M_REF)
        return top_ash(trial, x)

    return call


# C12：括号求根的手算例，T = 1。
def test_hand_examples_of_the_endpoint_and_illinois_search():
    seen = []

    def coarse(w):
        seen.append(w)
        return 2.0 + 10.0 * w

    with pytest.raises(InertBalanceError) as caught:
        fine_share(coarse, 0.1, 1.0)
    assert "粗灰带出超过目标" in str(caught.value)
    assert caught.value.info["stop"] == "coarse"

    # E = 1 + w：端点 w = 0 已是解；下一轮仍保留 0（第 2 步直接满足）。
    w, info = fine_share(lambda w: 1.0 + w, 0.116, 1.0)
    assert w == 0.0 and info["E"] == pytest.approx(1.0)
    again, _info = fine_share(lambda w: 1.0 + w, 0.0, 1.0)
    assert again == 0.0

    # E = 0,5 + w，起点 w_f = 0 未达目标，Illinois 得 0,5。
    w, info = fine_share(lambda w: 0.5 + w, 0.0, 1.0)
    assert w == pytest.approx(0.5, abs=1.0e-9) and info["stop"] == "root"

    with pytest.raises(InertBalanceError) as caught:
        fine_share(lambda w: 0.2 + 0.5 * w, 0.1, 1.0)
    assert "细灰不足" in str(caught.value)

    # E ≡ 1：任意 w_f 已满足，保留当前值。
    for start in (0.0, 0.37, 1.0):
        kept, info = fine_share(lambda w: 1.0, start, 1.0)
        assert kept == start and info["stop"] == "satisfied"

    # 不连续函数使 50 次内不达容差：报告求根失败，评价次数达上限。
    def step(w):
        return 0.9 if w < 1.0 / 3.0 else 1.1

    with pytest.raises(InertBalanceError) as caught:
        fine_share(step, 0.5, 1.0, limit=50)
    assert "求根失败" in str(caught.value)
    assert caught.value.info["stop"] == "failed" and caught.value.info["evaluations"] >= 50


# C12b：端点响应只求一次（缓存命中不重算）。
def test_response_is_evaluated_once_per_distinct_share():
    calls = []

    def response(w):
        calls.append(w)
        return 0.5 + w

    w, info = fine_share(response, 0.0, 1.0)
    assert w == pytest.approx(0.5)
    assert len(calls) == len(set(calls)) == info["evaluations"]

    calls.clear()
    kept, info = fine_share(response, 1.0, 1.5)
    assert kept == 1.0 and info["stop"] == "satisfied"
    # 初值就是端点：该端点只求一次。
    assert calls == [0.0, 1.0]


# C13：φ、Ṁ_A 的输入有效域。
def test_invalid_phi_and_empty_ash_inputs_stop():
    released = SimpleNamespace(
        fuel_up_kg_s=np.zeros((1, 1)),
        undried_exit_kg_s=np.zeros(1), dried_exit_kg_s=np.zeros(1),
    )
    state = SimpleNamespace(bed_char_kg=np.array([[0.01]]))
    physical = SimpleNamespace(moisture_mass_fraction=0.0, ash_mass_fraction_wf=1.0e-3)
    for phi in (0.0, -0.5, 1.5):
        case = SimpleNamespace(physical=physical, fuel_kg_s=1.0, n_free=0, fly_ash_fraction=phi)
        with pytest.raises(PrecalcError, match="飞灰份额"):
            _drain_terms(case, state, released, np.zeros(1), np.array([0.4]))
    # Ṁ_A = 0：带入的灰全部随未转化燃料离开。
    case = SimpleNamespace(physical=physical, fuel_kg_s=1.0, n_free=0, fly_ash_fraction=0.5)
    all_out = SimpleNamespace(
        fuel_up_kg_s=np.zeros((1, 1)),
        undried_exit_kg_s=np.array([1.0]), dried_exit_kg_s=np.zeros(1),
    )
    with pytest.raises(PrecalcError, match="不是正数"):
        _drain_terms(case, state, all_out, np.zeros(1), np.array([0.4]))


# C15：床灰只进 k = 1 与 k = K；反应颗粒数与焦炭、气体行逐位不变。
def test_bed_ash_goes_only_to_the_smallest_and_largest_class():
    split = bed_ash_by_class(np.array([0.4, 0.6]), 0.25, 4)
    assert split[:, 0] == pytest.approx([0.1, 0.15])
    assert split[:, -1] == pytest.approx([0.3, 0.45])
    assert np.all(split[:, 1:-1] == 0.0)
    assert np.allclose(split.sum(axis=1), [0.4, 0.6])
    # 单粒级时两箱落到同一粒级，总量仍是 m_inert。
    assert bed_ash_by_class(np.array([0.4]), 0.25, 1) == pytest.approx(np.array([[0.4]]))


def test_bed_ash_split_leaves_char_and_gas_rows_bitwise_unchanged():
    """ADR 0044 的反应颗粒数只按焦炭质量计，细灰份额不改它。"""
    frozen = _ash_frozen(share=0.0)
    state = _ash_state(frozen)
    left = evaluate_residual(pack_state(state, frozen), frozen)
    right_frozen = replace(frozen, bed_ash_fine_share=1.0)
    right = evaluate_residual(pack_state(state, right_frozen), right_frozen)
    assert np.array_equal(left.char_consumption_bed_kg_s, right.char_consumption_bed_kg_s)
    assert np.array_equal(left.char_consumption_free_kg_s, right.char_consumption_free_kg_s)
    assert int(np.count_nonzero(left.values - right.values)) > 0


# C14：E(w) 按定义重新求自由空间灰稳态。
def test_e_of_w_resolves_the_freeboard_ash_steady_state():
    frozen = _ash_frozen(share=0.0)
    state = _ash_state(frozen)
    response = _response(frozen, state)
    e0, e1 = response(0.0), response(1.0)
    assert e0 != e1
    trial = replace(frozen, bed_ash_fine_share=0.3)
    values, iterations = free_ash_state(trial, state, SolverConfig(), _M_REF)
    assert iterations >= 1
    residual = evaluate_residual(values, trial)
    rows = [slot.index for slot in unknown_slots(trial) if slot.kind == "ash"]
    assert float(np.max(np.abs(residual.values[rows]))) / _M_REF <= 1.0e-6
    # 子求解的解确实随 w 变化，不是保持顶部灰存量不变。
    other, _count = free_ash_state(replace(frozen, bed_ash_fine_share=1.0), state, SolverConfig(), _M_REF)
    top = unpack_state(values, trial).freeboard_ash_free_kg[0]
    shifted = unpack_state(other, trial).freeboard_ash_free_kg[0]
    assert not np.allclose(top, shifted)


# C11：不动点上飞灰 = φ·Ṁ_A、底灰 = (1 − φ)·Ṁ_A，全炉灰平衡闭合；φ = 1 时排灰为 0。
def test_fixed_point_ash_balance_closes_with_the_prescribed_split():
    frozen = _ash_frozen(phi=0.7, share=0.0, target=7.0e-4)
    state = _ash_state(frozen)
    w, info = update_fine_share(frozen, state, SolverConfig(), _M_REF)
    fly = info["E"]
    assert 0.0 < w < 1.0
    assert fly == pytest.approx(info["target"], rel=1.0e-6)
    physical = SimpleNamespace(moisture_mass_fraction=0.0, ash_mass_fraction_wf=1.0e-3)
    case = SimpleNamespace(physical=physical, fuel_kg_s=1.0, n_free=2, fly_ash_fraction=0.7)
    released = SimpleNamespace(
        fuel_up_kg_s=np.zeros((1, 2)),
        undried_exit_kg_s=np.zeros((3, 2)), dried_exit_kg_s=np.zeros((3, 2)),
    )
    terms = _drain_terms(case, state, released, np.array([0.38]), np.array([0.4]))
    bottom = terms["ash"]
    assert terms["available"] == pytest.approx(1.0e-3)
    assert fly == pytest.approx(0.7 * terms["available"])
    assert bottom == pytest.approx(0.3 * terms["available"])
    assert fly + bottom == pytest.approx(terms["available"], rel=1.0e-6)
    # φ = 1：排灰（灰与焦炭）为 0。
    zero_case = SimpleNamespace(physical=physical, fuel_kg_s=1.0, n_free=2, fly_ash_fraction=1.0)
    zero = _drain_terms(zero_case, state, released, np.array([0.38]), np.array([0.4]))
    assert zero["ash"] == 0.0 and np.all(zero["char"] == 0.0)


# C16：外层收敛门与 (b) 复核。
def test_outer_gate_requires_the_ash_balance():
    monitor = OuterMonitor(np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), ash_relative=1.0e-3)
    assert not ash_met(monitor, ASH_RELATIVE_TOLERANCE)
    assert ash_met(replace(monitor, ash_relative=1.0e-5), ASH_RELATIVE_TOLERANCE)
    assert ash_met(replace(monitor, ash_relative=None), ASH_RELATIVE_TOLERANCE)
    frozen = _ash_frozen(phi=0.5, share=0.3, target=1.0e-3)
    state = make_state(frozen, char=1.0e-2, ash=1.0e-3)
    assembly = evaluate_residual(pack_state(state, frozen), frozen)
    relative = driver._ash_relative(state, frozen, assembly)
    assert relative is not None and relative >= ASH_RELATIVE_TOLERANCE
    # 把顶部带出灰调到目标上：复核通过。
    exact = replace(state, freeboard_ash_free_kg=state.freeboard_ash_free_kg.copy())
    exact.freeboard_ash_free_kg[0] = 1.0e-3 / float(np.sum(frozen.a_free_per_s[0]))
    near = driver._ash_relative(exact, frozen, assembly)
    assert near is not None and near < ASH_RELATIVE_TOLERANCE


# C17：S18 的 overlay 只改 φ。
def test_s18_overlay_changes_only_the_fly_ash_fraction(tmp_path):
    specs = {spec.name: spec for spec in campaign.run_specs()}
    assert "S18" in specs and "φ = 0.3" in specs["S18"].change
    dest = tmp_path / "S18.json"
    adiabatic = campaign.overlay(str(_CASE), str(dest), specs["S18"])
    assert adiabatic is False
    original = json.loads(_CASE.read_text(encoding="utf-8"))
    changed = json.loads(dest.read_text(encoding="utf-8"))
    original["cases"]["1"]["fly_ash_fraction"]["value"] = 0.3
    assert changed == original
    assert campaign.config_for("S18") == SolverConfig()
