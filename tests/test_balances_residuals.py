# -*- coding: utf-8 -*-
"""式 (2.4)–(2.8) 的残差。期望值是手算和守恒，不用旧代码输出。"""
from __future__ import annotations

import inspect

import numpy as np
import pytest

from hamel.balances import (
    ResidualInputError,
    energy_residual,
    gas_mole_residual,
    solid_mass_flow_residual,
    stream_enthalpy_flow,
)

def _zeros(shape):
    return np.zeros(shape)


def test_gas_hand_calculation_three_cells():
    # 下标 0 是最上层。只给组分 0 一组完整数字；组分 1 只放交换项，用来核对正负号。
    # 悬浮相，组分 0：
    #   外部加入 [1, 0, 4]，再循环 [0, 2, 0]，反应 [0, 1, 0]，出流 [8, 3, 5]
    #   来自下一层的入流是下一层出流：[3, 5, 0]
    #   交换 [1, 0, -3]，悬浮相减去它
    #   cell 0：1 + 0 + 3 + 0 − 8 − 1 = −5
    #   cell 1：0 + 2 + 5 + 1 − 3 − 0 = 5
    #   cell 2：4 + 0 + 0 + 0 − 5 − (−3) = 2
    # 组分 1 的交换是 [2, 0, −3]，其余为 0，所以悬浮相残差是 [−2, 0, 3]。
    # 气泡相，组分 0：
    #   外部加入 [0, 0, 6]，再循环 [1, 0, 0]，反应 [0, 0, 1]，出流 [2, 4, 7]
    #   来自下一层：[4, 7, 0]；交换仍是 [1, 0, −3]，气泡相加上它
    #   cell 0：0 + 1 + 4 + 0 − 2 + 1 = 4
    #   cell 1：0 + 0 + 7 + 0 − 4 + 0 = 3
    #   cell 2：6 + 0 + 0 + 1 − 7 + (−3) = −3
    # 组分 1 的气泡相残差是交换本身：[2, 0, −3]。
    external_suspension = np.array([[1.0, 0.0], [0.0, 0.0], [4.0, 0.0]])
    recycle_suspension = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 0.0]])
    reaction_suspension = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    out_suspension = np.array([[8.0, 0.0], [3.0, 0.0], [5.0, 0.0]])
    external_bubble = np.array([[0.0, 0.0], [0.0, 0.0], [6.0, 0.0]])
    recycle_bubble = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    reaction_bubble = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
    out_bubble = np.array([[2.0, 0.0], [4.0, 0.0], [7.0, 0.0]])
    exchange = np.array([[1.0, 2.0], [0.0, 0.0], [-3.0, -3.0]])

    result = gas_mole_residual(
        external_suspension,
        recycle_suspension,
        reaction_suspension,
        out_suspension,
        external_bubble,
        recycle_bubble,
        reaction_bubble,
        out_bubble,
        exchange,
    )

    assert np.array_equal(result.suspension_mol_s, np.array([[-5.0, -2.0], [5.0, 0.0], [2.0, 3.0]]))
    assert np.array_equal(result.bubble_mol_s, np.array([[4.0, 2.0], [3.0, 0.0], [-3.0, -3.0]]))
    assert np.array_equal(result.top_suspension_out_mol_s, np.array([8.0, 0.0]))
    assert np.array_equal(result.top_bubble_out_mol_s, np.array([2.0, 0.0]))


def test_gas_conservation_cancels_internal_flows_and_exchange():
    # 把所有 cell、两相残差相加：相邻的上行出流两两抵消，交换项也抵消。
    # 剩下外部加入、再循环、反应源，再减去离开最上层的两相出流。
    external_suspension = np.array([[1.0, 0.0], [0.0, 0.0], [4.0, 0.0]])
    recycle_suspension = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 0.0]])
    reaction_suspension = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    out_suspension = np.array([[8.0, 0.0], [3.0, 0.0], [5.0, 0.0]])
    external_bubble = np.array([[0.0, 0.0], [0.0, 0.0], [6.0, 0.0]])
    recycle_bubble = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    reaction_bubble = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
    out_bubble = np.array([[2.0, 0.0], [4.0, 0.0], [7.0, 0.0]])
    exchange = np.array([[1.0, 2.0], [0.0, 0.0], [-3.0, -3.0]])
    result = gas_mole_residual(
        external_suspension,
        recycle_suspension,
        reaction_suspension,
        out_suspension,
        external_bubble,
        recycle_bubble,
        reaction_bubble,
        out_bubble,
        exchange,
    )
    sources = (
        external_suspension
        + recycle_suspension
        + reaction_suspension
        + external_bubble
        + recycle_bubble
        + reaction_bubble
    )
    summed = result.suspension_mol_s.sum(axis=0) + result.bubble_mol_s.sum(axis=0)
    assert np.array_equal(summed + result.top_suspension_out_mol_s + result.top_bubble_out_mol_s, sources.sum(axis=0))


def test_gas_external_feed_enters_only_the_bottom_then_the_cell_above():
    # 只在最下层悬浮相加入 9 mol/s。其它流都是 0 时，只有该 cell 的残差是 9。
    # 再把这 9 mol/s 写成该 cell 的向上出流：本 cell 残差回到 0，上一层残差变为 9。
    external = _zeros((3, 1))
    external[2, 0] = 9.0
    zeros = _zeros((3, 1))
    before = gas_mole_residual(external, zeros, zeros, zeros, zeros, zeros, zeros, zeros, zeros)
    assert np.array_equal(before.suspension_mol_s, np.array([[0.0], [0.0], [9.0]]))
    assert np.array_equal(before.bubble_mol_s, zeros)

    out = _zeros((3, 1))
    out[2, 0] = 9.0
    after = gas_mole_residual(external, zeros, zeros, out, zeros, zeros, zeros, zeros, zeros)
    assert np.array_equal(after.suspension_mol_s, np.array([[0.0], [9.0], [0.0]]))


def test_solid_hand_calculation_matches_ten_flow_terms():
    # 一种固体，3 个 cell、2 个粒级。中间 cell、粒级 0 的十项是：
    #   外部加入 1，再循环 2，离开体系 3，本 cell 上行 4，本 cell 下行 5，
    #   下方 cell 同一粒级的上行 6，上方 cell 同一粒级的下行 7，
    #   反应 8，自其它粒级进入 9，进入其它粒级 10。
    #   1 + 2 − 3 − 4 − 5 + 6 + 7 + 8 + 9 − 10 = 11
    # 粒级 1 不混入这一项：下方上行只加到上一层的同一粒级。
    # cell 0 粒级 0：1 + 0 − 0 − 11 − 7 + 4 + 0 + 0 + 0 − 0 = −13
    # cell 2 粒级 0：0 + 0 − 0 − 6 − 9 + 0 + 5 + 0 + 0 − 0 = −10
    # cell 0 粒级 1：10 − 2 = 8
    # cell 1 粒级 1：下方上行 1 + 上方下行 2 = 3
    # cell 2 粒级 1：−1
    external = np.array([[1.0, 10.0], [1.0, 0.0], [0.0, 0.0]])
    recycle = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 0.0]])
    leaving = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 0.0]])
    reaction = np.array([[0.0, 0.0], [8.0, 0.0], [0.0, 0.0]])
    from_class = np.array([[0.0, 0.0], [9.0, 0.0], [0.0, 0.0]])
    to_class = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 0.0]])
    upward = np.array([[11.0, 0.0], [4.0, 0.0], [6.0, 1.0]])
    downward = np.array([[7.0, 2.0], [5.0, 0.0], [9.0, 0.0]])

    ten_terms = 1.0 + 2.0 - 3.0 - 4.0 - 5.0 + 6.0 + 7.0 + 8.0 + 9.0 - 10.0
    assert ten_terms == 11.0

    result = solid_mass_flow_residual(
        external, recycle, leaving, reaction, from_class, to_class, upward, downward
    )
    assert np.array_equal(result.residual_kg_s, np.array([[-13.0, 8.0], [11.0, 3.0], [-10.0, -1.0]]))
    assert result.residual_kg_s[1, 0] == ten_terms
    assert np.array_equal(result.top_upward_out_kg_s, np.array([11.0, 0.0]))
    assert np.array_equal(result.bottom_downward_out_kg_s, np.array([9.0, 0.0]))


def test_solid_conservation_leaves_only_external_terms_and_boundary_flows():
    # 内部的上行、下行逐对抵消。全床残差之和加上最上层上行、最下层下行，
    # 等于外部加入 + 再循环 − 离开体系 + 反应 + 粒级迁入 − 粒级迁出。
    external = np.array([[1.0, 10.0], [1.0, 0.0], [0.0, 0.0]])
    recycle = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 0.0]])
    leaving = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 0.0]])
    reaction = np.array([[0.0, 0.0], [8.0, 0.0], [0.0, 0.0]])
    from_class = np.array([[0.0, 0.0], [9.0, 0.0], [0.0, 0.0]])
    to_class = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 0.0]])
    upward = np.array([[11.0, 0.0], [4.0, 0.0], [6.0, 1.0]])
    downward = np.array([[7.0, 2.0], [5.0, 0.0], [9.0, 0.0]])
    result = solid_mass_flow_residual(
        external, recycle, leaving, reaction, from_class, to_class, upward, downward
    )
    sources = external + recycle - leaving + reaction + from_class - to_class
    summed = result.residual_kg_s.sum(axis=0) + result.top_upward_out_kg_s + result.bottom_downward_out_kg_s
    assert np.array_equal(summed, sources.sum(axis=0))


def test_solid_bottom_feed_reaches_the_cell_above_on_the_same_size_class():
    # 外部加入只出现在最下层、粒级 1。写成该粒级的上行出流后，进入上一层的同一粒级。
    external = _zeros((3, 2))
    external[2, 1] = 4.0
    zeros = _zeros((3, 2))
    before = solid_mass_flow_residual(external, zeros, zeros, zeros, zeros, zeros, zeros, zeros)
    assert np.array_equal(before.residual_kg_s[:, 0], np.array([0.0, 0.0, 0.0]))
    assert np.array_equal(before.residual_kg_s[:, 1], np.array([0.0, 0.0, 4.0]))

    upward = _zeros((3, 2))
    upward[2, 1] = 4.0
    after = solid_mass_flow_residual(external, zeros, zeros, zeros, zeros, zeros, upward, zeros)
    assert np.array_equal(after.residual_kg_s[:, 1], np.array([0.0, 4.0, 0.0]))
    assert np.array_equal(after.residual_kg_s[:, 0], np.zeros(3))


def test_solid_residual_has_no_inventory_argument():
    names = set(inspect.signature(solid_mass_flow_residual).parameters)
    assert names.isdisjoint({"m", "K", "k", "holdup", "solid_mass_in_cell"})


def test_energy_hand_calculation_and_conservation():
    # cell 0：1 + 0 + 3 + 0 − 8 − 1 − 0 − 0 − 1 = −6
    #   其中 3 是下一层的上行能量流，最上层没有来自上方的下行入流。
    # cell 1：0 + 2 + 5 + 1 − 3 − 7 − 1 − 0 − 0 = −3
    # cell 2：4 + 0 + 0 + 7 − 5 − 6 − 0 − 2 − 0 = −2
    #   最下层没有来自更下方的上行入流；7 是上一层的下行能量流。
    external = np.array([1.0, 0.0, 4.0])
    recycle = np.array([0.0, 2.0, 0.0])
    upward = np.array([8.0, 3.0, 5.0])
    downward = np.array([1.0, 7.0, 6.0])
    leaving = np.array([0.0, 1.0, 0.0])
    wall = np.array([0.0, 0.0, 2.0])
    exchanger = np.array([1.0, 0.0, 0.0])
    result = energy_residual(external, recycle, upward, downward, leaving, wall, exchanger)
    assert np.array_equal(result.residual_w, np.array([-6.0, -3.0, -2.0]))
    assert result.top_upward_out_w == 8.0
    assert result.bottom_downward_out_w == 6.0

    sources = external + recycle - leaving - wall - exchanger
    assert result.residual_w.sum() + result.top_upward_out_w + result.bottom_downward_out_w == sources.sum()


def test_stream_enthalpy_flow_sums_phases_classes_and_required_s_term():
    # 气相：1×10 + 2×20 + 3×30 + 4×40 = 300
    # 第一种固体两个粒级：0.5×100 + 1.5×200 = 350
    # 第二种固体一个粒级：2×5 = 10
    # s 项由调用方给出，这里是 7 W。合计 667。
    gas_flow = np.array([[1.0, 2.0], [3.0, 4.0]])
    gas_enthalpy = np.array([[10.0, 20.0], [30.0, 40.0]])
    mass = [np.array([0.5, 1.5]), np.array([2.0])]
    enthalpy = [np.array([100.0, 200.0]), np.array([5.0])]
    assert stream_enthalpy_flow(gas_flow, gas_enthalpy, mass, enthalpy, 7.0) == 667.0


def test_stream_enthalpy_flow_requires_s_term_and_accepts_one_solid_vector():
    gas_flow = np.zeros((1, 1))
    gas_enthalpy = np.zeros((1, 1))
    with pytest.raises(TypeError):
        stream_enthalpy_flow(gas_flow, gas_enthalpy, [], [])
    # 2×4 + 3×5 + (−1) = 22
    assert stream_enthalpy_flow(gas_flow, gas_enthalpy, np.array([2.0, 3.0]), np.array([4.0, 5.0]), -1.0) == 22.0


def test_negative_flows_enter_the_residual_unchanged():
    value = np.array([[-4.0]])
    zeros = np.zeros((1, 1))
    gas = gas_mole_residual(value, zeros, zeros, zeros, zeros, zeros, zeros, zeros, zeros)
    assert gas.suspension_mol_s[0, 0] == -4.0
    solid = solid_mass_flow_residual(value, zeros, zeros, zeros, zeros, zeros, zeros, zeros)
    assert solid.residual_kg_s[0, 0] == -4.0
    energy = energy_residual(
        np.array([-4.0]),
        np.zeros(1),
        np.zeros(1),
        np.zeros(1),
        np.zeros(1),
        np.zeros(1),
        np.zeros(1),
    )
    assert energy.residual_w[0] == -4.0


def test_shape_mismatch_and_nan_raise():
    ok = np.zeros((2, 2))
    wrong = np.zeros((2, 1))
    with pytest.raises(ResidualInputError):
        gas_mole_residual(ok, ok, ok, ok, ok, ok, ok, wrong, ok)
    dirty = ok.copy()
    dirty[0, 1] = np.nan
    with pytest.raises(ResidualInputError):
        gas_mole_residual(dirty, ok, ok, ok, ok, ok, ok, ok, ok)

    with pytest.raises(ResidualInputError):
        solid_mass_flow_residual(ok, ok, ok, ok, ok, ok, wrong, ok)
    with pytest.raises(ResidualInputError):
        energy_residual(np.zeros((2, 1)), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2))
    wall = np.zeros(2)
    wall[1] = np.nan
    with pytest.raises(ResidualInputError):
        energy_residual(np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), wall, np.zeros(2))
    with pytest.raises(ResidualInputError):
        stream_enthalpy_flow(np.zeros((1, 1)), np.array([[np.nan]]), [], [], 0.0)
    with pytest.raises(ResidualInputError):
        stream_enthalpy_flow(np.zeros((1, 1)), np.zeros((1, 1)), [], [], np.nan)


def test_infinite_inputs_raise_and_finite_negatives_stay():
    # 0 × ∞ 会变成 NaN。无穷大在入口拒绝，不进入残差，也不截断有限负值。
    zeros_gas = np.zeros((1, 1))
    for bad in (np.inf, -np.inf):
        enthalpy = np.array([[bad]])
        with pytest.raises(ResidualInputError):
            stream_enthalpy_flow(zeros_gas, enthalpy, [], [], 0.0)
        with pytest.raises(ResidualInputError):
            stream_enthalpy_flow(zeros_gas, zeros_gas, [], [], bad)
        upward = np.array([0.0, bad])
        with pytest.raises(ResidualInputError):
            energy_residual(
                np.zeros(2), np.zeros(2), upward, np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)
            )
    assert stream_enthalpy_flow(zeros_gas, zeros_gas, [], [], -3.0) == -3.0


def test_conservation_holds_for_negative_and_single_cell_inputs():
    rng = np.random.default_rng(11)
    for n_cell in (1, 3):
        gas_shape = (n_cell, 2)
        flows = [rng.normal(size=gas_shape) for _ in range(9)]
        gas = gas_mole_residual(*flows)
        sources = flows[0] + flows[1] + flows[2] + flows[4] + flows[5] + flows[6]
        summed = gas.suspension_mol_s.sum(axis=0) + gas.bubble_mol_s.sum(axis=0)
        assert np.allclose(summed + gas.top_suspension_out_mol_s + gas.top_bubble_out_mol_s, sources.sum(axis=0))

        solid_shape = (n_cell, 3)
        solid_flows = [rng.normal(size=solid_shape) for _ in range(8)]
        solid = solid_mass_flow_residual(*solid_flows)
        solid_sources = (
            solid_flows[0] + solid_flows[1] - solid_flows[2] + solid_flows[3] + solid_flows[4] - solid_flows[5]
        )
        solid_sum = solid.residual_kg_s.sum(axis=0) + solid.top_upward_out_kg_s + solid.bottom_downward_out_kg_s
        assert np.allclose(solid_sum, solid_sources.sum(axis=0))

        energy_flows = [rng.normal(size=n_cell) for _ in range(7)]
        energy = energy_residual(*energy_flows)
        energy_sources = energy_flows[0] + energy_flows[1] - energy_flows[4] - energy_flows[5] - energy_flows[6]
        assert np.isclose(
            energy.residual_w.sum() + energy.top_upward_out_w + energy.bottom_downward_out_w,
            energy_sources.sum(),
        )
