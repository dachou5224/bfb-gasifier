"""衡算装配：打包、守恒、两状态、灰与粒级迁移。"""
from __future__ import annotations

import math

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.properties.data import GAS_SPECIES, atomic_weight, atoms, molar_mass_kg_mol
from hamel.solver import (
    CASE1_BED_CELLS, CASE1_FREEBOARD_CELLS, AssemblyError, FreeboardDomainError,
    TrialNotEvaluable, balance_index, evaluate_residual, geometric_index, pack_state,
    secondary_air_balance_index, unknown_slots, unpack_state,
)
from hamel.solver.frozen import FrozenQuantities
from hamel.solver.vector import ReactorState

_SPECIES = {name: index for index, name in enumerate(GAS_SPECIES)}
_TOL = 1.0e-8


def _zeros(shape):
    return np.zeros(shape, dtype=float)


def make_frozen(n_bed=2, n_free=2, n_class=2, **overrides) -> FrozenQuantities:
    """测试用冻结量。未改写的系数是小的正数，不是工况 1 的预计算。"""
    diameters = np.linspace(1.0e-4, 2.0e-4, n_class)
    species = len(GAS_SPECIES)
    feed = _zeros(species)
    feed[_SPECIES["N2"]] = 0.04
    feed[_SPECIES["O2"]] = 0.008
    feed[_SPECIES["H2O"]] = 0.008
    feed[_SPECIES["CO2"]] = 0.002
    feed[_SPECIES["H2"]] = 0.001
    feed[_SPECIES["CO"]] = 0.001
    fields = dict(
        diameter_m=diameters,
        bed_pressure_pa=np.full(n_bed, 5.0e5),
        freeboard_pressure_pa=np.full(n_free, 5.0e5),
        k_auf_per_s=np.full(n_bed, 0.05),
        m_dot_ab_kg_s=np.full(n_bed, 0.01),
        geometric_mass_kg=np.full(n_bed, 0.4),
        k_bd_per_s=np.full(n_bed, 0.2),
        v_bubble_exchange_m3=np.full(n_bed, 1.0e-3),
        v_suspension_m3=np.full(n_bed, 2.0e-3),
        v_bubble_m3=np.full(n_bed, 1.0e-3),
        v_freeboard_m3=np.full(n_free, 4.0e-3),
        primary_feed_mol_s=feed,
        feed_to_suspension=0.7,
        primary_temperature_k=458.15,
        secondary_feed_mol_s=feed * 0.5,
        secondary_temperature_k=458.15,
        secondary_balance_index=secondary_air_balance_index(n_free),
        pyro_bed_mol_s=_zeros((n_bed, species)),
        pyro_freeboard_mol_s=_zeros((n_free, species)),
        char_gen_bed_kg_s=_zeros((n_bed, n_class)),
        char_gen_free_kg_s=_zeros((n_free, n_class)),
        char_gen_hover_kg_s=_zeros((n_free, n_class)),
        ash_gen_bed_kg_s=_zeros((n_bed, n_class)),
        ash_gen_free_kg_s=_zeros((n_free, n_class)),
        ash_gen_hover_kg_s=_zeros((n_free, n_class)),
        fuel_bed_kg=_zeros((n_bed, n_class)),
        fuel_up_bed_kg_s=_zeros((n_bed, n_class)),
        fuel_down_bed_kg_s=_zeros((n_bed, n_class)),
        fuel_exit_bed_kg_s=_zeros((n_bed, n_class)),
        fuel_up_free_kg_s=_zeros((n_free, n_class)),
        fuel_down_free_kg_s=_zeros((n_free, n_class)),
        fuel_feed_kg_s=_zeros(n_class),
        fuel_formation_j_per_kg=-1.0e6,
        fuel_heat_capacity_j_per_kg_k=1500.0,
        fuel_feed_temperature_k=298.15,
        hover=np.zeros((n_free, n_class), dtype=bool),
        free_absent=np.zeros((n_free, n_class), dtype=bool),
        a_free_per_s=np.full((n_free, n_class), 0.1),
        a_hover_per_s=_zeros((n_free, n_class)),
        b_free_per_s=np.full((n_free, n_class), 0.02),
        b_hover_per_s=_zeros((n_free, n_class)),
        theta=_zeros((n_free, n_class)),
        ejection_to_free=np.ones(n_class),
        ejection_to_hover=_zeros(n_class),
        u_r_bed_m_s=np.full(n_bed, 0.4),
        u_r_freeboard_m_s=np.full((n_free, n_class), 0.8),
        u_r_hover_m_s=np.full((n_free, n_class), 0.5),
        rho_p_kg_m3=800.0,
        drain_char_kg_s=_zeros(n_class),
        drain_ash_kg_s=0.0,
        k_w_w_m2_k=2.0,
        bed_height_m=0.15,
        freeboard_height_m=0.2875,
        bed_diameter_m=0.15,
        freeboard_diameter_m=0.25,
        epsilon_mf=0.46,
        particle_tolerance_w=SolverConfig().particle_tolerance_w,
        particle_max_iterations=SolverConfig().particle_max_iterations,
        coupling_tolerance_k=SolverConfig().coupling_tolerance_k,
        coupling_max_iterations=SolverConfig().coupling_max_iterations,
    )
    fields["m_dot_ab_kg_s"][-1] = 0.0
    fields.update(overrides)
    if "m_dot_ab_kg_s" in overrides:
        fields["m_dot_ab_kg_s"] = np.asarray(overrides["m_dot_ab_kg_s"], dtype=float)
        fields["m_dot_ab_kg_s"][-1] = 0.0
    return FrozenQuantities(**fields)


def make_state(frozen: FrozenQuantities, temperature=1100.0, char=2.0e-3, ash=1.0e-4) -> ReactorState:
    flow = _zeros(len(GAS_SPECIES))
    flow[_SPECIES["N2"]] = 0.03
    flow[_SPECIES["O2"]] = 0.004
    flow[_SPECIES["H2O"]] = 0.006
    flow[_SPECIES["CO2"]] = 0.003
    flow[_SPECIES["H2"]] = 0.002
    flow[_SPECIES["CO"]] = 0.002
    flow[_SPECIES["CH4"]] = 0.0005
    bed_gas = np.zeros((frozen.n_bed, 2, len(GAS_SPECIES)))
    bed_gas[:] = flow
    free_gas = np.zeros((frozen.n_freeboard, len(GAS_SPECIES)))
    free_gas[:] = flow * 1.4
    char_bed = np.full((frozen.n_bed, frozen.n_class), char)
    char_free = np.full((frozen.n_freeboard, frozen.n_class), char)
    char_hover = np.zeros_like(char_free)
    ash_free = np.full_like(char_free, ash)
    ash_hover = np.zeros_like(char_free)
    char_free[frozen.free_absent] = 0.0
    ash_free[frozen.free_absent] = 0.0
    char_hover[frozen.hover] = char
    ash_hover[frozen.hover] = ash
    return ReactorState(
        bed_gas, free_gas, np.full(frozen.n_bed, temperature), np.full(frozen.n_freeboard, temperature),
        char_bed, char_free, char_hover, ash_free, ash_hover,
    )


def _evaluate(frozen, state):
    return evaluate_residual(pack_state(state, frozen), frozen)


def test_balance_and_secondary_air_index():
    assert balance_index(0, 8) == 7
    assert geometric_index(0, 8) == 7
    assert balance_index(geometric_index(3, 8), 8) == 3
    # 1.4875 m 是自由空间自下而上第 2 格，8 格时衡算下标为 6。
    assert secondary_air_balance_index(CASE1_FREEBOARD_CELLS) == 6
    assert geometric_index(secondary_air_balance_index(), CASE1_FREEBOARD_CELLS) == 1
    assert CASE1_BED_CELLS == 8


def test_pack_unpack_round_trip():
    frozen = make_frozen()
    state = make_state(frozen)
    packed = pack_state(state, frozen)
    restored = unpack_state(packed, frozen)
    again = pack_state(restored, frozen)
    assert again == pytest.approx(packed)
    assert restored.bed_temperature_k == pytest.approx(state.bed_temperature_k)
    assert restored.freeboard_char_free_kg == pytest.approx(state.freeboard_char_free_kg)
    names = [slot.name for slot in unknown_slots(frozen)]
    assert len(names) == len(set(names)) == len(packed)
    assert names[0].startswith("freeboard[b=0,geo=1].gas_outflow.gas.CO")
    assert any(slot.zone == "bed" and slot.kind == "temperature" and slot.balance_index == frozen.n_bed - 1 for slot in unknown_slots(frozen))
    assert not any(slot.kind == "ash" and slot.zone == "bed" for slot in unknown_slots(frozen))


def test_free_state_absent_is_omitted():
    hover = np.zeros((2, 2), dtype=bool)
    absent = np.zeros((2, 2), dtype=bool)
    # 顶格：粒级 2 只有悬停状态（床面抛出只进自由状态，故最下层不能缺自由状态）。
    hover[0, 1] = True
    absent[0, 1] = True
    # ADR 0048 起床面抛出只进自由状态，故自由状态不存在的格子只能在最下层以上没有抛入时出现。
    theta_up = np.array([[0.0, 1.0], [0.0, 0.0]])
    frozen = make_frozen(
        hover=hover, free_absent=absent, theta_up=theta_up,
        a_free_per_s=_zeros((2, 2)), b_free_per_s=_zeros((2, 2)),
        c_up_per_s=_zeros((2, 2)), c_down_per_s=_zeros((2, 2)),
        k_auf_per_s=_zeros(2),
    )
    slots = unknown_slots(frozen)
    assert not any(slot.balance_index == 0 and slot.size_class == 1 and slot.phase == "free" for slot in slots)
    assert any(slot.balance_index == 0 and slot.size_class == 1 and slot.phase == "hover" and slot.kind == "char" for slot in slots)
    state = make_state(frozen)
    result = _evaluate(frozen, state)
    assert result.char_gap_kg_s == pytest.approx(0.0, abs=1.0e-9)


def _closure_scale(result) -> float:
    """大速率相减后，绝对残差只保证到与最大项成比例的舍入。"""
    return 1.0e-8 * (1.0 + float(np.max(np.abs(result.values))))


def test_conservation_of_transport_reaction_and_elements():
    frozen = make_frozen()
    state = make_state(frozen)
    result = _evaluate(frozen, state)
    scale = _closure_scale(result)
    assert np.max(np.abs(result.gas_gap_mol_s)) <= scale
    assert abs(result.char_gap_kg_s) <= 1.0e-9
    assert abs(result.ash_gap_kg_s) <= 1.0e-9
    assert abs(result.energy_gap_w) <= 1.0e-4
    for element, gap in result.element_gap_kg_s.items():
        assert abs(gap) <= scale, element
    # 试探点本身不是解，分块残差不必为零。
    assert np.max(np.abs(result.values)) > 1.0e-6


def test_residual_uses_the_supplied_trial_point():
    frozen = make_frozen()
    first = make_state(frozen, temperature=1050.0)
    second = make_state(frozen, temperature=1120.0)
    first_values = _evaluate(frozen, first).values
    second_values = _evaluate(frozen, second).values
    repeat = _evaluate(frozen, first).values
    assert first_values == pytest.approx(repeat, rel=0.0, abs=1.0e-9)
    assert not np.allclose(first_values, second_values)


def test_pure_ash_steady_state():
    """m_C = 0、灰入流 0.001 kg/s、出流系数 0.2 s^-1 时 m_A = 0.005 kg。"""
    # 床区最上层（衡算 0）灰份额 × K = 0.001。m_inert = 0.4 − 0.02 = 0.38。
    # 出流系数放在接受抛出的最下格。
    frozen = make_frozen(
        n_class=1,
        k_auf_per_s=np.array([0.001 / 0.38, 0.0]),
        a_free_per_s=np.array([[0.0], [0.2]]),
        b_free_per_s=np.array([[0.0], [0.0]]),
        m_dot_ab_kg_s=np.array([0.0, 0.0]),
    )
    state = make_state(frozen, char=0.02, ash=0.0)
    state.bed_char_kg[:] = 0.02
    state.freeboard_char_free_kg[:] = 0.0
    state.freeboard_ash_free_kg[:] = 0.0
    state.freeboard_ash_free_kg[-1, 0] = 0.005
    result = _evaluate(frozen, state)
    ash_slots = [slot for slot in result.slots if slot.kind == "ash" and slot.balance_index == 1]
    assert len(ash_slots) == 1
    assert result.values[ash_slots[0].index] == pytest.approx(0.0, abs=1.0e-12)


def test_empty_ash_inventory_has_zero_flow_and_large_ratio_conserves():
    frozen = make_frozen(a_free_per_s=np.full((2, 2), 0.2), b_free_per_s=np.full((2, 2), 0.0), k_auf_per_s=np.zeros(2))
    empty = make_state(frozen, char=0.01, ash=0.0)
    empty.freeboard_char_free_kg[:] = 0.0
    empty.freeboard_ash_free_kg[:] = 0.0
    empty_result = _evaluate(frozen, empty)
    ash_values = [empty_result.values[slot.index] for slot in empty_result.slots if slot.kind == "ash"]
    assert np.allclose(ash_values, 0.0, atol=1.0e-12)
    wide = make_state(frozen, char=1.0e-4, ash=0.02)
    wide_result = _evaluate(frozen, wide)
    assert wide_result.ash_gap_kg_s == pytest.approx(0.0, abs=1.0e-9)
    assert wide_result.element_gap_kg_s["C"] == pytest.approx(0.0, abs=1.0e-8)


def test_burnout_sends_ash_to_smallest_class():
    frozen = make_frozen(
        k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), b_free_per_s=_zeros((2, 2)),
    )
    state = make_state(frozen, char=0.0, ash=0.0)
    state.bed_char_kg[:] = 0.01
    state.freeboard_char_free_kg[-1, 1] = 0.004
    state.freeboard_ash_free_kg[-1, 1] = 0.001
    result = _evaluate(frozen, state)
    cell = frozen.n_freeboard - 1
    ratio = 0.001 / 0.004
    particle = frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m**3 / 6.0
    expected_factor = particle[1] / (particle[1] - particle[0])
    consumption = result.char_consumption_free_kg_s[cell, 1]
    assert consumption > 0.0
    assert result.char_migration_free_kg_s[cell, 1] / consumption == pytest.approx(particle[0] / (particle[1] - particle[0]), rel=1.0e-12)
    ash_slots = {slot.size_class: slot for slot in result.slots if slot.kind == "ash" and slot.balance_index == cell and slot.phase == "free"}
    assert result.values[ash_slots[1].index] / consumption == pytest.approx(-ratio * expected_factor, rel=1.0e-9)
    assert result.values[ash_slots[0].index] + result.values[ash_slots[1].index] == pytest.approx(0.0, abs=1.0e-9)


def test_no_ash_shrinking_steady_state_matches_linear_solution():
    frozen = make_frozen(
        k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), b_free_per_s=_zeros((2, 2)),
    )
    probe = make_state(frozen, char=0.0, ash=0.0)
    probe.bed_char_kg[:] = 0.01
    probe.freeboard_char_free_kg[-1] = np.array([0.0, 1.0])
    probed = _evaluate(frozen, probe)
    cell = frozen.n_freeboard - 1
    rate_sum = probed.char_consumption_free_kg_s[cell, 1] + probed.char_migration_free_kg_s[cell, 1]
    small = make_state(frozen, char=0.0, ash=0.0)
    small.bed_char_kg[:] = 0.01
    small.freeboard_char_free_kg[-1] = np.array([1.0, 0.0])
    small_rate = _evaluate(frozen, small).char_consumption_free_kg_s[cell, 0]
    mass = 0.2
    steady = make_state(frozen, char=0.0, ash=0.0)
    steady.bed_char_kg[:] = 0.01
    steady.freeboard_char_free_kg[-1, 1] = mass
    steady.freeboard_char_free_kg[-1, 0] = probed.char_migration_free_kg_s[cell, 1] * mass / small_rate
    frozen_gen = make_frozen(
        k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), b_free_per_s=_zeros((2, 2)),
        char_gen_free_kg_s=np.array([[0.0, 0.0], [0.0, rate_sum * mass]]),
    )
    result = _evaluate(frozen_gen, steady)
    slots = [slot for slot in result.slots if slot.kind == "char" and slot.zone == "freeboard" and slot.balance_index == cell]
    assert np.allclose([result.values[slot.index] for slot in slots], 0.0, atol=1.0e-8)


def test_two_ash_ratios_and_two_states_sum_to_equation_2_6():
    hover = np.ones((2, 2), dtype=bool)
    frozen = make_frozen(
        hover=hover, k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), a_hover_per_s=_zeros((2, 2)),
        b_free_per_s=_zeros((2, 2)), b_hover_per_s=_zeros((2, 2)),
    )
    state = make_state(frozen, char=0.0, ash=0.0)
    state.bed_char_kg[:] = 0.01
    cell = frozen.n_freeboard - 1
    state.freeboard_char_free_kg[cell] = np.array([0.002, 0.003])
    state.freeboard_ash_free_kg[cell] = np.array([0.001, 0.0005])
    state.freeboard_char_hover_kg[cell] = np.array([0.002, 0.003])
    state.freeboard_ash_hover_kg[cell] = np.array([0.004, 0.006])
    result = _evaluate(frozen, state)
    particle = frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m**3 / 6.0
    for klass, char, ash in ((0, 0.002, 0.001), (1, 0.003, 0.0005)):
        assert (char + ash) / particle[klass] > 0.0
    free_rate = result.char_consumption_free_kg_s[cell] / state.freeboard_char_free_kg[cell]
    hover_rate = result.char_consumption_hover_kg_s[cell] / state.freeboard_char_hover_kg[cell]
    assert not np.allclose(free_rate, hover_rate)
    for klass in range(2):
        inflow = 0.0
        if klass == 0:
            inflow = result.char_migration_free_kg_s[cell, 1] + result.char_migration_hover_kg_s[cell, 1]
        expected = (
            -result.char_consumption_free_kg_s[cell, klass] - result.char_migration_free_kg_s[cell, klass]
            - result.char_consumption_hover_kg_s[cell, klass] - result.char_migration_hover_kg_s[cell, klass]
            + inflow
        )
        slots = [
            slot for slot in result.slots
            if slot.kind == "char" and slot.zone == "freeboard" and slot.balance_index == cell and slot.size_class == klass
        ]
        assert sum(result.values[slot.index] for slot in slots) == pytest.approx(expected, abs=1.0e-9)
    carbon = result.element_gap_kg_s["C"]
    assert carbon == pytest.approx(0.0, abs=1.0e-8)
    hand_n = (0.003 + 0.006) / particle[1]
    assert hand_n == pytest.approx((0.003 + 0.006) / (frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m[1] ** 3 / 6.0))


def test_hover_without_reaction_and_captured_fines_stop():
    hover = np.zeros((2, 1), dtype=bool)
    hover[-1, 0] = True
    frozen = make_frozen(n_class=1, hover=hover, c_up_per_s=np.array([[0.0], [0.05]]))
    state = make_state(frozen, char=0.002, ash=1.0e-4)
    state.bed_gas_mol_s[:] = 0.0
    state.freeboard_gas_mol_s[:] = 0.0
    state.bed_gas_mol_s[..., _SPECIES["N2"]] = 0.05
    state.freeboard_gas_mol_s[..., _SPECIES["N2"]] = 0.05
    with pytest.raises(FreeboardDomainError, match="没有反应"):
        _evaluate(frozen, state)

    captured = make_frozen(
        n_class=1, hover=hover, c_up_per_s=np.array([[0.0], [0.05]]),
        k_auf_per_s=np.array([0.05, 0.05]),
    )
    flowing = make_state(captured, char=0.02, ash=1.0e-4)
    with pytest.raises(FreeboardDomainError, match="最小粒级"):
        _evaluate(captured, flowing)


def test_solid_volume_above_voidage_limit_stops():
    frozen = make_frozen()
    state = make_state(frozen, char=50.0, ash=50.0)
    with pytest.raises(FreeboardDomainError, match="固体体积"):
        _evaluate(frozen, state)


def test_corella_outside_its_interval_is_not_evaluable():
    frozen = make_frozen()
    state = make_state(frozen, temperature=1300.0)
    with pytest.raises(Exception, match="600"):
        _evaluate(frozen, state)


def test_zero_gas_flow_is_not_evaluable():
    frozen = make_frozen()
    state = make_state(frozen)
    state.freeboard_gas_mol_s[0] = 0.0
    with pytest.raises(TrialNotEvaluable, match="总摩尔流"):
        _evaluate(frozen, state)


def test_fuel_mass_reduces_inert_and_missing_layout_is_rejected():
    frozen = make_frozen(fuel_bed_kg=np.full((2, 2), 0.05))
    state = make_state(frozen, char=0.01, ash=0.0)
    result = _evaluate(frozen, state)
    expected = 0.4 - 0.02 - 0.1
    assert result.m_inert_bed_kg == pytest.approx(expected)
    with pytest.raises(AssemblyError, match="自由状态不存在"):
        make_frozen(free_absent=np.ones((2, 2), dtype=bool), hover=np.zeros((2, 2), dtype=bool))


def test_case1_smoke_prints_block_magnitudes(capsys):
    """工况 1 的临时冻结量和一个试探 x。临时系数不是预计算结果。"""
    n_bed, n_free, n_class = 8, 8, 7
    diameters = np.array([5.06e-05, 1.527e-4, 3.466e-4, 5.928e-4, 8.385e-4, 1.3863e-3, 2.4885e-3])
    feed = _case1_primary()
    secondary = _zeros(len(GAS_SPECIES))
    secondary[_SPECIES["O2"]] = 0.373030744673468 * 0.21
    secondary[_SPECIES["N2"]] = 0.373030744673468 * 0.79
    frozen = make_frozen(
        n_bed=n_bed, n_free=n_free, n_class=n_class, diameter_m=diameters,
        primary_feed_mol_s=feed, secondary_feed_mol_s=secondary,
        geometric_mass_kg=np.full(n_bed, 3.9629153924566762 / n_bed),
        k_auf_per_s=np.full(n_bed, 0.02),
        rho_p_kg_m3=700.0, k_w_w_m2_k=5.0,
    )
    state = make_state(frozen, temperature=1100.0, char=0.01, ash=0.0)
    state.freeboard_ash_free_kg[:] = 1.0e-4
    result = _evaluate(frozen, state)
    blocks = {kind: [] for kind in ("gas_outflow", "temperature", "char", "ash")}
    for slot in result.slots:
        blocks[slot.kind].append(result.values[slot.index])
    print("工况 1 临时冻结量，试探温度 1100 K")
    for kind, values in blocks.items():
        array = np.asarray(values)
        print(f"{kind}: n={array.size} max|.|={np.max(np.abs(array)):.6e}")
    print(
        "闭合",
        f"gas={np.max(np.abs(result.gas_gap_mol_s)):.3e}",
        f"char={result.char_gap_kg_s:.3e}",
        f"ash={result.ash_gap_kg_s:.3e}",
        f"energy={result.energy_gap_w:.3e}",
        " ".join(f"{name}={gap:.3e}" for name, gap in result.element_gap_kg_s.items()),
    )
    captured = capsys.readouterr().out
    assert "gas_outflow" in captured
    scale = _closure_scale(result)
    assert np.max(np.abs(result.gas_gap_mol_s)) <= scale
    assert abs(result.char_gap_kg_s) <= 1.0e-6
    assert abs(result.ash_gap_kg_s) <= 1.0e-6
    assert abs(result.energy_gap_w) <= max(1.0e-2, scale)
    for gap in result.element_gap_kg_s.values():
        assert abs(gap) <= scale


def _case1_primary():
    air = 0.6419598861822472
    purge = 0.07683689757393693
    steam = 0.001888888888888889 / molar_mass_kg_mol("H2O")
    feed = _zeros(len(GAS_SPECIES))
    feed[_SPECIES["O2"]] = air * 0.21
    feed[_SPECIES["N2"]] = air * 0.79 + purge
    feed[_SPECIES["H2O"]] = steam
    return feed


def test_pure_ash_does_not_react_as_char():
    """含 O2、H2O、H2 的纯灰不产生异相碳汇，灰仍按出流系数离开。"""
    frozen = make_frozen(
        n_class=1, k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=np.array([[0.0], [0.2]]), b_free_per_s=np.zeros((2, 1)),
    )
    state = make_state(frozen, char=0.01, ash=0.0)
    state.freeboard_char_free_kg[:] = 0.0
    state.freeboard_ash_free_kg[:] = 0.0
    state.freeboard_ash_free_kg[-1, 0] = 0.005
    result = _evaluate(frozen, state)
    assert np.allclose(result.char_consumption_free_kg_s, 0.0, atol=1.0e-15)
    assert np.allclose(result.char_migration_free_kg_s, 0.0, atol=1.0e-15)
    ash_slot = next(slot for slot in result.slots if slot.kind == "ash" and slot.balance_index == 1)
    assert result.values[ash_slot.index] == pytest.approx(-0.001, abs=1.0e-12)
    assert abs(result.element_gap_kg_s["C"]) <= _closure_scale(result)
    assert abs(result.ash_gap_kg_s) <= 1.0e-12


def test_pyrolysis_ash_into_smallest_hover_stops():
    hover = np.zeros((2, 1), dtype=bool)
    hover[-1, 0] = True
    ash_gen = np.zeros((2, 1))
    ash_gen[-1, 0] = 0.001
    frozen = make_frozen(n_class=1, hover=hover, k_auf_per_s=np.zeros(2), ash_gen_hover_kg_s=ash_gen)
    state = make_state(frozen, char=0.002, ash=0.0)
    with pytest.raises(FreeboardDomainError, match="最小粒级"):
        _evaluate(frozen, state)


def test_large_class_ash_transports_without_char():
    """ADR 0044 第 4 项：取消 m_C = 0 的切换，k > 1 的灰按同一线性式迁移与释放，不向 k = 1 改记。"""
    frozen = make_frozen(
        k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=np.array([[0.0, 0.0], [0.2, 0.2]]), b_free_per_s=_zeros((2, 2)),
    )
    cell = frozen.n_freeboard - 1
    particle = frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m**3 / 6.0
    q_two = particle[0] / (particle[1] - particle[0])
    probe = make_state(frozen, char=0.01, ash=0.0)
    probe.freeboard_char_free_kg[:] = 0.0
    probe.freeboard_char_free_kg[cell, 1] = 0.002
    probe.freeboard_ash_free_kg[cell, 1] = 0.0005
    alpha = _evaluate(frozen, probe).char_consumption_free_kg_s[cell, 1] / 0.002
    state = make_state(frozen, char=0.01, ash=0.0)
    state.freeboard_char_free_kg[:] = 0.0
    state.freeboard_ash_free_kg[:] = 0.0
    state.freeboard_ash_free_kg[cell, 1] = 0.004
    result = _evaluate(frozen, state)
    slot = next(
        item for item in result.slots
        if item.kind == "ash" and item.balance_index == cell and item.size_class == 1
    )
    expected = -0.2 * 0.004 - alpha * 0.004 * q_two - alpha * 0.004
    assert result.values[slot.index] == pytest.approx(expected, rel=1.0e-12)
    doubled = make_state(frozen, char=0.01, ash=0.0)
    doubled.freeboard_char_free_kg[:] = 0.0
    doubled.freeboard_ash_free_kg[:] = 0.0
    doubled.freeboard_ash_free_kg[cell, 1] = 0.008
    again = _evaluate(frozen, doubled)
    assert again.values[slot.index] == pytest.approx(2.0 * expected, rel=1.0e-12)


def test_char_and_ash_sinks_are_linear_in_the_inventories():
    """ADR 0044 第 3、5 项：消耗 α_s·M_C、迁移与灰释放都对存量线性，M_C = 0 处连续。"""
    frozen = make_frozen(
        k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), b_free_per_s=_zeros((2, 2)),
    )
    cell = frozen.n_freeboard - 1
    state = make_state(frozen, char=0.0, ash=0.0)
    state.bed_char_kg[:] = 0.01
    state.freeboard_char_free_kg[cell, 1] = 0.002
    state.freeboard_ash_free_kg[cell, 1] = 0.0005
    base = _evaluate(frozen, state)
    alpha = base.char_consumption_free_kg_s[cell, 1] / 0.002
    particle = frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m**3 / 6.0
    step = particle[0] / (particle[1] - particle[0])
    assert base.char_migration_free_kg_s[cell, 1] == pytest.approx(alpha * 0.002 * step, rel=1.0e-12)
    rows = {
        (item.size_class, item.kind): item.index
        for item in base.slots
        if item.balance_index == cell and item.zone == "freeboard"
    }
    # k = 2 的灰行：−α·M_A·q_2 − α·M_A；同一份流进入 k = 1 的灰行。
    assert base.values[rows[(1, "ash")]] == pytest.approx(-alpha * 0.0005 * (1.0 + step), rel=1.0e-12)
    assert base.values[rows[(0, "ash")]] == pytest.approx(alpha * 0.0005 * (1.0 + step), rel=1.0e-12)
    doubled = make_state(frozen, char=0.0, ash=0.0)
    doubled.bed_char_kg[:] = 0.01
    doubled.freeboard_char_free_kg[cell, 1] = 0.004
    doubled.freeboard_ash_free_kg[cell, 1] = 0.001
    twice = _evaluate(frozen, doubled)
    assert twice.char_consumption_free_kg_s[cell, 1] == pytest.approx(2.0 * base.char_consumption_free_kg_s[cell, 1])
    assert twice.values[rows[(1, "ash")]] == pytest.approx(2.0 * base.values[rows[(1, "ash")]])
    # 只有灰、没有焦炭：焦炭的消耗与迁移为 0；灰按同一线性式迁移与释放，不切换。
    ash_only = make_state(frozen, char=0.0, ash=0.0)
    ash_only.bed_char_kg[:] = 0.01
    ash_only.freeboard_ash_free_kg[cell, 1] = 0.001
    only = _evaluate(frozen, ash_only)
    assert only.char_consumption_free_kg_s[cell, 1] == 0.0
    assert only.char_migration_free_kg_s[cell, 1] == 0.0
    assert only.values[rows[(1, "ash")]] == pytest.approx(-alpha * 0.001 * (1.0 + step), rel=1.0e-12)
    assert only.values[rows[(0, "ash")]] == pytest.approx(alpha * 0.001 * (1.0 + step), rel=1.0e-12)


def test_free_and_hover_states_have_different_unit_char_rates():
    """ADR 0044 第 2 项：两个状态的 u_r 不同，α_s 一般不同。"""
    hover = np.ones((2, 2), dtype=bool)
    frozen = make_frozen(
        hover=hover, k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), a_hover_per_s=_zeros((2, 2)),
        b_free_per_s=_zeros((2, 2)), b_hover_per_s=_zeros((2, 2)),
    )
    state = make_state(frozen, char=0.0, ash=0.0)
    state.bed_char_kg[:] = 0.01
    cell = frozen.n_freeboard - 1
    state.freeboard_char_free_kg[cell, 1] = 0.002
    state.freeboard_char_hover_kg[cell, 1] = 0.002
    result = _evaluate(frozen, state)
    free_rate = result.char_consumption_free_kg_s[cell, 1] / 0.002
    hover_rate = result.char_consumption_hover_kg_s[cell, 1] / 0.002
    assert free_rate > 0.0 and hover_rate > 0.0
    assert free_rate != pytest.approx(hover_rate, rel=1.0e-6)


def test_single_cell_steady_state_with_migration_and_ash():
    """无输运时 k = 2 的稳态 m_C = G/(α_s(1+q_2))，m_A 由灰行给出。"""
    cell = 1
    inflow = 4.0e-5
    generation = np.zeros((2, 2))
    generation[cell, 1] = inflow
    frozen = make_frozen(
        k_auf_per_s=np.zeros(2), m_dot_ab_kg_s=np.zeros(2),
        a_free_per_s=_zeros((2, 2)), b_free_per_s=_zeros((2, 2)),
        char_gen_free_kg_s=generation,
    )
    probe = make_state(frozen, char=0.0, ash=0.0)
    probe.bed_char_kg[:] = 0.01
    probe.freeboard_char_free_kg[cell, 1] = 1.0
    probe.freeboard_ash_free_kg[cell, 1] = 0.0
    probed = _evaluate(frozen, probe)
    alpha = probed.char_consumption_free_kg_s[cell, 1]
    particle = frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m**3 / 6.0
    q_two = particle[0] / (particle[1] - particle[0])
    steady = make_state(frozen, char=0.0, ash=0.0)
    steady.bed_char_kg[:] = 0.01
    steady.freeboard_char_free_kg[cell, 1] = inflow / (alpha * (1.0 + q_two))
    result = _evaluate(frozen, steady)
    slot = next(
        item for item in result.slots
        if item.kind == "char" and item.balance_index == cell and item.size_class == 1 and item.phase == "free"
    )
    assert result.values[slot.index] == pytest.approx(0.0, abs=1.0e-12)


def test_element_weights_are_available_for_the_closure():
    assert atoms("CO")["C"] == 1
    assert atomic_weight("C") > 12.0


def test_public_equation_docstrings_cite_the_source():
    from hamel.solver.frozen import FrozenQuantities
    from hamel.solver.residual import evaluate_residual, wall_area_m2
    from hamel.solver.vector import balance_index, pack_state, unknown_slots

    for function in (
        evaluate_residual, wall_area_m2, balance_index, pack_state, unknown_slots,
        FrozenQuantities.n_bed.fget, FrozenQuantities.n_class.fget,
    ):
        assert "Source: Hamel Eq." in (function.__doc__ or ""), function.__name__


def test_theta_down_migration_carries_source_enthalpy_between_cells():
    """θ↓ 的焦炭与灰迁移流在两格能量方程中按源格温度分别扣除与接收。ADR 0048 第 2 项。

    只改变非零 θ↓，让源格与目标格的自由空间温度不同：源格（top，1150 K）应减少
    E = m_C·h_C(T_src) + m_A·h_A(T_src)，下方目标格应在同一 E 上增加，且两格之差为零。
    全炉求和会互相抵消，不能替代这里的逐格核对。
    """
    from hamel.properties import ash_specific_heat_j_per_kg_k, char_specific_enthalpy_j_per_kg

    hover = np.zeros((2, 3), dtype=bool)
    hover[:, 1] = True
    base = make_frozen(n_class=3, hover=hover, theta_up=_zeros((2, 3)), theta_down=_zeros((2, 3)))
    theta_down = _zeros((2, 3))
    theta_down[0, 1] = 1.0
    moved = make_frozen(n_class=3, hover=hover, theta_up=_zeros((2, 3)), theta_down=theta_down)
    state = make_state(base, char=2.0e-3, ash=1.0e-4)
    state.freeboard_temperature_k[:] = [1150.0, 1050.0]
    before = _evaluate(base, state)
    after = _evaluate(moved, state)

    def slot(kind, phase, balance, klass, result):
        for item in result.slots:
            if (item.kind, item.zone, item.balance_index, item.size_class) == (kind, "freeboard", balance, klass):
                if phase is not None and item.phase != phase:
                    continue
                return float(result.values[item.index])
        raise AssertionError((kind, phase, balance, klass))

    char_mass = slot("char", "hover", 1, 1, after) - slot("char", "hover", 1, 1, before)
    ash_mass = slot("ash", "hover", 1, 1, after) - slot("ash", "hover", 1, 1, before)
    assert char_mass > 0.0 and ash_mass > 0.0
    source_k = 1150.0
    expected = char_mass * char_specific_enthalpy_j_per_kg(source_k) + ash_mass * ash_specific_heat_j_per_kg_k() * (
        source_k - 298.15
    )
    source_delta = slot("temperature", None, 0, -1, after) - slot("temperature", None, 0, -1, before)
    target_delta = slot("temperature", None, 1, -1, after) - slot("temperature", None, 1, -1, before)
    assert source_delta == pytest.approx(-expected, rel=1.0e-10)
    assert target_delta == pytest.approx(expected, rel=1.0e-10)
    assert source_delta + target_delta == pytest.approx(0.0, abs=1.0e-9)

    # 若误用目标格温度（1050 K）计焓，源格扣除量会差出这一项。
    wrong = char_mass * char_specific_enthalpy_j_per_kg(1050.0) + ash_mass * ash_specific_heat_j_per_kg_k() * (
        1050.0 - 298.15
    )
    assert abs(source_delta + wrong) > 1.0e-3 * abs(expected)
