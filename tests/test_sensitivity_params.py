"""T25d：默认值逐位不变，改值进入计算，入口定位，悬停改送，诊断。"""
from __future__ import annotations

import math
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.hydrodynamics.bed import BedConditions, BedProfile
from hamel.properties.data import GAS_SPECIES
from hamel.pyrolysis.release import FuelTransport, feed_cell_index
from hamel.solver.campaign import config_for, overlay, run_specs
from hamel.solver.frozen import AssemblyError
from hamel.solver.precalc import EPSILON_MF, PHI_S, _bounds, _exchange_diffusivity, _load
from hamel.solver.residual import _bed_ash_by_class, evaluate_residual, wall_area_m2
from hamel.solver.sensitivity import (
    bubble_balance_ratio, edge_diameters, interface_balance_index, pack_reported, redirect_ledger,
)
from hamel.solver.vector import pack_state, unknown_slots, unpack_state
from test_assembly import make_frozen, make_state
from test_driver import _CASE
from test_precalc import _run, _transport

_SPECIES = {name: index for index, name in enumerate(GAS_SPECIES)}


def _pair():
    frozen = make_frozen()
    state = make_state(frozen)
    state.bed_gas_mol_s[:, 1, _SPECIES["C6H6"]] = 1.0e-4
    state.freeboard_gas_mol_s[:, _SPECIES["C6H6"]] = 1.0e-4
    return frozen, state


def test_defaults_match_the_frozen_residual_bit_for_bit():
    frozen, state = _pair()
    base = evaluate_residual(pack_state(state, frozen), frozen)
    explicit = make_frozen(
        particle_emissivity=0.85, hover_closure="A",
        serio_frequency_per_s=5.42e4, serio_activation_j_mol=100.5e3,
    )
    other = evaluate_residual(pack_state(state, explicit), explicit)
    assert np.array_equal(base.values, other.values)
    changed = make_frozen(particle_emissivity=0.79)
    shifted = evaluate_residual(pack_state(state, changed), changed)
    assert not np.array_equal(base.values, shifted.values)
    second = make_frozen(serio_frequency_per_s=9.96e2, serio_activation_j_mol=64.5e3)
    serio = evaluate_residual(pack_state(state, second), second)
    assert not np.array_equal(base.values, serio.values)


def test_case_file_voidage_matches_the_module_literal():
    case = _load(_CASE, "1")
    assert case.epsilon_mf == EPSILON_MF == 0.46
    assert case.phi_s == PHI_S == 0.8
    assert case.solver.particle_emissivity == 0.85
    assert case.solver.interphase_diffusivity == "hamel_5_20"


def test_fuller_changes_only_the_exchange_coefficient():
    conditions = BedConditions(1.0e-3, 0.46, 1600.0, 1.2, 2.8e-5, 2.8e-5 / 1.2, 101300.0, 2.0e-5, 0.15, "technical", "B", True)
    profile = BedProfile(
        0.2, np.array([0.2]), np.array([1.2]), np.array([0.04]), np.array([0.7]),
        np.array([3.0]), np.array([0.4]), np.array([0.01]),
    )
    flow = SimpleNamespace(bed_temperature=np.array([1100.0]), pressure=101300.0, u0=np.array([0.8]))
    case = SimpleNamespace(solver=SolverConfig())
    assert _exchange_diffusivity(case, flow, profile, conditions) is profile
    fuller = _exchange_diffusivity(
        SimpleNamespace(solver=replace(SolverConfig(), interphase_diffusivity="fuller_o2_n2")),
        flow, profile, conditions,
    )
    assert fuller.k_bd[0] != profile.k_bd[0]
    assert fuller.d_b[0] == profile.d_b[0]
    assert fuller.u_b[0] == profile.u_b[0]


def test_secondary_air_must_land_on_a_cell_boundary_and_fuel_stays_inside_a_cell():
    assert interface_balance_index(1.4875, 1.2, 8, 0.2875) == 6
    assert interface_balance_index(1.4875, 1.2, 16, 0.14375) == 13
    with pytest.raises(AssemblyError):
        interface_balance_index(1.50, 1.2, 16, 0.14375)
    edges = edge_diameters(8, 0.2875, 0.15, 0.25)
    assert edges[0] == 0.15 and edges[1] == 0.25 and np.all(edges[2:] == 0.25)
    fine = edge_diameters(16, 0.14375, 0.15, 0.25)
    assert fine[0] == 0.15
    assert fine[1] == pytest.approx(0.20)
    assert fine[2] == pytest.approx(0.25)
    radius_bed, radius_free = 0.075, 0.125
    slant = math.sqrt(0.2875**2 + (radius_free - radius_bed) ** 2)
    assert wall_area_m2("freeboard", 7, 8, 0.2875, 0.15, 0.25) == math.pi * (radius_bed + radius_free) * slant
    assert wall_area_m2("freeboard", 0, 8, 0.2875, 0.15, 0.25) == math.pi * 0.25 * 0.2875
    case = SimpleNamespace(n_bed=16, n_free=16, bed_cell_m=0.075, free_cell_m=0.14375)
    bottoms, tops = _bounds(case)
    index = feed_cell_index(0.12, bottoms, tops)
    assert bottoms[index] <= 0.12 <= tops[index]
    assert not math.isclose(bottoms[index], 1.4875)


def _closure_frozen(closure: str):
    n_free, n_class = 2, 2
    hover = np.zeros((n_free, n_class), dtype=bool)
    if closure == "A":
        hover[:, 1] = True
    theta = np.zeros((n_free, n_class))
    theta[:, 1] = 1.0
    a_hover = np.zeros((n_free, n_class))
    b_hover = np.zeros((n_free, n_class))
    a_hover[1, 1] = 0.2
    b_hover[0, 1] = 0.15
    fuel_up = np.zeros((n_free, n_class))
    fuel_down = np.zeros((n_free, n_class))
    fuel_up[1, 1] = 1.0e-4
    fuel_down[0, 1] = 2.0e-4
    ejection_to_free = np.ones(n_class)
    ejection_to_hover = np.zeros(n_class)
    ejection_to_free[1] = 0.75
    ejection_to_hover[1] = 0.25
    return make_frozen(
        n_bed=1, n_free=n_free, n_class=n_class, hover=hover, hover_closure=closure,
        a_hover_per_s=a_hover, b_hover_per_s=b_hover, theta=theta,
        ejection_to_free=ejection_to_free, ejection_to_hover=ejection_to_hover,
        fuel_up_free_kg_s=fuel_up, fuel_down_free_kg_s=fuel_down,
        fuel_up_bed_kg_s=np.full((1, n_class), 1.0e-4),
        char_gen_hover_kg_s=np.zeros((n_free, n_class)), ash_gen_hover_kg_s=np.zeros((n_free, n_class)),
    )


def test_redirect_ledger_cancels_item_by_item_and_the_gaps_stay_closed():
    frozen_a = _closure_frozen("A")
    state = make_state(frozen_a, char=2.0e-3, ash=1.0e-4)
    state.freeboard_char_hover_kg[:] = 0.0
    state.freeboard_ash_hover_kg[:] = 0.0
    zeros = np.zeros_like(state.freeboard_char_free_kg)
    inert = frozen_a.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen_a.fuel_bed_kg.sum(axis=1)
    share = _bed_ash_by_class(state.bed_char_kg, inert)
    items = redirect_ledger(frozen_a, state, zeros, zeros, zeros, zeros, share, zeros, zeros)
    assert items == ()
    for closure in ("B", "C"):
        frozen = _closure_frozen(closure)
        moved = redirect_ledger(frozen, state, zeros, zeros, zeros, zeros, share, zeros, zeros)
        assert moved
        assert sum(item["char"] for item in moved) > 0.0
        assert sum(item["enthalpy"] for item in moved) != 0.0
        counted = {}
        for item in moved:
            counted[item["kind"]] = counted.get(item["kind"], 0.0) + item["char"]
        assert counted["a_t"] == pytest.approx(float((frozen.a_hover_per_s[1] * state.freeboard_char_free_kg[1]).sum()))
        outflow = (frozen.a_free_per_s + frozen.a_hover_per_s + frozen.b_free_per_s + frozen.b_hover_per_s) * state.freeboard_char_free_kg
        assert float(outflow[1].sum()) == pytest.approx(float(((frozen.a_free_per_s[1] + frozen.a_hover_per_s[1] + frozen.b_free_per_s[1] + frozen.b_hover_per_s[1]) * state.freeboard_char_free_kg[1]).sum()))
        result = evaluate_residual(pack_state(state, frozen), frozen)
        assert abs(result.char_gap_kg_s) <= 1.0e-8
        assert abs(result.ash_gap_kg_s) <= 1.0e-8
        assert abs(result.energy_gap_w) <= 1.0e-2
        assert result.redirect_char_kg_s == pytest.approx(sum(item["char"] for item in moved), rel=1.0e-9, abs=1.0e-12)
        assert result.redirect_ash_kg_s == pytest.approx(sum(item["ash"] for item in moved), rel=1.0e-9, abs=1.0e-12)
        assert result.redirect_fuel_kg_s == pytest.approx(sum(item["fuel"] for item in moved), rel=1.0e-9, abs=1.0e-12)
        assert result.redirect_enthalpy_w == pytest.approx(sum(item["enthalpy"] for item in moved), rel=1.0e-9, abs=1.0e-9)
    state.freeboard_char_hover_kg[:, 1] = 2.0e-3
    state.freeboard_ash_hover_kg[:, 1] = 1.0e-4
    same = evaluate_residual(pack_state(state, frozen_a), frozen_a)
    assert abs(same.char_gap_kg_s) <= 1.0e-8
    assert abs(same.energy_gap_w) <= 1.0e-2
    assert same.redirect_char_kg_s == 0.0


def test_release_closures_solve_separately_and_keep_the_split():
    captured = _transport(False)
    base = _run(np.array([2.0]), np.array([3.0]), captured)
    bed = captured.bed_count
    for closure in ("B", "C"):
        transport = replace(captured, hover=np.zeros_like(captured.hover), closure=closure)
        solved = _run(np.array([2.0]), np.array([3.0]), transport)
        assert solved.char_hover_kg_s.sum() == 0.0
        assert abs(solved.mass_gap_kg_s).max() <= 1.0e-6
        if closure == "B":
            assert solved.undried_kg[bed].sum() + solved.dried_kg[bed].sum() != pytest.approx(
                base.undried_kg[bed].sum() + base.dried_kg[bed].sum()
            )


def test_char_solves_separately_and_report_only_the_difference_from_a():
    """异相消耗在固定气体上几乎不随存量变化。用 A 在参考存量上的焦炭格残差配平生成量，
    A 就是焦炭根；B、C 用同一床区与自由空间生成量各自求解。悬停生成只留给 A。
    与 A 的存量差只核对“不等于改送流”，不把差值当成验收等式。"""
    generation = _char_generation_that_roots_a()
    results = {}
    for closure in ("A", "B", "C"):
        frozen = replace(
            _closure_frozen(closure),
            char_gen_bed_kg_s=generation["bed"],
            char_gen_free_kg_s=generation["free"],
            char_gen_hover_kg_s=generation["hover"] if closure == "A" else np.zeros_like(generation["hover"]),
        )
        state = _closure_state(frozen, closure)
        solved = _solve_char(frozen, state)
        residual = evaluate_residual(pack_state(solved, frozen), frozen)
        char_indexes = [slot.index for slot in residual.slots if slot.kind == "char"]
        assert np.linalg.norm(residual.values[char_indexes]) <= 1.0e-8
        assert abs(residual.char_gap_kg_s) <= 1.0e-6
        assert abs(residual.ash_gap_kg_s) <= 1.0e-6
        assert abs(residual.energy_gap_w) <= 1.0e-2
        for element_gap in residual.element_gap_kg_s.values():
            assert abs(element_gap) <= 1.0e-5
        hover_slots = [slot for slot in unknown_slots(frozen) if slot.phase == "hover"]
        assert (closure == "A") == (len(hover_slots) > 0)
        results[closure] = (frozen, solved, residual)
    frozen_b, state_b, residual_b = results["B"]
    _frozen_a, state_a, residual_a = results["A"]
    _frozen_c, state_c, residual_c = results["C"]
    assert residual_a.redirect_char_kg_s == 0.0
    assert residual_b.redirect_char_kg_s > 0.0
    assert residual_c.redirect_char_kg_s > 0.0
    assert not np.allclose(state_b.freeboard_char_free_kg, state_a.freeboard_char_free_kg)
    assert not np.allclose(state_c.freeboard_char_free_kg, state_a.freeboard_char_free_kg)
    difference = float((state_b.freeboard_char_free_kg - state_a.freeboard_char_free_kg).sum())
    assert difference != pytest.approx(residual_b.redirect_char_kg_s)
    assert _elutriation_from_gap(frozen_b, residual_b) == pytest.approx(_free_exit(frozen_b, state_b))
    assert _elutriation_from_gap(_frozen_c, residual_c) == pytest.approx(
        _free_exit(_frozen_c, state_c) + residual_c.redirect_char_kg_s
    )


def _closure_state(frozen, closure: str):
    state = make_state(frozen, char=2.0e-3, ash=1.0e-4)
    if closure == "A":
        state.freeboard_char_hover_kg[:] = 2.0e-3
        state.freeboard_ash_hover_kg[:] = 1.0e-4
    else:
        state.freeboard_char_hover_kg[:] = 0.0
        state.freeboard_ash_hover_kg[:] = 0.0
    return state


def _char_generation_that_roots_a():
    frozen = _closure_frozen("A")
    residual = evaluate_residual(pack_state(_closure_state(frozen, "A"), frozen), frozen)
    bed = np.zeros_like(frozen.char_gen_bed_kg_s)
    free = np.zeros_like(frozen.char_gen_free_kg_s)
    hover = np.zeros_like(frozen.char_gen_hover_kg_s)
    for slot in residual.slots:
        if slot.kind != "char":
            continue
        amount = -float(residual.values[slot.index])
        if slot.zone == "bed":
            bed[slot.balance_index, slot.size_class] = amount
        elif slot.phase == "free":
            free[slot.balance_index, slot.size_class] = amount
        else:
            hover[slot.balance_index, slot.size_class] = amount
    return {"bed": bed, "free": free, "hover": hover}


def _free_exit(frozen, state) -> float:
    return float(((frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * state.freeboard_char_free_kg[0]).sum())


def _elutriation_from_gap(frozen, residual) -> float:
    """由已闭合的全炉焦炭账反推出口流。格残差之和里不含出口。"""
    char_sum = float(sum(residual.values[slot.index] for slot in residual.slots if slot.kind == "char"))
    generated = float(frozen.char_gen_bed_kg_s.sum() + frozen.char_gen_free_kg_s.sum() + frozen.char_gen_hover_kg_s.sum())
    consumed = float(
        residual.char_consumption_bed_kg_s.sum()
        + residual.char_consumption_free_kg_s.sum()
        + residual.char_consumption_hover_kg_s.sum()
    )
    return float(residual.char_gap_kg_s - char_sum + generated - consumed - frozen.drain_char_kg_s.sum())


def _solve_char(frozen, state, steps=8):
    values = pack_state(state, frozen)
    indexes = [slot.index for slot in unknown_slots(frozen) if slot.kind == "char"]

    def residual_at(trial):
        return evaluate_residual(trial, frozen)

    for _ in range(steps):
        current = residual_at(values)
        rhs = current.values[indexes]
        if np.linalg.norm(rhs) <= 1.0e-8:
            return unpack_state(values, frozen)
        jacobian = np.zeros((len(indexes), len(indexes)))
        for column, index in enumerate(indexes):
            step = 1.0e-7 * max(1.0, abs(values[index]))
            bumped = values.copy()
            bumped[index] += step
            jacobian[:, column] = (residual_at(bumped).values[indexes] - rhs) / step
        delta = np.linalg.solve(jacobian, -rhs)
        accepted = False
        for scale in (1.0, 0.5, 0.25):
            trial = values.copy()
            trial[indexes] = np.maximum(values[indexes] + scale * delta, 1.0e-9)
            if np.linalg.norm(residual_at(trial).values[indexes]) < np.linalg.norm(rhs):
                values = trial
                accepted = True
                break
        if not accepted:
            break
    residual = residual_at(values)
    char_norm = np.linalg.norm(residual.values[indexes])
    raise AssertionError(f"焦炭未解到格残差 1e-8 以内，当前范数 {char_norm}")


def test_diagnostics_publish_existing_intermediates():
    profile = SimpleNamespace(
        epsilon_b=np.array([0.2, 0.1]), u_b=np.array([1.5, 1.2]), d_b=np.array([0.2, 0.05]),
    )
    flow = SimpleNamespace(
        u0=np.array([0.8, 0.7]), free_mol_s=np.array([[0.02, 0.01]]), free_temperature=np.array([1000.0]),
        pressure=101300.0,
    )
    released = SimpleNamespace(
        mass_gap_kg_s=np.array([1.0e-6, -2.0e-6, 0.0]),
        element_gap_kg_s={"C": np.array([1.0e-7, 0.0, 3.0e-7])},
    )
    case = SimpleNamespace(bed_diameter_m=0.15, free_diameter_m=0.25, n_free=2)
    reported = pack_reported(profile, flow, case, released)
    assert reported["cells_db_exceeds_d"] == 1
    assert reported["kurkela_bed_m_s"] == 0.74
    assert reported["kurkela_freeboard_m_s"] == 0.67
    assert len(reported["freeboard_mass_gap_kg_s"]) == 2
    assert reported["freeboard_closed_form"] is None
    assert "不新增计算" in reported["freeboard_closed_form_note"]
    ratio = bubble_balance_ratio(0.01, 1000.0, 101300.0, reported["visible_bubble_m3_s"][0])
    assert ratio is not None and ratio > 0.0


def test_campaign_rows_are_wired_and_overlay_changes_the_copy(tmp_path):
    assert [item.name for item in run_specs() if not item.wired] == []
    specs = {item.name: item for item in run_specs()}
    assert config_for("S3").particle_emissivity == 0.79
    assert config_for("S5").interphase_diffusivity == "fuller_o2_n2"
    assert config_for("S13").hover_closure == "B"
    assert config_for("S14").hover_closure == "C"
    assert config_for("S15").serio_frequency_per_s == pytest.approx(9.96e2)
    assert config_for("baseline").hover_closure == "A"
    dest = tmp_path / "s9.json"
    overlay(_CASE, dest, specs["S9"])
    loaded = _load(dest, "1")
    assert loaded.epsilon_mf == pytest.approx(0.40)
    assert loaded.phi_s == pytest.approx(0.6)
    doubled = tmp_path / "s16.json"
    overlay(_CASE, doubled, specs["S16"])
    wide = _load(doubled, "1")
    assert wide.n_bed == 16 and wide.n_free == 16
    assert wide.bed_cell_m == pytest.approx(0.075)
    assert wide.free_cell_m == pytest.approx(0.14375)
    assert interface_balance_index(wide.secondary_height_m, wide.n_bed * wide.bed_cell_m, wide.n_free, wide.free_cell_m) == 13

