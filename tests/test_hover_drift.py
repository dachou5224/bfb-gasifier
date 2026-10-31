"""悬停状态的粒级内漂移（ADR 0051）：窗口、去向、δ^h、δ^f、μ 与衡算。

例子都用两到四个格与给定的常气速、逐格物性，不是 Hamel 的鬼泡剖面。
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.freeboard import ejection_velocity_samples, terminal_settling_velocity
from hamel.freeboard.accounting import (
    DRIFT_BURN,
    DRIFT_EMPTY,
    DRIFT_FREE,
    DRIFT_HOVER,
    DRIFT_MIGRATION,
    ClassWindow,
    account_freeboard,
    drift_windows,
)
from hamel.solver import evaluate_residual, pack_state, unknown_slots
from hamel.solver import newton as newton_module
from hamel.solver.frozen import AssemblyError, FreeboardDomainError
from hamel.solver.vector import jacobian_layout
from hamel.freeboard import stable_equilibria, trace_until_fate
from hamel.solver.sensitivity import solid_enthalpy_w
from tests.test_assembly import make_frozen, make_state

_RHO_P = 1000.0
_PHI = 0.8
# ADR 0051 审查给出的变物性反例：格 j+1、j+2 的密度与粘度不同。
_PROPS = ((5.0, 4.0e-5), (5.0, 4.0e-5), (3.0, 2.0e-5), (5.0, 4.0e-5))
_DIAMETERS_M = (0.45e-3, 0.40e-3, 0.35e-3, 0.30e-3)
_LOWER_M, _UPPER_M = 0.25e-3, 0.5e-3


def _config(**changes) -> SolverConfig:
    return replace(SolverConfig(), **changes)


def _settling(rho_g: float, mu_g: float, diameter_m: float, phi_s: float = _PHI) -> float:
    return terminal_settling_velocity(rho_g, mu_g, diameter_m, _RHO_P, phi_s)


def _cellwise(*speeds: float):
    def gas(cell, _height):
        return float(speeds[cell])

    return gas


def _column(count: int, length: float = 0.5):
    return np.linspace(0.0, length * count, count + 1)


def _interface_reference():
    """反例：ū_{j+1} = v_t(0,4 mm)、ū_{j+2} = v_t(0,3 mm)，ū_{j+1} < ū_{j+2}，界面 j+1 仍稳定。"""
    props = _PROPS[:3]
    rho = np.array([item[0] for item in props])
    mu = np.array([item[1] for item in props])
    speeds = [_settling(rho[0], mu[0], 0.6e-3), _settling(rho[1], mu[1], 0.4e-3), _settling(rho[2], mu[2], 0.3e-3)]
    return _column(3), _cellwise(*speeds), rho, mu, speeds


def _account(edges, gas, rho, mu, samples, window, first_class=False, bubble_velocity=1.0, sample_count=48):
    speeds, shares, direct = ejection_velocity_samples(bubble_velocity, sample_count)
    return account_freeboard(
        edges, gas, rho, mu, float(_DIAMETERS_M[2]), _RHO_P, _PHI, speeds, shares, direct,
        _config(), samples=samples, class_window=window, first_class=first_class,
    )


def test_adjacent_windows_join_and_empty_windows_do_not_contribute():
    """相邻界面的窗口首尾相接；气速单调时没有窗口，δ 与 μ 都为 0。ADR 0051 第 2 项。"""
    edges, gas, rho, mu, _speeds = _interface_reference()
    window = ClassWindow(0.10e-3, 1.0e-3, 1.0)
    table = drift_windows(edges, gas, rho, mu, _RHO_P, [_PHI], window)
    assert table[0][0].low_m == pytest.approx(table[1][0].high_m, rel=0.0, abs=0.0)
    assert table[0][0].high_m > table[0][0].low_m
    assert table[0][0].fate == DRIFT_HOVER
    result = _account(edges, gas, rho, mu, ((_DIAMETERS_M[2], _PHI, 1.0),), ClassWindow(_LOWER_M, _UPPER_M, 1.0))
    assert result.cells[0].drift_hover > 0.0 and result.cells[0].drift_free == 0.0
    flat = drift_windows(edges, _cellwise(2.0, 2.0, 2.0), rho, mu, _RHO_P, [_PHI], window)
    assert all(item.fate == DRIFT_EMPTY and item.exit_share == 0.0 for row in flat for item in row)
    flat_book = _account(edges, _cellwise(2.0, 2.0, 2.0), rho, mu, ((_DIAMETERS_M[2], _PHI, 1.0),), ClassWindow(_LOWER_M, _UPPER_M, 1.0))
    assert all(cell.drift_hover == 0.0 and cell.drift_free == 0.0 and cell.migration_mu == 0.0 for cell in flat_book.cells)


def test_variable_properties_keep_the_next_interface_stable():
    """只比较气速会判错：ū_{j+1} < ū_{j+2} 时界面 j+1 仍稳定，去向是悬停。ADR 0051 审查反例。"""
    edges, gas, rho, mu, speeds = _interface_reference()
    assert speeds[1] < speeds[2]
    windows = drift_windows(edges, gas, rho, mu, _RHO_P, [_PHI], ClassWindow(_LOWER_M, _UPPER_M, 1.0))
    assert windows[0][0].fate == DRIFT_HOVER
    # 离开点的粒径等于界面 j+1 窗口的上沿；下沿由格 j+2 自己的物性给出，低于离开点。
    assert windows[0][0].low_m == pytest.approx(0.4e-3, rel=1e-12)
    assert _settling(rho[1], mu[1], windows[0][0].low_m) == pytest.approx(speeds[1], rel=1e-12)
    assert _settling(rho[2], mu[2], windows[1][0].low_m) == pytest.approx(speeds[2], rel=1e-9)


def test_topmost_target_cell_sends_the_drift_to_the_free_state():
    """格 j+1 是最上一格时没有界面 j+1，漂移进入自由状态。ADR 0051 第 2 项。"""
    edges, gas, rho, mu, _speeds = _interface_reference()
    windows = drift_windows(edges[:3], gas, rho[:2], mu[:2], _RHO_P, [_PHI], ClassWindow(_LOWER_M, _UPPER_M, 1.0))
    assert windows[0][0].fate == DRIFT_FREE


def test_facets_with_different_fates_accumulate_separately():
    """同一（格，粒级）不同 φ_p 的去向不同：δ^h 与 δ^f 分别累计，不合并成一个去向。"""
    rho = np.array([item[0] for item in _PROPS])
    mu = np.array([item[1] for item in _PROPS])
    speeds = [_settling(rho[index], mu[index], diameter) for index, diameter in enumerate(_DIAMETERS_M)]
    edges = _column(4)
    faceted = [_PHI, 1.0]
    windows = drift_windows(edges, _cellwise(*speeds), rho, mu, _RHO_P, faceted, ClassWindow(0.2e-3, 0.6e-3, 1.0))
    # 格 1：φ = 0,8 漂移进悬停，φ = 1,0 漂移进自由。
    assert [item.fate for item in windows[1]] == [DRIFT_HOVER, DRIFT_FREE]
    samples = ((0.37e-3, _PHI, 0.5), (0.25e-3, 1.0, 0.5))
    book = _account(edges, _cellwise(*speeds), rho, mu, samples, ClassWindow(0.2e-3, 0.6e-3, 1.0))
    cell = book.cells[1]
    assert cell.drift_hover > 0.0 and cell.drift_free > 0.0 and cell.migration_mu == 0.0
    assert [window.weight for window in cell.windows] == pytest.approx([0.5, 0.5])
    expected = {fate: sum(part.exit_share * part.weight for part in cell.windows if part.fate == fate) for fate in (DRIFT_HOVER, DRIFT_FREE)}
    assert cell.drift_hover == pytest.approx(expected[DRIFT_HOVER], rel=1e-12)
    assert cell.drift_free == pytest.approx(expected[DRIFT_FREE], rel=1e-12)
    assert sum(part.weight for part in cell.windows) == pytest.approx(1.0)


def test_full_window_degrades_to_the_adr_0044_migration_factor():
    """窗口盖满粒级、下沿是 a_k 时 δ = 0、μ = q_k。ADR 0051 第 2 项。"""
    edges = _column(2)
    rho, mu = np.full(2, 1.2), np.full(2, 2.0e-5)
    previous_m, diameter_m = 0.2e-3, 0.35e-3
    window = ClassWindow(_LOWER_M, _UPPER_M, math.log(diameter_m / previous_m) / math.log(_UPPER_M / _LOWER_M))
    speeds = [2.0 * _settling(1.2, 2.0e-5, _UPPER_M), 0.5 * _settling(1.2, 2.0e-5, _LOWER_M)]
    book = _account(edges, _cellwise(*speeds), rho, mu, ((diameter_m, _PHI, 1.0),), window)
    cell = book.cells[0]
    assert cell.windows[0].fate == DRIFT_MIGRATION
    assert cell.drift_hover == 0.0 and cell.drift_free == 0.0
    cube = (previous_m / diameter_m) ** 3
    assert cell.migration_mu == pytest.approx(cube / (1.0 - cube), rel=1e-13)


def test_first_class_stays_and_burns_without_migration():
    """粒级 1 的窗口下沿为 a_1 的部分留在本状态燃尽，μ_j,1 = 0；下沿在 a_1 以上照常漂移。"""
    edges = _column(2)
    rho, mu = np.full(2, 1.2), np.full(2, 2.0e-5)
    window = ClassWindow(_LOWER_M, _UPPER_M, 1.0)
    full = [2.0 * _settling(1.2, 2.0e-5, _UPPER_M), 0.5 * _settling(1.2, 2.0e-5, _LOWER_M)]
    book = _account(edges, _cellwise(*full), rho, mu, ((0.35e-3, _PHI, 1.0),), window, first_class=True)
    assert book.cells[0].windows[0].fate == DRIFT_BURN
    assert book.cells[0].migration_mu == 0.0 and book.cells[0].drift_free == 0.0 and book.cells[0].drift_hover == 0.0
    partial = [2.0 * _settling(1.2, 2.0e-5, _UPPER_M), _settling(1.2, 2.0e-5, 0.35e-3)]
    drift = _account(edges, _cellwise(*partial), rho, mu, ((0.45e-3, _PHI, 1.0),), window, first_class=True)
    assert drift.cells[0].windows[0].fate == DRIFT_FREE
    assert drift.cells[0].migration_mu == 0.0 and drift.cells[0].drift_free > 0.0


def test_window_chain_path_length_stays_inside_the_class():
    """同一 φ_p 沿窗口链的路径长度之和不超过 ρ_k·ln(b_k/a_k)。ADR 0051 第 2 项。"""
    edges, gas, rho, mu, _speeds = _interface_reference()
    window = ClassWindow(0.28e-3, 0.48e-3, 1.0)
    table = drift_windows(edges, gas, rho, mu, _RHO_P, [_PHI], window)
    total = sum(row[0].length for row in table)
    assert 0.0 < total <= math.log(window.upper_m / window.lower_m) + 1.0e-12


def test_hand_computed_two_window_example():
    """两窗口：上格进入自由状态、下格进入悬停状态，δ 由 1/(e^{3L}−1) 手算。"""
    edges, gas, rho, mu, _speeds = _interface_reference()
    window = ClassWindow(_LOWER_M, _UPPER_M, 1.0)
    table = drift_windows(edges, gas, rho, mu, _RHO_P, [_PHI], window)
    assert [row[0].fate for row in table] == [DRIFT_HOVER, DRIFT_FREE]
    book = _account(edges, gas, rho, mu, ((_DIAMETERS_M[2], _PHI, 1.0),), window)
    for cell, row in zip(book.cells, table):
        item = row[0]
        length = math.log(item.high_m / item.low_m)
        assert item.length == pytest.approx(length, rel=1e-12)
        assert item.exit_share == pytest.approx(1.0 / math.expm1(3.0 * length), rel=1e-12)
    assert book.cells[0].drift_hover == pytest.approx(table[0][0].exit_share, rel=1e-12)
    assert book.cells[1].drift_free == pytest.approx(table[1][0].exit_share, rel=1e-12)
    assert book.cells[0].drift_free == 0.0 and book.cells[1].drift_hover == 0.0


# ω 归集的手算例：同一气速剖面上两个 φ 分别向上、向下被捕获，落在不同的界面。
_OMEGA_SPEEDS = (1.2, 1.1, 0.5)
_OMEGA_CLASS = ClassWindow(0.15e-3, 0.34e-3, 1.0)  # ρ_k = 1，故 L = ln(d′_hi/d′_lo)
_OMEGA_SAMPLES = ((0.23e-3, 1.0, 0.5), (0.25e-3, 0.7, 0.5))


def _diameter_at(speed: float, phi: float) -> float:
    """v_t(d, φ) = speed 的根，测试内独立二分解，不复用模块的求根。"""
    low, high = 5.0e-5, 2.0e-3
    for _ in range(200):
        middle = 0.5 * (low + high)
        if terminal_settling_velocity(1.2, 2.0e-5, middle, _RHO_P, phi) < speed:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _hand_window(cell: int, phi: float, lower: float, upper: float):
    """该格该 φ 的窗口 [d′_lo, d′_hi] 与出口因子，按 ADR 0051 第 2 项手算。"""
    high = min(_diameter_at(_OMEGA_SPEEDS[cell], phi), upper)
    low = max(_diameter_at(_OMEGA_SPEEDS[cell + 1], phi), lower)
    length = math.log(high / low)
    return low, high, length, 1.0 / math.expm1(3.0 * length)


def _omega_book():
    edges = _column(3, 0.5)
    return account_freeboard(
        edges, _cellwise(*_OMEGA_SPEEDS), np.full(3, 1.2), np.full(3, 2.0e-5), 0.245e-3, _RHO_P, 0.8,
        np.array([]), np.array([]), 1.0, _config(), samples=_OMEGA_SAMPLES, class_window=_OMEGA_CLASS,
    )


def test_omega_is_collected_at_the_capture_target_cell():
    """ω_p 归集到实际悬停目标格：向下进入时是源格下方相邻格，不按源格计。审查意见第 2 条。"""
    edges = _column(3, 0.5)
    rho, mu = np.full(3, 1.2), np.full(3, 2.0e-5)
    gas = _cellwise(*_OMEGA_SPEEDS)
    # 同一格中点起步：φ=1,0 的 0,21 mm 向下进入 R*_0（目标 0），φ=0,7 的 0,28 mm 向上（目标 1）。
    for diameter, phi, kind, target in ((0.23e-3, 1.0, "down", 0), (0.25e-3, 0.7, "up", 1)):
        terminal = terminal_settling_velocity(1.2, 2.0e-5, diameter, _RHO_P, phi)
        equilibria = stable_equilibria(edges, gas, rho, mu, diameter, _RHO_P, phi)
        middle = 0.5 * (edges[1] + edges[2])
        traced = trace_until_fate(
            middle, _OMEGA_SPEEDS[1] - terminal, edges, gas, rho, mu, diameter, _RHO_P, phi, _config(), equilibria,
        )
        assert traced.fate == "capture"
        assert traced.entry_kind == kind and traced.entered_cell == target
        assert traced.entry_cell == 1
    book = _omega_book()
    upper = {window.sphericity: window for window in book.cells[0].windows}
    lower = {window.sphericity: window for window in book.cells[1].windows}
    # 目标 0 只收到 φ=1,0，目标 1 只收到 φ=0,7；源格计数不改变这条归集。
    assert upper[1.0].weight == pytest.approx(1.0) and upper[0.7].weight == 0.0
    assert lower[0.7].weight == pytest.approx(1.0) and lower[1.0].weight == 0.0
    # 权重不得落到空窗口上（旧写法把捕获写回源格时会落到这里）。
    for cell in book.cells:
        for window in cell.windows:
            if window.weight > 0.0:
                assert window.fate != DRIFT_EMPTY


def test_omega_hand_computes_drift_factors_per_target_cell():
    """两个 φ 各自的目标窗口不同：δ^h、δ^f、μ 由 ω 乘各自的出口因子手算。"""
    book = _omega_book()
    upper = {window.sphericity: window for window in book.cells[0].windows}
    lower = {window.sphericity: window for window in book.cells[1].windows}
    low, high, length, share = _hand_window(0, 1.0, _OMEGA_CLASS.lower_m, _OMEGA_CLASS.upper_m)
    assert upper[1.0].low_m == pytest.approx(low, rel=1e-9)
    assert upper[1.0].high_m == pytest.approx(high, rel=1e-9)
    assert upper[1.0].length == pytest.approx(length, rel=1e-9)
    assert upper[1.0].exit_share == pytest.approx(share, rel=1e-9)
    assert book.cells[0].drift_hover == pytest.approx(share, rel=1e-9)
    assert book.cells[0].drift_free == 0.0 and book.cells[0].migration_mu == 0.0
    low, high, length, share = _hand_window(1, 0.7, _OMEGA_CLASS.lower_m, _OMEGA_CLASS.upper_m)
    assert lower[0.7].low_m == pytest.approx(low, rel=1e-9)
    assert lower[0.7].high_m == pytest.approx(high, rel=1e-9)
    assert lower[0.7].length == pytest.approx(length, rel=1e-9)
    assert lower[0.7].exit_share == pytest.approx(share, rel=1e-9)
    assert book.cells[1].drift_free == pytest.approx(share, rel=1e-9)
    assert book.cells[1].drift_hover == 0.0 and book.cells[1].migration_mu == 0.0


def test_drift_target_check_uses_the_upward_neighbour():
    """衡算下标 0 为顶，漂移从 i 到 i−1：源 [1:] 对目标 [:-1]。审查意见第 1 条的探针。"""
    # 审查探针：3 格 2 粒级，只有 hover[2, 1]，drift_hover_per_s[2, 1] = 1，目标 hover[1, 1] 为假。
    probe = np.zeros((3, 2), dtype=bool)
    probe[2, 1] = True
    drift = np.zeros((3, 2))
    drift[2, 1] = 1.0
    with pytest.raises(AssemblyError, match="进入悬停状态的漂移目标没有悬停状态"):
        make_frozen(n_free=3, n_class=2, hover=probe, drift_hover_per_s=drift)
    # 最下源格合法、目标存在时通过；目标上方那格没有悬停状态不影响。
    valid = probe.copy()
    valid[1, 1] = True
    frozen = make_frozen(n_free=3, n_class=2, hover=valid, drift_hover_per_s=drift)
    assert frozen.drift_hover_per_s[2, 1] == 1.0
    assert not frozen.hover[0, 1]


def test_drift_target_check_does_not_reject_a_legal_upward_step():
    """源在衡算下标 1、目标 0：目标有悬停状态即为合法，不因下格没有悬停状态被拒。"""
    hover = np.zeros((3, 2), dtype=bool)
    hover[0, 0] = True
    hover[1, 0] = True
    drift = np.zeros((3, 2))
    drift[1, 0] = 1.0
    frozen = make_frozen(n_free=3, n_class=2, hover=hover, drift_hover_per_s=drift)
    assert frozen.drift_hover_per_s[1, 0] == 1.0
    assert not frozen.hover[2, 0]


def test_drift_target_check_rejects_a_missing_free_state_and_top_source():
    """漂移要进入的自由状态不存在时拒绝；最上一格不能作为源。"""
    hover = np.ones((3, 2), dtype=bool)
    absent = np.zeros((3, 2), dtype=bool)
    absent[0, 0] = True
    theta = np.zeros((3, 2))
    theta[0, 0] = 1.0
    drift_free = np.zeros((3, 2))
    drift_free[1, 0] = 0.5
    with pytest.raises(AssemblyError, match="进入自由状态的漂移目标没有自由状态"):
        make_frozen(n_free=3, n_class=2, hover=hover, free_absent=absent, theta=theta, drift_free_per_s=drift_free)
    # 最上一格（衡算下标 0）没有上方相邻格，那里不能有漂移。
    top = np.zeros((3, 2))
    top[0, 1] = 1.0
    with pytest.raises(AssemblyError, match="最上一格没有上方相邻格"):
        make_frozen(n_free=3, n_class=2, hover=hover, drift_hover_per_s=top)


_DRIFT_HOVER = 0.3
_DRIFT_FREE = 0.2
# 冒烟格：粒级 2（k = 2）。粒级 1 的悬停状态收到灰流时按 ADR 0026 第 5 项报告并停止，另有测试。
_KLASS = 1


def _drift_frozen(drift=True, klass=_KLASS, hover=None, migration_mu=None):
    mask = np.ones((3, 2), dtype=bool) if hover is None else hover
    drift_hover = np.zeros((3, 2))
    drift_free = np.zeros((3, 2))
    if drift:
        drift_hover[2, klass], drift_free[2, klass] = _DRIFT_HOVER, _DRIFT_FREE
    return make_frozen(
        n_free=3, n_class=2, hover=mask,
        migration_mu_per_s=np.zeros((3, 2)) if migration_mu is None else migration_mu,
        drift_hover_per_s=drift_hover, drift_free_per_s=drift_free,
    )


def _solid_rows(frozen, state, kind: str, phase: str, klass: int) -> dict:
    assembly = evaluate_residual(pack_state(state, frozen), frozen)
    rows = {}
    for slot in unknown_slots(frozen):
        if slot.kind == kind and slot.zone == "freeboard" and slot.size_class == klass and slot.phase == phase:
            rows[slot.balance_index] = float(assembly.values[slot.index])
    return rows


def test_drift_moves_char_and_ash_between_the_two_cells_only():
    """漂移的焦炭与灰只在源格与上方相邻格之间转移，别的行的残差不变。ADR 0051 第 3 项。"""
    still, moving = _drift_frozen(drift=False), _drift_frozen()
    state = make_state(still)
    for kind in ("char", "ash"):
        hover_before = _solid_rows(still, state, kind, "hover", _KLASS)
        hover_after = _solid_rows(moving, state, kind, "hover", _KLASS)
        free_before = _solid_rows(still, state, kind, "free", _KLASS)
        free_after = _solid_rows(moving, state, kind, "free", _KLASS)
        moved = hover_before[2] - hover_after[2]
        assert moved > 0.0
        assert hover_after[0] == hover_before[0] and free_after[2] == free_before[2]
        share = _DRIFT_HOVER / (_DRIFT_HOVER + _DRIFT_FREE)
        assert hover_after[1] - hover_before[1] == pytest.approx(moved * share, rel=1e-12)
        assert free_after[1] - free_before[1] == pytest.approx(moved * (1.0 - share), rel=1e-12)
        total = sum(hover_after[cell] - hover_before[cell] + free_after[cell] - free_before[cell] for cell in hover_after)
        assert total == pytest.approx(0.0, abs=1e-12)


def test_first_class_drift_into_the_free_state_conserves_char_and_ash():
    """粒级 1 的漂移进入上方格的自由状态：焦炭与灰逐格扣除与接收，μ_j,1 = 0。"""
    hover = np.zeros((3, 2), dtype=bool)
    hover[2, 0] = True  # 目标格只有自由状态：k = 1 的悬停状态收灰时按 ADR 0026 第 5 项停止
    still = make_frozen(n_free=3, n_class=2, hover=hover)
    drift_free = np.zeros((3, 2))
    drift_free[2, 0] = 0.5
    moving = make_frozen(n_free=3, n_class=2, hover=hover, drift_free_per_s=drift_free)
    state = make_state(still)
    for kind in ("char", "ash"):
        hover_before = _solid_rows(still, state, kind, "hover", 0)
        hover_after = _solid_rows(moving, state, kind, "hover", 0)
        free_before = _solid_rows(still, state, kind, "free", 0)
        free_after = _solid_rows(moving, state, kind, "free", 0)
        moved = hover_before[2] - hover_after[2]
        assert moved > 0.0
        assert free_after[1] - free_before[1] == pytest.approx(moved, rel=1e-12)
        assert free_after[2] == free_before[2] and hover_after.get(1, 0.0) == 0.0
        assert moving.migration_mu_per_s is None
        assert moving.hover[1, 0] == False


def test_drift_into_the_first_class_hover_state_reports_the_domain_rule():
    """粒级 1 的悬停状态收到灰流时按 ADR 0026 第 5 项报告并停止，本 ADR 不改这条规则。"""
    frozen = _drift_frozen(klass=0)
    with pytest.raises(FreeboardDomainError, match="最小粒级"):
        evaluate_residual(pack_state(make_state(frozen), frozen), frozen)


def _energy_rows(frozen, state) -> dict:
    assembly = evaluate_residual(pack_state(state, frozen), frozen)
    rows = {}
    for slot in unknown_slots(frozen):
        if slot.kind == "temperature" and slot.zone == "freeboard":
            rows[slot.balance_index] = float(assembly.values[slot.index])
    return rows


def _temperatures(*values: float) -> np.ndarray:
    return np.asarray(values, dtype=float)


def test_drift_enthalpy_is_valued_at_the_source_cell_temperature():
    """漂移流按源格温度计焓，在两格的能量方程中逐格扣除与接收。ADR 0051 第 3 项。"""
    still, moving = _drift_frozen(drift=False), _drift_frozen()
    hot_source = make_state(still)
    hot_source.freeboard_temperature_k = _temperatures(1150.0, 1100.0, 1000.0)
    cool_target = make_state(still)
    cool_target.freeboard_temperature_k = _temperatures(1150.0, 900.0, 1000.0)
    for state in (hot_source, cool_target):
        before = _energy_rows(still, state)
        after = _energy_rows(moving, state)
        source_temperature = float(state.freeboard_temperature_k[2])
        char_mass = -(_solid_rows(moving, state, "char", "hover", _KLASS)[2] - _solid_rows(still, state, "char", "hover", _KLASS)[2])
        ash_mass = -(_solid_rows(moving, state, "ash", "hover", _KLASS)[2] - _solid_rows(still, state, "ash", "hover", _KLASS)[2])
        expected = solid_enthalpy_w(char_mass, ash_mass, 0.0, source_temperature, moving)
        assert after[1] - before[1] == pytest.approx(expected, rel=1e-10, abs=1e-9)
        assert after[2] - before[2] == pytest.approx(-expected, rel=1e-10, abs=1e-9)
        assert after[0] - before[0] == pytest.approx(0.0, abs=1e-9)
    hot = _energy_rows(moving, hot_source)
    cool = _energy_rows(moving, cool_target)
    still_hot, still_cool = _energy_rows(still, hot_source), _energy_rows(still, cool_target)
    # 只改接收格（格 1）的温度时，漂移带来的两格能量差不变；改源格温度才会改变。
    assert hot[1] - still_hot[1] == pytest.approx(cool[1] - still_cool[1], rel=1e-12)
    hotter = make_state(still)
    hotter.freeboard_temperature_k = _temperatures(1150.0, 1100.0, 1100.0)
    warm = _energy_rows(moving, hotter)
    assert warm[1] - hot[1] != pytest.approx(0.0, abs=1e-6)


def test_degenerate_mu_matches_the_adr_0044_residual_bitwise():
    """μ 取 ADR 0044 的 q_k、δ = 0 时，残差与 ADR 0044 逐位相同。ADR 0051 实现约束。"""
    bare = make_frozen(n_free=3, n_class=2, hover=np.ones((3, 2), dtype=bool))
    mass = bare.rho_p_kg_m3 * math.pi * bare.diameter_m**3 / 6.0
    step = np.zeros_like(mass)
    step[1:] = mass[:-1] / (mass[1:] - mass[:-1])
    degenerate = make_frozen(
        n_free=3, n_class=2, hover=np.ones((3, 2), dtype=bool),
        migration_mu_per_s=np.tile(step, (3, 1)),
    )
    state = make_state(bare)
    assert np.array_equal(
        evaluate_residual(pack_state(state, bare), bare).values,
        evaluate_residual(pack_state(state, degenerate), degenerate).values,
    )
    assert degenerate.migration_mu_per_s[0, 0] == 0.0


def test_drift_only_couples_adjacent_cells_in_the_grouped_jacobian():
    """漂移只连接相邻格，分组雅可比与逐列差分逐位相同。ADR 0046 第 3、5 项。"""
    frozen = _drift_frozen(migration_mu=np.array([[0.0, 0.02], [0.0, 0.02], [0.0, 0.02]]))
    state = make_state(frozen)
    x = pack_state(state, frozen)
    scale = np.ones(x.size)
    evaluate = lambda trial: evaluate_residual(trial, frozen).values
    columnwise = newton_module._jacobian(x, evaluate(x), evaluate, scale, lambda trial: True)
    grouped, _stats = newton_module._grouped_jacobian(x, evaluate(x), evaluate, scale, lambda trial: True, jacobian_layout(unknown_slots(frozen), frozen))
    assert np.array_equal(grouped, columnwise)
