"""自由空间参考轨迹的归宿与输运系数。

格内气速为常数，稳定点只在界面（ADR 0047）；粒级内取样与半格捕获区（ADR 0048）。
例子都用两三个格与给定的常气速，不是 Hamel 的鬼泡剖面。
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.freeboard import (
    FreeboardError,
    account_freeboard,
    cell_mean_velocity,
    ejection_velocity_samples,
    ghost_bubble_velocity,
    stable_equilibria,
    terminal_settling_velocity,
    trace_until_fate,
)
from hamel.freeboard.accounting import sphericity_points, within_class_samples

_GAS = 1.2
_MU = 2.0e-5
_DENSITY = 1200.0
_PHI = 0.8
_DIAMETER = 2.0e-4
# 界面上下气速：下格 1,0、上格 0,2；v_t ≈ 0,849 落在两者之间。
_LOWER = 1.0
_UPPER = 0.2


def _config(**changes) -> SolverConfig:
    return replace(SolverConfig(), **changes)


def _relative(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def _column(count: int, length: float, rho_g: float, mu: float):
    edges = np.linspace(0.0, length, count + 1)
    return edges, np.full(count, rho_g), np.full(count, mu)


def _cellwise(*speeds: float):
    def gas(cell, _height):
        return float(speeds[cell])

    return gas


def _account(edges, gas, rho, mu, diameter, rho_p, phi, bubble_velocity, samples, scale=1.0, **config):
    speeds, shares, direct = ejection_velocity_samples(bubble_velocity, samples)
    return account_freeboard(
        edges, gas, rho, mu, diameter, rho_p, phi, speeds, shares, direct, _config(**config), scale
    )


def test_within_class_samples_are_log_uniform_and_mean_to_eight_tenths():
    """质量在 ln d 上均匀的等分中点；球形度截断正态的等概率分位中点。ADR 0048 第 1 项。"""
    low, high = 0.25e-3, 0.5e-3
    points = within_class_samples(low, high)
    assert len(points) == 64
    assert sum(weight for _d, _phi, weight in points) == pytest.approx(1.0)
    diameters = sorted({round(diameter, 12) for diameter, _phi, _w in points})
    assert len(diameters) == 8
    step = (np.log(high) - np.log(low)) / 8.0
    expected = [float(np.exp(np.log(low) + (index + 0.5) * step)) for index in range(8)]
    assert diameters == pytest.approx(expected, rel=1.0e-12)
    phi = sphericity_points()
    assert len(phi) == 8
    assert all(0.3 <= value <= 1.0 for value in phi)
    assert float(np.mean(phi)) == pytest.approx(0.8, abs=2.0e-3)
    assert list(phi) == sorted(phi)


def test_interface_capture_is_recorded_at_the_cell_midpoint():
    """R* 为界面两侧各半格：进入捕获区的时刻是最后一次越过格中点。ADR 0048 第 2 项。"""
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    equilibria = stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI)
    assert len(equilibria) == 1
    assert equilibria[0].height_m == pytest.approx(1.0)
    assert equilibria[0].cells == (1, 2)
    assert equilibria[0].half_width_m == pytest.approx(0.25)
    traced = trace_until_fate(0.2, 2.0, edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, _config(), equilibria)
    assert traced.fate == "capture"
    assert traced.entered_cell in (1, 2)
    assert traced.capture_time_s is not None and traced.capture_time_s > 0.0
    assert traced.velocity_reversals > 2
    # 停留时间只算到进入捕获区为止。
    assert sum(traced.residence_s) == pytest.approx(traced.capture_time_s)
    # 起点在格 0 中点以下，第一次向上越过格中点记入自由状态。
    assert traced.upward_free[0] == 1.0
    assert all(value == 0.0 for value in traced.residence_s[2:]) or traced.entered_cell == 2


def test_downward_entry_goes_to_the_cell_below():
    """向下越过格中点、此后不离开 R*_(i-1) 的，悬停状态记在格 i−1。ADR 0048 第 2 项。"""
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_UPPER, _LOWER, _UPPER)
    equilibria = stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI)
    assert equilibria and equilibria[0].cells == (1, 2)
    # 格 2 的中点起步、初速向下：起点即 R*_1 的分界，按初速方向归下方界面。
    traced = trace_until_fate(1.25, -0.6, edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, _config(), equilibria)
    assert traced.fate == "capture"
    assert traced.entry_kind == "down"
    assert traced.entered_cell == 1
    assert traced.capture_time_s == pytest.approx(0.0)


def test_local_generation_splits_into_up_and_down_shares():
    """当地生成按初速方向归入上方或下方的 R*；两侧都没有时进自由状态。ADR 0048 第 2 项。"""
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_UPPER, _LOWER, _UPPER)
    result = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8)
    assert result.capture_cells == (1, 2)
    for cell in result.cells:
        assert cell.theta_up + cell.theta_down <= 1.0 + 1.0e-12
    # 格 1 的中点以上属于 R*_1（向上）；格 2 的中点以下属于 R*_1（向下）。
    assert result.cells[0].theta_up == 0.0 and result.cells[0].theta_down == 0.0
    assert result.cells[1].theta_up > 0.0 and result.cells[2].theta_down > 0.0


def test_c_terms_and_theta_conserve_mass_per_class():
    """c↑、c↓、θ↑、θ↓ 逐格装配后，参考存量的进出与抛出份额一致。ADR 0048 第 2 项。"""
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    result = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 16)
    direct = ejection_velocity_samples(1.0, 16)[2]
    assert result.ejection_hover == 0.0
    assert result.ejection_free == pytest.approx(1.0 - direct)
    assert result.capture_share + result.return_share + result.elutriation_share == pytest.approx(1.0)
    for cell in result.cells:
        if cell.m_hat is None:
            continue
        assert cell.c_up >= 0.0 and cell.c_down >= 0.0
    assert result.cells[0].c_down == 0.0
    # 自由状态每单位参考存量的出流 = 上行 + 下行 + 格内转入；转入的份额等于捕获份额。
    taken = sum((cell.c_up + cell.c_down) * cell.m_hat for cell in result.cells if cell.m_hat)
    assert taken == pytest.approx(result.capture_share, rel=1.0e-9)


def test_halving_delta_leaves_the_books_unchanged():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    full = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8)
    half = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8, scale=0.5)
    assert half.cells == full.cells
    assert half.return_share == full.return_share
    assert half.elutriation_share == full.elutriation_share
    assert half.capture_share == full.capture_share


def test_doubling_samples_keeps_the_split_coefficients_within_1e_4():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    gas = _cellwise(_LOWER, _LOWER, _UPPER)
    base = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 48)
    fine = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 96)
    assert fine.capture_share == pytest.approx(base.capture_share, rel=1e-4, abs=0.0)
    for left, right in zip(base.cells, fine.cells):
        assert left.free_state == right.free_state
        for name in ("a_f", "b_f", "c_up", "c_down"):
            first, second = getattr(left, name), getattr(right, name)
            if first is None:
                assert second is None
                continue
            assert _relative(first, second) < 1.0e-4


def test_monotonic_profiles_are_not_captured():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    for speed, expect in ((0.2, "return"), (1.0 + 4.0, "elutriate")):
        gas = _cellwise(speed, speed, speed)
        assert stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI) == ()
        result = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, speed, 8)
        direct = ejection_velocity_samples(speed, 8)[2]
        assert result.capture_cells == ()
        assert result.capture_share == 0.0
        assert result.ejection_hover == 0.0
        assert all(cell.c_up == 0.0 and cell.c_down == 0.0 for cell in result.cells if cell.c_up is not None)
        if expect == "return":
            assert result.return_share == pytest.approx(1.0)
            assert result.elutriation_share == pytest.approx(0.0)
        else:
            assert result.elutriation_share == pytest.approx(1.0 - direct)


def test_neutral_cell_reports_and_stops():
    edges, rho, mu = _column(2, 1.0, _GAS, _MU)
    settling = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI)
    gas = _cellwise(settling, settling)
    with pytest.raises(FreeboardError, match="中性"):
        stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI)


def test_reaction_free_solution_matches_the_reference_books():
    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    speed = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI) + 4.0
    gas = _cellwise(speed, speed, speed)
    result = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 48)
    assert result.capture_share == 0.0
    count = len(result.cells)
    matrix = np.zeros((count, count))
    right = np.zeros(count)
    right[0] = result.ejection_free
    for index, cell in enumerate(result.cells):
        assert cell.free_state == "存在"
        matrix[index, index] = cell.a_f + cell.b_f
        if index:
            matrix[index, index - 1] = -result.cells[index - 1].a_f
        if index + 1 < count:
            matrix[index, index + 1] = -result.cells[index + 1].b_f
    inventory = np.linalg.solve(matrix, right)
    for solved, cell in zip(inventory, result.cells):
        assert solved == pytest.approx(cell.m_hat)
    elutriated = result.cells[-1].a_f * inventory[-1]
    assert elutriated == pytest.approx(result.elutriation_share)


def test_bed_samples_reach_hover_only_through_the_c_terms():
    """床面抛出不在任何 R* 内，F^ej,t = 0；捕获经 c↑、c↓ 记账。ADR 0048 第 2 项。"""
    edges, rho, mu = _column(2, 1.0, _GAS, _MU)
    gas = _cellwise(_LOWER, _UPPER)
    result = _account(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI, 1.0, 8)
    assert result.capture_cells == (0, 1)
    assert result.ejection_hover == 0.0
    assert result.capture_share > 0.0
    assert result.cells[0].c_up > 0.0


def test_oversized_gas_region_has_no_interface_even_if_speeds_differ():
    """上下两格都高于 v_t（或都低于）时没有稳定界面。"""
    edges, rho, mu = _column(2, 1.0, _GAS, _MU)
    slow = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI)
    for lower, upper in ((slow + 0.4, slow + 0.2), (0.3 * slow, 0.1 * slow)):
        gas = _cellwise(lower, upper)
        assert stable_equilibria(edges, gas, rho, mu, _DIAMETER, _DENSITY, _PHI) == ()


def test_ghost_bubble_velocity_feeds_the_cell_mean():
    edges = np.array([0.0, 0.2875, 2.3])
    diameters = np.array([0.15, 0.25, 0.25])
    flows = np.array([0.04, 0.07])

    def gas(cell, height):
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, 1.6, 0.05)

    mean = cell_mean_velocity(gas, 0, 0.0, 0.2875)
    samples = [gas(0, 0.2875 * (index + 0.5) / 200) for index in range(200)]
    assert mean == pytest.approx(float(np.mean(samples)), rel=1e-5)


_TRACK = Path(__file__).resolve().parents[1] / "reports" / "t26h-trajectories" / "precalc_old.json"

# 回归快照，非正确性证据：T26h 代码（整格 R*、单一代表粒径、φ_S = 0,8）在工况 1
# 第 1 轮预计算的自由空间气速剖面上给出的系数。按衡算下标（0 为顶格）。
_OLD_BOOKS = {
    0: {"a_f": [2.638029649799, 2.744255674413, 2.874699467031, 3.034890480081, 3.231624675540, 3.473256793634, 3.769916846954, 4.163938172065],
        "b_f": [0.0] * 8, "theta": [0.0] * 8},
    2: {"a_f": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, None],
        "b_f": [2.327001483620, 1.418286528098, 1.122601192794, 0.887542493932, 0.644668108542, 0.0, None, None],
        "theta": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]},
    6: {"a_f": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.254639862742, 0.862158852254],
        "b_f": [26.465546002428, 17.583038328474, 15.734527804972, 14.872782465341, 14.320886542151, 4.023643376474, 3.564958078382, 2.955973207729],
        "theta": [0.0] * 8},
}


@pytest.mark.skipif(not _TRACK.exists(), reason="T26h 的旧预计算快照不在本地")
def test_old_representative_size_and_whole_cell_rule_reproduce_the_old_books():
    """兼容测试：注入旧代表粒径、单一 φ_S = 0,8 与整格 R*，逐位复现 T26h 的系数。

    生产取样（N_d = N_φ = 8、半格 R*）不复现这些数；这里只核对重构没有改变旧规则的读数。
    """
    old = json.loads(_TRACK.read_text(encoding="utf-8"))
    column = old["trajectory"]
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float)
    rho = np.asarray(column["rho_g"], dtype=float)
    mu = np.asarray(column["mu_g"], dtype=float)

    def gas(cell, height):
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, column["bubble_m_s"], column["bubble_diameter_m"])

    speeds, shares, direct = ejection_velocity_samples(column["bubble_m_s"], 48)
    for klass, expected in _OLD_BOOKS.items():
        book = account_freeboard(
            edges, gas, rho, mu, float(column["diameter_m"][klass]), column["rho_p"], 0.8,
            speeds, shares, direct, _config(), region="cell",
        )
        for geometric, cell in enumerate(book.cells):
            balance = len(book.cells) - 1 - geometric
            first, second = expected["a_f"][balance], expected["b_f"][balance]
            if cell.m_hat is None:
                assert first is None and second is None
                continue
            assert cell.a_f == pytest.approx(first, rel=1.0e-10, abs=1.0e-14)
            assert cell.b_f == pytest.approx(second, rel=1.0e-10, abs=1.0e-14)
            assert cell.theta_up == pytest.approx(expected["theta"][balance], abs=1.0e-12)
