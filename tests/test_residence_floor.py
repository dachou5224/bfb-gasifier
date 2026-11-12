"""近中性停留时间下限速度。解释性选择：ADR 0057。"""
from __future__ import annotations

import math

import numpy as np
import pytest
from dataclasses import replace

from hamel.cases.config import SolverConfig
from hamel.freeboard import FreeboardError, account_freeboard, terminal_settling_velocity, trace_until_fate
from hamel.freeboard.accounting import _classify

RHO_G, MU, RHO_P, PHI = 1.2, 2.0e-5, 1200.0, 0.8


def _config(**changes):
    return replace(SolverConfig(), **changes)


def _segment(cell, start, end, slip, h0, h1):
    return (cell, start, end, slip, h0, h1)


def _fate(segments, *, u_min, capture_time=None, capture_height=None, slip_piece_check=True):
    midpoints = [] if capture_time is None else [(capture_time, 0, "up")]
    fate = "capture" if capture_time is not None else "elutriate"
    return _classify(
        1, fate, segments[-1][2], segments[-1][5], 0.0, 0, 0, list(segments), [], midpoints,
        frozenset({0}) if capture_time is not None else frozenset(), "half", 1, None, u_min, capture_height,
    )


def test_near_neutral_residence_is_height_over_u_min():
    """|ε| 远小于 u_min 但不中性：计入停留为 Δh/u_min；u_min = 0 恢复解析时间。"""
    terminal = terminal_settling_velocity(RHO_G, MU, 2.0e-4, RHO_P, PHI)
    edges = np.array([0.0, 1.0])
    rho, mu = np.array([RHO_G]), np.array([MU])

    def gas(_cell, _height, terminal=terminal):
        return terminal + 1.0e-6

    plain = trace_until_fate(0.0, 0.0, edges, gas, rho, mu, 2.0e-4, RHO_P, PHI, _config(freeboard_residence_floor_m_s=0.0))
    capped = trace_until_fate(0.0, 0.0, edges, gas, rho, mu, 2.0e-4, RHO_P, PHI, _config(freeboard_residence_floor_m_s=1.0e-2))
    assert plain.fate == capped.fate == "elutriate"
    assert plain.time_s > 80.0 and plain.time_s == pytest.approx(capped.time_s)
    assert sum(plain.residence_s) == pytest.approx(plain.time_s)
    assert sum(capped.residence_s) == pytest.approx(1.0 / 1.0e-2)
    assert plain.entry_count == capped.entry_count


def test_fast_segment_is_bitwise_unchanged():
    """|Δh|/t ≥ u_min 时计入值与 u_min = 0 逐位相同。"""
    segment = [_segment(0, 0.0, 1.0, 0.4, 0.0, 0.5)]
    plain = _fate(segment, u_min=0.0)
    capped = _fate(segment, u_min=1.0e-2)
    assert plain.residence_s == capped.residence_s
    assert plain.slip_m == capped.slip_m
    assert plain.residence_s[0] == 1.0


def test_capture_inside_segment_uses_the_free_portion():
    """整段 20 s、位移 0.30 m；捕获前 10 s、位移 0.05 m；u_min = 0.01 计入 5 s，不是 10 s。"""
    segment = [_segment(0, 0.0, 20.0, 1.0, 0.0, 0.30)]
    counted = _fate(segment, u_min=0.01, capture_time=10.0, capture_height=0.05)
    plain = _fate(segment, u_min=0.0, capture_time=10.0, capture_height=0.05)
    assert counted.residence_s[0] == pytest.approx(5.0)
    assert plain.residence_s[0] == pytest.approx(10.0)
    assert counted.slip_m[0] == pytest.approx(1.0 * 5.0 / 20.0)
    assert plain.slip_m[0] == pytest.approx(1.0 * 10.0 / 20.0)
    assert counted.slip_m[0] / plain.slip_m[0] == pytest.approx(counted.residence_s[0] / plain.residence_s[0])


def test_negative_or_non_finite_floor_is_rejected():
    edges = np.array([0.0, 1.0])
    rho, mu = np.array([RHO_G]), np.array([MU])

    def gas(_cell, _height):
        return 2.0

    for value in (-1.0e-3, math.nan, math.inf):
        with pytest.raises(FreeboardError, match="freeboard_residence_floor_m_s"):
            trace_until_fate(0.0, 1.0, edges, gas, rho, mu, 2.0e-4, RHO_P, PHI, _config(freeboard_residence_floor_m_s=value))


def test_fate_crossings_and_theta_do_not_depend_on_the_floor():
    edges = np.array([0.0, 0.5, 1.0])
    rho, mu = np.full(2, RHO_G), np.full(2, MU)
    diameter = 2.0e-4
    terminal = terminal_settling_velocity(RHO_G, MU, diameter, RHO_P, PHI)

    def gas(cell, _height, terminal=terminal):
        return (terminal + 1.0, terminal - 0.2)[cell]

    speeds = np.array([1.5, 2.0])
    shares = np.array([0.4, 0.6])
    books = []
    for u_min in (0.0, 1.0e-2):
        books.append(account_freeboard(
            edges, gas, rho, mu, diameter, RHO_P, PHI, speeds, shares, 0.0, _config(freeboard_residence_floor_m_s=u_min),
        ))
    plain, capped = books
    assert (plain.return_share, plain.elutriation_share, plain.capture_share) == (capped.return_share, capped.elutriation_share, capped.capture_share)
    for left, right in zip(plain.cells, capped.cells):
        assert (left.theta_up, left.theta_down, left.free_state) == (right.theta_up, right.theta_down, right.free_state)
