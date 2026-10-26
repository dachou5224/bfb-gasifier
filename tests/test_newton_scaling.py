"""差分步长与行列缩放。求解手段：ADR 0041。阻尼、可行条件与收敛门的测试仍在 test_driver.py。"""
import math

import numpy as np

from hamel.solver.gates import ResidualBlocks
from hamel.solver.newton import _SQRT_EPS, _difference_basis, _jacobian, _newton_step


def test_relative_step_matches_the_power_derivative_at_a_tiny_value():
    """f(c) = −c − K·c^0.3。c = 1e-20 时新步长对准解析导数，旧步长不准。"""
    concentration = 1.0e-20
    coefficient = 1.0
    scale = np.array([1.0])
    point = np.array([concentration])

    def residual(trial):
        value = float(trial[0])
        return np.array([-value - coefficient * value ** 0.3])

    base = residual(point)
    matrix = _jacobian(point, base, residual, scale, lambda trial: True)
    analytic = -1.0 - coefficient * 0.3 * concentration ** (-0.7)
    assert abs(matrix[0, 0] - analytic) / abs(analytic) < 1.0e-6

    old_step = _SQRT_EPS * max(abs(concentration), float(scale[0]))
    moved = point.copy()
    moved[0] += old_step
    old_derivative = (residual(moved)[0] - base[0]) / old_step
    assert abs(old_derivative - analytic) / abs(analytic) >= 1.0e-6
    assert _difference_basis(concentration, scale[0]) == concentration
    assert _difference_basis(0.0, 3.5) == 3.5
    assert math.isclose(_SQRT_EPS * _difference_basis(0.0, 3.5), _SQRT_EPS * 3.5)


def test_scaled_newton_step_matches_the_unscaled_step_when_well_conditioned():
    """良态线性问题。行列缩放不改变牛顿步。"""
    state = np.array([1.5, 0.5, 2.5])
    target = np.array([1.0, 2.0, 3.0])
    matrix_a = np.array([[2.0, 0.1, 0.0], [0.1, 1.5, 0.2], [0.0, 0.2, 3.0]])
    scales = np.array([1.0, 1.0, 1.0])
    blocks = ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 4.0, 2.0, 8.0)

    def residual(trial):
        return matrix_a @ np.asarray(trial, dtype=float) - target

    base = residual(state)
    matrix = _jacobian(state, base, residual, scales, lambda trial: True)
    plain = np.linalg.solve(matrix, -base)
    scaled = _newton_step(matrix, base, state, scales, blocks)
    assert np.linalg.norm(scaled - plain) / np.linalg.norm(plain) < 1.0e-10
