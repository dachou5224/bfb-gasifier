"""差分步长与行列缩放。求解手段：ADR 0043。阻尼、可行条件与收敛门的测试仍在 test_driver.py。"""
import numpy as np

from hamel.solver.gates import ResidualBlocks
from hamel.solver.newton import _SQRT_EPS, _column_scale, _difference_basis, _jacobian, _newton_step


def _quadratic(trial):
    value = float(trial[0])
    return np.array([value * value + 0.3 * value])


def test_difference_step_matches_adr_0017():
    """步长为 √ε·max(|x|, s)，与 ADR 0017 第 6 项一致。"""
    samples = ((0.25, 1.0), (4.0, 1.0), (0.0, 3.5), (1.0e-20, 1.197))
    for value, scale in samples:
        point = np.array([value])
        scales = np.array([scale])
        step = _SQRT_EPS * max(abs(value), scale)
        base = _quadratic(point)
        matrix = _jacobian(point, base, _quadratic, scales, lambda trial: True)
        moved = point.copy()
        moved[0] += step
        expected = (_quadratic(moved)[0] - base[0]) / step
        assert matrix[0, 0] == expected
        assert _difference_basis(value, scale) == max(abs(value), scale)


def test_column_scale_is_max_of_abs_value_and_scale():
    """列缩放 C = diag(max(|x|, s))。|x| 小于 s 时取 s。"""
    column = _column_scale(np.array([1.0e-8, 2.0, 0.0]), np.array([1.0, 0.5, 3.5]))
    assert column.tolist() == [1.0, 2.0, 3.5]


def test_scaled_newton_step_matches_the_unscaled_step_when_well_conditioned():
    """良态线性问题。行列缩放不改变牛顿步。含 |x| < s 的列。"""
    _assert_scaled_step(
        np.array([1.5, 0.5, 2.5]), np.array([1.0, 1.0, 1.0]),
        np.array([[2.0, 0.1, 0.0], [0.1, 1.5, 0.2], [0.0, 0.2, 3.0]]),
        np.array([1.0, 2.0, 3.0]),
        ResidualBlocks(np.array([0]), np.array([1]), np.array([2]), 4.0, 2.0, 8.0),
    )
    _assert_scaled_step(
        np.array([1.0e-8, 2.0]), np.array([1.0, 0.5]),
        np.array([[1.2, 0.1], [0.05, 0.8]]), np.array([0.3, 0.4]),
        ResidualBlocks(np.array([0]), np.array([1]), np.array([], dtype=int), 4.0, 2.0, 8.0),
    )


def _assert_scaled_step(state, scales, matrix_a, target, blocks):
    def residual(trial):
        return matrix_a @ np.asarray(trial, dtype=float) - target

    base = residual(state)
    matrix = _jacobian(state, base, residual, scales, lambda trial: True)
    plain = np.linalg.solve(matrix, -base)
    scaled = _newton_step(matrix, base, state, scales, blocks)
    assert np.linalg.norm(scaled - plain) / np.linalg.norm(plain) < 1.0e-10
