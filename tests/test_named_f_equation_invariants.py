"""H1–H3 的独立恒等式。预期值用手算，不从被测残差反推。"""

import numpy as np
import pytest

from src.core.constants import Rg
from src.solvers.wirsum_1997_named_f import (
    _conc_major,
    c_tar_components,
    path_char_consumption_kg_s,
    recycle_bins_from_leave,
)


def test_h1_major_plus_tar_exceeds_p_over_rt_by_the_tar_fraction() -> None:
    T = 1000.0  # [K]
    P = 101325.0  # [Pa]
    n_major = np.array([1.0, 0.0, 0.0, 2.0, 0.0, 0.5, 6.5])
    n_tar1, n_tar2 = 0.05, 0.05
    c_star = P / (Rg * T)
    y_tar = 0.10 / (float(np.sum(n_major)) + 0.10)
    c1, c2, c_sum = c_tar_components(n_major, n_tar1, n_tar2, T, P)
    major = _conc_major(n_major, T, P)
    assert float(np.sum(major)) == pytest.approx(c_star)
    assert c_sum == pytest.approx(y_tar * c_star)
    assert c1 + c2 == pytest.approx(c_sum)
    assert float(np.sum(major)) + c_sum == pytest.approx(c_star * (1.0 + y_tar))
    z1, z2, zsum = c_tar_components(n_major, 0.0, 0.0, T, P)
    assert (z1, z2, zsum) == (0.0, 0.0, 0.0)
    assert float(np.sum(_conc_major(n_major, T, P))) == pytest.approx(c_star)


def test_h2_recycle_follows_post_reaction_leave_when_both_classes_remain() -> None:
    # 反应前离开 0.30 与 0.70 kg/s，回料总量 0.40 kg/s。
    # 反应后离开 0.10 与 0.65 kg/s。手算比例 0.10/0.75 与 0.65/0.75。
    post = np.array([0.10, 0.65])
    pre_split = np.array([0.40 * 0.30, 0.40 * 0.70])
    out = recycle_bins_from_leave(post, 0.40, pre_split)
    assert out[0] == pytest.approx(0.40 * 0.10 / 0.75)
    assert out[1] == pytest.approx(0.40 * 0.65 / 0.75)
    assert float(np.sum(out)) == pytest.approx(0.40)
    assert out[0] != pytest.approx(pre_split[0])


def test_h2_zero_leave_does_not_erase_a_defined_inlet_split() -> None:
    inlet = np.array([0.25, 0.15])
    out = recycle_bins_from_leave(np.zeros(2), 0.40, inlet)
    assert out == pytest.approx(inlet)


def test_h3_time_weighted_consumption_can_exceed_the_return_stream() -> None:
    # tau 是 w·t 的和。回落占时间 10/10.1，质量只占 0.20。
    # 消耗 1 kg/s 时，分给回落的是 10/10.1 kg/s，大于 0.20 kg/s。
    on_return, on_leave = path_char_consumption_kg_s(
        np.array([1.0, 0.0]),
        np.array([10.0, 0.0]),
        np.array([0.1, 0.0]),
    )
    assert on_return[0] == pytest.approx(10.0 / 10.1)
    assert on_leave[0] == pytest.approx(0.1 / 10.1)
    m_return = 0.20 * 1.0 - float(on_return[0])
    assert m_return < 0.0
    assert float(on_return[0] + on_leave[0]) == pytest.approx(1.0)
