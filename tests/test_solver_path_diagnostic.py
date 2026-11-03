"""限时诊断工具的收敛解检查（#145 审查 P2）。"""

import numpy as np

from tools.solver_path_diagnostic import final_violations


def _check(x, inert):
    names = ["char", "ash", "gas", "T"]
    stocks = np.array([True, True, False, False])
    gas = np.array([False, False, True, False])
    return final_violations(np.asarray(x, dtype=float), names, stocks, gas, np.asarray(inert, dtype=float), 1e-12, 1e-12)


def test_positive_stocks_with_negative_inert_are_rejected():
    # 焦炭为正，但总焦炭超过几何存量，惰性物为负：不能只看 x 中的负值。
    out = _check([0.6, 0.0, 0.1, 1000.0], [-0.1])
    assert not out["ok"]
    assert out["inert_negative_cells"] == [0]
    assert out["stocks_below_tau"] == {}


def test_roundoff_stocks_are_reported_but_accepted():
    out = _check([-5e-20, 0.01, 0.1, 1000.0], [0.4])
    assert out["ok"]
    assert out["stocks_within_tau"] == {"char": -5e-20}


def test_stock_and_gas_below_tolerance_are_rejected():
    assert not _check([-1e-9, 0.01, 0.1, 1000.0], [0.4])["ok"]
    assert not _check([0.0, 0.01, -1e-9, 1000.0], [0.4])["ok"]
