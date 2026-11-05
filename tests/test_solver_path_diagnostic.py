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


def _fake_problem():
    from types import SimpleNamespace
    names = [
        "bed[b=1,geo=0].gas_outflow.suspension.O2",
        "bed[b=0,geo=1].gas_outflow.suspension.O2",
        "freeboard[b=0,geo=1].gas_outflow.gas.CO",
    ]
    return SimpleNamespace(
        o2=np.array([0, 1]), names=names, n_ref=1.0,
        blocks=SimpleNamespace(gas=np.array([0, 1, 2])),
        frozen=SimpleNamespace(n_bed=2, secondary_balance_index=5),
    )


def test_o2inlet_power_columns_exclude_logarithmic_o2():
    from tools.solver_path_diagnostic import transforms
    x = np.array([1e-7, 1e-20, 0.3])
    to_x, from_x, scales, pow_idx = transforms(_fake_problem(), {"o2inlet", "gaslog"}, x, np.ones(3))
    assert pow_idx.tolist() == [0]                     # 只有进风格（床区最底格）取幂变量
    y = from_x(x)
    assert y[1] == np.log(1e-20)                       # 非进风格 O₂ 是对数列
    assert np.allclose(to_x(y), x)


def test_projection_keeps_logarithmic_columns():
    from tools.solver_path_diagnostic import _projected_damp
    y = np.array([0.01, np.log(1e-29), np.log(0.3)])
    delta = np.array([-0.05, -1.0, 0.0])
    point, _, _ = _projected_damp(
        y, delta, np.array([0]), 1.0, lambda v: np.zeros(3), lambda v: True, lambda r: 0.0, 30, None,
    )
    assert point[0] == 0.0                              # 幂变量投影到 s = 0
    assert point[1] == np.log(1e-29) - 1.0              # 对数列保持原试探值，不截到 0


def test_outer_coupled_hard_timeout_covers_initial_construction(monkeypatch):
    import time
    from pathlib import Path
    import hamel.solver.precalc as precalc
    import tools.outer_coupled as outer
    original = precalc._drain
    monkeypatch.setattr(outer, "construct", lambda: time.sleep(10))
    started = time.monotonic()
    log = outer.run(2, 1.0, {"gaslog"}, 60.0, Path(".scratch"), total_limit_s=1.0)  # 超时前不写文件
    assert log["stop"] == "TIMEOUT"
    assert time.monotonic() - started < 5.0
    assert precalc._drain is original                   # 替换函数已恢复
