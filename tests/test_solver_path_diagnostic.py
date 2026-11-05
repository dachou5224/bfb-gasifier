"""第 3 层收敛解检查（#145 审查 P2；ADR 0054 第 2 项第 3 层）。"""

import numpy as np

from hamel.solver.newton import IDENTITY, Transform, TrialRule

_NAMES = ("char", "ash", "gas", "T")


def _check(x, inert, geometric=0.5):
    rule = TrialRule(
        IDENTITY, _NAMES, np.array([0, 1]), np.array([2]), np.array([3], dtype=int),
        ((0,),), np.array([geometric], dtype=float), np.zeros(1), 1.0e-12, 1.0e-12,
    )
    return rule.final(np.asarray(x, dtype=float))


def test_positive_stocks_with_negative_inert_are_rejected():
    # 焦炭 0,6 kg 与几何存量 0,5 kg 给出惰性物 −0,1 kg：不能只看 x 中的负值。
    out = _check([0.6, 0.0, 0.1, 1000.0], -0.1)
    assert not out["ok"]
    assert out["inert_negative_cells"] == [0]
    assert out["stocks_below_tau"] == []


def test_roundoff_stocks_are_reported_but_accepted():
    out = _check([-5e-20, 0.01, 0.1, 1000.0], 0.4, geometric=1.0)
    assert out["ok"]
    assert out["stocks_within_tau"] == {"char": -5e-20}


def test_stock_and_gas_below_tolerance_are_rejected():
    assert not _check([-1e-9, 0.01, 0.1, 1000.0], 0.4, geometric=1.0)["ok"]
    assert not _check([0.0, 0.01, -1e-9, 1000.0], 0.4, geometric=1.0)["ok"]


def test_projection_clamps_only_the_power_columns():
    """o2proj（仅诊断）：只把幂变量 s 截到 0，对数列保持原试探值。ADR 0054 审查记录。"""
    from tools.solver_path_diagnostic import projected_damp
    transform = Transform(np.array([1]), np.array([0]), np.array([0.5]), np.array([1.0]))
    y = np.array([0.01, np.log(1.0e-29)])
    point, _, _ = projected_damp(
        y, np.array([-0.05, -1.0]), 1.0, lambda v: np.zeros(2), lambda v: True,
        lambda r: 0.0, 30, None, transform=transform,
    )
    assert point[0] == 0.0
    assert point[1] == np.log(1.0e-29) - 1.0


def test_outer_coupled_hard_timeout_covers_initial_construction(monkeypatch):
    import time
    from pathlib import Path
    import hamel.solver.precalc as precalc
    import tools.outer_coupled as outer
    original = precalc._drain_terms
    monkeypatch.setattr(outer, "construct", lambda: time.sleep(10))
    started = time.monotonic()
    log = outer.run(2, 1.0, {"corellaX"}, 60.0, Path(".scratch"), total_limit_s=1.0)  # 超时前不写文件
    assert log["stop"] == "TIMEOUT"
    assert time.monotonic() - started < 5.0
    assert precalc._drain_terms is original
