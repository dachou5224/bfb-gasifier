"""研究工具 tools/initial_solids_diagnostic.py 的期限保护：回溯中的超时必须向外传播。"""

from types import SimpleNamespace

import numpy as np
import pytest

import tools.initial_solids_diagnostic as tool


def _one_dimensional(monkeypatch, fail_on_call: int):
    calls = {"n": 0}

    def evaluate(x, frozen):
        calls["n"] += 1
        if calls["n"] == fail_on_call:
            raise tool._Timeout()
        return SimpleNamespace(values=np.array([x[0] - 1.0]))

    monkeypatch.setattr(tool, "_consistent", lambda x, frozen, case, recompute: (frozen, None))
    monkeypatch.setattr(tool.residual, "evaluate_residual", evaluate)
    monkeypatch.setattr(tool, "unknown_slots", lambda frozen: [SimpleNamespace(kind="char", index=0)])
    return SimpleNamespace(geometric_mass_kg=np.array([1.0]))


def test_timeout_during_line_search_propagates(monkeypatch):
    # 第 1 次评价是当前点，第 2 次是雅可比列，第 3 次是第一个回溯试探点。
    frozen = _one_dimensional(monkeypatch, fail_on_call=3)
    with pytest.raises(tool._Timeout):
        tool._newton(np.array([0.0]), frozen, None, 1.0, False, limit=1)


def test_timeout_is_not_an_ordinary_exception():
    assert not issubclass(tool._Timeout, Exception)


def test_main_records_timeout_and_returns_nonzero(monkeypatch, tmp_path):
    monkeypatch.setattr(tool, "_OUT", tmp_path)

    def capture(context):
        raise tool._Timeout()

    monkeypatch.setattr(tool, "_capture", capture)
    monkeypatch.setattr("sys.argv", ["initial_solids_diagnostic", "--timeout-s", "60"])
    assert tool.main() == 1
    assert '"TIMEOUT"' in (tmp_path / "study.json").read_text(encoding="utf-8")
