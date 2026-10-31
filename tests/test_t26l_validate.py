"""T26l 验证脚本的期限保护与历史研究工具的适用性守卫。

PR #139 审查 P2 两条：验证入口须有覆盖整个 solve_case 的硬超时；hover_drift_closure
在 ADR 0052 删除预启动之后不能再把零起点当作 A 平衡对照点。
"""

from __future__ import annotations

import pytest

from tools import hover_drift_closure, t26l_validate


def test_timeout_is_not_an_ordinary_exception():
    assert not issubclass(t26l_validate._Timeout, Exception)


def test_solve_timeout_is_reported_as_timeout(monkeypatch, tmp_path):
    """期限在 solve_case 内部抛 BaseException 子类时，报告记 TIMEOUT 而不是 DOMAIN_STOP。"""
    def raise_timeout(*args, **kwargs):
        raise t26l_validate._Timeout()

    monkeypatch.setattr(t26l_validate, "solve_case", raise_timeout)
    monkeypatch.setattr(t26l_validate, "_OUT", tmp_path)
    t26l_validate.main(str(tmp_path))
    assert '"TIMEOUT"' in (tmp_path / "validation.json").read_text(encoding="utf-8")


def test_hover_drift_closure_refuses_on_the_current_model():
    """ADR 0052 删除预启动、ADR 0051 改用 α_s·μ：历史工具必须明确拒绝运行。"""
    with pytest.raises(RuntimeError, match="ADR 0052"):
        hover_drift_closure._capture()
