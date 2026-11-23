"""T26v 验证入口：硬墙钟记 TIMEOUT 并以非零状态退出（CHARTER §7）。不跑长工况。"""
import time

from hamel.solver.newton import SolveTimeout


def test_blocking_solve_is_interrupted_and_recorded_as_timeout(tmp_path, monkeypatch):
    from tools import t26v_validate as mod

    monkeypatch.setattr(mod, "_SCRATCH", tmp_path)

    def blocking(*args, **kwargs):
        time.sleep(2.0)  # 不轮询 deadline 的阻塞替身
        raise AssertionError("硬墙钟没有中断求解")

    monkeypatch.setattr(mod.driver, "solve_case", blocking)
    summary = mod.run("baseline", False, 0.2)
    assert summary["status"] == "TIMEOUT"
    assert "墙钟上限" in summary["reason"]
    # 中断后返回正式实现，不把替身泄漏给后续调用。
    assert mod.driver.evaluate_residual.__name__ == "evaluate_residual"
    assert mod.initial.evaluate_residual.__name__ == "evaluate_residual"


def test_main_exits_nonzero_on_timeout(monkeypatch):
    from tools import t26v_validate as mod

    monkeypatch.setattr(mod, "run", lambda label, truncate, limit: {"status": "TIMEOUT", "reason": "墙钟"})
    assert mod.main(["baseline", "--limit-s", "1"]) == 2


def test_solve_timeout_is_also_recorded_as_timeout(tmp_path, monkeypatch):
    from tools import t26v_validate as mod

    monkeypatch.setattr(mod, "_SCRATCH", tmp_path)

    def raising(*args, **kwargs):
        raise SolveTimeout("超过期限")

    monkeypatch.setattr(mod.driver, "solve_case", raising)
    assert mod.run("s20", True, 5.0)["status"] == "TIMEOUT"


def test_exception_is_not_reported_as_timeout(tmp_path, monkeypatch):
    from tools import t26v_validate as mod

    monkeypatch.setattr(mod, "_SCRATCH", tmp_path)

    def failing(*args, **kwargs):
        raise RuntimeError("求解实现出错")

    monkeypatch.setattr(mod.driver, "solve_case", failing)
    summary = mod.run("baseline", False, 5.0)
    assert summary["status"] == "EXCEPTION" and "求解实现出错" in summary["reason"]
