"""T26e 验证脚本的监督入口与奇异报告。不重跑长工况。"""
import json
import subprocess

import numpy as np


def test_point_worker_prints_json_and_singular_direction_is_reported(tmp_path, monkeypatch, capsys):
    from tools import t26e_validate as mod

    monkeypatch.setattr(mod, "_OUT", tmp_path)
    monkeypatch.setattr(mod, "evaluate_residual", lambda trial, frozen: type("Residual", (), {"values": np.zeros_like(trial)})())

    def fake_column(x, residual, evaluate, scales, feasible, blocks, names, kinds):
        unresolved = [{"name": names[0], "x": float(x[0]), "step": 1.0}]
        return np.zeros(x.size, dtype=bool), np.eye(x.size), None, unresolved, [{"name": names[1], "x": float(x[1]), "start": 1.0, "step": 2.0, "cap": 3.0}]

    monkeypatch.setattr(mod, "_column", fake_column)
    mod.point()
    summary = json.loads(capsys.readouterr().out)
    assert summary["step_norm"] is None
    assert summary["unresolved"] == 1
    text = (tmp_path / "validation.md").read_text(encoding="utf-8")
    assert "未形成牛顿步，无方向探针。" in text
    assert "被放大 1 列，到上限仍不可分辨 1 列。" in text
    saved = json.loads((tmp_path / "point.json").read_text(encoding="utf-8"))
    assert saved["direction"] is None
    assert saved["unresolved"][0]["step"] == 1.0


def test_spawn_parses_a_successful_json_summary_and_keeps_stderr(tmp_path, monkeypatch, capsys):
    from tools import t26e_validate as mod

    monkeypatch.setattr(
        mod, "_process",
        lambda args, timeout: subprocess.CompletedProcess(args, 0, '{"step_norm": 3.0, "unresolved": 2}', "残差评价 50\n"),
    )
    folder = tmp_path / "logs"
    mod._spawn("point", 5, folder)
    assert json.loads(capsys.readouterr().out)["step_norm"] == 3.0
    assert (folder / "stderr.log").read_text(encoding="utf-8") == "残差评价 50\n"
