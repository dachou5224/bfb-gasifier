"""预检的监督入口和由残差项作的除法。不跑长求解。"""
import json
import subprocess

import numpy as np

from tools import t26f_precheck as mod


def test_root_specific_rate_recovers_the_inventory_from_the_balance_terms():
    """残差为 0 且消耗与存量成正比时，稳态存量就是当前 m_C。"""
    mass = 2.0e-6
    consumption = 4.0e-4
    row = mod._hover_terms(mass, consumption, 0.0, 0.0, 0.0, 0.0)
    assert row["specific_consumption_per_s"] == consumption / mass
    assert row["inflow_kg_s"] == consumption
    assert row["steady_kg"] == mass
    assert mod._hover_terms(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)["steady_kg"] is None


def test_volume_ratio_uses_the_voidage_limit():
    parts = mod._volume_parts(solid_kg=5.0, rho_p=1000.0, cell_volume=2.0, epsilon_mf=0.5)
    assert parts["solid_m3"] == 0.005
    assert parts["limit_m3"] == 1.0
    assert parts["ratio"] == 0.005


def test_supervise_records_timeout_and_the_report_keeps_the_scope_sentence(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mod, "_OUT", tmp_path)
    monkeypatch.setattr(mod, "_process", lambda args, timeout: None)
    record = mod.supervise("initial", 1.0)
    assert record["status"] == "TIMEOUT"
    assert record["reason"] == "墙钟上限内没有返回"
    saved = json.loads((tmp_path / "initial" / "record.json").read_text(encoding="utf-8"))
    assert saved["scope"] == mod._SCOPE
    mod.write_report({"initial": record, "failure": None})
    text = (tmp_path / "precheck.md").read_text(encoding="utf-8")
    assert mod._SCOPE in text
    assert "不能证明有效域内无正解" in text
    assert capsys.readouterr().out == ""


def test_supervise_keeps_stderr_when_the_child_exits(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "_OUT", tmp_path)
    monkeypatch.setattr(
        mod, "_process",
        lambda args, timeout: subprocess.CompletedProcess(args, 1, "", "预计算失败\n"),
    )
    record = mod.supervise("failure", 1.0)
    assert record["status"] == "ERROR"
    assert (tmp_path / "failure" / "stderr.log").read_text(encoding="utf-8") == "预计算失败\n"
    assert np.isfinite(record["supervised_wall_clock_s"])
