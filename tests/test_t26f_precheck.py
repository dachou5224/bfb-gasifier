"""预检的监督入口和由残差项作的除法。不跑长求解。"""
import json
import subprocess

import numpy as np

from tools import t26f_precheck as mod


def test_frozen_rate_estimate_equals_the_inventory_at_a_linear_root():
    """残差为 0 且消耗与存量成正比时，冻结比速率估算就是当前 m_C。"""
    mass = 2.0e-6
    consumption = 4.0e-4
    row = mod._hover_terms(mass, consumption, 0.0, 0.0, 0.0, 0.0)
    assert row["specific_consumption_per_s"] == consumption / mass
    assert row["inflow_kg_s"] == consumption
    assert row["frozen_rate_estimate_kg"] == mass
    assert mod._hover_terms(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)["frozen_rate_estimate_kg"] is None


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


def test_frozen_rate_estimate_is_not_the_root_when_the_specific_rate_varies():
    """消耗 C(M) = M + 9（灰使比速率随存量变化），入流 20。真实根 M = 11，估算从 M = 1 给出 2。"""
    def consumption(m):
        return m + 9.0
    current = 1.0
    inflow = 20.0
    row = mod._hover_terms(current, consumption(current), 0.0, inflow - consumption(current), 0.0, 0.0)
    assert row["inflow_kg_s"] == inflow
    assert row["frozen_rate_estimate_kg"] == 2.0
    root = 11.0
    assert inflow - consumption(root) == 0.0
    assert inflow - consumption(row["frozen_rate_estimate_kg"]) != 0.0


class _Frozen:
    n_bed = 3
    n_freeboard = 2
    n_class = 2
    geometric_mass_kg = np.array([1.0, 2.0, 4.0])
    k_auf_per_s = np.array([0.5, 0.25, 0.125])
    m_dot_ab_kg_s = np.array([0.3, 0.2, 0.0])
    b_free_per_s = np.array([[0.0, 0.0], [1.0, 2.0]])
    b_hover_per_s = np.zeros((2, 2))
    ejection_to_free = np.array([1.0, 1.0])
    ejection_to_hover = np.array([0.0, 0.0])


class _State:
    bed_char_kg = np.array([[0.2, 0.3], [0.4, 0.6], [1.0, 1.0]])
    freeboard_char_free_kg = np.array([[0.0, 0.0], [0.01, 0.02]])


def test_interface_pairs_the_upward_flow_of_the_lower_cell_with_the_downward_flow_of_the_upper_cell():
    face = mod._interface(_State(), _Frozen())
    first, second, surface = face["faces"]
    assert first["face"] == "床区 0/1"
    assert first["up_kg_s"] == 0.25 * 2.0
    assert first["down_kg_s"] == 0.3
    assert first["net_up_kg_s"] == 0.25 * 2.0 - 0.3
    assert first["char_up_kg_s"] == [0.25 * 0.4, 0.25 * 0.6]
    assert first["char_down_kg_s"] == [0.3 * 0.2 / 1.0, 0.3 * 0.3 / 1.0]
    assert second["up_kg_s"] == 0.125 * 4.0
    assert second["down_kg_s"] == 0.2
    assert surface["face"] == "床面"
    assert surface["up_kg_s"] == 0.5 * 1.0
    assert surface["char_down_kg_s"] == [1.0 * 0.01, 2.0 * 0.02]
    assert face["m_dot_ab_0_kg_s"] == 0.3
