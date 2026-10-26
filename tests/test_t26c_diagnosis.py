"""T26c 诊断里的线段求根与依赖列。不跑工况。"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np

from hamel.solver.vector import unknown_slots
from test_assembly import make_frozen
from tools.t26c_diagnosis import _char_switches, _sign_switches, dependent_columns, segment_root


def test_segment_root_finds_the_open_interval_and_keeps_a_zero_start():
    assert segment_root(1.0, -1.0) == 0.5
    assert segment_root(0.0, 1.0) == 0.0
    assert segment_root(-1.0, -2.0) is None
    assert segment_root(1.0, 2.0) is None
    assert segment_root(0.0, 0.0) == 0.0


def test_switch_lists_keep_crossings_and_skip_absent_cells():
    start = np.array([[1.0, 1.0], [0.2, -0.1]])
    end = np.array([[1.0, -1.0], [0.2, -0.2]])
    mask = np.array([[True, True], [False, True]])
    rows = _char_switches("自由", start, end, mask)
    assert [item["where"] for item in rows] == ["自由 衡算格 0 k=2", "自由 衡算格 1 k=2"]
    assert rows[0]["h"] == 0.5
    assert rows[1]["h"] is None
    bounds = _sign_switches("气相", np.array([0.2, 0.4]), np.array([-0.2, 0.4]), "Ṅ < 0")
    assert len(bounds) == 1 and bounds[0]["h"] == 0.5


def test_frozen_pickle_roundtrip_accepts_tighter_particle_tolerances(tmp_path):
    from dataclasses import replace

    frozen = make_frozen()
    path = tmp_path / "frozen.pkl"
    with path.open("wb") as handle:
        pickle.dump(frozen, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with path.open("rb") as handle:
        loaded = pickle.load(handle)
    tight = replace(
        loaded,
        particle_tolerance_w=loaded.particle_tolerance_w / 10.0,
        coupling_tolerance_k=loaded.coupling_tolerance_k / 10.0,
        particle_max_iterations=loaded.particle_max_iterations * 2,
        coupling_max_iterations=loaded.coupling_max_iterations * 2,
    )
    assert tight.particle_tolerance_w == frozen.particle_tolerance_w / 10.0
    assert tight.coupling_max_iterations == frozen.coupling_max_iterations * 2


def test_char_equation_depends_on_other_indices():
    frozen = make_frozen()
    slots = unknown_slots(frozen)
    row = next(
        slot for slot in slots
        if slot.zone == "bed" and slot.kind == "char" and slot.balance_index == frozen.n_bed - 1 and slot.size_class == 0
    )
    columns = dependent_columns(frozen, slots, row.index)
    assert row.index in columns
    assert any(
        slots[index].kind == "char" and slots[index].balance_index == row.balance_index and slots[index].size_class == 1
        for index in columns
    )
    assert any(slots[index].species == "O2" and slots[index].phase == "suspension" for index in columns)
    assert any(
        slots[index].kind == "char" and slots[index].balance_index == row.balance_index - 1 and slots[index].size_class == 0
        for index in columns
    )
    assert not any(slots[index].phase == "bubble" for index in columns)
    quiet = make_frozen(m_dot_ab_kg_s=np.zeros(frozen.n_bed))
    assert not any(
        slots[index].balance_index == row.balance_index - 1
        for index in dependent_columns(quiet, slots, row.index)
    )


def test_saved_baseline_uses_the_structural_column_union():
    root = Path(__file__).resolve().parents[1]
    record = json.loads((root / "reports/t26c-diagnosis/baseline/record.json").read_text(encoding="utf-8"))
    rows = [item["index"] for item in record["diagnostics"]["rounds"][0]["newton"][-1]["largest"]]
    with (root / "reports/t26c-diagnosis/baseline/frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    slots = unknown_slots(frozen)
    columns = set()
    for row in rows:
        columns.update(dependent_columns(frozen, slots, row))
    assert len(columns) == 62
    assert any(slots[index].species == "O2" and index not in rows for index in columns)
    assert set(rows) < columns


def test_supervised_run_kills_the_process_group(tmp_path, monkeypatch):
    from dataclasses import replace

    from hamel.solver.campaign import config_for
    import tools.t26c_run as runner

    monkeypatch.setattr(runner, "_OUT", tmp_path)
    monkeypatch.setattr(runner, "config_for", lambda name: replace(config_for(name), wall_timeout_s=0.4))
    monkeypatch.setattr(runner, "_worker_command", lambda name, outer: [sys.executable, "-c", "import time; time.sleep(30)"])
    record = runner.run("baseline", None)
    assert record["status"] == "TIMEOUT"
    assert record["reason"] == "墙钟上限内没有返回"
    saved = json.loads((tmp_path / "baseline" / "record.json").read_text(encoding="utf-8"))
    assert saved["status"] == "TIMEOUT"
