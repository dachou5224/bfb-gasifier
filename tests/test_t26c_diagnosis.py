"""T26c 诊断里的线段求根。不跑工况。"""
import pickle

import numpy as np

from test_assembly import make_frozen
from tools.t26c_diagnosis import _char_switches, _sign_switches, segment_root


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
