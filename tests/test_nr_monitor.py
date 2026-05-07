import numpy as np
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import _nr_monitor
from src.core.reactor import Reactor, ReactorConfig


def test_build_vorab_x0_sanity_forwards_major_gibbs_solver_mode(monkeypatch):
    cfg = ReactorConfig(n_cells=1)
    cfg.thesis_mode = True
    cfg.major_gibbs_solver_mode = "shadow_compare"
    reactor = Reactor(cfg)

    captured: dict[str, str] = {}

    def _fake_t_profile(**kwargs):
        return np.array([float(kwargs["T_inlet"])], dtype=float)

    def _fake_generate_initial_x0(**kwargs):
        captured["major_gibbs_solver_mode"] = str(kwargs.get("major_gibbs_solver_mode"))

    monkeypatch.setattr(_nr_monitor, "estimate_axial_T_profile", _fake_t_profile)
    monkeypatch.setattr(_nr_monitor, "generate_initial_x0", _fake_generate_initial_x0)

    _nr_monitor.build_vorab_x0_sanity(reactor)
    assert captured["major_gibbs_solver_mode"] == "shadow_compare"
