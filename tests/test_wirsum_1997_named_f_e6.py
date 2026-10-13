"""F_E6 焦油出口未知量。不改速率式，不宣称炉温。"""

from __future__ import annotations

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f_e5 import freeze_e5_q, solve_named_e5_from_e4
from src.solvers.wirsum_1997_named_f_e6 import (
    e6_source_snapshot,
    freeze_e6_q,
    solve_named_e6_from_e5,
)
from src.solvers.wirsum_1997_named_f_registry import omitted_named

_E5_CACHE: dict | None = None


def _e5() -> dict:
    global _E5_CACHE
    if _E5_CACHE is None:
        _E5_CACHE = solve_named_e5_from_e4(Q_W=0.0, max_iter=20)
    return _E5_CACHE


def test_e6_embed_is_already_root() -> None:
    e5 = _e5()
    snap = e6_source_snapshot(e5["x"], freeze_e5_q(Q_W=0.0), freeze_e6_q(Q_W=0.0))
    assert abs(float(snap["F6_tar1"])) < 1.0e-15
    assert abs(float(snap["F6_tar2"])) < 1.0e-15
    assert abs(float(snap["d_energy"])) < 1.0e-8
    assert float(snap["F6_l2"]) < 1.0e-6
    assert float(snap["N_tar1"]) > 0.0
    assert "tar_transport" not in omitted_named(6)
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_e6_solve_does_not_change_rate_or_claim_furnace() -> None:
    e5 = _e5()
    sol = solve_named_e6_from_e5(Q_W=0.0, max_iter=20, e5=e5)
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E6"
    assert sol["n_unknowns"] == 43
    assert sol["r10_tar_oxidation"] is False
    assert sol["Q_W_is_fitted_to_1100K"] is False
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["e5_physical_convergence"] is True
    assert sol["physical_convergence"] is True
    assert sol["residual_l2"] < 1.0e-6
    assert sol["equilibrated_rank"] == 43
    assert abs(float(sol["T_K"]) - float(sol["T_e5_K"])) < 1.0e-6
    assert float(sol["N_tar1"]) > 0.0
    assert abs(float(sol["T_K"]) - 1100.0) > 50.0
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
