"""F_E5 壁面热 Q_W。具名 q，禁止拟合 1100 K，不宣称炉温。"""

from __future__ import annotations

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f_e4 import freeze_e4_q, solve_named_e4_from_e3
from src.solvers.wirsum_1997_named_f_e5 import (
    NAMED_Q_W_EXAMPLE_W,
    e5_source_snapshot,
    freeze_e5_q,
    solve_named_e5_from_e4,
)
from src.solvers.wirsum_1997_named_f_registry import omitted_named

_E4_CACHE: dict | None = None


def _e4() -> dict:
    global _E4_CACHE
    if _E4_CACHE is None:
        _E4_CACHE = solve_named_e4_from_e3(max_iter=40)
    return _E4_CACHE


def test_e5_qw0_matches_e4_residual() -> None:
    e4 = _e4()
    snap = e5_source_snapshot(e4["x"], freeze_e4_q(), freeze_e5_q(Q_W=0.0))
    assert snap["dF_l2"] < 1.0e-12
    assert abs(float(snap["dF_energy"])) < 1.0e-12
    assert "energy_wall_heat" not in omitted_named(5)
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_e5_energy_shifts_by_minus_qw() -> None:
    e4 = _e4()
    qw = NAMED_Q_W_EXAMPLE_W
    snap = e5_source_snapshot(e4["x"], freeze_e4_q(), freeze_e5_q(Q_W=qw))
    assert float(snap["dF_energy"]) == -qw
    assert abs(float(snap["F5_l2"]) - qw) < 1.0e-6


def test_e5_solve_does_not_fit_or_claim_furnace() -> None:
    e4 = _e4()
    sol0 = solve_named_e5_from_e4(Q_W=0.0, max_iter=20, e4=e4)
    assert sol0["rate_formula_changed"] is False
    assert sol0["f_name"] == "F_E5"
    assert sol0["n_unknowns"] == 41
    assert sol0["Q_W"] == 0.0
    assert sol0["Q_W_is_fitted_to_1100K"] is False
    assert sol0["not_a_furnace_convergence_claim"] is True
    assert sol0["production_hamel_outer_refresh"] is False
    assert sol0["physical_convergence"] is True
    assert sol0["residual_l2"] < 1.0e-6
    assert sol0["equilibrated_rank"] == 41
    assert abs(float(sol0["T_K"]) - float(sol0["T_e4_K"])) < 1.0e-6

    solq = solve_named_e5_from_e4(Q_W=NAMED_Q_W_EXAMPLE_W, max_iter=20, e4=e4)
    assert solq["Q_W"] == NAMED_Q_W_EXAMPLE_W
    assert solq["Q_W_is_fitted_to_1100K"] is False
    assert solq["not_a_furnace_convergence_claim"] is True
    assert solq["physical_convergence"] is True
    assert solq["residual_l2"] < 1.0e-6
    assert solq["equilibrated_rank"] == 41
    # 壁面散热，T 下降。禁止把根写成 1100 K。
    assert float(solq["T_K"]) < float(solq["T_e4_K"]) - 4.0
    assert abs(float(solq["T_K"]) - 1100.0) > 50.0
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
