"""F_E6+R10 焦油氧化。不改 rate_R10，不宣称炉温。"""

from __future__ import annotations

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q, solve_named_e6_from_e5
from src.solvers.wirsum_1997_named_f_r10 import (
    freeze_e6_r10_q,
    r10_source_snapshot,
    solve_named_e6_then_r10,
)
from src.solvers.wirsum_1997_named_f_registry import omitted_named

_E6_CACHE: dict | None = None


def _e6() -> dict:
    global _E6_CACHE
    if _E6_CACHE is None:
        _E6_CACHE = solve_named_e6_from_e5(Q_W=0.0, max_iter=20)
    return _E6_CACHE


def test_r10_at_e6_root_is_hard_source() -> None:
    e6 = _e6()
    snap = r10_source_snapshot(e6["x"], freeze_e6_q(Q_W=0.0), freeze_e6_r10_q(Q_W=0.0))
    assert snap["dF_l2"] > 1.0e3
    assert float(snap["dF_O2_d"]) < -1.0
    assert omitted_named(6)["freeboard_reaction"] == "omitted, named"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r10_solve_does_not_change_rate_or_claim_furnace() -> None:
    e6 = _e6()
    sol = solve_named_e6_then_r10(Q_W=0.0, max_iter=20, e6=e6)
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E6+R10"
    assert sol["n_unknowns"] == 43
    assert sol["r10_tar_oxidation"] is True
    assert sol["r11_tar_conversion"] is False
    assert sol["Q_W_is_fitted_to_1100K"] is False
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["e6_physical_convergence"] is True
    assert sol["apply_startwert"] is False
    assert sol["x_nonneg_projection"] is False
    assert sol["advance_mode"] == "heuristic_search"
    assert sol["coupling_mode"] == "composite_P"
    assert sol["endpoint_cert"]["clip_applied"] is False
    assert "solver_stop_reason" in sol
    # 去掉残差 clip 后启发式路径不保证仍到旧端点；认证走未投影 F。
    if sol["endpoint_certified"]:
        assert sol["residual_l2"] < 1.0e-6
        assert sol["equilibrated_rank"] == 43
        assert float(sol["N_tar1"]) + float(sol["N_tar2"]) > 0.0
    assert abs(float(sol["T_K"]) - 1100.0) > 50.0
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
