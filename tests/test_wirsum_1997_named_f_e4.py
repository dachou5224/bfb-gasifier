"""F_E4 DAEM X(T)。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import daem_X_at_T
from src.solvers.wirsum_1997_named_f_e3 import freeze_e3_q, solve_named_e3_from_r234
from src.solvers.wirsum_1997_named_f_e4 import (
    apply_energy_t_startwert_e4,
    e4_source_snapshot,
    freeze_e4_q,
    solve_named_e4_from_e3,
)
from src.solvers.wirsum_1997_named_f_registry import omitted_named
from src.solvers.wirsum_1997_stiffness_full_gas import IDX

_E3: dict | None = None


def _e3() -> dict:
    global _E3
    if _E3 is None:
        _E3 = solve_named_e3_from_r234(max_iter=20)
    return _E3


def test_e4_daem_at_e3_root_reassigns_pyrolysis() -> None:
    e3 = _e3()
    snap = e4_source_snapshot(e3["x"], freeze_e3_q(), freeze_e4_q())
    assert snap["X_vm"] > 0.99
    assert snap["dF_l2"] > 100.0
    assert abs(float(snap["dF_energy"])) > 100.0
    assert snap["n_py_sum_daem"] > snap["n_py_sum_frozen"]
    sw = apply_energy_t_startwert_e4(e3["x"], freeze_e3_q(), freeze_e4_q(), daem_scale=1.0)
    assert sw["applied"] is True
    assert abs(float(sw["E_after_W"])) < 1.0
    assert 0.0 < float(sw["dT_K"]) < 20.0
    assert omitted_named(4)["energy_wall_heat"] == "omitted, named"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_e4_x_vm_increases_with_t() -> None:
    nq = freeze_e4_q()
    tau = float(nq.tau_daem_s)
    assert daem_X_at_T(900.0, tau) > daem_X_at_T(700.0, tau)
    assert daem_X_at_T(1312.0, tau) > 0.99


def test_e4_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_e4_from_e3(max_iter=40, e3=_e3())
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E4"
    assert sol["n_unknowns"] == 41
    assert sol["scaled_linear_solve"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["e3_physical_convergence"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["coupling_mode"] == "composite_P"
    assert sol["physical_convergence"] is True
    assert sol["residual_l2"] < 1.0e-6
    assert sol["equilibrated_rank"] == 41
    assert float(sol["X_vm"]) > 0.99
    # 密相出口残氧不够烧掉全部重分配 H2；反向 WGSR 接到平衡。仍不是 HTW。
    assert float(sol["T_K"]) > float(sol["T_e3_K"])
    assert 1200.0 < float(sol["T_K"]) < 1500.0
    assert 0.0 < float(sol["x"][IDX["H2"]]) < 0.02
    assert 0.0 < float(sol["x"][IDX["O2"]]) < 0.015
    assert sol["pyro_reverse_wgsr_startwert"] is True
    assert sol["energy_t_startwert"]["applied"] is True
    assert abs(float(sol["energy_W"])) < 1.0
    assert sol["bound_feasible_direction"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
