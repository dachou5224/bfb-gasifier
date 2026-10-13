"""F_E3 灰库存。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import coupled_residual_named
from src.solvers.wirsum_1997_named_f_e3 import (
    apply_energy_t_startwert,
    embed_e2_root_in_e3,
    freeze_e3_q,
    mash_eq26_kg,
    residual_named_e3,
    solve_named_e3_from_r234,
)
from src.solvers.wirsum_1997_named_f_r234 import freeze_homog_r234_q, solve_named_homog_then_r234
from src.solvers.wirsum_1997_named_f_registry import omitted_named
from src.solvers.wirsum_1997_stiffness_full_gas import IDX

_R234: dict | None = None


def _r234():
    global _R234
    if _R234 is None:
        _R234 = solve_named_homog_then_r234(max_iter=20)
    return _R234


def test_e3_eq26_startwert_zeros_ash_row() -> None:
    r234 = _r234()
    nq3 = freeze_e3_q()
    x41 = embed_e2_root_in_e3(r234["x"], nq3)
    F = coupled_residual_named(x41, nq3)
    assert x41.shape == (41,)
    assert abs(float(F[39])) < 1.0e-12
    assert float(x41[39]) == mash_eq26_kg(nq3)
    assert omitted_named(3)["energy_wall_heat"] == "omitted, named"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_e3_energy_t_startwert_closes_fly_ash_enthalpy() -> None:
    r234 = _r234()
    nq3 = freeze_e3_q()
    x41 = embed_e2_root_in_e3(r234["x"], nq3)
    f0 = residual_named_e3(x41, nq3, rho_scale=0.0)
    assert abs(float(f0[-1])) > 1.0
    assert float(np.linalg.norm(f0[:39])) < 1.0e-8
    sw = apply_energy_t_startwert(x41, nq3)
    assert sw["applied"] is True
    assert abs(float(sw["E_after_W"])) < 1.0
    assert abs(float(sw["dT_K"])) < 1.0


def test_e3_embed_does_not_use_frozen_rho75() -> None:
    r234 = _r234()
    nq2 = freeze_homog_r234_q()
    nq3 = freeze_e3_q()
    x41 = embed_e2_root_in_e3(r234["x"], nq3)
    V_d = float(nq3.energy.species.geo_ex.V_d)
    rho_star = float(x41[39]) / V_d
    assert rho_star < 0.1 * float(nq2.rho_ash)
    F = coupled_residual_named(x41, nq3)
    assert abs(float(F[-1])) > 1.0


def test_e3_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_e3_from_r234(max_iter=20, r234=_r234())
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E3"
    assert sol["n_unknowns"] == 41
    assert sol["scaled_linear_solve"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["r234_physical_convergence"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["physical_convergence"] is True
    assert sol["equilibrated_rank"] == 41
    assert abs(float(sol["M_ash"]) - float(sol["M_ash_eq26"])) < 1.0e-8
    assert abs(float(sol["T_K"]) - float(sol["T_r234_K"])) < 20.0
    assert float(sol["x"][IDX["H2"]]) < 1.0e-8
    assert sol["energy_t_startwert"]["applied"] is True
    assert abs(float(sol["energy_W"])) < 1.0
    assert sol["bound_feasible_direction"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
