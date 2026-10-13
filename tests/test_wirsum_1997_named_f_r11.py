"""F_E6+R10+R11 焦油转化。不改 rate_R11_*，不宣称炉温，不升 F_E7。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry
from src.solvers.wirsum_1997_named_f import R10_FUEL_TYPE, default_x0_named
from src.solvers.wirsum_1997_named_f_e6 import solve_named_e6_from_e5
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q, solve_named_e6_then_r10
from src.solvers.wirsum_1997_named_f_r11 import (
    freeze_e6_r10_r11_q,
    r11_analytic_jacobian,
    r11_source_snapshot,
    solve_named_r10_then_r11,
)
from src.solvers.wirsum_1997_named_f_registry import omitted_named
from src.solvers.wirsum_1997_startwert_r11 import (
    apply_same_phase_r11_startwert,
    apply_scaled_bubble_r11_predictor,
)

_E6_CACHE: dict | None = None
_R10_CACHE: dict | None = None


def _e6() -> dict:
    global _E6_CACHE
    if _E6_CACHE is None:
        _E6_CACHE = solve_named_e6_from_e5(Q_W=0.0, max_iter=20)
    return _E6_CACHE


def _r10() -> dict:
    global _R10_CACHE
    if _R10_CACHE is None:
        _R10_CACHE = solve_named_e6_then_r10(Q_W=0.0, max_iter=20, e6=_e6())
    return _R10_CACHE


def _vec(**kwargs: float) -> np.ndarray:
    n = np.zeros(N_GAS, dtype=np.float64)
    for name, value in kwargs.items():
        n[GAS_SPECIES_INDEX[name]] = value
    return n


def test_r11_startwert_consumes_limiting_tar() -> None:
    sto = get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)
    n = _vec(H2O=20.0, CO=0.1, H2=0.1, CH4=0.0)
    diag = apply_same_phase_r11_startwert(n, 0.4, 0.6, sto)
    assert diag["xi11_mol_s"] > 0.0
    assert diag["N_tar1_after"] + diag["N_tar2_after"] < 1.0
    assert n[GAS_SPECIES_INDEX["H2O"]] < 20.0


def test_r11_bubble_predictor_keeps_dense_h2() -> None:
    r10 = _r10()
    x0 = np.asarray(r10["x"], dtype=np.float64)
    packed = apply_scaled_bubble_r11_predictor(x0, ds=0.05)
    z = packed["x"]
    assert packed["xi11_mol_s"] > 0.0
    assert float(z[40] + z[41]) < float(x0[40] + x0[41])
    assert abs(float(z[2]) - float(x0[2])) < 1.0e-15


def test_r11_zero_tar_keeps_finite_right_derivative() -> None:
    nq = freeze_e6_r10_r11_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[40] = 0.0
    x[41] = 0.0
    J = r11_analytic_jacobian(x, nq, rate_scale=1.0)
    assert float(np.linalg.norm(J[:, 40])) > 0.0
    assert float(np.linalg.norm(J[:, 41])) > 0.0


def test_r11_at_r10_root_consumes_tar_and_steam() -> None:
    r10 = _r10()
    snap = r11_source_snapshot(
        r10["x"], freeze_e6_r10_q(Q_W=0.0), freeze_e6_r10_r11_q(Q_W=0.0)
    )
    assert snap["dF_l2"] > 0.0
    assert float(snap["dF_H2O_d"]) < 0.0
    assert float(snap["dF_tar1"]) + float(snap["dF_tar2"]) != 0.0
    assert omitted_named(6)["freeboard_reaction"] == "omitted, named"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r11_solve_does_not_change_rate_or_claim_furnace() -> None:
    r10 = _r10()
    sol = solve_named_r10_then_r11(Q_W=0.0, max_iter=20, r10=r10)
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E6+R10+R11"
    assert sol["n_unknowns"] == 43
    assert sol["r11_tar_conversion"] is True
    assert sol["r10_tar_oxidation"] is True
    assert sol["freeboard_reaction"] is False
    assert sol["Q_W_is_fitted_to_1100K"] is False
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["x_nonneg_projection"] is False
    assert sol["apply_startwert"] is True
    assert sol["advance_mode"] == "root_correcting"
    assert sol["coupling_mode"] == "composite_P"
    assert sol["endpoint_cert"]["clip_applied"] is False
    after = sol["r11_after_startwert"]
    before = sol["r11_at_r10_root"]
    assert after["dF_l2"] < before["dF_l2"]
    if sol["endpoint_certified"]:
        assert sol["residual_l2"] < 1.0e-6
        assert sol["equilibrated_rank"] == 43
        tar = float(sol["N_tar1"]) + float(sol["N_tar2"])
        tar0 = float(r10["N_tar1"]) + float(r10["N_tar2"])
        assert tar <= tar0 + 1.0e-12
    assert abs(float(sol["T_K"]) - 1100.0) > 50.0
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
