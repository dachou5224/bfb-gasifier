"""F_E7 自由板均相。不改 rate_*，不宣称炉温，不升 F_E8。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import default_x0_named
from src.solvers.wirsum_1997_named_f_e6 import solve_named_e6_from_e5
from src.solvers.wirsum_1997_named_f_e7 import (
    e7_source_snapshot,
    freeze_e7_q,
    residual_named_e7,
    residual_r11_plus_scaled_fb,
    solve_named_e7_from_r11,
)
from src.solvers.wirsum_1997_named_f_r10 import solve_named_e6_then_r10
from src.solvers.wirsum_1997_named_f_r11 import freeze_e6_r10_r11_q, solve_named_r10_then_r11
from src.solvers.wirsum_1997_named_f_registry import omitted_named
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_startwert_e7 import apply_named_e7_fb_startwert

_E6_CACHE: dict | None = None
_R10_CACHE: dict | None = None
_R11_CACHE: dict | None = None


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


def _r11() -> dict:
    global _R11_CACHE
    if _R11_CACHE is None:
        _R11_CACHE = solve_named_r10_then_r11(Q_W=0.0, max_iter=20, r10=_r10())
    return _R11_CACHE


def _vec(**kwargs: float) -> np.ndarray:
    n = np.zeros(N_MAJ, dtype=np.float64)
    for name, value in kwargs.items():
        n[IDX[name]] = value
    return n


def test_e7_fb_startwert_consumes_mixed_overlap() -> None:
    nq = freeze_e7_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[0:7] = _vec(H2=2.0, CO=0.5, CH4=0.4, N2=8.0)
    x[7:14] = _vec(O2=1.5, N2=3.0)
    x[16:23] = x[0:7] + x[7:14]
    x[23:30] = x[16:23]
    x[30:37] = x[16:23]
    packed = apply_named_e7_fb_startwert(x)
    assert packed["xi12_mol_s"] > 0.0
    n_fb = packed["x"][16:23]
    assert float(n_fb[IDX["O2"]]) < 1.5
    np.testing.assert_allclose(packed["x"][23:30], n_fb)
    np.testing.assert_allclose(packed["x"][0:7], x[0:7])


def test_e7_scale_zero_matches_r11_residual() -> None:
    nq7 = freeze_e7_q(Q_W=0.0)
    nq11 = freeze_e6_r10_r11_q(Q_W=0.0)
    x = default_x0_named(nq7)
    f0 = residual_r11_plus_scaled_fb(x, nq11, nq7, fb_scale=0.0)
    f11 = residual_r11_plus_scaled_fb(x, nq11, nq7, fb_scale=0.0)
    f_named = residual_named_e7(x, nq7)
    np.testing.assert_allclose(f0, f11)
    assert not np.allclose(f0[16:23], f_named[16:23])
    assert omitted_named(7).get("freeboard_reaction") is None
    assert omitted_named(7)["K_eject_in_this_P"] == "omitted, named"
    assert omitted_named(7)["production_hamel_outer_refresh"] == "omitted, named"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_e7_at_r11_root_changes_fb_rows() -> None:
    r11 = _r11()
    nq11 = freeze_e6_r10_r11_q(Q_W=0.0)
    nq7 = freeze_e7_q(Q_W=0.0)
    snap = e7_source_snapshot(r11["x"], nq11, nq7)
    assert snap["dF_fb_l2"] > 1.0
    packed = apply_named_e7_fb_startwert(r11["x"])
    after = e7_source_snapshot(packed["x"], nq11, nq7)
    assert after["dF_l2"] < snap["dF_l2"]
    assert after["src_O2"] == 0.0
    cfg = ReactorConfig()
    assert cfg.vorab_outer_source_policy != "hamel_outer_refresh"


def test_e7_solve_does_not_change_rate_or_claim_furnace() -> None:
    r11 = _r11()
    sol = solve_named_e7_from_r11(Q_W=0.0, max_iter=20, r11=r11)
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E7"
    assert sol["n_unknowns"] == 43
    assert sol["freeboard_reaction"] is True
    assert sol["freeboard_r10_r11"] is False
    assert sol["production_freeboard_step"] is False
    assert sol["r10_tar_oxidation"] is True
    assert sol["r11_tar_conversion"] is True
    assert sol["Q_W_is_fitted_to_1100K"] is False
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["x_nonneg_projection"] is False
    assert sol["apply_startwert"] is True
    assert sol["advance_mode"] == "root_correcting"
    assert sol["coupling_mode"] == "composite_P"
    assert sol["endpoint_cert"]["clip_applied"] is False
    assert sol["K_eject_in_this_P"] == "omitted, named"
    assert sol["endpoint_certified"] is True
    assert sol["endpoint_cert"]["gas_mol_s_l2"] < 1.0e-6
    assert sol["endpoint_cert"]["tar_mol_s_l2"] < 1.0e-6
    assert sol["endpoint_cert"]["energy_W_abs"] < 1.0
    assert sol["equilibrated_rank"] == 43
    assert float(sol["N_fb_O2"]) >= -1.0e-15
    after = sol["e7_after_startwert"]
    before = sol["e7_at_r11_root"]
    assert after["dF_l2"] < before["dF_l2"]
    assert abs(float(sol["T_K"]) - 1100.0) > 50.0
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    assert cfg.vorab_outer_source_policy != "hamel_outer_refresh"
