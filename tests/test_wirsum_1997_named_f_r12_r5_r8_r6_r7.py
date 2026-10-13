"""F_E0+R12+R5+R8+R6 根上接入 R7。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import coupled_residual_named
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import (
    freeze_r12_r5_r8_r6_q,
    solve_named_r12_r5_r8_then_r6,
)
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import (
    freeze_r12_r5_r8_r6_r7_q,
    r7_analytic_jacobian,
    r7_source_snapshot,
    residual_r12_r5_r8_r6_plus_scaled_r7,
    solve_named_r12_r5_r8_r6_then_r7,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def test_r7_source_at_r6_root_is_tiny() -> None:
    r6 = solve_named_r12_r5_r8_then_r6(max_iter=20, apply_startwert=True)
    nq6 = freeze_r12_r5_r8_r6_q()
    nq7 = freeze_r12_r5_r8_r6_r7_q()
    snap = r7_source_snapshot(r6["x"], nq6, nq7)
    assert snap["T_K"] > 1200.0
    assert snap["N_d_CH4"] < 1.0e-4
    assert snap["dF_l2"] < 1.0e-10
    assert abs(float(snap["dF_energy"])) < 1.0e-6
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r7_scale_zero_matches_r6_residual() -> None:
    r6 = solve_named_r12_r5_r8_then_r6(max_iter=20, apply_startwert=True)
    nq6 = freeze_r12_r5_r8_r6_q()
    nq7 = freeze_r12_r5_r8_r6_r7_q()
    f0 = residual_r12_r5_r8_r6_plus_scaled_r7(r6["x"], nq6, nq7, r7_scale=0.0)
    f6 = coupled_residual_named(r6["x"], nq6, rate_scale=1.0)
    assert np.allclose(f0, f6)
    J = r7_analytic_jacobian(r6["x"], nq7, rate_scale=1.0)
    assert J.shape == (40, 40)
    z = np.asarray(r6["x"], dtype=np.float64).copy()
    z[IDX["CH4"]] = 0.0
    J0 = r7_analytic_jacobian(z, nq7, rate_scale=1.0)
    # CH4=0：右导数必须非零。
    assert float(np.linalg.norm(J0[:, IDX["CH4"]])) > 0.0


def test_r12_r5_r8_r6_r7_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_r12_r5_r8_r6_then_r7(max_iter=20)
    assert sol["rate_formula_changed"] is False
    assert sol["enabled_rxns"] == ("R12", "R5", "R8", "R6", "R7")
    assert sol["scaled_linear_solve"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["r6_physical_convergence"] is True
    hist = sol["continuation"]
    assert hist[0]["r7_scale"] == 0.0
    assert hist[0]["failure_reason"] == "newton_merit"
    assert hist[-1]["r7_scale"] == 1.0
    assert hist[-1]["equilibrated_rank"] == 40
    assert hist[-1]["failure_reason"] == "newton_merit"
    assert sol["physical_convergence"] is True
    assert abs(float(sol["T_K"]) - float(sol["T_r6_K"])) < 1.0
    assert float(sol["x"][IDX["CH4"]]) < 1.0e-4
    assert float(sol["x"][IDX["H2"]]) < 1.0e-8
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
