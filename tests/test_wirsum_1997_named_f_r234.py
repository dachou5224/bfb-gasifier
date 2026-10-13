"""均相根上接入 R2/R3/R4。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import coupled_residual_named
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import freeze_r12_r5_r8_r6_r7_q
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import solve_named_r12_r5_r8_r6_then_r7
from src.solvers.wirsum_1997_named_f_r234 import (
    freeze_homog_r234_q,
    r234_analytic_jacobian,
    r234_source_snapshot,
    residual_homog_plus_scaled_r234,
    solve_named_homog_then_r234,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def test_r234_source_at_homog_root_is_dense_char() -> None:
    homog = solve_named_r12_r5_r8_r6_then_r7(max_iter=20)
    nqh = freeze_r12_r5_r8_r6_r7_q()
    nq = freeze_homog_r234_q()
    snap = r234_source_snapshot(homog["x"], nqh, nq)
    assert snap["T_K"] > 1200.0
    # R2 主导：密相产 H2/CO、耗 H2O；气泡无炭。
    assert snap["dF_l2"] > 0.5
    assert snap["dF_H2_d"] > 0.5
    assert snap["dF_CO_d"] > 0.5
    assert snap["dF_H2O_d"] < 0.0
    assert snap["dF_M_fine"] < 0.0
    assert snap["dF_M_coarse"] < 0.0
    assert abs(float(snap["dF_energy"])) < 1.0e-6
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r234_scale_zero_matches_homog_residual() -> None:
    homog = solve_named_r12_r5_r8_r6_then_r7(max_iter=20)
    nqh = freeze_r12_r5_r8_r6_r7_q()
    nq = freeze_homog_r234_q()
    f0 = residual_homog_plus_scaled_r234(homog["x"], nqh, nq, r234_scale=0.0)
    fh = coupled_residual_named(homog["x"], nqh, rate_scale=1.0)
    assert np.allclose(f0, fh)
    J = r234_analytic_jacobian(homog["x"], nq, rate_scale=1.0)
    assert J.shape == (40, 40)
    assert float(np.linalg.norm(J[0:7, 0:7])) > 0.0
    assert float(np.linalg.norm(J[7:14, 7:14])) == 0.0
    assert float(np.linalg.norm(J[14:16, 14:16])) > 0.0
    z = np.asarray(homog["x"], dtype=np.float64).copy()
    z[IDX["H2"]] = 0.0
    J0 = r234_analytic_jacobian(z, nq, rate_scale=1.0)
    # H2=0：R3 右导数必须非零。
    assert float(np.linalg.norm(J0[:, IDX["H2"]])) > 0.0


def test_homog_r234_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_homog_then_r234(max_iter=20)
    assert sol["rate_formula_changed"] is False
    assert sol["enabled_rxns"] == ("R12", "R5", "R8", "R6", "R7", "R2", "R3", "R4")
    assert sol["scaled_linear_solve"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["homog_physical_convergence"] is True
    hist = sol["continuation"]
    assert hist[0]["r234_scale"] == 0.0
    assert hist[0]["failure_reason"] == "newton_merit"
    assert hist[-1]["r234_scale"] == 1.0
    assert hist[-1]["equilibrated_rank"] == 40
    assert hist[-1]["failure_reason"] == "newton_merit"
    assert sol["physical_convergence"] is True
    # R2 产的 H2 仍被剩余 O2 经 R12 烧掉。T 几乎不动。不是 HTW。
    assert abs(float(sol["T_K"]) - float(sol["T_homog_K"])) < 10.0
    assert float(sol["x"][IDX["H2"]]) < 1.0e-8
    assert float(sol["x"][14]) < float(sol["r234_at_homog_root"]["M_fine"])
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
