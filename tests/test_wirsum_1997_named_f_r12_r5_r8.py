"""F_E0+R12+R5 根上接入 R8。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import coupled_residual_named
from src.solvers.wirsum_1997_named_f_r12_r5 import freeze_r12_r5_q, solve_named_r12_then_r5
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import (
    freeze_r12_r5_r8_q,
    r8_analytic_jacobian,
    r8_source_snapshot,
    residual_r12_r5_plus_scaled_r8,
    solve_named_r12_r5_then_r8,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def test_r8_source_at_r5_root_is_dense_wgsr() -> None:
    r5 = solve_named_r12_then_r5(max_iter=20, apply_startwert=True)
    nq5 = freeze_r12_r5_q()
    nq8 = freeze_r12_r5_r8_q()
    snap = r8_source_snapshot(r5["x"], nq5, nq8)
    assert snap["T_K"] > 1000.0
    assert snap["dF_l2"] > 0.0
    # 密相产 H2、耗 CO；气泡无灰，R8 源为 0。
    assert snap["dF_H2_d"] > 0.0
    assert snap["dF_CO_d"] < 0.0
    assert abs(snap["dF_H2_b"]) < 1.0e-12
    assert abs(snap["dF_energy"]) < 1.0e-6
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r8_scale_zero_matches_r5_residual() -> None:
    r5 = solve_named_r12_then_r5(max_iter=20, apply_startwert=True)
    nq5 = freeze_r12_r5_q()
    nq8 = freeze_r12_r5_r8_q()
    f0 = residual_r12_r5_plus_scaled_r8(r5["x"], nq5, nq8, r8_scale=0.0)
    f5 = coupled_residual_named(r5["x"], nq5, rate_scale=1.0)
    assert np.allclose(f0, f5)
    J = r8_analytic_jacobian(r5["x"], nq8, rate_scale=1.0)
    assert J.shape == (40, 40)
    assert float(np.linalg.norm(J[0:7, 0:7])) > 0.0
    assert float(np.linalg.norm(J[7:14, 7:14])) == 0.0


def test_r12_r5_r8_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_r12_r5_then_r8(max_iter=20)
    assert sol["rate_formula_changed"] is False
    assert sol["enabled_rxns"] == ("R12", "R5", "R8")
    assert sol["scaled_linear_solve"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["r5_physical_convergence"] is True
    hist = sol["continuation"]
    assert hist[0]["r8_scale"] == 0.0
    assert hist[0]["failure_reason"] == "newton_merit"
    assert hist[-1]["r8_scale"] == 1.0
    assert hist[-1]["equilibrated_rank"] == 40
    assert hist[-1]["failure_reason"] == "newton_merit"
    # 数学根：R8 源 ~0.04 mol/s，H2 仍被 R12 烧掉。T 几乎不动。不是 HTW。
    assert sol["physical_convergence"] is True
    assert abs(float(sol["T_K"]) - float(sol["T_r5_K"])) < 5.0
    assert float(sol["x"][IDX["H2"]]) < 1.0e-8
    assert float(sol["x"][IDX["CO"]]) < float(sol["r8_at_r5_root"]["N_d_CO"])
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
