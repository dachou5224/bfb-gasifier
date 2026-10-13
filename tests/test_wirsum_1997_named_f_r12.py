"""仅 R12：Hamel 5.38–5.39 不改式，处理同相重叠与解析雅可比。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import coupled_residual_named, freeze_stage_q
from src.solvers.wirsum_1997_named_f_continuation import certify_named_full_rate_endpoint
from src.solvers.wirsum_1997_named_f_r12 import (
    apply_named_r12_startwert,
    freeze_r12_only_q,
    jacobian_named_r12,
    r12_x0_explosion_snapshot,
    solve_named_r12_only,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def test_r12_x0_explosion_is_same_phase_overlap_not_formula() -> None:
    snap = r12_x0_explosion_snapshot()
    assert snap["x0_l2"] > 1.0e12
    assert snap["N_d_H2"] > 0.0
    assert snap["N_d_O2"] > 0.0
    assert snap["dense_xi12_mol_s"] > 0.0
    assert snap["after_startwert_l2"] < 2.0e5
    assert snap["not_a_furnace_convergence_claim"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r12_analytic_jac_full_rank_after_startwert() -> None:
    nq = freeze_r12_only_q()
    e0 = freeze_stage_q(0)
    from src.solvers.wirsum_1997_named_f import default_x0_named, solve_inner_named

    root = solve_inner_named(e0, max_iter=40)
    packed = apply_named_r12_startwert(root["x"])
    J = jacobian_named_r12(packed["x"], nq, rate_scale=1.0)
    rank = int(np.linalg.matrix_rank(J))
    assert rank == 40
    F = coupled_residual_named(packed["x"], nq)
    # Startwert 改了 N，能量行会动；物种重叠已消，不应再是 10^14。
    assert float(np.linalg.norm(F)) < 1.0e5
    assert abs(float(F[IDX["H2"]])) < 1.0
    assert packed["x"][IDX["H2"]] == 0.0 or packed["x"][IDX["H2"]] < packed["x"][IDX["O2"]]


def test_r12_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_r12_only(max_iter=20)
    assert sol["rate_formula_changed"] is False
    assert sol["enabled_rxns"] == ("R12",)
    assert sol["scaled_linear_solve"] is True
    # 数学根：H2 耗尽流形上 ||F|| 可压到阈值以下。不是 HTW / 1100 K。
    assert sol["physical_convergence"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert abs(float(sol["T_K"]) - 1100.0) > 200.0
    assert float(sol["T_K"]) < 1000.0
    assert float(sol["x"][IDX["H2"]]) < 1.0e-8
    hist = sol["continuation"]
    assert hist[0]["rate_scale"] == 0.0
    h2 = [row["N_d_H2"] for row in hist]
    assert min(h2) < h2[0]
    assert sol["residual_l2"] < 1.0e-6
    cert = certify_named_full_rate_endpoint(
        sol["x"], freeze_r12_only_q(), coupling="frozen_q",
    )
    assert cert["endpoint_certified"] is True
    assert cert["clip_applied"] is False
    assert cert["gas_mol_s_l2"] < 1.0e-6
    assert cert["energy_W_abs"] < 1.0
    full = hist[-1]
    assert full["rate_scale"] == 1.0
    assert full["rank_last"] == 40
    assert full["failure_reason"] == "newton_merit"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r12_perturbed_x0_stays_on_same_branch() -> None:
    sol = solve_named_r12_only(max_iter=20)
    x_p = np.asarray(sol["x"], dtype=np.float64).copy()
    x_p[-1] = float(x_p[-1]) + 0.5
    sol_p = solve_named_r12_only(
        x0=x_p, from_e0_root=False, continuation=False,
        apply_startwert=False, max_iter=20,
    )
    cert = certify_named_full_rate_endpoint(
        sol_p["x"], freeze_r12_only_q(), coupling="frozen_q",
    )
    assert cert["endpoint_certified"] is True
    dT = abs(float(sol_p["T_K"]) - float(sol["T_K"]))
    dH2 = abs(float(sol_p["x"][IDX["H2"]]) - float(sol["x"][IDX["H2"]]))
    # 不同根则报分支；不以离 1100 K 判对错。
    assert dT < 2.0, f"r12 branch: dT={dT} K"
    assert dH2 < 1.0e-6, f"r12 branch: dH2={dH2} mol/s"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    nq = freeze_r12_only_q()
    sol = solve_named_r12_only(max_iter=12)
    x = sol["x"]
    F = coupled_residual_named(x, nq, rate_scale=1.0)
    J = jacobian_named_r12(x, nq, rate_scale=1.0)
    from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank, named_r12_newton_step

    info = named_r12_scaled_rank(J, F, x)
    assert info["raw_rank"] < 40
    assert info["equilibrated_rank"] == 40
    dx, rank, _cond = named_r12_newton_step(J, F, x)
    assert rank == 40
    assert dx.shape == x.shape
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
