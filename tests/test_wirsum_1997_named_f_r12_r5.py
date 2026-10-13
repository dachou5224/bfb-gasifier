"""F_E0+R12 根上接入 R5。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.wirsum_1997_named_f import coupled_residual_named
from src.solvers.wirsum_1997_named_f_r12 import freeze_r12_only_q, solve_named_r12_only
from src.solvers.wirsum_1997_named_f_r12_r5 import (
    freeze_r12_r5_q,
    r5_analytic_jacobian,
    r5_source_snapshot,
    residual_r12_plus_scaled_r5,
    solve_named_r12_then_r5,
)
from src.solvers.wirsum_1997_startwert_r5 import (
    apply_named_r5_startwert,
    apply_same_phase_r5_startwert,
    same_phase_r5_overlap_mol_s,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def _vec(**kwargs: float) -> np.ndarray:
    n = np.zeros(N_GAS, dtype=np.float64)
    for name, value in kwargs.items():
        n[GAS_SPECIES_INDEX[name]] = value
    return n


def test_r5_startwert_consumes_limiting_co() -> None:
    n = _vec(CO=1.0, O2=1.0, CO2=0.2)
    diag = apply_same_phase_r5_startwert(n)
    assert diag["xi5_mol_s"] == pytest.approx(1.0)
    assert n[GAS_SPECIES_INDEX["CO"]] == pytest.approx(0.0)
    assert n[GAS_SPECIES_INDEX["O2"]] == pytest.approx(0.5)
    assert n[GAS_SPECIES_INDEX["CO2"]] == pytest.approx(1.2)
    assert same_phase_r5_overlap_mol_s(n) == pytest.approx(0.0)


def test_r5_startwert_consumes_limiting_o2() -> None:
    n = _vec(CO=4.0, O2=1.0, CO2=0.0)
    apply_same_phase_r5_startwert(n)
    assert n[GAS_SPECIES_INDEX["CO"]] == pytest.approx(2.0)
    assert n[GAS_SPECIES_INDEX["O2"]] == pytest.approx(0.0)
    assert n[GAS_SPECIES_INDEX["CO2"]] == pytest.approx(2.0)
    assert same_phase_r5_overlap_mol_s(n) == pytest.approx(0.0)


def test_r5_startwert_conserves_c_and_o_atoms() -> None:
    n = _vec(CO=1.25, O2=5.94, CO2=3.0, H2=0.4, H2O=1.0)
    c_before = n[GAS_SPECIES_INDEX["CO"]] + n[GAS_SPECIES_INDEX["CO2"]]
    o_before = (
        n[GAS_SPECIES_INDEX["CO"]]
        + 2.0 * n[GAS_SPECIES_INDEX["CO2"]]
        + 2.0 * n[GAS_SPECIES_INDEX["O2"]]
        + n[GAS_SPECIES_INDEX["H2O"]]
    )
    apply_same_phase_r5_startwert(n)
    c_after = n[GAS_SPECIES_INDEX["CO"]] + n[GAS_SPECIES_INDEX["CO2"]]
    o_after = (
        n[GAS_SPECIES_INDEX["CO"]]
        + 2.0 * n[GAS_SPECIES_INDEX["CO2"]]
        + 2.0 * n[GAS_SPECIES_INDEX["O2"]]
        + n[GAS_SPECIES_INDEX["H2O"]]
    )
    assert c_after == pytest.approx(c_before)
    assert o_after == pytest.approx(o_before)
    assert n[GAS_SPECIES_INDEX["H2"]] == pytest.approx(0.4)


def test_r5_source_nonzero_at_r12_root() -> None:
    r12 = solve_named_r12_only(max_iter=20)
    nq12 = freeze_r12_only_q()
    nq_both = freeze_r12_r5_q()
    snap = r5_source_snapshot(r12["x"], nq12, nq_both)
    assert snap["N_d_O2"] > 0.0
    assert snap["N_d_CO"] > 0.0
    assert snap["N_d_H2O"] > 0.0
    assert snap["dF_l2"] > 0.0
    assert snap["dF_CO_d"] < 0.0
    assert abs(snap["dF_energy"]) < 1.0e-6
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r5_startwert_drops_bubble_source_at_r12_root() -> None:
    r12 = solve_named_r12_only(max_iter=20)
    nq12 = freeze_r12_only_q()
    nq_both = freeze_r12_r5_q()
    before = r5_source_snapshot(r12["x"], nq12, nq_both)
    packed = apply_named_r5_startwert(r12["x"])
    after = r5_source_snapshot(packed["x"], nq12, nq_both)
    assert abs(float(before["dF_CO_b"])) > 100.0
    assert packed["bubble"]["overlap_after_mol_s"] == pytest.approx(0.0, abs=1e-15)
    assert packed["dense"]["overlap_after_mol_s"] == pytest.approx(0.0, abs=1e-15)
    assert after["dF_l2"] < 0.01 * before["dF_l2"]
    assert abs(after["dF_CO_b"]) < 1.0e-6
    J = r5_analytic_jacobian(packed["x"], nq_both, rate_scale=1.0)
    # CO=0 但 O2/H2O>0：右导数必须非零，否则 Startwert 流形亏秩。
    assert float(np.linalg.norm(J[:, IDX["CO"]])) > 0.0
    assert float(np.linalg.norm(J[:, 7 + IDX["CO"]])) > 0.0


def test_r5_scale_zero_matches_r12_residual() -> None:
    r12 = solve_named_r12_only(max_iter=20)
    nq12 = freeze_r12_only_q()
    nq_both = freeze_r12_r5_q()
    f0 = residual_r12_plus_scaled_r5(r12["x"], nq12, nq_both, r5_scale=0.0)
    f12 = coupled_residual_named(r12["x"], nq12, rate_scale=1.0)
    assert np.allclose(f0, f12)
    J = r5_analytic_jacobian(r12["x"], nq_both, rate_scale=1.0)
    assert J.shape == (40, 40)
    assert float(np.linalg.norm(J)) > 0.0


def test_r12_r5_without_startwert_still_not_furnace() -> None:
    sol = solve_named_r12_then_r5(
        max_iter=20, apply_startwert=False, advance_mode="heuristic_search",
    )
    assert sol["apply_startwert"] is False
    snap = sol["r5_at_r12_root"]
    assert abs(float(snap["dF_CO_b"])) > 100.0 * abs(float(snap["dF_CO_d"]))
    hist = sol["continuation"]
    assert hist[0]["r5_scale"] == 0.0
    tiny = next(row for row in hist if abs(row["r5_scale"] - 1.0e-6) < 1.0e-18)
    assert tiny["failure_reason"] == "newton_merit"
    assert sol["physical_convergence"] is False
    assert sol["not_a_furnace_convergence_claim"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r12_r5_startwert_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_r12_then_r5(max_iter=20, apply_startwert=True)
    assert sol["rate_formula_changed"] is False
    assert sol["enabled_rxns"] == ("R12", "R5")
    assert sol["scaled_linear_solve"] is True
    assert sol["apply_startwert"] is True
    assert sol["advance_mode"] == "root_correcting"
    assert sol["coupling_mode"] == "frozen_q"
    assert sol["endpoint_certified"] is True
    assert sol["startwert"] is not None
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["r12_physical_convergence"] is True
    after = sol["r5_after_startwert"]
    before = sol["r5_at_r12_root"]
    assert after["dF_l2"] < 0.01 * before["dF_l2"]
    hist = sol["continuation"]
    assert hist[0]["r5_scale"] == 1.0e-4
    assert hist[-1]["r5_scale"] == 1.0
    assert hist[-1]["equilibrated_rank"] == 40
    assert hist[-1]["failure_reason"] == "newton_merit"
    # 数学根：烧掉剩余 CO 后 T 升到 ~1100 K 量级。仍是 F_E0+R12+R5，不是 HTW。
    assert sol["physical_convergence"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    assert float(sol["T_K"]) > float(sol["T_r12_K"]) + 200.0
    assert 1000.0 < float(sol["T_K"]) < 1200.0
    assert float(sol["x"][7 + IDX["CO"]]) < float(sol["x"][IDX["CO"]])
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
