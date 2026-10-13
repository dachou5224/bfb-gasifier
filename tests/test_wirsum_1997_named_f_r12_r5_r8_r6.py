"""F_E0+R12+R5+R8 根上接入 R6。不改速率式，不宣称炉温。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.solvers.wirsum_1997_named_f import coupled_residual_named
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import freeze_r12_r5_r8_q, solve_named_r12_r5_then_r8
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import (
    freeze_r12_r5_r8_r6_q,
    r6_analytic_jacobian,
    r6_source_snapshot,
    residual_r12_r5_r8_plus_scaled_r6,
    solve_named_r12_r5_r8_then_r6,
)
from src.solvers.wirsum_1997_startwert_r5 import apply_named_r5_startwert, same_phase_r5_overlap_mol_s
from src.solvers.wirsum_1997_startwert_r6 import (
    apply_named_r6_startwert,
    apply_same_phase_r6_startwert,
    same_phase_r6_overlap_mol_s,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def _vec(**kwargs: float) -> np.ndarray:
    n = np.zeros(N_GAS, dtype=np.float64)
    for name, value in kwargs.items():
        n[GAS_SPECIES_INDEX[name]] = value
    return n


def test_r6_startwert_consumes_limiting_ch4() -> None:
    n = _vec(CH4=1.0, O2=3.0, CO=0.2, H2O=0.1)
    diag = apply_same_phase_r6_startwert(n)
    assert diag["xi6_mol_s"] == pytest.approx(1.0)
    assert n[GAS_SPECIES_INDEX["CH4"]] == pytest.approx(0.0)
    assert n[GAS_SPECIES_INDEX["O2"]] == pytest.approx(1.5)
    assert n[GAS_SPECIES_INDEX["CO"]] == pytest.approx(1.2)
    assert n[GAS_SPECIES_INDEX["H2O"]] == pytest.approx(2.1)
    assert same_phase_r6_overlap_mol_s(n) == pytest.approx(0.0)


def test_r6_startwert_consumes_limiting_o2() -> None:
    n = _vec(CH4=2.0, O2=1.5, CO=0.0, H2O=0.0)
    apply_same_phase_r6_startwert(n)
    assert n[GAS_SPECIES_INDEX["CH4"]] == pytest.approx(1.0)
    assert n[GAS_SPECIES_INDEX["O2"]] == pytest.approx(0.0)
    assert n[GAS_SPECIES_INDEX["CO"]] == pytest.approx(1.0)
    assert n[GAS_SPECIES_INDEX["H2O"]] == pytest.approx(2.0)
    assert same_phase_r6_overlap_mol_s(n) == pytest.approx(0.0)


def test_r6_startwert_conserves_c_h_o_atoms() -> None:
    n = _vec(CH4=0.8, O2=2.0, CO=0.3, CO2=1.0, H2=0.1, H2O=0.4)
    c_before = n[GAS_SPECIES_INDEX["CH4"]] + n[GAS_SPECIES_INDEX["CO"]] + n[GAS_SPECIES_INDEX["CO2"]]
    h_before = (
        4.0 * n[GAS_SPECIES_INDEX["CH4"]]
        + 2.0 * n[GAS_SPECIES_INDEX["H2"]]
        + 2.0 * n[GAS_SPECIES_INDEX["H2O"]]
    )
    o_before = (
        n[GAS_SPECIES_INDEX["CO"]]
        + 2.0 * n[GAS_SPECIES_INDEX["CO2"]]
        + n[GAS_SPECIES_INDEX["H2O"]]
        + 2.0 * n[GAS_SPECIES_INDEX["O2"]]
    )
    apply_same_phase_r6_startwert(n)
    c_after = n[GAS_SPECIES_INDEX["CH4"]] + n[GAS_SPECIES_INDEX["CO"]] + n[GAS_SPECIES_INDEX["CO2"]]
    h_after = (
        4.0 * n[GAS_SPECIES_INDEX["CH4"]]
        + 2.0 * n[GAS_SPECIES_INDEX["H2"]]
        + 2.0 * n[GAS_SPECIES_INDEX["H2O"]]
    )
    o_after = (
        n[GAS_SPECIES_INDEX["CO"]]
        + 2.0 * n[GAS_SPECIES_INDEX["CO2"]]
        + n[GAS_SPECIES_INDEX["H2O"]]
        + 2.0 * n[GAS_SPECIES_INDEX["O2"]]
    )
    assert c_after == pytest.approx(c_before)
    assert h_after == pytest.approx(h_before)
    assert o_after == pytest.approx(o_before)
    assert n[GAS_SPECIES_INDEX["H2"]] == pytest.approx(0.1)


def test_r6_source_at_r8_root_is_hard_both_phases() -> None:
    r8 = solve_named_r12_r5_then_r8(max_iter=20)
    nq8 = freeze_r12_r5_r8_q()
    nq6 = freeze_r12_r5_r8_r6_q()
    snap = r6_source_snapshot(r8["x"], nq8, nq6)
    assert snap["N_d_CH4"] > 0.0
    assert snap["N_d_O2"] > 0.0
    assert snap["N_b_CH4"] > 0.0
    assert snap["N_b_O2"] > 0.0
    assert abs(float(snap["dF_CH4_d"])) > 100.0
    assert abs(float(snap["dF_CH4_b"])) > 100.0
    assert abs(float(snap["dF_energy"])) < 1.0e-6
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r6_then_r5_startwert_drops_source_at_r8_root() -> None:
    r8 = solve_named_r12_r5_then_r8(max_iter=20)
    nq8 = freeze_r12_r5_r8_q()
    nq6 = freeze_r12_r5_r8_r6_q()
    before = r6_source_snapshot(r8["x"], nq8, nq6)
    packed6 = apply_named_r6_startwert(r8["x"])
    packed5 = apply_named_r5_startwert(packed6["x"])
    after = r6_source_snapshot(packed5["x"], nq8, nq6)
    assert packed6["dense"]["overlap_after_mol_s"] == pytest.approx(0.0, abs=1e-15)
    assert packed6["bubble"]["overlap_after_mol_s"] == pytest.approx(0.0, abs=1e-15)
    assert same_phase_r5_overlap_mol_s(packed5["x"][0:7]) == pytest.approx(0.0, abs=1e-15)
    assert same_phase_r5_overlap_mol_s(packed5["x"][7:14]) == pytest.approx(0.0, abs=1e-15)
    assert after["dF_l2"] < 0.01 * before["dF_l2"]
    J = r6_analytic_jacobian(packed5["x"], nq6, rate_scale=1.0)
    o2_d = float(packed5["x"][IDX["O2"]])
    o2_b = float(packed5["x"][7 + IDX["O2"]])
    # CH4=0 且同相仍有 O2 时，右导数必须非零。两反应物都被后续 R5 Startwert
    # 抽干时，R6 对该相 CH4 列可以为 0。
    if o2_d > 1.0e-12:
        assert float(np.linalg.norm(J[:, IDX["CH4"]])) > 0.0
    if o2_b > 1.0e-12:
        assert float(np.linalg.norm(J[:, 7 + IDX["CH4"]])) > 0.0


def test_r6_scale_zero_matches_r8_residual() -> None:
    r8 = solve_named_r12_r5_then_r8(max_iter=20)
    nq8 = freeze_r12_r5_r8_q()
    nq6 = freeze_r12_r5_r8_r6_q()
    f0 = residual_r12_r5_r8_plus_scaled_r6(r8["x"], nq8, nq6, r6_scale=0.0)
    f8 = coupled_residual_named(r8["x"], nq8, rate_scale=1.0)
    assert np.allclose(f0, f8)
    J = r6_analytic_jacobian(r8["x"], nq6, rate_scale=1.0)
    assert J.shape == (40, 40)
    assert float(np.linalg.norm(J[0:7, 0:7])) > 0.0
    assert float(np.linalg.norm(J[7:14, 7:14])) > 0.0


def test_r12_r5_r8_r6_without_startwert_still_not_furnace() -> None:
    sol = solve_named_r12_r5_r8_then_r6(max_iter=20, apply_startwert=False)
    assert sol["apply_startwert"] is False
    snap = sol["r6_at_r8_root"]
    assert abs(float(snap["dF_CH4_d"])) > 100.0
    hist = sol["continuation"]
    assert hist[0]["r6_scale"] == 0.0
    tiny = next(row for row in hist if abs(row["r6_scale"] - 1.0e-8) < 1.0e-18)
    assert tiny["residual_l2"] < 1.0e-6
    last = hist[-1]
    assert last["r6_scale"] == 1.0
    assert sol["physical_convergence"] is True
    assert sol["not_a_furnace_convergence_claim"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_r12_r5_r8_r6_startwert_solve_does_not_change_rate_or_claim_furnace() -> None:
    sol = solve_named_r12_r5_r8_then_r6(max_iter=20, apply_startwert=True)
    assert sol["rate_formula_changed"] is False
    assert sol["enabled_rxns"] == ("R12", "R5", "R8", "R6")
    assert sol["scaled_linear_solve"] is True
    assert sol["apply_startwert"] is True
    assert sol["startwert_r6"] is not None
    assert sol["startwert_r5_after_r6"] is not None
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["r8_physical_convergence"] is True
    after = sol["r6_after_startwert"]
    before = sol["r6_at_r8_root"]
    assert after["dF_l2"] < 0.01 * before["dF_l2"]
    hist = sol["continuation"]
    if sol.get("startwert_abandoned"):
        assert hist[0]["r6_scale"] == 0.0
        assert hist[-1]["r6_scale"] == 1.0
    else:
        assert hist[0]["r6_scale"] == 1.0e-4
        assert hist[-1]["r6_scale"] == 1.0
        assert hist[-1]["failure_reason"] == "newton_merit"
    assert hist[-1]["equilibrated_rank"] == 40
    assert sol["physical_convergence"] is True
    # CH4 氧化热进这一格。仍不是 HTW。
    assert float(sol["T_K"]) > float(sol["T_r8_K"]) + 50.0
    assert 1200.0 < float(sol["T_K"]) < 1500.0
    assert float(sol["x"][IDX["CH4"]]) < 1.0e-4
    assert float(sol["x"][7 + IDX["CH4"]]) < 1.0e-4
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
