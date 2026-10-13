"""F_E0 门禁与 F_E2–F_E8 一次一项省略。不改生产默认。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_coupled_f_energy import (
    F_NAME,
    dump_e0_E_of_T_json,
    freeze_energy_q,
    coupled_residual as residual_e0,
    default_x0,
)
from src.solvers.wirsum_1997_named_f import (
    ash_feed_kg_s,
    coupled_residual_named,
    daem_X_at_T,
    default_x0_named,
    freeze_stage_q,
    homog_sources,
    homogeneous_full_rate_gate,
    n_unknowns,
    outer_P_named,
    rho_ash_used,
    stage_payload,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_registry import f_name, omitted_named
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def test_e0_contract_and_scan_json(tmp_path: Path) -> None:
    assert F_NAME == "F_E0"
    path = tmp_path / "named_f_e0_E_of_T.json"
    payload = dump_e0_E_of_T_json(str(path))
    assert payload["f_name"] == "F_E0"
    assert payload["n_unknowns"] == 40
    assert payload["Q_W"] == 0.0
    assert payload["outer_P_updates_X_vm"] is False
    assert payload["not_a_furnace_convergence_claim"] is True
    assert payload["sign_change"] is True
    assert payload["energy_W"][0] > 0.0
    assert payload["energy_W"][1] > 0.0
    assert payload["energy_W"][2] < 0.0
    nq = freeze_stage_q(0)
    x = default_x0(nq.energy)
    np.testing.assert_allclose(
        coupled_residual_named(x, nq), residual_e0(x, nq.energy), atol=1e-12
    )


def test_homogeneous_gate_does_not_claim_physical_root() -> None:
    gate = homogeneous_full_rate_gate()
    assert gate["physical_convergence"] is False
    assert gate["T_retained"] is True
    assert gate["fictitious_Q_W"] is False
    assert gate["full_raw_rank"] <= 16
    assert gate["fe2_may_claim_physical_convergence"] is False


def test_e2_bubble_r8_zero_and_no_physical_claim() -> None:
    nq = freeze_stage_q(2)
    assert f_name(2) == "F_E2"
    assert n_unknowns(2) == 40
    assert "R5_R8_R12" not in omitted_named(2)
    assert omitted_named(2)["ash_inventory"] == "omitted, named"
    x = default_x0_named(nq)
    qT = nq.energy.species
    _sd, sb = homog_sources(x[:39], qT, rate_scale=1.0, rho_ash=75.0)
    # 气泡 R8 关闭：灰催化不进 bubble
    _sd0, sb0 = homog_sources(x[:39], qT, rate_scale=1.0, rho_ash=0.0)
    np.testing.assert_allclose(sb, sb0)
    payload = stage_payload(2)
    assert payload["physical_convergence"] is False
    assert payload["homogeneous_gate"]["fe2_may_claim_physical_convergence"] is False


def test_e3_ash_inventory_and_cyclone_split() -> None:
    nq = freeze_stage_q(3)
    assert n_unknowns(3) == 41
    assert "ash_inventory" not in omitted_named(3)
    k = nq.energy.species.k
    fly = float(k.K_auf) / (float(k.K_auf) + float(k.K_ab))
    K_e = float(nq.energy.species.K_eject[1])
    m_ash = ash_feed_kg_s(nq)
    M = m_ash / max(fly * K_e, 1.0e-30)
    x = default_x0_named(nq)
    x[39] = M
    F = coupled_residual_named(x, nq)
    assert abs(float(F[39])) < 1e-12
    from dataclasses import replace

    k2 = replace(k, K_auf=2.0 * float(k.K_auf), K_ab=2.0 * float(k.K_ab))
    nq2 = freeze_stage_q(3)
    nq2 = replace(nq2, energy=replace(nq2.energy, species=replace(nq2.energy.species, k=k2)))
    F2 = coupled_residual_named(x, nq2)
    assert abs(float(F2[39])) < 1e-12
    u = unpack_named(x, 3)
    rho = rho_ash_used(nq, u)
    assert rho == pytest.approx(u["M_ash"] / nq.energy.species.geo_ex.V_d)


def test_e4_daem_x_increases_with_T_and_P_is_repeatable() -> None:
    nq = freeze_stage_q(4)
    assert "DAEM_X_from_T" not in omitted_named(4)
    x_lo = daem_X_at_T(700.0, 1.0)
    x_hi = daem_X_at_T(1100.0, 1.0)
    assert x_hi >= x_lo
    x = default_x0_named(nq)
    x_copy = x.copy()
    q1 = outer_P_named(x, nq)
    q2 = outer_P_named(x, q1)
    np.testing.assert_allclose(x, x_copy)
    np.testing.assert_allclose(q1.energy.species.n_py_d, q2.energy.species.n_py_d)
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_e5_Q_W_is_named_q_not_fitted() -> None:
    nq0 = freeze_stage_q(5, Q_W=0.0)
    nqQ = freeze_stage_q(5, Q_W=1000.0)
    assert "energy_wall_heat" not in omitted_named(5)
    x = default_x0_named(nq0)
    E0 = float(coupled_residual_named(x, nq0)[-1])
    EQ = float(coupled_residual_named(x, nqQ)[-1])
    assert EQ == pytest.approx(E0 - 1000.0, rel=0.0, abs=1e-6)
    payload = stage_payload(5, Q_W=1000.0)
    assert payload["Q_W"] == 1000.0
    assert payload["Q_W_is_fitted_to_1100K"] is False


def test_e6_tar_is_in_x() -> None:
    nq = freeze_stage_q(6)
    assert n_unknowns(6) == 43
    assert "tar_transport" not in omitted_named(6)
    x = default_x0_named(nq)
    u = unpack_named(x, 6)
    F = coupled_residual_named(x, nq)
    assert abs(float(F[40])) < 1e-12
    assert abs(float(F[41])) < 1e-12
    assert u["N_tar1"] == pytest.approx(nq.energy.n_tar["TAR1"])


def test_e7_freeboard_gas_source_changes_fb_rows() -> None:
    nq6 = freeze_stage_q(6)
    nq7 = freeze_stage_q(7)
    x = default_x0_named(nq7)
    F6 = coupled_residual_named(x, nq6)
    F7 = coupled_residual_named(x, nq7)
    assert not np.allclose(F6[16:23], F7[16:23])
    assert omitted_named(7).get("freeboard_reaction") is None


def test_e8_updates_K_eject_not_cyclone_K() -> None:
    nq = freeze_stage_q(8)
    assert "K_eject_in_this_P" not in omitted_named(8)
    assert omitted_named(8)["cyclone_K_update_P"] == "omitted, named"
    assert omitted_named(8)["production_hamel_outer_refresh"] == "omitted, named"
    x = default_x0_named(nq)
    x[-1] = 1500.0
    x_copy = x.copy()
    q2 = outer_P_named(x, nq)
    np.testing.assert_allclose(x, x_copy)
    assert q2.energy.species.k == nq.energy.species.k
    # 抛射 K=m(M=1)/M_geom 不随 T；P 仍更新轨迹气密，不更新旋风 K。
    assert q2.energy.species.geo_fb.rho_g != nq.energy.species.geo_fb.rho_g
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
