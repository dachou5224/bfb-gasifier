"""C+D 接入联立 F：两级分迁移 + 热解 q 源。不改生产默认。"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.core.species import MW_C_KG_PER_MOL
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import (
    N_X,
    OMITTED,
    atom_gap_from_F,
    build_stage_cd_coupled_payload,
    coupled_residual,
    default_x0,
    edges_from_x,
    extents_from_x,
    freeze_q,
    omitted_named,
    pack,
    q_diagnostics,
    recycle_changes_fine_root,
    solve_coupled,
    unpack,
)
from src.solvers.wirsum_1997_cyclone_cells import CycloneCellK


def test_unknowns_match_equations() -> None:
    q = freeze_q()
    x = default_x0(q)
    F = coupled_residual(x, q)
    assert x.size == N_X == 39
    assert F.size == 39
    assert set(unpack(x)) == {
        "N_d",
        "N_b",
        "M_fine",
        "M_coarse",
        "N_fb",
        "N_c1",
        "N_c2",
        "M_c1",
        "M_c2",
    }


def test_migrate_cancels_between_bins() -> None:
    q = freeze_q()
    x = default_x0(q)
    F = coupled_residual(x, q)
    ed = edges_from_x(x, q)
    ext = extents_from_x(x, q)
    mig = float(ext["m_migrate_kg_s"])
    assert mig > 0.0
    F_sum = float(F[14]) + float(F[15])
    expect = (
        float(np.sum(q.m_feed))
        + float(np.sum(ed["m_return"]))
        + float(ed["m_rez"])
        - float(np.sum(ed["m_eject"]))
        - float(np.sum(ext["consume_kg_s"]))
    )
    assert F_sum == pytest.approx(expect, abs=1e-12)


def test_pyrolysis_is_q_and_does_not_rewrite_N() -> None:
    q = freeze_q()
    x = default_x0(q)
    F1 = coupled_residual(x, q)
    F1b = coupled_residual(x, q)
    np.testing.assert_allclose(F1, F1b)
    q2 = replace(q, n_py_d=q.n_py_d * 1.25)
    F2 = coupled_residual(x, q2)
    assert F1[0:7] != pytest.approx(F2[0:7])
    u = unpack(x)
    assert u["N_d"].shape == (7,)
    d = q_diagnostics(q)
    assert d["rewrites_unknown_N"] is False
    assert d["vm_pool_basis"] == d["expected_vm_pool_basis"]
    assert d["hides_H_over_umf"] is False
    assert float(np.sum(q.n_py_d)) > 0.0
    assert q.nC_tar_omitted > 0.0
    assert q.m_feed[0] == pytest.approx(0.0)
    assert q.m_feed[1] == pytest.approx(q.nC_koks * MW_C_KG_PER_MOL)
    assert q.nC_koks != pytest.approx(q.nC_fixed)


def test_r1_consume_matches_gas_and_solid() -> None:
    q = freeze_q()
    x = default_x0(q)
    ext = extents_from_x(x, q)
    F = coupled_residual(x, q)
    ed = edges_from_x(x, q)
    F_Mc = float(F[15])
    expect = (
        float(q.m_feed[1])
        + float(ed["m_return"][1])
        + float(ed["m_rez_bin"][1])
        - float(ext["m_migrate_kg_s"])
        - float(ed["m_eject"][1])
        - float(ext["consume_kg_s"][1])
    )
    assert F_Mc == pytest.approx(expect, abs=1e-12)
    assert ext["consume_kg_s"][0] == pytest.approx(
        MW_C_KG_PER_MOL * ext["e1"][0], abs=1e-12
    )
    assert ext["model"] == "two_bin_cohort_representative"


def test_atom_and_mass_identity_at_arbitrary_x() -> None:
    q = freeze_q()
    x = default_x0(q) * 1.13
    x[14] = 0.3
    x[15] = 0.7
    gap = atom_gap_from_F(x, q)
    for key in ("C", "H", "O", "N", "mass_kg_s"):
        assert gap[key] < 1.0e-9, key
    nC_q = q.nC_koks + q.nC_tar_omitted + float(
        q.n_py_d[0] + q.n_py_d[1] + q.n_py_d[4]
    ) + float(q.nC_unreleased)
    assert nC_q == pytest.approx(q.nC_fuel, abs=1e-9)


def test_omitted_are_named() -> None:
    names = omitted_named()
    assert set(names) == set(OMITTED)
    assert "tar_transport" in names
    assert "ash_inventory" in names
    assert "energy_wall_heat" in names
    assert "pyrolysis_drying" not in names
    assert "multi_size_migration" not in names
    d = q_diagnostics(freeze_q())
    assert d["engineering_K_fb_leave"] == "not_used"
    assert d["cyclone"]["not_a_separation_law"] is True
    assert d["fate"][0] == "reaches_next_cell"
    assert d["fate"][1] == "returns"


def test_pack_unpack_roundtrip() -> None:
    q = freeze_q()
    x = default_x0(q)
    y = pack(**unpack(x))
    np.testing.assert_allclose(x, y)


def test_newton_finds_physical_root() -> None:
    q = freeze_q()
    out = solve_coupled(q)
    assert out["physical_convergence"] is True
    assert out["residual_l2"] < 1.0e-10
    assert out["rank_last"] == 39
    assert out["q_purity_gap"] < 1.0e-14
    assert float(np.min(out["x"])) > 0.0
    for key in ("C", "H", "O", "N", "mass_kg_s"):
        assert out["atom_gap"][key] < 1.0e-9, key
    F_again = coupled_residual(out["x"], q)
    np.testing.assert_allclose(F_again, out["F"], atol=1e-14)
    ext = out["extent"]
    assert ext["m_migrate_kg_s"] > 0.0


def test_recycle_k_changes_fine_root() -> None:
    rec = recycle_changes_fine_root(
        k_lo=CycloneCellK(K_auf=0.25, K_ab=0.75),
        k_hi=CycloneCellK(K_auf=0.25, K_ab=2.25),
    )
    assert rec["both_physical"] is True
    assert rec["gap_kg"] > 1.0e-6


def test_payload_and_production_defaults_unchanged() -> None:
    payload = build_stage_cd_coupled_payload()
    assert payload["solved"] is False
    assert payload["not_a_furnace_convergence_claim"] is True
    assert payload["n_unknowns"] == 39
    assert payload["physical_convergence"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    assert cfg.solid_transport_inventory_mode == "engineering_bottom_threshold"
    assert cfg.freeboard_trajectory_model == "analytical_wirsum"
    assert "y_CO_target" not in str(payload)
