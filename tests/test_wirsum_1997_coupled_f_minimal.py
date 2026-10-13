"""Q6：等温单粒径联立 F。一张残差，不是四套根拼接。"""

from __future__ import annotations

import numpy as np
import pytest

from src.solvers.wirsum_1997_coupled_f_minimal import (
    N_X,
    OMITTED,
    atom_gap_from_F,
    build_stage_q6_payload,
    coupled_residual,
    default_x0,
    edges_from_x,
    freeze_q,
    omitted_named,
    pack,
    q_diagnostics,
    r1_extent_from_x,
    recycle_changes_bed_root,
    solve_coupled,
    unpack,
)
from src.solvers.wirsum_1997_cyclone_cells import CycloneCellK


def test_unknowns_match_equations() -> None:
    q = freeze_q()
    x = default_x0(q)
    F = coupled_residual(x, q)
    assert x.size == N_X == 38
    assert F.size == 38
    assert set(unpack(x)) == {"N_d", "N_b", "M", "N_fb", "N_c1", "N_c2", "M_c1", "M_c2"}
    assert q_diagnostics(q)["n_unknowns"] == 38


def test_r1_extent_is_unique() -> None:
    q = freeze_q()
    x = default_x0(q)
    a = r1_extent_from_x(x, q)
    b = r1_extent_from_x(x, q)
    assert a["e1"] == pytest.approx(b["e1"])
    F = coupled_residual(x, q)
    src_from_e1 = a["e1"]
    # 床固残差含同一 e1
    ed = edges_from_x(x, q)
    from src.core.species import MW_C_KG_PER_MOL

    F_M = float(F[14])
    expect = q.m_feed + ed["m_return"] + ed["m_rez"] - ed["m_eject"] - MW_C_KG_PER_MOL * src_from_e1
    assert F_M == pytest.approx(expect, abs=1e-12)


def test_edges_cancel_on_same_cut() -> None:
    q = freeze_q()
    x = default_x0(q)
    ed = edges_from_x(x, q)
    assert ed["m_eject"] == pytest.approx(ed["m_leave"] + ed["m_return"])
    assert ed["m_rez"] == pytest.approx(ed["m_ab0"])


def test_recycle_k_changes_bed_residual() -> None:
    q1 = freeze_q(k=CycloneCellK(K_auf=0.25, K_ab=0.75))
    q2 = freeze_q(k=CycloneCellK(K_auf=0.25, K_ab=2.25))
    x = default_x0(q1)
    F1 = coupled_residual(x, q1)
    F2 = coupled_residual(x, q2)
    assert F1[14] != pytest.approx(F2[14])
    assert F1[37] != pytest.approx(F2[37])


def test_atom_and_mass_identity_at_arbitrary_x() -> None:
    q = freeze_q()
    x = default_x0(q)
    x = x * 1.17
    x[14] = 0.4
    gap = atom_gap_from_F(x, q)
    for key, val in gap.items():
        assert val < 1.0e-9, key


def test_omitted_are_named() -> None:
    names = omitted_named()
    assert set(names) == set(OMITTED)
    assert all(v == "omitted, named" for v in names.values())
    assert "energy_wall_heat" in names
    assert "cyclone_K_update_P" in names
    d = q_diagnostics(freeze_q())
    assert d["engineering_K_fb_leave"] == "not_used"
    assert d["cyclone"]["not_a_separation_law"] is True


def test_pack_unpack_roundtrip() -> None:
    q = freeze_q()
    x = default_x0(q)
    u = unpack(x)
    y = pack(**u)
    np.testing.assert_allclose(x, y)


def test_newton_finds_physical_root() -> None:
    q = freeze_q()
    out = solve_coupled(q)
    assert out["physical_convergence"] is True
    assert out["residual_l2"] < 1.0e-10
    assert out["rank_last"] == 38
    assert out["q_purity_gap"] < 1.0e-14
    assert float(np.min(out["x"])) > 0.0
    for key, val in out["atom_gap"].items():
        assert val < 1.0e-9, key
    F_again = coupled_residual(out["x"], q)
    np.testing.assert_allclose(F_again, out["F"], atol=1e-14)


def test_recycle_k_changes_bed_root() -> None:
    rec = recycle_changes_bed_root(
        k_lo=CycloneCellK(K_auf=0.25, K_ab=0.75),
        k_hi=CycloneCellK(K_auf=0.25, K_ab=2.25),
    )
    assert rec["both_physical"] is True
    assert rec["gap_kg"] > 1.0e-6


def test_payload_does_not_claim_furnace() -> None:
    payload = build_stage_q6_payload()
    assert payload["solved"] is False
    assert payload["not_a_furnace_convergence_claim"] is True
    assert payload["physical_convergence"] is True
    assert payload["n_unknowns"] == 38
    assert "y_CO_target" not in str(payload)
