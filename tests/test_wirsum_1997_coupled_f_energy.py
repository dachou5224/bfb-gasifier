"""能量行 + 约化外层 P。不改生产默认，不把 Q_W 当闭合。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.core.species import SOLID_H_MODE_INTEGRAL
from src.solvers.wirsum_1997_coupled_f_energy import (
    N_X,
    OMITTED,
    Q_WALL_W,
    build_stage_energy_outer_payload,
    coupled_residual,
    default_x0,
    energy_residual_W,
    freeze_energy_q,
    omitted_named,
    outer_P,
    pack,
    q_at_T,
    solve_inner,
    unpack,
)
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import (
    atom_gap_from_F,
    coupled_residual as species_residual,
    solve_coupled,
)


def test_unknowns_match_equations() -> None:
    qe = freeze_energy_q()
    x = default_x0(qe)
    F = coupled_residual(x, qe)
    assert x.size == N_X == 40
    assert F.size == 40
    assert "T" in unpack(x)
    assert Q_WALL_W == 0.0
    assert omitted_named()["energy_wall_heat"] == "omitted, named"
    assert "production_hamel_outer_refresh" in OMITTED


def test_energy_uses_integral_solid_and_counts_gas_hf_once() -> None:
    qe = freeze_energy_q()
    assert qe.solid_enthalpy_mode == SOLID_H_MODE_INTEGRAL
    x = default_x0(qe)
    E = energy_residual_W(x, qe)
    assert np.isfinite(E)
    x2 = x.copy()
    x2[-1] = x[-1] + 20.0
    dE = energy_residual_W(x2, qe) - E
    assert dE != 0.0


def test_F_pure_given_q() -> None:
    qe = freeze_energy_q()
    x = default_x0(qe)
    F1 = coupled_residual(x, qe)
    F2 = coupled_residual(x, qe)
    np.testing.assert_allclose(F1, F2)
    y = x.copy()
    y[0] += 1.0e-6
    _ = coupled_residual(y, qe)
    F3 = coupled_residual(x, qe)
    np.testing.assert_allclose(F1, F3)


def test_outer_P_does_not_rewrite_x() -> None:
    qe = freeze_energy_q()
    spec = solve_coupled(qe.species)
    x = pack(spec["x"], 1200.0)
    x_copy = x.copy()
    q2 = outer_P(x, qe)
    np.testing.assert_allclose(x, x_copy)
    assert q2.species.K_eject.tolist() == qe.species.K_eject.tolist()
    assert not np.allclose(q2.species.n_py_d, qe.species.n_py_d)
    F_a = coupled_residual(x, qe)
    F_b = coupled_residual(x, q2)
    assert np.linalg.norm(F_a - F_b) > 0.0


def test_species_block_still_has_atom_identity() -> None:
    qe = freeze_energy_q()
    spec = solve_coupled(qe.species)
    gap = atom_gap_from_F(spec["x"], qe.species)
    for key in ("C", "H", "O", "N", "mass_kg_s"):
        assert gap[key] < 1.0e-9, key


def test_species_root_then_energy_row_does_not_claim_furnace() -> None:
    qe = freeze_energy_q()
    inner = solve_inner(qe)
    assert inner["species_physical"] is True
    assert inner["rewrites_unknown_N"] is False
    assert np.isfinite(inner["energy_W"])
    assert inner["T_K"] > 0.0
    # Hegermann 焦炭给料后绝热 40 维仍有根；T≈743 K 不是气化炉工况。
    assert inner["physical_convergence"] is True
    assert inner["rank_last"] == 40
    assert 700.0 < inner["T_K"] < 850.0
    assert abs(inner["T_K"] - 1100.0) > 200.0
    assert abs(inner["energy_W"]) < 1.0e-6
    qT = q_at_T(qe.species, inner["T_K"])
    Fs = species_residual(inner["x"][:39], qT)
    assert inner["species_residual_l2"] < 1.0e-10
    assert float(np.linalg.norm(Fs)) < 1.0e-8


def test_payload_and_production_defaults() -> None:
    payload = build_stage_energy_outer_payload()
    assert payload["solved"] is False
    assert payload["not_a_furnace_convergence_claim"] is True
    assert payload["n_unknowns"] == 40
    assert payload["f_name"] == "F_E0"
    assert payload["Q_W_meaning"] == "omitted_named"
    assert payload["outer_P_X_vm_frozen"] is True
    assert payload["Q_W"] == 0.0
    assert payload["gas_hf_counted_once"] is True
    assert payload["outer_rewrites_N"] is False
    assert payload["outer_rewrites_T"] is False
    assert payload["species_physical"] is True
    hist = payload["outer_history"]
    assert len(hist) == 2
    assert hist[0]["x_unchanged"] is True
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    assert cfg.solid_enthalpy_mode == "hamel_eq28_component"
    assert "y_CO_target" not in str(payload)


def test_dry_mass_matches_cho_atoms() -> None:
    from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import default_char_feed
    from src.solvers.wirsum_1997_pyrolysis_ledger import dry_mass_from_cho_atoms

    feed = default_char_feed(m_char_kg_s=0.002)
    assert feed.m_dry_kg_s == pytest.approx(dry_mass_from_cho_atoms(feed), rel=1e-12)
    assert feed.m_dry_kg_s < 0.02


def test_E_of_T_scan_buckets_and_sign_change() -> None:
    from src.solvers.wirsum_1997_coupled_f_energy import (
        enthalpy_buckets_W,
        scan_E_at_species_roots,
    )

    qe = freeze_energy_q()
    b = enthalpy_buckets_W(default_x0(qe), qe)
    assert b["sum_W"] == pytest.approx(b["energy_W"], abs=1e-6)
    assert b["Q_W"] == 0.0
    T = np.array([500.0, 700.0, 1100.0, 1500.0])
    out = scan_E_at_species_roots(qe, T, refresh_pyrolysis=False)
    assert bool(np.all(out["species_physical"]))
    assert out["sign_change"] is True
    assert out["not_a_furnace_convergence_claim"] is True
    E500 = out["rows"][0]["energy_W"]
    E1500 = out["rows"][-1]["energy_W"]
    assert E500 > 0.0
    assert E1500 < 0.0
