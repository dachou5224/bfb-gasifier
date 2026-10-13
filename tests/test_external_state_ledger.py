"""工作包 C：未解状态也能出外部账，适配器用量和合同对照。"""

import numpy as np
import pytest
from dataclasses import replace

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wesseling_two_zone import wesseling1_input_reconciliation


def _p3_q_and_x():
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    return outer_P_named(x, nq), x


def test_certified_endpoint_ledger_is_built_from_streams() -> None:
    q, x = _p3_q_and_x()
    book = external_imbalance(x, q)
    inlet = book["inlet_mol_s"]["total"]
    outlet = book["outlet_mol_s"]["total"]
    for el in ("C", "H", "O"):
        assert book["gap_mol_s"][el] == pytest.approx(inlet[el] - outlet[el])
    assert book["full_element_certification"] is False
    assert "NH3" in book["unmodeled"]["N"]
    assert "H2S" in book["unmodeled"]["S"]
    assert book["interface_char"]["excluded_from_external_gap"] is True
    assert book["interface_char"]["up_must_equal_down"] is False
    assert book["interface_char"]["up_stream_paired"] is True
    assert book["interface_char"]["bed_exit_up_kg_s"] == pytest.approx(
        book["interface_char"]["freeboard_inlet_up_kg_s"]
    )
    assert book["enthalpy_W"]["Q_W"] == 0.0
    assert book["enthalpy_W"]["fuel_formation"] != book["enthalpy_W"]["Q_W"]
    trial = x.copy()
    trial[-1] += 40.0
    other = external_imbalance(trial, q)
    assert other["enthalpy_W"]["gap"] != pytest.approx(book["enthalpy_W"]["gap"])


def test_upward_and_downward_char_are_different_streams() -> None:
    q, x = _p3_q_and_x()
    sp = replace(q.energy.species, w_leave=np.ones(2), w_return=np.zeros(2))
    q2 = replace(q, energy=replace(q.energy, species=sp))
    face = external_imbalance(x, q2)["interface_char"]
    assert face["up_must_equal_down"] is False
    assert face["bed_exit_up_kg_s"] == pytest.approx(face["freeboard_inlet_up_kg_s"])
    assert face["freeboard_return_after_reaction_kg_s"] == pytest.approx(face["bed_inlet_down_kg_s"])
    assert face["leave_after_reaction_kg_s"] > 0.0


def test_wesseling1_adapter_matches_applied_inputs_and_keeps_size_separate() -> None:
    report = wesseling1_input_reconciliation()
    assert report["all_applied_inputs_match"] is True
    assert report["secondary_air_Nm3_h"] is None
    assert report["secondary_air_Nm3_h_historical"] == pytest.approx(30.1)
    assert report["secondary_air_xi"] is None
    assert report["secondary_air_xi_digitized"] == pytest.approx(0.6)
    assert report["secondary_air_placed"] is False
    assert report["steam_status"] == "missing"
    assert report["fuel_size_applied_to_bins"] is False
    assert report["fuel_size_mm"]["n_size_classes"] == 10
    assert report["fuel_size_mm"]["feed_psd_measured"] is False
    assert report["adapter_char_diameter_m"][0] < 1.0e-3
    assert report["Q_W"] == 0.0
    assert report["formation_enthalpy_J_kg"] == pytest.approx(-3.96e6)


def test_nh3_h2s_bypass_is_in_the_atom_stream_and_the_energy_row() -> None:
    from src.solvers.wirsum_1997_named_f import energy_zones_named, minor_bypass_enthalpy_W
    from src.solvers.wirsum_1997_named_f import default_x0_named
    from src.solvers.wesseling_two_zone import q_on_wesseling1_path

    q = q_on_wesseling1_path(0.22)
    sp = q.energy.species
    n_nh3 = float(sp.n_NH3_bypass_mol_s)
    n_h2s = float(sp.n_H2S_bypass_mol_s)
    assert n_nh3 > 0.0 and n_h2s > 0.0
    assert sp.minor_bypass_in_major_density is False
    x = default_x0_named(q)
    x[-2:] = 1100.0
    book = external_imbalance(x, q)
    assert book["outlet_mol_s"]["minor_bypass"]["H"] == pytest.approx(3.0 * n_nh3 + 2.0 * n_h2s)
    assert book["minor_bypass"]["covers_subsequent_reactions"] is False
    assert book["minor_bypass"]["full_NS_model"] is False
    assert book["h_o_gap_not_waived"] is True
    bare = replace(sp, n_NH3_bypass_mol_s=0.0, n_H2S_bypass_mol_s=0.0)
    q_bare = replace(q, energy=replace(q.energy, species=bare))
    gap_h = external_imbalance(x, q)["gap_mol_s"]["H"]
    gap_h_bare = external_imbalance(x, q_bare)["gap_mol_s"]["H"]
    assert gap_h_bare - gap_h == pytest.approx(3.0 * n_nh3 + 2.0 * n_h2s)
    e_with = sum(energy_zones_named(x, q, 1100.0, 1100.0))
    e_bare = sum(energy_zones_named(x, q_bare, 1100.0, 1100.0))
    assert e_bare - e_with == pytest.approx(minor_bypass_enthalpy_W(sp, 1100.0))
