"""Wesseling 1 输入合同：缺项保持缺失，循环不进外部账。"""

import pytest

from src.solvers.wesseling_case_contract import (
    air_kmol_h,
    external_element_balance_kmol_h,
    external_enthalpy_w,
    fuel_element_kmol_h,
    map_wesseling_inputs,
    wesseling1_external_inlet_kmol_h,
    aa3_exit,
    wesseling12_conditional_definitions,
    wesseling12_input_ledger,
)


def test_wesseling1_keeps_missing_inputs_and_separates_tolerances() -> None:
    mapped = map_wesseling_inputs("CASE_HTW_WESSELING_1")
    fields = {row["name"]: row for row in mapped["fields"]}
    assert fields["bed_material"]["status"] == "missing"
    assert fields["primary_air_Nm3_h"]["status"] == "missing"
    assert fields["steam_kg_h"]["status"] == "missing"
    assert mapped["conditional"] is True
    obs = mapped["observation"]
    assert obs["paper_simulation_temperature_K"] == 1100
    assert obs["measured_temperature_K"] == 1120
    assert obs["digitization_rel"] == 0.05
    assert obs["model_allowance"] == "not_set_from_results"
    fuel = fuel_element_kmol_h("CASE_HTW_WESSELING_1")
    assert fuel["dry_fuel_kg_h"] == pytest.approx(2806.6, abs=0.2)
    assert fuel["C"] > 0.0
    air = air_kmol_h(5536.9, 1162.7)
    assert air["O2_kmol_h"] / (air["O2_kmol_h"] + air["N2_kmol_h"]) == pytest.approx(0.210, abs=0.002)


def test_wesseling12_ledger_keeps_heating_value_and_size_mass_empty() -> None:
    led = wesseling12_input_ledger()
    heat = led["heating_value"]
    size = led["fuel_size"]
    assert heat["value"] is None
    assert heat["hhv_or_lhv"] is None
    assert heat["dry_or_as_received"] is None
    assert heat["fitted_from_exit"] is False
    assert size["wesseling_feed_class_mass_kg_h"] is None
    assert size["wesseling_feed_psd_in_table_7_1"] is False
    assert size["applied_to_char_bins"] is False
    assert size["model_discretization_n_classes"] == 10
    air = led["cases"]["CASE_HTW_WESSELING_1"]
    oxygen = led["cases"]["CASE_HTW_WESSELING_2"]
    assert air["agent"] == "Air/Steam (Luft/Dampf)"
    assert oxygen["agent"] == "O2/Steam (O2/Dampf)"
    closure = led["agent_closure"]
    assert closure["steam_kg_h_in_table_7_1_ocr"] is None
    assert closure["secondary_Nm3_h_in_table_7_1_ocr"] is None
    assert closure["secondary_Nm3_h_in_case_json"]["CASE_HTW_WESSELING_1"] == pytest.approx(30.1)
    assert closure["secondary_Nm3_h_used_as_measurement"] is False
    assert "128-129" in closure["pdf_visual_check"]
    assert air["secondary_Nm3_h_measured"] is None
    assert closure["complete_case"] is False
    assert closure["experimental_comparison_blocked"] is True
    assert closure["public_background_sources_supply_the_missing_numbers"] is False
    assert closure["numerical_target_steam_is_an_assumption"] is True
    assert closure["two_zone_can_place_freeboard_injection"] is False
    assert oxygen["steam_status"] == "missing"
    assert heat["present_in_table_7_1_ocr"] is False
    assert heat["thermal_input_W"] is None
    assert led["experimental_comparison_done"] is False
    assert led["bypass_accounting"]["enters_seven_species_molar_total"] is False


def test_wesseling1_definition_keeps_missing_inputs_conditional() -> None:
    definition = wesseling12_conditional_definitions()
    air = definition["table_fields"]["CASE_HTW_WESSELING_1"]
    oxygen = definition["table_fields"]["CASE_HTW_WESSELING_2"]
    assert air["agent"] == "Air/Steam (Luft/Dampf)"
    assert air["lambda_O2"] == pytest.approx(0.337)
    assert oxygen["agent"] == "O2/Steam (O2/Dampf)"
    assert oxygen["lambda_O2"] == pytest.approx(0.268)
    assert "steam_kg_h" in definition["unavailable_not_filled"]
    assert definition["numerical_assumptions"]["steam_kg_h"] == 0.0
    assert definition["sources_may_still_exist"] is True
    assert definition["experimental_comparison_done"] is False
    assert definition["two_zone_certification_is_not_the_multicell_map"] is True


def test_aa3_keeps_each_experimental_input_open() -> None:
    exit_ = aa3_exit()
    assert exit_["milestone_reached"]["conditional_s1"] is False
    assert exit_["milestone_reached"]["paired_experiment"] is False
    assert exit_["named_kernel_s1_is_only_the_first_milestone"] is True
    assert exit_["fitted_from_exit"] is False
    for item in exit_["inputs"].values():
        assert item["assumption_is_a_measurement"] is False
    assert exit_["inputs"]["steam"]["measured"] is None
    assert exit_["inputs"]["steam"]["numerical_assumption"]["steam_kg_h"] == 0.0
    assert exit_["inputs"]["heat_loss"]["measured_W"] is None
    assert exit_["inputs"]["heat_loss"]["numerical_assumption"]["Q_W"] == 0.0
    assert exit_["inputs"]["feed_size"]["measured_sieve_mass_kg_h"] is None
    assert exit_["inputs"]["staged_injection"]["measured_secondary_Nm3_h"] is None
    assert exit_["inputs"]["staged_injection"]["two_zone_can_place_it"] is False
    assert exit_["experimental_comparison_done"] is False


def test_wesseling2_steam_stays_missing() -> None:
    mapped = map_wesseling_inputs("CASE_HTW_WESSELING_2")
    fields = {row["name"]: row for row in mapped["fields"]}
    assert fields["steam_kg_h"]["status"] == "missing"
    assert fields["primary_O2_Nm3_h"]["status"] == "missing"
    assert mapped["conditional"] is True


def test_recycle_cancels_and_is_not_an_external_flow() -> None:
    balance = external_element_balance_kmol_h(
        {"C": 10.0},
        {"C": 4.0},
        recycle={"up": 3.0, "down": -3.0},
    )
    assert balance["C"] == pytest.approx(6.0)
    assert external_enthalpy_w(1000.0, 400.0, recycle_up_w=50.0, recycle_down_w=-50.0) == pytest.approx(600.0)


def test_wesseling1_inlet_is_larger_than_the_named_kernel_and_omits_steam() -> None:
    from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

    inlet = wesseling1_external_inlet_kmol_h()
    assert inlet["steam_included"] is False
    assert inlet["primary_split_included"] is False
    named = freeze_e8_q(Q_W=0.0)
    o2_named = float(named.energy.species.N_in_b[5])
    o2_case = float(inlet["O2_kmol_h"]) * 1000.0 / 3600.0
    assert o2_case > 10.0 * o2_named
    closed = external_element_balance_kmol_h(
        {"C": float(inlet["C"])},
        {"C": 0.0},
        recycle={"cyclone_to_bed": 1.0, "bed_to_cyclone": -1.0},
    )
    assert closed["C"] == pytest.approx(float(inlet["C"]))


def test_recycle_must_cancel_each_element_not_only_the_sum() -> None:
    from src.solvers.wesseling_case_contract import element_pair_cancels

    assert element_pair_cancels({"C": 2.0, "H": 1.0}, {"C": -2.0, "H": -1.0})
    assert element_pair_cancels({"C": 2.0, "H": 0.0}, {"C": 0.0, "H": -2.0}) is False


def test_wesseling1_ash_bed_and_heat_loss_stay_on_separate_bases() -> None:
    from src.solvers.wesseling_case_contract import wesseling1_inventory_and_energy

    book = wesseling1_inventory_and_energy()
    assert book["reactive_ash_from_fuel_kg_h"] == pytest.approx(2806.6 * 0.1141, rel=0.002)
    assert book["fixed_bed_inventory_kg"] is None
    assert book["heat_loss_W"] is None
    assert book["recycle_mass_flow_kg_h"] is None
    assert book["exit_mole_fractions_are_not_an_element_balance"] is True
    assert book["heat_loss_fraction_range"] == (0.05, 0.10)


def test_wesseling1_conditional_boundary_uses_the_pilot_geometry() -> None:
    from src.solvers.wesseling_case_contract import wesseling1_conditional_boundary
    from src.solvers.wesseling_two_zone import build_wesseling1_two_zone_q

    bound = wesseling1_conditional_boundary()
    assert bound["steam_status"] == "missing"
    assert bound["secondary_air_placed_at_xi"] is False
    assert bound["two_zone_scores_secondary_temperature_peak"] is False
    assert bound["model_allowance"] == "not_set_from_results"
    assert bound["O2_mol_s"] == pytest.approx(1162.7 / 22.414 / 3.6, rel=0.002)
    assert bound["area_m2"] == pytest.approx(3.141592653589793 * 0.09, rel=1.0e-6)
    assert bound["H_freeboard_m"] == pytest.approx(8.5)
    nq = build_wesseling1_two_zone_q()
    assert nq.stage == 11
    assert float(nq.reference.area) == pytest.approx(bound["area_m2"])
    assert float(nq.reference.H_bed) == pytest.approx(6.0)
    assert float(nq.energy.species.N_in_b[5]) == pytest.approx(bound["O2_mol_s"])
    assert float(nq.energy.T_in) == pytest.approx(293.0)
    assert float(nq.Q_W) == 0.0
    assert abs(float(nq.energy.feed.h_f_dry_J_kg) - (-3.96e6)) < 1.0


def test_path_step_is_between_the_kernels_and_is_not_named_wesseling() -> None:
    from src.solvers.wesseling_two_zone import q_on_wesseling1_path, wesseling1_input_path

    named = q_on_wesseling1_path(0.0)
    mid = q_on_wesseling1_path(0.02)
    target = q_on_wesseling1_path(1.0)
    o2_named = float(named.energy.species.N_in_b[5])
    o2_mid = float(mid.energy.species.N_in_b[5])
    o2_target = float(target.energy.species.N_in_b[5])
    assert o2_named < o2_mid < o2_target
    from src.core.elemental_ledger import H2O_MOLAR_MASS_KG_PER_MOL

    assert float(mid.reference.trajectory_continue_s) == pytest.approx(0.02 * 15.0)
    assert float(mid.reference.trajectory_extra_physical_s) == pytest.approx(20.0)
    assert int(mid.reference.trajectory_max_work_steps) == 20000
    assert float(mid.reference.trajectory_wall_s) == pytest.approx(5.0)
    assert float(named.reference.trajectory_extra_physical_s) == 0.0
    assert float(mid.energy.species.n_moist[3]) == pytest.approx(
        float(mid.energy.feed.m_moist_kg_s) / H2O_MOLAR_MASS_KG_PER_MOL
    )
    assert float(mid.energy.feed.h_f_dry_J_kg) == pytest.approx(-3.96e6)
    assert tuple(float(v) for v in mid.energy.species.d_p) == tuple(float(v) for v in named.energy.species.d_p)
    path = wesseling1_input_path()
    assert path["intermediate_name"] == "not CASE_HTW_WESSELING_1"
    assert path["target_reached"] is False
    assert path["fuel_size_applied_to_bins"] is False
    assert path["fuel_heating_value"] is None


def test_path_limit_matches_the_target_inlets() -> None:
    import numpy as np

    from src.core.elemental_ledger import H2O_MOLAR_MASS_KG_PER_MOL
    from src.solvers.wesseling_two_zone import (
        build_wesseling1_two_zone_q,
        q_on_wesseling1_path,
        wesseling1_input_path,
    )

    def water(q):
        return float(q.energy.feed.m_moist_kg_s) / H2O_MOLAR_MASS_KG_PER_MOL

    for s in (0.0, 0.02, 0.22, 0.999999, 1.0):
        q = q_on_wesseling1_path(s)
        assert float(q.energy.species.n_moist[3]) == pytest.approx(water(q))
    end = q_on_wesseling1_path(1.0)
    near = q_on_wesseling1_path(0.999999)
    assert np.allclose(end.energy.species.N_in_d, 0.0)
    assert np.max(np.abs(near.energy.species.N_in_d - end.energy.species.N_in_d)) < 1.0e-4
    assert np.allclose(near.energy.species.N_in_b, end.energy.species.N_in_b, rtol=1.0e-4, atol=1.0e-6)
    assert float(end.reference.trajectory_continue_s) == pytest.approx(15.0)
    assert abs(float(near.reference.trajectory_continue_s) - 15.0) < 1.0e-3
    built = build_wesseling1_two_zone_q()
    assert np.allclose(end.energy.species.n_py_d, built.energy.species.n_py_d, rtol=1.0e-8, atol=1.0e-10)
    assert np.allclose(near.energy.species.n_py_d, end.energy.species.n_py_d, rtol=1.0e-4, atol=1.0e-6)
    assert end.energy.species.nC_unreleased == pytest.approx(built.energy.species.nC_unreleased)
    assert end.energy.species.nH_unreleased == pytest.approx(built.energy.species.nH_unreleased)
    assert end.energy.species.nO_unreleased == pytest.approx(built.energy.species.nO_unreleased)
    for s in (0.02, 0.22, 1.0):
        q = q_on_wesseling1_path(s)
        feed = q.energy.feed
        ash_frac = float(feed.ash_dry_wt) / 100.0 if float(feed.ash_dry_wt) > 1.0 else float(feed.ash_dry_wt)
        from src.solvers.wirsum_1997_named_f import ash_feed_kg_s
        assert ash_feed_kg_s(q) == pytest.approx(float(feed.m_dry_kg_s) * ash_frac)
        assert float(feed.m_dry_kg_s) * float(feed.h_f_dry_J_kg) == pytest.approx(
            float(feed.m_dry_kg_s) * (-3.96e6)
        )
        assert np.linalg.norm(q.energy.species.n_py_d - q_on_wesseling1_path(0.0).energy.species.n_py_d) > 0.0
    path = wesseling1_input_path()
    assert "not a physical certification frontier" in path["historical_s_0_22"]


def test_parameter_continuation_stops_short_of_an_uncertified_target() -> None:
    from src.solvers.wirsum_1997_named_f_control import advance_parameter_factor

    current, step, stop = advance_parameter_factor(1.0, 0.95, 0.05, passed=False)
    assert current == pytest.approx(1.0)
    assert step == pytest.approx(0.025)
    assert stop is False
    current, step, stop = advance_parameter_factor(current, 0.95, 0.012, passed=False)
    assert current == pytest.approx(1.0)
    assert stop is True
    reached, _, stop = advance_parameter_factor(1.0, 0.95, 0.05, passed=True)
    assert reached == pytest.approx(0.95)
    assert stop is False
