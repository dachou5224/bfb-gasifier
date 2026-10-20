# -*- coding: utf-8 -*-
"""T9：工况装配、验收目标和外部元素平衡。"""
from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from hamel.cases.assemble import (
    PHYSICAL_FIELDS,
    CaseAssemblyError,
    CitedFieldNotApproved,
    MissingFieldError,
    RecordOnlyFieldError,
    assemble_case,
    load_case_file,
    read_case_field,
    read_field,
)
from hamel.cases.balance import (
    BalanceInputError,
    ElementStream,
    GasStream,
    SolidStream,
    external_balance,
)
from hamel.cases.config import PhysicalConfig, SolverConfig
from hamel.cases.judge import (
    ABSENT,
    BORDERLINE,
    FAIL,
    PASS,
    PENDING,
    JudgementError,
    Measurement,
    Premises,
    judge,
)
from hamel.cases.targets import TargetError, acceptance_targets

ROOT = Path(__file__).resolve().parents[1]
VTT = ROOT / "data" / "cases" / "vtt_pressurized.json"
WSV = ROOT / "data" / "cases" / "wsv400.json"
DIGITIZED = ROOT / "data" / "digitized"
V_M = 0.022414

# 摘自 origin/task/T1-probe-wsv400a:probe/out/thesis_default.json 的 balance 段。
# 方法字符串：inlet_atoms - (bed_top gas atoms + bed_top solid upflow atoms)。
# 未复制整个文件。碳缺口 gap_mol_atom_s.C = -0.5602918220004856。
_T1_INLET = {"C": 0.2738387515148151, "H": 0.715256443614865, "O": 0.5909509026926438}
_T1_BED_TOP_GAS = {"C": 0.833952794302232, "H": 0.1939489547969353, "O": 1.0514418843153517}
_T1_BED_TOP_SOLID = {"C": 0.0001777792130686537, "H": 0.0, "O": 0.0}
_T1_CARBON_GAP = -0.5602918220004856

_OK = Premises(residuals_met=True, balance_closed=True, state_finite_nonnegative=True)

# docs/STATUS.md「WSV400 A–F 验收目标」表中的 Hamel 模拟值。
_WSV_SIM = {
    "A": {"T_FRB_K": 1048, "CO": 0.169, "CO2": 0.122, "H2": 0.081, "CH4": 0.022},
    "B": {"T_FRB_K": 997, "CO": 0.198, "CO2": 0.115, "H2": 0.090, "CH4": 0.033},
    "C": {"T_FRB_K": 1027, "CO": 0.183, "CO2": 0.124, "H2": 0.097, "CH4": 0.026},
    "D": {"T_FRB_K": 1067, "CO": 0.164, "CO2": 0.129, "H2": 0.090, "CH4": 0.021},
    "E": {"T_FRB_K": 1094, "CO": 0.149, "CO2": 0.134, "H2": 0.078, "CH4": 0.017},
    "F": {"T_FRB_K": 1044, "CO": 0.182, "CO2": 0.116, "H2": 0.096, "CH4": 0.027},
}


def _mol_s(volume_m3_per_h: float) -> float:
    return volume_m3_per_h / V_M / 3600.0


def _by_name(items):
    return {item.quantity: item for item in items}


def test_physical_config_has_no_defaults_and_solver_fields_match_their_adrs():
    assert len(PHYSICAL_FIELDS) <= 30
    assert {field.name for field in dataclasses.fields(PhysicalConfig)} == set(PHYSICAL_FIELDS)
    solver = SolverConfig()
    assert {field.name for field in dataclasses.fields(SolverConfig)} == {
        "pyrolysis_radial_intervals",
        "pyrolysis_time_rtol",
        "pyrolysis_activation_nodes",
        "pyrolysis_tail_delta_t_k",
        "pyrolysis_initial_shell_fraction",
        "freeboard_ejection_samples",
        "freeboard_max_step_s",
        "freeboard_time_limit_s",
        "freeboard_delta_h_m",
        "freeboard_delta_u_m_s",
    }
    assert solver.pyrolysis_radial_intervals == 24
    assert solver.pyrolysis_time_rtol == 1.0e-6
    assert solver.pyrolysis_activation_nodes == 16
    assert solver.pyrolysis_tail_delta_t_k == 0.02
    assert solver.pyrolysis_initial_shell_fraction == 1.0e-3
    assert solver.freeboard_ejection_samples == 48
    assert solver.freeboard_max_step_s == 2.0e-4
    assert solver.freeboard_time_limit_s == 80.0
    assert solver.freeboard_delta_h_m == 1.0e-3
    assert solver.freeboard_delta_u_m_s == 1.0e-2
    with pytest.raises(TypeError):
        PhysicalConfig()
    with pytest.raises(TypeError):
        Premises()


def test_assemble_vtt_pressurized_1():
    assembled = assemble_case(VTT, "CASE_VTT_PRESSURIZED_1")
    assert assembled.plant_id == "VTT_PRESSURIZED"
    assert assembled.case_key == "1"
    assert assembled.case_id == "CASE_VTT_PRESSURIZED_1"
    physical = assembled.physical
    assert physical.bed_diameter == 0.15
    assert physical.freeboard_diameter == 0.25
    assert physical.reactor_height == 3.5
    assert physical.bed_zone_height == 1.2
    assert physical.freeboard_height == 2.3
    assert physical.fuel_feed_height == 0.12
    assert physical.cyclone_cell_count == 2
    assert physical.fuel_feed_count == 1
    assert physical.agent_feed_count == 2
    assert physical.carbon_mass_fraction_wf == 0.638
    assert physical.fixed_carbon_mass_fraction_wf == 0.427
    assert physical.particle_size_max == 0.004
    assert physical.lambda_o2 == 0.33
    assert physical.pressure == 500_000.0
    assert physical.inlet_temperature == 458.15
    assert physical.heating_value == 21_000_000.0
    assert physical.recirculation is False
    assert physical.moisture_mass_fraction == pytest.approx(0.118)
    assert physical.fuel_mass_flow == pytest.approx(46.9 / 3600.0)
    assert physical.steam_mass_flow == pytest.approx(6.8 / 3600.0)
    assert physical.primary_air_amount_flow == pytest.approx(_mol_s(51.8))
    assert physical.secondary_air_amount_flow == pytest.approx(_mol_s(30.1))
    assert physical.nitrogen_purge_volume_flow == pytest.approx(_mol_s(6.2))
    assert read_case_field(VTT, "1", "bed_diameter") == 0.15


def test_missing_field_is_an_error_on_both_plants():
    with pytest.raises(MissingFieldError, match="bed_inventory_diameter"):
        read_case_field(VTT, "1", "bed_inventory_diameter")
    with pytest.raises(MissingFieldError, match="inlet_gas_mole_fractions"):
        read_case_field(VTT, "2", "inlet_gas_mole_fractions")
    with pytest.raises(MissingFieldError, match="inlet_temperature"):
        assemble_case(VTT, "2")
    with pytest.raises(MissingFieldError, match="bed_height"):
        read_case_field(WSV, "A", "bed_height")


def test_cited_field_needs_an_adr_and_record_only_stays_out():
    with pytest.raises(CitedFieldNotApproved, match="nitrogen_purge_volume_flow"):
        read_case_field(VTT, "3", "nitrogen_purge_volume_flow")
    with pytest.raises(RecordOnlyFieldError, match="thermal_design_power"):
        read_case_field(VTT, "1", "thermal_design_power")
    assert read_case_field(WSV, "A", "fuel_feed_height") == pytest.approx(0.2)


def test_vtt1_targets_follow_the_charter_formula():
    targets = _by_name(acceptance_targets("CASE_VTT_PRESSURIZED_1", DIGITIZED))
    assert list(targets) == ["T_FRB_K", "X_C", "CO", "CO2", "H2", "CH4", "H2O"]
    expected = {
        "T_FRB_K": (1164.0, 20.0, 13.0, None),
        "X_C": (78.1, 3.0, 1.3, None),
        "CO": (0.164, 0.10 * 0.164, 0.009, "wet"),
        "CO2": (0.091, 0.015, 0.004, "wet"),
        "H2": (0.109, 0.015, 0.004, "wet"),
        "CH4": (0.013, 0.005, 0.004, "wet"),
        "H2O": (0.122, 0.015, 0.013, "wet"),
    }
    for name, (sim, width, error, basis) in expected.items():
        target = targets[name]
        assert target.present
        assert target.sim == sim
        assert target.lower == pytest.approx(sim - width)
        assert target.upper == pytest.approx(sim + width)
        assert target.est_error == error
        assert target.basis == basis
        assert target.source == "fig_7_10.json"
    # STATUS 表里的 CO 区间是舍入显示；代码必须用精确区间。
    assert targets["CO"].lower == pytest.approx(0.1476)
    assert targets["CO"].upper == pytest.approx(0.1804)
    assert targets["CO"].lower != pytest.approx(0.148)


def test_wsv400_targets_match_status_table():
    for letter, sims in _WSV_SIM.items():
        case_id = f"CASE_WSV400_{letter}"
        targets = _by_name(acceptance_targets(case_id, DIGITIZED))
        for name, sim in sims.items():
            assert targets[name].sim == sim
            assert targets[name].present
            assert targets[name].basis == ("dry" if name != "T_FRB_K" else None)
        assert targets["H2O"].present is False
        assert targets["X_C"].present is False
    case_a = _by_name(acceptance_targets("CASE_WSV400_A", DIGITIZED))
    assert case_a["T_FRB_K"].lower == 1028
    assert case_a["T_FRB_K"].upper == 1068
    assert case_a["T_FRB_K"].est_error == 8
    assert case_a["CO2"].lower == pytest.approx(0.107)
    assert case_a["CO2"].upper == pytest.approx(0.137)
    assert case_a["H2"].lower == pytest.approx(0.066)
    assert case_a["H2"].upper == pytest.approx(0.096)
    assert case_a["CH4"].lower == pytest.approx(0.017)
    assert case_a["CH4"].upper == pytest.approx(0.027)
    assert case_a["CO"].lower == pytest.approx(0.169 - 0.0169)
    assert case_a["CO"].upper == pytest.approx(0.169 + 0.0169)
    assert case_a["T_exit_K"].sim == 1000
    assert case_a["T_exit_K"].lower == 980
    assert case_a["T_exit_K"].upper == 1020
    assert case_a["T_exit_K"].est_error == 8
    assert case_a["T_exit_K"].source == "fig_7_16.json"
    assert "T_exit_K" not in _by_name(acceptance_targets("CASE_WSV400_B", DIGITIZED))
    with pytest.raises(TargetError, match="CASE_VTT_PRESSURIZED_2"):
        acceptance_targets("CASE_VTT_PRESSURIZED_2", DIGITIZED)


def test_balance_is_zero_for_a_closed_pair_and_reports_added_carbon():
    gas = GasStream({"CO": 1.0, "H2": 0.5, "H2O": 0.25}, enthalpy_flow_w=100.0)
    solid = SolidStream(
        mass_flow_kg_s=0.012011,
        element_mass_fraction={"C": 1.0, "H": 0.0, "O": 0.0},
        enthalpy_flow_w=5.0,
    )
    closed = external_balance([gas, solid], [gas, solid])
    assert closed.C_mol_s == pytest.approx(0.0)
    assert closed.H_mol_s == pytest.approx(0.0)
    assert closed.O_mol_s == pytest.approx(0.0)
    assert closed.energy_w == pytest.approx(0.0)

    extra = ElementStream({"C": 0.25, "H": 0.0, "O": 0.0}, enthalpy_flow_w=0.0)
    opened = external_balance([gas], [gas, extra])
    assert opened.C_mol_s == pytest.approx(-0.25)
    assert opened.H_mol_s == pytest.approx(0.0)
    assert opened.O_mol_s == pytest.approx(0.0)

    hotter = GasStream({"CO": 1.0}, enthalpy_flow_w=400.0)
    colder = GasStream({"CO": 1.0}, enthalpy_flow_w=150.0)
    assert external_balance([hotter], [colder]).energy_w == pytest.approx(250.0)
    with pytest.raises(BalanceInputError, match="未知组分"):
        external_balance([GasStream({"C2H4": 1.0}, enthalpy_flow_w=0.0)], [])


def test_t1_bed_top_carbon_gap():
    # 能量数不在这三组原子里。两侧显式写 0，只核对碳；不把它当成 T1 的能量残差。
    gap = external_balance(
        [ElementStream(_T1_INLET, enthalpy_flow_w=0.0)],
        [
            ElementStream(_T1_BED_TOP_GAS, enthalpy_flow_w=0.0),
            ElementStream(_T1_BED_TOP_SOLID, enthalpy_flow_w=0.0),
        ],
    )
    assert gap.C_mol_s == pytest.approx(_T1_CARBON_GAP)
    assert gap.C_mol_s == pytest.approx(-0.56, abs=0.01)


def _measured(targets, **overrides):
    measured = {
        item.quantity: Measurement(item.sim, item.unit, item.basis)
        for item in targets
        if item.present
    }
    for key, value in overrides.items():
        current = measured[key]
        measured[key] = (
            value if isinstance(value, Measurement) else Measurement(value, current.unit, current.basis)
        )
    return measured


def _vtt_measured(**overrides):
    targets = acceptance_targets("CASE_VTT_PRESSURIZED_1", DIGITIZED)
    return targets, _measured(targets, **overrides)


def test_judgement_covers_pass_fail_and_both_sides_of_the_boundary():
    targets, measured = _vtt_measured()
    passed = judge("CASE_VTT_PRESSURIZED_1", measured, targets, _OK)
    assert passed.status == PASS
    # CO2 合格区间 [0.076, 0.106]，读数误差 0.004。中心点远离两端。
    assert _by_name(passed.quantities)["CO2"].status == PASS

    outside = _measured(targets, CO2=0.076 - 0.02)
    failed = judge("CASE_VTT_PRESSURIZED_1", outside, targets, _OK)
    assert _by_name(failed.quantities)["CO2"].status == FAIL
    assert failed.status == FAIL

    inside_near = _measured(targets, CO2=0.076 + 0.001)
    near_in = judge("CASE_VTT_PRESSURIZED_1", inside_near, targets, _OK)
    assert _by_name(near_in.quantities)["CO2"].status == BORDERLINE
    assert near_in.status == PENDING

    outside_near = _measured(targets, CO2=0.076 - 0.001)
    near_out = judge("CASE_VTT_PRESSURIZED_1", outside_near, targets, _OK)
    assert _by_name(near_out.quantities)["CO2"].status == BORDERLINE
    assert near_out.status == PENDING

    on_error = _measured(targets, CO2=0.076 + 0.004)
    assert _by_name(judge("CASE_VTT_PRESSURIZED_1", on_error, targets, _OK).quantities)["CO2"].status == PASS
    beyond = _measured(targets, CO2=0.076 - 0.004)
    assert _by_name(judge("CASE_VTT_PRESSURIZED_1", beyond, targets, _OK).quantities)["CO2"].status == FAIL


def test_failed_premise_skips_the_boundary():
    targets, measured = _vtt_measured(CO2=0.076 + 0.001)
    verdict = judge(
        "CASE_VTT_PRESSURIZED_1",
        measured,
        targets,
        Premises(residuals_met=False, balance_closed=True, state_finite_nonnegative=True),
    )
    carbon = _by_name(verdict.quantities)["CO2"]
    assert verdict.status == FAIL
    assert verdict.premises_met is False
    assert carbon.status == FAIL
    assert carbon.status != BORDERLINE
    assert carbon.distance_to_boundary is None

    unbalanced = judge(
        "CASE_VTT_PRESSURIZED_1",
        _vtt_measured()[1],
        targets,
        Premises(residuals_met=True, balance_closed=False, state_finite_nonnegative=True),
    )
    assert unbalanced.status == FAIL
    assert unbalanced.balance_closed is False
    assert _by_name(unbalanced.quantities)["T_FRB_K"].status == FAIL


def test_absent_quantity_and_a_single_failure():
    targets = acceptance_targets("CASE_WSV400_A", DIGITIZED)
    measured = _measured(targets)
    verdict = judge("CASE_WSV400_A", measured, targets, _OK)
    assert _by_name(verdict.quantities)["H2O"].status == ABSENT
    assert _by_name(verdict.quantities)["X_C"].status == ABSENT
    assert verdict.status == PASS

    measured = _measured(targets, H2=0.2)
    failed = judge("CASE_WSV400_A", measured, targets, _OK)
    assert failed.status == FAIL
    assert _by_name(failed.quantities)["H2"].status == FAIL
    assert _by_name(failed.quantities)["CO"].status == PASS


def test_incomplete_or_foreign_targets_are_not_a_pass():
    targets, measured = _vtt_measured()
    with pytest.raises(JudgementError):
        judge("CASE_VTT_PRESSURIZED_1", _measured(targets[:1]), targets[:1], _OK)
    duplicated = targets + (targets[0],)
    with pytest.raises(JudgementError, match="重复"):
        judge("CASE_VTT_PRESSURIZED_1", measured, duplicated, _OK)
    with pytest.raises(JudgementError, match="CASE_VTT_PRESSURIZED_2"):
        judge("CASE_VTT_PRESSURIZED_2", measured, targets, _OK)


def test_rejects_mixed_basis_and_units():
    targets, measured = _vtt_measured()
    dry = dict(measured)
    dry["CO"] = Measurement(measured["CO"].value, "mol/mol", "dry")
    with pytest.raises(JudgementError, match="CO"):
        judge("CASE_VTT_PRESSURIZED_1", dry, targets, _OK)
    as_mole_fraction = dict(measured)
    as_mole_fraction["X_C"] = Measurement(78.1, "mol/mol", None)
    with pytest.raises(JudgementError, match="X_C"):
        judge("CASE_VTT_PRESSURIZED_1", as_mole_fraction, targets, _OK)
    bare = {key: item.value for key, item in measured.items()}
    with pytest.raises(JudgementError, match="T_FRB_K"):
        judge("CASE_VTT_PRESSURIZED_1", bare, targets, _OK)


def test_cited_approval_stays_on_vtt_and_dimensions_are_checked():
    doc = load_case_file(VTT)
    doc["plant_id"] = "OTHER"
    with pytest.raises(CitedFieldNotApproved, match="inlet_temperature"):
        read_field(doc, "1", "inlet_temperature")
    with pytest.raises(CitedFieldNotApproved, match="nitrogen_purge_volume_flow"):
        read_field(doc, "1", "nitrogen_purge_volume_flow")
    peat = load_case_file(VTT)
    peat["fuels"]["finn_torf"]["particle_size_max"] = peat["fuels"]["rhein_braunkohle"]["particle_size_max"]
    with pytest.raises(CitedFieldNotApproved, match="particle_size_max"):
        read_field(peat, "5", "particle_size_max")
    assert read_case_field(VTT, "5", "particle_size_range") == (0.0, 0.004)

    wrong_unit = load_case_file(VTT)
    wrong_unit["shared"]["bed_diameter"]["unit"] = "K"
    with pytest.raises(CaseAssemblyError, match="bed_diameter"):
        read_field(wrong_unit, "1", "bed_diameter")
    text = load_case_file(VTT)
    text["shared"]["bed_diameter"]["value"] = "0.15"
    with pytest.raises(CaseAssemblyError, match="bed_diameter"):
        read_field(text, "1", "bed_diameter")
    not_a_flow = load_case_file(VTT)
    not_a_flow["shared"]["bed_diameter"]["unit"] = "m3(VN)/h"
    with pytest.raises(CaseAssemblyError, match="bed_diameter"):
        read_field(not_a_flow, "1", "bed_diameter")
