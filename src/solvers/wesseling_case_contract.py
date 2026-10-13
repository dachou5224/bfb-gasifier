"""Wesseling 工况输入合同与外部元素账。不读取求解残差。

缺失量保持缺失。数字化 5% 只描述读图误差，不是模型验收容差。

Source: data/validation_cases.json；docs/audit_reports/plan_ab_diagnosis_2026-09-28.md 工作包 3
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_ROOT = Path(__file__).resolve().parents[2]
_CASES = _ROOT / "data" / "validation_cases.json"
_STP_M3_PER_KMOL = 22.414  # [Nm³/kmol]
_KG_PER_G = 1.0e-3


def _status(value: Any) -> str:
    text = "" if value is None else str(value)
    if text.startswith("__UNVERIFIED__"):
        return "unverified"
    if text.startswith("__MISSING__") or text.startswith("__NOT_TABULATED__") or value is None:
        return "missing"
    if text.startswith("__CALCULATED__"):
        return "calculated"
    return "tabulated"


def load_case(case_id: str) -> dict[str, Any]:
    data = json.loads(_CASES.read_text(encoding="utf-8"))
    if case_id not in data:
        raise KeyError(case_id)
    return data[case_id]


def map_wesseling_inputs(case_id: str) -> dict[str, Any]:
    """把工况 JSON 收成带状态的输入。不补默认床料、蒸汽或一次风份额。"""
    case = load_case(case_id)
    fuel = case["inputs"]["fuel"]
    agent = case["inputs"]["gasification_agent"]
    reactor = case["inputs"]["reactor"]
    op = case["inputs"]["operating_conditions"]
    outputs = case["outputs"]
    steam = agent.get("steam_flow_kg_h", agent.get("steam_feed_kg_h"))
    rows = {
        "fuel_as_received_kg_h": fuel["feed_rate_kg_h"],
        "moisture_wt_pct": fuel["proximate_analysis"]["moisture_wt_pct"],
        "ash_dry_wt_pct": fuel["proximate_analysis"]["ash_dry_wt_pct"],
        "vm_daf_pct": fuel["proximate_analysis"]["VM_daf_pct"],
        "ultimate_dry_wt_pct": fuel["ultimate_analysis_dry_wt_pct"],
        "particle_diameter_mm": fuel["particle_diameter_mm"]["mean_range"],
        "pressure_MPa": op["pressure_MPa"],
        "bed_height_m": reactor["bed_height_m"],
        "diameter_m": reactor["diameter_m"],
        "bed_material": reactor.get("bed_material"),
        "bed_particle_diameter_mm": reactor.get("bed_particle_diameter_mm"),
        "heat_loss_note": case["inputs"].get("heat_loss_note"),
        "steam_kg_h": steam,
        "recycle": op.get("recirculation"),
    }
    if "total_air_Nm3_h" in agent:
        rows["total_air_Nm3_h"] = agent["total_air_Nm3_h"]
        rows["secondary_air_Nm3_h"] = agent.get("secondary_air_Nm3_h")
        rows["primary_air_Nm3_h"] = agent.get("primary_air_Nm3_h")
    if "total_O2_Nm3_h" in agent:
        rows["total_O2_Nm3_h"] = agent["total_O2_Nm3_h"]
        rows["secondary_O2_Nm3_h"] = agent.get("secondary_O2_Nm3_h")
        rows["primary_O2_Nm3_h"] = agent.get("primary_O2_Nm3_h")
    mapped = []
    for key, value in rows.items():
        mapped.append({"name": key, "value": value, "status": _status(value)})
    return {
        "case_id": case_id,
        "conditional": any(item["status"] in ("missing", "unverified") for item in mapped),
        "fields": mapped,
        "observation": {
            "basis": "wet mole fraction for the tabulated exit gas",
            "paper_simulation_temperature_K": outputs.get("exit_temperature_K"),
            "measured_temperature_K": outputs.get("exit_temperature_measured_K"),
            "digitization_rel": 0.05,
            "digitization_source": outputs.get("_precision"),
            "measurement_uncertainty": "not_declared_in_table_7_1",
            "model_allowance": "not_set_from_results",
        },
    }


def fuel_element_kmol_h(case_id: str) -> dict[str, float]:
    """干基元素加水分中的 H、O。不含未给定的蒸汽。"""
    case = load_case(case_id)
    fuel = case["inputs"]["fuel"]
    wet = float(fuel["feed_rate_kg_h"])
    moisture = float(fuel["proximate_analysis"]["moisture_wt_pct"]) / 100.0
    dry = wet * (1.0 - moisture)
    ult = fuel["ultimate_analysis_dry_wt_pct"]
    molar = {"C": 12.011, "H": 1.00794, "O": 15.999, "N": 14.0067, "S": 32.065}
    out = {el: dry * float(ult[el]) / 100.0 / molar[el] for el in molar}
    water = wet * moisture / 18.01528
    out["H"] += 2.0 * water
    out["O"] += water
    out["dry_fuel_kg_h"] = dry
    return out


def air_kmol_h(total_nm3_h: float, o2_nm3_h: float) -> dict[str, float]:
    """由总体积和氧体积得到 O2、N2。不把二次风再加一次。"""
    o2 = float(o2_nm3_h) / _STP_M3_PER_KMOL
    n2 = (float(total_nm3_h) - float(o2_nm3_h)) / _STP_M3_PER_KMOL
    return {"O2_kmol_h": o2, "N2_kmol_h": n2}


def external_element_balance_kmol_h(
    inlet: dict[str, float],
    outlet: dict[str, float],
    recycle: dict[str, float] | None = None,
) -> dict[str, float]:
    """外部进出之差。循环项必须成对抵消，不进入这个差。"""
    loop = recycle or {}
    if loop and abs(sum(loop.values())) > 1.0e-9:
        raise ValueError("recycle streams must sum to zero")
    keys = set(inlet) | set(outlet)
    return {key: float(inlet.get(key, 0.0)) - float(outlet.get(key, 0.0)) for key in sorted(keys)}


def external_enthalpy_w(
    inlet_w: float,
    outlet_w: float,
    recycle_up_w: float = 0.0,
    recycle_down_w: float = 0.0,
) -> float:
    """外部焓差 [W]。向上和向下的循环焓必须等大反号，不留下净项。"""
    if abs(float(recycle_up_w) + float(recycle_down_w)) > 1.0e-6:
        raise ValueError("recycle enthalpy must cancel")
    return float(inlet_w) - float(outlet_w)


def wesseling1_external_inlet_kmol_h() -> dict[str, float | bool]:
    """燃料加总空气的外部进料。蒸汽和一次/二次分配没有加进去。"""
    fuel = fuel_element_kmol_h("CASE_HTW_WESSELING_1")
    case = load_case("CASE_HTW_WESSELING_1")
    agent = case["inputs"]["gasification_agent"]
    calc = agent["calculated_flows"]
    air = air_kmol_h(float(agent["total_air_Nm3_h"]), float(calc["O2_actual_Nm3_h"]))
    inlet = {el: float(fuel[el]) for el in ("C", "H", "O", "N", "S")}
    inlet["O"] += 2.0 * air["O2_kmol_h"]
    inlet["N"] += 2.0 * air["N2_kmol_h"]
    inlet["steam_included"] = False
    inlet["primary_split_included"] = False
    inlet["O2_kmol_h"] = air["O2_kmol_h"]
    inlet["N2_kmol_h"] = air["N2_kmol_h"]
    return inlet


def wesseling1_inventory_and_energy() -> dict[str, Any]:
    """固定床料和燃料灰是两套库存。热值缺失时不算热损失功率。"""
    case = load_case("CASE_HTW_WESSELING_1")
    fuel = fuel_element_kmol_h("CASE_HTW_WESSELING_1")
    ash_frac = float(case["inputs"]["fuel"]["proximate_analysis"]["ash_dry_wt_pct"]) / 100.0
    return {
        "reactive_ash_from_fuel_kg_h": float(fuel["dry_fuel_kg_h"]) * ash_frac,
        "reactive_ash_basis": "dry fuel ash fraction; this ash can leave with the gas and solids",
        "fixed_bed_inventory_kg": None,
        "fixed_bed_basis": "not fuel ash and not the named-kernel geometric gap",
        "heat_loss_fraction_range": (0.05, 0.10),
        "fuel_heating_value": None,
        "heat_loss_W": None,
        "recycle_topology": "cyclone return to bed",
        "recycle_mass_flow_kg_h": None,
        "exit_mole_fractions_are_not_an_element_balance": True,
    }


def bed_phase_volumes_m3(
    area_m2: float,
    h_bed_m: float,
    eps_b: float,
    eps_d: float,
) -> dict[str, float]:
    """床体积按相含率拆开。``V_d`` 是悬浮相里的气体体积，不是悬浮相总体积。

    ``eps_d`` 是悬浮相内空隙率，``eps_b`` 是床层气泡相含率。
    恒等式：气泡 + 悬浮相气体 + 悬浮相固体 = A H。
    """
    v_bed = float(area_m2) * float(h_bed_m)
    v_b = float(eps_b) * v_bed
    v_emulsion = (1.0 - float(eps_b)) * v_bed
    v_d = float(eps_d) * v_emulsion
    v_solids = (1.0 - float(eps_d)) * v_emulsion
    return {
        "V_bed_m3": v_bed,
        "V_b_m3": v_b,
        "V_emulsion_m3": v_emulsion,
        "V_d_m3": v_d,
        "V_solids_m3": v_solids,
    }


def wesseling1_conditional_boundary() -> dict[str, Any]:
    """Wesseling 1 两区条件边界。缺的蒸汽、一次风和床料不从出口回填。

    全部表载空气进气泡入口。二次风仍留在这股总风里，两区没有 ξ≈0.6 的喷口。
    相体积由 Wesseling 的 A、H 和相含率拆出。含率本身不是表 7.1 的测量值。
    """
    case = load_case("CASE_HTW_WESSELING_1")
    fuel = case["inputs"]["fuel"]
    wet = float(fuel["feed_rate_kg_h"])
    moisture = float(fuel["proximate_analysis"]["moisture_wt_pct"]) / 100.0
    dry = wet * (1.0 - moisture)
    ult = fuel["ultimate_analysis_dry_wt_pct"]
    molar = {"C": 12.011, "H": 1.00794, "O": 15.999, "N": 14.0067, "S": 32.065}
    dry_kmol_h = {el: dry * float(ult[el]) / 100.0 / molar[el] for el in molar}
    inlet = wesseling1_external_inlet_kmol_h()
    diameter = float(case["inputs"]["reactor"]["diameter_m"])
    area = 3.141592653589793 * (0.5 * diameter) ** 2
    h_bed = float(case["inputs"]["reactor"]["bed_height_m"])
    h_total = float(case["inputs"]["reactor"]["total_height_m"])
    eps_b = 0.13  # [-] 具名核冻结气泡含率，表 7.1 没有
    eps_d = 0.45  # [-] 具名核冻结悬浮相空隙率，表 7.1 没有
    volumes = bed_phase_volumes_m3(area, h_bed, eps_b, eps_d)
    vm = float(fuel["proximate_analysis"]["VM_daf_pct"]) / 100.0
    return {
        "conditional": True,
        "steam_kg_h": 0.0,
        "primary_split_included": False,
        "secondary_air_placed_at_xi": False,
        "two_zone_scores_secondary_temperature_peak": False,
        "dry_element_kmol_h": dry_kmol_h,
        "moisture_kg_h": wet * moisture,
        "O2_mol_s": float(inlet["O2_kmol_h"]) / 3.6,
        "N2_mol_s": float(inlet["N2_kmol_h"]) / 3.6,
        "diameter_m": diameter,
        "area_m2": area,
        "H_bed_m": h_bed,
        "H_bed_status": "section_7_1_2_note",
        "H_freeboard_m": h_total - h_bed,
        "pressure_Pa": float(case["inputs"]["operating_conditions"]["pressure_MPa"]) * 1.0e6,
        "T_in_K": float(case["inputs"]["operating_conditions"]["T_inlet_K"]),
        "eps_b": eps_b,
        "eps_d": eps_d,
        "voidage_status": "kernel_default_not_tabulated",
        "V_bed_m3": volumes["V_bed_m3"],
        "V_d_m3": volumes["V_d_m3"],
        "V_b_m3": volumes["V_b_m3"],
        "V_solids_m3": volumes["V_solids_m3"],
        "char_diameter_m": (2.0e-4, 3.0e-4),
        "char_diameter_basis": "named-kernel bins, not the 1.5–3 mm fuel and not the missing bed material",
        "formation_enthalpy_basis": "named-kernel schematic dry-fuel formation enthalpy, not an HHV or LHV",
        "Q_W": 0.0,
        "launch_window_s": 5.0,
        "paper_simulation_temperature_K": 1100,
        "measured_temperature_K": 1120,
        "digitization_rel": 0.05,
        "model_allowance": "not_set_from_results",
        "nC_fixed_share_of_fuel_carbon": 1.0 - vm,
        "nC_vm_share_of_fuel_carbon": vm,
        "ash_dry_wt_pct": float(fuel["proximate_analysis"]["ash_dry_wt_pct"]),
        "VM_daf_pct": float(fuel["proximate_analysis"]["VM_daf_pct"]),
        "secondary_air_Nm3_h": None,
        "secondary_air_Nm3_h_historical": 30.1,
        "secondary_air_xi": None,
        "secondary_air_xi_digitized": 0.6,
        "fuel_size_mm": {
            "mean_range": list(fuel["particle_diameter_mm"]["mean_range"]),
            "max_feed": float(fuel["particle_diameter_mm"]["max_feed"]),
            "n_size_classes": int(fuel["particle_diameter_mm"]["n_size_classes"]),
            "feed_psd_measured": False,
        },
        "steam_status": "missing",
        "steam_kg_h_from_table": None,
        "steam_temperature_K": None,
        "numerical_target_steam_kg_h": 0.0,
        "complete_air_steam_case": False,
        "steam_source": "Tabelle 7.1 lists Luft/Dampf and lambda_O2, not steam kg/h",
        "bed_material_status": "missing",
        "recycle_mass_flow_status": "missing",
    }


def wesseling12_input_ledger() -> dict[str, Any]:
    """Wesseling 1/2 的共同设备、不同气化剂，以及仍空着的热值和粒级质量。

    数值直接来自 ``data/validation_cases.json``。热值和各档质量保持空，
    不用出口温度或组成去填。
    """
    rows = {}
    for case_id in ("CASE_HTW_WESSELING_1", "CASE_HTW_WESSELING_2"):
        case = load_case(case_id)
        fuel = case["inputs"]["fuel"]
        agent = case["inputs"]["gasification_agent"]
        reactor = case["inputs"]["reactor"]
        op = case["inputs"]["operating_conditions"]
        size = fuel["particle_diameter_mm"]
        steam = agent.get("steam_flow_kg_h", agent.get("steam_feed_kg_h"))
        rows[case_id] = {
            "simulation": case.get("hamel_simulation_id"),
            "table": case.get("hamel_table_ref"),
            "fuel": fuel["type"],
            "feed_rate_kg_h": float(fuel["feed_rate_kg_h"]),
            "moisture_wt_pct": float(fuel["proximate_analysis"]["moisture_wt_pct"]),
            "ash_dry_wt_pct": float(fuel["proximate_analysis"]["ash_dry_wt_pct"]),
            "VM_daf_pct": float(fuel["proximate_analysis"]["VM_daf_pct"]),
            "ultimate_dry_wt_pct": {
                el: float(fuel["ultimate_analysis_dry_wt_pct"][el]) for el in ("C", "H", "O", "N", "S")
            },
            "agent": agent["type"],
            "lambda_O2": float(op["ER"]),
            "pressure_MPa": float(op["pressure_MPa"]),
            "T_inlet_K": float(op["T_inlet_K"]),
            "diameter_m": float(reactor["diameter_m"]),
            "bed_height_m": float(reactor["bed_height_m"]),
            "total_height_m": float(reactor["total_height_m"]),
            "secondary_xi_measured": None,
            "secondary_xi_digitized": agent.get("secondary_injection_xi_digitized"),
            "secondary_agent": agent["secondary_injection_agent"],
            "secondary_Nm3_h_measured": None,
            "secondary_Nm3_h_historical": agent.get("secondary_air_Nm3_h_historical", agent.get("secondary_O2_Nm3_h_historical")),
            "steam_status": _status(steam),
            "primary_split_status": "missing",
            "paper_simulation_temperature_K": case["outputs"].get("exit_temperature_K"),
            "measured_temperature_K": case["outputs"].get("exit_temperature_measured_K"),
        }
    return {
        "plant": "HTW Wesseling",
        "common_equipment": {
            "diameter_m": 0.6,
            "bed_height_m": 6.0,
            "total_height_m": 14.5,
            "pressure_MPa": 2.5,
            "T_inlet_K": 293,
            "recycle": "cyclone return to bed",
        },
        "cases": rows,
        "heating_value": {
            "dry_or_as_received": None,
            "hhv_or_lhv": None,
            "value": None,
            "unit": None,
            "present_in_table_7_1_ocr": False,
            "required_by": "Hamel Tabelle 2.2 lists Heizwert as an operating input",
            "formation_enthalpy_needs_dry_HHV": True,
            "heat_loss_frac_of_thermal_input": (0.05, 0.10),
            "thermal_input_W": None,
            "heat_loss_W": None,
            "named_kernel_h_f_dry_J_kg": -3.96e6,
            "named_kernel_h_f_scope": "schematic dry-fuel formation enthalpy, not an HHV or LHV, and not a substitute for Q_W",
            "fitted_from_exit": False,
        },
        "fuel_size": {
            "wesseling_feed_class_mass_kg_h": None,
            "wesseling_feed_psd_in_table_7_1": False,
            "model_discretization_n_classes": 10,
            "model_discretization_source": "Hamel §2.2 adopts the ten char size classes of Dersch (1994); this is not a measured Wesseling feed",
            "table_A_2_1_1": "laboratory combustion batches of dry Rhenish brown coal at single diameters, including 1.5 mm and 3.0 mm",
            "lab_coal_max_mm": 6.0,
            "lab_coal_source": "Hamel §6.1.2, the laboratory bed, not the Wesseling pilot feed",
            "applied_to_char_bins": False,
            "char_bins_m": [2.0e-4, 3.0e-4],
            "char_bins_scope": "named-kernel diameters",
        },
        "agent_closure": {
            "steam_kg_h_in_table_7_1_ocr": None,
            "steam_temperature_K": None,
            "injection_levels_footnote": 10,
            "injection_levels_source": "Tabelle 7.1 footnote: the agent enters on up to 10 levels",
            "freeboard_agent_mentioned": True,
            "freeboard_agent_source": "Hamel §7.1.2 narrative for Bild 7.4: an air/steam mixture is added in the freeboard",
            "secondary_Nm3_h_in_table_7_1_ocr": None,
            "secondary_Nm3_h_in_case_json": {
                "CASE_HTW_WESSELING_1": 30.1,
                "CASE_HTW_WESSELING_2": 37.5,
            },
            "secondary_Nm3_h_used_as_measurement": False,
            "pdf_visual_check": "Tabelle 7.1 on printed pages 128-129 has no steam kg/h, heating value, sieve masses, or 30.1/37.5 Nm3/h",
            "pdf_path": "docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf",
            "table_caption_citation": "Hamel et al. (1999), 4th ISCS Beijing, pp. 411-420; full text is not in this workspace",
            "primary_measurement_source": "Rheinbraun AG and Krupp Uhde GmbH, internal report 1997; not public",
            "public_background_sources_supply_the_missing_numbers": False,
            "xi_0_6_meaning": "historical digitized value; not re-digitized from Bild 7.4 in this pass and not a nozzle coordinate",
            "two_zone_can_place_freeboard_injection": False,
            "numerical_target_puts_lambda_air_at_bubble_inlet": True,
            "numerical_target_steam_kg_h": 0.0,
            "numerical_target_steam_is_an_assumption": True,
            "complete_case": False,
            "experimental_comparison_blocked": True,
            "experimental_order": "conditional Wesseling 1 numerical root, then real steam/heating-value/sieve/staged-injection boundaries and the multi-cell map, then the experiment; Wesseling 2 uses the same model afterwards",
        },
        "bypass_accounting": {
            "species": ("NH3", "H2S"),
            "enters_external_atoms_and_enthalpy": True,
            "enters_seven_species_molar_total": False,
            "enters_gas_velocity_density_and_mole_fraction": False,
            "enters_particle_gas_field": False,
            "subsequent_reactions": False,
            "nitrogen_sulfur_reaction_network": False,
        },
        "experimental_comparison_done": False,
        "experimental_comparison_withheld_until": "Wesseling 1 target root, independent conservation, and input closure",
    }


def wesseling12_conditional_definitions() -> dict[str, Any]:
    """空气/蒸汽 Wesseling 1 的输入定义，以及同模型下的 Wesseling 2 对照。

    表载量来自印刷页 128。缺项保持空，数值假设单独列出。
    """
    ledger = wesseling12_input_ledger()
    air = ledger["cases"]["CASE_HTW_WESSELING_1"]
    oxygen = ledger["cases"]["CASE_HTW_WESSELING_2"]
    return {
        "priority_case": "CASE_HTW_WESSELING_1",
        "comparison_case": "CASE_HTW_WESSELING_2",
        "same_model_and_solver": True,
        "two_zone_certification_is_not_the_multicell_map": True,
        "table_fields": {
            "CASE_HTW_WESSELING_1": {
                "agent": air["agent"],
                "pressure_MPa": air["pressure_MPa"],
                "feed_rate_kg_h": air["feed_rate_kg_h"],
                "lambda_O2": air["lambda_O2"],
                "moisture_wt_pct": air["moisture_wt_pct"],
                "ultimate_dry_wt_pct": air["ultimate_dry_wt_pct"],
            },
            "CASE_HTW_WESSELING_2": {
                "agent": oxygen["agent"],
                "pressure_MPa": oxygen["pressure_MPa"],
                "feed_rate_kg_h": oxygen["feed_rate_kg_h"],
                "lambda_O2": oxygen["lambda_O2"],
                "moisture_wt_pct": oxygen["moisture_wt_pct"],
                "ultimate_dry_wt_pct": oxygen["ultimate_dry_wt_pct"],
            },
        },
        "calculated_not_metered": {
            "CASE_HTW_WESSELING_1_total_air_Nm3_h": 5536.9,
            "CASE_HTW_WESSELING_2_total_O2_Nm3_h": 1921.6,
            "basis": "lambda_O2 times stoichiometric oxygen, not a flowmeter reading",
        },
        "unavailable_not_filled": (
            "steam_kg_h",
            "steam_temperature",
            "heating_value_basis_and_value",
            "sieve_masses",
            "secondary_flow_Nm3_h",
            "nozzle_coordinate",
        ),
        "numerical_assumptions": {
            "steam_kg_h": 0.0,
            "Q_W": 0.0,
            "h_f_dry_J_kg": -3.96e6,
            "all_lambda_air_at_bubble_inlet": True,
            "char_bins_m": (2.0e-4, 3.0e-4),
        },
        "sources_checked_this_pass": {
            "table_7_1_printed_pages": (128, 129),
            "hamel_et_al_1999_iscs_full_text": "not_in_workspace",
            "rheinbraun_krupp_uhde_1997": "internal_report_not_public",
            "schumacher_1993_braunkohle_abstract": "no steam kg/h, heating value, or sieve for simulation 1",
            "netl_gtc9808p": "Wesseling pilot capacity and air or oxygen/steam operation; not the missing table columns",
        },
        "sources_may_still_exist": True,
        "experimental_comparison_done": False,
    }


def aa3_exit() -> dict[str, Any]:
    """数值路线的出口，以及五项实验输入的来源和缺项。

    不从出口反推，不把数值假设写成测量。
    """
    definition = wesseling12_conditional_definitions()
    ledger = wesseling12_input_ledger()
    heat = ledger["heating_value"]
    size = ledger["fuel_size"]
    agent = ledger["agent_closure"]
    return {
        "last_joint_pass_s": 0.418,
        "last_joint_pass_file": "docs/audit_reports/z1_z3_review_2026-10-08/aa1_aa2.json",
        "last_joint_pass_is_a_homotopy_root": True,
        "route": (
            "conditional s=1",
            "repeatable convergence and cost at fixed inputs",
            "full-furnace discretization with transport, element and zone-energy ledgers",
            "paired Wesseling 1 air/steam and Wesseling 2 oxygen/steam on the same model",
        ),
        "named_kernel_s1_is_only_the_first_milestone": True,
        "milestone_reached": {
            "conditional_s1": False,
            "fixed_input_cost": False,
            "full_furnace_ledger": False,
            "paired_experiment": False,
        },
        "same_model_and_solver": True,
        "inputs": {
            "steam": {
                "source": "Tabelle 7.1 names Luft/Dampf or O2/Dampf; printed pages 128-129 give no kg/h or temperature",
                "measured": None,
                "numerical_assumption": {"steam_kg_h": 0.0},
                "assumption_is_a_measurement": False,
            },
            "heating_value_and_formation_enthalpy": {
                "source": heat["required_by"],
                "measured": None,
                "basis": None,
                "numerical_assumption": {
                    "h_f_dry_J_kg": heat["named_kernel_h_f_dry_J_kg"],
                    "scope": heat["named_kernel_h_f_scope"],
                },
                "assumption_is_a_measurement": False,
            },
            "feed_size": {
                "source": size["model_discretization_source"],
                "measured_sieve_mass_kg_h": size["wesseling_feed_class_mass_kg_h"],
                "numerical_assumption": {"char_bins_m": list(size["char_bins_m"])},
                "assumption_is_a_measurement": False,
            },
            "staged_injection": {
                "source": agent["injection_levels_source"],
                "freeboard_note": agent["freeboard_agent_source"],
                "measured_secondary_Nm3_h": None,
                "historical_not_used": {
                    "CASE_HTW_WESSELING_1_Nm3_h": 30.1,
                    "CASE_HTW_WESSELING_2_Nm3_h": 37.5,
                    "xi_digitized": 0.6,
                },
                "two_zone_can_place_it": False,
                "numerical_assumption": "all lambda-derived gas at the bubble inlet",
                "assumption_is_a_measurement": False,
            },
            "heat_loss": {
                "source": "Hamel Tabelle 2.2 gives 0.05 to 0.10 of thermal input; thermal input needs the missing heating value",
                "measured_W": None,
                "fraction_range": heat["heat_loss_frac_of_thermal_input"],
                "numerical_assumption": {"Q_W": 0.0},
                "assumption_is_a_measurement": False,
            },
        },
        "fitted_from_exit": False,
        "sources_may_still_exist": definition["sources_may_still_exist"],
        "experimental_comparison_done": False,
    }


def ejection_and_wesseling_open_items() -> dict[str, Any]:
    """抛射链和 Wesseling 输入的来源与未决项。不改速率，不把缺项写成默认值。"""
    bound = wesseling1_conditional_boundary()
    return {
        "eq_3_62": "P''=(1-cos xi_j) * integral of f(dl/d_b); printed page 39",
        "eq_3_63_literal": "P''_k=(1-P'')*P''_{k-1}",
        "eq_3_63_rewritten_as_geometric_power": False,
        "P_double_prime_1": None,
        "approach_angle_deg": {"clift_grace": 30.0, "briens": 48.0, "selected": None},
        "n_B_cap": 5,
        "zeta_w_pemberton": 0.4,
        "george_grace_ejected_wake_fraction": 0.4,
        "live_hamel_extended_uses_eq_367": False,
        "fuel_formation_enthalpy_basis": bound["formation_enthalpy_basis"],
        "fuel_size_applied_to_bins": False,
        "steam_status": bound["steam_status"],
        "secondary_air_placed": False,
        "secondary_air_Nm3_h": bound["secondary_air_Nm3_h"],
        "secondary_air_Nm3_h_historical": bound["secondary_air_Nm3_h_historical"],
        "secondary_air_xi": bound["secondary_air_xi"],
        "secondary_air_xi_digitized": bound["secondary_air_xi_digitized"],
    }


def element_pair_cancels(
    up: dict[str, float],
    down: dict[str, float],
    *,
    tol: float = 1.0e-9,
) -> bool:
    """同一界面的上下流必须按元素分别抵消。元素和为零不够。"""
    keys = set(up) | set(down)
    return all(abs(float(up.get(key, 0.0)) + float(down.get(key, 0.0))) <= tol for key in keys)
