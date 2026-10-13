"""LU 验证 JSON 加载、单位换算与 inlet 拓扑抽取。

Source: ``data/validation_cases.json``；``docs/hamel_submodels/00_readme_and_citation_rules.md``
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from src.core.feed_inlet import compute_gas_feeds_mol_s

CASE_LU_VALIDATION_KEY = "CASE_HTW_WESSELING_1"

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _validation_cases_path() -> Path:
    return _REPO_ROOT / "data" / "validation_cases.json"


def _numeric_or_none(val: object) -> float | None:
    if isinstance(val, (int, float)):
        return float(val)
    return None


def nm3_h_to_mol_s_stp(flow_nm3_h: float) -> float:
    """Normvolumenstrom [Nm3/h] -> [mol/s]，按 STP 22.414 Nm3/kmol。"""
    return float(flow_nm3_h) * (1000.0 / 22.414) / 3600.0


def air_nm3_h_to_species_mol_s(flow_nm3_h: float) -> dict[str, float]:
    """空气 Normvolumenstrom 转成 `O2/N2` 摩尔流。"""
    n_tot = nm3_h_to_mol_s_stp(flow_nm3_h)
    return {
        "O2_mol_s": 0.21 * n_tot,
        "N2_mol_s": 0.79 * n_tot,
    }


def o2_nm3_h_to_species_mol_s(flow_nm3_h: float) -> dict[str, float]:
    """纯氧 Normvolumenstrom 转成 `O2/N2` 摩尔流。"""
    n_tot = nm3_h_to_mol_s_stp(flow_nm3_h)
    return {
        "O2_mol_s": n_tot,
        "N2_mol_s": 0.0,
    }


def _resolve_shared_ref(section: dict, root: dict | None, section_name: str | None = None) -> dict:
    if not isinstance(section, dict):
        return {}
    ref = section.get("_ref")
    if not ref:
        return dict(section)
    if root is None:
        root = load_validation_json_root()
    shared = root.get(str(ref), {})
    if not isinstance(shared, dict):
        return dict(section)
    if section_name and isinstance(shared.get(section_name), dict):
        shared = shared[section_name]
    merged = dict(shared)
    merged.update({k: v for k, v in section.items() if k != "_ref"})
    return merged


def extract_case_inlet_topology_hints(node: dict, root: dict | None = None) -> dict[str, object]:
    """从任意 validation case 节点中抽取 inlet / recycle / geometry 提示。"""
    ins = node.get("inputs", {})
    ga = ins.get("gasification_agent", {})
    react = _resolve_shared_ref(ins.get("reactor", {}), root, section_name="reactor")
    op = ins.get("operating_conditions", {})

    primary_air_nm3_h = _numeric_or_none(ga.get("primary_air_Nm3_h"))
    total_air_nm3_h = _numeric_or_none(ga.get("total_air_Nm3_h"))
    total_o2_nm3_h = _numeric_or_none(ga.get("total_O2_Nm3_h"))
    secondary_air_nm3_h = _numeric_or_none(ga.get("secondary_air_Nm3_h"))
    secondary_o2_nm3_h = _numeric_or_none(ga.get("secondary_O2_Nm3_h"))
    steam_feed_kg_h = _numeric_or_none(ga.get("steam_feed_kg_h"))
    secondary_xi = _numeric_or_none(ga.get("secondary_injection_xi_reactor"))
    secondary_agent = str(ga.get("secondary_injection_agent", "")).strip()

    if primary_air_nm3_h is None and total_air_nm3_h is not None and secondary_air_nm3_h is not None:
        primary_air_nm3_h = max(total_air_nm3_h - secondary_air_nm3_h, 0.0)

    primary_o2_nm3_h = _numeric_or_none(ga.get("primary_O2_Nm3_h"))
    if primary_o2_nm3_h is None and total_o2_nm3_h is not None and secondary_o2_nm3_h is not None:
        primary_o2_nm3_h = max(total_o2_nm3_h - secondary_o2_nm3_h, 0.0)

    primary_air_species = None if primary_air_nm3_h is None else air_nm3_h_to_species_mol_s(primary_air_nm3_h)
    primary_o2_species = None if primary_o2_nm3_h is None else o2_nm3_h_to_species_mol_s(primary_o2_nm3_h)
    secondary_air_species = None if secondary_air_nm3_h is None else air_nm3_h_to_species_mol_s(secondary_air_nm3_h)
    secondary_o2_species = None if secondary_o2_nm3_h is None else o2_nm3_h_to_species_mol_s(secondary_o2_nm3_h)

    secondary_agent_nm3_h = secondary_air_nm3_h if secondary_air_nm3_h is not None else secondary_o2_nm3_h
    secondary_agent_species = secondary_air_species if secondary_air_species is not None else secondary_o2_species

    bed_diameter = _numeric_or_none(react.get("bed_diameter_m"))
    freeboard_diameter = _numeric_or_none(react.get("freeboard_diameter_m"))
    legacy_diameter = _numeric_or_none(react.get("diameter_m"))
    if bed_diameter is None:
        bed_diameter = legacy_diameter
    if freeboard_diameter is None:
        freeboard_diameter = legacy_diameter

    bed_height = _numeric_or_none(react.get("bed_height_m"))
    reactor_height = _numeric_or_none(react.get("height_m")) or _numeric_or_none(react.get("total_height_m"))

    recirculation = bool(op.get("recirculation", False))
    recycle_topology = str(op.get("recycle_topology", "")).strip()
    recirculation_note = str(op.get("recirculation_note", "")).strip()
    cyclone_present = bool(op.get("cyclone_present", recirculation))
    if not recycle_topology and recirculation and "cyclone" in recirculation_note.lower():
        recycle_topology = "cyclone return"

    return {
        "plant": str(node.get("plant", "")),
        "secondary_injection_xi": secondary_xi,
        "secondary_injection_agent": secondary_agent,
        "secondary_injection_note": str(ga.get("secondary_injection_note", "")).strip(),
        "total_air_Nm3_h": total_air_nm3_h,
        "total_O2_Nm3_h": total_o2_nm3_h,
        "secondary_air_Nm3_h": secondary_air_nm3_h,
        "secondary_O2_Nm3_h": secondary_o2_nm3_h,
        "secondary_agent_Nm3_h": secondary_agent_nm3_h,
        "secondary_air_O2_mol_s": None if secondary_air_species is None else secondary_air_species["O2_mol_s"],
        "secondary_air_N2_mol_s": None if secondary_air_species is None else secondary_air_species["N2_mol_s"],
        "secondary_O2_O2_mol_s": None if secondary_o2_species is None else secondary_o2_species["O2_mol_s"],
        "secondary_O2_N2_mol_s": None if secondary_o2_species is None else secondary_o2_species["N2_mol_s"],
        "secondary_agent_O2_mol_s": None if secondary_agent_species is None else secondary_agent_species["O2_mol_s"],
        "secondary_agent_N2_mol_s": None if secondary_agent_species is None else secondary_agent_species["N2_mol_s"],
        "primary_air_Nm3_h": primary_air_nm3_h,
        "primary_O2_Nm3_h": primary_o2_nm3_h,
        "primary_air_O2_mol_s": None if primary_air_species is None else primary_air_species["O2_mol_s"],
        "primary_air_N2_mol_s": None if primary_air_species is None else primary_air_species["N2_mol_s"],
        "primary_O2_O2_mol_s": None if primary_o2_species is None else primary_o2_species["O2_mol_s"],
        "primary_O2_N2_mol_s": None if primary_o2_species is None else primary_o2_species["N2_mol_s"],
        "steam_feed_kg_h": steam_feed_kg_h,
        "recirculation": recirculation,
        "recirculation_note": recirculation_note,
        "recycle_topology": recycle_topology,
        "cyclone_present": cyclone_present,
        "reactor_height_m": reactor_height,
        "bed_height_m": bed_height,
        "bed_diameter_m": bed_diameter,
        "freeboard_diameter_m": freeboard_diameter,
        "secondary_air_port_mm_above_nozzle": _numeric_or_none(react.get("secondary_air_port_mm_above_nozzle")),
        "fuel_inlet_mm_above_nozzle": _numeric_or_none(react.get("fuel_inlet_mm_above_nozzle")),
    }


CASE_FLATTEN_FALLBACKS: dict[str, dict[str, float]] = {
    "CASE_VTT_PRESSURIZED_PEAT_6": {
        "moisture_wt_pct": 14.8,
        "T_inlet_K": 293.0,
        "bed_height_m": 1.2,
        "secondary_injection_xi_reactor": 0.6,
    },
}

MULTI_CASE_MATRIX_KEYS: tuple[str, ...] = (
    "CASE_HTW_WESSELING_1",
    "CASE_HTW_WESSELING_2",
    "CASE_VTT_PRESSURIZED_PEAT_6",
)


def _json_numeric_or_default(val: object, default: float) -> float:
    parsed = _numeric_or_none(val)
    return float(default) if parsed is None else float(parsed)


def _infer_primary_agent(agent_type: str) -> str:
    text = str(agent_type)
    if "O2" in text and "Steam" in text and "Air" not in text:
        return "o2_steam"
    if "Air" in text and "Steam" not in text:
        return "air"
    return "air_steam"


def _infer_fuel_type(fuel_label: str) -> str:
    text = str(fuel_label).lower()
    if any(token in text for token in ("peat", "torf", "wood", "sawdust", "biomass")):
        return "biomass"
    return "coal"


def _infer_reactor_family(plant: str, scale: str) -> str:
    plant_l = str(plant).lower()
    if "htw" in plant_l:
        return "HTW_pressurised"
    if "vtt" in plant_l:
        return "VTT_pressurised"
    if "wsv" in plant_l:
        return "WSV_atmospheric"
    return str(scale or "unknown")


def _infer_secondary_xi_from_axial(node: dict) -> float | None:
    axial = node.get("outputs", {}).get("axial_profiles", {})
    coord = str(axial.get("_coordinate", "")).lower()
    if "secondary" not in coord:
        return None
    for token in ("0.6", "0.60", "ξ≈0.6", "xi≈0.6"):
        if token.replace("ξ", "xi") in coord or token in coord:
            return 0.6
    return None


def flatten_validation_case(node: dict, *, case_key: str) -> dict:
    """将 validation_cases.json 中嵌套结构展平为 Reactor 测试所用扁平字段。"""
    ins = node["inputs"]
    fuel = ins["fuel"]
    pa = fuel["proximate_analysis"]
    ua = fuel["ultimate_analysis_dry_wt_pct"]
    react = _resolve_shared_ref(ins.get("reactor", {}), None, section_name="reactor")
    op = ins["operating_conditions"]
    ga = ins["gasification_agent"]
    hints = extract_case_inlet_topology_hints(node)
    fallbacks = CASE_FLATTEN_FALLBACKS.get(case_key, {})

    agent_type = str(ga.get("type", ""))
    primary_agent = _infer_primary_agent(agent_type)

    pd = fuel.get("particle_diameter_mm", {})
    if isinstance(pd, dict):
        particle_range = tuple(pd.get("mean_range", (1.5, 3.0)))
        n_size_classes = int(pd.get("n_size_classes", 10))
    else:
        particle_range = (1.5, 3.0)
        n_size_classes = 10

    bed_height_m = _json_numeric_or_default(
        react.get("bed_height_m"),
        float(fallbacks.get("bed_height_m", 6.0)),
    )
    reactor_height_m = _json_numeric_or_default(react.get("height_m"), bed_height_m)
    diameter_m = _json_numeric_or_default(
        react.get("diameter_m"),
        float(hints["freeboard_diameter_m"] or hints["bed_diameter_m"] or 0.6),
    )
    secondary_xi = hints["secondary_injection_xi"]
    if secondary_xi is None:
        secondary_xi = _json_numeric_or_default(
            fallbacks.get("secondary_injection_xi_reactor"),
            float("nan"),
        )
        if not np.isfinite(secondary_xi):
            secondary_xi = _infer_secondary_xi_from_axial(node)

    return {
        "_comment": node.get("_comment", ""),
        "_source_case_key": case_key,
        "reactor": _infer_reactor_family(str(node.get("plant", "")), str(node.get("scale", ""))),
        "fuel": str(fuel.get("type", "")),
        "fuel_type": _infer_fuel_type(str(fuel.get("type", ""))),
        "P": _json_numeric_or_default(op.get("pressure_MPa"), 0.101325) * 1e6,
        "T_inlet": _json_numeric_or_default(op.get("T_inlet_K"), float(fallbacks.get("T_inlet_K", 293.0))),
        "fuel_feed": float(fuel["feed_rate_kg_h"]),
        "primary_agent": primary_agent,
        "primary_agent_note": agent_type,
        "ER": float(op["ER"]),
        "recirculation": bool(hints["recirculation"]),
        "recirculation_note": str(hints["recirculation_note"]),
        "recycle_topology": str(hints["recycle_topology"]),
        "cyclone_present": bool(hints["cyclone_present"]),
        "H_bed": bed_height_m,
        "D_freeboard": diameter_m,
        "bed_diameter_m": float(hints["bed_diameter_m"] or diameter_m),
        "freeboard_diameter_m": float(hints["freeboard_diameter_m"] or diameter_m),
        "bed_height_m": bed_height_m,
        "reactor_height_m": reactor_height_m,
        "primary_air_Nm3_h": hints["primary_air_Nm3_h"],
        "total_air_Nm3_h": hints["total_air_Nm3_h"],
        "secondary_air_Nm3_h": hints["secondary_air_Nm3_h"],
        "primary_O2_Nm3_h": hints["primary_O2_Nm3_h"],
        "total_O2_Nm3_h": hints["total_O2_Nm3_h"],
        "secondary_O2_Nm3_h": hints["secondary_O2_Nm3_h"],
        "secondary_injection_xi": secondary_xi,
        "secondary_injection_agent": str(hints["secondary_injection_agent"]),
        "secondary_injection_note": str(hints["secondary_injection_note"]),
        "secondary_agent_Nm3_h": hints["secondary_agent_Nm3_h"],
        "secondary_agent_O2_mol_s": hints["secondary_agent_O2_mol_s"],
        "secondary_agent_N2_mol_s": hints["secondary_agent_N2_mol_s"],
        "steam_feed_kg_h": _numeric_or_none(ga.get("steam_feed_kg_h")),
        "particle_diameter_mm_range": particle_range,
        "n_size_classes": n_size_classes,
        "moisture_wt": _json_numeric_or_default(
            pa.get("moisture_wt_pct"),
            float(fallbacks.get("moisture_wt_pct", 15.0)),
        ),
        "ash_dry_wt": float(pa["ash_dry_wt_pct"]),
        "VM_daf": float(pa["VM_daf_pct"]),
        "C_dry": float(ua["C"]),
        "H_dry": float(ua["H"]),
        "O_dry": float(ua["O"]),
        "N_dry": float(ua["N"]),
        "S_dry": float(ua["S"]),
    }


def load_validation_case_flat(case_key: str = CASE_LU_VALIDATION_KEY) -> dict:
    """读取并展平任意 validation case（含 LU）。"""
    node = load_validation_case_node(case_key)
    return flatten_validation_case(node, case_key=case_key)


def load_case_LU() -> dict:
    return load_validation_case_flat(CASE_LU_VALIDATION_KEY)


def load_validation_json_root() -> dict:
    """读取 validation_cases.json 根对象（含 validation_performance_summary）。"""
    data_path = _validation_cases_path()
    with open(data_path, encoding="utf-8") as f:
        return json.load(f)


def load_validation_case_node(key: str = CASE_LU_VALIDATION_KEY) -> dict:
    """读取 validation_cases.json 中完整 CASE 节点（含 inputs / outputs）。"""
    data = load_validation_json_root()
    if key not in data:
        raise KeyError(f"{_validation_cases_path()} 中缺少键 {key!r}")
    return data[key]


def json_numeric_or_none(val: object) -> float | None:
    if isinstance(val, (int, float)):
        return float(val)
    return None


def validation_numeric_tolerances() -> dict[str, float]:
    """与根节点 validation_performance_summary 对齐的相对容差（解析为标量）。"""
    return {
        "rtol_T": 0.15,
        "rtol_CO_CO2_H2": 0.15,
        "rtol_CH4": 0.20,
        "rtol_carbon_conv": 0.15,
    }


def strict_validation_gate(
    result: dict[str, Any],
    *,
    rms_max: float = 0.15,
) -> tuple[bool, list[str]]:
    """Return whether a solved state is trustworthy enough for strict validation."""
    reasons: list[str] = []
    converged = bool(
        result.get(
            "converged_fully",
            result.get("converged_outer", result.get("converged", False)),
        )
    )
    if not converged:
        reasons.append("solver_not_converged")

    rms_raw = result.get("rms_scaled_final")
    if rms_raw is None:
        reasons.append("missing_rms_scaled_final")
    else:
        rms = float(rms_raw)
        if not np.isfinite(rms):
            reasons.append("nonfinite_rms_scaled_final")
        elif rms > float(rms_max):
            reasons.append(f"rms_scaled_final>{float(rms_max):.3f}")

    component_raw = result.get("rms_scaled_component_max_final")
    if component_raw is not None:
        component = float(component_raw)
        if not np.isfinite(component):
            reasons.append("nonfinite_rms_scaled_component_max_final")
        elif component > float(rms_max):
            reasons.append(f"rms_scaled_component_max_final>{float(rms_max):.3f}")

    max_abs_raw = result.get("max_abs_scaled_final")
    if max_abs_raw is not None:
        max_abs = float(max_abs_raw)
        if not np.isfinite(max_abs):
            reasons.append("nonfinite_max_abs_scaled_final")
        elif max_abs > float(rms_max):
            reasons.append(f"max_abs_scaled_final>{float(rms_max):.3f}")

    return len(reasons) == 0, reasons


def estimate_gas_feeds(case: dict) -> dict[str, float]:
    """与 ``ReactorConfig(ER=..., primary_agent=...)`` / ``compute_gas_feeds_mol_s`` 一致。"""
    fuel_kg_s = case["fuel_feed"] / 3600.0
    o2, h2o, n2 = compute_gas_feeds_mol_s(
        fuel_feed_kg_s=fuel_kg_s,
        moisture_wt=case["moisture_wt"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
    )
    return {"O2_feed": o2, "H2O_feed": h2o, "N2_feed": n2}
