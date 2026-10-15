"""通用地装配 ``data/cases/*.json``。

同一套键名与 ``status`` 语义。``status=missing`` 报错并指出字段名，不填默认值。
``status=cited`` 只有 ADR 0005、0006 批准进入模型的字段才可装配。
体积流量按 ADR 0004 换算为 mol/s。本模块不写物理方程。
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from hamel.cases.config import PhysicalConfig, SolverConfig

# 解释性选择：ADR 0004。DIN 1343，V_M = 0,022414 m³/mol。
V_M_M3_PER_MOL = 0.022414
_HOUR_S = 3600.0

# ADR 0005：批准为模型输入的 cited 字段。
_CITED_MODEL = {
    "fuel_feed_height",
    "bed_zone_height",
    "freeboard_height",
    "fixed_carbon_mass_fraction_wf",
    "particle_size_max",
    "particle_size_range",
    "inlet_temperature",
    "heating_value",
}
# ADR 0005：只这些工况的进气温度与收到基低位热值。
_CITED_CASES = {
    "inlet_temperature": {"1", "3", "7", "11", "14"},
    "heating_value": {"1", "3", "7", "11", "14"},
}
# ADR 0006：仅工况 1、2 把吹扫氮计入模型。
_NITROGEN_CASES = {"1", "2"}
# ADR 0005：文献里有数，但不作模型输入。
_RECORD_ONLY = {
    "distributor_to_outlet_height",
    "pressure_vessel_height",
    "pressure_vessel_inside_diameter",
    "ash_sample_pipe_diameter",
    "thermal_design_power",
}
_AMOUNT_FROM_VOLUME = {
    "primary_air_amount_flow": "primary_air_volume_flow",
    "secondary_air_amount_flow": "secondary_air_volume_flow",
}
_INT_FIELDS = {"cyclone_cell_count", "fuel_feed_count", "agent_feed_count"}
_SI_UNITS = {"m", "kg/s", "kg/kg", "mol/s", "Pa", "K", "J/kg", "1", "W", "mol/mol"}
_NORMAL_VOLUME = "m3(VN)/h"

# 与 PhysicalConfig 的字段一致。第 1 关 VTT 加压 1 都能取到。
PHYSICAL_FIELDS = (
    "bed_diameter",
    "freeboard_diameter",
    "reactor_height",
    "bed_zone_height",
    "freeboard_height",
    "fuel_feed_height",
    "cyclone_cell_count",
    "fuel_feed_count",
    "agent_feed_count",
    "ash_mass_fraction_wf",
    "volatiles_mass_fraction_waf",
    "carbon_mass_fraction_wf",
    "hydrogen_mass_fraction_wf",
    "oxygen_mass_fraction_wf",
    "nitrogen_mass_fraction_wf",
    "sulfur_mass_fraction_wf",
    "fixed_carbon_mass_fraction_wf",
    "particle_size_max",
    "fuel_mass_flow",
    "moisture_mass_fraction",
    "lambda_o2",
    "primary_air_amount_flow",
    "secondary_air_amount_flow",
    "nitrogen_purge_volume_flow",
    "steam_mass_flow",
    "inlet_temperature",
    "pressure",
    "recirculation",
    "heating_value",
)

class CaseAssemblyError(Exception):
    """装配失败。``field`` 是工况文件里的键名。"""

    def __init__(self, field: str, message: str):
        self.field = field
        super().__init__(f"{field}：{message}")

class MissingFieldError(CaseAssemblyError):
    """``status=missing``，或值为空。不填默认值。"""

class CitedFieldNotApproved(CaseAssemblyError):
    """``status=cited``，但没有 ADR 批准把它配进模型。"""

class RecordOnlyFieldError(CaseAssemblyError):
    """ADR 0005 写明只作记录、不作模型输入。"""

class FieldNotInPlant(CaseAssemblyError):
    """这个装置的工况文件里没有该键。"""

@dataclass(frozen=True)
class CaseConfig:
    """一次装配的结果。物理参数与求解器参数分开。"""

    plant_id: str
    case_key: str
    case_id: str
    physical: PhysicalConfig
    solver: SolverConfig

def normal_volume_to_mol_s(volume_m3_per_h: float) -> float:
    """把 m³(V_N)/h 换成 mol/s。

    解释性选择：ADR 0004。Ṅ = V̇_N / V_M / 3600。
    """
    return float(volume_m3_per_h) / V_M_M3_PER_MOL / _HOUR_S

def load_case_file(path: Path | str) -> dict:
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)

def resolve_case_key(doc: dict, case_id: str) -> str:
    """接受工况编号（``\"1\"``）或 ``case_id``。"""
    cases = doc["cases"]
    if case_id in cases:
        return str(case_id)
    for key, case in cases.items():
        if case.get("case_id") == case_id:
            return str(key)
    raise KeyError(f"找不到工况 {case_id}")

def read_field(doc: dict, case_id: str, name: str):
    """读取一个字段并换成 SI。``status=missing`` 时报错。"""
    case_key = resolve_case_key(doc, case_id)
    return _to_si(doc, case_key, name, _field_object(doc, case_key, name))

def read_case_field(path: Path | str, case_id: str, name: str):
    return read_field(load_case_file(path), case_id, name)

def assemble_case(path: Path | str, case_id: str) -> CaseConfig:
    """按 ``PHYSICAL_FIELDS`` 装配。任一字段 missing、未批准或对不上，整案失败。"""
    doc = load_case_file(path)
    case_key = resolve_case_key(doc, case_id)
    case = doc["cases"][case_key]
    values = {name: read_field(doc, case_key, name) for name in PHYSICAL_FIELDS}
    return CaseConfig(
        plant_id=str(doc["plant_id"]),
        case_key=case_key,
        case_id=str(case["case_id"]),
        physical=PhysicalConfig(**values),
        solver=SolverConfig(),
    )

def _field_object(doc: dict, case_key: str, name: str) -> dict:
    case = doc["cases"][case_key]
    fuel_id = case.get("fuel_id")
    fuel = doc.get("fuels", {}).get(fuel_id, {}) if fuel_id else {}
    for container in (case, fuel, doc.get("shared", {})):
        obj = container.get(name) if isinstance(container, dict) else None
        if isinstance(obj, dict) and "status" in obj:
            return obj
    raise FieldNotInPlant(name, "工况文件中没有这个字段")

def _to_si(doc: dict, case_key: str, name: str, obj: dict):
    status = obj.get("status")
    if name in _RECORD_ONLY and status == "cited":
        raise RecordOnlyFieldError(name, "status=cited，ADR 0005 只作记录，不作模型输入")
    if status == "missing" or obj.get("value") is None:
        raise MissingFieldError(name, "status=missing，原文未给出，不能填默认值")
    if status == "cited":
        _require_cited(name, case_key)
    elif status not in {"source", "derived"}:
        raise CaseAssemblyError(name, f"无法识别的 status={status}")

    value = obj["value"]
    unit = obj.get("unit")
    if name == "nitrogen_purge_volume_flow":
        _require_unit(name, unit, _NORMAL_VOLUME)
        return normal_volume_to_mol_s(value)
    if unit == _NORMAL_VOLUME:
        return normal_volume_to_mol_s(value)
    if name in _AMOUNT_FROM_VOLUME:
        _require_unit(name, unit, "mol/s")
        volume_name = _AMOUNT_FROM_VOLUME[name]
        converted = read_field(doc, case_key, volume_name)
        if not math.isclose(converted, float(value), rel_tol=1e-12, abs_tol=1e-12):
            raise CaseAssemblyError(
                name,
                f"按 ADR 0004 由体积流量换算得 {converted} mol/s，与文件中的 {value} 不一致",
            )
        return float(value)
    if unit not in _SI_UNITS:
        raise CaseAssemblyError(name, f"单位 {unit} 不是已登记的 SI 单位，不能静默换算")
    if name == "recirculation":
        if not isinstance(value, bool):
            raise CaseAssemblyError(name, "recirculation 必须是布尔值")
        return value
    if name in _INT_FIELDS:
        if isinstance(value, bool) or not isinstance(value, int):
            raise CaseAssemblyError(name, "计数必须是整数")
        return value
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    return float(value)

def _require_cited(name: str, case_key: str) -> None:
    if name == "nitrogen_purge_volume_flow":
        if case_key not in _NITROGEN_CASES:
            raise CitedFieldNotApproved(
                name,
                f"status=cited，ADR 0006 未批准工况 {case_key} 计入模型",
            )
        return
    allowed_cases = _CITED_CASES.get(name)
    if name in _CITED_MODEL and (allowed_cases is None or case_key in allowed_cases):
        return
    if allowed_cases is not None and case_key not in allowed_cases:
        raise CitedFieldNotApproved(
            name,
            f"status=cited，ADR 0005 未批准工况 {case_key} 使用该字段",
        )
    raise CitedFieldNotApproved(name, "status=cited，没有 ADR 批准进入模型")

def _require_unit(name: str, unit: object, expected: str) -> None:
    if unit != expected:
        raise CaseAssemblyError(name, f"单位是 {unit}，按登记应为 {expected}")
