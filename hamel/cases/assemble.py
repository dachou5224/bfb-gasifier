"""通用地装配 data/cases/*.json。missing 报错；cited 只在 ADR 批准的装置内可用。"""
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
# ADR 0006：仅 VTT 加压工况 1、2 把吹扫氮计入模型。
_NITROGEN_CASES = {"1", "2"}
_VTT_PLANT = "VTT_PRESSURIZED"
# ADR 0005：这两项分别只属于一种燃料。
_CITED_FUELS = {
    "particle_size_max": {"rhein_braunkohle"},
    "particle_size_range": {"finn_torf"},
}
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
_NORMAL_VOLUME = "m3(VN)/h"

def _field_spec() -> dict[str, tuple[str, str]]:
    """字段名 → (类型, 原文单位)。只有 volume 才做 ADR 0004 换算。"""
    groups = (
        ("float", "m", (
            "bed_diameter", "freeboard_diameter", "reactor_height", "bed_zone_height",
            "freeboard_height", "fuel_feed_height", "particle_size_max",
        )),
        ("float", "kg/kg", (
            "ash_mass_fraction_wf", "volatiles_mass_fraction_waf", "carbon_mass_fraction_wf",
            "hydrogen_mass_fraction_wf", "oxygen_mass_fraction_wf", "nitrogen_mass_fraction_wf",
            "sulfur_mass_fraction_wf", "fixed_carbon_mass_fraction_wf", "moisture_mass_fraction",
        )),
        ("float", "kg/s", ("fuel_mass_flow", "steam_mass_flow")),
        ("float", "mol/s", ("primary_air_amount_flow", "secondary_air_amount_flow")),
        ("float", "K", ("inlet_temperature",)),
        ("float", "Pa", ("pressure",)),
        ("float", "J/kg", ("heating_value",)),
        ("float", "1", ("lambda_o2",)),
        ("int", "1", ("cyclone_cell_count", "fuel_feed_count", "agent_feed_count")),
        ("bool", "1", ("recirculation",)),
        ("volume", _NORMAL_VOLUME, (
            "nitrogen_purge_volume_flow", "primary_air_volume_flow", "secondary_air_volume_flow",
        )),
        ("range", "m", ("particle_size_range",)),
    )
    spec: dict[str, tuple[str, str]] = {}
    for kind, unit, names in groups:
        for name in names:
            spec[name] = (kind, unit)
    return spec

_FIELD_SPEC = _field_spec()
# 体积流量原字段和粒径范围不进 PhysicalConfig，但读取时仍要核量纲。
_NOT_CONFIG = {"primary_air_volume_flow", "secondary_air_volume_flow", "particle_size_range"}
PHYSICAL_FIELDS = tuple(name for name in _FIELD_SPEC if name not in _NOT_CONFIG)

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
    if status == "missing" or obj.get("value") is None:
        raise MissingFieldError(name, "status=missing，原文未给出，不能填默认值")
    fuel_id = doc["cases"][case_key].get("fuel_id")
    if status == "cited":
        _require_cited(str(doc.get("plant_id")), name, case_key, fuel_id)
    elif status not in {"source", "derived"}:
        raise CaseAssemblyError(name, f"无法识别的 status={status}")
    if name not in _FIELD_SPEC:
        raise CaseAssemblyError(name, "没有登记量纲，不能装配")
    kind, expected = _FIELD_SPEC[name]
    _require_unit(name, obj.get("unit"), expected)
    value = obj["value"]
    if kind == "volume":
        return normal_volume_to_mol_s(_require_finite(name, value))
    if name in _AMOUNT_FROM_VOLUME:
        number = _require_finite(name, value)
        converted = read_field(doc, case_key, _AMOUNT_FROM_VOLUME[name])
        if not math.isclose(converted, number, rel_tol=1e-12, abs_tol=1e-12):
            raise CaseAssemblyError(
                name,
                f"按 ADR 0004 由体积流量换算得 {converted} mol/s，与文件中的 {number} 不一致",
            )
        return number
    if kind == "bool":
        if not isinstance(value, bool):
            raise CaseAssemblyError(name, "必须是布尔值")
        return value
    if kind == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            raise CaseAssemblyError(name, "必须是整数")
        return value
    if kind == "range":
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise CaseAssemblyError(name, "必须是两个有限数值")
        return tuple(_require_finite(name, item) for item in value)
    return _require_finite(name, value)

def _require_finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CaseAssemblyError(name, "必须是有限数值")
    number = float(value)
    if not math.isfinite(number):
        raise CaseAssemblyError(name, "必须是有限数值")
    return number

def _require_cited(plant_id: str, name: str, case_key: str, fuel_id: str | None) -> None:
    if plant_id != _VTT_PLANT:
        raise CitedFieldNotApproved(name, f"status=cited，ADR 未批准装置 {plant_id}")
    if name in _RECORD_ONLY:
        raise RecordOnlyFieldError(name, "status=cited，ADR 0005 只作记录，不作模型输入")
    fuels = _CITED_FUELS.get(name)
    if fuels is not None and fuel_id not in fuels:
        raise CitedFieldNotApproved(name, f"status=cited，ADR 0005 未批准燃料 {fuel_id}")
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
