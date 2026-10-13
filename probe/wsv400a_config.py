"""把 CASE_WSV400_A 与 _WSV400_SHARED 映到现有 ReactorConfig。

只写已有字段。论文模式开关保持类默认值；三组设置只动初值或离散预算。
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from src.core.feed_inlet import X_O2_AIR
from src.core.reactor import ReactorConfig

_ROOT = Path(__file__).resolve().parents[1]
_CASES = _ROOT / "data" / "validation_cases.json"
_PARITY = _ROOT / "data" / "hamel_fig7_parity_digitized.json"

# 与 validation JSON 的 Nm³ 口径一致：22.414 Nm³/kmol。
_NM3_PER_MOL = 22.414e-3

CASE_ID = "CASE_WSV400_A"
SHARED_ID = "_WSV400_SHARED"

# 三组上限。键必须是 ReactorConfig 已有字段。
SETTINGS: dict[str, dict[str, Any]] = {
    "thesis_default": {},
    "freeboard_n4": {"n_freeboard_cells": 4},
    "seed_budget_progress": {
        "vorab_nr_temperature_seed_mode_thesis": "budget_progress",
    },
}

SETTING_NOTES: dict[str, str] = {
    "thesis_default": (
        "thesis_mode=True，求解开关与温度种子保持 ReactorConfig 类默认值；"
        "max_global_iter=20，即 Reactor.solve 的默认参数。"
        "n_freeboard_cells 保持 0，自由空间高度写入但不进入方程图。"
    ),
    "freeboard_n4": (
        "与 thesis_default 相同，只把 n_freeboard_cells 从默认 0 改为 4。"
        "1.8 m 自由空间分成 4 格（每格 0.45 m）。格子数不在 JSON 里，"
        "这是离散预算，不是按出口调出来的。"
    ),
    "seed_budget_progress": (
        "与 thesis_default 相同，只把 vorab_nr_temperature_seed_mode_thesis "
        "从默认 adiabatic_anchor 改为已有模式 budget_progress。"
        "种子温度数值仍用类默认，不用出口温度。"
    ),
}


def load_validation() -> dict[str, Any]:
    return json.loads(_CASES.read_text(encoding="utf-8"))


def load_parity_wsv400_a() -> dict[str, Any]:
    parity = json.loads(_PARITY.read_text(encoding="utf-8"))
    return parity["fig_7_17_wsv400"]["cases"]["A"]


def _record(
    *,
    config_field: str | None,
    value: Any,
    source: str,
    status: str,
    note: str,
) -> dict[str, Any]:
    return {
        "config_field": config_field,
        "value": value,
        "source": source,
        "status": status,
        "note": note,
    }


def base_mapping_records(data: dict[str, Any]) -> list[dict[str, Any]]:
    """列出写入、沿用默认、以及无法写入的量。不改数值。"""
    case = data[CASE_ID]
    shared = data[SHARED_ID]
    fuel = case["inputs"]["fuel"]
    prox = fuel["proximate_analysis"]
    ult = fuel["ultimate_analysis_dry_wt_pct"]
    agent = case["inputs"]["gasification_agent"]
    reactor = case["inputs"]["reactor"]
    op = case["inputs"]["operating_conditions"]
    shared_reactor = shared["reactor"]
    feed_kg_s = float(fuel["feed_rate_kg_h"]) / 3600.0
    h_bed = float(reactor["bed_height_m"])
    h_total = float(reactor["total_height_m"])
    records = [
        _record(
            config_field="thesis_mode",
            value=True,
            source="任务卡：现有论文模式",
            status="task",
            note="打开后由 ReactorConfig.__post_init__ 套论文模式默认，本探针不再改那些开关。",
        ),
        _record(
            config_field="fuel_feed",
            value=feed_kg_s,
            source="CASE_WSV400_A.inputs.fuel.feed_rate_kg_h = 35.2",
            status="tabulated",
            note="35.2 kg/h ÷ 3600 → kg/s。",
        ),
        _record(
            config_field="moisture_wt",
            value=float(prox["moisture_wt_pct"]),
            source="CASE_WSV400_A.inputs.fuel.proximate_analysis.moisture_wt_pct",
            status="tabulated",
            note="湿基 wt%。",
        ),
        _record(
            config_field="ash_dry_wt",
            value=float(prox["ash_dry_wt_pct"]),
            source="CASE_WSV400_A.inputs.fuel.proximate_analysis.ash_dry_wt_pct",
            status="tabulated",
            note="干基 wt%。",
        ),
        _record(
            config_field="VM_daf",
            value=float(prox["VM_daf_pct"]),
            source="CASE_WSV400_A.inputs.fuel.proximate_analysis.VM_daf_pct",
            status="tabulated",
            note="无水无灰基 wt%。",
        ),
        _record(
            config_field="C_dry",
            value=float(ult["C"]),
            source="CASE_WSV400_A.inputs.fuel.ultimate_analysis_dry_wt_pct.C",
            status="tabulated",
            note="与 _WSV400_SHARED.fuel_properties_case_A_to_E 的 C 相同。",
        ),
        _record(
            config_field="H_dry",
            value=float(ult["H"]),
            source="CASE_WSV400_A.inputs.fuel.ultimate_analysis_dry_wt_pct.H",
            status="tabulated",
            note="干基 wt%。",
        ),
        _record(
            config_field="O_dry",
            value=float(ult["O"]),
            source="CASE_WSV400_A.inputs.fuel.ultimate_analysis_dry_wt_pct.O",
            status="tabulated",
            note="干基 wt%。",
        ),
        _record(
            config_field="S_dry",
            value=float(ult["S"]),
            source="CASE_WSV400_A.inputs.fuel.ultimate_analysis_dry_wt_pct.S",
            status="tabulated",
            note="只进入 ER 化学计量。Reactor 构造不把 S_dry 抄进 SolidProps.sulfur_fraction。",
        ),
        _record(
            config_field="nitrogen_fraction",
            value=float(ult["N"]),
            source="CASE_WSV400_A.inputs.fuel.ultimate_analysis_dry_wt_pct.N",
            status="tabulated",
            note="干基 wt%。elemental_ledger 里按 /100 使用，与 C_dry 同一百分数口径。",
        ),
        _record(
            config_field="fuel_type",
            value="biomass",
            source="CASE_WSV400_A.inputs.fuel.type = Natural Wood (Holz)",
            status="inferred",
            note="TarFuelType 只有 coal 与 biomass，没有 wood。木材取 biomass。",
        ),
        _record(
            config_field="D_bed",
            value=float(reactor["bed_diameter_m"]),
            source="CASE_WSV400_A.inputs.reactor.bed_diameter_m；_WSV400_SHARED.reactor.bed_diameter_m",
            status="tabulated",
            note="两处都是 0.4 m。",
        ),
        _record(
            config_field="H_bed",
            value=h_bed,
            source="CASE_WSV400_A.inputs.reactor.bed_height_m",
            status="inferred",
            note=(
                "JSON 自带 bed_height_source：由图 7.16 的 ξ≈0.4 推成约 40%×3 m。"
                "_WSV400_SHARED.reactor.bed_height_m 仍是 __MISSING__。"
            ),
        ),
        _record(
            config_field="H_freeboard",
            value=h_total - h_bed,
            source="total_height_m 3.0 − bed_height_m 1.2",
            status="inferred",
            note="总高在工况与共享块里都有。床高本身是图上推断，自由空间高度跟着减出来。",
        ),
        _record(
            config_field="P",
            value=float(op["pressure_MPa"]) * 1.0e6,
            source="CASE_WSV400_A.inputs.operating_conditions.pressure_MPa = 0.1",
            status="tabulated",
            note="0.1 MPa = 1.0e5 Pa。",
        ),
        _record(
            config_field="T_inlet",
            value=float(op["T_inlet_K"]),
            source="CASE_WSV400_A.inputs.operating_conditions.T_inlet_K",
            status="tabulated",
            note="293 K，不是 293.15。",
        ),
        _record(
            config_field="ER",
            value=float(op["ER"]),
            source="CASE_WSV400_A.inputs.operating_conditions.ER",
            status="tabulated",
            note="0.4532。ReactorConfig.__post_init__ 会据此重写 O2/H2O/N2。",
        ),
        _record(
            config_field="primary_agent",
            value="air",
            source="CASE_WSV400_A.inputs.gasification_agent.type = Air (only, no steam)",
            status="tabulated",
            note="已有枚举 air。",
        ),
        _record(
            config_field="steam_to_o2_molar",
            value=0.0,
            source="CASE_WSV400_A.inputs.gasification_agent.steam_feed_kg_h = 0",
            status="tabulated",
            note="类默认 0.8 会在纯空气工况里加入蒸汽，故显式写 0。",
        ),
        _record(
            config_field="recirculation_frac",
            value=0.0,
            source="CASE_WSV400_A.inputs.operating_conditions.recirculation = false",
            status="tabulated",
            note="类默认 0.1。工况写明不循环。",
        ),
        _record(
            config_field="recycle_gas",
            value=False,
            source="CASE_WSV400_A.inputs.operating_conditions.recirculation = false",
            status="tabulated",
            note="类默认 True。与 recirculation false 对齐。",
        ),
        _record(
            config_field="n_cells",
            value=10,
            source="ReactorConfig 类默认",
            status="class_default",
            note="JSON 没有轴向格子数。",
        ),
        _record(
            config_field="n_freeboard_cells",
            value=0,
            source="ReactorConfig 类默认",
            status="class_default",
            note="默认不求解自由空间。freeboard_n4 那一组才改成 4。",
        ),
        _record(
            config_field="d_p",
            value=0.00225,
            source="ReactorConfig 类默认",
            status="class_default",
            note=(
                "燃料给出的是等效粒径范围 2–26 mm、主导档 14–17 mm（质量分数 0.22）。"
                "没有单一粒径。床料 __MISSING__。默认分布模式 coarsest_class 会把质量放到最粗档，"
                "所以不把范围写进 d_p_min/d_p_max。"
            ),
        ),
        _record(
            config_field="rho_s",
            value=1400.0,
            source="ReactorConfig 类默认",
            status="class_default",
            note="_WSV400_SHARED.reactor.bed_material = __MISSING__，不另设密度。",
        ),
        _record(
            config_field="phi_s",
            value=0.86,
            source="ReactorConfig 类默认",
            status="class_default",
            note="JSON 无形状系数。",
        ),
        _record(
            config_field="eps_mf",
            value=0.45,
            source="ReactorConfig 类默认",
            status="class_default",
            note="JSON 无最小流化空隙率。",
        ),
        _record(
            config_field="HHV_dry",
            value=22.0,
            source="ReactorConfig 类默认",
            status="class_default",
            note="工况输入没有热值。hamel_plant_engineering.thermal_input_kW = 100 未换算成 HHV。",
        ),
        _record(
            config_field="heat_loss_frac",
            value=0.1,
            source="ReactorConfig 类默认",
            status="class_default",
            note="JSON 无热损失分数。不按出口温度调整。",
        ),
        _record(
            config_field=None,
            value=shared_reactor.get("bed_material"),
            source="_WSV400_SHARED.reactor.bed_material",
            status="unmapped",
            note="值是 __MISSING__，且没有对应的床料字段可填。",
        ),
        _record(
            config_field=None,
            value=reactor.get("freeboard_diameter_m"),
            source="CASE_WSV400_A.inputs.reactor.freeboard_diameter_m；共享块同值 0.6 m",
            status="unmapped",
            note="ReactorConfig 只有 D_bed。自由空间计算沿用床径。不新增字段。",
        ),
        _record(
            config_field=None,
            value=reactor.get("fuel_inlet_mm_above_nozzle"),
            source="CASE_WSV400_A.inputs.reactor.fuel_inlet_mm_above_nozzle",
            status="unmapped",
            note="喷嘴上方 200 mm。没有轴向进料位置字段，不改 solid_lower_inlet_frac。",
        ),
        _record(
            config_field=None,
            value=fuel["particle_diameter_mm"],
            source="CASE_WSV400_A.inputs.fuel.particle_diameter_mm",
            status="unmapped",
            note="范围与主导档都在，单一 d_p 的含义不明，保持类默认。",
        ),
        _record(
            config_field="freeboard_secondary_O2_mol_s",
            value=0.0,
            source="CASE_WSV400_A.inputs.gasification_agent.secondary_air_Nm3_h = 0",
            status="tabulated",
            note="类默认已是 0，与二次风为 0 一致，不另改。",
        ),
        _record(
            config_field=None,
            value=agent.get("primary_air_Nm3_h"),
            source="CASE_WSV400_A.inputs.gasification_agent.primary_air_Nm3_h",
            status="crosscheck",
            note="进料走 ER，不把 50 Nm³/h 再写成一套摩尔流。构造后与计算值对照。",
        ),
    ]
    return records


def build_reactor_config(setting: str) -> tuple[ReactorConfig, dict[str, Any]]:
    if setting not in SETTINGS:
        known = ", ".join(SETTINGS)
        raise KeyError(f"未知设置 {setting!r}，可选 {known}")
    data = load_validation()
    records = base_mapping_records(data)
    applied = {
        rec["config_field"]: rec["value"]
        for rec in records
        if rec["config_field"] is not None and rec["status"] != "crosscheck"
    }
    # 二次风默认已是 0，不重复写入。
    applied.pop("freeboard_secondary_O2_mol_s", None)
    overrides = dict(SETTINGS[setting])
    applied.update(overrides)
    cfg = ReactorConfig(**applied)
    for key, value in overrides.items():
        setattr(cfg, key, value)
    o2_nm3_h = float(cfg.O2_feed) * 3600.0 * _NM3_PER_MOL
    air_nm3_h = o2_nm3_h / X_O2_AIR
    meta = {
        "case_id": CASE_ID,
        "shared_id": SHARED_ID,
        "setting": setting,
        "setting_note": SETTING_NOTES[setting],
        "setting_overrides": overrides,
        "mapping": records,
        "computed_feeds_mol_s": {
            "O2_feed": float(cfg.O2_feed),
            "H2O_feed": float(cfg.H2O_feed),
            "N2_feed": float(cfg.N2_feed),
        },
        "computed_air_Nm3_h": {
            "O2": o2_nm3_h,
            "air_from_O2_over_0_21": air_nm3_h,
            "json_primary_air_Nm3_h": 50.0,
            "json_total_air_Nm3_h": 50.0,
            "json_calculated_air_actual_Nm3_h": 50.3,
            "nm3_per_kmol": 22.414,
        },
        "solve_kwargs": {
            "max_global_iter": 20,
            "tol_global": 1.0e-4,
            "solver": "global_nr",
        },
    }
    return cfg, meta
