#!/usr/bin/env python3
"""Build 15-case master summary table (JSON + Markdown) from validation_cases.json.

Usage:
    python3 scripts/build_validation_cases_master_table.py
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
CASES_PATH = ROOT / "data" / "validation_cases.json"
LITERATURE_PATH = ROOT / "data" / "validation_plant_literature.json"
JSON_OUT = ROOT / "data" / "validation_cases_master_table.json"
MD_OUT = ROOT / "docs/validation_cases_master_table.md"

CASE_ORDER: tuple[str, ...] = (
    "CASE_HTW_WESSELING_1",
    "CASE_HTW_WESSELING_2",
    "CASE_HTW_BERRENRATH",
    "CASE_VTT_PRESSURIZED_COAL_1",
    "CASE_VTT_PRESSURIZED_PEAT_6",
    "CASE_VTT_PRESSURIZED_PEAT_13",
    "CASE_VTT_PRESSURIZED_SAWDUST_14",
    "CASE_VTT_ATMOSPHERIC_PEAT_D",
    "CASE_VTT_ATMOSPHERIC_PEAT_E",
    "CASE_WSV400_A",
    "CASE_WSV400_B",
    "CASE_WSV400_C",
    "CASE_WSV400_D",
    "CASE_WSV400_E",
    "CASE_WSV400_F",
)

PLANT_ORDER: tuple[str, ...] = (
    "PLANT_HTW_WESSELING",
    "PLANT_HTW_BERRENRATH",
    "PLANT_VTT",
    "PLANT_WSV400",
)


def _load_plants(root: dict) -> dict[str, dict[str, Any]]:
    plants = root.get("validation_plants", {})
    return {k: v for k, v in plants.items() if not k.startswith("_") and isinstance(v, dict)}


def _case_to_plant_id(plants: dict[str, dict[str, Any]]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for plant_id, plant in plants.items():
        for case_id in plant.get("case_ids", []):
            mapping[str(case_id)] = plant_id
    return mapping


def _plant_summary_for_case(
    case_id: str,
    node: dict[str, Any],
    plants: dict[str, dict[str, Any]],
    case_to_plant: dict[str, str],
) -> dict[str, Any]:
    plant_id = node.get("plant_group_id") or case_to_plant.get(case_id)
    if not plant_id or plant_id not in plants:
        return {"plant_group_id": plant_id}
    plant = plants[plant_id]
    return {
        "plant_group_id": plant_id,
        "plant_group_name_zh": plant.get("display_name_zh"),
        "plant_group_name": plant.get("display_name"),
        "validation_tier": plant.get("validation_tier"),
        "validation_tier_zh": plant.get("validation_tier_zh"),
        "location": plant.get("location"),
        "operation_period": plant.get("operation_period"),
        "hamel_section": plant.get("hamel_section"),
        "data_provenance": plant.get("data_provenance"),
        "validation_grouping": plant.get("validation_grouping"),
    }


def _is_missing(val: object) -> bool:
    if val is None:
        return True
    if isinstance(val, str) and val.startswith("__"):
        return True
    return False


def _fmt(val: object, *, digits: int = 3) -> str | None:
    if _is_missing(val):
        return None
    if isinstance(val, bool):
        return "是" if val else "否"
    if isinstance(val, (int, float)):
        if isinstance(val, float) and abs(val - round(val)) > 1e-9:
            return f"{val:.{digits}f}".rstrip("0").rstrip(".")
        return str(int(val)) if isinstance(val, float) and val == int(val) else str(val)
    return str(val)


def _resolve_reactor(reactor: dict, root: dict) -> dict:
    if not isinstance(reactor, dict):
        return {}
    ref = reactor.get("_ref")
    base: dict = {}
    if ref and isinstance(root.get(ref), dict):
        shared = root[ref]
        if isinstance(shared.get("reactor"), dict):
            base = dict(shared["reactor"])
    merged = dict(base)
    merged.update({k: v for k, v in reactor.items() if not str(k).startswith("_") and k != "_ref"})
    return merged


def _gas_frac(out: dict, species: str) -> tuple[str | None, str | None]:
    """Return (value, basis) where basis is 'wet', 'dry', or None."""
    wet = out.get("exit_gas_wet_mol_frac")
    dry = out.get("exit_gas_dry_mol_frac")
    if species == "H2O" and isinstance(wet, dict):
        v = _fmt(wet.get("H2O"))
        return v, "wet" if v is not None else (None, None)
    if isinstance(wet, dict) and not _is_missing(wet.get(species)):
        v = _fmt(wet.get(species))
        if v is not None:
            return v, "wet"
    if isinstance(dry, dict) and not _is_missing(dry.get(species)):
        v = _fmt(dry.get(species))
        if v is not None:
            return v, "dry"
    return None, None


def _build_case_record(
    case_id: str,
    node: dict,
    root: dict,
    plants: dict[str, dict[str, Any]],
    case_to_plant: dict[str, str],
) -> dict[str, Any]:
    ins = node.get("inputs", {})
    fuel = ins.get("fuel", {})
    ga = ins.get("gasification_agent", {})
    react = _resolve_reactor(ins.get("reactor", {}), root)
    op = ins.get("operating_conditions", {})
    out = node.get("outputs", {})
    pa = fuel.get("proximate_analysis", {})
    ua = fuel.get("ultimate_analysis_dry_wt_pct", {})
    pd = fuel.get("particle_diameter_mm", {})

    wet_basis = isinstance(out.get("exit_gas_wet_mol_frac"), dict)
    dry_basis = isinstance(out.get("exit_gas_dry_mol_frac"), dict)

    co, co_basis = _gas_frac(out, "CO")
    co2, _ = _gas_frac(out, "CO2")
    h2, _ = _gas_frac(out, "H2")
    ch4, _ = _gas_frac(out, "CH4")
    h2o, _ = _gas_frac(out, "H2O")
    gas_basis = co_basis or ("wet" if wet_basis else ("dry" if dry_basis else None))

    record: dict[str, Any] = {
        "case_id": case_id,
        "hamel_simulation_id": node.get("hamel_simulation_id"),
        "hamel_table_ref": node.get("hamel_table_ref"),
        "comment": node.get("_comment"),
        "plant": node.get("plant"),
        "plant_group": _plant_summary_for_case(case_id, node, plants, case_to_plant),
        "country": node.get("country"),
        "scale": node.get("scale"),
        "source_input": node.get("_source_input"),
        "source_output": node.get("_source_output"),
        "reactor": {
            "bed_diameter_m": _fmt(react.get("bed_diameter_m") or react.get("diameter_m")),
            "freeboard_diameter_m": _fmt(react.get("freeboard_diameter_m") or react.get("diameter_m")),
            "total_height_m": _fmt(react.get("total_height_m") or react.get("height_m")),
            "bed_height_m": _fmt(react.get("bed_height_m")),
            "expanded_bed_height_m": _fmt(react.get("expanded_bed_height_m")),
            "reactor_height_to_cyclone_m": _fmt(react.get("reactor_height_to_cyclone_m")),
            "secondary_air_port_mm": _fmt(react.get("secondary_air_port_mm_above_nozzle")),
            "fuel_inlet_mm_above_nozzle": _fmt(react.get("fuel_inlet_mm_above_nozzle")),
            "bed_material": None if _is_missing(react.get("bed_material")) else react.get("bed_material"),
            "bed_particle_diameter_mm": None if _is_missing(react.get("bed_particle_diameter_mm")) else react.get("bed_particle_diameter_mm"),
        },
        "fuel": {
            "type": fuel.get("type"),
            "feed_rate_kg_h": _fmt(fuel.get("feed_rate_kg_h")),
            "proximate_analysis": {
                "moisture_wt_pct": _fmt(pa.get("moisture_wt_pct")) if isinstance(pa, dict) else None,
                "ash_dry_wt_pct": _fmt(pa.get("ash_dry_wt_pct")) if isinstance(pa, dict) else None,
                "VM_daf_pct": _fmt(pa.get("VM_daf_pct")) if isinstance(pa, dict) else None,
            },
            "ultimate_analysis_dry_wt_pct": {
                "C": _fmt(ua.get("C")) if isinstance(ua, dict) else None,
                "H": _fmt(ua.get("H")) if isinstance(ua, dict) else None,
                "O": _fmt(ua.get("O")) if isinstance(ua, dict) else None,
                "N": _fmt(ua.get("N")) if isinstance(ua, dict) else None,
                "S": _fmt(ua.get("S")) if isinstance(ua, dict) else None,
            },
            "particle_diameter_mm": pd if isinstance(pd, dict) else None,
        },
        "inlets": {
            "gasification_agent_type": ga.get("type"),
            "pressure_MPa": _fmt(op.get("pressure_MPa")),
            "ER_lambda_O2": _fmt(op.get("ER")),
            "T_inlet_K": _fmt(op.get("T_inlet_K")),
            "recirculation": op.get("recirculation") if not _is_missing(op.get("recirculation")) else None,
            "primary_air_Nm3_h": _fmt(ga.get("primary_air_Nm3_h")),
            "secondary_air_Nm3_h": _fmt(ga.get("secondary_air_Nm3_h")),
            "total_air_Nm3_h": _fmt(ga.get("total_air_Nm3_h")),
            "primary_O2_Nm3_h": _fmt(ga.get("primary_O2_Nm3_h")),
            "secondary_O2_Nm3_h": _fmt(ga.get("secondary_O2_Nm3_h")),
            "total_O2_Nm3_h": _fmt(ga.get("total_O2_Nm3_h")),
            "steam_flow_kg_h": _fmt(ga.get("steam_feed_kg_h")),
            "secondary_injection_xi": _fmt(ga.get("secondary_injection_xi_reactor")),
            "secondary_injection_agent": ga.get("secondary_injection_agent"),
        },
        "outlets": {
            "gas_basis": gas_basis,
            "exit_temperature_K": _fmt(out.get("exit_temperature_K")),
            "exit_temperature_note": out.get("exit_temperature_note"),
            "carbon_conversion_pct": _fmt(out.get("carbon_conversion_pct")),
            "CO": co,
            "CO2": co2,
            "H2": h2,
            "CH4": ch4,
            "H2O": h2o,
            "NH3_ppm": _fmt(out.get("nitrogen_compounds", {}).get("NH3_ppm"))
            if isinstance(out.get("nitrogen_compounds"), dict)
            else None,
            "HCN_ppm": _fmt(out.get("nitrogen_compounds", {}).get("HCN_ppm"))
            if isinstance(out.get("nitrogen_compounds"), dict)
            else None,
            "parity_plot_available": bool(out.get("parity_plot_data", {}).get("available")),
            "axial_profile_available": bool(out.get("axial_profiles", {}).get("available")),
            "output_precision": out.get("_precision"),
        },
    }

    # strip None-only nested empties for cleaner JSON
    for section in ("reactor", "fuel", "inlets", "outlets"):
        record[section] = {k: v for k, v in record[section].items() if v is not None}

    prox = record["fuel"].get("proximate_analysis", {})
    if prox and all(v is None for v in prox.values()):
        record["fuel"].pop("proximate_analysis", None)
    ult = record["fuel"].get("ultimate_analysis_dry_wt_pct", {})
    if ult and all(v is None for v in ult.values()):
        record["fuel"].pop("ultimate_analysis_dry_wt_pct", None)

    return record


def _md_cell(val: object) -> str:
    if val is None or val == "":
        return "—"
    text = str(val).replace("|", "\\|").replace("\n", " ")
    return text


def _numeric_values(values: list[object]) -> list[float]:
    result: list[float] = []
    for value in values:
        if value is None or isinstance(value, bool):
            continue
        try:
            result.append(float(value))
        except (TypeError, ValueError):
            continue
    return result


def _range_text(
    values: list[object],
    *,
    scale: float = 1.0,
    digits: int = 3,
    thousands: bool = False,
) -> str:
    numbers = [value * scale for value in _numeric_values(values)]
    if not numbers:
        return "—"

    def format_number(value: float) -> str:
        text = f"{value:,.{digits}f}" if thousands else f"{value:.{digits}f}"
        return text.rstrip("0").rstrip(".") if "." in text else text

    return f"{format_number(min(numbers))}–{format_number(max(numbers))}"


def _load_literature() -> dict[str, Any]:
    return json.loads(LITERATURE_PATH.read_text(encoding="utf-8"))


def _format_ref_bibliography(ref: dict[str, Any]) -> str:
    parts: list[str] = [ref.get("authors", "")]
    if ref.get("title"):
        parts.append(f"*{ref['title']}*")
    venue = ref.get("venue", "")
    if ref.get("year"):
        venue = f"{venue}, {ref['year']}" if venue else str(ref["year"])
    if venue:
        parts.append(venue)
    if ref.get("volume") and ref.get("pages"):
        parts.append(f"{ref['volume']}({ref.get('issue', '?')}): {ref['pages']}")
    elif ref.get("pages"):
        parts.append(f"pp. {ref['pages']}")
    if ref.get("isbn"):
        parts.append(f"ISBN {ref['isbn']}")
    return ". ".join(p for p in parts if p)


def _format_ref_links(ref: dict[str, Any]) -> str:
    links: list[str] = []
    if ref.get("doi"):
        links.append(f"[DOI](https://doi.org/{ref['doi']})")
    for url in ref.get("urls") or []:
        if url and (not ref.get("doi") or url != f"https://doi.org/{ref['doi']}"):
            links.append(f"[链接]({url})")
    if not links:
        note = ref.get("access_note") or "无公开链接"
        return note.split("；")[0].split(";")[0]
    return " · ".join(links)


def _render_literature_section(
    literature: dict[str, Any],
    plants: dict[str, dict[str, Any]],
) -> list[str]:
    refs = literature.get("refs", {})
    plant_map = literature.get("plant_literature_map", {})
    role_labels = {
        "primary_data": "实测数据（primary）",
        "secondary_compilation": "二次汇编（Hamel）",
        "engineering_background": "装置工程背景",
    }

    lines: list[str] = [
        "",
        "## 7. 原始文献（按装置分组）",
        "",
        f"> SSOT：[`data/validation_plant_literature.json`](../data/validation_plant_literature.json)",
        "",
        "下列为 Hamel Ch.7 各装置实测数据的**原始出处**；Hamel (2001) VDI 469 为 Table 7.1–7.4 的二次汇编。",
        "无公开 URL 的条目见「获取方式」。",
        "",
    ]

    for plant_id in PLANT_ORDER:
        if plant_id not in plant_map:
            continue
        mapping = plant_map[plant_id]
        plant = plants.get(plant_id, {})
        lines.append(f"### 7.{PLANT_ORDER.index(plant_id) + 1} `{plant_id}` — {plant.get('display_name_zh', plant_id)}")
        lines.append("")
        if mapping.get("data_provider_note"):
            lines.append(f"**数据提供**：{mapping['data_provider_note']}")
            lines.append("")
        if mapping.get("hamel_tables_figures"):
            lines.append(
                "**Hamel 引用**："
                + "；".join(str(x) for x in mapping["hamel_tables_figures"])
            )
            lines.append("")
        lines.extend(
            [
                "| 角色 | ref_id | 文献 | 链接 | 数据用途 |",
                "|------|--------|------|------|----------|",
            ]
        )
        for role_key, role_zh in role_labels.items():
            for ref_id in mapping.get(role_key, []):
                ref = refs.get(ref_id, {})
                lines.append(
                    f"| {role_zh} | `{ref_id}` | {_md_cell(_format_ref_bibliography(ref))} | "
                    f"{_md_cell(_format_ref_links(ref))} | {_md_cell(ref.get('data_use'))} |"
                )
        lines.append("")

    lines.extend(
        [
            "### 7.5 全文献书目索引",
            "",
            "| ref_id | 短引 | DOI | 链接 | 获取说明 |",
            "|--------|------|-----|------|----------|",
        ]
    )
    for ref_id, ref in refs.items():
        doi_cell = f"[{ref['doi']}](https://doi.org/{ref['doi']})" if ref.get("doi") else "—"
        url_cell = _format_ref_links(ref)
        lines.append(
            f"| `{ref_id}` | {_md_cell(ref.get('short_key'))} | {doi_cell} | "
            f"{_md_cell(url_cell)} | {_md_cell(ref.get('access_note'))} |"
        )
    lines.append("")
    return lines


def _render_markdown(
    cases: list[dict[str, Any]],
    plants: dict[str, dict[str, Any]],
    literature: dict[str, Any],
) -> str:
    lines: list[str] = [
        "# Validation Cases 主表（15 Case / 4 套装置）",
        "",
        f"> 自动生成自 [`data/validation_cases.json`](../data/validation_cases.json)，"
        f"脚本 [`scripts/build_validation_cases_master_table.py`](../scripts/build_validation_cases_master_table.py)。"
        f" 更新日期：{date.today().isoformat()}",
        "",
        "15 个 Hamel 验证工况对应 **4 套物理装置**：Wesseling 中试 2 组、Berrenrath 示范 1 组、"
        "VTT 实验台 6 组（常压/加压为同一装置的不同压力工况）、WSV 400 研究炉 6 组。",
        "",
        "单位：压力 MPa；流量 Nm³/h（STP）；燃料分析 wt%；出口摩尔分数为 **湿基** 或 **干基**（见「气基」列）。",
        "",
        "## 0. 装置分组（Validation 用）",
        "",
        "| 装置组 ID | 名称 | 规模 | Case 数 | 运行期 | Hamel | 数据出处 | 分组说明 |",
        "|-----------|------|------|---------|--------|-------|----------|----------|",
    ]

    for plant_id in PLANT_ORDER:
        if plant_id not in plants:
            continue
        p = plants[plant_id]
        lines.append(
            f"| `{plant_id}` | {_md_cell(p.get('display_name_zh'))} | "
            f"{_md_cell(p.get('validation_tier_zh'))} | {_md_cell(p.get('case_count'))} | "
            f"{_md_cell(p.get('operation_period'))} | {_md_cell(p.get('hamel_section'))} | "
            f"{_md_cell(p.get('data_provenance'))} | {_md_cell(p.get('validation_grouping'))} |"
        )

    lines.extend(
        [
            "",
            "### 0.1 Case → 装置组映射",
            "",
            "| Case ID | 装置组 | 规模 | 装置实例名 | Hamel Sim |",
            "|---------|--------|------|------------|-----------|",
        ]
    )
    for c in cases:
        pg = c.get("plant_group", {})
        lines.append(
            f"| {_md_cell(c['case_id'])} | `{_md_cell(pg.get('plant_group_id'))}` | "
            f"{_md_cell(pg.get('validation_tier_zh'))} | {_md_cell(c.get('plant'))} | "
            f"{_md_cell(c.get('hamel_simulation_id'))} |"
        )

    wet_cases = [c for c in cases if c.get("outlets", {}).get("gas_basis") == "wet"]
    dry_cases = [c for c in cases if c.get("outlets", {}).get("gas_basis") == "dry"]
    temperature_K = _numeric_values([c.get("outlets", {}).get("exit_temperature_K") for c in cases])
    temperature_C = [value - 273.15 for value in temperature_K]
    carbon_conversion = _numeric_values(
        [c.get("outlets", {}).get("carbon_conversion_pct") for c in cases]
    )

    lines.extend(
        [
            "",
            "### 0.2 跨工况实验覆盖范围",
            "",
            "以下范围由主表收录的工况汇总而成，不对应某一个单独 Case。出口组成来自论文图的数字化读数，主要用于说明验证窗口，不能替代原始测量表。",
            "",
            "| 指标 | 汇总范围 | 说明 |",
            "|------|----------|------|",
            f"| 反应器压力 | {_range_text([c.get('inlets', {}).get('pressure_MPa') for c in cases])} MPa | 从常压实验台到加压 HTW 装置 |",
            f"| 燃料给料量 | {_range_text([c.get('fuel', {}).get('feed_rate_kg_h') for c in cases], digits=1, thousands=True)} kg/h | 按各试验表给出的入炉燃料流量 |",
            f"| 化学计量氧比 | {_range_text([c.get('inlets', {}).get('ER_lambda_O2') for c in cases])} | 空气、空气／水蒸气和氧气／水蒸气工况 |",
            f"| 燃料含水量 | {_range_text([c.get('fuel', {}).get('proximate_analysis', {}).get('moisture_wt_pct') for c in cases], digits=1)} wt% | 褐煤、泥炭、锯末与木片 |",
            f"| 代表性出口／平均自由空间温度 | {_range_text(temperature_K, digits=0)} K（{_range_text(temperature_C, digits=0)} ℃） | 共 {len(temperature_K)} 组温度数据 |",
            f"| 总碳转化率 | {_range_text(carbon_conversion, digits=1)}% | 仅统计有实验值的 {len(carbon_conversion)} 组工况 |",
            "",
            "产品气组成必须按湿基和干基分开汇总：",
            "",
            "| 气基 | CO | CO₂ | H₂ | CH₄ | H₂O |",
            "|------|----|-----|----|-----|-----|",
            f"| 湿基（{len(wet_cases)} 组） | {_range_text([c['outlets'].get('CO') for c in wet_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('CO2') for c in wet_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('H2') for c in wet_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('CH4') for c in wet_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('H2O') for c in wet_cases], scale=100, digits=1)} mol% |",
            f"| 干基（{len(dry_cases)} 组） | {_range_text([c['outlets'].get('CO') for c in dry_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('CO2') for c in dry_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('H2') for c in dry_cases], scale=100, digits=1)} mol% | {_range_text([c['outlets'].get('CH4') for c in dry_cases], scale=100, digits=1)} mol% | — |",
            "",
            "湿基与干基浓度不能直接合并比较；干基组没有水蒸气浓度。各图读数的数字化误差约为 ±5%，范围端点不宜作高精度设计值使用。",
        ]
    )

    lines.extend(
        [
            "",
            "## 1. 总览表",
            "",
            "| Case ID | 装置组 | Hamel | 装置 | 燃料 | P [MPa] | ER | 循环 | T_out [K] | Xc [%] | CO | CO₂ | H₂ | CH₄ | H₂O | 气基 |",
            "|---------|--------|-------|------|------|---------|-----|------|-----------|--------|-----|-----|-----|-----|-----|------|",
        ]
    )

    for c in cases:
        inl = c["inlets"]
        ou = c["outlets"]
        pg = c.get("plant_group", {})
        lines.append(
            "| {case_id} | {group} | {hamel} | {plant} | {fuel} | {p} | {er} | {rec} | {t} | {xc} | {co} | {co2} | {h2} | {ch4} | {h2o} | {basis} |".format(
                case_id=_md_cell(c["case_id"]),
                group=_md_cell(pg.get("plant_group_id")),
                hamel=_md_cell(c.get("hamel_simulation_id")),
                plant=_md_cell(c.get("plant")),
                fuel=_md_cell(c["fuel"].get("type")),
                p=_md_cell(inl.get("pressure_MPa")),
                er=_md_cell(inl.get("ER_lambda_O2")),
                rec=_md_cell("是" if inl.get("recirculation") else "否" if inl.get("recirculation") is not None else None),
                t=_md_cell(ou.get("exit_temperature_K")),
                xc=_md_cell(ou.get("carbon_conversion_pct")),
                co=_md_cell(ou.get("CO")),
                co2=_md_cell(ou.get("CO2")),
                h2=_md_cell(ou.get("H2")),
                ch4=_md_cell(ou.get("CH4")),
                h2o=_md_cell(ou.get("H2O")),
                basis=_md_cell(ou.get("gas_basis")),
            )
        )

    lines.extend(
        [
            "",
            "## 2. 反应器尺寸",
            "",
            "| Case ID | D_bed [m] | D_fb [m] | H_total [m] | H_bed [m] | H_exp [m] | 二级风高度 [mm] | 进料高度 [mm] |",
            "|---------|-----------|----------|-------------|-----------|-----------|-----------------|---------------|",
        ]
    )
    for c in cases:
        r = c["reactor"]
        lines.append(
            f"| {_md_cell(c['case_id'])} | {_md_cell(r.get('bed_diameter_m'))} | "
            f"{_md_cell(r.get('freeboard_diameter_m'))} | {_md_cell(r.get('total_height_m'))} | "
            f"{_md_cell(r.get('bed_height_m'))} | {_md_cell(r.get('expanded_bed_height_m'))} | "
            f"{_md_cell(r.get('secondary_air_port_mm'))} | {_md_cell(r.get('fuel_inlet_mm_above_nozzle'))} |"
        )

    lines.extend(
        [
            "",
            "## 3. 进料与边界",
            "",
            "| Case ID | 燃料 [kg/h] | 氧化剂 | 主空气 | 次空气 | 总空气 | 主O₂ | 次O₂ | 总O₂ | 蒸汽 [kg/h] | ξ_sec | T_in [K] |",
            "|---------|-------------|--------|--------|--------|--------|------|------|------|-------------|-------|----------|",
        ]
    )
    for c in cases:
        inl = c["inlets"]
        lines.append(
            f"| {_md_cell(c['case_id'])} | {_md_cell(c['fuel'].get('feed_rate_kg_h'))} | "
            f"{_md_cell(inl.get('gasification_agent_type'))} | "
            f"{_md_cell(inl.get('primary_air_Nm3_h'))} | {_md_cell(inl.get('secondary_air_Nm3_h'))} | "
            f"{_md_cell(inl.get('total_air_Nm3_h'))} | {_md_cell(inl.get('primary_O2_Nm3_h'))} | "
            f"{_md_cell(inl.get('secondary_O2_Nm3_h'))} | {_md_cell(inl.get('total_O2_Nm3_h'))} | "
            f"{_md_cell(inl.get('steam_flow_kg_h'))} | {_md_cell(inl.get('secondary_injection_xi'))} | "
            f"{_md_cell(inl.get('T_inlet_K'))} |"
        )

    lines.extend(
        [
            "",
            "## 4. 燃料工业/元素分析",
            "",
            "| Case ID | M [%] | Ash_d [%] | VM_daf [%] | C_d [%] | H_d [%] | O_d [%] | N_d [%] | S_d [%] |",
            "|---------|-------|-----------|------------|---------|---------|---------|---------|---------|",
        ]
    )
    for c in cases:
        pa = c["fuel"].get("proximate_analysis", {})
        ua = c["fuel"].get("ultimate_analysis_dry_wt_pct", {})
        lines.append(
            f"| {_md_cell(c['case_id'])} | {_md_cell(pa.get('moisture_wt_pct'))} | "
            f"{_md_cell(pa.get('ash_dry_wt_pct'))} | {_md_cell(pa.get('VM_daf_pct'))} | "
            f"{_md_cell(ua.get('C'))} | {_md_cell(ua.get('H'))} | {_md_cell(ua.get('O'))} | "
            f"{_md_cell(ua.get('N'))} | {_md_cell(ua.get('S'))} |"
        )

    lines.extend(
        [
            "",
            "## 5. 数据来源",
            "",
            "| Case ID | 输入来源 | 输出来源 | 轴向剖面 | Parity |",
            "|---------|----------|----------|----------|--------|",
        ]
    )
    for c in cases:
        ou = c["outlets"]
        lines.append(
            f"| {_md_cell(c['case_id'])} | {_md_cell(c.get('source_input'))} | "
            f"{_md_cell(c.get('source_output'))} | "
            f"{_md_cell('有' if ou.get('axial_profile_available') else '无')} | "
            f"{_md_cell('有' if ou.get('parity_plot_available') else '无')} |"
        )

    lines.extend(
        [
            "",
            "## 6. JSON SSOT",
            "",
            "结构化全字段副本：[`data/validation_cases_master_table.json`](../data/validation_cases_master_table.json)",
            "",
            "原始嵌套 SSOT：[`data/validation_cases.json`](../data/validation_cases.json)",
            "",
            "装置原始文献 SSOT：[`data/validation_plant_literature.json`](../data/validation_plant_literature.json)",
            "",
        ]
    )
    lines.extend(_render_literature_section(literature, plants))
    return "\n".join(lines)


def main() -> None:
    root = json.loads(CASES_PATH.read_text(encoding="utf-8"))
    literature = _load_literature()
    plants = _load_plants(root)
    case_to_plant = _case_to_plant_id(plants)
    cases: list[dict[str, Any]] = []
    for case_id in CASE_ORDER:
        if case_id not in root:
            raise KeyError(f"missing {case_id} in validation_cases.json")
        cases.append(_build_case_record(case_id, root[case_id], root, plants, case_to_plant))

    plant_payload = {pid: plants[pid] for pid in PLANT_ORDER if pid in plants}

    payload = {
        "_meta": {
            "description": "15 validation cases master summary — flattened from validation_cases.json",
            "source": str(CASES_PATH.relative_to(ROOT)),
            "literature_source": str(LITERATURE_PATH.relative_to(ROOT)),
            "generated_by": "scripts/build_validation_cases_master_table.py",
            "generated_on": date.today().isoformat(),
            "case_count": len(cases),
            "plant_group_count": len(plant_payload),
            "literature_ref_count": len(literature.get("refs", {})),
            "notes": [
                "15 case 对应 4 套物理装置；见 validation_plants 与每 case 的 plant_group",
                "出口气组成优先取 validation_cases.json outputs；湿基/干基见 outlets.gas_basis",
                "缺失字段在 JSON 中为 null；MD 表中以 — 表示",
                "原始文献书目与链接见 validation_literature 与 MD §7",
            ],
        },
        "validation_plants": plant_payload,
        "validation_literature": {
            "_meta": literature.get("_meta"),
            "refs": literature.get("refs"),
            "plant_literature_map": literature.get("plant_literature_map"),
        },
        "cases": cases,
    }

    JSON_OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD_OUT.write_text(_render_markdown(cases, plant_payload, literature), encoding="utf-8")
    print(f"Wrote {JSON_OUT} ({len(cases)} cases)")
    print(f"Wrote {MD_OUT}")


if __name__ == "__main__":
    main()
