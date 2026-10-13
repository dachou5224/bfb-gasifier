#!/usr/bin/env python3
"""Merge Hamel Table 7.1–7.4 + plant engineering into validation_cases.json.

Source SSOT: data/hamel_table7_extractions.json (from Hamel 1999 Ch.7 OCR/text).

Usage:
    python3 scripts/enrich_validation_cases_from_hamel.py
    python3 scripts/enrich_validation_cases_from_hamel.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
from copy import deepcopy
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
CASES_PATH = ROOT / "data" / "validation_cases.json"
HAMEL_PATH = ROOT / "data" / "hamel_table7_extractions.json"


def _is_missing(val: object) -> bool:
    if val is None:
        return True
    if isinstance(val, str) and val.startswith("__"):
        return True
    return False


def _deep_set_if_missing(target: dict, key: str, value: object) -> None:
    if value is None:
        return
    if key not in target or _is_missing(target.get(key)):
        target[key] = value


def _fuel_analysis_block(row: dict, *, coal: dict | None, peat: dict | None, sawdust: dict | None) -> dict:
    fuel_name = str(row.get("fuel", "")).lower()
    if "braunk" in fuel_name or "coal" in fuel_name:
        base = dict(coal or {})
    elif "torf" in fuel_name or "peat" in fuel_name:
        base = dict(peat or {})
    elif "säge" in fuel_name or "saw" in fuel_name or "mehl" in fuel_name:
        base = dict(sawdust or {})
    else:
        base = {}
    for k in ("moisture_wt_pct", "ash_dry_wt_pct", "VM_daf_pct", "C_dry", "H_dry", "O_dry", "N_dry", "S_dry"):
        if k in row and row[k] is not None:
            base[k] = row[k]
    return base


def _apply_table_row(case: dict, row: dict, *, plant_key: str, hamel: dict) -> None:
    inputs = case.setdefault("inputs", {})
    fuel = inputs.setdefault("fuel", {})
    ga = inputs.setdefault("gasification_agent", {})
    react = inputs.setdefault("reactor", {})
    op = inputs.setdefault("operating_conditions", {})

    eng = hamel["plant_engineering"].get(plant_key, {})
    case["hamel_table_ref"] = case.get("hamel_table_ref") or {}
    case["hamel_plant_engineering"] = eng

    # Fuel feed & analysis
    _deep_set_if_missing(fuel, "feed_rate_kg_h", row.get("fuel_feed_kg_h"))
    fa = row.get("fuel_analysis")
    if fa is None:
        fa = _fuel_analysis_block(
            row,
            coal=hamel.get("table_7_2_vtt_pressurized", {}).get("coal_fuel_analysis"),
            peat=hamel.get("table_7_2_vtt_pressurized", {}).get("peat_fuel_analysis")
            or hamel.get("table_7_3_vtt_atmospheric", {}).get("peat_fuel_analysis"),
            sawdust=hamel.get("table_7_2_vtt_pressurized", {}).get("sawdust_fuel_analysis"),
        )
    pa = fuel.get("proximate_analysis")
    if not isinstance(pa, dict):
        pa = {}
        fuel["proximate_analysis"] = pa
    ua = fuel.get("ultimate_analysis_dry_wt_pct")
    if not isinstance(ua, dict):
        ua = {}
        fuel["ultimate_analysis_dry_wt_pct"] = ua

    for k in ("moisture_wt_pct", "ash_dry_wt_pct", "VM_daf_pct"):
        if fa.get(k) is not None:
            _deep_set_if_missing(pa, k, fa[k])
    elem_map = {"C_dry": "C", "H_dry": "H", "O_dry": "O", "N_dry": "N", "S_dry": "S"}
    for src, dst_key in elem_map.items():
        val = fa.get(src)
        if val is not None:
            _deep_set_if_missing(ua, dst_key, val)
    # Clean erroneous duplicate keys from prior partial fills
    for bad in ("C_dry", "H_dry", "O_dry", "N_dry", "S_dry"):
        ua.pop(bad, None)

    if _is_missing(fuel.get("particle_diameter_mm")):
        fuel["particle_diameter_mm"] = {
            "mean_range": [1.0, 4.0],
            "_source": "Hamel Ch.7 VTT/Kurkela typical feed; exact sieve not tabulated in Table 7.2",
        }

    # Gasification agent
    if row.get("primary_air_Nm3_h") is not None:
        _deep_set_if_missing(ga, "primary_air_Nm3_h", row["primary_air_Nm3_h"])
        total = float(row["primary_air_Nm3_h"]) + float(row.get("secondary_air_Nm3_h") or 0.0)
        _deep_set_if_missing(ga, "total_air_Nm3_h", total)
    if row.get("secondary_air_Nm3_h") is not None:
        _deep_set_if_missing(ga, "secondary_air_Nm3_h", row["secondary_air_Nm3_h"])
    if row.get("steam_kg_h") is not None:
        _deep_set_if_missing(ga, "steam_feed_kg_h", row["steam_kg_h"])
    if row.get("steam_to_dry_fuel_ratio") is not None:
        ga["steam_to_dry_fuel_ratio"] = row["steam_to_dry_fuel_ratio"]

    # Operating
    if row.get("pressure_MPa") is not None:
        _deep_set_if_missing(op, "pressure_MPa", row["pressure_MPa"])
    if row.get("ER") is not None:
        _deep_set_if_missing(op, "ER", row["ER"])
    if row.get("recirculation") is not None:
        _deep_set_if_missing(op, "recirculation", bool(row["recirculation"]))
        if row["recirculation"]:
            _deep_set_if_missing(op, "recycle_topology", "cyclone return to bed")
            _deep_set_if_missing(op, "cyclone_present", True)

    _deep_set_if_missing(op, "T_inlet_K", eng.get("T_inlet_K", 293))

    # Reactor geometry from plant engineering
    if eng:
        _deep_set_if_missing(react, "bed_diameter_m", eng.get("bed_diameter_m"))
        _deep_set_if_missing(react, "freeboard_diameter_m", eng.get("freeboard_diameter_m"))
        h_tot = eng.get("total_height_m") or eng.get("reactor_height_to_cyclone_m") or eng.get("reactor_height_m")
        if h_tot is not None:
            _deep_set_if_missing(react, "height_m", h_tot)
            _deep_set_if_missing(react, "total_height_m", h_tot)
        h_bed = eng.get("expanded_bed_height_m") or eng.get("bed_height_m")
        if h_bed is not None:
            _deep_set_if_missing(react, "bed_height_m", h_bed)
        if eng.get("fuel_inlet_mm_above_nozzle") is not None:
            _deep_set_if_missing(react, "fuel_inlet_mm_above_nozzle", eng["fuel_inlet_mm_above_nozzle"])
        if eng.get("secondary_air_port_mm_above_nozzle") is not None:
            _deep_set_if_missing(react, "secondary_air_port_mm_above_nozzle", eng["secondary_air_port_mm_above_nozzle"])
        react["bed_height_source"] = eng.get("bed_height_source", react.get("bed_height_source"))
        react["hamel_engineering_note"] = eng.get("source")


def _plant_key_for_case(case_key: str) -> str:
    if "HTW_WESSELING" in case_key:
        return "HTW_WESSELING"
    if "HTW_BERRENRATH" in case_key:
        return "HTW_BERRENRATH"
    if "VTT_PRESSURIZED" in case_key:
        return "VTT_PRESSURIZED"
    if "VTT_ATMOSPHERIC" in case_key:
        return "VTT_ATMOSPHERIC"
    if "WSV400" in case_key:
        return "WSV400"
    return ""


def _get_table_row(hamel: dict, mapping: dict) -> dict | None:
    table_id = mapping["table"]
    sim = mapping["sim"]
    if table_id == "7.1":
        row = hamel["table_7_1_htw"]["rows"][sim]
        row = dict(row)
        row["fuel_analysis"] = row.get("fuel_analysis") or hamel["table_7_1_htw"]["rows"][sim]["fuel_analysis"]
        return row
    if table_id == "7.2":
        row = dict(hamel["table_7_2_vtt_pressurized"]["rows"][sim])
        fuel = str(row.get("fuel", "")).lower()
        if "braunk" in fuel:
            row.update({k: v for k, v in hamel["table_7_2_vtt_pressurized"]["coal_fuel_analysis"].items()})
        elif "torf" in fuel:
            row.update({k: v for k, v in hamel["table_7_2_vtt_pressurized"]["peat_fuel_analysis"].items()})
        else:
            row.update({k: v for k, v in hamel["table_7_2_vtt_pressurized"]["sawdust_fuel_analysis"].items()})
        if row.get("moisture_wt_pct") is not None:
            row["fuel_analysis"] = {
                "moisture_wt_pct": row["moisture_wt_pct"],
                "ash_dry_wt_pct": row.get("ash_dry_wt_pct"),
                "VM_daf_pct": row.get("VM_daf_pct"),
                "C_dry": row.get("C_dry"),
                "H_dry": row.get("H_dry"),
                "O_dry": row.get("O_dry"),
                "N_dry": row.get("N_dry"),
                "S_dry": row.get("S_dry"),
            }
        return row
    if table_id == "7.3":
        row = dict(hamel["table_7_3_vtt_atmospheric"]["rows"][sim])
        row.update(hamel["table_7_3_vtt_atmospheric"]["peat_fuel_analysis"])
        row["fuel_analysis"] = {
            "moisture_wt_pct": row.get("moisture_wt_pct"),
            **hamel["table_7_3_vtt_atmospheric"]["peat_fuel_analysis"],
        }
        return row
    if table_id == "7.4":
        return dict(hamel["table_7_4_wsv400"]["rows"][sim])
    return None


def enrich_cases(cases: dict, hamel: dict) -> tuple[dict, list[str]]:
    log: list[str] = []
    out = deepcopy(cases)
    for case_key, mapping in hamel["case_mapping"].items():
        if case_key not in out:
            log.append(f"SKIP missing case {case_key}")
            continue
        row = _get_table_row(hamel, mapping)
        if row is None:
            log.append(f"SKIP no row for {case_key}")
            continue
        plant_key = _plant_key_for_case(case_key)
        case = out[case_key]
        case["hamel_simulation_id"] = mapping["sim"]
        case["hamel_table_ref"] = f"Tabelle {mapping['table']}, Hamel (1999) Ch.7"
        _apply_table_row(case, row, plant_key=plant_key, hamel=hamel)
        log.append(f"OK {case_key} ← Table {mapping['table']} sim {mapping['sim']}")

    # HTW secondary injection notes (not in Table 7.1 numeric columns)
    for key, sec in (
        ("CASE_HTW_WESSELING_1", {"secondary_air_Nm3_h": 30.1, "agent": "Air", "xi": 0.6}),
        ("CASE_HTW_WESSELING_2", {"secondary_O2_Nm3_h": 37.5, "agent": "O2", "xi": 0.6}),
    ):
        ga = out[key]["inputs"]["gasification_agent"]
        if key.endswith("_1"):
            _deep_set_if_missing(ga, "secondary_air_Nm3_h", sec["secondary_air_Nm3_h"])
            _deep_set_if_missing(ga, "secondary_injection_agent", sec["agent"])
        else:
            _deep_set_if_missing(ga, "secondary_O2_Nm3_h", sec["secondary_O2_Nm3_h"])
            _deep_set_if_missing(ga, "secondary_injection_agent", sec["agent"])
        _deep_set_if_missing(ga, "secondary_injection_xi_reactor", sec["xi"])

    # VTT pressurized secondary xi from Fig 7.7 coordinate (0.6)
    for k in out:
        if k.startswith("CASE_VTT_PRESSURIZED") and "gasification_agent" in out[k].get("inputs", {}):
            ga = out[k]["inputs"]["gasification_agent"]
            if ga.get("secondary_air_Nm3_h"):
                _deep_set_if_missing(ga, "secondary_injection_xi_reactor", 0.6)

    # Update meta
    out["_meta"]["hamel_enrichment"] = {
        "script": "scripts/enrich_validation_cases_from_hamel.py",
        "hamel_table_ssot": "data/hamel_table7_extractions.json",
        "last_run_note": "Filled Table 7.1–7.4 inputs + §7.1 plant engineering; HTW steam abs. flow still not tabulated",
    }
    return out, log


def _count_missing(obj: object) -> int:
    if isinstance(obj, dict):
        return sum(_count_missing(v) for v in obj.values())
    if isinstance(obj, list):
        return sum(_count_missing(v) for v in obj)
    if isinstance(obj, str) and obj.startswith("__"):
        return 1
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cases = json.loads(CASES_PATH.read_text(encoding="utf-8"))
    hamel = json.loads(HAMEL_PATH.read_text(encoding="utf-8"))
    before = _count_missing(cases)
    enriched, log = enrich_cases(cases, hamel)
    after = _count_missing(enriched)

    for line in log:
        print(line)
    print(f"Missing markers: {before} → {after}")

    if not args.dry_run:
        CASES_PATH.write_text(json.dumps(enriched, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"Wrote {CASES_PATH}")


if __name__ == "__main__":
    main()
