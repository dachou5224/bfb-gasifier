#!/usr/bin/env python3
"""Merge Fig 7.10 / 7.13 / 7.17 parity digitization into validation_cases.json outputs.

SSOT: data/hamel_fig7_parity_digitized.json

Usage:
    python3 scripts/merge_fig7_parity_outputs.py
    python3 scripts/merge_fig7_parity_outputs.py --dry-run
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
CASES_PATH = ROOT / "data" / "validation_cases.json"
PARITY_PATH = ROOT / "data" / "hamel_fig7_parity_digitized.json"

SPECIES = ("CO", "CO2", "H2", "CH4", "H2O")


def _is_missing(val: object) -> bool:
    if val is None:
        return True
    if isinstance(val, str) and val.startswith("__"):
        return True
    return False


def _wet_to_dry(wet: dict[str, float]) -> dict[str, float]:
    h2o = wet.get("H2O", 0.0)
    dry_factor = 1.0 - h2o
    if dry_factor <= 0.0:
        raise ValueError(f"invalid H2O={h2o}")
    dry = {sp: wet[sp] / dry_factor for sp in ("CO", "CO2", "H2", "CH4") if sp in wet}
    return dry


def _build_outputs_from_parity(
    parity_case: dict[str, Any],
    *,
    basis: str,
    fig: str,
    overwrite_missing_only: bool,
    existing: dict[str, Any] | None,
) -> dict[str, Any]:
    existing = existing or {}
    out: dict[str, Any] = dict(existing)

    t = parity_case.get("T_FRB_K", {})
    t_meas = t.get("meas")
    if t_meas is not None and (overwrite_missing_only and _is_missing(out.get("exit_temperature_K")) or not overwrite_missing_only):
        if _is_missing(out.get("exit_temperature_K")):
            out["exit_temperature_K"] = int(round(t_meas))
            out["exit_temperature_note"] = (
                f"Mean freeboard T (T_FRB) from Fig. {fig} parity Messung axis; "
                "may differ from axial outlet profile at ξ=1"
            )

    xc = parity_case.get("carbon_conversion_pct", {})
    if xc.get("meas") is not None and _is_missing(out.get("carbon_conversion_pct")):
        out["carbon_conversion_pct"] = float(xc["meas"])

    wet: dict[str, float] = {}
    dry: dict[str, float] = {}
    for sp in SPECIES:
        pt = parity_case.get(sp, {})
        if pt.get("meas") is None:
            continue
        val = float(pt["meas"])
        if basis == "wet_mol_frac":
            wet[sp] = val
        elif sp != "H2O":
            dry[sp] = val

    if basis == "wet_mol_frac" and wet:
        if _is_missing(out.get("exit_gas_wet_mol_frac")):
            wet_out = {k: round(v, 3) for k, v in wet.items()}
            wet_out["N2"] = "__CALCULATED__"
            out["exit_gas_wet_mol_frac"] = wet_out
        if _is_missing(out.get("exit_gas_dry_mol_frac")) and "H2O" in wet:
            dry_calc = _wet_to_dry(wet)
            out["exit_gas_dry_mol_frac"] = {k: round(v, 3) for k, v in dry_calc.items()}
            out["exit_gas_dry_mol_frac"]["_note"] = "dry = wet / (1 - H2O) from parity Messung"

    if basis == "dry_mol_frac" and dry:
        if _is_missing(out.get("exit_gas_dry_mol_frac")):
            dry_out = {k: round(v, 3) for k, v in dry.items()}
            dry_out["N2"] = "__CALCULATED__  (balance for air gasification)"
            dry_out["_note"] = "Fig. 7.17 reports dry basis (mol/mol trocken)"
            out["exit_gas_dry_mol_frac"] = dry_out

    parity_block: dict[str, Any] = {
        "available": True,
        "source": f"Fig. {fig}, digitized 2026-06-04",
        "messung_axis": {},
        "sim_axis": {},
    }
    if t:
        parity_block["T_FRB_K"] = t
    if xc:
        parity_block["carbon_conversion_pct"] = xc
    for sp in SPECIES:
        if sp in parity_case:
            parity_block["messung_axis"][sp] = parity_case[sp].get("meas")
            parity_block["sim_axis"][sp] = parity_case[sp].get("sim")

    out["parity_plot_data"] = parity_block
    out["_precision"] = out.get("_precision") or f"Digitized from Fig. {fig} parity (Messung axis), ±5% reading error"
    if str(out.get("_precision", "")).startswith("__"):
        out["_precision"] = f"Digitized from Fig. {fig} parity (Messung axis), ±5% reading error"

    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    parity = json.loads(PARITY_PATH.read_text(encoding="utf-8"))
    cases = json.loads(CASES_PATH.read_text(encoding="utf-8"))

    fig710 = parity["fig_7_10_vtt_pressurized"]
    fig713 = parity["fig_7_13_vtt_atmospheric"]
    fig717 = parity["fig_7_17_wsv400"]
    mapping = parity["case_id_map"]

    fig_sections = {
        "7.10": fig710,
        "7.13": fig713,
        "7.17": fig717,
    }

    updated: list[str] = []
    for case_id, ref in mapping.items():
        if case_id not in cases:
            continue
        fig = ref["fig"]
        sim = ref["sim"]
        section = fig_sections[fig]
        pdata = section["cases"][sim]
        basis = section["species_basis"]

        existing_out = cases[case_id].get("outputs", {})
        new_out = _build_outputs_from_parity(
            pdata,
            basis=basis,
            fig=fig,
            overwrite_missing_only=True,
            existing=deepcopy(existing_out) if isinstance(existing_out, dict) else {},
        )
        cases[case_id]["outputs"] = new_out
        updated.append(case_id)
        print(f"OK {case_id} ← Fig. {fig} sim {sim}")

    if not args.dry_run:
        CASES_PATH.write_text(json.dumps(cases, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"Wrote {CASES_PATH} ({len(updated)} cases)")


if __name__ == "__main__":
    main()
