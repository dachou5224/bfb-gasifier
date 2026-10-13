#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 JSON 里的 calibration 复算 WSV400 出口读数，并与文件中的数值对照。

像素和刻度只记在 data/digitized/fig_7_17.json、fig_7_16.json 的 calibration 中。
本脚本做线性拟合和四舍五入（0.5 进位），不读取 PDF。
"""

from __future__ import annotations

import json
import math
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def round_half_up(value: float, ndigits: int):
    quantum = Decimal("1") if ndigits == 0 else Decimal("0." + "0" * (ndigits - 1) + "1")
    rounded = Decimal(str(value)).quantize(quantum, rounding=ROUND_HALF_UP)
    if ndigits == 0:
        return int(rounded)
    return float(rounded)


def fit_axis(samples):
    pixels = np.array([p for p, _ in samples], dtype=float)
    values = np.array([v for _, v in samples], dtype=float)
    slope, intercept = np.polyfit(pixels, values, 1)
    return float(slope), float(intercept)


def apply_axis(slope: float, intercept: float, pixel: float, ndigits: int):
    return round_half_up(slope * pixel + intercept, ndigits)


def values_717(doc):
    cal = doc["calibration"]
    out = {}
    for variable, grid in cal["grids"].items():
        sx, ix = fit_axis(grid["x"])
        sy, iy = fit_axis(grid["y"])
        nd = cal["rounding_digits"][variable]
        out[variable] = {
            letter: [apply_axis(sx, ix, x, nd), apply_axis(sy, iy, y, nd)]
            for letter, (x, y) in cal["centers"][variable].items()
        }
    return out


def values_716(doc):
    cal = doc["calibration"]
    sx, ix = fit_axis(cal["grids"]["CO_x"])
    tx, it = fit_axis(cal["grids"]["T_x"])
    points = cal["points"]
    nd = cal["rounding_digits"]
    return {
        "CO_sim": apply_axis(sx, ix, points["CO_sim"]["x"], nd["CO"]),
        "CO_meas": apply_axis(sx, ix, points["CO_meas"]["x"], nd["CO"]),
        "T_sim": apply_axis(tx, it, points["T_sim"]["x"], nd["T"]),
        "T_meas": apply_axis(tx, it, points["T_meas"]["x"], nd["T"]),
    }


def recorded_717(doc):
    out = {key: {} for key in doc["calibration"]["grids"]}
    for case in doc["cases"]:
        letter = case["case_label_in_figure"]
        for key in out:
            out[key][letter] = case[key]
    return out


OUTLET_716 = (("sim", "T"), ("meas", "T"), ("sim", "CO"), ("meas", "CO"))


def is_error_number(value):
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
        and value >= 0
    )


def series_variables(doc):
    keys = []
    for item in doc.get("series", []):
        variable = item.get("variable")
        if variable and variable not in keys:
            keys.append(variable)
    return keys


def quantity_keys_717(doc):
    keys = []
    for case in doc.get("cases", []):
        for key, value in case.items():
            if isinstance(value, list) and key not in keys:
                keys.append(key)
    return keys


def check_error_map(doc, label, required):
    mismatches = []
    errors = doc.get("est_error")
    units = doc.get("units")
    if not isinstance(errors, dict):
        mismatches.append(f"{label} 缺少 est_error")
        errors = {}
    if not isinstance(units, dict):
        mismatches.append(f"{label} 缺少 units")
        units = {}
    for key in required:
        if key not in errors:
            mismatches.append(f"{label} est_error 缺少 {key}")
        elif not is_error_number(errors[key]):
            mismatches.append(f"{label} est_error.{key} 不是有限非负数值")
        unit = units.get(key)
        if not isinstance(unit, str) or not unit:
            mismatches.append(f"{label} units 缺少 {key}")
    return mismatches


def check_716_series(doc):
    mismatches = []
    seen = set()
    for item in doc.get("series", []):
        key = (item.get("kind"), item.get("variable"))
        points = item.get("points")
        if key in seen:
            mismatches.append(f"7.16 序列 {key[0]}/{key[1]} 重复")
            continue
        seen.add(key)
        outlet = (
            isinstance(points, list)
            and len(points) == 1
            and isinstance(points[0], (list, tuple))
            and len(points[0]) == 2
            and points[0][0] == 1
        )
        if not outlet:
            mismatches.append(f"7.16 序列 {key[0]}/{key[1]} 必须恰有一个 ξ = 1 的点")
    for kind, variable in OUTLET_716:
        if (kind, variable) not in seen:
            mismatches.append(f"7.16 缺少序列 {kind}/{variable}")
    return mismatches


def check_716_units(doc):
    units = doc.get("units") if isinstance(doc.get("units"), dict) else {}
    mismatches = []
    for item in doc.get("series", []):
        variable = item.get("variable")
        if variable in units and item.get("unit") != units[variable]:
            mismatches.append(f"7.16 {variable} 的序列单位与 units 不一致")
    return mismatches


def recorded_716(doc):
    series = {(item["kind"], item["variable"]): item["points"][0][1] for item in doc["series"]}
    return {
        "CO_sim": series[("sim", "CO")],
        "CO_meas": series[("meas", "CO")],
        "T_sim": series[("sim", "T")],
        "T_meas": series[("meas", "T")],
    }


def pixels_of(doc717, doc716):
    pixels = []
    cal = doc717["calibration"]
    for grid in cal["grids"].values():
        pixels.extend(p for p, _ in grid["x"])
        pixels.extend(p for p, _ in grid["y"])
    for centers in cal["centers"].values():
        for x, y in centers.values():
            pixels.extend((x, y))
    cal16 = doc716["calibration"]
    for sample in cal16["grids"]["xi_y"] + cal16["grids"]["CO_x"] + cal16["grids"]["T_x"]:
        pixels.append(sample[0])
    for point in cal16["points"].values():
        pixels.append(point["x"])
    return pixels


def check(doc717=None, doc716=None):
    if doc717 is None:
        doc717 = json.loads((ROOT / "data/digitized/fig_7_17.json").read_text())
    if doc716 is None:
        doc716 = json.loads((ROOT / "data/digitized/fig_7_16.json").read_text())
    mismatches = []
    mismatches.extend(check_error_map(doc717, "7.17", quantity_keys_717(doc717)))
    mismatches.extend(check_error_map(doc716, "7.16", series_variables(doc716)))
    mismatches.extend(check_716_units(doc716))
    series_errors = check_716_series(doc716)
    mismatches.extend(series_errors)
    got = values_717(doc717)
    recorded = recorded_717(doc717)
    for variable, cases in got.items():
        for letter, pair in cases.items():
            if pair != recorded[variable][letter]:
                mismatches.append(f"7.17 {variable} {letter}: 复算 {pair}，JSON {recorded[variable][letter]}")
    if not series_errors:
        got716 = values_716(doc716)
        recorded716 = recorded_716(doc716)
        for key, value in got716.items():
            if value != recorded716[key]:
                mismatches.append(f"7.16 {key}: 复算 {value}，JSON {recorded716[key]}")
    report = (ROOT / "data/digitized/REPORT.md").read_text()
    for pixel in pixels_of(doc717, doc716):
        token = f"{pixel:.2f}".rstrip("0").rstrip(".")
        if re.search(rf"(?<![\d.]){re.escape(token)}0*(?!\d)", report) is None:
            mismatches.append(f"REPORT 缺少像素 {pixel}")
    return mismatches


def main():
    mismatches = check()
    if mismatches:
        print("\n".join(mismatches), file=sys.stderr)
        return 1
    print("复算与 fig_7_17.json、fig_7_16.json 的 calibration 一致，REPORT 含全部像素。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
