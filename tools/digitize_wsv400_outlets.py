#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 JSON 里的 calibration 复算 WSV400 出口读数，并与文件中的数值对照。

像素和刻度只记在 data/digitized/fig_7_17.json、fig_7_16.json 的 calibration 中。
本脚本做线性拟合和四舍五入（0.5 进位），不读取 PDF。
"""

from __future__ import annotations

import json
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


def check():
    doc717 = json.loads((ROOT / "data/digitized/fig_7_17.json").read_text())
    doc716 = json.loads((ROOT / "data/digitized/fig_7_16.json").read_text())
    mismatches = []
    for key, value in doc717["est_error"].items():
        if not isinstance(value, (int, float)):
            mismatches.append(f"7.17 est_error.{key} 不是数值")
        if key not in doc717["units"]:
            mismatches.append(f"7.17 est_error.{key} 没有对应的 units")
    for key, value in doc716["est_error"].items():
        if not isinstance(value, (int, float)):
            mismatches.append(f"7.16 est_error.{key} 不是数值")
    got = values_717(doc717)
    recorded = recorded_717(doc717)
    for variable, cases in got.items():
        for letter, pair in cases.items():
            if pair != recorded[variable][letter]:
                mismatches.append(f"7.17 {variable} {letter}: 复算 {pair}，JSON {recorded[variable][letter]}")
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
