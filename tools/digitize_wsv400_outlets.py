#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由记录的像素坐标复算 WSV400 出口读数，并与 data/digitized 中的 JSON 对照。

渲染：PyMuPDF 1.26.5，zoom = 4（4 px/pt）。坐标原点是裁图左上角，x 向右，y 向下。
网格线像素是该线墨迹的质心。Bild 7.17 的点是空心方框四边描边中线的中点，
横轴为实验值，纵轴为模拟值。Bild 7.16 只取 ξ = 1 紧下方的曲线或符号中心。

不读取 PDF。换算只依赖本文件里的像素和刻度，便于离线复算。
"""

from __future__ import annotations

import json
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Bild 7.17：PDF 第 83 页左半，clip (40, 20, 420, 490) pt。
# 每条网格： (像素, 刻度)。x 为竖直网格（实验轴），y 为水平网格（模拟轴）。
FIG_717 = {
    "T_FRB": {
        "x": [(161.0, 1000), (319.0, 1100), (476.0, 1200), (631.5, 1300)],
        "y": [(70.0, 1300), (224.6, 1200), (379.0, 1100), (533.0, 1000)],
        "ndigits": 0,
        # 方框中心 (x 实验, y 模拟)
        "centers": {
            "A": (264.5, 459.5),
            "B": (252.5, 538.5),
            "C": (204.0, 491.0),
            "D": (235.25, 429.5),
            "E": (279.0, 388.0),
            "F": (169.0, 465.75),
        },
    },
    "CO": {
        "x": [(162.9, 0.1), (397.5, 0.2), (631.0, 0.3)],
        "y": [(678.6, 0.3), (907.5, 0.2), (1138.5, 0.1)],
        "ndigits": 3,
        "centers": {
            "A": (330.5, 980.5),
            "B": (292.0, 913.5),
            "C": (319.0, 948.0),
            "D": (345.75, 990.5),
            "E": (367.5, 1025.5),
            "F": (269.5, 948.5),
        },
    },
    "CO2": {
        "x": [(865.5, 0.0), (1099.5, 0.1), (1326.9, 0.2)],
        "y": [(673.0, 0.2), (903.9, 0.1), (1132.1, 0.0)],
        "ndigits": 3,
        "centers": {
            "A": (1206.5, 851.5),
            "B": (1272.25, 869.0),
            "C": (1258.5, 847.25),
            "D": (1239.5, 835.5),
            "E": (1217.25, 825.5),
            "F": (1235.25, 866.25),
        },
    },
    "H2": {
        "x": [(167.5, 0.0), (404.0, 0.1), (640.0, 0.2)],
        "y": [(1295.4, 0.2), (1525.6, 0.1), (1756.0, 0.0)],
        "ndigits": 3,
        "centers": {
            "A": (373.5, 1568.5),
            "B": (426.75, 1549.5),
            "C": (443.0, 1533.0),
            "D": (457.0, 1549.5),
            "E": (455.5, 1576.5),
            "F": (331.0, 1534.5),
        },
    },
    "CH4": {
        "x": [(868.5, 0.00), (1106.0, 0.05), (1337.0, 0.10)],
        "y": [(1291.0, 0.10), (1521.9, 0.05), (1751.5, 0.00)],
        "ndigits": 3,
        "centers": {
            "A": (1092.5, 1652.0),
            "B": (1077.0, 1602.0),
            "C": (1090.0, 1630.5),
            "D": (1051.5, 1655.0),
            "E": (1025.0, 1672.0),
            "F": (1058.5, 1628.0),
        },
    },
}

# Bild 7.16：PDF 第 82 页右半，clip (420, 20, 820, 500) pt。出口在顶边 ξ = 1。
FIG_716 = {
    "xi_y": [(879.00, 1.0), (1025.95, 0.8), (1174.03, 0.6), (1321.91, 0.4), (1470.90, 0.2), (1617.99, 0.0)],
    "gas_x": [(530.10, 0.0), (661.91, 0.1), (792.61, 0.2), (925.88, 0.3)],
    "T_x": [(1027.52, 800), (1158.99, 1000), (1291.12, 1200), (1422.37, 1400)],
    # 出口符号的 x。y 用来说明它紧贴哪条 ξ 线，不参与温度或 CO 的换算。
    "outlets": {
        "CO_sim": {"x": 754.5, "y_span": (882, 889), "ndigits": 3},
        "CO_meas": {"x": 761.25, "y_span": (889, 901), "ndigits": 3},
        "T_sim": {"x": 1159.0, "y_span": (884, 888), "ndigits": 0},
        "T_meas": {"x": 1166.5, "y_span": (890, 904), "ndigits": 0},
    },
}


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


def values_717():
    out = {}
    for variable, panel in FIG_717.items():
        sx, ix = fit_axis(panel["x"])
        sy, iy = fit_axis(panel["y"])
        nd = panel["ndigits"]
        out[variable] = {
            letter: [apply_axis(sx, ix, x, nd), apply_axis(sy, iy, y, nd)]
            for letter, (x, y) in panel["centers"].items()
        }
    return out


def values_716():
    sx, ix = fit_axis(FIG_716["gas_x"])
    tx, it = fit_axis(FIG_716["T_x"])
    outlets = FIG_716["outlets"]
    return {
        "CO_sim": apply_axis(sx, ix, outlets["CO_sim"]["x"], 3),
        "CO_meas": apply_axis(sx, ix, outlets["CO_meas"]["x"], 3),
        "T_sim": apply_axis(tx, it, outlets["T_sim"]["x"], 0),
        "T_meas": apply_axis(tx, it, outlets["T_meas"]["x"], 0),
    }


def load_recorded():
    fig717 = json.loads((ROOT / "data/digitized/fig_7_17.json").read_text())
    fig716 = json.loads((ROOT / "data/digitized/fig_7_16.json").read_text())
    recorded = {key: {} for key in FIG_717}
    for case in fig717["cases"]:
        letter = case["case_label_in_figure"]
        for key in FIG_717:
            recorded[key][letter] = case[key if key != "T_FRB" else "T_FRB"]
    series = {(item["kind"], item["variable"]): item["points"][0][1] for item in fig716["series"]}
    recorded716 = {
        "CO_sim": series[("sim", "CO")],
        "CO_meas": series[("meas", "CO")],
        "T_sim": series[("sim", "T")],
        "T_meas": series[("meas", "T")],
    }
    return recorded, recorded716


def check():
    got = values_717()
    got716 = values_716()
    recorded, recorded716 = load_recorded()
    mismatches = []
    for variable, cases in got.items():
        for letter, pair in cases.items():
            if pair != recorded[variable][letter]:
                mismatches.append(f"7.17 {variable} {letter}: 复算 {pair}，JSON {recorded[variable][letter]}")
    for key, value in got716.items():
        if value != recorded716[key]:
            mismatches.append(f"7.16 {key}: 复算 {value}，JSON {recorded716[key]}")
    report = (ROOT / "data/digitized/REPORT.md").read_text()
    pixels = []
    for panel in FIG_717.values():
        pixels.extend(p for p, _ in panel["x"])
        pixels.extend(p for p, _ in panel["y"])
        for x, y in panel["centers"].values():
            pixels.extend((x, y))
    for sample in FIG_716["gas_x"] + FIG_716["T_x"] + FIG_716["xi_y"]:
        pixels.append(sample[0])
    for outlet in FIG_716["outlets"].values():
        pixels.append(outlet["x"])
    for pixel in pixels:
        token = f"{pixel:.2f}".rstrip("0").rstrip(".")
        if re.search(rf"(?<![\d.]){re.escape(token)}0*(?!\d)", report) is None:
            mismatches.append(f"REPORT 缺少像素 {pixel}")
    return mismatches


def main():
    mismatches = check()
    if mismatches:
        print("\n".join(mismatches), file=sys.stderr)
        return 1
    print("复算与 fig_7_17.json、fig_7_16.json 一致，REPORT 含全部像素。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
