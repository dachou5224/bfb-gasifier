#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""KURKELA_1992 Fig. 1（p. 4）褐煤累积粒径曲线的数字化与复算。

两种用法：
- 默认（不读 PDF）：用 data/digitized/kurkela1992_fig1.json 里 calibration 记录的坐标框像素与圆心像素，
  复算累积分布与粒级质量分数，并与文件中的数值对照。测试用这一模式。
- ``--trace <kurkela1992.pdf>``：在本地渲染页面，检测坐标框，用环形模板定位褐煤的空心圆，打印
  calibration 与 points，供写入 JSON。PDF 只在本地，不入库（CHARTER §11）。

横轴为对数坐标（坐标框左端 0,1 mm、右端 5 mm），纵轴为线性坐标（下端 0 %、上端 100 %）。
只读褐煤（BROWN COAL，空心圆、点线），不读泥炭与锯末。
"""

from __future__ import annotations

import json
import math
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data" / "digitized" / "kurkela1992_fig1.json"

PAGE_INDEX = 3  # p. 4，零起
ZOOM = 8
CLIP_PT = (60, 80, 320, 265)  # 图区（pt），只为裁切
X_FRAME_MM = (0.1, 5.0)
Y_FRAME_PCT = (0.0, 100.0)
SIEVES_MM = [0.1, 0.25, 0.5, 0.71, 1.0, 2.0, 3.15]  # 圆心横坐标应落在的标准筛孔
SIEVE_TOLERANCE = 0.06  # 圆心换算的粒径与筛孔的相对偏差上限；0,1 mm 的圆压在左边框上，圆心偏右约 5 %
DIGITS = 1


def round_half_up(value: float, ndigits: int) -> float:
    quantum = Decimal("1") if ndigits == 0 else Decimal("0." + "0" * (ndigits - 1) + "1")
    return float(Decimal(str(value)).quantize(quantum, rounding=ROUND_HALF_UP))


def axes(cal: dict):
    (xl, xr), (yb, yt) = cal["x_frame_px"], cal["y_frame_px"]
    decade = (xr - xl) / math.log10(X_FRAME_MM[1] / X_FRAME_MM[0])

    def size_mm(px: float) -> float:
        return X_FRAME_MM[0] * 10 ** ((px - xl) / decade)

    def cum_pct(py: float) -> float:
        return (yb - py) / (yb - yt) * (Y_FRAME_PCT[1] - Y_FRAME_PCT[0])

    return size_mm, cum_pct


def recompute(doc: dict) -> tuple[list, list]:
    """由 calibration 复算 [筛孔 mm, 累积 %] 与粒级质量分数 %。"""
    cal = doc["calibration"]
    size_mm, cum_pct = axes(cal)
    points = []
    for c, sieve in zip(cal["circle_centers_px"], SIEVES_MM):
        points.append([sieve, round_half_up(min(cum_pct(c[1]), 100.0), DIGITS)])
    fractions = []
    previous = 0.0
    for sieve, cum in points:
        fractions.append(round_half_up(cum - previous, DIGITS))
        previous = cum
    return points, fractions


def valid_error(value) -> bool:
    """误差须为有限、非负的实数；布尔值不算（bool 是 int 的子类）。"""
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value >= 0)


def check_frame(cal: dict) -> list:
    """坐标框：两轴各两个有限实数端点，正跨度，方向正确（右 > 左，下 > 上）。"""
    problems = []
    for key in ("x_frame_px", "y_frame_px"):
        pair = cal.get(key)
        if (not isinstance(pair, (list, tuple)) or len(pair) != 2
                or not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in pair)):
            problems.append(f"{key} 须为两个有限实数")
    if problems:
        return problems
    (xl, xr), (yb, yt) = cal["x_frame_px"], cal["y_frame_px"]
    if not xr > xl:
        problems.append("坐标框横向跨度须为正（右端像素大于左端）")
    if not yb > yt:
        problems.append("坐标框纵向跨度须为正（下端像素大于上端）")
    return problems


def check(doc: dict | None = None) -> list:
    doc = doc if doc is not None else json.loads(JSON_PATH.read_text(encoding="utf-8"))
    problems = []
    for key in ("size_mm", "cumulative_pct"):
        if not valid_error(doc.get("est_error", {}).get(key)):
            problems.append(f"est_error.{key} 须为有限、非负的实数")
    cal = doc["calibration"]
    frame_problems = check_frame(cal)
    if frame_problems:
        # 坐标框无效时不构造坐标变换，直接返回诊断
        return problems + frame_problems
    size_mm, _ = axes(cal)
    if len(cal["circle_centers_px"]) != len(SIEVES_MM):
        problems.append("圆心个数与筛孔个数不一致")
    for c, sieve in zip(cal["circle_centers_px"], SIEVES_MM):
        d = size_mm(c[0])
        if abs(d / sieve - 1) > SIEVE_TOLERANCE:
            problems.append(f"筛孔 {sieve} mm：圆心换算为 {d:.3f} mm")
    points, fractions = recompute(doc)
    if points != doc["points"]:
        problems.append(f"points 与 calibration 复算不一致：{points} != {doc['points']}")
    if fractions != doc["class_mass_fractions_pct"]:
        problems.append(f"粒级质量分数与复算不一致：{fractions} != {doc['class_mass_fractions_pct']}")
    if any(b < a for (_, a), (_, b) in zip(points, points[1:])):
        problems.append("累积分布不单调")
    if abs(sum(fractions) - 100.0) > 0.2:
        problems.append(f"粒级质量分数之和为 {sum(fractions):.1f} %")
    return problems


def _groups(indices, gap=3):
    out = []
    for t in indices:
        if out and t - out[-1][-1] <= gap:
            out[-1].append(t)
        else:
            out.append([t])
    return [float(np.mean(g)) for g in out]


def trace(pdf_path: str) -> dict:
    import fitz  # PyMuPDF，仅本地追踪时需要
    from scipy.signal import fftconvolve

    page = fitz.open(pdf_path)[PAGE_INDEX]
    pix = page.get_pixmap(matrix=fitz.Matrix(ZOOM, ZOOM), clip=fitz.Rect(*CLIP_PT), colorspace=fitz.csGRAY)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width).astype(float)
    h, w = img.shape
    dark = img < 120
    rows = _groups([y for y in range(h) if dark[y, 300:1850].sum() > 0.9 * 1550])
    cols = _groups([x for x in range(w) if dark[100:1300, x].sum() > 0.9 * 1200])
    y_frame = (max(rows), min(rows))
    x_frame = (min(cols), max(cols))
    cal = {"render": {"page": PAGE_INDEX + 1, "zoom": ZOOM, "clip_pt": list(CLIP_PT), "dark_threshold": 120},
           "x_frame_px": [round(x_frame[0], 1), round(x_frame[1], 1)],
           "y_frame_px": [round(y_frame[0], 1), round(y_frame[1], 1)]}
    size_mm, cum_pct = axes(cal)
    decade = (x_frame[1] - x_frame[0]) / math.log10(X_FRAME_MM[1] / X_FRAME_MM[0])

    # 环形模板定位空心圆；在每个筛孔的横坐标附近搜索
    ring = (img < 140).astype(float)
    r = 17
    yy, xx = np.mgrid[-r - 4:r + 5, -r - 4:r + 5]
    rr = np.hypot(yy, xx)
    template = np.where(np.abs(rr - r) <= 1.5, 1.0, 0.0) - np.where(rr < r - 4, 0.6, 0.0)
    score = fftconvolve(ring, template[::-1, ::-1], mode="same")
    # 累积曲线单调上升：每个圆只在上一个圆的上方（像素 y 更小）搜索
    centers = []
    y_floor = int(y_frame[0])
    for sieve in SIEVES_MM:
        x_expect = x_frame[0] + decade * math.log10(sieve / X_FRAME_MM[0])
        x_lo, x_hi = int(x_expect - 30), int(x_expect + 30)
        y_top = int(y_frame[1]) - 10
        window = score[y_top:y_floor, x_lo:x_hi]
        j = np.unravel_index(int(np.argmax(window)), window.shape)
        cy, cx = y_top + j[0], x_lo + j[1]
        centers.append([float(cx), float(cy)])
        y_floor = cy + 5
    cal["circle_radius_px"] = r
    cal["circle_centers_px"] = centers
    doc = {"calibration": cal}
    doc["points"], doc["class_mass_fractions_pct"] = recompute(doc)
    return doc


def main(argv: list) -> int:
    if len(argv) >= 2 and argv[0] == "--trace":
        print(json.dumps(trace(argv[1]), ensure_ascii=False, indent=2))
        return 0
    problems = check()
    for p in problems:
        print("问题：" + p)
    print("通过" if not problems else f"{len(problems)} 个问题")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
