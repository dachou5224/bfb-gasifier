#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bild 5.11（PDF 56 右 / S. 93）K_eq 曲线的数字化与复算。

两种用法：
- 默认（不读 PDF）：用 data/digitized/fig_5_11.json 里 calibration 记录的网格线像素和曲线中心像素，
  做线性拟合，复算 points，并与文件中的数值对照。测试用这一模式。
- ``--trace <Hamel PDF>``：在本地渲染页面，检测网格线，沿实线中心逐列追踪，打印 calibration 与 points，
  供写入 JSON。PDF 只在本地，不入库（CHARTER §11）。

只读图中实线 K_eq（按式 5.43 计算的平衡常数），不读实验比值 K_exp 的散点。
"""

from __future__ import annotations

import json
import math
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data" / "digitized" / "fig_5_11.json"

PAGE_INDEX = 55  # PDF 56，零起
ZOOM = 6
DARK = 110  # 灰度阈值：实线与竖网格线
GRID_LIGHT = 160  # 横网格线较淡
X_LABELS = [1050, 1100, 1150, 1200, 1250]  # 标了数字的竖网格线；最左的断轴不用
Y_LABELS = [1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
T_POINTS = [1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250]
WINDOW = 3  # 取点列两侧各 3 列的中位数
DIGITS = 3


def round_half_up(value: float, ndigits: int) -> float:
    quantum = Decimal("0." + "0" * (ndigits - 1) + "1")
    return float(Decimal(str(value)).quantize(quantum, rounding=ROUND_HALF_UP))


def fit(pixels, values):
    slope, intercept = np.polyfit(np.asarray(pixels, float), np.asarray(values, float), 1)
    return float(slope), float(intercept)


def recompute(doc: dict) -> list:
    """由 calibration 复算 [T_K, K_eq]。"""
    cal = doc["calibration"]
    sy, iy = fit([p for p, _ in cal["y_grid"]], [v for _, v in cal["y_grid"]])
    return [[c["T_K"], round_half_up(sy * c["pixel_y_center"] + iy, DIGITS)] for c in cal["centerline"]]


MAX_GRID_RESIDUAL_PX = 2.0  # 当前标定实测最大残差：横轴 1,0 px，纵轴 1,25 px
MAX_SPACING_DEVIATION = 0.02  # 相邻网格线间距与平均间距之差不超过 2 %（当前实测在 1 % 以内）


def check_grid(name: str, grid: list, pixel_rises_with_value: bool) -> list:
    """网格标定本身：标签与像素严格单调且方向正确、线性拟合残差、等间距。"""
    if len(grid) < 3:
        return [f"{name}：网格线少于 3 条"]
    problems = []
    pixels = np.asarray([p for p, _ in grid], float)
    values = np.asarray([v for _, v in grid], float)
    if not (np.all(np.isfinite(pixels)) and np.all(np.isfinite(values))):
        return [f"{name}：像素或标签不是有限实数"]
    if np.ptp(pixels) == 0 or np.ptp(values) == 0:
        return [f"{name}：像素或标签的跨度为 0"]
    dv = np.diff(values)
    if not (np.all(dv > 0) or np.all(dv < 0)):
        problems.append(f"{name}：标签不严格单调")
    steps = np.diff(pixels) / np.where(dv == 0, np.nan, dv)
    if not (np.all(steps > 0) if pixel_rises_with_value else np.all(steps < 0)):
        problems.append(f"{name}：像素与标签的方向不一致或不单调")
    slope, intercept = fit(values, pixels)
    residual = float(np.abs(pixels - (slope * values + intercept)).max())
    if residual > MAX_GRID_RESIDUAL_PX:
        problems.append(f"{name}：拟合残差 {residual:.2f} px 超过 {MAX_GRID_RESIDUAL_PX} px")
    gaps = np.abs(np.diff(pixels))
    if float(np.abs(gaps / gaps.mean() - 1).max()) > MAX_SPACING_DEVIATION:
        problems.append(f"{name}：网格间距不均匀 {np.round(gaps, 1).tolist()}")
    return problems


def valid_error(value) -> bool:
    """误差须为有限、非负的实数；布尔值不算（bool 是 int 的子类）。"""
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value >= 0)


def check(doc: dict | None = None) -> list:
    doc = doc if doc is not None else json.loads(JSON_PATH.read_text(encoding="utf-8"))
    problems = []
    err = doc.get("est_error", {})
    for key in ("T_K", "K_eq"):
        if not valid_error(err.get(key)):
            problems.append(f"est_error.{key} 须为有限、非负的实数")
    cal = doc["calibration"]
    grid_problems = (check_grid("x_grid", cal["x_grid"], pixel_rises_with_value=True)
                     + check_grid("y_grid", cal["y_grid"], pixel_rises_with_value=False))
    if grid_problems:
        # 网格标定无效时不再拟合与复算，直接返回诊断
        return problems + grid_problems
    sx, ix = fit([p for p, _ in cal["x_grid"]], [v for _, v in cal["x_grid"]])
    for c in cal["centerline"]:
        t_from_pixel = sx * c["pixel_x"] + ix
        if abs(t_from_pixel - c["T_K"]) > 1.0:  # 1250 K 取框内列，偏 0,6 K
            problems.append(f"T={c['T_K']}：取点列换算为 {t_from_pixel:.2f} K")
    if recompute(doc) != doc["points"]:
        problems.append(f"points 与 calibration 复算不一致：{recompute(doc)} != {doc['points']}")
    return problems


def _groups(indices, gap=3):
    out = []
    for t in indices:
        if out and t - out[-1][-1] <= gap:
            out[-1].append(t)
        else:
            out.append([t])
    return [float(np.mean(g)) for g in out]


def _runs(column, offset):
    ys = np.where(column)[0]
    out = []
    for y in ys:
        if out and y - out[-1][-1] <= 1:
            out[-1].append(y)
        else:
            out.append([y])
    return [(r[0] + offset, r[-1] + offset) for r in out]


def trace(pdf_path: str) -> dict:
    import fitz  # PyMuPDF，仅本地追踪时需要

    page = fitz.open(pdf_path)[PAGE_INDEX]
    pix = page.get_pixmap(matrix=fitz.Matrix(ZOOM, ZOOM), colorspace=fitz.csGRAY)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    right = img[:, img.shape[1] // 2:]
    dark = right < DARK
    h, w = right.shape

    # 竖网格线：在图区高度内几乎连续的深色列
    y0, y1 = int(h * 0.18), int(h * 0.62)
    col = dark[y0:y1].sum(0)
    vgrid = _groups([x for x in range(w) if col[x] > 0.8 * (y1 - y0)])
    # 左端为断轴边框，标数字的是其后 5 条
    vgrid = vgrid[-5:]
    # 横网格线：在两条外侧竖网格线之间几乎连续的较淡行
    xl, xr = int(vgrid[0]) - 200, int(vgrid[-1]) - 50
    row = (right < GRID_LIGHT)[:, xl:xr].sum(1)
    hgrid = [y for y in _groups([y for y in range(h) if row[y] > 0.6 * (xr - xl)]) if 640 < y < 1900]
    hgrid = hgrid[-8:]

    sx, ix = fit(X_LABELS, vgrid)  # T → 像素
    sy, iy = fit(Y_LABELS, hgrid)  # K_eq → 像素

    # 从 1100 K（曲线穿过 1,0 网格线处）向两侧追踪实线中心；穿过网格线或散点时按斜率外推
    centers = {}

    def follow(x, y, step, slope):
        gap, hist = 0, []
        while vgrid[0] - 260 < x < vgrid[-1] + 2:
            runs = _runs(dark[640:1850, x], 640)
            near = [(abs((a + b) / 2 - y), (a + b) / 2) for a, b in runs if 6 <= b - a + 1 <= 16]
            near = [c for c in near if c[0] <= 6]
            if near:
                _, c = min(near)
                centers[x] = c
                hist = (hist + [(x, c)])[-40:]
                if len(hist) > 10:
                    slope = np.polyfit([p for p, _ in hist], [q for _, q in hist], 1)[0]
                y, gap = c + slope * step, 0
            else:
                gap += 1
                y += slope * step
                if gap > 40:
                    break
            x += step

    start = int(round(sx * 1100 + ix))
    follow(start, sy * 1.0 + iy, 1, 0.6)
    follow(start - 1, sy * 1.0 + iy, -1, 0.6)

    centerline = []
    for t in T_POINTS:
        x = int(round(sx * t + ix))
        if t == X_LABELS[-1]:
            x -= 3  # 1250 K 处与右边框重合，取框内最后几列
        win = [centers[k] for k in range(x - WINDOW, x + WINDOW + 1) if k in centers]
        centerline.append({"T_K": t, "pixel_x": x, "pixel_y_center": round(float(np.median(win)), 1)})

    doc = {
        "calibration": {
            "render": {"pdf_page": PAGE_INDEX + 1, "zoom": ZOOM, "half": "right",
                       "dark_threshold": DARK, "grid_threshold": GRID_LIGHT},
            "x_grid": [[round(p, 1), v] for p, v in zip(vgrid, X_LABELS)],
            "y_grid": [[round(p, 1), v] for p, v in zip(hgrid, Y_LABELS)],
            "centerline": centerline,
        }
    }
    doc["points"] = recompute(doc)
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
