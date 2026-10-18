#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bild 3.3 与 Bild 3.5（PDF 26）五条 Berechnung 实线的数字化与复算。

两种用法：
- 默认（不读 PDF、不读页面图像）：用 data/digitized/fig_3_3.json 与 fig_3_5.json
  里记录的图框像素和每个取点的曲线中心像素换算纵坐标，并与文件中的 series 对照。
  测试用这一模式。
- ``--trace <页面 PNG 或 Hamel PDF>``：沿实线中心追踪。加上 ``--write`` 时把 JSON
  和叠图写回仓库。页面图像与 PDF 只在本地，不入库（CHARTER §11）。

只读图注 „Berechnung“ 的五条实线，不读测量符号。
压力点按 ADR 0013：0,2、0,5、1,0、1,5、2,0、2,5 MPa。
Bild 3.3 自上而下为 0,1、0,2、0,3、0,4、0,6 m/s（图注引线）。
Bild 3.5 自上而下为 0,6、0,4、0,3、0,2、0,1 m/s（曲线旁的标注）。
两条以上的线落在同一条没有亮隙的墨段里时，这些点记 null，并写明原因。
"""

from __future__ import annotations

import json
import math
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PATH_33 = ROOT / "data" / "digitized" / "fig_3_3.json"
PATH_35 = ROOT / "data" / "digitized" / "fig_3_5.json"
OVERLAY_DIR = ROOT / "data" / "digitized" / "overlays"

PRESSURES = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5]
# 图中自上而下的 u0−umf（m/s）。
ORDER_33 = [0.1, 0.2, 0.3, 0.4, 0.6]
ORDER_35 = [0.6, 0.4, 0.3, 0.2, 0.1]
PDF_PAGE_INDEX = 25  # PDF 26，零起
DARK = 45
DIGITS = 3
# 实线宽约 6 px。竖列切开一条斜线时可以到约 16 px，仍是一条线。


def round_half_up(value: float, ndigits: int) -> float:
    quantum = Decimal("0." + "0" * (ndigits - 1) + "1")
    return float(Decimal(str(value)).quantize(quantum, rounding=ROUND_HALF_UP))


def _finite_positive(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value > 0


def _x_edge(frame: dict, y: float, side: str) -> float:
    """竖轴随高度略斜，左右端的 x 按 y 在顶边与底边之间线性插值。"""
    top = frame[f"x_{side}"]
    bot = frame.get(f"x_{side}_bottom", top)
    y0 = frame["y_top_left"]
    y1 = frame["y_bottom_left"]
    return top + (y - y0) / (y1 - y0) * (bot - top)


def _p_of(frame: dict, x: float, y: float) -> float:
    left = _x_edge(frame, y, "left")
    right = _x_edge(frame, y, "right")
    return frame["p_left"] + (x - left) / (right - left) * (frame["p_right"] - frame["p_left"])


def _x_of(frame: dict, p: float, y: float) -> float:
    left = _x_edge(frame, y, "left")
    right = _x_edge(frame, y, "right")
    return left + (p - frame["p_left"]) / (frame["p_right"] - frame["p_left"]) * (right - left)


def _y_edge(frame: dict, x: float, which: str) -> float:
    t = (x - frame["x_left"]) / (frame["x_right"] - frame["x_left"])
    return frame[f"y_{which}_left"] + t * (frame[f"y_{which}_right"] - frame[f"y_{which}_left"])


def _value_of(frame: dict, x: float, y: float) -> float:
    top = _y_edge(frame, x, "top")
    bot = _y_edge(frame, x, "bottom")
    return frame["y_bottom_value"] + (bot - y) / (bot - top) * (
        frame["y_top_value"] - frame["y_bottom_value"]
    )


def _fit_residual(pixels, labels) -> float:
    px = np.array(pixels, dtype=float)
    val = np.array(labels, dtype=float)
    slope, intercept = np.polyfit(px, val, 1)
    return float(np.max(np.abs(val - (slope * px + intercept))))


def recompute_doc(doc: dict) -> list:
    """由 calibration 复算每个系列的 [P, 值或 null]。"""
    frame = doc["calibration"]["frame"]
    out = []
    for series in doc["calibration"]["series"]:
        points = []
        for sample in series["samples"]:
            y = sample["pixel_y_center"]
            if y is None:
                points.append([sample["P_MPa"], None])
            else:
                points.append([
                    sample["P_MPa"],
                    round_half_up(_value_of(frame, sample["pixel_x"], y), DIGITS),
                ])
        out.append({"u0_minus_umf_m_s": series["u0_minus_umf_m_s"], "points": points})
    return out


def check(doc: dict) -> list:
    problems = []
    frame = doc["calibration"]["frame"]
    if not frame["x_left"] < frame["x_right"]:
        problems.append("图框左右像素顺序不对")
    if not frame["y_top_left"] < frame["y_bottom_left"]:
        problems.append("图框上下像素顺序不对")
    if not frame["y_top_value"] > frame["y_bottom_value"]:
        problems.append("纵轴数值上大下小")
    order = ORDER_33 if doc["figure"] == "3.3" else ORDER_35
    labels = [s["u0_minus_umf_m_s"] for s in doc["calibration"]["series"]]
    if labels != order:
        problems.append(f"系列顺序 {labels} 不是自上而下的 {order}")
    if recompute_doc(doc) != doc["series"]:
        problems.append("series 与 calibration 复算不一致")
    err_key = "u_b_r" if doc["figure"] == "3.3" else "d_b"
    if not _finite_positive(doc.get("est_error", {}).get(err_key)):
        problems.append(f"est_error.{err_key} 须为正的有限数值")
    if not _finite_positive(doc.get("est_error", {}).get("P_MPa")):
        problems.append("est_error.P_MPa 须为正的有限数值")
    y_pixels = [row[0] for row in doc["calibration"]["y_grid"]]
    y_labels = [row[1] for row in doc["calibration"]["y_grid"]]
    if y_pixels != sorted(y_pixels) or y_labels != sorted(y_labels, reverse=True):
        problems.append("y_grid 须上大下小")
    elif _fit_residual(y_pixels, y_labels) > 0.015 * (y_labels[0] - y_labels[-1]):
        problems.append("y_grid 偏离直线")
    x_pixels = [row[0] for row in doc["calibration"]["x_grid"]]
    x_labels = [row[1] for row in doc["calibration"]["x_grid"]]
    if x_pixels != sorted(x_pixels) or x_labels != sorted(x_labels):
        problems.append("x_grid 须从左到右增大")
    elif _fit_residual(x_pixels, x_labels) > 0.02:
        problems.append("x_grid 偏离直线")
    problems.extend(_grid_against_frame(doc))
    for series in doc["calibration"]["series"]:
        got = [s["P_MPa"] for s in series["samples"]]
        if got != PRESSURES:
            problems.append(f"{series['u0_minus_umf_m_s']} 的压力点不是 {PRESSURES}")
        for sample in series["samples"]:
            if sample["pixel_y_center"] is None and not sample.get("note"):
                problems.append(f"P={sample['P_MPa']} 的空读数没有说明原因")
            y_for_p = sample["pixel_y_center"]
            if y_for_p is None:
                y_for_p = 0.5 * (frame["y_bottom_left"] + frame["y_bottom_right"])
            p_back = _p_of(frame, sample["pixel_x"], y_for_p)
            if abs(p_back - sample["P_MPa"]) > 0.03:
                problems.append(
                    f"P={sample['P_MPa']} 的 pixel_x {sample['pixel_x']} 换算为 {p_back:.3f} MPa"
                )
    at_1 = {}
    for series in doc["series"]:
        at_1[series["u0_minus_umf_m_s"]] = dict(series["points"]).get(1.0)
    vals = [at_1.get(v) for v in (0.1, 0.2, 0.3, 0.4, 0.6)]
    if any(v is None for v in vals):
        problems.append("1,0 MPa 处有空读数")
    elif doc["figure"] == "3.3" and not all(a > b for a, b in zip(vals, vals[1:])):
        problems.append("Bild 3.3 在 1,0 MPa 未呈现 0,1 在上、0,6 在下")
    elif doc["figure"] == "3.5" and not all(a < b for a, b in zip(vals, vals[1:])):
        problems.append("Bild 3.5 在 1,0 MPa 未呈现 0,6 在上、0,1 在下")
    return problems


# 刻度线本身有几像素宽，带心又是在一段宽度上取的，和单条参考线可以差这么多。
GRID_PX = 6


def _grid_against_frame(doc: dict) -> list:
    """用登记的参考行、列，把刻度像素经图框换算后与登记刻度比较。

    只检查网格自身的直线拟合时，图框比例错了、再按错的图框重算 series，检查仍会通过。
    """
    frame = doc["calibration"]["frame"]
    ref = doc["calibration"].get("grid_reference") or {}
    if "y_grid_x" not in ref or "x_grid_y" not in ref:
        return ["grid_reference 未标明横刻度所在的参考列和竖刻度所在的参考行"]
    problems = []
    y_px = abs(frame["y_bottom_left"] - frame["y_top_left"])
    y_val = abs(frame["y_top_value"] - frame["y_bottom_value"])
    y_tol = GRID_PX / y_px * y_val
    x_ref = ref["y_grid_x"]
    for pixel_y, label in doc["calibration"]["y_grid"]:
        got = _value_of(frame, x_ref, pixel_y)
        if abs(got - label) > y_tol:
            problems.append(
                f"y_grid 在参考列 x={x_ref} 处把 y={pixel_y} 换算为 {got:.4f}，登记刻度是 {label}"
            )
    x_px = abs(frame["x_right"] - frame["x_left"])
    x_val = abs(frame["p_right"] - frame["p_left"])
    x_tol = GRID_PX / x_px * x_val
    y_ref = ref["x_grid_y"]
    for pixel_x, label in doc["calibration"]["x_grid"]:
        got = _p_of(frame, pixel_x, y_ref)
        if abs(got - label) > x_tol:
            problems.append(
                f"x_grid 在参考行 y={y_ref} 处把 x={pixel_x} 换算为 {got:.4f} MPa，登记刻度是 {label}"
            )
    return problems


def check_files() -> list:
    problems = []
    for path in (PATH_33, PATH_35):
        doc = json.loads(path.read_text(encoding="utf-8"))
        problems += [f"{path.name}：{p}" for p in check(doc)]
    return problems


# --- 追踪（只在 --trace 时用到页面） ---


def _load_gray(path: str) -> np.ndarray:
    src = Path(path)
    if src.suffix.lower() == ".pdf":
        import fitz

        page = fitz.open(src)[PDF_PAGE_INDEX]
        pix = page.get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), colorspace=fitz.csGRAY)
        return np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    from PIL import Image

    return np.array(Image.open(src).convert("L"))


def _darkest_col(img, x0, x1, y0, y1, thresh=40) -> int:
    best = None
    for x in range(x0, x1):
        col = img[y0:y1, x]
        score = int((col < thresh).sum())
        if best is None or score > best[0]:
            best = (score, x)
    return best[1]


def _border_y(img, x, y0, y1) -> int:
    return y0 + int(np.argmin(img[y0:y1, x]))


def _frame_33(img) -> dict:
    x_left = float(_darkest_col(img, 368, 392, 1250, 1800))
    x_right = float(_darkest_col(img, 1430, 1455, 1250, 1800))
    return {
        "x_left": x_left,
        "x_right": x_right,
        "x_left_bottom": float(_darkest_col(img, 368, 392, 1800, 1860)),
        "x_right_bottom": float(_darkest_col(img, 1430, 1455, 1800, 1860)),
        "y_top_left": float(np.median([_border_y(img, x, 1168, 1192) for x in (420, 500, 600)])),
        "y_top_right": float(np.median([_border_y(img, x, 1165, 1190) for x in (1200, 1300, 1380)])),
        "y_bottom_left": float(np.median([_border_y(img, x, 1848, 1872) for x in (420, 500, 600)])),
        "y_bottom_right": float(np.median([_border_y(img, x, 1845, 1870) for x in (1200, 1300, 1380)])),
        "p_left": 0.0,
        "p_right": 2.5,
        "y_top_value": 2.5,
        "y_bottom_value": 0.0,
    }


def _frame_35(img) -> dict:
    return {
        "x_left": float(_darkest_col(img, 2070, 2095, 1280, 1400)),
        "x_right": float(_darkest_col(img, 3138, 3160, 1280, 1400)),
        "x_left_bottom": float(_darkest_col(img, 2060, 2085, 1800, 1910)),
        "x_right_bottom": float(_darkest_col(img, 3125, 3150, 1800, 1910)),
        "y_top_left": float(_border_y(img, 2140, 1230, 1255)),
        "y_top_right": float(_border_y(img, 3080, 1240, 1265)),
        "y_bottom_left": float(_border_y(img, 2140, 1910, 1935)),
        "y_bottom_right": float(_border_y(img, 3080, 1925, 1950)),
        "p_left": 0.0,
        "p_right": 2.5,
        "y_top_value": 0.4,
        "y_bottom_value": 0.0,
    }


def _runs(column: np.ndarray, y0: int, thresh: int = DARK):
    runs = []
    start = None
    for i, dark in enumerate(column < thresh):
        if dark and start is None:
            start = i
        elif not dark and start is not None:
            runs.append((y0 + start, y0 + i - 1))
            start = None
    if start is not None:
        runs.append((y0 + start, y0 + len(column) - 1))
    return runs


def _center(img, x: int, ya: int, yb: int) -> float:
    ys = np.arange(ya, yb + 1)
    weight = (255 - img[ys, x]).astype(float)
    return float(np.average(ys, weights=np.maximum(weight, 1.0)))


def _band_center(img, x0, x1, y0, y1) -> float:
    """一条横网格或横边框在 y0..y1 内的灰度质心。"""
    weights = []
    ys = []
    for y in range(y0, y1):
        dark = int((img[y, x0:x1] < 80).sum())
        if dark:
            weights.append(dark)
            ys.append(y)
    return float(np.average(ys, weights=weights))


def _grid_reference(figure: str, frame: dict) -> dict:
    """横刻度在图宽中段取带心，竖刻度在固定高度取最暗列。检查用同一参考。"""
    if figure == "3.3":
        x0, x1 = frame["x_left"] + 40, frame["x_right"] - 80
        y_ref = 1525.0  # 竖网格搜索带 y=1300..1750 的中线
    else:
        x0, x1 = frame["x_left"] + 50, frame["x_right"] - 50
        y_ref = 1400.0  # 竖网格搜索带 y=1300..1500 的中线
    return {
        "y_grid_x": round(0.5 * (x0 + x1), 1),
        "x_grid_y": y_ref,
        "note": "横刻度像素是参考列上的带心；竖刻度像素是参考行上的竖线。复算检查按这两处经图框换算。",
    }


def _grids_33(img, frame) -> tuple[list, list]:
    x0, x1 = int(frame["x_left"]) + 40, int(frame["x_right"]) - 80
    y_rows = [
        (_band_center(img, x0, x1, 1172, 1186), 2.5),
        (_band_center(img, x0, x1, 1310, 1324), 2.0),
        (_band_center(img, x0, x1, 1444, 1456), 1.5),
        (_band_center(img, x0, x1, 1583, 1595), 1.0),
        (_band_center(img, x0, x1, 1719, 1731), 0.5),
        (_band_center(img, x0, x1, 1855, 1866), 0.0),
    ]
    y_ref0, y_ref1 = 1300, 1750
    x_rows = []
    for p in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5):
        guess = int(round(_x_of(frame, p, 1600)))
        x_rows.append((float(_darkest_col(img, guess - 8, guess + 10, y_ref0, y_ref1)), p))
    return [[round(px, 1), p] for px, p in x_rows], [[round(py, 1), v] for py, v in y_rows]


def _grids_35(img, frame) -> tuple[list, list]:
    x0, x1 = int(frame["x_left"]) + 50, int(frame["x_right"]) - 50
    # 0,2 的横网格在 y≈1582；0,4 与 0 是顶边和底边。0,1 与 0,3 的网格较淡，不拿来标定。
    y_rows = [
        (_band_center(img, x0, x1, 1236, 1252), 0.4),
        (_band_center(img, x0, x1, 1576, 1594), 0.2),
        (_band_center(img, x0, x1, 1918, 1938), 0.0),
    ]
    y_ref0, y_ref1 = 1300, 1500
    x_rows = []
    for p in (0.0, 2.5):
        guess = int(round(_x_of(frame, p, 1400)))
        lo, hi = (guess - 6, guess + 12) if p == 0 else (guess - 12, guess + 8)
        x_rows.append((float(_darkest_col(img, lo, hi, y_ref0, y_ref1)), p))
    return [[round(px, 1), p] for px, p in x_rows], [[round(py, 1), v] for py, v in y_rows]


def _track(img, frame, seeds, x0, x_lo, x_hi, window, max_w, slope_limit):
    """从 seeds 向两侧走。返回每条线的 {x: (y, 是否单独占一条暗段)}。

    同一条暗段里落入两条以上追踪线时，这一列记为分不开，y 留在上一次分开的位置。
    线距相对种子处的间距突变时，拒绝这一步，避免被斜向引线或标注吸走。
    """
    y_top = int(min(frame["y_top_left"], frame["y_top_right"])) + 8
    y_bot = int(max(frame["y_bottom_left"], frame["y_bottom_right"])) - 6
    series = [{x0: (float(y), True)} for y in seeds]
    seed_gaps = np.diff(sorted(seeds))

    def step(x_from, x_to):
        direction = 1 if x_to > x_from else -1
        x = x_from
        while x != x_to:
            nxt = x + direction
            preds = []
            for s in series:
                xs = sorted(xx for xx in s if (xx - x) * direction <= 0)
                ok = [(xx, s[xx][0]) for xx in xs if s[xx][1]][-12:]
                if len(ok) >= 4:
                    slope = (ok[-1][1] - ok[0][1]) / (ok[-1][0] - ok[0][0])
                else:
                    slope = 0.0
                preds.append(s[x][0] + float(np.clip(slope, -slope_limit, slope_limit)) * direction)
            if float((img[y_top:y_bot, nxt] < 40).mean()) > 0.35:
                for s in series:
                    s[nxt] = (s[x][0], False)
                x = nxt
                continue
            found = [(a, b) for a, b in _runs(img[y_top:y_bot, nxt], y_top) if 3 <= b - a + 1 <= 30]
            host = []
            for pred in preds:
                hit = None
                for r, (a, b) in enumerate(found):
                    if a - 2 <= pred <= b + 2:
                        hit = r
                        break
                host.append(hit)
            taken = set()
            chosen = []
            for i, s in enumerate(series):
                h = host[i]
                width = None if h is None else found[h][1] - found[h][0] + 1
                if h is not None and host.count(h) == 1 and width <= max_w and h not in taken:
                    chosen.append((_center(img, nxt, *found[h]), True, h))
                    taken.add(h)
                elif h is not None and host.count(h) > 1 and width <= 28:
                    chosen.append((s[x][0], False, h))
                    taken.add(h)
                else:
                    best = None
                    for r, (a, b) in enumerate(found):
                        if r in taken or b - a + 1 > max_w:
                            continue
                        cy = _center(img, nxt, a, b)
                        if abs(cy - preds[i]) <= window and (best is None or abs(cy - preds[i]) < best[0]):
                            best = (abs(cy - preds[i]), cy, r)
                    if best is None:
                        chosen.append((preds[i], False, None))
                    else:
                        chosen.append((best[1], True, best[2]))
                        taken.add(best[2])
            # 跟丢的线用两侧仍分开的邻线，按种子处的间距比例找回
            ys = [c[0] for c in chosen]
            for i in range(len(chosen)):
                if chosen[i][1]:
                    continue
                above = next((j for j in range(i - 1, -1, -1) if chosen[j][1]), None)
                below = next((j for j in range(i + 1, len(chosen)) if chosen[j][1]), None)
                if above is None or below is None or seeds[below] == seeds[above]:
                    continue
                frac = (seeds[i] - seeds[above]) / (seeds[below] - seeds[above])
                y_hat = ys[above] + frac * (ys[below] - ys[above])
                best = None
                for r, (a, b) in enumerate(found):
                    if r in taken or b - a + 1 > max_w:
                        continue
                    cy = _center(img, nxt, a, b)
                    if abs(cy - y_hat) <= 8 and (best is None or abs(cy - y_hat) < best[0]):
                        best = (abs(cy - y_hat), cy, r)
                if best is not None:
                    chosen[i] = (best[1], True, best[2])
                    taken.add(best[2])
                    ys[i] = best[1]
            gaps = np.diff([c[0] for c in chosen])
            # 只在五条都声称分开时检查间距，避免一条被吸到邻线上
            if all(c[1] for c in chosen) and len(gaps) == len(seed_gaps):
                if np.any(gaps < 8) or np.any(gaps > seed_gaps * 1.8 + 12) or np.any(gaps < seed_gaps * 0.45):
                    chosen = [(s[x][0], False, None) for s in series]
            for s, (y, ok, _h) in zip(series, chosen):
                s[nxt] = (float(y), ok)
            x = nxt

    step(x0, x_lo)
    step(x0, x_hi)
    return series


def _expected_y(series: dict, x: int) -> float:
    ok = sorted((xx, yy) for xx, (yy, good) in series.items() if good)
    left = [p for p in ok if p[0] <= x]
    right = [p for p in ok if p[0] >= x]
    if left and right:
        (x0, y0), (x1, y1) = left[-1], right[0]
        if x1 == x0:
            return y0
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    if left:
        seg = [p for p in left if p[0] >= left[-1][0] - 80] or left[-2:]
        slope = 0.0 if len(seg) < 2 else (seg[-1][1] - seg[0][1]) / (seg[-1][0] - seg[0][0])
        return seg[-1][1] + slope * (x - seg[-1][0])
    seg = [p for p in right if p[0] <= right[0][0] + 80] or right[:2]
    slope = 0.0 if len(seg) < 2 else (seg[-1][1] - seg[0][1]) / (seg[-1][0] - seg[0][0])
    return seg[0][1] + slope * (x - seg[0][0])


def _consensus_runs(img, x, y0, y1, radius=2):
    """邻列里多数为墨的行才算墨。单列上 1～2 px 的缺口不把一条墨带拆开。"""
    dark = np.zeros(y1 - y0, dtype=int)
    for dx in range(-radius, radius + 1):
        xx = x + dx
        if 0 <= xx < img.shape[1]:
            dark += img[y0:y1, xx] < DARK
    mask = dark >= radius + 1
    runs = []
    start = None
    for i, on in enumerate(mask):
        if on and start is None:
            start = i
        elif not on and start is not None:
            if i - start >= 3:
                runs.append((y0 + start, y0 + i - 1))
            start = None
    if start is not None and len(mask) - start >= 3:
        runs.append((y0 + start, y0 + len(mask) - 1))
    return runs


def _sample(img, frame, tracked, pressures, search=10):
    """在每个压力上取同一列。两条以上追踪线落入同一墨段则这些点为 null。"""
    y_top = int(min(frame["y_top_left"], frame["y_top_right"])) + 4
    y_bot = int(max(frame["y_bottom_left"], frame["y_bottom_right"])) - 4
    n = len(tracked)
    columns = []
    for p in pressures:
        y_guess = float(np.median([_expected_y(s, int(round(_x_of(frame, p, 1600)))) for s in tracked]))
        nominal = int(round(_x_of(frame, p, y_guess)))
        if nominal >= int(frame["x_right"]) - 3:
            nominal = int(frame["x_right"]) - 4
        # 2,5 MPa 的名义列贴着右边框，优先向框内找
        offsets = list(range(0, -search - 1, -1)) + list(range(1, search + 1))
        if p >= 2.5:
            offsets = list(range(0, -16, -1))
        best = None
        for dx in offsets:
            x = nominal + dx
            if x <= frame["x_left"] + 4 or x >= frame["x_right"] - 1:
                continue
            if float((img[y_top:y_bot, x] < 40).mean()) > 0.35:
                continue
            found = _consensus_runs(img, x, y_top, y_bot)
            preds = [_expected_y(s, x) for s in tracked]
            host = []
            for pred in preds:
                hit = None
                for r, (a, b) in enumerate(found):
                    if a - 3 <= pred <= b + 3:
                        hit = r
                        break
                if hit is None:
                    # 预测落在宽墨带外侧几个像素时，仍算这条墨带，不把分不开的线记成“没找到”
                    near = None
                    for r, (a, b) in enumerate(found):
                        if b - a + 1 < 12:
                            continue
                        dist = 0 if a <= pred <= b else min(abs(pred - a), abs(pred - b))
                        if dist <= 24 and (near is None or dist < near[0]):
                            near = (dist, r)
                    if near is not None:
                        hit = near[1]
                host.append(hit)
            # 单独一条时只要线宽量级的墨段；几条挤在一起时允许更宽的共用墨段
            clean = True
            for h in host:
                if h is None:
                    clean = False
                    continue
                width = found[h][1] - found[h][0] + 1
                share = host.count(h)
                if not ((share == 1 and width <= 16) or (share > 1 and width <= 36)):
                    clean = False
            # 取离目标压力最近、且暗段宽度说得通的一列。不奖励“刚好裂开”的列。
            rank = (clean, -abs(dx))
            if best is None or rank > best[0]:
                best = (rank, x, found, host, preds)
        picked = None if best is None else best[1:]
        if picked is None:
            columns.append((nominal, None, [None] * n, [_expected_y(s, nominal) for s in tracked]))
        else:
            columns.append(picked)
    samples = [[] for _ in range(n)]
    for p, (x, found, host, preds) in zip(pressures, columns):
        for i in range(n):
            if found is None or host[i] is None:
                samples[i].append({
                    "P_MPa": p,
                    "pixel_x": float(x),
                    "pixel_y_center": None,
                    "note": "该列追踪位置没有落在实线上",
                })
                continue
            share = [j for j, h in enumerate(host) if h == host[i]]
            a, b = found[host[i]]
            width = b - a + 1
            if len(share) > 1 and width <= 36:
                samples[i].append({
                    "P_MPa": p,
                    "pixel_x": float(x),
                    "pixel_y_center": None,
                    "note": (
                        f"P={p} MPa 处有 {len(share)} 条 Berechnung 实线落在同一条"
                        f"宽 {width} px 的墨段里，中间没有亮隙，线距小于线宽，分不开"
                    ),
                })
            elif len(share) == 1 and width <= 16:
                samples[i].append({
                    "P_MPa": p,
                    "pixel_x": float(x),
                    "pixel_y_center": round(_center(img, x, a, b), 1),
                    "note": None,
                })
            else:
                # 没有合格列时也不把过宽的单独墨段写成正常读数
                samples[i].append({
                    "P_MPa": p,
                    "pixel_x": float(x),
                    "pixel_y_center": None,
                    "note": (
                        f"P={p} MPa 附近没有线宽合格的列："
                        f"该墨段宽 {width} px，单独一条不得超过 16 px，共用不得超过 36 px"
                    ),
                })
    return samples


def _seed_33(img, frame):
    """1,5 MPa 处五条线已分开，自上而下即 0,1…0,6。"""
    x = int(round(_x_of(frame, 1.5, frame["y_bottom_left"])))
    centers = [
        _center(img, x, a, b)
        for a, b in _runs(img[1760:1840, x], 1760)
        if 4 <= b - a + 1 <= 10
    ]
    if len(centers) < 5:
        raise RuntimeError(f"Bild 3.3 在 x={x} 只找到 {len(centers)} 条线")
    return centers[:5], x


def _seed_35(img, frame):
    """约 1,7 MPa 处五条线分开且远离标注，自上而下即 0,6…0,1。"""
    x = 2800
    y_top = int(_y_edge(frame, x, "top")) + 12
    y_bot = int(_y_edge(frame, x, "bottom")) - 20
    runs = []
    for a, b in _runs(img[y_top:y_bot, x], y_top):
        if 4 <= b - a + 1 <= 10:
            c = _center(img, x, a, b)
            if 0.07 <= _value_of(frame, x, c) <= 0.36:
                runs.append(c)
    if len(runs) < 5:
        raise RuntimeError(f"Bild 3.5 在 x={x} 只找到 {len(runs)} 条线")
    return runs[:5], x


def _doc(figure, frame, order, samples, x_grid, y_grid, quantity, unit, error, caption, side, page, evidence):
    series_cal = []
    series = []
    for label, sample in zip(order, samples):
        points = []
        for item in sample:
            if item["pixel_y_center"] is None:
                points.append([item["P_MPa"], None])
            else:
                points.append([
                    item["P_MPa"],
                    round_half_up(_value_of(frame, item["pixel_x"], item["pixel_y_center"]), DIGITS),
                ])
        series_cal.append({"u0_minus_umf_m_s": label, "samples": sample})
        series.append({"u0_minus_umf_m_s": label, "points": points})
    return {
        "figure": figure,
        "pdf_page": 26,
        "printed_page": page,
        "side": side,
        "caption_de": caption,
        "read_date": "2026-10-17",
        "tool": (
            "refs/hamel/pages/026.png（300 dpi）。图框取轴边墨迹最暗的行列；"
            "Berechnung 实线从五条已经分开的一列向两侧追踪，取点列取该暗段的灰度质心。"
            "同一暗段里有两条以上追踪线时记 null。"
            "复算：python3 tools/digitize_fig_3_5.py；重新追踪：同命令 --trace <页面 PNG 或 PDF> --write。"
        ),
        "units": {"P_MPa": "MPa", quantity: unit},
        "est_error": {"P_MPa": 0.01, quantity: error},
        "est_error_note": (
            "压力按图框左右端 0 与 2,5 MPa 线性换算，轴端定位约 ±2 px（约 0,005 MPa），登记 0,01 MPa。"
            f"{quantity} 的线宽约 6 px，质心约 ±2 px，再加图框倾斜的标定误差，分开的点用这一误差。"
            "分不开的点是 null，不另给数值。"
        ),
        "quantity_definition": {"P_MPa": "横轴 Druck", quantity: evidence},
        "points_format": f"[P_MPa, {quantity}]",
        "conditions": {
            "T_C": 25,
            "d_P_mm": 0.7,
            "rho_P_kg_m3": 2600,
            "h_m": 0.55,
            "source": "图中标注",
        },
        "series_order_from_top": order,
        "series": series,
        "use": "ADR 0013：在首次运行前用 Bild 3.5（主）与 Bild 3.3（分不开时）判别式 (3.43) 与 α_b。",
        "calibration": {
            "render": {"pdf_page": 26, "source": "refs/hamel/pages/026.png", "dpi": 300, "dark_threshold": DARK},
            "frame": frame,
            "grid_reference": _grid_reference(figure, frame),
            "x_grid": x_grid,
            "y_grid": y_grid,
            "series": series_cal,
        },
    }


def trace(path: str):
    img = _load_gray(path)
    frame33 = _frame_33(img)
    frame35 = _frame_35(img)
    seed33, x33 = _seed_33(img, frame33)
    seed35, x35 = _seed_35(img, frame35)
    tr33 = _track(img, frame33, seed33, x33, 450, int(frame33["x_right"]) - 4, window=8, max_w=16, slope_limit=4)
    tr35 = _track(
        img, frame35, seed35, x35,
        int(frame35["x_left"]) + 20, int(frame35["x_right"]) - 5,
        window=9, max_w=11, slope_limit=3,
    )
    s33 = _sample(img, frame33, tr33, PRESSURES)
    s35 = _sample(img, frame35, tr35, PRESSURES, search=14)
    xg33, yg33 = _grids_33(img, frame33)
    xg35, yg35 = _grids_35(img, frame35)
    doc33 = _doc(
        "3.3", frame33, ORDER_33, s33, xg33, yg33, "u_b_r", "m/s", 0.02,
        "Bild 3.3: Geschwindigkeit der Blasendurchströmung in Abhängigkeit des Drucks. "
        "Vergleich von Messwerten nach Olowson und Almstedt (1990) mit den nach Gleichung (3.44) berechneten Ergebnissen.",
        "left", "32",
        "纵轴 Gasgeschwindigkeit u_b,r relativ zur Blase。只读 Berechnung 实线。"
        "引线：标注 0,1 m/s 指向最上方曲线，0,6 m/s 指向最下方曲线；自上而下为 0,1、0,2、0,3、0,4、0,6 m/s。",
    )
    doc35 = _doc(
        "3.5", frame35, ORDER_35, s35, xg35, yg35, "d_b", "m", 0.005,
        "Bild 3.5: Blasendurchmesser in Abhängigkeit des Drucks. "
        "Vergleich von Messwerten nach Olowson und Almstedt (1990) mit den nach Gleichung (3.14) bis (3.44) berechneten Ergebnissen.",
        "right", "33",
        "纵轴 Blasendurchmesser d_b。只读 Berechnung 实线。"
        "曲线旁标注 0,6 m/s 在最上方，0,1 m/s 在最下方；自上而下为 0,6、0,4、0,3、0,2、0,1 m/s。",
    )
    return doc33, doc35, {"33": tr33, "35": tr35, "img": img}


def write_overlay(img, doc, tracked, dest: Path):
    from PIL import Image, ImageDraw, ImageFont

    rgb = np.stack([img, img, img], axis=-1)
    colors = {
        0.1: (210, 30, 30),
        0.2: (20, 140, 40),
        0.3: (20, 50, 200),
        0.4: (200, 110, 0),
        0.6: (140, 0, 170),
    }
    frame = doc["calibration"]["frame"]
    for series, tr in zip(doc["calibration"]["series"], tracked):
        color = colors[series["u0_minus_umf_m_s"]]
        for x, (y, ok) in tr.items():
            if not ok:
                continue
            yy, xx = int(round(y)), int(x)
            if 1 <= yy < rgb.shape[0] - 1 and 0 <= xx < rgb.shape[1]:
                rgb[yy - 1:yy + 2, xx] = color
    pad = 36
    x0 = max(0, int(frame["x_left"]) - 80)
    x1 = min(rgb.shape[1], int(max(frame["x_right"], frame.get("x_right_bottom", frame["x_right"]))) + 30)
    y0 = max(0, int(min(frame["y_top_left"], frame["y_top_right"])) - pad)
    y1 = min(rgb.shape[0], int(max(frame["y_bottom_left"], frame["y_bottom_right"])) + pad)
    crop = Image.fromarray(rgb[y0:y1, x0:x1])
    draw = ImageDraw.Draw(crop)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 20)
    except OSError:
        font = ImageFont.load_default()
    # 标签放在线距还比较开的压力上，并向右错开，避免叠在一起
    p_lab = 0.45 if doc["figure"] == "3.3" else 1.25
    y_ref = 0.5 * (frame["y_top_left"] + frame["y_bottom_left"])
    x_lab = _x_of(frame, p_lab, y_ref)
    for i, (series, tr) in enumerate(zip(doc["calibration"]["series"], tracked)):
        ys = [y for x, (y, ok) in tr.items() if ok and abs(x - x_lab) <= 4]
        if not ys:
            continue
        y = float(np.median(ys))
        label = f"{series['u0_minus_umf_m_s']:.1f}".replace(".", ",") + " m/s"
        tx = x_lab - x0 + 36 + i * 62
        ty = y - y0 - 11
        color = colors[series["u0_minus_umf_m_s"]]
        draw.line([(x_lab - x0, y - y0), (tx - 4, y - y0)], fill=color, width=2)
        draw.text((tx, ty), label, fill=color, font=font)
    dest.parent.mkdir(parents=True, exist_ok=True)
    crop.save(dest)


def _print_doc(doc):
    print(doc["figure"], "自上而下", doc["series_order_from_top"])
    for series in doc["series"]:
        print(" ", series["u0_minus_umf_m_s"], series["points"])


def main(argv: list) -> int:
    if argv and argv[0] == "--trace":
        if len(argv) < 2:
            print("用法：--trace <页面 PNG 或 PDF> [--write]", file=sys.stderr)
            return 2
        doc33, doc35, extra = trace(argv[1])
        _print_doc(doc33)
        _print_doc(doc35)
        if "--write" in argv:
            PATH_33.write_text(json.dumps(doc33, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            PATH_35.write_text(json.dumps(doc35, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            write_overlay(extra["img"], doc33, extra["33"], OVERLAY_DIR / "fig_3_3.png")
            write_overlay(extra["img"], doc35, extra["35"], OVERLAY_DIR / "fig_3_5.png")
            problems = check_files()
            for item in problems:
                print("问题：" + item)
            print("已写入" if not problems else f"已写入，但有 {len(problems)} 个问题")
            return 1 if problems else 0
        return 0
    problems = check_files()
    for item in problems:
        print("问题：" + item)
    print("通过" if not problems else f"{len(problems)} 个问题")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
