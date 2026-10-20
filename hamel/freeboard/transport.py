"""自由空间气速。悬停判定和输运系数不在本模块。"""

from __future__ import annotations

import math

import numpy as np

from hamel.freeboard._common import FreeboardError, _positive


def cross_section_m2(cell: int, height_m: float, edges_m, diameters_m) -> float:
    """格内直径线性插值，截面为圆。

    Source: Hamel Eq.3.74, S. 43
    """
    edges = np.asarray(edges_m, dtype=float)
    diameters = np.asarray(diameters_m, dtype=float)
    span = float(edges[cell + 1] - edges[cell])
    if span <= 0.0:
        raise FreeboardError("格边界必须严格递增")
    clamped = min(max(float(height_m), float(edges[cell])), float(edges[cell + 1]))
    fraction = (clamped - float(edges[cell])) / span
    diameter = float(diameters[cell]) + fraction * float(diameters[cell + 1] - diameters[cell])
    if diameter <= 0.0:
        raise FreeboardError("直径必须为正")
    return math.pi * diameter * diameter / 4.0


def superficial_velocity(cell: int, height_m: float, edges_m, diameters_m, volume_flow_m3_s) -> float:
    """格内气体体积流除以当地截面。二次风由调用方写进该格的体积流。

    解释性选择：ADR 0016 第 7 条。
    Source: Hamel Eq.3.74, S. 43
    """
    flows = np.asarray(volume_flow_m3_s, dtype=float)
    flow = float(flows[cell])
    if flow < 0.0:
        raise FreeboardError("气体体积流不能为负")
    return flow / cross_section_m2(cell, height_m, edges_m, diameters_m)


def ghost_bubble_velocity(
    cell: int,
    height_m: float,
    edges_m,
    diameters_m,
    volume_flow_m3_s,
    bubble_velocity_m_s: float,
    bubble_diameter_m: float,
) -> float:
    """u_gb(h) = u_0(h) + (u_b,WS − u_0(h)) exp(−h / (3 d_b,WS))。

    h 从床面起算。解释性选择：ADR 0016 第 7 条，u_g = u_gb(h)。
    Source: Hamel Eq.3.74, S. 43
    Source: Hamel Eq.3.75, S. 43
    Source: Hamel Eq.3.76, S. 43
    """
    approach = superficial_velocity(cell, height_m, edges_m, diameters_m, volume_flow_m3_s)
    bubble_velocity = _positive("bubble_velocity_m_s", bubble_velocity_m_s)
    bubble_diameter = _positive("bubble_diameter_m", bubble_diameter_m)
    beta = 1.0 / (3.0 * bubble_diameter)
    return approach + (bubble_velocity - approach) * math.exp(-beta * float(height_m))


