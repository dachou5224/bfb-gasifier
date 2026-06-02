from __future__ import annotations


def estimate_drag_coefficient(re_s: float) -> float:
    """简化阻力系数关联式（占位，用于统一接口）。"""
    re = max(float(re_s), 1e-12)
    if re < 1.0:
        return 24.0 / re
    if re < 1000.0:
        return 24.0 / re * (1.0 + 0.15 * re**0.687)
    return 0.44
