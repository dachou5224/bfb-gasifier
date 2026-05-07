from __future__ import annotations


def effective_rate_multiplier(requested_multiplier: float, *, lower: float = 0.0, upper: float = 1.0e6) -> float:
    """统一动力学速率倍率裁剪，避免极端输入。"""
    return float(min(max(float(requested_multiplier), float(lower)), float(upper)))
