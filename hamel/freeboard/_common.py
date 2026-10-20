"""自由空间模块共用的异常与重力常数。"""

from __future__ import annotations

import math

import numpy as np

# 与 ADR 0014 的床料质量换算、hamel.hydrodynamics 相同。
G = 9.81

# 式 (3.67) 的尾涡份额。解释性选择：ADR 0016 第 5 条。
ZETA_W = 0.4
F_W = 0.25
N_BUBBLES = 5
# 式 (3.55) 的 σ，Δl/d_b 的无量纲标准差。ADR 0016 第 3 条。
SPACING_SIGMA = 0.4
# Briens–Baron 半角。解释性选择：ADR 0016 第 2 条。
APPROACH_ANGLE_RAD = math.radians(48.0)
# 式 (3.73)。解释性选择：ADR 0016 第 6 条。
EJECTION_MEAN_FACTOR = 1.53
EJECTION_SIGMA_FACTOR = 0.6


class FreeboardError(ValueError):
    """原方程下不能继续：输入不成立，或到时间上限仍未确定归宿。"""


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or isinstance(value, np.ndarray):
        raise FreeboardError(f"{name} 必须是有限标量")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise FreeboardError(f"{name} 必须是有限标量") from exc
    if not math.isfinite(number):
        raise FreeboardError(f"{name} 不是有限值")
    return number


def _positive(name: str, value: object) -> float:
    number = _finite(name, value)
    if number <= 0.0:
        raise FreeboardError(f"{name} 必须大于 0")
    return number


