from __future__ import annotations

import math
import numpy as np


def compound_heat_loss_fraction(total_frac: float, weight: float) -> float:
    """将总热损失比例映射到某一轴向权重分段。"""
    total = float(np.clip(total_frac, 0.0, 1.0))
    seg_weight = float(max(weight, 0.0))
    if total <= 0.0 or seg_weight <= 0.0:
        return 0.0
    return float(1.0 - math.pow(max(1.0 - total, 0.0), seg_weight))
