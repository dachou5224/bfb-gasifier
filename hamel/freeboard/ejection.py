"""床面抛出通量与抛出初速分布。

平均泡间距用式 (3.57)。式 (3.58) 是 Baron 的轨道修正，还要 h_WS 和 A_0。
冻结稿写明正文没有在两式之间选定，ADR 0016 也没有选。任务卡给定的床面气泡量
是 ε_b、u_b、d_b，不含这两项，因此抛出通量按式 (3.57) 计算。
"""

from __future__ import annotations

import math

import numpy as np

from hamel.freeboard._common import (
    APPROACH_ANGLE_RAD,
    EJECTION_MEAN_FACTOR,
    EJECTION_SIGMA_FACTOR,
    F_W,
    N_BUBBLES,
    SPACING_SIGMA,
    ZETA_W,
    FreeboardError,
    _finite,
    _positive,
)


def _normal_cdf(x: float, mean: float, sigma: float) -> float:
    return 0.5 * (1.0 + math.erf((x - mean) / (sigma * math.sqrt(2.0))))


def _normal_ppf(probability: float, mean: float, sigma: float) -> float:
    low = mean - 12.0 * sigma
    high = mean + 12.0 * sigma
    for _ in range(70):
        mid = 0.5 * (low + high)
        if _normal_cdf(mid, mean, sigma) < probability:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def spacing_ratio(bubble_fraction: float) -> float:
    """Δl_m / d_b。气泡在床截面上均匀分布。

    Source: Hamel Eq.3.57, S. 38
    """
    fraction = _positive("bubble_fraction", bubble_fraction)
    if fraction >= 1.0:
        raise FreeboardError("ε_b 必须小于 1")
    return (6.0 * fraction / math.pi) ** (-1.0 / 3.0)


def approach_factor() -> float:
    """C = 1 − cos ξ，ξ = 48°。解释性选择：ADR 0016 第 1、2 条。

    Source: Hamel Eq.3.61, S. 39
    """
    return 1.0 - math.cos(APPROACH_ANGLE_RAD)


def pair_coalescence_probability(bubble_fraction: float) -> float:
    """两泡聚并概率 P″ = C · P′。P′ 是式 (3.55) 在 Δl/d_b < 1 上的积分。

    σ = 0,4 是 Δl/d_b 的标准差。ADR 0016 第 3 条。
    Source: Hamel Eq.3.55, S. 37
    Source: Hamel Eq.3.56, S. 38
    Source: Hamel Eq.3.62, S. 39
    """
    ratio = spacing_ratio(bubble_fraction)
    probability = _normal_cdf(1.0, ratio, SPACING_SIGMA)
    return approach_factor() * probability


def coalescence_weight(bubble_fraction: float) -> float:
    """Σ_{k=2}^{n_B} (k−1) P″_k。P″_k = (1−P″)^{k−1} P″，n_B = 5。

    解释性选择：ADR 0016 第 4 条。尾部 k > 5 不重新归一。
    Source: Hamel Eq.3.63, S. 39
    Source: Hamel Eq.3.64, S. 39
    """
    probability = pair_coalescence_probability(bubble_fraction)
    if probability < 0.0 or probability > 1.0:
        raise FreeboardError("P″ 必须在 [0, 1] 内")
    total = 0.0
    for count in range(2, N_BUBBLES + 1):
        term = (1.0 - probability) ** (count - 1) * probability
        total += (count - 1) * term
    return total


def ejection_flux(
    bubble_fraction: float,
    bubble_velocity_m_s: float,
    bubble_diameter_m: float,
    particle_density: float,
    suspension_voidage: float,
    mean_diameter_m: float,
) -> float:
    """床面抛出质量通量 F₀，kg/(m²·s)。空隙率照式面用 (1−ε_d)。

    ζ_w = 0,4，f_w = 0,25。解释性选择：ADR 0016 第 5 条。
    Source: Hamel Eq.3.66, S. 40
    Source: Hamel Eq.3.67, S. 40
    """
    fraction = _positive("bubble_fraction", bubble_fraction)
    if fraction >= 1.0:
        raise FreeboardError("ε_b 必须小于 1")
    velocity = _positive("bubble_velocity_m_s", bubble_velocity_m_s)
    diameter = _positive("bubble_diameter_m", bubble_diameter_m)
    density = _positive("particle_density", particle_density)
    voidage = _finite("suspension_voidage", suspension_voidage)
    if not 0.0 < voidage < 1.0:
        raise FreeboardError("ε_d 必须在 (0, 1) 内")
    mean_diameter = _positive("mean_diameter_m", mean_diameter_m)
    roof = 3.0 * mean_diameter / diameter
    weight = coalescence_weight(fraction)
    bracket = roof + (ZETA_W * F_W - roof) * weight
    return fraction * velocity * density * (1.0 - voidage) * bracket


def ejection_mass_flow(
    bubble_fraction: float,
    bubble_velocity_m_s: float,
    bubble_diameter_m: float,
    particle_density: float,
    suspension_voidage: float,
    mean_diameter_m: float,
    bed_area_m2: float,
) -> float:
    """床面抛出总质量流 F₀·A，kg/s。按粒级的分配不在这里做。

    Source: Hamel Eq.3.67, S. 40
    """
    area = _positive("bed_area_m2", bed_area_m2)
    return ejection_flux(
        bubble_fraction,
        bubble_velocity_m_s,
        bubble_diameter_m,
        particle_density,
        suspension_voidage,
        mean_diameter_m,
    ) * area


def ejection_velocity_samples(bubble_velocity_m_s: float, count: int) -> tuple[np.ndarray, np.ndarray, float]:
    """正初速样本。返回 (u_P,0, dF, 直接回落份额)。

    u_P,0 ~ N(1,53 u_b, (0,6 u_b)²)。正半轴按等概率切成 count 份，
    代表点取该份的概率中位。u_P,0 ≤ 0 的份额用分布函数一次算完，计入床面回落，
    不发射。各发射份额与直接回落份额之和为 1。解释性选择：ADR 0016 第 6 条。
    Source: Hamel Eq.3.72, S. 41
    Source: Hamel Eq.3.73, S. 42
    """
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise FreeboardError("样本数必须是正整数")
    velocity = _positive("bubble_velocity_m_s", bubble_velocity_m_s)
    mean = EJECTION_MEAN_FACTOR * velocity
    sigma = EJECTION_SIGMA_FACTOR * velocity
    direct = _normal_cdf(0.0, mean, sigma)
    launched = 1.0 - direct
    speeds = np.empty(count)
    shares = np.full(count, launched / count)
    for index in range(count):
        left = direct + launched * index / count
        right = direct + launched * (index + 1) / count
        speeds[index] = _normal_ppf(0.5 * (left + right), mean, sigma)
        if speeds[index] <= 0.0:
            raise FreeboardError("正半轴样本的代表初速不是正数")
    return speeds, shares, direct
