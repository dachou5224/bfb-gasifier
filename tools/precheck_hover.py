#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""自由空间内部悬停预检。不进入 hamel/，不改闭合。

气速取 u_gb(h)（式 3.74），u_g = u_gb（ADR 0016 第 7 条）。终端沉降速度由式 (3.82)
在加速度为零时与式 (3.85) 联立，z 取 Re_∞ → ∞（ADR 0016 第 8 条）。

两类参考轨迹分开报告「会 / 不会 / 未检验」：
- 床面抛出样本：式 (3.72) 的起始速度，σ 按 ADR 0016 第 6 条；
- 当地生成：床面样本到不了的 cell，从中点以 u_g(h_i) − v_t 起步（ADR 0016 第 9、10 条）。

只有两类都检查过、都不趋近稳定点，综合结论才是「不会」。任一为「会」则综合为「会」。
没有稳定高度时，两类都不可能趋近，综合为「不会」，并写明已检查。

趋近判据（预检用，不是悬停闭合方案）：
- 与某个稳定高度相距 ≤ 5 mm 且 |u_P| ≤ 0,01 m/s；
- 或积分满 40 s 时仍朝该高度移动，距离 < 5 cm 且 |u_P| < 0,05 m/s
  （停留时间发散，不把截断当成落回或带出）；
- 落到床面、到达出口，或当地轨迹离开所在 cell：不趋近；
- 时间用尽且不满足上列趋近条件：未检验。
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from hamel.hydrodynamics.bed import (
    G,
    BedConditions,
    HydroSearch,
    HydrodynamicsError,
    bed_profile,
    minimum_fluidization_velocity,
    particle_density_from_bed_mass,
)

# 趋近判据。写在这里，测试直接用同一组数。
HEIGHT_TOL_M = 5.0e-3
SPEED_TOL_M_S = 1.0e-2
TIME_CAP_S = 40.0
APPROACH_HEIGHT_M = 5.0e-2
APPROACH_SPEED_M_S = 5.0e-2
STEP_CAP = 80_000

# 床面样本离散数。临时：预检固定 17 个等概率正初速，不是 SolverConfig。
N_EJECTION_SAMPLES = 17

# ε_mf、φ_S 的扫描端点与中点。临时区间，待 T28（ADR 0022 C6；任务卡）。
EPSILON_MF_GRID = (0.40, 0.475, 0.55)
PHI_S_GRID = (0.60, 0.75, 0.90)

# ADR 0020 登记的搜索参数。该 ADR 仍是待审查；T19 函数不设默认值，预检作为调用方传入。
HYDRO_SEARCH = HydroSearch(
    1.0e-12,
    1.0,
    80,
    1.0e-8,
    2.0,
    1.35,
    36,
    8000.0,
    60,
    1.0e-6,
)


@dataclass(frozen=True)
class Particle:
    """自由空间单颗粒。物性在一条轨迹内不变。"""

    diameter_m: float
    rho_p: float
    rho_g: float
    mu_g: float
    phi_s: float

    @property
    def stokes(self) -> float:
        """18 μ_g / (ρ_P d_P²)，式 (3.82) 中 24/Re 项化成的线性系数。"""
        return 18.0 * self.mu_g / (self.rho_p * self.diameter_m**2)

    @property
    def quadratic(self) -> float:
        """3 ρ_g z / (4 ρ_P d_P)。解释性选择：ADR 0016 第 8 项。"""
        return 3.0 * self.rho_g * drag_limit_z(self.phi_s) / (4.0 * self.rho_p * self.diameter_m)

    @property
    def g_eff(self) -> float:
        return (self.rho_p - self.rho_g) * G / self.rho_p


@dataclass
class Flight:
    outcome: str
    height_peak_m: float
    time_s: float


@dataclass
class FamilyVerdict:
    """一类轨迹的结论。status 只取 会、不会、未检验。"""

    status: str
    detail: str


@dataclass
class HoverAssessment:
    v_t_m_s: float
    stable_heights_m: list[float]
    bed: FamilyVerdict
    local: FamilyVerdict
    overall: str
    unreachable_cells: list[int] = field(default_factory=list)


def drag_limit_z(phi_s: float) -> float:
    """z = 73,69·exp(−5,0784·φ_S)。24/Re_∞ 在 Re_∞ → ∞ 时消失。

    解释性选择：ADR 0016 第 8 项。
    Source: Hamel Eq.3.86, S. 45
    """
    if not 0.0 < phi_s <= 1.0:
        raise ValueError("φ_S 必须在 (0, 1] 内")
    return 73.69 * math.exp(-5.0784 * phi_s)


def terminal_settling_velocity(particle: Particle) -> float:
    """静气中的终端沉降速度 v_t > 0，向上为正时平衡相对速度为 −v_t。

    令式 (3.82) 的加速度为零，c_W 用式 (3.85)。该阻力下 v_t 有二次方程的正根。
    Source: Hamel Eq.3.82, S. 44
    Source: Hamel Eq.3.85, S. 45
    Source: Hamel Eq.3.87, S. 45
    """
    stokes = particle.stokes
    quadratic = particle.quadratic
    weight = particle.g_eff
    if quadratic <= 0.0:
        return weight / stokes
    return (-stokes + math.sqrt(stokes**2 + 4.0 * quadratic * weight)) / (2.0 * quadratic)


def motion_residual(particle: Particle, relative_m_s: float) -> float:
    """式 (3.82) 在给定相对速度 u_r = u_P − u_g 下的加速度。静气沉降时 u_r = −v_t。"""
    stokes = particle.stokes
    quadratic = particle.quadratic
    return -stokes * relative_m_s - quadratic * relative_m_s * abs(relative_m_s) - particle.g_eff


def stable_heights(height_m, velocity_m_s, v_t: float) -> list[float]:
    """u_g = v_t 且沿高度穿过时气速在下降的位置。上升穿过的是不稳定点，不收。

    height_m、velocity_m_s 为同一条光滑段上的采样，高度严格递增。
    """
    height = np.asarray(height_m, dtype=float)
    velocity = np.asarray(velocity_m_s, dtype=float)
    if height.ndim != 1 or velocity.shape != height.shape or height.size < 2:
        raise ValueError("稳定高度需要至少两个采样点")
    gap = velocity - v_t
    found = []
    for index in range(height.size - 1):
        left = gap[index]
        right = gap[index + 1]
        if left == 0.0 and velocity[index + 1] < velocity[index]:
            found.append(float(height[index]))
            continue
        if left > 0.0 and right <= 0.0 and velocity[index + 1] < velocity[index]:
            span = left - right
            fraction = 0.0 if span == 0.0 else left / span
            found.append(float(height[index] + fraction * (height[index + 1] - height[index])))
    return found


def _captured(height: float, velocity: float, anchors: list[float]) -> bool:
    return any(abs(height - anchor) <= HEIGHT_TOL_M and abs(velocity) <= SPEED_TOL_M_S for anchor in anchors)


def _still_approaching(height: float, velocity: float, anchors: list[float], low: float, high: float) -> bool:
    for anchor in anchors:
        if not (low < anchor < high):
            continue
        toward = (anchor - height) * velocity > 0.0
        if toward and abs(height - anchor) < APPROACH_HEIGHT_M and abs(velocity) < APPROACH_SPEED_M_S:
            return True
    return False


def _advance(height: float, velocity: float, step: float, u_of_h, particle: Particle):
    gas = float(u_of_h(height))
    relative = velocity - gas
    stokes = particle.stokes
    explicit = -particle.quadratic * relative * abs(relative) - particle.g_eff
    if stokes * step < 1.0e-8:
        updated = velocity + step * (-stokes * relative + explicit)
    else:
        decay = math.exp(-stokes * step)
        updated = gas + (velocity - gas) * decay + explicit * (1.0 - decay) / stokes
    return height + 0.5 * (velocity + updated) * step, updated


def fly(
    height_0: float,
    velocity_0: float,
    u_of_h,
    particle: Particle,
    low: float,
    high: float,
    anchors: list[float],
) -> Flight:
    """积分 dh/dt = u_P 与式 (3.82)，直到离开 [low, high]、趋近稳定点或时间用尽。

    起点落在下边界且向上时，不把起点本身当成已经落回。
    Source: Hamel Eq.3.82, S. 44
    """
    height = height_0
    velocity = velocity_0
    elapsed = 0.0
    peak = height_0
    inside = [anchor for anchor in anchors if low < anchor < high]
    if _captured(height, velocity, inside):
        return Flight("hover", peak, elapsed)
    for _ in range(STEP_CAP):
        speed = max(abs(velocity), 1.0e-3)
        step = min(0.02, 1.0e-3 / speed)
        height, velocity = _advance(height, velocity, step, u_of_h, particle)
        elapsed += step
        if not math.isfinite(height) or not math.isfinite(velocity):
            return Flight("unresolved", peak, elapsed)
        peak = max(peak, height)
        if height < low:
            return Flight("left_down", peak, elapsed)
        if height > high:
            return Flight("left_up", peak, elapsed)
        if _captured(height, velocity, inside):
            return Flight("hover", peak, elapsed)
        if elapsed >= TIME_CAP_S:
            if _still_approaching(height, velocity, inside, low, high):
                return Flight("hover", peak, elapsed)
            return Flight("unresolved", peak, elapsed)
    return Flight("unresolved", peak, elapsed)


def _family(outcomes: list[str], empty_detail: str) -> FamilyVerdict:
    if not outcomes:
        return FamilyVerdict("不会", empty_detail)
    if any(item == "hover" for item in outcomes):
        return FamilyVerdict("会", "有轨迹趋近稳定高度")
    if any(item == "unresolved" for item in outcomes):
        return FamilyVerdict("未检验", "积分时间用尽，且没有样本进入趋近判据")
    return FamilyVerdict("不会", "已积分的样本都在有限时间内离开，没有趋近稳定高度")


def combine_verdicts(bed: str, local: str) -> str:
    """两类都是「不会」才综合为「不会」。任一为「会」则综合为「会」。"""
    if bed == "会" or local == "会":
        return "会"
    if bed == "不会" and local == "不会":
        return "不会"
    return "未检验"


def assess_hover(
    u_of_h,
    height_samples_m,
    velocity_samples_m_s,
    cell_edges_m,
    particle: Particle,
    samples: list[tuple[float, float]],
) -> HoverAssessment:
    """对一条气速剖面做两类轨迹预检。

    samples 为 (u_P,0, dF_s)。dF_s ≤ 0 的不发射。cell_edges_m 从床面向上。
    气速采样可以分段传入：这里用一条已拼好的采样；稳定高度只在采样的相邻点之间找，
    调用方不要把间断两侧连成一段。
    """
    v_t = terminal_settling_velocity(particle)
    anchors = stable_heights(height_samples_m, velocity_samples_m_s, v_t)
    edges = np.asarray(cell_edges_m, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("cell 边界必须严格递增")
    low = float(edges[0])
    high = float(edges[-1])
    if not anchors:
        detail = "已检查：这段自由空间里没有 du_g/dh < 0 的 u_g = v_t，两类轨迹都不可能趋近"
        bed = FamilyVerdict("不会", detail)
        local = FamilyVerdict("不会", detail)
        return HoverAssessment(v_t, [], bed, local, "不会", [])

    emitted = [(speed, share) for speed, share in samples if share > 0.0]
    peaks = []
    bed_outcomes = []
    for speed, _share in emitted:
        flown = fly(low, speed, u_of_h, particle, low, high, anchors)
        peaks.append(flown.height_peak_m)
        bed_outcomes.append(flown.outcome)
    bed = _family(bed_outcomes, "没有 dF_s > 0 的床面样本")

    unreachable = []
    local_outcomes = []
    for index in range(edges.size - 1):
        cell_low = float(edges[index])
        cell_high = float(edges[index + 1])
        entered = any(peak > cell_low + 1.0e-9 for peak in peaks)
        if entered:
            continue
        unreachable.append(index)
        midpoint = 0.5 * (cell_low + cell_high)
        start = float(u_of_h(midpoint)) - v_t
        flown = fly(midpoint, start, u_of_h, particle, cell_low, cell_high, anchors)
        local_outcomes.append(flown.outcome)
    if unreachable:
        local = _family(local_outcomes, "没有床面不可达的 cell")
    else:
        local = FamilyVerdict("不会", "床面样本到达了全部自由空间 cell，没有当地生成轨迹要积分")
    return HoverAssessment(v_t, anchors, bed, local, combine_verdicts(bed.status, local.status), unreachable)


def ejection_samples(bubble_velocity_m_s: float, count: int = N_EJECTION_SAMPLES) -> list[tuple[float, float]]:
    """床面抛出初速。均值 1,53 u_b，σ = 0,6 u_b，只发射 u_P,0 > 0 的部分，份额和为 1。

    解释性选择：ADR 0016 第 6 条。正半轴按等概率切成 count 份，代表点取区间中位初速。
    Source: Hamel Eq.3.72, S. 41
    Source: Hamel Eq.3.73, S. 42
    """
    if count < 1:
        raise ValueError("样本数至少为 1")
    mean = 1.53 * bubble_velocity_m_s
    sigma = 0.6 * bubble_velocity_m_s
    if sigma <= 0.0:
        raise ValueError("u_b 必须为正，否则式 (3.72) 的标准差没有定义")

    def cdf(speed: float) -> float:
        return 0.5 * (1.0 + math.erf((speed - mean) / (sigma * math.sqrt(2.0))))

    def inverse(probability: float) -> float:
        low = min(0.0, mean - 8.0 * sigma)
        high = mean + 8.0 * sigma
        for _ in range(80):
            mid = 0.5 * (low + high)
            if cdf(mid) < probability:
                low = mid
            else:
                high = mid
        return 0.5 * (low + high)

    floor = cdf(0.0)
    mass = 1.0 - floor
    if mass <= 0.0:
        return []
    samples = []
    share = 1.0 / count
    for index in range(count):
        left = floor + mass * index / count
        right = floor + mass * (index + 1) / count
        samples.append((inverse(0.5 * (left + right)), share))
    return samples


def assess_piecewise(u_of_h, pieces, cell_edges_m, particle: Particle, samples) -> HoverAssessment:
    """在各光滑段内找稳定高度，再沿整段自由空间积分。段与段之间不连线。"""
    v_t = terminal_settling_velocity(particle)
    anchors = []
    for start, stop in pieces:
        grid = np.linspace(start, stop, 500)
        speed = np.array([float(u_of_h(item)) for item in grid])
        anchors.extend(stable_heights(grid, speed, v_t))
    return _assess_with_anchors(u_of_h, anchors, cell_edges_m, particle, samples, v_t)


def _assess_with_anchors(u_of_h, anchors, cell_edges_m, particle, samples, v_t) -> HoverAssessment:
    edges = np.asarray(cell_edges_m, dtype=float)
    low = float(edges[0])
    high = float(edges[-1])
    if not anchors:
        detail = "已检查：这段自由空间里没有 du_g/dh < 0 的 u_g = v_t，两类轨迹都不可能趋近"
        bed = FamilyVerdict("不会", detail)
        local = FamilyVerdict("不会", detail)
        return HoverAssessment(v_t, [], bed, local, "不会", [])
    emitted = [(speed, share) for speed, share in samples if share > 0.0]
    peaks = []
    bed_outcomes = []
    for speed, _share in emitted:
        flown = fly(low, speed, u_of_h, particle, low, high, anchors)
        peaks.append(flown.height_peak_m)
        bed_outcomes.append(flown.outcome)
    bed = _family(bed_outcomes, "没有 dF_s > 0 的床面样本")
    unreachable = []
    local_outcomes = []
    for index in range(edges.size - 1):
        cell_low = float(edges[index])
        cell_high = float(edges[index + 1])
        if any(peak > cell_low + 1.0e-9 for peak in peaks):
            continue
        unreachable.append(index)
        midpoint = 0.5 * (cell_low + cell_high)
        start = float(u_of_h(midpoint)) - v_t
        flown = fly(midpoint, start, u_of_h, particle, cell_low, cell_high, anchors)
        local_outcomes.append(flown.outcome)
    if unreachable:
        local = _family(local_outcomes, "没有床面不可达的 cell")
    else:
        local = FamilyVerdict("不会", "床面样本到达了全部自由空间 cell，没有当地生成轨迹要积分")
    return HoverAssessment(v_t, anchors, bed, local, combine_verdicts(bed.status, local.status), unreachable)


# --- 工况 1 的估算输入。每一项在 INPUT_NOTES 里写出处或「临时」。 ---

BED_DIAMETER_M = 0.15  # Hamel §7.1.3，S. 132
FREEBOARD_DIAMETER_M = 0.25  # Hamel §7.1.3，S. 132
BED_HEIGHT_M = 1.2  # ADR 0014 A2
FREEBOARD_HEIGHT_M = 2.3  # KURKELA_1992 表 3，经 ADR 0005
SECONDARY_HEIGHT_M = 1.4875  # 床面 1,2 m + 自由空间第一格 0,2875 m，ADR 0014 B
CELL_LENGTH_M = 0.2875  # ADR 0014 B，自由空间 2,3 m 分 8 格
PRESSURE_PA = 500_000.0  # 表 7.2，S. 133，0,5 MPa
BED_TEMPERATURE_K = 819.0 + 273.15  # KURKELA_1992 表 5 工况 1 床温，仅作预检
BED_PRESSURE_DROP_PA = 2200.0  # KURKELA_1992 表 5，ADR 0014 A3
PRIMARY_AIR_MOL_S = 0.6419598861822472  # 表 7.2，ADR 0004
SECONDARY_AIR_MOL_S = 0.373030744673468  # 表 7.2，ADR 0004
STEAM_KG_S = 0.001888888888888889  # 表 7.2，6,8 kg/h
PURGE_N2_MOL_S = 0.07683689757393693  # KURKELA_1992 表 5，ADR 0006 暂定计入
N_ORIFICE = 139  # ADR 0019 为 139,5；T19 只接受整数。临时取整，见 INPUT_NOTES
H_A_M = 0.0  # ADR 0019 第 1 项
GAS_CONSTANT = 8.314462618  # 临时：摩尔气体常数，J/(mol·K)
AIR_MOLAR_MASS = 0.0289647  # 临时：干空气
WATER_MOLAR_MASS = 0.01801528  # 临时
N2_MOLAR_MASS = 0.0280134  # 临时
FINE_DIAMETER_M = 0.05e-3  # 临时：<0,1 mm 没有下限，床料代表粒径用筛孔上限的一半

INPUT_NOTES = (
    "床径 0,15 m、自由空间直径 0,25 m：Hamel §7.1.3，S. 132。",
    "床高 1,2 m：ADR 0014 A2。自由空间高 2,3 m：KURKELA_1992 表 3（ADR 0005）。",
    "二次风口 1,4875 m：ADR 0014 B，自由空间第一格上沿；以上计入二次风。",
    "床温 819 °C：KURKELA_1992 表 5 工况 1，只作预检。全床与自由空间气体体积都按这个温度，临时"
    "（二次风加入温度是 458,15 K，混合后的温度本预检不求）。",
    "压力 0,5 MPa：表 7.2。一次风、二次风、蒸汽：表 7.2，摩尔流量按 ADR 0004。"
    "吹扫氮按 ADR 0006 暂定计入。",
    "干空气摩尔质量 0,0289647 kg/mol、水 0,01801528 kg/mol、氮 0,0280134 kg/mol、"
    "R = 8,314462618 J/(mol·K)：临时。",
    "气体粘度用空气 Sutherland 公式（1,716·10⁻⁵ Pa·s，S = 110,4 K）：临时，T21 未合入。",
    "n_0：ADR 0019 为 139,5。T19 的式 (3.37) 只接受正整数，临时取 139。"
    "h_A = 0：ADR 0019 第 1 项。",
    "Geldart 取 B：临时，待 T28。原文没有印出分类边界。",
    "ε_mf ∈ {0,40, 0,475, 0,55}，φ_S ∈ {0,60, 0,75, 0,90}：临时区间，待 T28。"
    "ρ_P 用 T19 的 particle_density_from_bed_mass，床料质量按 ADR 0014 A3。",
    "床料代表粒径：有上下限的筛分级取几何平均；<0,1 mm 没有下限，临时取 0,05 mm。"
    "只用于反求 ρ_P 和 u_b,WS，不作为悬停粒级。",
    "悬停粒级用 kurkela1992_fig1.json 的筛孔上限与下限。<0,1 mm 的下限未检验。",
    "床区积分每格 6 段，与 T19 冒烟相同。搜索参数用 ADR 0020 登记值（该 ADR 待审查）。",
    "鬼泡坐标 h 从床面起算：u_gb(0) = u_b,WS。这是式 (3.76) 的读法，任务卡没有另给原点。",
    "床面样本 17 个等概率正初速：临时离散，不是求解器参数。",
)


def sutherland_air_viscosity(temperature_k: float) -> float:
    """空气动力粘度。临时：Sutherland，参考 273,15 K、1,716·10⁻⁵ Pa·s、S = 110,4 K。"""
    reference = 273.15
    return 1.716e-5 * (temperature_k / reference) ** 1.5 * (reference + 110.4) / (temperature_k + 110.4)


def case1_gas():
    """工况 1 预检用的气体状态。组成是床层气体（一次风 + 蒸汽 + 吹扫氮），自由空间不另算二次风稀释。"""
    steam_mol_s = STEAM_KG_S / WATER_MOLAR_MASS
    amounts = {
        "air": PRIMARY_AIR_MOL_S,
        "steam": steam_mol_s,
        "n2": PURGE_N2_MOL_S,
    }
    total = sum(amounts.values())
    mass = (
        amounts["air"] * AIR_MOLAR_MASS
        + amounts["steam"] * WATER_MOLAR_MASS
        + amounts["n2"] * N2_MOLAR_MASS
    )
    molar_mass = mass / total
    rho_g = PRESSURE_PA * molar_mass / (GAS_CONSTANT * BED_TEMPERATURE_K)
    mu_g = sutherland_air_viscosity(BED_TEMPERATURE_K)
    return {
        "primary_mol_s": total,
        "secondary_mol_s": SECONDARY_AIR_MOL_S,
        "molar_mass": molar_mass,
        "rho_g": rho_g,
        "mu_g": mu_g,
        "bed_area_m2": math.pi * (BED_DIAMETER_M / 2.0) ** 2,
        "freeboard_area_m2": math.pi * (FREEBOARD_DIAMETER_M / 2.0) ** 2,
    }


def superficial_velocity(height_above_distributor_m: float, gas: dict) -> float:
    """所在高度的表观气速。二次风口以上才加入二次风。锥段直径从 0,15 m 线性到 0,25 m。"""
    amount = gas["primary_mol_s"]
    if height_above_distributor_m > SECONDARY_HEIGHT_M:
        amount += gas["secondary_mol_s"]
    volume = amount * GAS_CONSTANT * BED_TEMPERATURE_K / PRESSURE_PA
    return volume / cross_section(height_above_distributor_m, gas)


def cross_section(height_above_distributor_m: float, gas: dict) -> float:
    if height_above_distributor_m <= BED_HEIGHT_M:
        return gas["bed_area_m2"]
    if height_above_distributor_m >= SECONDARY_HEIGHT_M:
        return gas["freeboard_area_m2"]
    fraction = (height_above_distributor_m - BED_HEIGHT_M) / CELL_LENGTH_M
    diameter = BED_DIAMETER_M + (FREEBOARD_DIAMETER_M - BED_DIAMETER_M) * fraction
    return math.pi * (diameter / 2.0) ** 2


def size_bounds():
    """筛孔上下限。最细一档没有下限。质量分数来自同一 JSON，只用于床料代表粒径。"""
    path = Path(__file__).resolve().parents[1] / "data" / "digitized" / "kurkela1992_fig1.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    fractions = [item / 100.0 for item in data["class_mass_fractions_pct"]]
    bounds_mm = [None, 0.1, 0.25, 0.5, 0.71, 1.0, 2.0, 3.15]
    classes = []
    for index, fraction in enumerate(fractions):
        classes.append((bounds_mm[index], bounds_mm[index + 1], fraction))
    return classes


def bed_representative_diameter_m(classes) -> float:
    total = 0.0
    moment = 0.0
    for low_mm, high_mm, fraction in classes:
        if low_mm is None:
            diameter = FINE_DIAMETER_M
        else:
            diameter = math.sqrt(low_mm * high_mm) * 1.0e-3
        total += fraction
        moment += fraction * diameter
    return moment / total


def ghost_bubble_velocity(
    height_above_bed_m: float,
    bubble_velocity: float,
    bubble_diameter: float,
    gas: dict,
    *,
    secondary: bool | None = None,
) -> float:
    """u_gb(h) = u_0(h) + (u_b,WS − u_0(h)) exp(−β_A h)，h 从床面起算。

    β_A = 1/(3 d_b,WS)。u_0 是该绝对高度的表观气速。secondary 为 None 时，
    绝对高度高于 1,4875 m 才计入二次风。
    Source: Hamel Eq.3.74, S. 43
    Source: Hamel Eq.3.75, S. 43
    Source: Hamel Eq.3.76, S. 43
    """
    absolute = BED_HEIGHT_M + height_above_bed_m
    if secondary is None:
        approach = superficial_velocity(absolute, gas)
    else:
        amount = gas["primary_mol_s"] + (gas["secondary_mol_s"] if secondary else 0.0)
        volume = amount * GAS_CONSTANT * BED_TEMPERATURE_K / PRESSURE_PA
        approach = volume / cross_section(absolute, gas)
    beta = 1.0 / (3.0 * bubble_diameter)
    return approach + (bubble_velocity - approach) * math.exp(-beta * height_above_bed_m)


def surface_state(epsilon_mf: float, phi_s: float, diameter_m: float, gas: dict):
    """用 T19 由床料质量反求 ρ_P，并取床面上边界的 u_b、d_b。"""
    area = gas["bed_area_m2"]
    u_bed = superficial_velocity(BED_HEIGHT_M / 2.0, gas)
    rho_seed = gas["rho_g"] * 2.0
    conditions = BedConditions(
        diameter_m,
        epsilon_mf,
        rho_seed,
        gas["rho_g"],
        gas["mu_g"],
        gas["mu_g"] / gas["rho_g"],
        PRESSURE_PA,
        1.0e-5,
        BED_DIAMETER_M,
        "technical",
        "B",
        True,
    )
    edges = np.linspace(0.0, BED_HEIGHT_M, 9)
    velocities = np.full(8, u_bed)
    rho_p = particle_density_from_bed_mass(
        BED_PRESSURE_DROP_PA * area / G,
        edges,
        velocities,
        area,
        N_ORIFICE,
        H_A_M,
        6,
        phi_s,
        conditions,
        HYDRO_SEARCH,
    )
    u_mf = minimum_fluidization_velocity(
        epsilon_mf, gas["mu_g"] / gas["rho_g"], phi_s, diameter_m, rho_p, gas["rho_g"]
    )
    updated = BedConditions(
        diameter_m,
        epsilon_mf,
        rho_p,
        gas["rho_g"],
        gas["mu_g"],
        gas["mu_g"] / gas["rho_g"],
        PRESSURE_PA,
        1.0e-5,
        BED_DIAMETER_M,
        "technical",
        "B",
        True,
    )
    profile = bed_profile(edges, velocities, u_mf, area, N_ORIFICE, H_A_M, 6, updated, HYDRO_SEARCH)
    return {
        "rho_p": float(rho_p),
        "u_mf": float(u_mf),
        "u_b": float(profile.u_b[0]),
        "d_b": float(profile.d_b[0]),
        "u_bed": float(u_bed),
    }


def freeboard_edges() -> np.ndarray:
    return np.linspace(0.0, FREEBOARD_HEIGHT_M, 9)


def case1_pieces() -> list[tuple[float, float]]:
    """两段光滑气速。第二段从风口之上开始，避免把二次风间断采样进同一段。"""
    port = SECONDARY_HEIGHT_M - BED_HEIGHT_M
    return [(0.0, port), (port + 1.0e-4, FREEBOARD_HEIGHT_M)]


def case1_scan() -> dict:
    """扫 ε_mf、φ_S 与各筛孔边界。失败的密度角点记为未检验，不换一组参数重试。"""
    gas = case1_gas()
    classes = size_bounds()
    bed_diameter = bed_representative_diameter_m(classes)
    diameters = [(None, "未检验：kurkela1992_fig1 最细一档没有下限")]
    seen = set()
    for low_mm, high_mm, _fraction in classes:
        for bound in (low_mm, high_mm):
            if bound is None or bound in seen:
                continue
            seen.add(bound)
            diameters.append((bound * 1.0e-3, f"筛孔 {bound} mm"))
    rows = []
    states = []
    for epsilon_mf in EPSILON_MF_GRID:
        for phi_s in PHI_S_GRID:
            try:
                state = surface_state(epsilon_mf, phi_s, bed_diameter, gas)
                state["error"] = ""
            except HydrodynamicsError as exc:
                state = {
                    "rho_p": math.nan,
                    "u_mf": math.nan,
                    "u_b": math.nan,
                    "d_b": math.nan,
                    "u_bed": math.nan,
                    "error": str(exc),
                }
            state["epsilon_mf"] = epsilon_mf
            state["phi_s"] = phi_s
            states.append(state)
    for diameter_m, label in diameters:
        for state in states:
            row = {
                "label": label,
                "diameter_m": diameter_m,
                "epsilon_mf": state["epsilon_mf"],
                "phi_s": state["phi_s"],
                "rho_p": state["rho_p"],
                "u_b": state["u_b"],
                "d_b": state["d_b"],
                "v_t": math.nan,
                "stable_m": [],
                "bed": "未检验",
                "local": "未检验",
                "overall": "未检验",
                "note": state["error"] or label,
            }
            if diameter_m is None or state["error"]:
                rows.append(row)
                continue
            particle = Particle(diameter_m, state["rho_p"], gas["rho_g"], gas["mu_g"], state["phi_s"])

            def u_of_h(height, bubble=state["u_b"], diameter=state["d_b"], used=gas):
                return ghost_bubble_velocity(height, bubble, diameter, used)

            samples = ejection_samples(state["u_b"])
            pieces = case1_pieces()
            assessed = assess_piecewise(u_of_h, pieces, freeboard_edges(), particle, samples)
            row["v_t"] = assessed.v_t_m_s
            row["stable_m"] = [BED_HEIGHT_M + item for item in assessed.stable_heights_m]
            row["bed"] = assessed.bed.status
            row["local"] = assessed.local.status
            row["overall"] = assessed.overall
            row["note"] = f"{label}；床面：{assessed.bed.detail}；当地：{assessed.local.detail}"
            rows.append(row)
    return {
        "gas": gas,
        "bed_diameter_m": bed_diameter,
        "states": states,
        "rows": rows,
        "u_bed": superficial_velocity(0.6, gas),
        "u_freeboard": superficial_velocity(2.0, gas),
    }


def _same(values: list[str]) -> str:
    unique = sorted(set(values))
    if unique == ["不会"]:
        return "不会"
    if unique == ["会"]:
        return "会"
    if "会" in unique and "不会" in unique:
        return "取决于 ρ_P 区间"
    if "会" in unique:
        return "会"
    if "未检验" in unique and "不会" in unique:
        return "未检验"
    return unique[0] if len(unique) == 1 else "未检验"


def summary_table(report: dict) -> list[dict]:
    """每个筛孔边界一行。密度扫描内部结论不一致时写成「取决于 ρ_P 区间」。"""
    grouped = []
    current = None
    for row in report["rows"]:
        key = (row["label"], row["diameter_m"])
        if current is None or current["key"] != key:
            if current is not None:
                grouped.append(current)
            current = {"key": key, "rows": []}
        current["rows"].append(row)
    if current is not None:
        grouped.append(current)
    lines = []
    for group in grouped:
        usable = [item for item in group["rows"] if item["diameter_m"] is not None and math.isfinite(item["rho_p"])]
        label = group["key"][0]
        diameter = group["key"][1]
        if diameter is None or not usable:
            lines.append(
                {
                    "粒级": label,
                    "v_t": "—",
                    "h*": "—",
                    "床面抛出": "未检验",
                    "当地生成": "未检验",
                    "综合": "未检验",
                }
            )
            continue
        velocities = [item["v_t"] for item in usable]
        heights = [height for item in usable for height in item["stable_m"]]
        lines.append(
            {
                "粒级": label,
                "v_t": f"{min(velocities):.3f}–{max(velocities):.3f} m/s",
                "h*": "无" if not heights else f"{min(heights):.3f}–{max(heights):.3f} m（自布风器）",
                "床面抛出": _same([item["bed"] for item in usable]),
                "当地生成": _same([item["local"] for item in usable]),
                "综合": _same([item["overall"] for item in usable]),
            }
        )
    return lines


def hover_sentence(lines: list[dict]) -> str:
    """一句预检结论，不裁定闭合。"""
    if any(item["综合"] == "取决于 ρ_P 区间" for item in lines):
        head = "取决于 ρ_P 区间"
    elif any(item["综合"] == "会" for item in lines):
        head = "会"
    elif all(item["综合"] == "不会" for item in lines):
        head = "不会"
    else:
        head = "未检验"
    return (
        f"悬停预检的综合结果是「{head}」。"
        f"这只是所列临时输入下两类参考轨迹的计算记录，不是悬停闭合的裁定。"
    )


def combination_lines(rows: list[dict]) -> list[dict]:
    """每个（筛孔、ε_mf、φ_S）一行。最细档没有下限时只留一行。"""
    lines = []
    missing = False
    for row in rows:
        if row["diameter_m"] is None:
            if missing:
                continue
            missing = True
            lines.append(
                {
                    "粒级": row["label"],
                    "ε_mf": "—",
                    "φ_S": "—",
                    "ρ_P": "—",
                    "v_t": "—",
                    "h*": "—",
                    "床面抛出": "未检验",
                    "当地生成": "未检验",
                    "综合": "未检验",
                }
            )
            continue
        if not math.isfinite(row["v_t"]):
            lines.append(
                {
                    "粒级": row["label"],
                    "ε_mf": f"{row['epsilon_mf']:.3f}",
                    "φ_S": f"{row['phi_s']:.2f}",
                    "ρ_P": "—",
                    "v_t": "—",
                    "h*": "—",
                    "床面抛出": "未检验",
                    "当地生成": "未检验",
                    "综合": "未检验",
                }
            )
            continue
        heights = row["stable_m"]
        lines.append(
            {
                "粒级": row["label"],
                "ε_mf": f"{row['epsilon_mf']:.3f}",
                "φ_S": f"{row['phi_s']:.2f}",
                "ρ_P": f"{row['rho_p']:.0f}",
                "v_t": f"{row['v_t']:.3f}",
                "h*": "无" if not heights else "、".join(f"{height:.3f}" for height in heights),
                "床面抛出": row["bed"],
                "当地生成": row["local"],
                "综合": row["overall"],
            }
        )
    return lines


def format_combinations(lines: list[dict]) -> str:
    header = "| 粒级 | ε_mf | φ_S | ρ_P / kg/m³ | v_t / m/s | h* / m | 床面抛出 | 当地生成 | 综合 |"
    rule = "|---|---|---|---|---|---|---|---|---|"
    body = [
        "| {粒级} | {ε_mf} | {φ_S} | {ρ_P} | {v_t} | {h*} | {床面抛出} | {当地生成} | {综合} |".format(**item)
        for item in lines
    ]
    return "\n".join([header, rule, *body])


def format_table(lines: list[dict]) -> str:
    header = "| 粒级 | v_t | 稳定高度 h* | 床面抛出 | 当地生成 | 综合 |"
    rule = "|---|---|---|---|---|---|"
    body = [
        f"| {item['粒级']} | {item['v_t']} | {item['h*']} | {item['床面抛出']} | {item['当地生成']} | {item['综合']} |"
        for item in lines
    ]
    return "\n".join([header, rule, *body])


def main() -> None:
    report = case1_scan()
    lines = summary_table(report)
    print("工况 1 悬停预检")
    print(f"床层表观气速 {report['u_bed']:.3f} m/s，自由空间直筒表观气速 {report['u_freeboard']:.3f} m/s")
    print("对照：KURKELA_1992 表 5 给出床层 0,74 m/s、自由空间 0,67 m/s。这里不向这两个数调整输入。")
    print(f"床料代表粒径 {report['bed_diameter_m'] * 1e3:.3f} mm（临时）")
    print()
    finite = [item for item in report["states"] if math.isfinite(item["rho_p"])]
    if finite:
        print(
            "ρ_P {:.0f}–{:.0f} kg/m³，u_b,WS {:.3f}–{:.3f} m/s，d_b,WS {:.4f}–{:.4f} m".format(
                min(item["rho_p"] for item in finite),
                max(item["rho_p"] for item in finite),
                min(item["u_b"] for item in finite),
                max(item["u_b"] for item in finite),
                min(item["d_b"] for item in finite),
                max(item["d_b"] for item in finite),
            )
        )
    failed = [item for item in report["states"] if item["error"]]
    for item in failed:
        print(f"未反求到 ρ_P：ε_mf={item['epsilon_mf']}，φ_S={item['phi_s']}，{item['error']}")
    print()
    print(format_combinations(combination_lines(report["rows"])))
    print()
    print("按筛孔归并（密度角点不一致时写「取决于 ρ_P 区间」）：")
    print(format_table(lines))
    print()
    print(hover_sentence(lines))
    print()
    print("输入：")
    for note in INPUT_NOTES:
        print(f"- {note}")


if __name__ == "__main__":
    main()
