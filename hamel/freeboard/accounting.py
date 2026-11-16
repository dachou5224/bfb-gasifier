"""自由空间参考轨迹的归宿、记账与输运系数。

模型边界选择：ADR 0026（目标所有者 2026-10-19）。
格内常系数解析推进与界面稳定点：ADR 0047（目标所有者 2026-10-28；第 1 项为解释性选择，第 2、3 项为求解手段）。
粒级内取样、半格捕获区与源格→目标格：ADR 0048（解释性选择）。
悬停状态的粒级内漂移：解释性选择：ADR 0051 第 1、2 项。
捕获区取整格、进入即在格界，只是 ADR 0048 兼容测试用的旧规则，不是生产取样。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard._common import FreeboardError, _finite, _positive
from hamel.freeboard.motion import (
    _boundary_time,
    advance_cell,
    cell_height_at_velocity,
    cell_mean_velocity,
    drag_groups,
    motion_coefficients,
    terminal_settling_velocity,
)

_FREE_ABSENT = "自由状态不存在"
# |ū_g − v_t| 不超过这个倍数乘 v_t 时就认为该格中性。求解手段：ADR 0047 第 3 项。
_NEUTRAL = 1.0e-12
# 粒级内取样：ADR 0048 第 1 项。粒径在 ln d 上等分，球形度取截断正态的等概率分位中点。
SAMPLE_DIAMETERS = 8
SAMPLE_SPHERICITY = 8
SPHERICITY_SIGMA = 0.1
SPHERICITY_BOUNDS = (0.3, 1.0)
SPHERICITY_MEAN = 0.8
DEFAULT_REGION = "half"
# 悬停状态的下沿去向。解释性选择：ADR 0051 第 2 项。
DRIFT_MIGRATION = "粒级迁移"
DRIFT_HOVER = "漂移-悬停"
DRIFT_FREE = "漂移-自由"
DRIFT_BURN = "本级燃尽"
DRIFT_EMPTY = "空窗口"
# v_t(d, φ) = ū_g 二分求根的 d 区间。s 与 q 都随 d 单调减小，v_t 随 d 单调增大。
_ROOT_BOUNDS_M = (1.0e-7, 1.0e-2)


@dataclass(frozen=True)
class ClassWindow:
    """粒级边界 [a_k, b_k] 与 ADR 0051 第 2 项的 ρ_k = ln(d_k/d_(k−1))/ln(b_k/a_k)。"""

    lower_m: float
    upper_m: float
    rho_log: float


@dataclass(frozen=True)
class DriftWindow:
    """一个（格，粒级，φ_p）的悬停窗口、下沿去向与出口因子。ADR 0051 第 2 项。"""

    sphericity: float
    low_m: float
    high_m: float
    fate: str
    length: float
    exit_share: float
    weight: float = 0.0


@dataclass(frozen=True)
class StableEquilibrium:
    """一个稳定界面和它的捕获区：R* 为界面两侧各半格。"""

    height_m: float
    cells: tuple[int, ...]
    half_width_m: float


@dataclass(frozen=True)
class TrajectoryFate:
    """一条轨迹的归宿，以及进入捕获区之前的停留和穿出。"""

    fate: str
    time_s: float
    height_m: float
    velocity_m_s: float
    velocity_reversals: int
    capture_time_s: float | None
    captured_at_start: bool
    entry_cell: int | None
    entered_cell: int | None
    entry_kind: str | None
    entry_count: int
    # 本条的解析段数（advance_cell 调用次数）。求解器参数：ADR 0050。
    segments: int
    residence_s: tuple[float, ...]
    # 进入捕获区前 ∫|u_g−u_P|dt，与停留时间同一段。解释性选择：ADR 0035。
    slip_m: tuple[float, ...]
    upward_free: tuple[float, ...]
    downward_free: tuple[float, ...]
    # 进入记账的边界穿越序列（源格, 目标格或 None, 方向）。解释性选择：ADR 0058 第 1 项。
    crossings: tuple = ()


@dataclass(frozen=True)
class CellTransport:
    """一格一个粒级。自由状态不存在时，四个系数和 m̂ 都是 None。"""

    theta_up: float
    theta_down: float
    free_state: str
    a_f: float | None
    b_f: float | None
    c_up: float | None
    c_down: float | None
    m_hat: float | None
    mean_slip_m_s: float | None = None
    hover_slip_m_s: float | None = None
    # 悬停状态的粒级内漂移：进入上方格悬停、自由状态的因子与粒级迁移因子。解释性选择：ADR 0051。
    drift_hover: float = 0.0
    drift_free: float = 0.0
    migration_mu: float = 0.0
    windows: tuple = ()


@dataclass(frozen=True)
class ReferenceAccounting:
    """一个粒级的参考记账。份额对单位抛出，当地轨迹按每格权重 1。"""

    stable_heights_m: tuple[float, ...]
    capture_cells: tuple[int, ...]
    return_share: float
    elutriation_share: float
    capture_share: float
    ejection_free: float
    ejection_hover: float
    local_return: float
    local_elutriation: float
    cells: tuple[CellTransport, ...]
    roots: tuple = ()
    # ADR 0050 第 5 项：全部参考轨迹的段数与归宿时刻统计。只报告。
    max_segments: int = 0
    max_time_s: float = 0.0
    p99_time_s: float = 0.0
    longest: tuple | None = None
    # ADR 0058 第 1 项：(端点轨迹数, 定位用中间轨迹数, 第三种签名的区间数)。只报告。
    boundary: tuple = (0, 0, 0)


def _edges_and_gas(edges_m, rho_g, mu_g) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = np.asarray(edges_m, dtype=float)
    densities = np.asarray(rho_g, dtype=float)
    viscosities = np.asarray(mu_g, dtype=float)
    count = edges.size - 1
    if edges.ndim != 1 or count < 1 or np.any(np.diff(edges) <= 0.0):
        raise FreeboardError("自由空间边界必须严格递增")
    if densities.shape != (count,) or viscosities.shape != (count,):
        raise FreeboardError("气体密度和粘度必须每个格一个数")
    return edges, densities, viscosities


def stable_equilibria(
    edges_m,
    gas_velocity,
    rho_g,
    mu_g,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
) -> tuple[StableEquilibrium, ...]:
    """稳定点只在界面：下格 ū_g > v_t、上格 ū_g < v_t；R* 为界面两侧各半格。

    格内气速为常数时格内没有平衡点。v_t 用各自格的气体物性。
    求解手段：ADR 0047 第 1、3 项；捕获区取半格：解释性选择，ADR 0048 第 2 项。
    Source: Hamel Eq.3.82, S. 44
    """
    edges, densities, viscosities = _edges_and_gas(edges_m, rho_g, mu_g)
    count = edges.size - 1
    diameter = _positive("diameter_m", diameter_m)
    density_p = _positive("rho_p", rho_p)
    settling = [
        terminal_settling_velocity(float(densities[index]), float(viscosities[index]), diameter, density_p, phi_s)
        for index in range(count)
    ]
    speeds = [
        cell_mean_velocity(gas_velocity, index, float(edges[index]), float(edges[index + 1]))
        for index in range(count)
    ]
    for index in range(count):
        if abs(speeds[index] - settling[index]) <= _NEUTRAL * settling[index]:
            raise FreeboardError(
                f"格 {index} 中性（ū_g={speeds[index]:.6g} m/s，v_t={settling[index]:.6g} m/s），不能判定悬停"
            )
    found = []
    for face in range(1, count):
        lower, upper = face - 1, face
        if speeds[lower] > settling[lower] and speeds[upper] < settling[upper]:
            half = 0.5 * min(float(edges[face] - edges[lower]), float(edges[upper + 1] - edges[face]))
            found.append(StableEquilibrium(float(edges[face]), (lower, upper), float(half)))
    return tuple(found)


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def _normal_density(value: float) -> float:
    return math.exp(-0.5 * value * value) / math.sqrt(2.0 * math.pi)


def _truncated_mass(location: float, sigma: float = SPHERICITY_SIGMA) -> float:
    low, high = SPHERICITY_BOUNDS
    return _normal_cdf((high - location) / sigma) - _normal_cdf((low - location) / sigma)


def _truncated_mean(location: float, sigma: float = SPHERICITY_SIGMA) -> float:
    low, high = SPHERICITY_BOUNDS
    alpha, beta = (low - location) / sigma, (high - location) / sigma
    mass = _truncated_mass(location, sigma)
    return location + sigma * (_normal_density(alpha) - _normal_density(beta)) / mass


def sphericity_points(count: int = SAMPLE_SPHERICITY, sigma: float = SPHERICITY_SIGMA) -> tuple[float, ...]:
    """[0,3; 1,0] 截断正态的等概率分位中点，位置参数取使截断后均值为 0,8。解释性选择：ADR 0048 第 1 项。

    φ_S 的范围与平均值取自 S. 45（原文只给范围与平均值）。
    Source: Hamel Eq.3.86, S. 45
    """
    low, high = SPHERICITY_BOUNDS
    lower, upper = low, high
    for _ in range(200):
        middle = 0.5 * (lower + upper)
        if _truncated_mean(middle, sigma) < SPHERICITY_MEAN:
            lower = middle
        else:
            upper = middle
    location = 0.5 * (lower + upper)
    mass, floor = _truncated_mass(location, sigma), _normal_cdf((low - location) / sigma)
    points = []
    for index in range(int(count)):
        target = (index + 0.5) / int(count)
        lower, upper = low, high
        for _ in range(80):
            middle = 0.5 * (lower + upper)
            if (_normal_cdf((middle - location) / sigma) - floor) / mass < target:
                lower = middle
            else:
                upper = middle
        points.append(0.5 * (lower + upper))
    return tuple(points)


def within_class_samples(lower_m, upper_m, diameters: int = SAMPLE_DIAMETERS, sphericities: int = SAMPLE_SPHERICITY, sigma: float = SPHERICITY_SIGMA) -> tuple[tuple[float, float, float], ...]:
    """粒级内 N_d × N_φ 个取样点与权重。质量在 ln d 上均匀。解释性选择：ADR 0048 第 1 项。

    取样细化使类别分得更细，见 S. 49；这些样本喂给式 (3.88) 的轨迹。
    Source: Hamel Eq.3.88, S. 46
    """
    low, high = _positive("lower_m", lower_m), _positive("upper_m", upper_m)
    if high <= low:
        raise FreeboardError("粒级上边界必须高于下边界")
    count_d = int(diameters)
    step = (math.log(high) - math.log(low)) / count_d
    points = [math.exp(math.log(low) + (index + 0.5) * step) for index in range(count_d)]
    phi = sphericity_points(int(sphericities), float(sigma))
    weight = 1.0 / (count_d * len(phi))
    return tuple((diameter, value, weight) for diameter in points for value in phi)


def _settling_diameter(rho_g: float, mu_g: float, rho_p: float, phi_s: float, speed: float) -> float | None:
    """v_t(d, φ_S) = ū_g 的根；v_t 随 d 单调，二分。ADR 0051 第 2 项。

    根落在区间下界以外时返回 0（窗口下沿低于取样下界），落在上界以外时返回 None（无根，窗口为空）。
    Source: Hamel Eq.3.85, S. 45
    """
    low, high = _ROOT_BOUNDS_M

    def gap(diameter: float) -> float:
        return terminal_settling_velocity(rho_g, mu_g, diameter, rho_p, phi_s) - float(speed)

    if gap(low) >= 0.0:
        return 0.0
    if gap(high) <= 0.0:
        return None
    for _ in range(80):
        middle = 0.5 * (low + high)
        if middle <= low or middle >= high:
            break
        if gap(middle) < 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def drift_windows(
    edges_m,
    gas_velocity,
    rho_g,
    mu_g,
    rho_p: float,
    sphericities,
    class_window: ClassWindow,
    first_class: bool = False,
) -> tuple[tuple[DriftWindow, ...], ...]:
    """逐（格，φ_p）的悬停窗口 [d′_lo, d′_hi]、下沿去向与出口因子。ADR 0051 第 2 项。

    格 j 的窗口上沿由格 j 自己的物性与 ū_g,j 定义、下沿由格 j+1 的定义，故相邻界面的窗口
    首尾相接；下沿为 a_k 时是粒级迁移（k = 1 留在本状态燃尽），其余按格 j+1 自己的窗口是否
    非空判定去向。返回的每一项对应该格该 φ_p，空窗口也保留（fate = 空窗口）。
    Source: Hamel Eq.5.11, S. 86
    """
    edges, densities, viscosities = _edges_and_gas(edges_m, rho_g, mu_g)
    count = edges.size - 1
    density_p = _positive("rho_p", rho_p)
    speeds = [
        cell_mean_velocity(gas_velocity, index, float(edges[index]), float(edges[index + 1]))
        for index in range(count)
    ]
    lower, upper = _positive("class_lower_m", class_window.lower_m), _positive("class_upper_m", class_window.upper_m)
    if upper <= lower:
        raise FreeboardError("粒级上边界必须高于下边界")
    table = []
    for phi in sphericities:
        facet = [
            _settling_diameter(float(densities[index]), float(viscosities[index]), density_p, float(phi), speeds[index])
            for index in range(count)
        ]
        row = []
        for cell in range(count - 1):
            top, bottom = facet[cell], facet[cell + 1]
            if top is None or bottom is None or top <= bottom:
                row.append(DriftWindow(float(phi), 0.0, 0.0, DRIFT_EMPTY, 0.0, 0.0))
                continue
            high, low = min(top, upper), max(bottom, lower)
            if low >= high:
                row.append(DriftWindow(float(phi), low, high, DRIFT_EMPTY, 0.0, 0.0))
                continue
            if bottom <= lower:
                fate = DRIFT_BURN if first_class else DRIFT_MIGRATION
            elif cell + 2 < count and facet[cell + 2] is not None and facet[cell + 2] < facet[cell + 1]:
                fate = DRIFT_HOVER
            else:
                fate = DRIFT_FREE
            length = class_window.rho_log * math.log(high / low)
            share = 0.0 if length <= 0.0 else 1.0 / math.expm1(3.0 * length)
            row.append(DriftWindow(float(phi), low, high, fate, length, share))
        table.append(tuple(row))
    # 转成按格优先：table[cell][p]。
    return tuple(tuple(table[p][cell] for p in range(len(table))) for cell in range(count - 1))


def _window_factors(table, captured: dict, cell: int, first_class: bool) -> dict:
    """一格一个粒级的 ω_p、δ^h、δ^f、μ；没有给出粒级边界时全为 0。ADR 0051 第 2 项。"""
    if not table or cell >= len(table):
        return {"windows": (), "drift_hover": 0.0, "drift_free": 0.0, "migration_mu": 0.0}
    total = float(sum(captured.values()))
    open_count = sum(1 for item in table[cell] if item.fate != DRIFT_EMPTY)
    share = 1.0 / open_count if total == 0.0 and open_count else 0.0
    weighted = [
        replace(item, weight=(captured.get(item.sphericity, 0.0) / total) if total > 0.0 else (share if item.fate != DRIFT_EMPTY else 0.0))
        for item in table[cell]
    ]
    factors = {DRIFT_MIGRATION: 0.0, DRIFT_HOVER: 0.0, DRIFT_FREE: 0.0}
    for item in weighted:
        if item.fate in factors:
            factors[item.fate] += item.weight * item.exit_share
    if first_class:
        factors[DRIFT_MIGRATION] = 0.0
    return {
        "windows": tuple(weighted),
        "drift_hover": factors[DRIFT_HOVER],
        "drift_free": factors[DRIFT_FREE],
        "migration_mu": factors[DRIFT_MIGRATION],
    }


def _classify(
    count: int,
    fate: str,
    elapsed: float,
    height: float,
    velocity: float,
    reversals: int,
    start_cell: int,
    segments: list[tuple],
    crossings: list[tuple[float, int, int | None, str]],
    midpoints: list[tuple[float, int, str]],
    capture_cells: frozenset[int],
    region: str,
    segments_used: int = 0,
    captured_face: tuple[int, ...] | None = None,
    u_min: float = 0.0,
    capture_height: float | None = None,
) -> TrajectoryFate:
    entry_index = None
    source = None
    entered = None
    capture_time = None
    kind = None
    if region == "cell":
        source = start_cell if start_cell in capture_cells else None
        entered = start_cell if start_cell in capture_cells else None
        capture_time = 0.0 if entered is not None else None
        kind = "up" if entered is not None else None
        for index, (time, from_cell, to_cell, direction) in enumerate(crossings):
            enters = to_cell is not None and to_cell in capture_cells and from_cell not in capture_cells
            leaves = from_cell in capture_cells and (to_cell is None or to_cell not in capture_cells)
            if enters:
                entry_index = index
                source, entered, capture_time, kind = from_cell, to_cell, time, direction
            elif leaves:
                source, entered, capture_time, kind = None, None, None, None
    elif midpoints:
        # 半格捕获区：进入 R* 即越过格中点；取最后一次进入。
        capture_time, source, kind = midpoints[-1]
        entered = source if kind == "up" else source - 1
    elif fate == "capture" and captured_face is not None:
        # 起点本来就在 R* 内且此后没有离开：进入时刻为 0，悬停格取该界面的下方格。
        capture_time, source, kind = 0.0, captured_face[0], "up"
        entered = captured_face[0]
    captured = fate == "capture"
    if captured and (capture_time is None or entered is None or entered < 0):
        raise FreeboardError("判定为捕获，但轨迹没有确定的来源格")
    if not captured:
        capture_time, source, entered, kind, entry_index = None, None, None, None, None
    at_start = captured and entry_index is None and capture_time == 0.0
    residence = [0.0] * count
    slip = [0.0] * count
    upward_free = [0.0] * count
    downward_free = [0.0] * count
    for item in segments:
        cell, start, end, piece, h0, h1 = item
        if capture_time is None or end <= capture_time:
            span, h_end = end - start, h1
        elif start < capture_time:
            span = capture_time - start
            h_end = h1 if capture_height is None else capture_height
        else:
            continue
        # 自由部分的滑移与停留乘同一比例。解释性选择：ADR 0057。
        counted = span if not u_min else min(span, abs(h_end - h0) / u_min)
        residence[cell] += counted
        slip[cell] += piece * counted / (end - start)
    accounted = []
    for time, from_cell, to_cell, direction in crossings:
        if captured and capture_time is not None and time >= capture_time:
            continue
        if direction == "up":
            upward_free[from_cell] += 1.0
        else:
            downward_free[from_cell] += 1.0
        accounted.append((int(from_cell), None if to_cell is None else int(to_cell), str(direction)))
    return TrajectoryFate(
        fate,
        elapsed,
        height,
        velocity,
        reversals,
        capture_time if captured else None,
        at_start,
        source if captured else None,
        entered if captured else None,
        kind if captured else None,
        len(midpoints) if region == "half" else len(crossings),
        int(segments_used),
        tuple(residence),
        tuple(slip),
        tuple(upward_free),
        tuple(downward_free),
        tuple(accounted),
    )


def trace_until_fate(
    height_0: float,
    velocity_0: float,
    edges_m,
    gas_velocity,
    rho_g,
    mu_g,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
    config: SolverConfig,
    equilibria: tuple[StableEquilibrium, ...] | None = None,
    scale: float = 1.0,
    region: str = DEFAULT_REGION,
) -> TrajectoryFate:
    """逐格常系数解析推进到回落、带出，或穿过稳定界面时判为捕获。

    格内气速取 u_gb(h) 的高度平均。穿过稳定界面时 |u_P| ≤ δ_u，且按闭式 h(u)
    求得的越界距离 ≤ δ_h，即判为捕获；进入捕获区的时刻按最后一次越过格中点记账，
    只有捕获区之前的穿出才计入自由状态。解释性选择：ADR 0026 第 2 项、ADR 0047 第 1 项、
    ADR 0048 第 2 项、ADR 0057；求解手段：ADR 0047 第 2、3 项。
    Source: Hamel Eq.3.88, S. 46
    Source: Hamel Eq.3.93, S. 47
    """
    edges, densities, viscosities = _edges_and_gas(edges_m, rho_g, mu_g)
    count = edges.size - 1
    bounds = tuple(float(edge) for edge in edges)
    height = _finite("height_0", height_0)
    velocity = _finite("velocity_0", velocity_0)
    if not bounds[0] <= height <= bounds[-1]:
        raise FreeboardError("起步高度不在自由空间内")
    speeds = tuple(cell_mean_velocity(gas_velocity, index, bounds[index], bounds[index + 1]) for index in range(count))
    if equilibria is None:
        equilibria = stable_equilibria(edges_m, gas_velocity, rho_g, mu_g, diameter_m, rho_p, phi_s)
    interfaces = {item.cells[1]: item for item in equilibria}
    delta_h = scale * _positive("freeboard_delta_h_m", config.freeboard_delta_h_m)
    delta_u = scale * _positive("freeboard_delta_u_m_s", config.freeboard_delta_u_m_s)
    # 解析段数保护代替物理时间截止。求解器参数：ADR 0050。
    # 停留下限速度只改记账，不改轨迹。解释性选择：ADR 0057。
    u_min = float(config.freeboard_residence_floor_m_s)
    if not math.isfinite(u_min) or u_min < 0.0:
        raise FreeboardError("freeboard_residence_floor_m_s 必须是有限非负数")
    max_segments = int(config.freeboard_max_segments)
    if max_segments < 1:
        raise FreeboardError("freeboard_max_segments 必须是不小于 1 的整数")
    diameter = _positive("diameter_m", diameter_m)
    density_p = _positive("rho_p", rho_p)
    groups = tuple(
        drag_groups(float(densities[index]), float(viscosities[index]), diameter, density_p, phi_s)
        for index in range(count)
    )

    def acceleration(index: int, velocity_now: float) -> float:
        stokes, quadratic, weight = groups[index]
        relative = velocity_now - speeds[index]
        return -stokes * relative - quadratic * relative * abs(relative) - weight

    def coefficients(index: int, velocity_now: float) -> tuple[float, float, float]:
        return motion_coefficients(
            speeds[index], bool(velocity_now > speeds[index]), float(densities[index]), float(viscosities[index]),
            diameter, density_p, phi_s,
        )

    def crossing_captures(face: int, velocity_now: float) -> StableEquilibrium | None:
        """穿过稳定界面时的捕获判据：速度与按闭式 h(u) 算的越界距离。"""
        item = interfaces.get(face)
        if item is None or abs(velocity_now) > delta_u:
            return None
        entered = item.cells[1] if velocity_now > 0.0 else item.cells[0]
        a, b, c = coefficients(entered, velocity_now)
        turn = cell_height_at_velocity(a, b, c, velocity_now, bounds[face], 0.0)
        if turn is None or abs(turn - bounds[face]) > delta_h:
            return None
        return item

    cell = min(max(int(np.searchsorted(edges, height, side="right") - 1), 0), count - 1)
    start_cell = cell
    elapsed = 0.0
    open_time = 0.0
    reversals = 0
    segments_used = 0
    segments: list[tuple[int, float, float, float]] = []
    crossings: list[tuple[float, int, int | None, str]] = []
    midpoints: list[tuple[float, int, str]] = []
    pending = 0.0
    open_height = height
    outcome = None
    chosen = None
    span = bounds[-1] - bounds[0]
    for index in range(count):
        if abs(height - 0.5 * (bounds[index] + bounds[index + 1])) <= 1.0e-12 * span:
            midpoints.append((0.0, index, "up" if velocity > 0.0 else "down"))

    def close(end: float) -> None:
        nonlocal open_time, pending, open_height
        if end > open_time:
            segments.append((cell, open_time, end, pending, open_height, height))
            open_time, open_height = end, height
        pending = 0.0

    while segments_used < max_segments:
        if height <= bounds[0] and (velocity < 0.0 or (velocity == 0.0 and acceleration(0, 0.0) < 0.0)):
            close(elapsed)
            crossings.append((elapsed, cell, None, "down"))
            outcome = "return"
            break
        if height >= bounds[-1] and (velocity > 0.0 or (velocity == 0.0 and acceleration(count - 1, 0.0) > 0.0)):
            close(elapsed)
            crossings.append((elapsed, cell, None, "up"))
            outcome = "elutriate"
            break
        a, b, c = coefficients(cell, velocity)
        kind, tau, height_1, velocity_1 = advance_cell(
            height, velocity, bounds[cell], bounds[cell + 1], a, b, c, speeds[cell],
        )
        segments_used += 1
        if kind == "limit":
            break
        # 段内 u 不穿过 ū_g 与 0，故 ∫|u_g−u_P|dτ = |ū_g·τ − Δh|，按段解析。
        pending += abs(speeds[cell] * tau - (height_1 - height))
        if kind == "zero" or velocity * velocity_1 < 0.0:
            reversals += 1
        middle = 0.5 * (bounds[cell] + bounds[cell + 1])
        if tau > 0.0 and (height - middle) * (height_1 - middle) < 0.0:
            reached = _boundary_time(a, b, c, velocity, height, middle, tau)
            if reached is not None:
                midpoints.append((elapsed + reached, cell, "up" if height_1 > height else "down"))
        elapsed += tau
        height, velocity = height_1, velocity_1
        close(elapsed)
        if kind in ("gas", "zero"):
            continue
        step = 1 if kind == "upper" else -1
        face = cell + 1 if step > 0 else cell
        target = cell + step
        captured = crossing_captures(face, velocity)
        if target < 0:
            crossings.append((elapsed, cell, None, "down"))
            outcome = "return"
            break
        if target >= count:
            crossings.append((elapsed, cell, None, "up"))
            outcome = "elutriate"
            break
        crossings.append((elapsed, cell, target, "up" if step > 0 else "down"))
        cell = target
        if captured is not None:
            chosen = captured
            outcome = "capture"
            break
    if outcome is None:
        # 段数保护：报告粒径、球形度、起点、所在格、已走段数与已用物理时间。求解器参数：ADR 0050。
        raise FreeboardError(
            f"参考轨迹达到段数上限 {max_segments}（d={diameter:.4g} m，φ_S={phi_s:.4g}，"
            f"起点 h0={height_0:.4g} m、u0={velocity_0:.4g} m/s，当前格 {cell}，已用 {elapsed:.4g} s）"
        )
    capture_cells = frozenset(chosen.cells) if chosen is not None else frozenset()
    cap_h = None
    if region != "cell" and midpoints:
        index = midpoints[-1][1]
        cap_h = 0.5 * (bounds[index] + bounds[index + 1])
    return _classify(
        count, outcome, elapsed, height, velocity, reversals, start_cell, segments, crossings, midpoints,
        capture_cells, region, segments_used, None if chosen is None else chosen.cells, u_min, cap_h,
    )


def _sum_ledgers(ledgers, count: int, weight: float):
    """把若干条轨迹的停留、滑移与自由状态穿出按（样本权重 × 轨迹权重）加总。"""
    parts = ([0.0] * count, [0.0] * count, [0.0] * count, [0.0] * count)
    for ledger, share in ledgers:
        scale = weight * share
        for part, values in zip(parts, (ledger.residence_s, ledger.slip_m, ledger.upward_free, ledger.downward_free)):
            for cell in range(count):
                part[cell] += scale * values[cell]
    return parts


def _touch_max(touched: dict, cell: int, value: float) -> None:
    """按节点累加：同一节点在同一格的权重取各轨迹位权重的最大值。ADR 0058 第 1 项。"""
    touched[cell] = max(touched.get(cell, 0.0), value)

def _signature(ledger: TrajectoryFate) -> tuple:
    """离散记账签名：归宿、起点即在捕获区、进入格与来源格、方向、次数、穿越序列。ADR 0058 第 1 项。"""
    return (ledger.fate, ledger.captured_at_start, ledger.entry_cell, ledger.entered_cell, ledger.entry_kind,
            ledger.entry_count, tuple(map(int, ledger.upward_free)), tuple(map(int, ledger.downward_free)), ledger.crossings)


def boundary_weights(xs, signatures, probe, delta: float, bisections: int):
    """按 ADR 0058 第 1 项把 [x_0, x_(N_d+1)] 分给各节点的质量权重（单位 w，合计 N_d）。

    probe(x) 返回 (新节点下标, 该处签名)，由调用方登记。解释性选择：ADR 0058 第 1 项。
    Source: Hamel Eq.3.88, S. 46
    """
    tolerance = delta * 2.0 ** (-int(bisections))
    weights: dict[int, float] = {}
    counts = {"bisections": 0, "thirds": 0}
    if len(set(signatures)) == 1:  # 处处相同：与现行的等权求和逐位相同。
        return dict.fromkeys(range(1, len(xs) - 1), 1.0), counts

    def add(index, value):
        weights[index] = weights.get(index, 0.0) + value

    def locate(left, right, index_left, index_right, signature_left, signature_right):
        if right - left <= tolerance:
            middle = 0.5 * (left + right)
            add(index_left, (middle - left) / delta)
            add(index_right, (right - middle) / delta)
            return
        middle = 0.5 * (left + right)
        counts["bisections"] += 1
        index_middle, signature_middle = probe(middle)
        if signature_middle == signature_left:
            add(index_left, (middle - left) / delta)
            locate(middle, right, index_left, index_right, signature_left, signature_right)
        elif signature_middle == signature_right:
            add(index_right, (right - middle) / delta)
            locate(left, middle, index_left, index_right, signature_left, signature_right)
        else:
            counts["thirds"] += 1
            locate(left, middle, index_left, index_middle, signature_left, signature_middle)
            locate(middle, right, index_middle, index_right, signature_middle, signature_right)

    for position in range(len(xs) - 1):
        left, right = float(xs[position]), float(xs[position + 1])
        width = (right - left) / delta
        if signatures[position] == signatures[position + 1]:
            if position == 0 and len(xs) > 2:
                add(1, width)
            elif position == len(xs) - 2 and len(xs) > 2:
                add(position, width)
            else:
                add(position, 0.5 * width)
                add(position + 1, 0.5 * width)
        else:
            locate(left, right, position, position + 1, signatures[position], signatures[position + 1])
    return weights, counts


def account_freeboard(
    edges_m,
    gas_velocity,
    rho_g,
    mu_g,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
    ejection_speeds,
    ejection_shares,
    direct_return_share: float,
    config: SolverConfig,
    scale: float = 1.0,
    samples=None,
    region: str = DEFAULT_REGION,
    class_window: ClassWindow | None = None,
    first_class: bool = False,
) -> ReferenceAccounting:
    """由床面样本和各格当地轨迹得到 a^f、b^f、c↑、c↓、θ↑、θ↓ 与抛出拆分。

    每个粒级按 ADR 0048 第 1 项的取样点分别做参考轨迹，按权重汇总；m̂ = 0 时自由状态不存在。
    给出粒级边界时同时按 ADR 0051 第 2 项求悬停窗口与 δ^h、δ^f、μ；不给时不计漂移（只用于界面测试）。
    解释性选择：ADR 0016 第 9、10 条，ADR 0026 第 2、3 项，ADR 0048 第 1、2 项。
    解释性选择：ADR 0051 第 2 项。
    Source: Hamel Eq.2.6, S. 17
    Source: Hamel Eq.3.82, S. 44
    Source: Hamel Eq.5.11, S. 86
    """
    edges, densities, viscosities = _edges_and_gas(edges_m, rho_g, mu_g)
    count = edges.size - 1
    speeds_in = np.asarray(ejection_speeds, dtype=float)
    shares = np.asarray(ejection_shares, dtype=float)
    if speeds_in.ndim != 1 or speeds_in.shape != shares.shape:
        raise FreeboardError("抛出样本的速度和份额必须是同样长的一维数组")
    if speeds_in.size and (np.any(speeds_in <= 0.0) or np.any(shares < 0.0)):
        raise FreeboardError("抛出样本的初速必须为正，份额不能为负")
    direct = _finite("direct_return_share", direct_return_share)
    if direct < 0.0 or abs(float(shares.sum()) + direct - 1.0) > 1.0e-8:
        raise FreeboardError("抛出份额与直接回落之和必须为 1，直接回落不能为负")
    if region not in {"half", "cell"}:
        raise FreeboardError(f"不能识别的捕获区规则 {region}")
    entries = ((float(diameter_m), float(phi_s), 1.0),) if samples is None else tuple(samples)
    if not entries:
        raise FreeboardError("粒级内取样不能为空")
    bed_surface = float(edges[0])
    names = ("res", "slip", "up", "down", "c_up", "c_down", "slip_t", "w_t")
    totals = {name: [0.0] * count for name in names}
    local_totals = {name: [0.0] * count for name in names}
    theta_up, theta_down = [0.0] * count, [0.0] * count
    returned, elutriated, captured, local_return, local_elutriation = direct, 0.0, 0.0, 0.0, 0.0
    statistics: list[tuple] = []
    captured_phi: list[dict] = [{} for _ in range(count)]
    facets = tuple(dict.fromkeys(float(phi) for _d, phi, _w in entries))
    windows = () if class_window is None else drift_windows(edges, gas_velocity, densities, viscosities, rho_p, facets, class_window, bool(first_class))
    first = stable_equilibria(edges, gas_velocity, densities, viscosities, float(entries[0][0]), rho_p, float(entries[0][1]))
    bisections = int(config.freeboard_boundary_bisections)  # 解释性选择：ADR 0058 第 1 项。
    diameters = tuple(dict.fromkeys(float(item[0]) for item in entries))
    phis = tuple(dict.fromkeys(float(item[1]) for item in entries))
    grid = samples is not None and len(diameters) * len(phis) == len(entries) and len(diameters) >= 2
    if grid:
        step = math.log(diameters[1] / diameters[0])
        xs = [math.log(diameters[0]) - 0.5 * step] + [math.log(item) for item in diameters] + [math.log(diameters[-1]) + 0.5 * step]
        weight, wanted = float(entries[0][2]), 1.0 / len(phis)
    positions = [("bed", i, float(shares[i])) for i in range(speeds_in.size)] + [("local", c, 1.0) for c in range(count)]
    counts = {"endpoints": 2 * len(phis) * len(positions) if grid else 0, "bisections": 0, "thirds": 0}
    # 第一阶段：逐 φ_S 行、逐轨迹位沿 x = ln d 定位，得到该轨迹位的“轨迹–权重”表。ADR 0058 第 1 项。
    plans, lookup = [], {}
    for phi, single in ([(item, None) for item in phis] if grid else [(float(item[1]), item) for item in entries]):
        nodes, by_key, cache = [], {}, {}

        def add_node(x, diameter=None, phi=phi):
            key = round(float(x), 12)
            if key not in by_key:
                value = math.exp(float(x)) if diameter is None else float(diameter)
                terminal = [terminal_settling_velocity(float(densities[c]), float(viscosities[c]), value, rho_p, phi) for c in range(count)]
                nodes.append({"diameter": value, "terminal": terminal,
                              "equilibria": stable_equilibria(edges, gas_velocity, densities, viscosities, value, rho_p, phi)})
                by_key[key] = len(nodes) - 1
            return by_key[key]

        base = [add_node(x, d) for x, d in zip(xs, [math.exp(xs[0])] + list(diameters) + [math.exp(xs[-1])])] if grid else [add_node(math.log(single[0]), single[0])]
        plan = {}
        for kind, where, _share in positions:
            def ledger_at(index, kind=kind, where=where):
                key = (kind, where, index)
                if key not in cache:
                    node = nodes[index]
                    if kind == "bed":
                        height, speed = bed_surface, float(speeds_in[where])
                    else:
                        height = 0.5 * (edges[where] + edges[where + 1])
                        speed = cell_mean_velocity(gas_velocity, where, edges[where], edges[where + 1]) - node["terminal"][where]
                    cache[key] = trace_until_fate(height, speed, edges, gas_velocity, densities, viscosities, node["diameter"], rho_p, phi,
                                                  config, node["equilibria"], scale, region)
                    # ADR 0050 第 5 项：段数与归宿时刻，只报告。
                    statistics.append((cache[key].segments, cache[key].time_s, node["diameter"], phi, "床面抛出" if kind == "bed" else f"格 {where} 中点"))
                return cache[key]

            if grid:
                signatures = [_signature(ledger_at(index)) for index in base]

                def probe(x):
                    index = add_node(x)
                    return index, _signature(ledger_at(index))
                weights, used = boundary_weights(xs, signatures, probe, step, bisections)
                counts["bisections"] += used["bisections"]
                counts["thirds"] += used["thirds"]
                weights = {index: value * weight for index, value in weights.items()}
                if abs(sum(weights.values()) - wanted) > 1.0e-12 * wanted:
                    raise FreeboardError("每行每轨迹位的权重之和不是 1/N_φ")
            else:
                ledger_at(0)
                weights = {0: float(single[2])}
            plan[(kind, where)] = weights
        plans.append({"nodes": nodes, "plan": plan, "cache": cache})
        for position, index in enumerate(base):
            lookup[(round(phi, 15), round(nodes[index]["diameter"], 15))] = (len(plans) - 1, position)
    # 第二阶段：按取样点顺序累加；签名处处相同时与现行逐位相同。ADR 0058 第 1 项。
    for diameter, phi, _weight in entries:
        row, node = lookup[(round(float(phi), 15), round(float(diameter), 15))]
        plan, cache, terminal = plans[row]["plan"], plans[row]["cache"], plans[row]["nodes"][node]["terminal"]
        touched: dict[int, float] = {}
        for kind, where, share in positions:
            value, ledger = plan[(kind, where)].get(node, 0.0), cache[(kind, where, node)]
            if kind == "bed":
                mass = value * share
                if ledger.fate == "return":
                    returned += mass
                elif ledger.fate == "elutriate":
                    elutriated += mass
                else:
                    captured += mass
                    if ledger.entered_cell is None:
                        raise FreeboardError("捕获轨迹没有进入的格子")
                    totals["c_down" if ledger.entry_kind == "down" else "c_up"][ledger.entry_cell] += mass
                    _touch_max(touched, ledger.entered_cell, value)
                    captured_phi[ledger.entered_cell][phi] = captured_phi[ledger.entered_cell].get(phi, 0.0) + mass
                continue
            if ledger.fate == "return":
                local_return += value
            elif ledger.fate == "elutriate":
                local_elutriation += value
            if ledger.fate == "capture" and ledger.entered_cell is not None:
                captured_phi[ledger.entered_cell][phi] = captured_phi[ledger.entered_cell].get(phi, 0.0) + value
            entered = ledger.captured_at_start if region == "cell" else ledger.entry_count == 1
            if ledger.fate == "capture" and entered and ((ledger.entry_kind == "up" and ledger.entered_cell == where)
                                                         or (ledger.entry_kind == "down" and ledger.entered_cell == where - 1)):
                (theta_up if ledger.entry_kind == "up" else theta_down)[where] += value
                _touch_max(touched, where, value)
        for cell in range(count):
            if cell in touched:
                totals["w_t"][cell] += touched[cell]
                totals["slip_t"][cell] += touched[cell] * terminal[cell]
        for name, values in zip(("res", "slip", "up", "down"), _sum_ledgers(
                [(cache[("bed", j, node)], plan[("bed", j)].get(node, 0.0) * float(shares[j])) for j in range(speeds_in.size)], count, 1.0)):
            for cell in range(count):
                totals[name][cell] += values[cell]
        for name, values in zip(("res", "slip", "up", "down"), _sum_ledgers(
                [(cache[("local", c, node)], plan[("local", c)].get(node, 0.0)) for c in range(count)], count, 1.0)):
            for cell in range(count):
                local_totals[name][cell] += values[cell]
    cells = []
    for cell in range(count):
        bed_mass = totals["res"][cell]
        use_bed = bed_mass > 0.0 or totals["up"][cell] + totals["down"][cell] > 0.0
        source = totals if use_bed else local_totals
        mass = source["res"][cell]
        outflow = source["up"][cell] + source["down"][cell]
        if mass == 0.0:
            if outflow != 0.0: raise FreeboardError("自由段有穿出通量但停留时间为 0")
            cells.append(CellTransport(theta_up[cell], theta_down[cell], _FREE_ABSENT, None, None, None, None, None, **_window_factors(windows, captured_phi[cell], cell, first_class)))
            continue
        hover = None if totals["w_t"][cell] == 0.0 else totals["slip_t"][cell] / totals["w_t"][cell]
        cells.append(CellTransport(
            theta_up[cell], theta_down[cell], "存在",
            source["up"][cell] / mass, source["down"][cell] / mass,
            totals["c_up"][cell] / mass, totals["c_down"][cell] / mass,
            mass, source["slip"][cell] / mass, hover,
            **_window_factors(windows, captured_phi[cell], cell, first_class),
        ))
    return ReferenceAccounting(
        tuple(item.height_m for item in first),
        tuple(sorted({cell for item in first for cell in item.cells})),
        returned, elutriated, captured, float(shares.sum()), 0.0,
        local_return, local_elutriation, tuple(cells),
        tuple((item.cells, item.height_m, item.half_width_m) for item in first),
        max((item[0] for item in statistics), default=0),
        max((item[1] for item in statistics), default=0.0),
        float(np.percentile([item[1] for item in statistics], 99)) if statistics else 0.0,
        max(statistics, key=lambda item: item[1]) if statistics else None,
        tuple(counts.values()),
    )
