"""自由空间参考轨迹的归宿、记账与输运系数。

模型边界选择：ADR 0026（目标所有者 2026-10-19）。
格内常系数解析推进与界面稳定点：ADR 0047（目标所有者 2026-10-28；第 1 项为解释性选择，第 2、3 项为求解手段）。
粒级内取样、半格捕获区与源格→目标格：ADR 0048（解释性选择）。
捕获区取整格、进入即在格界，只是 ADR 0048 兼容测试用的旧规则，不是生产取样。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

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


def _classify(
    count: int,
    fate: str,
    elapsed: float,
    height: float,
    velocity: float,
    reversals: int,
    start_cell: int,
    segments: list[tuple[int, float, float, float]],
    crossings: list[tuple[float, int, int | None, str]],
    midpoints: list[tuple[float, int, str]],
    capture_cells: frozenset[int],
    region: str,
    segments_used: int = 0,
    captured_face: tuple[int, ...] | None = None,
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
    for cell, start, end, piece in segments:
        if capture_time is None or end <= capture_time:
            kept = end - start
        elif start < capture_time:
            kept = capture_time - start
        else:
            continue
        residence[cell] += kept
        slip[cell] += piece * kept / (end - start)
    for time, from_cell, _to_cell, direction in crossings:
        if captured and capture_time is not None and time >= capture_time:
            continue
        if direction == "up":
            upward_free[from_cell] += 1.0
        else:
            downward_free[from_cell] += 1.0
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
    ADR 0048 第 2 项；求解手段：ADR 0047 第 2、3 项。
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
    outcome = None
    chosen = None
    span = bounds[-1] - bounds[0]
    for index in range(count):
        if abs(height - 0.5 * (bounds[index] + bounds[index + 1])) <= 1.0e-12 * span:
            midpoints.append((0.0, index, "up" if velocity > 0.0 else "down"))

    def close(end: float) -> None:
        nonlocal open_time, pending
        if end > open_time:
            segments.append((cell, open_time, end, pending))
            open_time = end
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
    return _classify(
        count, outcome, elapsed, height, velocity, reversals, start_cell, segments, crossings, midpoints,
        capture_cells, region, segments_used, None if chosen is None else chosen.cells,
    )


def _sum_ledgers(ledgers, count: int, weight: float):
    """把若干条轨迹的停留、滑移与自由状态穿出按（样本权重 × 轨迹权重）加总。"""
    residence = [0.0] * count
    slip = [0.0] * count
    upward_free = [0.0] * count
    downward_free = [0.0] * count
    for ledger, share in ledgers:
        scale = weight * share
        for cell in range(count):
            residence[cell] += scale * ledger.residence_s[cell]
            slip[cell] += scale * ledger.slip_m[cell]
            upward_free[cell] += scale * ledger.upward_free[cell]
            downward_free[cell] += scale * ledger.downward_free[cell]
    return residence, slip, upward_free, downward_free


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
) -> ReferenceAccounting:
    """由床面样本和各格当地轨迹得到 a^f、b^f、c↑、c↓、θ↑、θ↓ 与抛出拆分。

    每个粒级按 ADR 0048 第 1 项的取样点分别做参考轨迹，按权重汇总；m̂ = 0 时自由状态不存在。
    解释性选择：ADR 0016 第 9、10 条，ADR 0026 第 2、3 项，ADR 0048 第 1、2 项。
    Source: Hamel Eq.2.6, S. 17
    Source: Hamel Eq.3.82, S. 44
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
    returned, elutriated, captured = direct, 0.0, 0.0
    local_return, local_elutriation = 0.0, 0.0
    statistics: list[tuple[int, float, float, float, str]] = []
    first: tuple[StableEquilibrium, ...] = ()
    for index, (diameter, phi, weight) in enumerate(entries):
        equilibria = stable_equilibria(edges, gas_velocity, densities, viscosities, diameter, rho_p, phi)
        if index == 0:
            first = equilibria
        capture_set = frozenset(cell for item in equilibria for cell in item.cells)

        def follow(height: float, speed: float, diameter=diameter, phi=phi, equilibria=equilibria) -> TrajectoryFate:
            return trace_until_fate(
                height, speed, edges, gas_velocity, densities, viscosities, diameter, rho_p, phi,
                config, equilibria, scale, region,
            )

        bed = [(follow(bed_surface, float(speeds_in[i])), float(shares[i])) for i in range(speeds_in.size)]
        local, settling = [], []
        for cell in range(count):
            terminal = terminal_settling_velocity(float(densities[cell]), float(viscosities[cell]), diameter, rho_p, phi)
            settling.append(terminal)
            middle = 0.5 * (edges[cell] + edges[cell + 1])
            local.append(follow(middle, cell_mean_velocity(gas_velocity, cell, edges[cell], edges[cell + 1]) - terminal))
        touched = [False] * count
        for ledger, share in bed:
            mass = weight * share
            if ledger.fate == "return":
                returned += mass
            elif ledger.fate == "elutriate":
                elutriated += mass
            else:
                captured += mass
                if ledger.entered_cell is None:
                    raise FreeboardError("捕获轨迹没有进入的格子")
                key = "c_down" if ledger.entry_kind == "down" else "c_up"
                totals[key][ledger.entry_cell] += mass
                touched[ledger.entry_cell] = True
        for cell, ledger in enumerate(local):
            if ledger.fate == "return":
                local_return += weight
            elif ledger.fate == "elutriate":
                local_elutriation += weight
            entered = ledger.captured_at_start if region == "cell" else ledger.entry_count == 1
            if ledger.fate != "capture" or not entered:
                continue
            if ledger.entry_kind == "up" and ledger.entered_cell == cell:
                theta_up[cell] += weight
                touched[cell] = True
            elif ledger.entry_kind == "down" and ledger.entered_cell == cell - 1:
                theta_down[cell] += weight
                touched[cell] = True
        for cell in range(count):
            if touched[cell]:
                totals["w_t"][cell] += weight
                totals["slip_t"][cell] += weight * settling[cell]
        bed_sum = dict(zip(("res", "slip", "up", "down"), _sum_ledgers(bed, count, weight)))
        for name, values in bed_sum.items():
            for cell in range(count):
                totals[name][cell] += values[cell]
        local_sum = dict(zip(("res", "slip", "up", "down"), _sum_ledgers([(ledger, 1.0) for ledger in local], count, weight)))
        for name, values in local_sum.items():
            for cell in range(count):
                local_totals[name][cell] += values[cell]
        # ADR 0050 第 5 项：记录段数与归宿时刻，只报告。
        for ledger, start in [(item, "床面抛出") for item, _share in bed] + [(item, f"格 {cell} 中点") for cell, item in enumerate(local)]:
            statistics.append((ledger.segments, ledger.time_s, diameter, phi, start))
    cells = []
    for cell in range(count):
        bed_mass, local_mass = totals["res"][cell], local_totals["res"][cell]
        use_bed = bed_mass > 0.0 or totals["up"][cell] + totals["down"][cell] > 0.0
        source = totals if use_bed else local_totals
        mass = source["res"][cell]
        outflow = source["up"][cell] + source["down"][cell]
        if mass == 0.0:
            if outflow != 0.0:
                raise FreeboardError("自由段有穿出通量但停留时间为 0")
            cells.append(CellTransport(theta_up[cell], theta_down[cell], _FREE_ABSENT, None, None, None, None, None))
            continue
        hover = None if totals["w_t"][cell] == 0.0 else totals["slip_t"][cell] / totals["w_t"][cell]
        cells.append(CellTransport(
            theta_up[cell], theta_down[cell], "存在",
            source["up"][cell] / mass, source["down"][cell] / mass,
            totals["c_up"][cell] / mass, totals["c_down"][cell] / mass,
            mass, source["slip"][cell] / mass, hover,
        ))
    return ReferenceAccounting(
        tuple(item.height_m for item in first),
        tuple(sorted({cell for item in first for cell in item.cells})),
        returned,
        elutriated,
        captured,
        float(shares.sum()),
        0.0,
        local_return,
        local_elutriation,
        tuple(cells),
        tuple((item.cells, item.height_m, item.half_width_m) for item in first),
        max((item[0] for item in statistics), default=0),
        max((item[1] for item in statistics), default=0.0),
        float(np.percentile([item[1] for item in statistics], 99)) if statistics else 0.0,
        max(statistics, key=lambda item: item[1]) if statistics else None,
    )
