"""自由空间参考轨迹的归宿、记账与输运系数。

模型边界选择：ADR 0026（目标所有者 2026-10-19）。
邻域只用来确定归宿。解释性选择：ADR 0033。
积分步与格边界处理同 integrate_freeboard，这里另记停留和穿越。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard._common import FreeboardError, _finite, _positive
from hamel.freeboard.motion import _step_velocity, drag_groups, terminal_settling_velocity

# 界面上 |u_g − v_t| 不超过这个值，就认为根在界面上。解释性选择：ADR 0033。
_FACE_RESIDUAL_M_S = 1.0e-8
# 与界面根相距不超过这个值的内部根并入界面根。解释性选择：ADR 0033。
_FACE_MATCH_M = 1.0e-6
_FREE_ABSENT = "自由状态不存在"


@dataclass(frozen=True)
class StableEquilibrium:
    """一个稳定平衡和它的捕获区。"""

    height_m: float
    cells: tuple[int, ...]
    half_width_m: float


@dataclass(frozen=True)
class TrajectoryFate:
    """一条轨迹的归宿，以及捕获时刻之前的停留和穿出。"""

    fate: str
    time_s: float
    height_m: float
    velocity_m_s: float
    velocity_reversals: int
    capture_time_s: float | None
    captured_at_start: bool
    entered_cell: int | None
    residence_s: tuple[float, ...]
    upward_free: tuple[float, ...]
    upward_hover: tuple[float, ...]
    downward_free: tuple[float, ...]
    downward_hover: tuple[float, ...]


@dataclass(frozen=True)
class CellTransport:
    """一格一个粒级。自由状态不存在时，四个系数和 m̂ 都是 None。"""

    theta: int
    free_state: str
    a_f: float | None
    a_t: float | None
    b_f: float | None
    b_t: float | None
    m_hat: float | None


@dataclass(frozen=True)
class ReferenceAccounting:
    """一个粒级的参考记账。份额对单位抛出，当地轨迹的悬停入流按每格权重 1。"""

    stable_heights_m: tuple[float, ...]
    capture_cells: tuple[int, ...]
    return_share: float
    elutriation_share: float
    capture_share: float
    ejection_free: float
    ejection_hover: float
    hover_inflow: tuple[float, ...]
    local_hover_inflow: tuple[float, ...]
    local_return: float
    local_elutriation: float
    cells: tuple[CellTransport, ...]


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


def _slope(gas_velocity, cell: int, height: float, left: float, right: float, side: str) -> float:
    """side 为 down / up 时用单侧差商，为 center 时在格内取中心差商。"""
    width = right - left
    dh = min(1.0e-4, width * 0.25)
    if side == "down":
        sample = max(left, height - dh)
        return (float(gas_velocity(cell, height)) - float(gas_velocity(cell, sample))) / (height - sample)
    if side == "up":
        sample = min(right, height + dh)
        return (float(gas_velocity(cell, sample)) - float(gas_velocity(cell, height))) / (sample - height)
    if height - left < dh:
        return _slope(gas_velocity, cell, height, left, right, "up")
    if right - height < dh:
        return _slope(gas_velocity, cell, height, left, right, "down")
    return (float(gas_velocity(cell, height + dh)) - float(gas_velocity(cell, height - dh))) / (2.0 * dh)


def _append_root(found: list[tuple[float, tuple[int, ...]]], height: float, cells: tuple[int, ...]) -> None:
    for index, (kept_height, kept_cells) in enumerate(found):
        if abs(height - kept_height) <= _FACE_MATCH_M:
            if len(cells) > len(kept_cells):
                found[index] = (height, cells)
            return
    found.append((height, cells))


def stable_equilibria(
    edges_m,
    gas_velocity,
    rho_g,
    mu_g,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
) -> tuple[StableEquilibrium, ...]:
    """u_g(h) = v_t 且 du_g/dh < 0 的根。界面根的捕获区是相邻两格。

    v_t 用该格的气体密度和粘度。解释性选择：ADR 0026 第 2 项；界面判据见 ADR 0033。
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
    found: list[tuple[float, tuple[int, ...]]] = []
    for face in range(1, count):
        height = float(edges[face])
        lower, upper = face - 1, face
        residual_lower = float(gas_velocity(lower, height)) - settling[lower]
        residual_upper = float(gas_velocity(upper, height)) - settling[upper]
        if abs(residual_lower) > _FACE_RESIDUAL_M_S or abs(residual_upper) > _FACE_RESIDUAL_M_S:
            continue
        slope_lower = _slope(gas_velocity, lower, height, float(edges[lower]), height, "down")
        slope_upper = _slope(gas_velocity, upper, height, height, float(edges[upper + 1]), "up")
        if slope_lower < 0.0 and slope_upper < 0.0:
            _append_root(found, height, (lower, upper))
        elif not (slope_lower > 0.0 and slope_upper > 0.0):
            raise FreeboardError("界面两侧的气速斜率符号不一致，不能判定这个根稳不稳定")
    samples = 32
    for cell in range(count):
        left = float(edges[cell])
        right = float(edges[cell + 1])
        previous_h = left
        previous_r = float(gas_velocity(cell, left)) - settling[cell]
        for step in range(1, samples + 1):
            height = left + (right - left) * step / samples
            residual = float(gas_velocity(cell, height)) - settling[cell]
            bracket = previous_r * residual < 0.0
            exact = abs(residual) <= _FACE_RESIDUAL_M_S and left < height < right
            if bracket or exact:
                lo, hi = previous_h, height
                r_lo = previous_r
                root = height
                if bracket:
                    for _ in range(60):
                        mid = 0.5 * (lo + hi)
                        r_mid = float(gas_velocity(cell, mid)) - settling[cell]
                        if r_lo * r_mid <= 0.0:
                            hi = mid
                        else:
                            lo = mid
                            r_lo = r_mid
                    root = 0.5 * (lo + hi)
                if left < root < right and _slope(gas_velocity, cell, root, left, right, "center") < 0.0:
                    _append_root(found, root, (cell,))
            previous_h, previous_r = height, residual
    equilibria = []
    for height, cells in sorted(found, key=lambda item: item[0]):
        outer_left = float(edges[cells[0]])
        outer_right = float(edges[cells[-1] + 1])
        half = 0.5 * min(height - outer_left, outer_right - height)
        if half <= 0.0:
            raise FreeboardError("稳定平衡到捕获区外边界的距离不是正数")
        equilibria.append(StableEquilibrium(height, cells, half))
    return tuple(equilibria)


def _check_delta(equilibria: tuple[StableEquilibrium, ...], config: SolverConfig) -> tuple[float, float]:
    delta_h = _positive("freeboard_delta_h_m", config.freeboard_delta_h_m)
    delta_u = _positive("freeboard_delta_u_m_s", config.freeboard_delta_u_m_s)
    step = _positive("freeboard_max_step_s", config.freeboard_max_step_s)
    if not equilibria:
        return delta_h, delta_u
    limit = min(item.half_width_m for item in equilibria)
    if delta_h > limit:
        raise FreeboardError(f"δ_h = {delta_h} m 超过 h* 到捕获区外边界距离的一半 {limit} m")
    if delta_u * step >= delta_h:
        raise FreeboardError("δ_u 与步长的乘积必须小于 δ_h，否则一步可能跨过邻域")
    return delta_h, delta_u


def _inside(height: float, velocity: float, equilibria: tuple[StableEquilibrium, ...], delta_h: float, delta_u: float) -> bool:
    return any(abs(height - item.height_m) <= delta_h and abs(velocity) <= delta_u for item in equilibria)


def _classify(
    count: int,
    fate: str,
    elapsed: float,
    height: float,
    velocity: float,
    reversals: int,
    start_cell: int,
    segments: list[tuple[int, float, float]],
    crossings: list[tuple[float, int, int | None, str]],
    capture_cells: frozenset[int],
) -> TrajectoryFate:
    entry_index = None
    entered = start_cell if start_cell in capture_cells else None
    capture_time = 0.0 if entered is not None else None
    for index, (time, from_cell, to_cell, _direction) in enumerate(crossings):
        enters = to_cell is not None and to_cell in capture_cells and from_cell not in capture_cells
        leaves = from_cell in capture_cells and (to_cell is None or to_cell not in capture_cells)
        if enters:
            entry_index = index
            entered = to_cell
            capture_time = time
        elif leaves:
            entered = None
            capture_time = None
    captured = fate == "capture"
    if captured and (capture_time is None or entered is None):
        raise FreeboardError("判定为捕获，但轨迹没有进入捕获区")
    if not captured:
        capture_time = None
        entered = None
        entry_index = None
    at_start = captured and entry_index is None
    residence = [0.0] * count
    upward_free = [0.0] * count
    upward_hover = [0.0] * count
    downward_free = [0.0] * count
    downward_hover = [0.0] * count
    for cell, start, end in segments:
        if capture_time is None:
            residence[cell] += end - start
        elif end <= capture_time:
            residence[cell] += end - start
        elif start < capture_time:
            residence[cell] += capture_time - start
    for index, (time, from_cell, _to_cell, direction) in enumerate(crossings):
        hover = captured and entry_index is not None and index == entry_index
        free = (not captured) or (capture_time is not None and time < capture_time)
        if not hover and not free:
            continue
        if direction == "up":
            (upward_hover if hover else upward_free)[from_cell] += 1.0
        else:
            (downward_hover if hover else downward_free)[from_cell] += 1.0
    return TrajectoryFate(
        fate,
        elapsed,
        height,
        velocity,
        reversals,
        capture_time if captured else None,
        at_start,
        entered if captured else None,
        tuple(residence),
        tuple(upward_free),
        tuple(upward_hover),
        tuple(downward_free),
        tuple(downward_hover),
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
    equilibria: tuple[StableEquilibrium, ...],
) -> TrajectoryFate:
    """把式 (3.82) 积到回落、带出或进入稳定平衡的邻域。

    邻域只决定归宿。捕获时刻是最后一次进入捕获区的时刻，起步就在里面且没有离开则为 0。
    时间上限到了仍未确定归宿就报错，不改判。解释性选择：ADR 0026 第 2 项，ADR 0033。
    Source: Hamel Eq.3.82, S. 44
    """
    edges, densities, viscosities = _edges_and_gas(edges_m, rho_g, mu_g)
    count = edges.size - 1
    height = _finite("height_0", height_0)
    velocity = _finite("velocity_0", velocity_0)
    if not float(edges[0]) <= height <= float(edges[-1]):
        raise FreeboardError("起步高度不在自由空间内")
    delta_h, delta_u = _check_delta(equilibria, config)
    step_cap = _positive("freeboard_max_step_s", config.freeboard_max_step_s)
    time_limit = _positive("freeboard_time_limit_s", config.freeboard_time_limit_s)
    diameter = _positive("diameter_m", diameter_m)
    density_p = _positive("rho_p", rho_p)
    groups = tuple(
        drag_groups(float(densities[index]), float(viscosities[index]), diameter, density_p, phi_s)
        for index in range(count)
    )
    bounds = tuple(float(edge) for edge in edges)
    capture_cells = frozenset(cell for item in equilibria for cell in item.cells)
    cell = int(np.searchsorted(edges, height, side="right") - 1)
    cell = min(max(cell, 0), count - 1)
    start_cell = cell
    elapsed = 0.0
    open_time = 0.0
    reversals = 0
    segments: list[tuple[int, float, float]] = []
    crossings: list[tuple[float, int, int | None, str]] = []

    def local_acceleration(index: int, height_now: float, velocity_now: float) -> float:
        gas_now = float(gas_velocity(index, height_now))
        stokes_now, quadratic_now, weight_now = groups[index]
        relative = velocity_now - gas_now
        return -stokes_now * relative - quadratic_now * relative * abs(relative) - weight_now

    def resolve_face(cell_now: int) -> tuple[int, str | None]:
        if height <= bounds[0]:
            accel = local_acceleration(0, height, velocity)
            if velocity < 0.0 or (velocity == 0.0 and accel < 0.0):
                return 0, "return"
            return 0, None
        if height >= bounds[-1]:
            last = count - 1
            accel = local_acceleration(last, height, velocity)
            if velocity > 0.0 or (velocity == 0.0 and accel > 0.0):
                return last, "elutriate"
            return last, None
        face = None
        for index in range(1, count):
            if height == bounds[index]:
                face = index
                break
        if face is None:
            return cell_now, None
        lower = face - 1
        upper = face
        if velocity > 0.0:
            return upper, None
        if velocity < 0.0:
            return lower, None
        accel_upper = local_acceleration(upper, height, 0.0)
        accel_lower = local_acceleration(lower, height, 0.0)
        leaves_upper = accel_upper < 0.0
        leaves_lower = accel_lower > 0.0
        if leaves_upper and leaves_lower:
            raise FreeboardError("格面两侧的加速度都指向对面，不能确定进入哪一格")
        if leaves_upper or (accel_upper == 0.0 and accel_lower < 0.0):
            return lower, None
        return upper, None

    def close(end: float) -> None:
        nonlocal open_time
        if end > open_time:
            segments.append((cell, open_time, end))
            open_time = end

    stalled = 0
    outcome = None
    while elapsed < time_limit:
        resolved, domain = resolve_face(cell)
        if resolved != cell:
            close(elapsed)
            direction = "up" if resolved > cell else "down"
            crossings.append((elapsed, cell, resolved, direction))
            cell = resolved
            open_time = elapsed
        if domain is not None:
            close(elapsed)
            crossings.append((elapsed, cell, None, "down" if domain == "return" else "up"))
            outcome = domain
            break
        if _inside(height, velocity, equilibria, delta_h, delta_u):
            if cell not in capture_cells:
                raise FreeboardError("进入邻域时粒子不在捕获区里")
            close(elapsed)
            outcome = "capture"
            break
        step = min(step_cap, time_limit - elapsed)
        stokes, quadratic, weight = groups[cell]

        def gas_at(height_now, index=cell):
            return gas_velocity(index, height_now)

        velocity_before = velocity
        height_1, velocity_1 = _step_velocity(height, velocity, step, gas_at, stokes, quadratic, weight)
        shrinks = 0
        while velocity * velocity_1 < 0.0 and step > step_cap * 1.0e-4 and shrinks < 16:
            step *= 0.5
            height_1, velocity_1 = _step_velocity(height, velocity, step, gas_at, stokes, quadratic, weight)
            shrinks += 1
        low = bounds[cell]
        high = bounds[cell + 1]
        boundary = None
        fraction = None
        if height_1 > high and height_1 != height:
            fraction = (high - height) / (height_1 - height)
            boundary = high if 0.0 < fraction < 1.0 else None
        elif height_1 < low and height_1 != height:
            fraction = (low - height) / (height_1 - height)
            boundary = low if 0.0 < fraction < 1.0 else None
        if fraction is not None and fraction <= 0.0:
            stalled += 1
            if stalled > 2:
                raise FreeboardError("粒子停在格边界上，这一步没有进入相邻格")
            continue
        stalled = 0
        if boundary is None:
            height, velocity = height_1, velocity_1
            elapsed += step
        else:
            used = fraction * step
            height, velocity = _step_velocity(height, velocity, used, gas_at, stokes, quadratic, weight)
            elapsed += used
            height = boundary
        if velocity_before * velocity < 0.0:
            reversals += 1
        if not math.isfinite(height) or not math.isfinite(velocity):
            raise FreeboardError("轨迹出现非有限的高度或速度")
    if outcome is None:
        raise FreeboardError(f"积分时间达到上限 {time_limit} s，归宿仍未确定")
    # 第一步 resolve 之前的格子才是起步格。跨过界面的那一步已经写进 crossings。
    return _classify(count, outcome, elapsed, height, velocity, reversals, start_cell, segments, crossings, capture_cells)


def _sum_ledgers(ledgers: list[tuple[TrajectoryFate, float]], count: int) -> tuple[list[float], list[float], list[float], list[float], list[float]]:
    residence = [0.0] * count
    upward_free = [0.0] * count
    upward_hover = [0.0] * count
    downward_free = [0.0] * count
    downward_hover = [0.0] * count
    for ledger, weight in ledgers:
        for cell in range(count):
            residence[cell] += weight * ledger.residence_s[cell]
            upward_free[cell] += weight * ledger.upward_free[cell]
            upward_hover[cell] += weight * ledger.upward_hover[cell]
            downward_free[cell] += weight * ledger.downward_free[cell]
            downward_hover[cell] += weight * ledger.downward_hover[cell]
    return residence, upward_free, upward_hover, downward_free, downward_hover


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
) -> ReferenceAccounting:
    """由床面样本和各格当地轨迹得到 a^f、a^t、b^f、b^t、抛出拆分和 θ。

    床面样本在该格有自由段时用床面样本，否则用全部当地轨迹在该格的自由段。
    m̂ = 0 时自由状态不存在，不给系数。解释性选择：ADR 0016 第 9、10 条，ADR 0026 第 3 项。
    Source: Hamel Eq.2.6, S. 17
    Source: Hamel Eq.3.82, S. 44
    """
    edges, densities, viscosities = _edges_and_gas(edges_m, rho_g, mu_g)
    count = edges.size - 1
    speeds = np.asarray(ejection_speeds, dtype=float)
    shares = np.asarray(ejection_shares, dtype=float)
    if speeds.ndim != 1 or speeds.shape != shares.shape:
        raise FreeboardError("抛出样本的速度和份额必须是同样长的一维数组")
    if speeds.size and (np.any(speeds <= 0.0) or np.any(shares < 0.0)):
        raise FreeboardError("抛出样本的初速必须为正，份额不能为负")
    direct = _finite("direct_return_share", direct_return_share)
    if direct < 0.0 or abs(float(shares.sum()) + direct - 1.0) > 1.0e-8:
        raise FreeboardError("抛出份额与直接回落之和必须为 1，直接回落不能为负")
    equilibria = stable_equilibria(edges, gas_velocity, densities, viscosities, diameter_m, rho_p, phi_s)
    _check_delta(equilibria, config)
    capture = tuple(sorted({cell for item in equilibria for cell in item.cells}))
    capture_set = frozenset(capture)
    bed_surface = float(edges[0])

    def follow(height: float, speed: float) -> TrajectoryFate:
        return trace_until_fate(
            height,
            speed,
            edges,
            gas_velocity,
            densities,
            viscosities,
            diameter_m,
            rho_p,
            phi_s,
            config,
            equilibria,
        )

    bed = [(follow(bed_surface, float(speeds[index])), float(shares[index])) for index in range(speeds.size)]
    local = []
    for cell in range(count):
        midpoint = 0.5 * (float(edges[cell]) + float(edges[cell + 1]))
        settling = terminal_settling_velocity(
            float(densities[cell]), float(viscosities[cell]), diameter_m, rho_p, phi_s
        )
        local.append(follow(midpoint, float(gas_velocity(cell, midpoint)) - settling))
    ejection_free = 0.0
    ejection_hover = 0.0
    hover_inflow = [0.0] * count
    returned = direct
    elutriated = 0.0
    for ledger, weight in bed:
        if ledger.fate == "return":
            returned += weight
            ejection_free += weight
        elif ledger.fate == "elutriate":
            elutriated += weight
            ejection_free += weight
        elif ledger.captured_at_start:
            if ledger.entered_cell != 0 or 0 not in capture_set:
                raise FreeboardError("床面样本起步即捕获时，应进入最下一格的悬停状态")
            ejection_hover += weight
            hover_inflow[0] += weight
        else:
            ejection_free += weight
            if ledger.entered_cell is None:
                raise FreeboardError("捕获轨迹没有进入的格子")
            hover_inflow[ledger.entered_cell] += weight
    local_hover = [0.0] * count
    local_return = 0.0
    local_elutriation = 0.0
    for cell, ledger in enumerate(local):
        if ledger.fate == "return":
            local_return += 1.0
        elif ledger.fate == "elutriate":
            local_elutriation += 1.0
        else:
            if ledger.entered_cell is None:
                raise FreeboardError("当地轨迹捕获后没有进入的格子")
            local_hover[ledger.entered_cell] += 1.0
    bed_residence, bed_up_f, bed_up_t, bed_down_f, bed_down_t = _sum_ledgers(bed, count)
    local_weighted = [(ledger, 1.0) for ledger in local]
    loc_residence, loc_up_f, loc_up_t, loc_down_f, loc_down_t = _sum_ledgers(local_weighted, count)
    cells = []
    for cell in range(count):
        theta = 1 if local[cell].captured_at_start else 0
        use_bed = bed_residence[cell] > 0.0 or bed_up_f[cell] + bed_up_t[cell] + bed_down_f[cell] + bed_down_t[cell] > 0.0
        if use_bed:
            mass, up_f, up_t, down_f, down_t = (
                bed_residence[cell],
                bed_up_f[cell],
                bed_up_t[cell],
                bed_down_f[cell],
                bed_down_t[cell],
            )
        else:
            mass, up_f, up_t, down_f, down_t = (
                loc_residence[cell],
                loc_up_f[cell],
                loc_up_t[cell],
                loc_down_f[cell],
                loc_down_t[cell],
            )
        outflow = up_f + up_t + down_f + down_t
        if mass == 0.0:
            if outflow != 0.0:
                raise FreeboardError("自由段有穿出通量但停留时间为 0")
            cells.append(CellTransport(theta, _FREE_ABSENT, None, None, None, None, None))
            continue
        cells.append(
            CellTransport(
                theta,
                "存在",
                up_f / mass,
                up_t / mass,
                down_f / mass,
                down_t / mass,
                mass,
            )
        )
    captured = sum(hover_inflow)
    return ReferenceAccounting(
        tuple(item.height_m for item in equilibria),
        capture,
        returned,
        elutriated,
        captured,
        ejection_free,
        ejection_hover,
        tuple(hover_inflow),
        tuple(local_hover),
        local_return,
        local_elutriation,
        tuple(cells),
    )
