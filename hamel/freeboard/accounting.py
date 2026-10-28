"""自由空间参考轨迹的归宿、记账与输运系数。

模型边界选择：ADR 0026（目标所有者 2026-10-19）。
格内常系数解析推进与界面稳定点：ADR 0047（目标所有者 2026-10-28；第 1 项为解释性选择，第 2、3 项为求解手段）。
格内根搜索、势能井、逐根 δ_h 与数值积分已按 ADR 0047 第 5 项撤销。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard._common import FreeboardError, _finite, _positive
from hamel.freeboard.motion import (
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


@dataclass(frozen=True)
class StableEquilibrium:
    """一个稳定界面和它的捕获区：R* 为相邻两格之并。"""

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
    # 捕获前 ∫|u_g−u_P|dt，与停留时间同一段。解释性选择：ADR 0035。
    slip_m: tuple[float, ...]
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
    mean_slip_m_s: float | None = None


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
    roots: tuple = ()


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
    """稳定点只在界面：下格 ū_g > v_t、上格 ū_g < v_t；R* 为两格之并。

    格内气速为常数时格内没有平衡点。v_t 用各自格的气体物性。
    求解手段：ADR 0047 第 1、3 项（细化 ADR 0026 第 2 项）。
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
    slip = [0.0] * count
    upward_free = [0.0] * count
    upward_hover = [0.0] * count
    downward_free = [0.0] * count
    downward_hover = [0.0] * count
    for cell, start, end, piece in segments:
        if capture_time is None or end <= capture_time:
            kept = end - start
        elif start < capture_time:
            kept = capture_time - start
        else:
            continue
        residence[cell] += kept
        slip[cell] += piece * kept / (end - start)
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
        tuple(slip),
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
    equilibria: tuple[StableEquilibrium, ...] | None = None,
    scale: float = 1.0,
) -> TrajectoryFate:
    """逐格常系数解析推进到回落、带出，或穿过稳定界面时判为捕获。

    格内气速取 u_gb(h) 的高度平均。穿过稳定界面时 |u_P| ≤ δ_u，且按闭式 h(u)
    求得的越界距离 ≤ δ_h，即判为捕获；捕获时刻按最后一次进入 R* 记账。
    解释性选择：ADR 0026 第 2 项、ADR 0047 第 1 项；求解手段：ADR 0047 第 2、3 项。
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
    time_limit = _positive("freeboard_time_limit_s", config.freeboard_time_limit_s)
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
    segments: list[tuple[int, float, float, float]] = []
    crossings: list[tuple[float, int, int | None, str]] = []
    pending = 0.0
    outcome = None
    chosen = None

    def close(end: float) -> None:
        nonlocal open_time, pending
        if end > open_time:
            segments.append((cell, open_time, end, pending))
            open_time = end
        pending = 0.0

    while elapsed < time_limit:
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
            height, velocity, bounds[cell], bounds[cell + 1], a, b, c, speeds[cell], time_limit - elapsed,
        )
        if kind == "limit":
            break
        # 段内 u 不穿过 ū_g 与 0，故 ∫|u_g−u_P|dτ = |ū_g·τ − Δh|，按段解析。
        pending += abs(speeds[cell] * tau - (height_1 - height))
        if kind == "zero" or velocity * velocity_1 < 0.0:
            reversals += 1
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
        raise FreeboardError(f"积分时间达到上限 {time_limit} s，归宿仍未确定")
    capture_cells = frozenset(chosen.cells) if chosen is not None else frozenset()
    return _classify(count, outcome, elapsed, height, velocity, reversals, start_cell, segments, crossings, capture_cells)


def _sum_ledgers(ledgers: list[tuple[TrajectoryFate, float]], count: int):
    residence = [0.0] * count
    slip = [0.0] * count
    upward_free = [0.0] * count
    upward_hover = [0.0] * count
    downward_free = [0.0] * count
    downward_hover = [0.0] * count
    for ledger, weight in ledgers:
        for cell in range(count):
            residence[cell] += weight * ledger.residence_s[cell]
            slip[cell] += weight * ledger.slip_m[cell]
            upward_free[cell] += weight * ledger.upward_free[cell]
            upward_hover[cell] += weight * ledger.upward_hover[cell]
            downward_free[cell] += weight * ledger.downward_free[cell]
            downward_hover[cell] += weight * ledger.downward_hover[cell]
    return residence, slip, upward_free, upward_hover, downward_free, downward_hover


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
) -> ReferenceAccounting:
    """由床面样本和各格当地轨迹得到 a^f、a^t、b^f、b^t、抛出拆分和 θ。

    床面样本在该格有自由段时用床面样本，否则用全部当地轨迹在该格的自由段。
    m̂ = 0 时自由状态不存在，不给系数。解释性选择：ADR 0016 第 9、10 条，ADR 0026 第 3 项。
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
    equilibria = stable_equilibria(edges, gas_velocity, densities, viscosities, diameter_m, rho_p, phi_s)
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
            scale,
        )

    bed = [(follow(bed_surface, float(speeds_in[index])), float(shares[index])) for index in range(speeds_in.size)]
    local = []
    for cell in range(count):
        midpoint = 0.5 * (float(edges[cell]) + float(edges[cell + 1]))
        settling = terminal_settling_velocity(
            float(densities[cell]), float(viscosities[cell]), diameter_m, rho_p, phi_s
        )
        mean = cell_mean_velocity(gas_velocity, cell, float(edges[cell]), float(edges[cell + 1]))
        local.append(follow(midpoint, mean - settling))
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
    bed_residence, bed_slip, bed_up_f, bed_up_t, bed_down_f, bed_down_t = _sum_ledgers(bed, count)
    local_weighted = [(ledger, 1.0) for ledger in local]
    loc_residence, loc_slip, loc_up_f, loc_up_t, loc_down_f, loc_down_t = _sum_ledgers(local_weighted, count)
    cells = []
    for cell in range(count):
        theta = 1 if local[cell].captured_at_start else 0
        use_bed = bed_residence[cell] > 0.0 or bed_up_f[cell] + bed_up_t[cell] + bed_down_f[cell] + bed_down_t[cell] > 0.0
        if use_bed:
            mass, moved, up_f, up_t, down_f, down_t = (
                bed_residence[cell], bed_slip[cell], bed_up_f[cell], bed_up_t[cell], bed_down_f[cell], bed_down_t[cell],
            )
        else:
            mass, moved, up_f, up_t, down_f, down_t = (
                loc_residence[cell], loc_slip[cell], loc_up_f[cell], loc_up_t[cell], loc_down_f[cell], loc_down_t[cell],
            )
        outflow = up_f + up_t + down_f + down_t
        if mass == 0.0:
            if outflow != 0.0:
                raise FreeboardError("自由段有穿出通量但停留时间为 0")
            cells.append(CellTransport(theta, _FREE_ABSENT, None, None, None, None, None))
            continue
        cells.append(CellTransport(theta, "存在", up_f / mass, up_t / mass, down_f / mass, down_t / mass, mass, moved / mass))
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
        tuple((item.cells, item.height_m, item.half_width_m) for item in equilibria),
    )
