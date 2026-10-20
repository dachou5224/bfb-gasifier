"""自由空间单颗粒的阻力、运动方程与数值积分。"""

from __future__ import annotations

import math
import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard._common import G, FreeboardError, _finite, _positive

# 式 (3.86) 取 Re_∞ → ∞。解释性选择：ADR 0016 第 8 条。
_Z_PREFACTOR = 73.69
_Z_EXPONENT = 5.0784


def drag_limit_z(phi_s: float) -> float:
    """z = 73,69·exp(−5,0784·φ_S)。

    解释性选择：ADR 0016 第 8 条。
    Source: Hamel Eq.3.86, S. 45
    """
    sphericity = _finite("phi_s", phi_s)
    if not 0.0 < sphericity <= 1.0:
        raise FreeboardError("φ_S 必须在 (0, 1] 内")
    return _Z_PREFACTOR * math.exp(-_Z_EXPONENT * sphericity)


def drag_groups(rho_g: float, mu_g: float, diameter_m: float, rho_p: float, phi_s: float) -> tuple[float, float, float]:
    """Stokes 系数、二次阻力系数、有效重力。

    线性项是 18 μ_g / (ρ_P d_P²)，二次项是 3 ρ_g z / (4 ρ_P d_P)。
    Source: Hamel Eq.3.85, S. 45
    Source: Hamel Eq.3.82, S. 44
    """
    density_g = _positive("rho_g", rho_g)
    viscosity = _positive("mu_g", mu_g)
    diameter = _positive("diameter_m", diameter_m)
    density_p = _positive("rho_p", rho_p)
    stokes = 18.0 * viscosity / (density_p * diameter**2)
    quadratic = 3.0 * density_g * drag_limit_z(phi_s) / (4.0 * density_p * diameter)
    weight = (density_p - density_g) * G / density_p
    if weight <= 0.0:
        raise FreeboardError("ρ_P 必须大于 ρ_g，否则终端沉降速度没有定义")
    return stokes, quadratic, weight


def drag_coefficient(reynolds: float, phi_s: float) -> float:
    """c_W = 24/Re + z。Re = 0 时 Stokes 项单独没有定义，调用方用加速度式。

    Source: Hamel Eq.3.85, S. 45
    """
    re = _positive("reynolds", reynolds)
    return 24.0 / re + drag_limit_z(phi_s)


def particle_acceleration(
    particle_velocity_m_s: float,
    gas_velocity_m_s: float,
    rho_g: float,
    mu_g: float,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
) -> float:
    """du_P/dt。相对速度换号时阻力方向跟着变，不另设分支。

    Source: Hamel Eq.3.82, S. 44
    Source: Hamel Eq.3.83, S. 44
    Source: Hamel Eq.3.85, S. 45
    """
    velocity = _finite("particle_velocity_m_s", particle_velocity_m_s)
    gas = _finite("gas_velocity_m_s", gas_velocity_m_s)
    stokes, quadratic, weight = drag_groups(rho_g, mu_g, diameter_m, rho_p, phi_s)
    relative = velocity - gas
    return -stokes * relative - quadratic * relative * abs(relative) - weight


def terminal_settling_velocity(rho_g: float, mu_g: float, diameter_m: float, rho_p: float, phi_s: float) -> float:
    """v_t > 0。加速度为零、相对速度为 −v_t 时的沉降速度。

    Source: Hamel Eq.3.82, S. 44
    Source: Hamel Eq.3.85, S. 45
    """
    stokes, quadratic, weight = drag_groups(rho_g, mu_g, diameter_m, rho_p, phi_s)
    if quadratic == 0.0:
        return weight / stokes
    return (-stokes + math.sqrt(stokes**2 + 4.0 * quadratic * weight)) / (2.0 * quadratic)


def motion_coefficients(
    gas_velocity_m_s: float,
    faster_than_gas: bool,
    rho_g: float,
    mu_g: float,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
) -> tuple[float, float, float]:
    """式 (3.88) 的 A、B、C。颗粒快于气体时用括号外的符号，否则用括号内的符号。

    Source: Hamel Eq.3.88, S. 46
    Source: Hamel Eq.3.89, S. 46
    Source: Hamel Eq.3.90, S. 46
    Source: Hamel Eq.3.91, S. 46
    """
    gas = _finite("gas_velocity_m_s", gas_velocity_m_s)
    if not isinstance(faster_than_gas, bool):
        raise FreeboardError("faster_than_gas 必须是布尔值")
    stokes, quadratic, weight = drag_groups(rho_g, mu_g, diameter_m, rho_p, phi_s)
    if faster_than_gas:
        linear = -stokes + 2.0 * quadratic * gas
        constant = stokes * gas - quadratic * gas**2 - weight
        return -quadratic, linear, constant
    linear = -stokes - 2.0 * quadratic * gas
    constant = stokes * gas + quadratic * gas**2 - weight
    return quadratic, linear, constant


def discriminant(a: float, b: float, c: float) -> float:
    """Δ = 4AC − B²。

    Source: Hamel Eq.3.92, S. 47
    """
    return 4.0 * _finite("a", a) * _finite("c", c) - _finite("b", b) ** 2


def closed_form_state(
    time_s: float,
    a: float,
    b: float,
    c: float,
    initial_velocity_m_s: float,
    initial_height_m: float,
    initial_time_s: float = 0.0,
) -> tuple[float, float]:
    """常系数且 Δ > 0 时的 u_P(t)、h_P(t)。只作测试基准，不用于求解。

    Source: Hamel Eq.3.93, S. 47
    Source: Hamel Eq.3.94, S. 48
    Source: Hamel Eq.3.95, S. 48
    """
    time = _finite("time_s", time_s)
    coefficient_a = _finite("a", a)
    coefficient_b = _finite("b", b)
    _finite("c", c)
    initial_velocity = _finite("initial_velocity_m_s", initial_velocity_m_s)
    initial_height = _finite("initial_height_m", initial_height_m)
    initial_time = _finite("initial_time_s", initial_time_s)
    if coefficient_a == 0.0:
        raise FreeboardError("A = 0 时式 (3.93) 没有定义")
    delta = discriminant(coefficient_a, coefficient_b, c)
    if delta <= 0.0:
        raise FreeboardError("式 (3.93) 只用于 Δ > 0")
    root = math.sqrt(delta)
    argument = (2.0 * coefficient_a * initial_velocity + coefficient_b) / root
    omega = 2.0 / root * math.atan(argument)
    elapsed = time - initial_time
    phase = 0.5 * omega * root + 0.5 * elapsed * root
    velocity = root / (2.0 * coefficient_a) * math.tan(phase) - coefficient_b / (2.0 * coefficient_a)
    cosine_now = math.cos(-0.5 * omega * root - 0.5 * elapsed * root)
    cosine_start = math.cos(0.5 * omega * root)
    if cosine_now <= 0.0 or cosine_start <= 0.0:
        raise FreeboardError("式 (3.95) 的余弦不是正数，这个时刻不在印刷式的实数域内")
    height = (
        initial_height
        - math.log(cosine_now) / coefficient_a
        - coefficient_b / (2.0 * coefficient_a) * elapsed
        + math.log(cosine_start) / coefficient_a
    )
    return velocity, height


def _forcing(velocity: float, gas: float, stokes: float, quadratic: float, weight: float) -> float:
    """u′ = −Stokes·u + f。f 含气体拖曳、二次阻力和重力。"""
    relative = velocity - gas
    return stokes * gas - quadratic * relative * abs(relative) - weight


def _advance(velocity: float, forcing_0: float, forcing_1: float, stokes: float, step_s: float) -> float:
    """线性阻力按梯形外力精确积分。外力在步内线性，步长趋于 0 时回到式 (3.82)。"""
    product = stokes * step_s
    if product < 1.0e-6:
        slope = -stokes * velocity + forcing_0
        average = 0.5 * (forcing_0 + forcing_1)
        return velocity + step_s * (-stokes * velocity + average) - 0.5 * stokes * slope * step_s * step_s
    decay = math.exp(-product)
    integral = forcing_0 * (1.0 - decay) / stokes + (forcing_1 - forcing_0) * (
        1.0 / stokes - (1.0 - decay) / (stokes * stokes * step_s)
    )
    return decay * velocity + integral


def _step_velocity(
    height: float,
    velocity: float,
    step_s: float,
    gas_at,
    stokes: float,
    quadratic: float,
    weight: float,
) -> tuple[float, float]:
    """预测－校正。气速和二次项都取步两端，线性项仍精确积分。"""
    gas_0 = float(gas_at(height))
    forcing_0 = _forcing(velocity, gas_0, stokes, quadratic, weight)
    predicted = _advance(velocity, forcing_0, forcing_0, stokes, step_s)
    gas_1 = float(gas_at(height + 0.5 * (velocity + predicted) * step_s))
    forcing_1 = _forcing(predicted, gas_1, stokes, quadratic, weight)
    updated = _advance(velocity, forcing_0, forcing_1, stokes, step_s)
    return height + 0.5 * (velocity + updated) * step_s, updated


def integrate_motion(
    height_0: float,
    velocity_0: float,
    duration_s: float,
    gas_velocity_m_s: float,
    rho_g: float,
    mu_g: float,
    diameter_m: float,
    rho_p: float,
    phi_s: float,
    step_s: float,
) -> tuple[float, float]:
    """给定常数气速，把式 (3.82) 积分到指定时刻。气速不随高度变化。

    Source: Hamel Eq.3.82, S. 44
    """
    height = _finite("height_0", height_0)
    velocity = _finite("velocity_0", velocity_0)
    duration = _positive("duration_s", duration_s)
    gas = _finite("gas_velocity_m_s", gas_velocity_m_s)
    step = _positive("step_s", step_s)
    diameter = _positive("diameter_m", diameter_m)
    density_p = _positive("rho_p", rho_p)
    stokes, quadratic, weight = drag_groups(rho_g, mu_g, diameter, density_p, phi_s)
    elapsed = 0.0
    while elapsed < duration:
        piece = min(step, duration - elapsed)
        height, velocity = _step_velocity(height, velocity, piece, lambda _height: gas, stokes, quadratic, weight)
        elapsed += piece
    return height, velocity



def integrate_freeboard(
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
) -> tuple[str, float, float, float]:
    """沿各格当地气速积分式 (3.82)，直到回到床面或飞出顶部。

    返回归宿、飞行时间、高度和速度。归宿只有回落和带出。
    时间上限到了仍未离开，就报错，不改判。
    Source: Hamel Eq.3.82, S. 44
    Source: Hamel Eq.3.88, S. 46
    """
    edges = np.asarray(edges_m, dtype=float)
    densities = np.asarray(rho_g, dtype=float)
    viscosities = np.asarray(mu_g, dtype=float)
    cell_count = edges.size - 1
    if densities.shape != (cell_count,) or viscosities.shape != (cell_count,):
        raise FreeboardError("气体密度和粘度必须每个格一个数")
    if cell_count < 1 or np.any(np.diff(edges) <= 0.0):
        raise FreeboardError("自由空间边界必须严格递增")
    height = _finite("height_0", height_0)
    velocity = _finite("velocity_0", velocity_0)
    if not edges[0] <= height <= edges[-1]:
        raise FreeboardError("起步高度不在自由空间内")
    step_cap = _positive("freeboard_max_step_s", config.freeboard_max_step_s)
    time_limit = _positive("freeboard_time_limit_s", config.freeboard_time_limit_s)
    diameter = _positive("diameter_m", diameter_m)
    density_p = _positive("rho_p", rho_p)
    groups = [
        drag_groups(float(densities[index]), float(viscosities[index]), diameter, density_p, phi_s)
        for index in range(cell_count)
    ]
    cell = int(np.searchsorted(edges, height, side="right") - 1)
    cell = min(max(cell, 0), cell_count - 1)
    elapsed = 0.0

    def local_acceleration(index: int, height_now: float, velocity_now: float) -> float:
        gas_now = float(gas_velocity(index, height_now))
        stokes_now, quadratic_now, weight_now = groups[index]
        relative = velocity_now - gas_now
        return -stokes_now * relative - quadratic_now * relative * abs(relative) - weight_now

    def resolve_face(cell_now: int) -> tuple[int, str | None]:
        """速度为零时用加速度判断是否正在离开。进入哪一格，就用那一格的物性。"""
        if height <= float(edges[0]):
            accel = local_acceleration(0, height, velocity)
            if velocity < 0.0 or (velocity == 0.0 and accel < 0.0):
                return 0, "return"
            return 0, None
        if height >= float(edges[-1]):
            last = cell_count - 1
            accel = local_acceleration(last, height, velocity)
            if velocity > 0.0 or (velocity == 0.0 and accel > 0.0):
                return last, "elutriate"
            return last, None
        face = None
        for index in range(1, cell_count):
            if height == float(edges[index]):
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

    stalled = 0
    while elapsed < time_limit:
        cell, fate = resolve_face(cell)
        if fate == "return":
            return fate, elapsed, float(edges[0]), velocity
        if fate == "elutriate":
            return fate, elapsed, float(edges[-1]), velocity
        step = min(step_cap, time_limit - elapsed)
        stokes, quadratic, weight = groups[cell]

        def gas_at(height_now, index=cell):
            return gas_velocity(index, height_now)

        height_1, velocity_1 = _step_velocity(height, velocity, step, gas_at, stokes, quadratic, weight)
        shrinks = 0
        while velocity * velocity_1 < 0.0 and step > step_cap * 1.0e-4 and shrinks < 16:
            step *= 0.5
            height_1, velocity_1 = _step_velocity(height, velocity, step, gas_at, stokes, quadratic, weight)
            shrinks += 1
        low = float(edges[cell])
        high = float(edges[cell + 1])
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
        if not math.isfinite(height) or not math.isfinite(velocity):
            raise FreeboardError("轨迹出现非有限的高度或速度")
    raise FreeboardError(f"积分时间达到上限 {time_limit} s，归宿仍未确定")
