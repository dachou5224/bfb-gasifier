"""自由空间单颗粒的阻力、运动方程与逐格常系数解析推进。"""

from __future__ import annotations

import math

import numpy as np

from hamel.freeboard._common import G, FreeboardError, _finite, _positive

# 式 (3.86) 取 Re_∞ → ∞。解释性选择：ADR 0016 第 8 条。
_Z_PREFACTOR = 73.69
_Z_EXPONENT = 5.0784
# |Δ| 不超过这个倍数乘 B² 时按 Δ = 0 的极限处理。求解手段：ADR 0047 第 2 项。
_NEAR_ZERO_DELTA = 1.0e-10
# 格内高度平均的高斯－勒让德节点与权重。解释性选择：ADR 0047 第 1 项。
_QUAD_X, _QUAD_W = np.polynomial.legendre.leggauss(16)


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


def cell_mean_velocity(gas_velocity, cell: int, left: float, right: float) -> float:
    """ū_g,i = 该格内 u_gb(h) 沿高度的平均值。

    解释性选择：ADR 0047 第 1 项。原文只写逐格应用，未说明格内取哪个值。
    Source: Hamel Eq.3.74, S. 43
    """
    left = _finite("cell_left_m", left)
    right = _finite("cell_right_m", right)
    if right <= left:
        raise FreeboardError("格边界必须严格递增")
    middle = 0.5 * (left + right)
    half = 0.5 * (right - left)
    total = sum(weight * float(gas_velocity(cell, middle + half * node)) for node, weight in zip(_QUAD_X, _QUAD_W))
    return float(0.5 * total)


def cell_coefficients(cell_gas_velocity: float, velocity: float, rho_g: float, mu_g: float, diameter_m: float, rho_p: float, phi_s: float) -> tuple[float, float, float]:
    """按颗粒是否快于该格气速选括号内外的符号。

    求解手段：ADR 0047 第 2 项。
    Source: Hamel Eq.3.89, S. 46
    """
    gas = _finite("cell_gas_velocity_m_s", cell_gas_velocity)
    return motion_coefficients(gas, bool(_finite("velocity", velocity) > gas), rho_g, mu_g, diameter_m, rho_p, phi_s)


def _is_near_zero_delta(a: float, b: float, c: float) -> bool:
    return abs(discriminant(a, b, c)) <= _NEAR_ZERO_DELTA * b * b


def cell_state(time_s: float, a: float, b: float, c: float, initial_velocity_m_s: float, initial_height_m: float) -> tuple[float, float]:
    """格内常系数下 u_P(τ)、h_P(τ) 的三种闭式。

    Δ > 0 用式 (3.93)–(3.95)；Δ < 0 与 Δ = 0 用自推的闭式（印刷式 3.98、3.99、3.103 有误）。
    求解手段：ADR 0047 第 2 项。
    Source: Hamel Eq.3.93, S. 47
    Source: Hamel Eq.3.95, S. 48
    """
    tau = _finite("time_s", time_s)
    coefficient_a = _finite("a", a)
    coefficient_b = _finite("b", b)
    coefficient_c = _finite("c", c)
    u0 = _finite("initial_velocity_m_s", initial_velocity_m_s)
    h0 = _finite("initial_height_m", initial_height_m)
    if coefficient_a == 0.0:
        raise FreeboardError("A = 0 时闭式没有定义")
    if _is_near_zero_delta(coefficient_a, coefficient_b, coefficient_c):
        rest = -coefficient_b / (2.0 * coefficient_a)
        drift = u0 - rest
        factor = 1.0 - coefficient_a * drift * tau
        if factor == 0.0:
            raise FreeboardError("Δ = 0 的闭式在这一点没有定义")
        return rest + drift / factor, h0 + rest * tau - math.log(abs(factor)) / coefficient_a
    delta = discriminant(coefficient_a, coefficient_b, coefficient_c)
    if delta > 0.0:
        root = math.sqrt(delta)
        theta = math.atan((2.0 * coefficient_a * u0 + coefficient_b) / root)
        angle = theta + 0.5 * root * tau
        cosine = math.cos(angle)
        if cosine <= 0.0:
            raise FreeboardError("式 (3.95) 的余弦不是正数，这个时刻不在实数域内")
        velocity = (root * math.tan(angle) - coefficient_b) / (2.0 * coefficient_a)
        height = h0 - math.log(cosine / math.cos(theta)) / coefficient_a - coefficient_b * tau / (2.0 * coefficient_a)
        return velocity, height
    rate = math.sqrt(-delta)
    psi = 2.0 * coefficient_a * u0 + coefficient_b
    ratio = (psi + rate) / (psi - rate)
    factor = ratio * math.exp(-rate * tau)
    if factor == 1.0:
        raise FreeboardError("Δ < 0 的闭式在这一点没有定义")
    velocity = (rate * (factor + 1.0) / (factor - 1.0) - coefficient_b) / (2.0 * coefficient_a)
    # 把 E₀ − e^{sτ} 拆成 −e^{sτ}(1 − E₀e^{−sτ})，避免大 τ 时指数溢出（同一式）。
    height = h0 + (
        -(rate + coefficient_b) * tau
        - 2.0 * math.log(abs(1.0 - ratio * math.exp(-rate * tau)))
        + 2.0 * math.log(abs(ratio - 1.0))
    ) / (2.0 * coefficient_a)
    return velocity, height


def time_to_velocity(a: float, b: float, c: float, initial_velocity_m_s: float, target_velocity_m_s: float) -> float | None:
    """格内速度由 u₀ 到 u 的解析时刻 τ(u) = ∫du/Q。同一分支内单调，无解时返回 None。

    求解手段：ADR 0047 第 2 项。
    Source: Hamel Eq.3.88, S. 46
    """
    coefficient_a = _finite("a", a)
    coefficient_b = _finite("b", b)
    coefficient_c = _finite("c", c)
    u0 = _finite("initial_velocity_m_s", initial_velocity_m_s)
    target = _finite("target_velocity_m_s", target_velocity_m_s)
    if target == u0:
        return 0.0
    if coefficient_a == 0.0:
        return None
    if _is_near_zero_delta(coefficient_a, coefficient_b, coefficient_c):
        rest = -coefficient_b / (2.0 * coefficient_a)
        drift = u0 - rest
        wanted = target - rest
        if drift == 0.0 or wanted == 0.0:
            return None
        tau = (wanted - drift) / (coefficient_a * drift * wanted)
        return tau if tau > 0.0 else None
    delta = discriminant(coefficient_a, coefficient_b, coefficient_c)
    if delta > 0.0:
        root = math.sqrt(delta)
        theta = math.atan((2.0 * coefficient_a * u0 + coefficient_b) / root)
        wanted = math.atan((2.0 * coefficient_a * target + coefficient_b) / root)
        if wanted <= theta:
            return None
        return 2.0 * (wanted - theta) / root
    rate = math.sqrt(-delta)
    psi0 = 2.0 * coefficient_a * u0 + coefficient_b
    psi = 2.0 * coefficient_a * target + coefficient_b
    if psi0 == rate or psi == rate or psi0 == -rate or psi == -rate:
        return None
    ratio = ((psi0 + rate) / (psi0 - rate)) / ((psi + rate) / (psi - rate))
    if ratio <= 0.0:
        return None
    tau = math.log(ratio) / rate
    return tau if tau > 0.0 else None


def cell_height_at_velocity(a: float, b: float, c: float, initial_velocity_m_s: float, initial_height_m: float, target_velocity_m_s: float) -> float | None:
    """闭式 h(u)：速度到 target 时的高度。无解返回 None。

    求解手段：ADR 0047 第 2、3 项。
    Source: Hamel Eq.3.88, S. 46
    """
    tau = time_to_velocity(a, b, c, initial_velocity_m_s, target_velocity_m_s)
    if tau is None:
        return None
    _velocity, height = cell_state(tau, a, b, c, initial_velocity_m_s, initial_height_m)
    return height


def _boundary_time(a: float, b: float, c: float, initial_velocity_m_s: float, initial_height_m: float, target_height_m: float, window_s: float | None = None) -> float | None:
    """在 [0, window_s]（None 为无期限）内对单调段上的 h(τ) 二分求根，绝对容差 1e-12 m。

    段内 u 不穿过 ū_g、也不变号，h(τ) 单调，故端点判定即可。求解手段：ADR 0047 第 2 项。
    """
    sign = 1.0 if target_height_m >= initial_height_m else -1.0

    def offset(tau: float) -> tuple[float, float]:
        _velocity, height = cell_state(tau, a, b, c, initial_velocity_m_s, initial_height_m)
        return (height - target_height_m) * sign, height

    if offset(0.0)[0] >= 0.0:
        return 0.0
    # 单调段上 h(τ) 单调无界，没有期限时用几何扩张找到跨越点的上界。
    if window_s is None:
        high = 1.0e-9
        for _ in range(200):
            if offset(high)[0] >= 0.0:
                break
            high *= 2.0
        else:
            return None
        low = 0.5 * high
    else:
        window = _positive("window_s", window_s)
        if offset(window)[0] < 0.0:
            return None
        low, high = 0.0, window
    for _ in range(200):
        middle = 0.5 * (low + high)
        value, height = offset(middle)
        if abs(height - target_height_m) <= 1.0e-12 or middle <= low or middle >= high:
            return middle
        if value >= 0.0:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def advance_cell(height_m: float, velocity_m_s: float, left_m: float, right_m: float, a: float, b: float, c: float, gas_velocity_m_s: float) -> tuple[str, float, float, float]:
    """格内推进到最早事件：到达格界、穿过 u_g 或转向。

    返回（事件，τ，h，u）。事件为 upper / lower / gas / zero。段内没有事件时返回
    limit，由调用方的段数保护报告。段数保护代替物理时间截止：求解器参数 ADR 0050。
    求解手段：ADR 0047 第 2 项。
    Source: Hamel Eq.3.88, S. 46
    """
    height = _finite("height_m", height_m)
    velocity = _finite("velocity_m_s", velocity_m_s)
    gas = _finite("gas_velocity_m_s", gas_velocity_m_s)
    events = []
    if velocity > gas:
        tau = time_to_velocity(a, b, c, velocity, gas)
        if tau is not None and tau > 0.0:
            events.append((tau, "gas"))
    if velocity != 0.0:
        tau = time_to_velocity(a, b, c, velocity, 0.0)
        if tau is not None and tau > 0.0:
            events.append((tau, "zero"))
    horizon = min(events)[0] if events else None
    if velocity != 0.0:
        direction = 1.0 if velocity > 0.0 else -1.0
    elif c != 0.0:
        direction = 1.0 if c > 0.0 else -1.0
    else:
        return "limit", 0.0, height, velocity
    if direction > 0.0 and height >= right_m:
        return "upper", 0.0, height, velocity
    if direction < 0.0 and height <= left_m:
        return "lower", 0.0, height, velocity
    target = right_m if direction > 0.0 else left_m
    tau = _boundary_time(a, b, c, velocity, height, target, horizon)
    if tau is not None:
        new_velocity, new_height = cell_state(tau, a, b, c, velocity, height)
        return ("upper" if direction > 0.0 else "lower"), tau, new_height, new_velocity
    if horizon is None:
        return "limit", 0.0, height, velocity
    kind = min(events)[1]
    new_velocity, new_height = cell_state(horizon, a, b, c, velocity, height)
    return kind, horizon, new_height, (0.0 if kind == "zero" else gas)
