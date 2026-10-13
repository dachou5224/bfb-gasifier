"""式 (3.74) 与现行两格均匀场的对照。不改认证模型的默认气场。"""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from src.physics.wirsum_1997_freeboard import (
    constant_property_boundary_event,
    drag_equilibrium_speed_m_s,
    hamel_eq386_asymptotic_z,
    particle_acceleration_cw,
    u_gb_hamel_374,
)

# 审查证据里 s=0.524 自由空间的冻结量
_U0 = 0.45507600092036804  # [m/s]
_UB = 1.3529554730444568  # [m/s]
_DB = 0.21  # [m]
_RHO = 5.515220076315504  # [kg/m³]
_RP = 1400.0  # [kg/m³]
_D = 2.0e-4  # [m]
_MU = 4.0e-5  # [Pa·s]
_H = 4.268799873299744  # [m]
_Z = float(hamel_eq386_asymptotic_z(0.8))
_U0_LAUNCH = 1.34  # [m/s]


def _slip_gas_speed() -> float:
    """使颗粒平衡速度为零的气体速度 [m/s]。"""
    u_eq = drag_equilibrium_speed_m_s(
        u_g=_U0, u_p=0.0, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
    )
    return float(_U0 - u_eq)


def _h_star() -> float:
    vt = _slip_gas_speed()
    return float(3.0 * _DB * math.log((_UB - _U0) / (vt - _U0)))


def _walk_to(height_m: float, n_cells: int, *, same_speed: bool) -> tuple[str, float]:
    """把自由空间分成 n 格。每格气速冻结在格心，或两格都用同一个表观速度。"""
    dh = _H / int(n_cells)
    target = float(height_m)
    i = 0
    h = 0.0
    u = _U0_LAUNCH
    t = 0.0
    for _visit in range(8 * int(n_cells) + 4):
        if same_speed:
            u_gas = _U0
        else:
            u_gas = u_gb_hamel_374((i + 0.5) * dh, _U0, _UB, _DB)
        local = h - i * dh
        if i * dh < target <= (i + 1) * dh:
            fate, dt, _z, _u, _reason = constant_property_boundary_event(
                u, u_g=u_gas, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
                cell_height=target - i * dh, height0=local,
            )
            return fate, float(t + dt)
        fate, dt, _z, u_end, reason = constant_property_boundary_event(
            u, u_g=u_gas, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
            cell_height=dh, height0=local,
        )
        t += float(dt)
        u = float(u_end)
        if fate == "undecided" or reason == "asymptotic_residence":
            return "undecided", float(t)
        if fate == "returns" and i == 0:
            return "returns", float(t)
        if fate == "returns":
            i -= 1
            h = (i + 1) * dh
            continue
        if fate == "reaches_next_cell":
            i += 1
            h = i * dh
            if i >= int(n_cells):
                return "reaches_next_cell", float(t)
            continue
        return fate, float(t)
    return "unresolved", float(t)


def test_stall_height_lies_inside_this_freeboard() -> None:
    h_star = _h_star()
    assert h_star == pytest.approx(4.067263686, abs=1.0e-6)
    assert h_star < _H
    u_eq_star = drag_equilibrium_speed_m_s(
        u_g=u_gb_hamel_374(h_star, _U0, _UB, _DB),
        u_p=0.0, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
    )
    u_eq_bed = drag_equilibrium_speed_m_s(
        u_g=u_gb_hamel_374(0.0, _U0, _UB, _DB),
        u_p=0.0, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
    )
    u_eq_roof = drag_equilibrium_speed_m_s(
        u_g=u_gb_hamel_374(_H, _U0, _UB, _DB),
        u_p=0.0, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
    )
    assert abs(u_eq_star) < 1.0e-8
    assert u_eq_bed > 0.5
    assert u_eq_roof < 0.0


def test_two_identical_cells_return_while_the_profile_is_still_rising() -> None:
    same = _walk_to(2.0, 2, same_speed=True)
    assert same[0] == "returns"
    assert same[1] == pytest.approx(41.336, abs=0.05)
    profile = _walk_to(2.0, 2, same_speed=False)
    assert profile[0] == "reaches_next_cell"


def test_piecewise_rise_converges_to_the_continuous_profile() -> None:
    def rhs(_t: float, y: np.ndarray) -> list[float]:
        gas = u_gb_hamel_374(float(y[0]), _U0, _UB, _DB)
        return [
            float(y[1]),
            particle_acceleration_cw(
                float(y[1]), u_g=gas, rho_g=_RHO, rho_p=_RP, d_p=_D, mu_g=_MU, z=_Z,
            ),
        ]

    def mark(_t: float, y: np.ndarray) -> float:
        return float(y[0]) - 2.0

    mark.terminal = True
    mark.direction = 1.0
    sol = solve_ivp(
        rhs, (0.0, 40.0), np.array([0.0, _U0_LAUNCH]),
        method="Radau", events=mark, rtol=1.0e-7, atol=1.0e-9,
    )
    assert sol.t_events[0].size == 1
    t_ref = float(sol.t_events[0][0])
    errors = []
    for n in (8, 16, 32, 64):
        fate, t_n = _walk_to(2.0, n, same_speed=False)
        assert fate == "reaches_next_cell"
        errors.append(abs(t_n - t_ref))
    assert errors[-1] < 0.02
    assert errors[-1] < errors[0]
