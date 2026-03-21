"""相体积分率：n_RZ、epsilon_d、epsilon_b。

Source: specs/02_hydrodynamics.md §3; Hamel (1999)
注意：n_RZ 必须按 Re_s 分段计算，禁止使用常数 4.65。
"""

from __future__ import annotations


def calc_n_rz(Re_s: float) -> float:
    """按 Re_s 分段计算 Richardson-Zaki 指数 n_RZ。

    分段经验式（按液固/气固流化通用形式）：
    - Re_s < 0.2: n = 4.65
    - 0.2 <= Re_s < 1: n = 4.4 * Re_s^(-0.03)
    - 1 <= Re_s < 500: n = 4.4 * Re_s^(-0.10)
    - Re_s >= 500: n = 2.4

    Source: specs/02_hydrodynamics.md §3（要求分段）
    """
    assert Re_s >= 0.0, "Re_s 必须非负"
    if Re_s < 0.2:
        return 4.65
    if Re_s < 1.0:
        return 4.4 * (Re_s ** -0.03)
    if Re_s < 500.0:
        return 4.4 * (Re_s ** -0.10)
    return 2.4


def calc_epsilon_b(u0: float, u_mf: float, u_b: float) -> float:
    """计算气泡相体积分率 epsilon_b。

    epsilon_b = (u0 - u_mf) / u_b

    Source: specs/02_hydrodynamics.md §2–§3
    """
    assert u0 >= 0.0 and u_mf >= 0.0 and u_b > 0.0, "速度输入需满足 u_b > 0"
    eps_b = (u0 - u_mf) / u_b
    return min(max(eps_b, 0.0), 0.95)


def calc_epsilon_d(epsilon_b: float) -> float:
    """计算悬浮相体积分率 epsilon_d = 1 - epsilon_b。

    Source: specs/02_hydrodynamics.md §3
    """
    assert 0.0 <= epsilon_b < 1.0, "epsilon_b 必须在 [0,1)"
    return 1.0 - epsilon_b
