"""相体积分率：n_RZ、床层级 dense fraction、epsilon_b。

Source: specs/02_hydrodynamics.md §3; Hamel (1999)
注意：n_RZ 必须按 Re_s 分段计算，禁止使用常数 4.65。
"""

from __future__ import annotations

import math


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
    """计算床层级 dense / suspension phase 体积分率 = 1 - epsilon_b。

    注意这不是 dense phase 内部孔隙率；后者应使用 ``calc_emulsion_porosity``。
    Source: specs/02_hydrodynamics.md §3
    """
    assert 0.0 <= epsilon_b < 1.0, "epsilon_b 必须在 [0,1)"
    return 1.0 - epsilon_b


def calc_visible_bubble_fraction(u0: float, u_b: float, u_d: float, n_b: float = 2.7) -> float:
    """按 Hamel Chapter 3 的 visible bubble 口径计算 epsilon_b。

    epsilon_b = (u0 - u_d) / (u_b + n_b*u_d - u_d)

    这里的 ``u_d`` 取 suspension / emulsion 气体特征速度，用于替代
    旧实现中直接以 ``u_mf`` 估算 excess-gas bubble hold-up 的简化式。
    """
    assert u0 >= 0.0 and u_b > 0.0 and u_d >= 0.0, "速度输入需非负且 u_b > 0"
    assert n_b > 0.0, "n_b 必须大于 0"
    denom = max(u_b + (n_b - 1.0) * u_d, 1e-12)
    eps_b = (u0 - u_d) / denom
    return min(max(eps_b, 0.0), 0.95)


def calc_emulsion_porosity(u_d: float, u_mf: float, eps_mf: float, n_rz: float) -> float:
    """按 Richardson-Zaki 关系计算 emulsion / suspension 孔隙率。

    epsilon_d,void = epsilon_mf * (u_d / u_mf)^(1 / n_RZ)

    注意这对应 Hamel 文中的乳化相孔隙率，不等同于床层内的相体积分率。
    """
    assert u_d >= 0.0 and u_mf > 0.0, "u_d/u_mf 必须满足 u_mf > 0"
    assert 0.0 < eps_mf < 1.0, "eps_mf 必须在 (0,1)"
    assert n_rz > 0.0, "n_rz 必须大于 0"
    ratio = max(u_d / u_mf, 1.0)
    eps = eps_mf * math.pow(ratio, 1.0 / n_rz)
    return min(max(eps, eps_mf), 0.99)


def calc_bulk_solid_holdup(eps_b: float, eps_d_voidage: float) -> float:
    """按 Hamel 两相层级计算整床固相体积分率 proxy。

    层级关系：
    - ``eps_b``: visible bubble volume fraction（床层级气泡体积分率）
    - ``1 - eps_b``: suspension / dense phase 占床层总体积的比例
    - ``eps_d_voidage``: suspension / dense phase 内部孔隙率

    因此床层平均空隙率为：
        eps_avg = eps_b + (1 - eps_b) * eps_d_voidage

    对应整床固相占比：
        bulk_solid = 1 - eps_avg = (1 - eps_b) * (1 - eps_d_voidage)

    这里故意不使用 ``eps_d`` 变量名，以避免与原文 ``epsilon_d``（悬浮相内部孔隙率）
    和代码中 ``eps_d = 1 - eps_b``（悬浮相总体积占比）之间的层级混淆。
    """
    assert 0.0 <= eps_b <= 1.0, "eps_b 必须在 [0,1]"
    assert 0.0 <= eps_d_voidage <= 1.0, "eps_d_voidage 必须在 [0,1]"
    return (1.0 - eps_b) * (1.0 - eps_d_voidage)
