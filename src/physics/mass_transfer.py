"""相间质量传递系数。

Source of truth:
- ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``
- Hamel (1999) Eq. 3.44, 3.50

Note
----
``03_hydrodynamics_core_chain.md`` 对 `(3.50)` 首项的 OCR 转录保留了量纲警告。
当前实现采用项目已锁定的一致性口径：

    K_bd = 3*u_br/(2*d_b) + sqrt(144*D_g*eps_mf*u_b/(pi*d_b^3))

其中 ``u_br`` 严格按 Eq. 3.44 含压力修正计算。
"""

from __future__ import annotations

import math

from src.core.constants import P0_HAMEL, n_b


def calc_u_br(u_d: float, P: float, n_b_factor: float = n_b, p_ref: float = P0_HAMEL) -> float:
    """Gleichung 3.44：论文符号 u_{b,r}（bubble through-flow），此处为 u_br。

    u_{b,r} = n_b * u_d * (P / P0)^(-0.15)

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.44
    """
    assert u_d >= 0.0, "u_d 必须非负 [m/s]"
    assert P > 0.0 and p_ref > 0.0, "P 与 P0 必须大于 0 [Pa]"
    assert n_b_factor > 0.0, "n_b 必须大于 0"
    return n_b_factor * u_d * ((P / p_ref) ** -0.15)


def calc_kbd(
    u_br: float,
    d_b: float,
    D_g: float,
    eps_mf: float,
    u_b: float,
) -> float:
    """计算气泡–悬浮相总传质系数 K_bd [1/s]（Hamel Gl. 3.50, p. 35）。

    K_bd = 3*u_{b,r}/(2*d_b) + sqrt(144*D_g*eps_mf*u_b/(pi*d_b^3))

    对流项：穿流（u_br 即论文符号 u_{b,r}，由 Gl. 3.44 含压力修正）。
    扩散项：Sit & Grace (1981) 渗透理论形式。

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.50
    """
    assert u_br >= 0.0, "u_br 必须非负 [m/s]"
    assert d_b > 0.0, "d_b 必须大于 0 [m]"
    assert D_g > 0.0, "D_g 必须大于 0 [m^2/s]"
    assert 0.0 < eps_mf < 1.0, "eps_mf 必须在 (0,1)"
    assert u_b >= 0.0, "u_b 必须非负 [m/s]"

    conv = 3.0 * u_br / (2.0 * d_b)
    diff = math.sqrt((144.0 * D_g * eps_mf * u_b) / (math.pi * (d_b**3)))
    return conv + diff
