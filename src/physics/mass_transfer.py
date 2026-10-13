"""相间质量传递系数。

Source of truth:
- ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``
- Hamel (1999) Eq. 3.44, 3.46–3.50

Note
----
Hamel PDF page 27 (``data/hamel_pages_png/page_27.png``) 给出 Eq. (3.50) 的等价形式：

    K_bd = V̇_{b,r}/V_b + sqrt(144 D_g ε_mf u_b / (π d_b³))
         = (3/2) u_{b,r}/d_b + sqrt(...)

其中 ``V_b = π d_b³/6``，故 ``V̇_{b,r}/V_b = (3/2) u_{b,r}/d_b``（量纲 [1/s]）。
分册旧转录中的 ``u_{b,r}/V_b`` 为 OCR 误识，**不得**用于实现。
"""

from __future__ import annotations

import math
from typing import Literal

from src.core.constants import P0_HAMEL, g, n_b

KbdModel = Literal["hamel_350", "kunii_levenspiel", "sit_grace_348", "preto_349"]


def calc_u_br(u_d: float, P: float, n_b_factor: float = n_b, p_ref: float = P0_HAMEL) -> float:
    """Gleichung 3.44：论文符号 u_{b,r}（bubble through-flow），此处为 u_br。

    u_{b,r} = n_b * u_d * (P / P0)^(-0.15)

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.44
    """
    assert u_d >= 0.0, "u_d 必须非负 [m/s]"
    assert P > 0.0 and p_ref > 0.0, "P 与 P0 必须大于 0 [Pa]"
    assert n_b_factor > 0.0, "n_b 必须大于 0"
    return n_b_factor * u_d * ((P / p_ref) ** -0.15)


def _kbd_diffusion_term(D_g: float, eps_mf: float, u_b: float, d_b: float) -> float:
    """Eq. (3.48)/(3.50) 第二项：Sit & Grace 扩散分量 [1/s]。"""
    return math.sqrt((144.0 * D_g * eps_mf * u_b) / (math.pi * (d_b**3)))


def calc_kbd_kunii_levenspiel(
    u_mf: float,
    d_b: float,
    D_g: float,
    u_b: float,
    eps_mf: float,
) -> float:
    """Kunii & Levenspiel (1977) 串联阻力模型 Eq. (3.46)–(3.47) [1/s]。

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.46–3.47
    """
    assert u_mf >= 0.0 and u_b >= 0.0 and d_b > 0.0 and D_g > 0.0
    assert 0.0 < eps_mf < 1.0
    k_bc = 4.5 * u_mf / d_b + 5.85 * (D_g**0.5) * (g**0.25) / (d_b**1.25)
    k_cd = 6.78 * math.sqrt(D_g * u_b * eps_mf / (d_b**3))
    return 1.0 / (1.0 / max(k_bc, 1e-30) + 1.0 / max(k_cd, 1e-30))


def calc_kbd_sit_grace_348(
    u_mf: float,
    d_b: float,
    D_g: float,
    eps_mf: float,
    u_b: float,
) -> float:
    """Sit & Grace (1981) 合并气泡模型 Eq. (3.48) [1/s]。

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.48
    """
    assert u_mf >= 0.0 and u_b >= 0.0 and d_b > 0.0 and D_g > 0.0
    assert 0.0 < eps_mf < 1.0
    return 2.0 * u_mf / d_b + _kbd_diffusion_term(D_g, eps_mf, u_b, d_b)


def calc_kbd_preto_349(u_mf: float, d_b: float, chi: float = 1.0) -> float:
    """Preto (1986) 对流主导模型 Eq. (3.49) [1/s]。

    χ=1 为慢速穿透气泡；快速气泡 χ→0.1（下限）。

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.49
    """
    assert u_mf >= 0.0 and d_b > 0.0
    chi_eff = float(max(min(chi, 1.0), 0.1))
    return 4.5 * chi_eff * u_mf / d_b


def calc_kbd(
    u_br: float,
    d_b: float,
    D_g: float,
    eps_mf: float,
    u_b: float,
    *,
    model: KbdModel = "hamel_350",
    u_mf: float | None = None,
    chi: float = 1.0,
) -> float:
    """计算气泡–悬浮相总传质系数 K_bd [1/s]。

    默认 ``model='hamel_350'`` 实现 Hamel 主链 Eq. (3.50)：

        K_bd = (3/2) u_{b,r}/d_b + sqrt(144 D_g ε_mf u_b / (π d_b³))

    可选 ``kunii_levenspiel`` / ``sit_grace_348`` / ``preto_349`` 对照文献模型 (3.46)–(3.49)。

    Source: ``docs/hamel_submodels/03_hydrodynamics_core_chain.md``; Hamel Eq. 3.50
    """
    assert u_br >= 0.0, "u_br 必须非负 [m/s]"
    assert d_b > 0.0, "d_b 必须大于 0 [m]"
    assert D_g > 0.0, "D_g 必须大于 0 [m^2/s]"
    assert 0.0 < eps_mf < 1.0, "eps_mf 必须在 (0,1)"
    assert u_b >= 0.0, "u_b 必须非负 [m/s]"

    if model == "hamel_350":
        conv = 1.5 * u_br / d_b
        return conv + _kbd_diffusion_term(D_g, eps_mf, u_b, d_b)
    if model == "kunii_levenspiel":
        u_mf_eff = 0.5 * u_br if u_mf is None else u_mf
        return calc_kbd_kunii_levenspiel(u_mf_eff, d_b, D_g, u_b, eps_mf)
    if model == "sit_grace_348":
        u_mf_eff = 0.5 * u_br if u_mf is None else u_mf
        return calc_kbd_sit_grace_348(u_mf_eff, d_b, D_g, eps_mf, u_b)
    if model == "preto_349":
        u_mf_eff = 0.5 * u_br if u_mf is None else u_mf
        return calc_kbd_preto_349(u_mf_eff, d_b, chi=chi)
    raise ValueError(f"Unknown K_bd model: {model}")
