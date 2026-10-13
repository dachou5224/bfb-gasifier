"""R1–R4 异相炭反应（床层/自由板固相）。

仓库历史编号保持为：
- R1: 炭燃烧
- R2: 异相水煤气反应
- R3: 加氢气化
- R4: Boudouard 反应

Hamel (1999) 原论文 5.1 的证据链为：
- 5.1.2 / Eq. (5.8)–(5.23) / Table 5.2：炭燃烧（Field + Hobbs，**模拟主链 n=1**）
- 5.1.3 / Eq. (5.24), (5.25)：异相水煤气反应（Hobbs）
- 5.1.4 / Eq. (5.26)–(5.33)：Boudouard 反应（Weeda）
- 5.1.5 / Eq. (5.34), (5.35)：加氢气化（Hobbs）

Ch.6 / Eq. (6.11)–(6.13) / Table 6.3 为**实验评价**半级动力学与灰层扩散，
不是 §5.1.2 模拟主链；仅通过 ``geometry="scm_ash"`` 可选启用灰层串联。

注意：论文中的 5.1.4/5.1.5 顺序与仓库历史函数名不同；本文件保留旧接口，
但在 docstring 中明确标注 thesis equation mapping。
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from src.core.constants import Rg
from src.kinetics.arrhenius import k_hobbs, k_standard

CharFuelForR1 = Literal["coal", "biomass"]
R1Geometry = Literal["spm", "scm_ash"]

# ── Hamel (1999) §5.1 动力学常数 ──
# Table 5.2：Rheinische Braunkohle -> Lignite A
R1_K0_LIGNITE_A: float = 1.2
R1_E_Rg_LIGNITE_A: float = 10_300.0
# Table 5.2 还给出 "alle Kohlen" 汇总行；本仓库在非 coal 分支上用它作保守回退。
R1_K0_ALL_COALS: float = 2.3
R1_E_Rg_ALL_COALS: float = 11_100.0

R2_k0: float = 3.42
R2_E_Rg: float = 15_600.0

# 论文 5.1.5 把加氢气化写作 k_ch,R4；仓库历史接口保留为 rate_R3。
R3_k0: float = 3.42e-3
R3_E_Rg: float = 15_600.0

R4_kf_k0: float = 2.72  # [(Pa·s)⁻¹] Hamel Eq. (5.27) k_{a,0}
R4_kf_E: float = 188_000.0  # [J/mol] Hamel Eq. (5.30)
R4_kb_k0: float = 8.31e-7  # [Pa⁻¹] Hamel Eq. (5.28)
R4_kb_E: float = -38_900.0  # [J/mol] 吸附焓，Hamel Eq. (5.31)
R4_kc_k0: float = 2.088e-14  # [Pa⁻¹] Hamel Eq. (5.29)
R4_kc_E: float = -24_500.0  # [J/mol] 吸附焓，Hamel Eq. (5.32)
# Eq. (5.33) 碳摩尔质量 M_C；与 elemental_ledger.ATOMIC_MASS_KG_PER_MOL["C"] 一致
R4_M_C: float = 12.011e-3  # [kg/mol]

# Ch.6 实验评价灰层扩散参考（Eq.6.13）；仅 scm_ash 路径使用。
D_A_ref: float = 1.0e-4  # [m²/s]


def r1_hamel_kinetic_constants(fuel: CharFuelForR1) -> tuple[float, float]:
    """Return Hobbs R1 constants in SI-ready form.

    Ref: Hamel (1999) Table 5.2, Eq. (5.23)

    `coal` follows Hamel's Rheinische Braunkohle -> Lignite A row.
    Hamel 5.1.2 does not provide a biomass-specific row; the non-coal branch
    therefore falls back to the table's generic "alle Kohlen" values.
    """
    if fuel == "coal":
        return R1_K0_LIGNITE_A, R1_E_Rg_LIGNITE_A * Rg
    return R1_K0_ALL_COALS, R1_E_Rg_ALL_COALS * Rg


def d_core_from_spm_char_conversion(char_conversion: float, d_p: float) -> float:
    """Shrinking-particle / core diameter proxy from conversion.

    Ref: Hamel (1999) §5.1.1.3 SPM；Ch.6 SCM 核径。

    对 §5.1.2 模拟主链（``geometry="spm"``），粒径类 ``d_p`` 已是当前外径，
    ``rate_R1`` 默认不用本函数收缩核径。保留供自由板/库存与 scm_ash 路径。
    """
    x_raw = float(char_conversion)
    if x_raw < 0.0:
        X = 0.0
    elif x_raw > 1.0 - 1e-9:
        X = 1.0 - 1e-9
    else:
        X = x_raw
    return float(d_p * (1.0 - X) ** (1.0 / 3.0))


def phi_c(d_core: float, T: float) -> float:
    """Combustion mechanism factor ``phi_c``.

    Ref: Hamel (1999) Eq. (5.9), (5.10)
    """
    d_c_mm = d_core * 1000.0
    p_c = k_standard(2500.0, 6240.0 * Rg, max(T, 300.0))
    if d_c_mm <= 0.05:
        return (2.0 * p_c + 2.0) / (p_c + 2.0)
    if d_c_mm >= 1.0:
        return 1.0
    numer = 2.0 * p_c + 2.0 - (p_c / 0.95) * (1000.0 * d_core - 0.05)
    return numer / (p_c + 2.0)


def D_d_A(D_A: float, Sh: float, d_p: float, d_c: float) -> float:
    """灰层扩散速度系数 [m/s]。

    Ref: Hamel (1999) Eq. (6.13) — Ch.6 实验评价路径
    """
    denom = 0.5 * (d_p + max(d_c, 1e-9))
    return Sh * D_A / max(denom, 1e-9)


def sherwood_number_eq517(
    d_p: float,
    D_g: float,
    u_slip: float | None = None,
    nu_g: float | None = None,
) -> float:
    """Sherwood number for film diffusion.

    Ref: Hamel (1999) Eq. (5.17)–(5.19)

    ``Sh = 2 + 0.6 Re^{0.5} Sc^{1/3}``，``Re = u_slip d_p / ν``，``Sc = ν / D_g``。
    缺滑移/粘度时取 Re→0 极限 ``Sh = 2``（仍是 Eq.5.17 的极限，非另设常数模型）。
    """
    d_p = max(float(d_p), 1e-9)
    D_g = max(float(D_g), 1e-30)
    if u_slip is None or nu_g is None:
        return 2.0
    u = max(float(u_slip), 0.0)
    nu = max(float(nu_g), 1e-30)
    if u <= 0.0:
        return 2.0
    Re = u * d_p / nu
    Sc = nu / D_g
    return 2.0 + 0.6 * (Re**0.5) * (Sc ** (1.0 / 3.0))


def _k_diff_sphere(D_g: float, d_p: float, Sh: float = 2.0) -> float:
    """外扩散传质系数 k_dg [m/s]。

    ``k_dg = Sh * D_g / d_p``

    Ref: Hamel (1999) Eq. (5.15)–(5.17) 语境下的 Field 气膜系数
    """
    return float(Sh) * float(D_g) / max(float(d_p), 1e-6)


def rate_R1(
    T: float,
    C_O2: float,
    d_p: float,
    D_g: float,
    d_core: float | None = None,
    fuel: CharFuelForR1 = "coal",
    *,
    u_slip: float | None = None,
    nu_g: float | None = None,
    geometry: R1Geometry = "spm",
) -> tuple[float, float]:
    """R1 char combustion surface rate [mol/(m²·s)] and ``alpha = 1/φ_c``.

    Ref: Hamel (1999) Eq. (5.8)–(5.23), Table 5.2；§5.1.1.3 SPM

    模拟主链（默认 ``geometry="spm"``）：
    - 反应级数 **n = 1**（Eq.5.13–5.14）
    - 浓度驱动力 ``C_O2`` [mol/m³]
    - 气膜 ``k_dg = Sh D_g / d_p``，``Sh`` 来自 Eq.5.17
    - 消去表面浓度：``r = C_O2 / (1/k_ch + 1/(φ_c k_dg))``
    - 无灰层（SPM）；``d_p`` 为当前粒径类外径；``d_core`` 忽略

    可选 ``geometry="scm_ash"``：叠加 Eq.6.13 灰层串联（Ch.6 评价路径），仍用 n=1。
    """
    C_O2 = max(float(C_O2), 0.0)
    d_p = max(float(d_p), 1e-9)
    D_g = max(float(D_g), 1e-30)
    T = max(float(T), 300.0)

    if geometry == "scm_ash" and d_core is not None:
        d_c = max(min(float(d_core), d_p), 1e-9)
    else:
        # SPM：粒径类外径即反应外径（§5.1.1.3）
        d_c = d_p

    phi = float(phi_c(d_c, T))
    alpha = 1.0 / max(phi, 1e-30)
    k0_r1, e_r1 = r1_hamel_kinetic_constants(fuel)
    k_ch = float(k_hobbs(k0_r1, e_r1, T))  # [m/s]

    Sh = sherwood_number_eq517(d_p, D_g, u_slip=u_slip, nu_g=nu_g)
    k_dg = _k_diff_sphere(D_g, d_p, Sh=Sh)  # [m/s]

    if C_O2 <= 1e-30 or k_ch <= 0.0:
        return 0.0, alpha

    if geometry == "scm_ash" and d_c < d_p - 1e-15:
        # Ch.6 对照：外表面基准串联（气膜 + 灰层 + 核表面反应面积比）
        # Ref: Hamel (1999) Eq. (6.10)–(6.13) 结构；级数仍用 §5.1.2 的 n=1
        D_d_A_val = D_d_A(D_A_ref, Sh, d_p, d_c)
        res_gas = 1.0 / max(k_dg, 1e-30)
        res_ash = ((d_p - d_c) / max(2.0 * D_d_A_val, 1e-30)) * (d_c / max(d_p, 1e-9))
        area_ratio = (d_p / max(d_c, 1e-9)) ** 2
        k_eff = 1.0 / (
            area_ratio / max(k_ch, 1e-30) + (1.0 / phi) * (res_gas + res_ash)
        )
    else:
        # Eq.5.12–5.14：Rc/Ap = C∞ / (1/k_ch + 1/(φ k_dg))
        k_eff = 1.0 / (1.0 / max(k_ch, 1e-30) + 1.0 / max(phi * k_dg, 1e-30))

    r = k_eff * C_O2
    return (float(r), float(alpha))


def rate_R2(T: float, C_H2O: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    """R2 heterogeneous water-gas reaction (Hamel surface rate).

    Ref: Hamel (1999) Eq. (5.24), (5.25)

    Hamel §5.1.3 仅给出 ``R_C = A_C * k_ch,R2 * C_H2O``，不含 R1/Ch.6 的气膜或
    灰层串联阻力。``d_p``、``D_g``、``d_core`` 保留以兼容既有调用签名，不参与计算。
    """
    _ = (d_p, D_g, d_core)
    C_H2O = max(C_H2O, 0.0)
    k_ch = k_hobbs(R2_k0, R2_E_Rg * Rg, T)
    return k_ch * C_H2O


def rate_R3(T: float, C_H2: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    """Repository R3 = thesis hydrogasification (Hamel R4 branch).

    Ref: Hamel (1999) Eq. (5.34), (5.35)

    Thesis §5.1.5 writes ``R_C ~ A_C * k_ch,R4 * C_H2`` with Hobbs ``k_hobbs`` form.
    ``d_p``、``D_g``、``d_core`` 保留以兼容调用签名，不参与计算。
    """
    _ = (d_p, D_g, d_core)
    C_H2 = max(C_H2, 0.0)
    k_ch = k_hobbs(R3_k0, R3_E_Rg * Rg, T)
    return k_ch * C_H2


def rate_R4(T: float, P_CO2: float, P_CO: float) -> float:
    """Repository R4 = thesis Boudouard branch (Weeda L-H).

    Source: docs/hamel_submodels（§5.1.4 待补）；Hamel (1999) Eq. (5.26)–(5.32)

    Returns
    -------
    R*_C,0 : float
        质量比速率 [1/s]，相对 Weeda 试样初始质量；分压 P_CO2/P_CO 单位 Pa。
        换算为摩尔速率须用 ``extent_R4`` / Eq. (5.33)：R = R* · m / M_C。
    """
    P_CO2 = max(float(P_CO2), 0.0)  # [Pa]
    P_CO = max(float(P_CO), 0.0)  # [Pa]
    kf = k_standard(R4_kf_k0, R4_kf_E, T)  # [(Pa·s)⁻¹]
    kb = k_standard(R4_kb_k0, R4_kb_E, T)  # [Pa⁻¹]
    kc = k_standard(R4_kc_k0, R4_kc_E, T)  # [Pa⁻¹]
    return kf * P_CO2 / (1.0 + kb * P_CO2 + kc * P_CO)


def extent_R4(T: float, P_CO2: float, P_CO: float, m_char_kg: float) -> float:
    """Boudouard 摩尔速率 [mol C/s]。

    Source: Hamel (1999) Eq. (5.33)
    ``R_C,0 = R*_C,0 · m_{j,k} / M_C``；粒径影响忽略，故直接用 ``rate_R4``。
    """
    m = max(float(m_char_kg), 0.0)  # [kg]
    return rate_R4(T, P_CO2, P_CO) * m / R4_M_C


def rate_R4_effective(
    T: float,
    P_CO2: float,
    P_CO: float,
    d_p: float,
    D_g: float,
    d_core: float | None = None,
) -> float:
    """Boudouard 质量比速率别名 [1/s]（Eq.5.33 忽略粒径/扩散）。

    Ref: Hamel (1999) Eq. (5.26)–(5.33)

    返回与 ``rate_R4`` 相同的 R* [1/s]。摩尔速率请用 ``extent_R4``，
    **不要**再乘外表面积。``d_p``、``D_g``、``d_core`` 仅保留调用签名兼容。
    """
    _ = (d_p, D_g, d_core)
    return rate_R4(T, P_CO2, P_CO)
