"""R1–R4 异相反应（炭气化/燃烧）。

R1: 炭燃烧（SCM，k_l）—— C + (1/Φc)O2 -> (2-2/Φc)CO + (2/Φc-1)CO2
R2: 炭 + H2O → CO + H2（水蒸气气化，SCM）
R3: 炭 + 2H2 → CH4（加氢气化，SCM）
R4: Boudouard C + CO2 → 2CO（Weeda，Langmuir-Hinshelwood，SCM）

所有反应速率以 [mol/(m²_ext·s)] 为单位（外表面基准），
在 cell 中乘以颗粒外表面积得到 [mol/s]。

注意：k_b、k_c（Boudouard）的 E 值为负（-E/(Rg*T) 为正，代表吸附热）。

灰层修正、核径收缩、Φc 机制因子：Hamel (1999) §5.1.2, §6.3.2, Eq.5.9–5.10, 6.11, 6.13

R1 采用 **Hamel 流化床有效动力学**（n=0.5 对 O2），k_l 串联式同 Eq.6.11，表面通量 r = k_l·C_O2^0.5。
未反应核径由 **SPM** 与炭转化率关联：d_core = d_p·(1−X)^(1/3)。

Source: specs/04_kinetics.md §R1–R4; docs/missing_parameters_summary.md §A;
        Weeda (1995); Hamel (1999) §5.1, Ch.6, Table 6.3
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from src.core.constants import Rg
from src.kinetics.arrhenius import k_hobbs

# -----------------------------------------------------------------------
# 参数来源：Hamel (1999) §5.1, Ch.6, Table 6.3; Hobbs et al. (1992); Weeda (1995)
# 详见 docs/missing_parameters_summary.md §A
# -----------------------------------------------------------------------

# R1 炭燃烧：Hamel 有效动力学（Shrinking Core / Particle，n=0.5）
# k0 与 Hobbs 同量纲 [m/(s·K)]，经 k_hobbs 得表面速率常数；通量 r = k_l·C_O2^0.5
# coal -> Braunkohle；biomass -> Wood/Holz
R1_HAMEL_K0_COAL: float = 39.0      # [m/(s·K)] Table 6.3 Braunkohle
R1_HAMEL_E_COAL: float = 50_000.0  # [J/mol]
R1_HAMEL_K0_BIOMASS: float = 48.0   # [m/(s·K)] Wood/Holz
R1_HAMEL_E_BIOMASS: float = 53_000.0  # [J/mol]

CharFuelForR1 = Literal["coal", "biomass"]


def r1_hamel_kinetic_constants(fuel: CharFuelForR1) -> tuple[float, float]:
    """R1 Hamel 有效动力学 (k0 [m/(s·K)], E [J/mol])。"""
    if fuel == "coal":
        return R1_HAMEL_K0_COAL, R1_HAMEL_E_COAL
    return R1_HAMEL_K0_BIOMASS, R1_HAMEL_E_BIOMASS


def d_core_from_spm_char_conversion(char_conversion: float, d_p: float) -> float:
    """SPM 未反应炭核直径：d_c = d_p · (1 − X)^(1/3)。

    X 为炭转化率 [0,1]（已反应炭占可反应炭的比例）。

    Source: Hamel (1999) SCM/SPM；球形颗粒、密度均匀假设
    """
    X = float(np.clip(char_conversion, 0.0, 1.0 - 1e-9))
    return float(d_p * (1.0 - X) ** (1.0 / 3.0))

# R2 炭+H2O: Hobbs Eq.5.24, k = 3.42·T·exp(-15600/T)
R2_k0: float = 3.42       # [m/(s·K)]
R2_E: float = 129_700.0   # [J/mol]

# R3 炭+2H2→CH4: 比 R2 慢约 3 个数量级，k = 3.42e-3·T·exp(-15600/T)
R3_k0: float = 3.42e-3    # [m/(s·K)]
R3_E: float = 129_700.0   # [J/mol]

# R4 Boudouard: Weeda (1995) Eq.5.26 氧交换模型
# R = k_a·P_CO2 / (1 + k_b·P_CO2 + k_c·P_CO)
R4_kf_k0: float = 2.72        # k_a [(Pa·s)⁻¹]
R4_kf_E: float = 188_000.0    # [J/mol]
R4_kb_k0: float = 8.31e-7     # k_b [Pa⁻¹] CO2 吸附
R4_kb_E: float = -38_900.0    # [J/mol] 吸附热
R4_kc_k0: float = 2.088e-14   # k_c [Pa⁻¹] CO 抑制
R4_kc_E: float = -24_500.0    # [J/mol] 吸附热

# 灰层扩散参考系数（Hamel Eq.6.13, p.113）
D_A_ref: float = 1.0e-4       # [m²/s] 灰层参考扩散系数


# -----------------------------------------------------------------------
# Φc 机制因子（CO/CO₂ 比，Arthur 1951 + Hamel Eq.5.9–5.10, p.86）
# -----------------------------------------------------------------------

def phi_c(d_core: float, T: float) -> float:
    """机制因子 Φc：定义 R1 化学计量 C + (1/Φc)O2 → (2-2/Φc)CO + (2/Φc-1)CO2。

    d_core [m], T [K]。公式中粒径以 mm 计：d_c_mm = d_core * 1000。

    - d_c <= 0.05 mm：偏 CO₂ 形成（边界层内二次氧化）
    - 0.05 < d_c < 1.0 mm：过渡区
    - d_c >= 1.0 mm：纯 CO 形成（Φc = 1）

    Source: Hamel (1999) Eq.5.9, 5.10; Arthur (1951)
    """
    d_c_mm = d_core * 1000.0  # [mm]
    p_c = 2500.0 * np.exp(-6240.0 / max(T, 300.0))  # Eq.5.10

    if d_c_mm <= 0.05:
        return (2.0 * p_c + 2.0) / (p_c + 2.0)
    if d_c_mm >= 1.0:
        return 1.0
    # 0.05 < d_c < 1.0 mm
    numer = 2.0 * p_c + 2.0 - (p_c / 0.95) * (1000.0 * d_core - 0.05)
    return numer / (p_c + 2.0)


# -----------------------------------------------------------------------
# 灰层有效扩散（Hamel Eq.6.13, p.113）
# -----------------------------------------------------------------------

def D_d_A(D_A: float, Sh: float, d_p: float, d_c: float) -> float:
    """灰层扩散速度 [m/s]。D_d,A = Sh·D_A / (0.5·(d_p + d_c))。

    用于 SCM 中氧通过多孔灰壳的扩散阻力。

    Source: Hamel (1999) Eq.6.13, p.113
    """
    denom = 0.5 * (d_p + max(d_c, 1e-9))
    return Sh * D_A / max(denom, 1e-9)


# -----------------------------------------------------------------------
# 外扩散速率系数（球形颗粒 Sherwood 关联）
# -----------------------------------------------------------------------

def _k_diff_sphere(D_g: float, d_p: float, Sh: float = 2.0) -> float:
    """外扩散传质系数 k_diff [m/s]。

    k_diff = Sh * D_g / d_p

    Sh = 2 (静止球 / Ranz-Marshall 低 Re 极限)

    Source: Hamel (1999) Eq.5.6
    """
    return Sh * D_g / max(d_p, 1e-6)


# -----------------------------------------------------------------------
# R1 炭燃烧（SCM，Hamel Eq.6.11）
# -----------------------------------------------------------------------

def rate_R1(
    T: float,
    C_O2: float,
    d_p: float,
    D_g: float,
    d_core: float | None = None,
    fuel: CharFuelForR1 = "coal",
) -> tuple[float, float]:
    """R1 炭燃烧速率 [mol/(m²·s)] 与化学计量因子 alpha = 1/Φc。

    C + (1/Φc)O2 -> (2-2/Φc)CO + (2/Φc-1)CO2
    返回 (r, alpha)，其中 alpha = 1/Φc 用于 O2 消耗与 CO/CO2 分配。

    **Hamel 有效动力学**：表面反应对 O2 为 **0.5 级**，r = k_l·C_O2^0.5；
    k_ch 为 Hamel Table 6.3（Braunkohle / Wood）与 k_hobbs 形式组合。

    总速率常数 k_l（Eq.6.11）：
    k_l = (1/Φc · [1/k_d,g · d_c²/d_p² + (d_p-d_c)/(2·D_d,A) · d_c/d_p] + 1/k_ch,R1)^(-1)

    Source: Hamel (1999) Eq.6.11, 6.13, Table 6.3; Arthur (1951) Eq.5.9–5.10
    """
    C_O2 = max(C_O2, 0.0)
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9)

    phi = phi_c(d_c, T)
    alpha = 1.0 / max(phi, 1e-6)

    k0_r1, e_r1 = r1_hamel_kinetic_constants(fuel)
    k_ch = k_hobbs(k0_r1, e_r1, T)
    k_d_g = _k_diff_sphere(D_g, d_p)
    D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)

    # Eq.6.11: 串联阻力
    term_gas = (1.0 / max(k_d_g, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    term_ash = ((d_p - d_c) / max(2.0 * D_d_A_val, 1e-30)) * (d_c / max(d_p, 1e-9))
    R_total = (1.0 / phi) * (term_gas + term_ash) + 1.0 / max(k_ch, 1e-30)
    k_l = 1.0 / max(R_total, 1e-30)

    r = k_l * (C_O2 ** 0.5)
    return (r, alpha)


# -----------------------------------------------------------------------
# R2 炭 + H2O → CO + H2（SCM）
# -----------------------------------------------------------------------

def rate_R2(
    T: float,
    C_H2O: float,
    d_p: float,
    D_g: float,
    d_core: float | None = None,
) -> float:
    """R2 水蒸气气化速率 [mol/(m²_ext·s)]，SCM 串联阻力。

    C + H2O -> CO + H2

    有效速率 = 1/(1/k_diff + 1/k_ash_diff + 1/k_chem) * C_H2O

    Source: Hobbs et al. (1992); Hamel (1999) Table 5.2
    """
    C_H2O = max(C_H2O, 0.0)
    k_chem = k_hobbs(R2_k0, R2_E, T)
    k_d = _k_diff_sphere(D_g, d_p)

    # 灰层扩散阻力（SCM，Hamel Eq.6.13 用 D_A_ref）
    if d_core is None:
        d_core = d_p * 0.8  # 默认由颗粒收缩状态计算
    d_c = max(min(d_core, d_p), 1e-9)
    r_p, r_c = d_p / 2.0, d_c / 2.0
    D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)
    # k_ash: 灰壳传质系数，k_ash = 2·D_d,A·d_p / (d_c·(d_p-d_c))
    k_ash = (2.0 * D_d_A_val * d_p / max(d_c * (d_p - d_c), 1e-12)) if d_p > d_c else 1e10

    k_eff_inv = 1.0 / max(k_d, 1e-30) + 1.0 / max(k_ash, 1e-30) + 1.0 / max(k_chem, 1e-30)
    return (1.0 / k_eff_inv) * C_H2O


# -----------------------------------------------------------------------
# R3 炭 + 2H2 → CH4（SCM）
# -----------------------------------------------------------------------

def rate_R3(
    T: float,
    C_H2: float,
    d_p: float,
    D_g: float,
    d_core: float | None = None,
) -> float:
    """R3 加氢气化速率 [mol/(m²_ext·s)]，SCM 串联阻力。

    C + 2H2 -> CH4

    Source: Hobbs et al. (1992); Hamel (1999) Table 5.2
    """
    C_H2 = max(C_H2, 0.0)
    k_chem = k_hobbs(R3_k0, R3_E, T)
    k_d = _k_diff_sphere(D_g, d_p)

    if d_core is None:
        d_core = d_p * 0.8
    d_c = max(min(d_core, d_p), 1e-9)
    r_p, r_c = d_p / 2.0, d_c / 2.0
    D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)
    k_ash = (2.0 * D_d_A_val * d_p / max(d_c * (d_p - d_c), 1e-12)) if d_p > d_c else 1e10

    k_eff_inv = 1.0 / max(k_d, 1e-30) + 1.0 / max(k_ash, 1e-30) + 1.0 / max(k_chem, 1e-30)
    return (1.0 / k_eff_inv) * C_H2


# -----------------------------------------------------------------------
# R4 Boudouard（Weeda 1995，Langmuir-Hinshelwood + SCM）
# -----------------------------------------------------------------------

def rate_R4(
    T: float,
    P_CO2: float,
    P_CO: float,
) -> float:
    """R4 Boudouard 速率 [mol/(m²·s)]，Langmuir-Hinshelwood 机制。

    C + CO2 -> 2CO

    r = k_f * P_CO2 / (1 + k_b * P_CO2 + k_c * P_CO)

    k_f, k_b, k_c 均使用 k_hobbs 形式。
    k_b, k_c 的 E < 0（-E/(Rg*T) > 0），代表吸附热。

    Source: Weeda (1995); Hamel (1999) Table 5.2
    """
    P_CO2 = max(P_CO2, 0.0)
    P_CO = max(P_CO, 0.0)

    kf = k_hobbs(R4_kf_k0, R4_kf_E, T)
    kb = k_hobbs(R4_kb_k0, R4_kb_E, T)
    kc = k_hobbs(R4_kc_k0, R4_kc_E, T)

    denom = 1.0 + kb * P_CO2 + kc * P_CO
    denom = max(denom, 1e-30)
    return kf * P_CO2 / denom


# -----------------------------------------------------------------------
# R4 含外扩散阻力的有效速率
# -----------------------------------------------------------------------

def rate_R4_effective(
    T: float,
    P_CO2: float,
    P_CO: float,
    d_p: float,
    D_g: float,
    d_core: float | None = None,
) -> float:
    """R4 有效速率 [mol/(m²_ext·s)]，含 SCM 扩散阻力。

    简化处理：用 L-H 表面速率 + 外扩散串联。
    严格实现需联立扩散-反应方程，此处为工程近似。
    TODO: 实现严格 SCM + L-H 联立求解

    Source: Weeda (1995); Hamel (1999) Table 5.2
    """
    r_surface = rate_R4(T, P_CO2, P_CO)
    k_d = _k_diff_sphere(D_g, d_p)

    # 将 r_surface 视为 "k_chem * C" 等效，外扩散限制
    C_CO2 = P_CO2 / (Rg * T) if T > 0 else 0.0
    r_diff_limited = k_d * max(C_CO2, 0.0)

    return min(r_surface, r_diff_limited) if r_diff_limited > 0 else r_surface
