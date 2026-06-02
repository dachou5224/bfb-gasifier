"""R1–R4 异相炭反应（床层/自由板固相）。

仓库历史编号保持为：
- R1: 炭燃烧
- R2: 异相水煤气反应
- R3: 加氢气化
- R4: Boudouard 反应

Hamel (1999) 原论文 5.1 的证据链为：
- 5.1.2 / Eq. (5.23) / Table 5.2：炭燃烧动力学（Field + Hobbs）
- 5.1.3 / Eq. (5.24), (5.25)：异相水煤气反应（Hobbs）
- 5.1.4 / Eq. (5.26)–(5.33)：Boudouard 反应（Weeda）
- 5.1.5 / Eq. (5.34), (5.35)：加氢气化（Hobbs）

注意：论文中的 5.1.4/5.1.5 顺序与仓库历史函数名不同；本文件保留旧接口，
但在 docstring 中明确标注 thesis equation mapping。
"""

from __future__ import annotations
from typing import Literal
import numpy as np
from src.core.constants import Rg
from src.kinetics.arrhenius import k_hobbs, k_jensen_r7, k_standard

CharFuelForR1 = Literal["coal", "biomass"]

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

R4_kf_k0: float = 2.72
R4_kf_E: float = 188_000.0
R4_kb_k0: float = 8.31e-7     
R4_kb_E: float = -38_900.0    
R4_kc_k0: float = 2.088e-14   
R4_kc_E: float = -24_500.0    

D_A_ref: float = 1.0e-4       

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
    """Shrinking-particle core diameter.

    Ref: Hamel (1999) uses the shrinking-particle formulation in §5.1.1.3 and p.84.
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
    if d_c_mm <= 0.05: return (2.0 * p_c + 2.0) / (p_c + 2.0)
    if d_c_mm >= 1.0: return 1.0
    numer = 2.0 * p_c + 2.0 - (p_c / 0.95) * (1000.0 * d_core - 0.05)
    return numer / (p_c + 2.0)

def D_d_A(D_A: float, Sh: float, d_p: float, d_c: float) -> float:
    denom = 0.5 * (d_p + max(d_c, 1e-9))
    return Sh * D_A / max(denom, 1e-9)

def _k_diff_sphere(D_g: float, d_p: float, Sh: float = 2.0) -> float:
    return Sh * D_g / max(d_p, 1e-6)

def rate_R1(T: float, C_O2: float, d_p: float, D_g: float, d_core: float | None = None, fuel: CharFuelForR1 = "coal") -> tuple[float, float]:
    """R1 char combustion rate with Field/Hobbs closure.

    Ref: Hamel (1999) Eq. (5.8)–(5.23), Table 5.2
    """
    C_O2 = max(C_O2, 0.0); P_O2 = C_O2 * Rg * T
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9)
    phi = phi_c(d_c, T); alpha = 1.0 / phi
    k0_r1, e_r1 = r1_hamel_kinetic_constants(fuel)
    k_ch = k_hobbs(k0_r1, e_r1, T); RT = Rg * T
    k_d_g_p = _k_diff_sphere(D_g, d_p) / RT; D_d_A_p = D_d_A(D_A_ref, 2.0, d_p, d_c) / RT
    res_gas = (1.0 / max(k_d_g_p, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_p, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_ext_p = 1.0 / max((1.0 / phi) * (res_gas + res_ash), 1e-30)
    a_quad = 1.0 / max(k_ch**2, 1e-30); b_quad = 1.0 / k_ext_p; c_quad = -P_O2
    if P_O2 <= 1e-10: r = 0.0
    else: r = (-b_quad + np.sqrt(max(b_quad**2 - 4.0 * a_quad * c_quad, 0.0))) / (2.0 * a_quad)
    return (r, alpha)

def rate_R2(T: float, C_H2O: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    """R2 heterogeneous water-gas reaction.

    Ref: Hamel (1999) Eq. (5.24), (5.25)
    """
    C_H2O = max(C_H2O, 0.0); k_chem = k_hobbs(R2_k0, R2_E_Rg * Rg, T); k_d = _k_diff_sphere(D_g, d_p)
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9); D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)
    res_gas = (1.0 / max(k_d, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_val, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_l = 1.0 / (res_gas + res_ash + 1.0/max(k_chem, 1e-30))
    return k_l * C_H2O

def rate_R3(T: float, C_H2: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    """Repository R3 = thesis hydrogasification branch.

    Ref: Hamel (1999) Eq. (5.34), (5.35)

    Thesis §5.1.5 writes the chemical constant as ``k_ch,R4``. The repository keeps the
    historical public function name `rate_R3`, but the rate law itself follows Eq. (5.35).
    """
    C_H2 = max(C_H2, 0.0); k_chem = k_jensen_r7(R3_k0, R3_E_Rg, T); k_d = _k_diff_sphere(D_g, d_p)
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9); D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)
    res_gas = (1.0 / max(k_d, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_val, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_l = 1.0 / (res_gas + res_ash + 1.0/max(k_chem, 1e-30))
    return k_l * C_H2

def rate_R4(T: float, P_CO2: float, P_CO: float) -> float:
    """Repository R4 = thesis Boudouard branch.

    Ref: Hamel (1999) Eq. (5.26)–(5.32)
    """
    kf = k_standard(R4_kf_k0, R4_kf_E, T); kb = k_standard(R4_kb_k0, R4_kb_E, T); kc = k_standard(R4_kc_k0, R4_kc_E, T)
    return kf * P_CO2 / (1.0 + kb * P_CO2 + kc * P_CO)

def rate_R4_effective(T: float, P_CO2: float, P_CO: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    """Boudouard rate with external/ash diffusion correction.

    Ref: Hamel (1999) Eq. (5.26)–(5.33)
    """
    r_kin = rate_R4(T, P_CO2, P_CO); RT = Rg * T; k_d_p = _k_diff_sphere(D_g, d_p) / RT
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9); D_d_A_p = D_d_A(D_A_ref, 2.0, d_p, d_c) / RT
    res_gas = (1.0 / max(k_d_p, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_p, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_ext_p = 1.0 / (res_gas + res_ash)
    return 1.0 / (1.0/max(r_kin, 1e-30) + 1.0/(k_ext_p * max(P_CO2, 1e-6)))
