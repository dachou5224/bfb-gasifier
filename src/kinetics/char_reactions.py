"""R1–R4 异相反应（炭气化/燃烧）。

R1: 炭燃烧（SCM，k_l）—— C + (1/Φc)O2 -> (2-2/Φc)CO + (2/Φc-1)CO2
R2: 炭 + H2O → CO + H2（水蒸气气化，SCM）
R3: 炭 + 2H2 → CH4（加氢气化，SCM）
R4: Boudouard C + CO2 → 2CO（Weeda，Langmuir-Hinshelwood，SCM）

Source: Hamel (1999) §5.1, Ch.6, Table 6.3; Hobbs et al. (1992); Weeda (1995)
"""

from __future__ import annotations
from typing import Literal
import numpy as np
from src.core.constants import Rg, g, P0
from src.kinetics.arrhenius import k_hobbs, k_standard

CharFuelForR1 = Literal["coal", "biomass"]

# ── 动力学常数重标定 (Calibration) ──
R1_HAMEL_K0_COAL: float = 0.5       
R1_HAMEL_E_COAL: float = 50_000.0  
R1_HAMEL_K0_BIOMASS: float = 0.6    
R1_HAMEL_E_BIOMASS: float = 53_000.0

R2_k0: float = 0.15       
R2_E: float = 129_700.0   

R3_k0: float = 3.42e-5    
R3_E: float = 129_700.0   

R4_kf_k0: float = 0.5        
R4_kf_E: float = 188_000.0    
R4_kb_k0: float = 8.31e-7     
R4_kb_E: float = -38_900.0    
R4_kc_k0: float = 2.088e-14   
R4_kc_E: float = -24_500.0    

D_A_ref: float = 1.0e-4       

def r1_hamel_kinetic_constants(fuel: CharFuelForR1) -> tuple[float, float]:
    if fuel == "coal": return R1_HAMEL_K0_COAL, R1_HAMEL_E_COAL
    return R1_HAMEL_K0_BIOMASS, R1_HAMEL_E_BIOMASS

def d_core_from_spm_char_conversion(char_conversion: float, d_p: float) -> float:
    X = float(np.clip(char_conversion, 0.0, 1.0 - 1e-9))
    return float(d_p * (1.0 - X) ** (1.0 / 3.0))

def phi_c(d_core: float, T: float) -> float:
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
    C_O2 = max(C_O2, 0.0); P_O2 = C_O2 * Rg * T
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9)
    # 强制工程常用值 phi = 1.8 (alpha = 0.55)，偏向 CO 生成
    phi = 1.8; alpha = 0.55
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
    C_H2O = max(C_H2O, 0.0); k_chem = k_hobbs(R2_k0, R2_E, T); k_d = _k_diff_sphere(D_g, d_p)
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9); D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)
    res_gas = (1.0 / max(k_d, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_val, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_l = 1.0 / (res_gas + res_ash + 1.0/max(k_chem, 1e-30))
    return k_l * C_H2O

def rate_R3(T: float, C_H2: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    C_H2 = max(C_H2, 0.0); k_chem = k_hobbs(R3_k0, R3_E, T); k_d = _k_diff_sphere(D_g, d_p)
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9); D_d_A_val = D_d_A(D_A_ref, 2.0, d_p, d_c)
    res_gas = (1.0 / max(k_d, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_val, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_l = 1.0 / (res_gas + res_ash + 1.0/max(k_chem, 1e-30))
    return k_l * C_H2

def rate_R4(T: float, P_CO2: float, P_CO: float) -> float:
    kf = k_standard(R4_kf_k0, R4_kf_E, T); kb = k_standard(R4_kb_k0, R4_kb_E, T); kc = k_standard(R4_kc_k0, R4_kc_E, T)
    return kf * P_CO2 / (1.0 + kb * P_CO2 + kc * P_CO)

def rate_R4_effective(T: float, P_CO2: float, P_CO: float, d_p: float, D_g: float, d_core: float | None = None) -> float:
    r_kin = rate_R4(T, P_CO2, P_CO); RT = Rg * T; k_d_p = _k_diff_sphere(D_g, d_p) / RT
    d_c = d_p if d_core is None else max(min(d_core, d_p), 1e-9); D_d_A_p = D_d_A(D_A_ref, 2.0, d_p, d_c) / RT
    res_gas = (1.0 / max(k_d_p, 1e-30)) * (d_c**2 / max(d_p**2, 1e-18))
    res_ash = ((d_p - d_c) / max(2.0 * D_d_A_p, 1e-30)) * (d_c / max(d_p, 1e-9))
    k_ext_p = 1.0 / (res_gas + res_ash)
    return 1.0 / (1.0/max(r_kin, 1e-30) + 1.0/(k_ext_p * max(P_CO2, 1e-6)))
