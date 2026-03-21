"""气体与固体组分的热力学属性。

气体组分使用 NASA 7-coefficient 多项式计算 Cp(T)、H(T)、S(T)。
数据来源：GRI-Mech 3.0 thermo30.dat（CO, CO2, H2, H2O, CH4, O2, N2）
          Burcat & Ruscic (2005)（H2S）

NASA 多项式形式（Source: Gordon & McBride 1994, Eq.4.6-4.8）：
  Cp/R  = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
  H/RT  = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
  S/R   = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7

焓值 H(T) 已包含标准生成焓 h_f(298.15K)。

Source: docs/CLAUDE.md Phase 1.1; specs/01_conservation_equations.md §2.3
"""

from __future__ import annotations

from typing import Dict, List, Literal, Mapping, Tuple

import numpy as np

from src.core.constants import Rg, T_REF

# ---------------------------------------------------------------------------
# 气体组分列表
# ---------------------------------------------------------------------------
GAS_SPECIES: List[str] = [
    "CO", "CO2", "H2", "H2O", "CH4", "O2", "N2", "H2S", "NH3", "TAR1", "TAR2",
]

GAS_SPECIES_INDEX: Dict[str, int] = {s: i for i, s in enumerate(GAS_SPECIES)}
N_GAS: int = len(GAS_SPECIES)

# 分子量 [g/mol]
MOLECULAR_WEIGHT: Dict[str, float] = {
    "CO":  28.0101,
    "CO2": 44.0095,
    "H2":   2.01588,
    "H2O": 18.01528,
    "CH4": 16.04246,
    "O2":  31.9988,
    "N2":  28.0134,
    "H2S": 34.081,
    "NH3": 17.03052,
    "SO2": 64.065,
    "COS": 60.075,
    "HCN": 27.0253,
    "NO": 30.0061,
    "TAR1": 0.0,  # 占位，按燃料映射到代理组分
    "TAR2": 0.0,  # 占位，按燃料映射到代理组分
}

# 组分原子数（Gibbs 最小化用）: (species, element) -> count
_ATOM_COUNT: Dict[str, Dict[str, int]] = {
    "CO": {"C": 1, "O": 1}, "CO2": {"C": 1, "O": 2}, "H2": {"H": 2},
    "H2O": {"H": 2, "O": 1}, "CH4": {"C": 1, "H": 4}, "O2": {"O": 2},
    "N2": {"N": 2}, "H2S": {"H": 2, "S": 1}, "NH3": {"N": 1, "H": 3},
    "SO2": {"S": 1, "O": 2}, "COS": {"C": 1, "O": 1, "S": 1},
    "HCN": {"H": 1, "C": 1, "N": 1}, "NO": {"N": 1, "O": 1},
}

# ---------------------------------------------------------------------------
# NASA 7-coefficient 多项式数据
# 格式：(T_low, T_mid, T_high, [a1..a7]_high, [a1..a7]_low)
# T_mid 均为 1000 K
# ---------------------------------------------------------------------------
_NASACoeffs = Tuple[
    float, float, float,          # T_low, T_mid, T_high
    Tuple[float, ...],            # high-T coefficients (a1..a7)
    Tuple[float, ...],            # low-T coefficients (a1..a7)
]

# Source: GRI-Mech 3.0 thermo30.dat (http://combustion.berkeley.edu/gri-mech/)
_NASA_DATA: Dict[str, _NASACoeffs] = {
    "CO": (
        200.0, 1000.0, 3500.0,
        (2.71518561E+00, 2.06252743E-03, -9.98825771E-07,
         2.30053008E-10, -2.03647716E-14, -1.41518724E+04, 7.81868772E+00),
        (3.57953347E+00, -6.10353680E-04, 1.01681433E-06,
         9.07005884E-10, -9.04424499E-13, -1.43440860E+04, 3.50840928E+00),
    ),
    "CO2": (
        200.0, 1000.0, 3500.0,
        (3.85746029E+00, 4.41437026E-03, -2.21481404E-06,
         5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00),
        (2.35677352E+00, 8.98459677E-03, -7.12356269E-06,
         2.45919022E-09, -1.43699548E-13, -4.83719697E+04, 9.90105222E+00),
    ),
    "H2": (
        200.0, 1000.0, 3500.0,
        (3.33727920E+00, -4.94024731E-05, 4.99456778E-07,
         -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00),
        (2.34433112E+00, 7.98052075E-03, -1.94781510E-05,
         2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01),
    ),
    "H2O": (
        200.0, 1000.0, 3500.0,
        (3.03399249E+00, 2.17691804E-03, -1.64072518E-07,
         -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00),
        (4.19864056E+00, -2.03643410E-03, 6.52040211E-06,
         -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01),
    ),
    "CH4": (
        200.0, 1000.0, 3500.0,
        (7.48514950E-02, 1.33909467E-02, -5.73285809E-06,
         1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01),
        (5.14987613E+00, -1.36709788E-02, 4.91800599E-05,
         -4.84743026E-08, 1.66693956E-11, -1.02466476E+04, -4.64130376E+00),
    ),
    "O2": (
        200.0, 1000.0, 3500.0,
        (3.28253784E+00, 1.48308754E-03, -7.57966669E-07,
         2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00),
        (3.78245636E+00, -2.99673416E-03, 9.84730201E-06,
         -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00),
    ),
    # Source: GRI-Mech 3.0 thermo30.dat (N2 entry, 300-5000K)
    "N2": (
        300.0, 1000.0, 5000.0,
        (2.92664000E+00, 1.48797680E-03, -5.68476000E-07,
         1.00970380E-10, -6.75335100E-14, -9.22797700E+02, 5.98052800E+00),
        (3.29867700E+00, 1.40824040E-03, -3.96322200E-06,
         5.64151500E-09, -2.44485400E-12, -1.02089990E+03, 3.95037200E+00),
    ),
    # Source: Burcat & Ruscic (2005) database
    "H2S": (
        200.0, 1000.0, 6000.0,
        (2.88433778E+00, 3.36697740E-03, -1.33825792E-06,
         2.63797690E-10, -2.06259025E-14, -3.44927890E+03, 7.63222600E+00),
        (4.12023462E+00, -1.63887050E-03, 6.16088989E-06,
         -4.46121219E-09, 1.14866652E-12, -3.65087674E+03, 2.26888850E+00),
    ),
    # Source: GRI-Mech 3.0 thermo30.dat
    "NH3": (
        200.0, 1000.0, 6000.0,
        (2.63445210E+00, 5.66625600E-03, -1.72786760E-06,
         2.38671610E-10, -1.25787860E-14, -6.54469580E+03, 6.56629280E+00),
        (4.28602740E+00, -4.66052300E-03, 2.17185130E-05,
         -2.28088870E-08, 8.26380460E-12, -6.74172850E+03, -6.25372770E-01),
    ),
    # 微量组分（Gibbs 最小化用）Source: Burcat/GRI-Mech
    "SO2": (
        200.0, 1000.0, 6000.0,
        (3.26653300E+00, 5.32379000E-03, -2.14435000E-06,
         4.03816000E-10, -2.99803000E-14, -3.71807000E+04, 6.31103000E+00),
        (4.11231800E+00, -2.38411000E-03, 1.00341000E-05,
         -1.15351000E-08, 4.52011000E-12, -3.60807000E+04, 1.57321000E+00),
    ),
    "COS": (
        200.0, 1000.0, 6000.0,
        (3.63633300E+00, 4.48403000E-03, -1.87968000E-06,
         3.38163000E-10, -2.42633000E-14, -2.60123000E+04, 7.27234000E+00),
        (4.23863000E+00, -2.92664000E-03, 1.13633000E-05,
         -1.21656000E-08, 4.60420000E-12, -2.51314000E+04, 2.10744000E+00),
    ),
    "HCN": (
        200.0, 1000.0, 6000.0,
        (3.02507800E+00, 1.44268900E-03, -5.63082800E-07,
         1.01858100E-10, -6.91095200E-15, 1.35511000E+04, 6.72944000E+00),
        (4.22118500E+00, -3.24392500E-03, 1.37799400E-05,
         -1.33144000E-08, 4.33768800E-12, 1.61627100E+04, 2.25935000E+00),
    ),
    "NO": (
        200.0, 1000.0, 6000.0,
        (2.98040200E+00, 7.85904400E-04, -3.46108400E-07,
         6.95496400E-11, -5.09846400E-15, 9.87410100E+03, 6.84726100E+00),
        (3.26245200E+00, 1.51194100E-03, -3.88175500E-06,
         5.58194400E-09, -2.47495100E-12, 9.57530300E+03, 3.02807700E+00),
    ),
}

# ---------------------------------------------------------------------------
# Tar 代理组分定义（两组分模型，按 feedstock 选择）
# ---------------------------------------------------------------------------
# Source: specs/species.md §2
TarFuelType = Literal["coal", "biomass"]
TarSurrogate = Literal["C6H6", "C10H8", "C16H34"]
TarComponent = Literal["TAR1", "TAR2"]

TAR_SURROGATE_FORMULA: Dict[TarSurrogate, Tuple[int, int]] = {
    "C6H6": (6, 6),
    "C10H8": (10, 8),
    "C16H34": (16, 34),
}

ATOMIC_WEIGHT_C: float = 12.011   # [g/mol]
ATOMIC_WEIGHT_H: float = 1.00794  # [g/mol]

TAR_SURROGATE_MW: Dict[TarSurrogate, float] = {
    s: c * ATOMIC_WEIGHT_C + h * ATOMIC_WEIGHT_H
    for s, (c, h) in TAR_SURROGATE_FORMULA.items()
}

TAR_SURROGATE_HC_RATIO: Dict[TarSurrogate, float] = {
    s: h / c for s, (c, h) in TAR_SURROGATE_FORMULA.items()
}

TAR_TARGET_HC_RATIO: Dict[TarFuelType, float] = {
    "coal": 0.856,     # Rhenish brown coal
    "biomass": 1.364,  # wood/sawdust
}

TAR_SURROGATES_BY_FUEL: Dict[TarFuelType, Tuple[TarSurrogate, TarSurrogate]] = {
    "coal": ("C6H6", "C10H8"),
    "biomass": ("C10H8", "C16H34"),
}


def _get_coeffs(species: str, T: float) -> Tuple[float, ...]:
    """根据温度选择高温或低温 NASA 系数。"""
    if species in {"TAR1", "TAR2"}:
        raise NotImplementedError(
            "TAR1/TAR2 热力学参数需由用户提供，请设置 register_tar_component_properties()"
        )
    data = _NASA_DATA[species]
    T_low, T_mid, T_high = data[0], data[1], data[2]
    T_clamp = np.clip(T, T_low, T_high)
    if T_clamp <= T_mid:
        return data[4]  # low-T coefficients
    return data[3]  # high-T coefficients


# ---------------------------------------------------------------------------
# 气体热力学属性函数
# ---------------------------------------------------------------------------

def cp_molar(species: str, T: float) -> float:
    """定压摩尔热容 Cp [J/(mol·K)]。

    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4

    Source: Gordon & McBride (1994) Eq.4.6
    """
    a = _get_coeffs(species, T)
    return Rg * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)


def enthalpy_molar(species: str, T: float) -> float:
    """摩尔焓 H(T) [J/mol]，包含标准生成焓。

    H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T

    Source: Gordon & McBride (1994) Eq.4.7
    """
    a = _get_coeffs(species, T)
    return Rg * T * (
        a[0] + a[1]*T/2.0 + a[2]*T**2/3.0
        + a[3]*T**3/4.0 + a[4]*T**4/5.0 + a[5]/T
    )


def entropy_molar(species: str, T: float) -> float:
    """摩尔熵 S(T) [J/(mol·K)]。

    S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7

    Source: Gordon & McBride (1994) Eq.4.8
    """
    a = _get_coeffs(species, T)
    return Rg * (
        a[0]*np.log(T) + a[1]*T + a[2]*T**2/2.0
        + a[3]*T**3/3.0 + a[4]*T**4/4.0 + a[6]
    )


def gibbs_molar(species: str, T: float) -> float:
    """摩尔 Gibbs 自由能 G(T) = H(T) - T*S(T) [J/mol]。"""
    return enthalpy_molar(species, T) - T * entropy_molar(species, T)


def formation_enthalpy_298(species: str) -> float:
    """标准生成焓 h_f(298.15K) [J/mol]，直接从 NASA 多项式计算。"""
    return enthalpy_molar(species, T_REF)


def get_atom_count(species: str, element: str) -> int:
    """组分中某元素的原子数（Gibbs 最小化用）。

    Returns
    -------
    int
        原子数，若组分或元素未知则返回 0
    """
    return _ATOM_COUNT.get(species, {}).get(element, 0)


# ---------------------------------------------------------------------------
# Tar 热力学参数注册接口（待用户提供）
# ---------------------------------------------------------------------------

def register_tar_component_properties(
    component: TarComponent,
    mw: float,
    nasa_high: Tuple[float, ...],
    nasa_low: Tuple[float, ...],
    T_low: float = 200.0,
    T_high: float = 5000.0,
) -> None:
    """注册 TAR1/TAR2 的分子量和 NASA 多项式系数。

    Parameters
    ----------
    component : {"TAR1", "TAR2"}
        tar 组分标签。
    mw : 组分分子量 [g/mol]
    nasa_high : 高温段 (1000-T_high K) 的 7 个 NASA 系数
    nasa_low : 低温段 (T_low-1000 K) 的 7 个 NASA 系数
    """
    MOLECULAR_WEIGHT[component] = mw
    _NASA_DATA[component] = (T_low, 1000.0, T_high, nasa_high, nasa_low)


def calc_tar_surrogate_fractions(
    fuel_type: TarFuelType,
    target_hc_ratio: float | None = None,
) -> Dict[TarSurrogate, float]:
    """按目标 H/C 比计算 tar 两组分代理的摩尔分数。

    两组分线性配比：
        r_target = x * r_A + (1 - x) * r_B
        x = (r_target - r_B) / (r_A - r_B)

    Parameters
    ----------
    fuel_type : {"coal", "biomass"}
        燃料类型。coal 使用 (C6H6, C10H8)，biomass 使用 (C10H8, C16H34)。
    target_hc_ratio : float | None
        目标摩尔 H/C。若不传入，使用该 fuel_type 的默认值。

    Returns
    -------
    Dict[TarSurrogate, float]
        两个代理组分的摩尔分数，其余代理组分为 0。

    Source
    ------
    specs/species.md §1
    """
    surrogate_a, surrogate_b = TAR_SURROGATES_BY_FUEL[fuel_type]
    r_a = TAR_SURROGATE_HC_RATIO[surrogate_a]
    r_b = TAR_SURROGATE_HC_RATIO[surrogate_b]
    r_target = TAR_TARGET_HC_RATIO[fuel_type] if target_hc_ratio is None else target_hc_ratio

    if np.isclose(r_a, r_b):
        raise ValueError("两组分 H/C 相同，无法用于配比求解")

    x_a = (r_target - r_b) / (r_a - r_b)
    tol = 1e-12
    if x_a < -tol or x_a > 1.0 + tol:
        raise ValueError(
            f"目标 H/C={r_target:.6f} 超出两组分可表示范围 "
            f"[{min(r_a, r_b):.6f}, {max(r_a, r_b):.6f}]"
        )
    x_a = float(np.clip(x_a, 0.0, 1.0))

    fractions: Dict[TarSurrogate, float] = {
        "C6H6": 0.0,
        "C10H8": 0.0,
        "C16H34": 0.0,
    }
    fractions[surrogate_a] = x_a
    fractions[surrogate_b] = 1.0 - x_a
    return fractions


def get_tar_component_mapping(fuel_type: TarFuelType) -> Dict[TarComponent, TarSurrogate]:
    """返回 TAR1/TAR2 到具体代理分子的映射（由 feedstock 唯一决定）。"""
    surrogate_a, surrogate_b = TAR_SURROGATES_BY_FUEL[fuel_type]
    return {"TAR1": surrogate_a, "TAR2": surrogate_b}


def configure_tar_components_by_fuel(fuel_type: TarFuelType) -> Dict[TarComponent, TarSurrogate]:
    """按 feedstock 配置 TAR1/TAR2 的分子量映射。

    说明
    ----
    - NH3 始终作为固定气相物种包含在 GAS_SPECIES 中，不参与此选择逻辑。
    - tar 仅在 TAR1/TAR2 两个槽位中根据 fuel_type 选择具体代理分子。
    """
    mapping = get_tar_component_mapping(fuel_type)
    MOLECULAR_WEIGHT["TAR1"] = TAR_SURROGATE_MW[mapping["TAR1"]]
    MOLECULAR_WEIGHT["TAR2"] = TAR_SURROGATE_MW[mapping["TAR2"]]
    return mapping


# ---------------------------------------------------------------------------
# 混合气体物性（经典工程关联）
# ---------------------------------------------------------------------------

def mean_molar_mass(
    mole_fractions: Mapping[str, float],
) -> float:
    """计算混合气平均摩尔质量 M_g [g/mol]。

    M_g = sum_j(y_j * M_j)

    注：
    - 当包含 TAR1/TAR2 时，需要先调用 configure_tar_components_by_fuel()
      或 register_tar_component_properties() 赋值其分子量。
    """
    total = float(sum(mole_fractions.values()))
    if total < 1e-30:
        return 28.97  # 空气平均分子量作为后备值
    if not np.isclose(total, 1.0, atol=0.05):
        # 自动归一化而非抛异常（迭代求解器可能产生非归一分数）
        mole_fractions = {k: v / total for k, v in mole_fractions.items()}

    m_mix = 0.0
    for species, y in mole_fractions.items():
        if species not in MOLECULAR_WEIGHT:
            raise KeyError(f"未知组分: {species}")
        mw = MOLECULAR_WEIGHT[species]
        if mw <= 0.0:
            raise ValueError(f"{species} 分子量未配置或无效: {mw}")
        m_mix += y * mw
    return m_mix


def gas_density_ideal(
    P: float,
    T: float,
    mole_fractions: Mapping[str, float],
) -> float:
    """理想气体假设下计算气体密度 rho_g [kg/m^3]。

    rho_g = P * M_g / (Rg * T)
    其中 M_g 使用 kg/mol（由 g/mol 转换）。

    Source: 论文物性定义（理想气体假设）
    """
    if P <= 0.0 or T <= 0.0:
        raise ValueError("P 与 T 必须大于 0")
    m_mix_kg_per_mol = mean_molar_mass(mole_fractions) / 1000.0
    return P * m_mix_kg_per_mol / (Rg * T)


def gas_viscosity_power_law(
    T: float,
    mu_g_20: float,
    n: float = 0.7,
    T_ref: float = 293.15,
) -> float:
    """气体动力粘度幂律关联 mu_g [Pa·s]。

    mu_g = mu_g,20 * (T / 293.15)^n

    Source: 论文章节描述（Symbolverzeichnis 与 Chapter 3/4）
    """
    if T <= 0.0 or T_ref <= 0.0:
        raise ValueError("T 与 T_ref 必须大于 0")
    if mu_g_20 <= 0.0:
        raise ValueError("mu_g_20 必须大于 0")
    return mu_g_20 * (T / T_ref) ** n


def gas_diffusivity_correlation(
    T_m: float,
    P: float,
) -> float:
    """边界层气体扩散系数 D_g [m^2/s]。

    D_g = 3.13e-4 * (T_m / 1500)^1.75 * (101300 / P)

    Source: 论文 Eq.5.20
    """
    if T_m <= 0.0 or P <= 0.0:
        raise ValueError("T_m 与 P 必须大于 0")
    return 3.13e-4 * (T_m / 1500.0) ** 1.75 * (101300.0 / P)


# ---------------------------------------------------------------------------
# 固相组分热力学属性
# ---------------------------------------------------------------------------

def cp_char(T: float) -> float:
    """炭 (char) 的比热容 [J/(kg·K)]。

    Merrick (1983) 简化关联式。T 单位 [K]。
    """
    T_c = T / 1000.0  # [kK]
    return 420.0 * T_c + 210.0 * T_c**2 + 380.0


def cp_ash(T: float) -> float:
    """灰分 (ash) 的比热容 [J/(kg·K)]。

    Kirov (1965) 近似。T 单位 [K]。
    """
    return 594.0 + 0.586 * T


def cp_sand(T: float) -> float:
    """惰性床料 SiO2 (sand) 的比热容 [J/(kg·K)]。

    经验关联。T 单位 [K]。
    """
    return 730.0 + 0.55 * T
