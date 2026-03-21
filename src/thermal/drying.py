"""干燥子模型（Agarwal + Crank-Nicolson）。

Hamel Kapitel 4（pp. 52–53）核心式与代码对应：
- **Gl. 4.2**：外表面 \(r=R_0\) 热流边界 \(\lambda_s\,\mathrm{d}T/\mathrm{d}r|_{R_0}=\alpha(T_{ws}-T_s)\)；\(\alpha\) 由 **Gl. 4.9** Nu 得。
- **Gl. 4.4**：修正蒸发焓 \(h_v' = h_v + (c_w + c_s/w_{0,tr})(T_e-T_0)\) → `corrected_evaporation_enthalpy`。
- **Gl. 4.6**：干壳温度边界 \(T|_{R_0}=T_s\)，\(T|_{r_e}=T_e\) → CN 径向求解边界结构。
- **Gl. 4.1 / 4.3**：壳层导热与前沿能量平衡（见论文）；数值用 CN 离散。

数值：Crank-Nicolson 径向温度场 + 蒸发前沿推进。
"""

from __future__ import annotations

import numpy as np


# -----------------------------------------------------------------------
# 物性参数（TODO: 待用户提供或从 feedstock 数据确定）
# -----------------------------------------------------------------------
LAMBDA_W: float = 0.6      # [W/(m·K)] 湿颗粒导热系数（褐煤典型值）
RHO_W: float = 1200.0      # [kg/m³]   湿颗粒密度
CP_W: float = 2000.0       # [J/(kg·K)] 湿颗粒比热容
H_EVAP: float = 2.26e6     # [J/kg]    水蒸发潜热（100°C）
T_EVAP: float = 373.15     # [K]       常压蒸发温度
C_WATER: float = 4_180.0   # [J/(kg·K)] 液态水比热


def thermal_diffusivity() -> float:
    """热扩散率 alpha_T [m²/s]。"""
    return LAMBDA_W / (RHO_W * CP_W)


def corrected_evaporation_enthalpy(
    h_v: float = H_EVAP,
    c_w: float = C_WATER,
    c_s: float = CP_W,
    w0_tr: float = 0.169,
    T_e: float = T_EVAP,
    T0: float = 300.0,
) -> float:
    """修正蒸发焓 h_v'（Hamel Gl. 4.4，p. 52–53）。

    h_v' = h_v + (c_w + c_s / w_{0,tr}) * (T_e - T_0)

    w_{0,tr}：干基初始含水率（参数 w0_tr）；c_w、c_s：水与固体比热。
    """
    w_eff = max(w0_tr, 1e-9)
    return h_v + (c_w + c_s / w_eff) * (T_e - T0)


def nusselt_particle(Re: float, Pr: float) -> float:
    """颗粒 Nu 关联式（Eq. 4.9）。"""
    return 2.0 + 1.2 * np.sqrt(max(Re, 0.0)) * (max(Pr, 0.0) ** (1.0 / 3.0))


def convective_htc_from_nusselt(
    d_p: float,
    u_rel: float,
    rho_g: float,
    mu_g: float,
    cp_g: float,
    lambda_g: float,
) -> float:
    """由 Eq. 4.9 计算颗粒表面对流换热系数 alpha。"""
    d_eff = max(d_p, 1e-9)
    Re = rho_g * max(u_rel, 0.0) * d_eff / max(mu_g, 1e-12)
    Pr = cp_g * max(mu_g, 1e-12) / max(lambda_g, 1e-12)
    Nu = nusselt_particle(Re, Pr)
    return Nu * lambda_g / d_eff


# -----------------------------------------------------------------------
# Crank-Nicolson FD 径向热传导求解器
# -----------------------------------------------------------------------

def solve_drying_CN(
    d_p: float,
    T_bed: float,
    T_init: float = 300.0,
    moisture_wt: float = 16.9,
    t_total: float = 10.0,
    Nr: int = 20,
    Nt: int = 200,
    h_conv: float | None = None,
    u_rel: float = 1.0,
    rho_g: float = 0.35,
    mu_g: float = 4.0e-5,
    cp_g: float = 1200.0,
    lambda_g: float = 0.08,
    return_history: bool = False,
) -> dict:
    """Crank-Nicolson 求解球形颗粒径向温度场与干燥进度。

    Parameters
    ----------
    d_p         : 颗粒直径 [m]
    T_bed       : 床层温度 [K]（外边界条件）
    T_init      : 颗粒初始温度 [K]
    moisture_wt : 初始含水率 [wt%]
    t_total     : 总模拟时间 [s]
    Nr          : 径向网格数
    Nt          : 时间步数
    h_conv      : 颗粒表面对流换热系数 [W/(m²·K)]，None 时由 Nu 关联式计算
    u_rel       : 气固相对速度 [m/s]（用于 Nu）
    rho_g       : 气体密度 [kg/m³]（用于 Nu）
    mu_g        : 气体动力黏度 [Pa·s]（用于 Nu）
    cp_g        : 气体比热 [J/(kg·K)]（用于 Nu）
    lambda_g    : 气体导热系数 [W/(m·K)]（用于 Nu）
    return_history : 是否返回径向温度史矩阵（供分层热解耦合）

    Returns
    -------
    dict with keys:
      't'        : 时间数组 [s]
      'X_dry'    : 干燥进度（0=湿，1=全干）随时间变化
      'T_center' : 颗粒中心温度 [K] 随时间
      'T_surface': 颗粒表面温度 [K] 随时间
      'r_evap'   : 蒸发前沿半径 [m] 随时间

    Source: docs/CLAUDE.md Phase 4.1; Crank-Nicolson FD for Stefan problem
    """
    R = d_p / 2.0     # [m] 颗粒半径
    alpha = thermal_diffusivity()
    dr = R / Nr
    dt = t_total / Nt
    if h_conv is None:
        h_conv = convective_htc_from_nusselt(
            d_p=d_p,
            u_rel=u_rel,
            rho_g=rho_g,
            mu_g=mu_g,
            cp_g=cp_g,
            lambda_g=lambda_g,
        )

    r = np.linspace(0, R, Nr + 1)  # 径向节点

    # 初始条件
    T = np.full(Nr + 1, T_init)
    t_arr = np.zeros(Nt + 1)
    X_dry = np.zeros(Nt + 1)
    T_center = np.zeros(Nt + 1)
    T_surface = np.zeros(Nt + 1)
    r_evap = np.full(Nt + 1, R)  # 蒸发前沿从外表面向内推进
    T_history = np.zeros((Nr + 1, Nt + 1)) if return_history else None

    T_center[0] = T_init
    T_surface[0] = T_init
    if return_history:
        T_history[:, 0] = T

    # 球坐标 Crank-Nicolson 系数
    # d(r²·dT/dr)/dr / r² = (1/alpha) dT/dt
    sigma = alpha * dt / (2.0 * dr**2)

    moisture_mass = moisture_wt / 100.0  # 初始质量分数
    h_evap_corr = corrected_evaporation_enthalpy(
        h_v=H_EVAP,
        c_w=C_WATER,
        c_s=CP_W,
        w0_tr=moisture_mass,
        T_e=T_EVAP,
        T0=T_init,
    )
    total_water = moisture_mass * RHO_W * (4.0 / 3.0 * np.pi * R**3)
    evaporated = 0.0

    for n in range(Nt):
        # 构建三对角矩阵（球坐标离散）
        N = Nr + 1
        A = np.zeros(N)
        B = np.zeros(N)
        C = np.zeros(N)
        D = np.zeros(N)

        # 中心对称 BC: dT/dr|_{r=0} = 0
        # 对 i=0 使用 L'Hôpital: ∂²T/∂r² + 2/r * ∂T/∂r -> 3 * ∂²T/∂r²
        B[0] = 1.0 + 6.0 * sigma
        C[0] = -6.0 * sigma
        D[0] = T[0] * (1.0 - 6.0 * sigma) + 6.0 * sigma * T[1]

        # 内部节点
        for i in range(1, Nr):
            ri = r[i]
            rp = ri + 0.5 * dr
            rm = ri - 0.5 * dr
            cp = sigma * rp**2 / ri**2
            cm = sigma * rm**2 / ri**2

            A[i] = -cm
            B[i] = 1.0 + cp + cm
            C[i] = -cp
            D[i] = cm * T[i - 1] + (1.0 - cp - cm) * T[i] + cp * T[i + 1]

        # 表面 BC: -lambda * dT/dr = h*(T_s - T_bed) (Robin BC)
        Bi = h_conv * dr / LAMBDA_W
        B[Nr] = 1.0 + Bi + 2.0 * sigma * (1.0 + 1.0 / Nr)
        A[Nr] = -2.0 * sigma * (1.0 + 1.0 / Nr)
        D[Nr] = (
            T[Nr] * (1.0 - Bi - 2.0 * sigma * (1.0 + 1.0 / Nr))
            + 2.0 * sigma * (1.0 + 1.0 / Nr) * T[Nr - 1]
            + 2.0 * Bi * T_bed
        )

        # Thomas 算法求解三对角系统
        T_new = _thomas_solve(A, B, C, D)

        # 蒸发处理：当 T > T_evap 时，该节点水分蒸发，温度锁定在 T_evap
        # 按 Eq. 4.3 使用修正蒸发焓 h_v' 计算等效蒸发量
        for i in range(Nr, -1, -1):
            if T_new[i] >= T_EVAP and r_evap[n] > r[i]:
                shell_vol = (4.0 / 3.0 * np.pi) * (r_evap[n]**3 - r[i]**3)
                dE = RHO_W * CP_W * shell_vol * (T_new[i] - T_EVAP)
                dm_evap = dE / max(h_evap_corr, 1e-12)
                evaporated += dm_evap
                T_new[i] = T_EVAP
                if total_water > 0:
                    r_evap[n + 1] = r[i]

        T = T_new.copy()
        if r_evap[n + 1] == R:
            r_evap[n + 1] = r_evap[n]

        t_arr[n + 1] = (n + 1) * dt
        X_dry[n + 1] = min(evaporated / max(total_water, 1e-30), 1.0)
        T_center[n + 1] = T[0]
        T_surface[n + 1] = T[Nr]
        if return_history:
            T_history[:, n + 1] = T

    out = {
        "t": t_arr,
        "X_dry": X_dry,
        "T_center": T_center,
        "T_surface": T_surface,
        "r_evap": r_evap,
        "r_nodes": r,
    }
    if return_history:
        out["T_history_rt"] = T_history
    return out


def _thomas_solve(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
) -> np.ndarray:
    """Thomas 算法（三对角矩阵追赶法）。"""
    n = len(b)
    c_ = np.zeros(n)
    d_ = np.zeros(n)
    x = np.zeros(n)

    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]

    for i in range(1, n):
        m = a[i] / (b[i] - a[i] * c_[i - 1])
        c_[i] = c[i] / (b[i] - a[i] * c_[i - 1]) if i < n - 1 else 0.0
        d_[i] = (d[i] - a[i] * d_[i - 1]) / (b[i] - a[i] * c_[i - 1])

    x[n - 1] = d_[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i + 1]

    return x


# -----------------------------------------------------------------------
# 稳态 cell 模型接口
# -----------------------------------------------------------------------

def drying_rate_for_cell(
    d_p: float,
    T_bed: float,
    tau_cell: float,
    moisture_wt: float = 16.9,
    T_init: float = 300.0,
) -> float:
    """计算颗粒在 cell 停留时间 tau_cell 内的干燥率 [kg_H2O/(kg_wet·s)]。

    Parameters
    ----------
    d_p         : 颗粒直径 [m]
    T_bed       : 床层温度 [K]
    tau_cell    : 颗粒在 cell 内的平均停留时间 [s]
    moisture_wt : 初始含水率 [wt%]
    T_init      : 颗粒初始温度 [K]

    Returns
    -------
    drying_rate : [kg_H2O/(kg_wet·s)] 平均干燥速率

    Source: docs/CLAUDE.md Phase 4.1
    """
    result = solve_drying_CN(
        d_p=d_p,
        T_bed=T_bed,
        T_init=T_init,
        moisture_wt=moisture_wt,
        t_total=max(tau_cell, 0.1),
        Nr=15,
        Nt=100,
        h_conv=None,  # 默认用 Eq. 4.9 Nu 关联式估算
    )
    X_final = result["X_dry"][-1]
    return X_final * (moisture_wt / 100.0) / max(tau_cell, 1e-10)
