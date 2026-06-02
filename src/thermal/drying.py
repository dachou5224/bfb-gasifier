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

from src.core.constants import MMHG_TO_PA, P0


# -----------------------------------------------------------------------
# 物性参数（Hamel Tabelle 4.1，Dersch 1994，褐煤干燥颗粒）
# -----------------------------------------------------------------------
RHO_W: float = 1250.0      # [kg/m³]   论文表 rho_T
CP_W: float = 1256.0       # [J/(kg·K)] 论文表 c_T
LAMBDA_W: float = 0.157    # [W/(m·K)] 由 a=0.1e-6 m²/s = lambda/(rho*c) 反算
H_EVAP: float = 2.26e6     # [J/kg]    水蒸发潜热（100°C）
T_EVAP: float = 373.15     # [K]       常压蒸发温度
C_WATER: float = 4_180.0   # [J/(kg·K)] 液态水比热


def saturation_temperature_water(pressure_pa: float) -> float:
    """水饱和温度 [K]（由压力决定）。

    采用 Antoine 经验式的两段参数，满足 Hamel 口径中“Te 取该压力下沸点，
    不再额外经验修正”的实现需求。
    """
    p = float(np.clip(pressure_pa, 611.0, 22.064e6))
    # Antoine 公式使用 mmHg、T[°C]
    p_mmhg = p / MMHG_TO_PA
    # 低温段/高温段切换（100°C 左右）
    if p <= P0:
        a, b, c = 8.07131, 1730.63, 233.426
    else:
        a, b, c = 8.14019, 1810.94, 244.485
    t_c = b / max(a - np.log10(max(p_mmhg, 1e-12)), 1e-12) - c
    return float(np.clip(t_c + 273.15, 273.15, 647.096))


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
    pressure_pa: float = P0,
    return_history: bool = False,
) -> dict:
    """Crank-Nicolson 球形颗粒径向干燥（Hamel Eq. 4.1-4.4, 4.6, Stefan 前沿）。

    **两阶段算法**（忠实于 Hamel Kapitel 4, Agarwal 1986 模型）：

    Phase 1 加热阶段：对整球求 CN，直到外表面温度达到 T_evap。
    Phase 2 干燥阶段（Stefan 前沿，Hamel Eq. 4.1 + 4.6 + 4.3）：
        * 仅在干壳 r_e ≤ r ≤ R 上求解热传导（Eq. 4.1）；
        * 内侧 Dirichlet BC：T|_{r_e} = T_e（Eq. 4.6）；
        * 外侧 Robin BC：λ dT/dr|_{R} = α(T_a - T_s)（Eq. 4.2）；
        * Stefan 条件驱动前沿推进（Eq. 4.3）：
              q_front = λ_s · dT/dr|_{r_e+}
              dr_e/dt = q_front / (ρ_s · w · h_v')
              dm/dt   = q_front · 4π r_e² / h_v'

    Parameters
    ----------
    d_p         : 颗粒直径 [m]
    T_bed       : 床层温度 [K]（外边界条件，即 T_a）
    T_init      : 颗粒初始温度 [K]
    moisture_wt : 初始含水率 [wt%]（湿基）
    t_total     : 总模拟时间 [s]
    Nr          : 径向网格数
    Nt          : 时间步数
    h_conv      : 颗粒表面对流换热系数 [W/(m²·K)]，None 时由 Eq. 4.9 计算
    u_rel       : 气固相对速度 [m/s]（用于 Eq. 4.9 Nu 关联式）
    rho_g       : 参考气体密度 [kg/m³]，计算 h_conv 时用压力修正
    mu_g        : 气体动力黏度 [Pa·s]
    cp_g        : 气体比热 [J/(kg·K)]
    lambda_g    : 气体导热系数 [W/(m·K)]
    pressure_pa : 操作压力 [Pa]（决定 T_evap）
    return_history : 返回径向温度矩阵（供热解 DAEM 耦合）

    Returns
    -------
    dict with keys:
      't'           : 时间数组 [s]
      'X_dry'       : 干燥进度（0=湿，1=全干）随时间
      'T_center'    : 颗粒中心温度 [K]
      'T_surface'   : 颗粒表面温度 [K]
      'r_evap'      : 蒸发前沿半径 [m]
      'T_history_rt': (Nr+1, Nt+1) 温度场矩阵（return_history=True 时）

    Sources: Hamel (1999) Eq. 4.1-4.4, 4.6, 4.9; Agarwal et al. (1986)
    """
    R = d_p / 2.0
    alpha = thermal_diffusivity()
    dr = R / max(Nr, 1)
    dt = t_total / max(Nt, 1)

    # h_conv：用压力修正 rho_g（理想气体，rho ∝ P）以还原高压效果
    if h_conv is None:
        rho_g_eff = rho_g * (pressure_pa / P0)
        h_conv = convective_htc_from_nusselt(
            d_p=d_p,
            u_rel=u_rel,
            rho_g=rho_g_eff,
            mu_g=mu_g,
            cp_g=cp_g,
            lambda_g=lambda_g,
        )

    r = np.linspace(0.0, R, Nr + 1)
    sigma = alpha * dt / (2.0 * dr**2)
    Bi_surf = h_conv * dr / LAMBDA_W  # 表面数值 Biot（Robin BC 系数）

    moisture_mass = moisture_wt / 100.0
    w0_tr = moisture_mass / max(1.0 - moisture_mass, 1e-12)
    t_evap_val = saturation_temperature_water(pressure_pa)
    h_evap_corr = corrected_evaporation_enthalpy(
        h_v=H_EVAP,
        c_w=C_WATER,
        c_s=CP_W,
        w0_tr=w0_tr,
        T_e=t_evap_val,
        T0=T_init,
    )
    total_water = moisture_mass * RHO_W * (4.0 / 3.0 * np.pi * R**3)
    evaporated = 0.0

    T = np.full(Nr + 1, T_init)
    t_arr = np.zeros(Nt + 1)
    X_dry_arr = np.zeros(Nt + 1)
    T_center_arr = np.zeros(Nt + 1)
    T_surface_arr = np.zeros(Nt + 1)
    r_evap_arr = np.full(Nt + 1, R)
    T_history = np.zeros((Nr + 1, Nt + 1)) if return_history else None

    T_center_arr[0] = T_init
    T_surface_arr[0] = T_init
    if return_history:
        T_history[:, 0] = T.copy()

    r_front = R
    drying_phase = T_init >= t_evap_val  # 初始已热则直接进入干燥阶段

    for n in range(Nt):
        N = Nr + 1
        Am = np.zeros(N)
        Bm = np.zeros(N)
        Cm = np.zeros(N)
        Dm = np.zeros(N)

        if not drying_phase:
            # ── Phase 1: 加热阶段，对整球求 CN ─────────────────────────
            Bm[0] = 1.0 + 6.0 * sigma
            Cm[0] = -6.0 * sigma
            Dm[0] = T[0] * (1.0 - 6.0 * sigma) + 6.0 * sigma * T[1]
            for i in range(1, Nr):
                ri = r[i]; rp = ri + 0.5 * dr; rm = ri - 0.5 * dr
                cp_i = sigma * rp**2 / ri**2
                cm_i = sigma * rm**2 / ri**2
                Am[i] = -cm_i
                Bm[i] = 1.0 + cp_i + cm_i
                Cm[i] = -cp_i
                Dm[i] = cm_i * T[i-1] + (1.0 - cp_i - cm_i) * T[i] + cp_i * T[i+1]
        else:
            # ── Phase 2: 干燥阶段 ────────────────────────────────────────
            fully_dry = evaporated >= total_water
            if fully_dry:
                # 全干：回到全球 CN（无湿核 Dirichlet BC），与 Phase 1 相同
                Bm[0] = 1.0 + 6.0 * sigma
                Cm[0] = -6.0 * sigma
                Dm[0] = T[0] * (1.0 - 6.0 * sigma) + 6.0 * sigma * T[1]
                for i in range(1, Nr):
                    ri = r[i]; rp = ri + 0.5 * dr; rm = ri - 0.5 * dr
                    cp_i = sigma * rp**2 / ri**2
                    cm_i = sigma * rm**2 / ri**2
                    Am[i] = -cm_i
                    Bm[i] = 1.0 + cp_i + cm_i
                    Cm[i] = -cp_i
                    Dm[i] = cm_i * T[i-1] + (1.0 - cp_i - cm_i) * T[i] + cp_i * T[i+1]
            else:
                # 仍在干燥：Stefan 前沿，最后一个湿核节点
                i_wet = max(0, min(int(r_front / dr), Nr - 1))
                # 湿核节点（0..i_wet）：Dirichlet T = T_evap（Eq. 4.6）
                for i in range(i_wet + 1):
                    Bm[i] = 1.0
                    Dm[i] = t_evap_val
                # 干壳内部节点（i_wet+1..Nr-1）：CN（Eq. 4.1）
                for i in range(i_wet + 1, Nr):
                    ri = r[i]; rp = ri + 0.5 * dr; rm = ri - 0.5 * dr
                    cp_i = sigma * rp**2 / ri**2
                    cm_i = sigma * rm**2 / ri**2
                    Am[i] = -cm_i
                    Bm[i] = 1.0 + cp_i + cm_i
                    Cm[i] = -cp_i
                    Dm[i] = cm_i * T[i-1] + (1.0 - cp_i - cm_i) * T[i] + cp_i * T[i+1]

        # 外表面 Robin BC（Eq. 4.2）——两阶段均使用
        Bm[Nr] = 1.0 + Bi_surf + 2.0 * sigma * (1.0 + 1.0 / Nr)
        Am[Nr] = -2.0 * sigma * (1.0 + 1.0 / Nr)
        Dm[Nr] = (
            T[Nr] * (1.0 - Bi_surf - 2.0 * sigma * (1.0 + 1.0 / Nr))
            + 2.0 * sigma * (1.0 + 1.0 / Nr) * T[Nr - 1]
            + 2.0 * Bi_surf * T_bed
        )

        T_new = _thomas_solve(Am, Bm, Cm, Dm)

        if not drying_phase:
            if T_new[Nr] >= t_evap_val:
                T_new[Nr] = t_evap_val
                drying_phase = True
                r_front = R
        else:
            # Stefan 条件（Eq. 4.3）：前沿热通量驱动蒸发与前沿推进
            if evaporated < total_water:
                i_wet = max(0, min(int(r_front / dr), Nr - 1))
                if i_wet < Nr - 1:
                    # 颗粒内部前沿：干壳导热热通量
                    q_front = LAMBDA_W * (T_new[i_wet + 1] - t_evap_val) / dr
                else:
                    # 前沿紧贴外表面（干壳厚度→0 极限）：外侧对流热通量
                    q_front = h_conv * max(T_bed - t_evap_val, 0.0)
                q_front = max(q_front, 0.0)

                # 蒸发质量（Eq. 4.3 积分）
                A_front = 4.0 * np.pi * max(r_front, dr * 0.5) ** 2
                dm = q_front * A_front * dt / max(h_evap_corr, 1e-12)
                dm = min(dm, total_water - evaporated)
                evaporated += max(dm, 0.0)

                # 前沿向内推进（Eq. 4.3 变形：dr_e/dt = q / (ρ_s w h_v')）
                denom = RHO_W * moisture_mass * max(h_evap_corr, 1e-12)
                r_front = max(r_front - q_front * dt / denom, 0.0)

            # 蒸发完成时强制将前沿归零（物理上颗粒全干，无湿核）
            if evaporated >= total_water:
                r_front = 0.0

            # 仅在仍有未蒸发水分时才锁定湿核温度（Eq. 4.6）
            if evaporated < total_water:
                i_wet_new = max(0, min(int(r_front / dr), Nr - 1))
                T_new[: i_wet_new + 1] = t_evap_val

        T = T_new
        t_arr[n + 1] = (n + 1) * dt
        X_dry_arr[n + 1] = min(evaporated / max(total_water, 1e-30), 1.0)
        T_center_arr[n + 1] = T[0]
        T_surface_arr[n + 1] = T[Nr]
        r_evap_arr[n + 1] = r_front
        if return_history:
            T_history[:, n + 1] = T

    out: dict = {
        "t": t_arr,
        "X_dry": X_dry_arr,
        "T_center": T_center_arr,
        "T_surface": T_surface_arr,
        "r_evap": r_evap_arr,
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
