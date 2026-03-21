"""热解模型（DAEM，Gauss-Hermite 数值积分）。

基于 Hamel 论文 Chapter 4（Eq. 4.10–4.12）：
1) 活化能服从 Gaussian 分布；
2) 每个活化能通道为一级反应；
3) 当颗粒内部存在径向温度梯度时，对半径做体积分加权。

验证指标：升温至 900°C，积分 100 s，挥发分释放率 > 80%
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from src.core.constants import Rg

# -----------------------------------------------------------------------
# DAEM 参数（Hamel Chapter 4, Table 4.1/4.2）
# -----------------------------------------------------------------------
@dataclass(frozen=True)
class DAEMKinetics:
    """DAEM 动力学参数。"""
    A: float        # [1/s] 频率因子 k0
    E0: float       # [J/mol] 平均活化能
    sigma: float    # [J/mol] 活化能标准差


DAEM_PARAMS = {
    # Table 4.1/4.2
    "brown_coal": DAEMKinetics(A=1.67e13, E0=192_000.0, sigma=40_000.0),
    "wood": DAEMKinetics(A=2_500.0, E0=67_500.0, sigma=13_500.0),
}

# 默认使用褐煤（与现有工况一致）
A_DAEM: float = DAEM_PARAMS["brown_coal"].A
E0_DAEM: float = DAEM_PARAMS["brown_coal"].E0
SIGMA_DAEM: float = DAEM_PARAMS["brown_coal"].sigma
VM_DAF_DEFAULT: float = 0.5342  # [-] daf 基挥发分含量（LU 工况）


# -----------------------------------------------------------------------
# Gauss-Hermite 积分节点与权重
# -----------------------------------------------------------------------

def _gauss_hermite_nodes(n: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """返回 n 阶 Gauss-Hermite 积分的节点和权重。

    Source: numpy.polynomial.hermite_e
    """
    xi, wi = np.polynomial.hermite_e.hermegauss(n)
    return xi, wi


# -----------------------------------------------------------------------
# DAEM 核心求解
# -----------------------------------------------------------------------

def _single_reaction_conversion(
    T_history: np.ndarray,
    t_history: np.ndarray,
    E: float,
    A: float = A_DAEM,
) -> float:
    """给定活化能 E 下，单个一级反应在 T(t) 温度史下的转化率。

    X(E) = 1 - exp(-A * integral_0^t exp(-E/(Rg*T(t'))) dt')

    用梯形法则数值积分。

    Source: DAEM 理论; Pitt (1962); Anthony & Howard (1976)
    """
    integrand = A * np.exp(np.clip(-E / (Rg * T_history), -200.0, 200.0))
    integral = np.trapezoid(integrand, t_history)
    return 1.0 - np.exp(-integral)


def gaussian_energy_pdf(E: np.ndarray | float, E0: float, sigma: float) -> np.ndarray | float:
    """Gaussian 活化能分布函数 f(E)（Eq. 4.11）。"""
    coef = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    return coef * np.exp(-((E - E0) ** 2) / (2.0 * sigma**2))


def daem_conversion(
    T_history: np.ndarray,
    t_history: np.ndarray,
    A: float = A_DAEM,
    E0: float = E0_DAEM,
    sigma: float = SIGMA_DAEM,
    n_quad: int = 10,
) -> float:
    """DAEM 整体挥发分转化率 X_VM [-]。

    X_VM = integral_{-inf}^{+inf} X(E) * f(E) dE

    其中 f(E) = 1/(sigma*sqrt(2*pi)) * exp(-(E-E0)^2/(2*sigma^2))
    用 Gauss-Hermite 积分：E_j = E0 + sqrt(2)*sigma*xi_j

    Parameters
    ----------
    T_history : 温度历史数组 [K]
    t_history : 时间历史数组 [s]
    A         : 频率因子 [1/s]
    E0        : 平均活化能 [J/mol]
    sigma     : 活化能标准差 [J/mol]
    n_quad    : Gauss-Hermite 积分阶数

    Returns
    -------
    X_VM : 0-1 之间的挥发分释放率

    Source: specs/03_drying_devolatilization.md §2; Anthony & Howard (1976)
    """
    xi, wi = _gauss_hermite_nodes(n_quad)

    # hermegauss 返回概率论 Hermite 节点/权重：
    # ∫ f(x) exp(-x²/2) dx ≈ Σ w_j f(x_j)
    # 令 E = E0 + sigma*x => f(E) dE = (1/sqrt(2π)) exp(-x²/2) dx
    X_total = 0.0
    for j in range(n_quad):
        E_j = E0 + sigma * xi[j]
        if E_j < 0:
            X_j = 1.0
        else:
            X_j = _single_reaction_conversion(T_history, t_history, E_j, A)
        X_total += wi[j] * X_j

    X_total /= np.sqrt(2.0 * np.pi)
    return float(np.clip(X_total, 0.0, 1.0))


def daem_conversion_radial(
    T_history_rt: np.ndarray,
    t_history: np.ndarray,
    r_nodes: np.ndarray,
    A: float = A_DAEM,
    E0: float = E0_DAEM,
    sigma: float = SIGMA_DAEM,
    n_quad: int = 10,
) -> float:
    """带径向体积分的 DAEM 挥发分释放率（Eq. 4.12）。

    Parameters
    ----------
    T_history_rt : shape=(Nr, Nt) 的温度历史矩阵 [K]
    t_history    : 时间数组 [s]，长度 Nt
    r_nodes      : 半径节点 [m]，长度 Nr
    """
    T_history_rt = np.asarray(T_history_rt, dtype=float)
    t_history = np.asarray(t_history, dtype=float)
    r_nodes = np.asarray(r_nodes, dtype=float)

    if T_history_rt.ndim != 2:
        raise ValueError("T_history_rt 必须为二维数组 (Nr, Nt)")
    if T_history_rt.shape[1] != t_history.size:
        raise ValueError("T_history_rt 的时间维必须与 t_history 一致")
    if T_history_rt.shape[0] != r_nodes.size:
        raise ValueError("T_history_rt 的径向维必须与 r_nodes 一致")
    if r_nodes.size < 2:
        raise ValueError("r_nodes 至少包含两个节点")

    # 逐半径点求局部 X_VM(r)
    X_r = np.array(
        [
            daem_conversion(
                T_history=T_history_rt[i, :],
                t_history=t_history,
                A=A,
                E0=E0,
                sigma=sigma,
                n_quad=n_quad,
            )
            for i in range(r_nodes.size)
        ]
    )

    R0 = float(np.max(r_nodes))
    if R0 <= 0:
        return 0.0

    # Eq. 4.12: 3/R0^3 ∫ X(r) r^2 dr
    integral = np.trapezoid(X_r * (r_nodes**2), r_nodes)
    X_vol = 3.0 * integral / (R0**3)
    return float(np.clip(X_vol, 0.0, 1.0))


# -----------------------------------------------------------------------
# 温度程序生成工具
# -----------------------------------------------------------------------

def _linear_heating_profile(
    T_start: float,
    T_end: float,
    t_total: float,
    n_points: int = 200,
) -> Tuple[np.ndarray, np.ndarray]:
    """生成线性升温温度历史。"""
    t = np.linspace(0, t_total, n_points)
    T = T_start + (T_end - T_start) * t / t_total
    return T, t


# -----------------------------------------------------------------------
# 稳态 cell 模型接口
# -----------------------------------------------------------------------

def devolatilization_rate_for_cell(
    T_bed: float,
    tau_cell: float,
    T_init: float = 400.0,
    VM_daf: float = VM_DAF_DEFAULT,
    A: float | None = None,
    E0: float | None = None,
    sigma: float | None = None,
    fuel_type: str = "brown_coal",
) -> Tuple[float, float]:
    """计算颗粒在 cell 停留时间 tau_cell 内的热解进度。

    假设颗粒从 T_init 线性升温至 T_bed（简化热传导）。

    Parameters
    ----------
    T_bed    : 床层温度 [K]
    tau_cell : 平均停留时间 [s]
    T_init   : 颗粒到达 cell 时的温度 [K]
    VM_daf   : daf 基挥发分质量分数 [-]
    A, E0, sigma : DAEM 参数；若为 None 则按 fuel_type 取 Table 4.1/4.2
    fuel_type    : "brown_coal" 或 "wood"

    Returns
    -------
    X_VM : 挥发分释放率 (0-1)
    devol_rate : [kg_VM/(kg_daf·s)] 平均热解速率

    Source: docs/CLAUDE.md Phase 4.2
    """
    if A is None or E0 is None or sigma is None:
        kin = DAEM_PARAMS.get(fuel_type, DAEM_PARAMS["brown_coal"])
        A = kin.A if A is None else A
        E0 = kin.E0 if E0 is None else E0
        sigma = kin.sigma if sigma is None else sigma

    T_history, t_history = _linear_heating_profile(
        T_start=T_init,
        T_end=T_bed,
        t_total=max(tau_cell, 0.01),
        n_points=200,
    )

    X_VM = daem_conversion(T_history, t_history, A, E0, sigma)
    devol_rate = X_VM * VM_daf / max(tau_cell, 1e-10)

    return X_VM, devol_rate
