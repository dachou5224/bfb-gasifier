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

from src.core.constants import Rg
from src.kinetics.arrhenius import k_standard

# ... (other code)

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
    # integrand = A * np.exp(np.clip(-E / (Rg * T_history), -200.0, 200.0))
    integrand = k_standard(A, E, T_history)
    integral = np.trapz(integrand, t_history)
    return 1.0 - np.exp(np.clip(-integral, -200.0, 0.0))


def gaussian_energy_pdf(E: np.ndarray | float, E0: float, sigma: float) -> np.ndarray | float:
    """Gaussian 活化能分布函数 f(E)（Eq. 4.11）。"""
    coef = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    return coef * np.exp(np.clip(-((E - E0) ** 2) / (2.0 * sigma**2), -200.0, 0.0))


def daem_conversion(
    T_history: np.ndarray,
    t_history: np.ndarray,
    A: float = A_DAEM,
    E0: float = E0_DAEM,
    sigma: float = SIGMA_DAEM,
    n_quad: int = 10,
) -> float | np.ndarray:
    """DAEM 整体挥发分转化率 X_VM [-]。支持向量化（T_history 为 2D 时）。"""
    xi, wi = _gauss_hermite_nodes(n_quad)
    
    T_history = np.asarray(T_history)
    t_history = np.asarray(t_history)
    
    # 维度处理：若 T_history 是 1D，扩展为 (1, Nt)
    is_1d = (T_history.ndim == 1)
    if is_1d:
        T_history = T_history[np.newaxis, :]
    
    Nr, Nt = T_history.shape
    X_total = np.zeros(Nr)
    
    for j in range(n_quad):
        E_j = E0 + sigma * xi[j]
        if E_j < 0:
            X_j = np.ones(Nr)
        else:
            # 向量化计算 k(t)
            # k = A * exp(-E_j / (Rg * T_history))
            k = k_standard(A, E_j, T_history) # shape (Nr, Nt)
            # 梯形积分
            integral = np.trapz(k, t_history, axis=1) # shape (Nr,)
            X_j = 1.0 - np.exp(np.clip(-integral, -200.0, 0.0))
        
        X_total += wi[j] * X_j

    X_total /= np.sqrt(2.0 * np.pi)
    res = np.clip(X_total, 0.0, 1.0)
    return float(res[0]) if is_1d else res


def daem_conversion_radial(
    T_history_rt: np.ndarray,
    t_history: np.ndarray,
    r_nodes: np.ndarray,
    A: float = A_DAEM,
    E0: float = E0_DAEM,
    sigma: float = SIGMA_DAEM,
    n_quad: int = 10,
) -> float:
    """带径向体积分的 DAEM 挥发分释放率（Eq. 4.12）。向量化版本。"""
    T_history_rt = np.asarray(T_history_rt, dtype=float)
    t_history = np.asarray(t_history, dtype=float)
    r_nodes = np.asarray(r_nodes, dtype=float)

    # 直接调用向量化的 daem_conversion
    X_r = daem_conversion(T_history_rt, t_history, A, E0, sigma, n_quad)
    assert isinstance(X_r, np.ndarray)

    R0 = float(np.max(r_nodes))
    if R0 <= 0:
        return 0.0

    # Eq. 4.12: 3/R0^3 ∫ X(r) r^2 dr
    integral = np.trapz(X_r * (r_nodes**2), r_nodes)
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
