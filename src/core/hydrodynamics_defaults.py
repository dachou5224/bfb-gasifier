"""Hamel (1999) 床层水力链默认配置（thesis / 非 thesis 共用）。

Source: docs/hamel_submodels/03_hydrodynamics_core_chain.md
       Hamel (1999) Eq.3.11–3.15, Eq.3.41–3.44, Eq.3.50
"""

from __future__ import annotations

from typing import Protocol

# Hamel 论文主链（Wein u_d + Hilligardt/Heinbockel 加压气泡 ODE + Heinbockel u_b）
HAMEL_HYDRO_U_D_CLOSURE: str = "wein_1992_eq312"
HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL: str = "hilligardt_ode"
HAMEL_HYDRO_PSI_B_STRATEGY: str = "wein_1992"
HAMEL_HYDRO_LAMBDA_STRATEGY: str = "hamel_280"
HAMEL_HYDRO_XI_STRATEGY: str = "hamel_regime"
HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY: str = "heinbockel_eq343"
HAMEL_HYDRO_BUBBLE_ODE_STRATEGY: str = "heinbockel_eq341"
HAMEL_HYDRO_V_B: float = 1.0  # [-] Eq. (3.43) 单泡速度权重 v_b

# 工程/审计对照分支（须显式 opt-in，不再作为默认）
LEGACY_HYDRO_U_D_CLOSURE: str = "current"
LEGACY_HYDRO_BUBBLE_DIAMETER_MODEL: str = "mori_wen"
LEGACY_HYDRO_BUBBLE_VELOCITY_STRATEGY: str = "hilligardt_eq313"
LEGACY_HYDRO_BUBBLE_ODE_STRATEGY: str = "hilligardt_eq333"


class _HydroCellLike(Protocol):
    u_d_closure: str
    bubble_diameter_model: str
    psi_b_strategy: str
    lambda_strategy: str
    xi_strategy: str
    bubble_velocity_strategy: str
    bubble_ode_strategy: str
    hydrodynamics_v_b: float


def apply_hamel_hydrodynamics_to_cell(cell: _HydroCellLike) -> None:
    """将 cell 水力策略设为 Hamel 主链（thesis_mode 安全锁与默认构建共用）。"""
    cell.u_d_closure = HAMEL_HYDRO_U_D_CLOSURE
    cell.bubble_diameter_model = HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL
    cell.psi_b_strategy = HAMEL_HYDRO_PSI_B_STRATEGY
    cell.lambda_strategy = HAMEL_HYDRO_LAMBDA_STRATEGY
    cell.xi_strategy = HAMEL_HYDRO_XI_STRATEGY
    cell.bubble_velocity_strategy = HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY
    cell.bubble_ode_strategy = HAMEL_HYDRO_BUBBLE_ODE_STRATEGY
    cell.hydrodynamics_v_b = float(HAMEL_HYDRO_V_B)
