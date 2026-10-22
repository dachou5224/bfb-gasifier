"""物理参数与求解器参数。字段名与 data/cases/*.json 的键一致。"""
from __future__ import annotations

from dataclasses import dataclass

@dataclass(frozen=True)
class PhysicalConfig:
    """已装配的物理参数，数值为 SI。不设默认值。"""

    bed_diameter: float
    freeboard_diameter: float
    # 解释性选择：ADR 0005。布风板至自由空间顶部，3.5 m。
    reactor_height: float
    # 解释性选择：ADR 0005。床区筒体高度，不是操作床高。
    bed_zone_height: float
    # 解释性选择：ADR 0005。
    freeboard_height: float
    # 解释性选择：ADR 0005。
    fuel_feed_height: float
    cyclone_cell_count: int
    fuel_feed_count: int
    agent_feed_count: int
    ash_mass_fraction_wf: float
    volatiles_mass_fraction_waf: float
    carbon_mass_fraction_wf: float
    hydrogen_mass_fraction_wf: float
    oxygen_mass_fraction_wf: float
    nitrogen_mass_fraction_wf: float
    sulfur_mass_fraction_wf: float
    # 解释性选择：ADR 0005。
    fixed_carbon_mass_fraction_wf: float
    # 解释性选择：ADR 0005。
    particle_size_max: float
    fuel_mass_flow: float
    moisture_mass_fraction: float
    lambda_o2: float
    # mol/s。由 primary_air_volume_flow 按 ADR 0004 换算。
    primary_air_amount_flow: float
    # mol/s。由 secondary_air_volume_flow 按 ADR 0004 换算。
    secondary_air_amount_flow: float
    # 值的单位是 mol/s。JSON 中此键为 m³(V_N)/h，装配时按 ADR 0004 换算。
    # 解释性选择：ADR 0006。仅工况 1、2 计入模型。
    nitrogen_purge_volume_flow: float
    steam_mass_flow: float
    # 解释性选择：ADR 0005。仅工况 1、3、7、11、14。
    inlet_temperature: float
    pressure: float
    recirculation: bool
    # 解释性选择：ADR 0005。收到基低位热值，J/kg。仅工况 1、3、7、11、14。
    heating_value: float

@dataclass(frozen=True)
class SolverConfig:
    """求解器参数。热解离散见 ADR 0027，自由空间轨迹见 ADR 0028，颗粒温度见 ADR 0033，悬停邻域见 ADR 0034，驱动字段见 ADR 0037，敏感性字段见 ADR 0039。"""

    pyrolysis_radial_intervals: int = 24
    pyrolysis_time_rtol: float = 1.0e-6
    pyrolysis_activation_nodes: int = 16
    pyrolysis_tail_delta_t_k: float = 0.02
    pyrolysis_initial_shell_fraction: float = 1.0e-3
    # 正半轴等概率份数。见 ADR 0028。邻域 δ 随悬停记账另登记。
    freeboard_ejection_samples: int = 48
    # 积分步长上限，s。见 ADR 0028。
    freeboard_max_step_s: float = 2.0e-4
    # 只用于报错。到时仍未离开自由空间就失败。见 ADR 0028。
    freeboard_time_limit_s: float = 80.0
    # 确定悬停归宿的高度邻域，m。记账与它无关。见 ADR 0034。
    freeboard_delta_h_m: float = 1.0e-3
    # 确定悬停归宿的速度邻域，m/s。见 ADR 0034。
    freeboard_delta_u_m_s: float = 1.0e-2
    particle_tolerance_w: float = 1.0e-4
    particle_max_iterations: int = 25
    coupling_tolerance_k: float = 1.0e-3
    coupling_max_iterations: int = 16
    # 床区几何积分的起始子步与上限。求解器参数：ADR 0036。
    bed_integration_intervals_per_cell: int = 8
    bed_integration_max_intervals_per_cell: int = 256
    # 失败判据，不是成功判据。ADR 0017 第 8 项；字段名见 ADR 0037。
    nu_max: int = 200
    outer_round_limit: int = 100
    # i ≤ 30。ADR 0017 第 7 项。
    damping_index_limit: int = 30
    # ADR 0012。元素 0,1 %，能量 0,5 %。
    residual_tolerance: float = 1.0e-6
    outer_temperature_tolerance_k: float = 0.1
    outer_relative_tolerance: float = 1.0e-4
    element_balance_tolerance: float = 1.0e-3
    energy_balance_tolerance: float = 5.0e-3
    # 单次运行的墙钟上限，s。依据见 ADR 0037。
    wall_timeout_s: float = 90_000.0
    # 只选择已批准的固定对，不接受任意数值。ADR 0039。
    particle_emissivity_case: str = "frozen"
    interphase_diffusivity: str = "hamel_5_20"
    hover_closure: str = "A"
    serio_pair: str = "first"


# 已批准的 ε_P 与 Serio 对。不是新的可调物理参数。ADR 0017 第 17 项；表 5.4，S. 97。
_EMISSIVITY = {"frozen": 0.85, "S3": 0.79, "S4": 0.91}
_SERIO = {"first": (5.42e4, 100.5e3), "second": (9.96e2, 64.5e3)}


def approved_emissivity(choice: str) -> float:
    """返回冻结值或 S3、S4 这一对里的一个。其它字符串拒绝。"""
    try:
        return _EMISSIVITY[choice]
    except KeyError as exc:
        raise ValueError("ε_P 只能取已批准的 frozen、S3 或 S4") from exc


def approved_serio(choice: str) -> tuple[float, float]:
    """返回表 5.4 的第一对或第二对，单位 s^-1 与 J/mol。"""
    try:
        return _SERIO[choice]
    except KeyError as exc:
        raise ValueError("Serio 参数只能取已批准的 first 或 second") from exc
