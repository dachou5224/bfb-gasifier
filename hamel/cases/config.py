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
    """求解器参数。阈值尚未登记 ADR，故不设字段。"""
