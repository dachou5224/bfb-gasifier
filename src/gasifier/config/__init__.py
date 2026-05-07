"""配置层兼容导出（当前复用 domain.config）。"""

from src.gasifier.domain.config import OperatingCondition, PlantData, SimulationConfig, SolidProperties

__all__ = ["SimulationConfig", "PlantData", "SolidProperties", "OperatingCondition"]
