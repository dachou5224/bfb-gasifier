"""域层：装置与运行配置 dataclass，以及轴向离散全局状态 GlobalState。"""

from src.gasifier.domain.config import (
    OperatingCondition,
    PlantData,
    SimulationConfig,
    SolidProperties,
)
from src.gasifier.domain.state import GlobalState

__all__ = [
    "GlobalState",
    "SimulationConfig",
    "PlantData",
    "SolidProperties",
    "OperatingCondition",
]
