"""域层：装置与运行配置 dataclass，以及 Phase0 ``GlobalState``（已弃用）。

配置 dataclass 仍供 ``gasifier.models.correlations`` 与 Phase0 演示引用。
生产状态与求解：``src.core.reactor.Reactor``。
"""

from __future__ import annotations

import importlib
import warnings

from src.gasifier.domain.config import (
    OperatingCondition,
    PlantData,
    SimulationConfig,
    SolidProperties,
)

__all__ = [
    "GlobalState",
    "SimulationConfig",
    "PlantData",
    "SolidProperties",
    "OperatingCondition",
]


def __getattr__(name: str):
    if name != "GlobalState":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    warnings.warn(
        "src.gasifier.domain.GlobalState is Phase0 demo only; "
        "use src.core.reactor.Reactor for production paths.",
        DeprecationWarning,
        stacklevel=2,
    )
    mod = importlib.import_module("src.gasifier.domain.state")
    return mod.GlobalState


def __dir__() -> list[str]:
    return sorted(__all__)
