"""预计算、外层与阻尼牛顿。本卡只装配残差。"""

from hamel.solver.frozen import (
    AssemblyError, FreeboardDomainError, FrozenQuantities, TrialNotEvaluable,
)
from hamel.solver.residual import AssemblyResidual, evaluate_residual, wall_area_m2
from hamel.solver.vector import (
    CASE1_BED_CELLS, CASE1_FREEBOARD_CELLS, CASE1_SECONDARY_AIR_M, ReactorState,
    UnknownSlot, balance_index, geometric_index, pack_state, secondary_air_balance_index,
    unknown_slots, unpack_state,
)
__all__ = [
    "AssemblyError", "AssemblyResidual", "CASE1_BED_CELLS", "CASE1_FREEBOARD_CELLS",
    "CASE1_SECONDARY_AIR_M", "FreeboardDomainError", "FrozenQuantities", "ReactorState",
    "TrialNotEvaluable", "UnknownSlot", "balance_index", "evaluate_residual",
    "geometric_index", "pack_state", "secondary_air_balance_index", "unknown_slots",
    "unpack_state", "wall_area_m2",
]
