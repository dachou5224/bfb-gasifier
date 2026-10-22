"""衡算装配、外层预计算与阻尼牛顿。"""

from hamel.solver.driver import solve_case
from hamel.solver.frozen import (
    AssemblyError, FreeboardDomainError, FrozenQuantities, TrialNotEvaluable,
)
from hamel.solver.precalc import (
    PrecalcError, PrecalcResult, downward_mass_flow, initial_precalculation, precalculate,
    require_geldart_b, wall_coefficient,
)
from hamel.solver.residual import AssemblyResidual, evaluate_residual, wall_area_m2
from hamel.solver.vector import (
    CASE1_BED_CELLS, CASE1_FREEBOARD_CELLS, CASE1_SECONDARY_AIR_M, ReactorState,
    UnknownSlot, balance_index, geometric_index, pack_state, secondary_air_balance_index,
    unknown_slots, unpack_state,
)
__all__ = [
    "AssemblyError", "AssemblyResidual", "CASE1_BED_CELLS", "CASE1_FREEBOARD_CELLS",
    "PrecalcError", "PrecalcResult", "downward_mass_flow", "initial_precalculation", "precalculate",
    "require_geldart_b", "solve_case", "wall_coefficient",
    "CASE1_SECONDARY_AIR_M", "FreeboardDomainError", "FrozenQuantities", "ReactorState",
    "TrialNotEvaluable", "UnknownSlot", "balance_index", "evaluate_residual",
    "geometric_index", "pack_state", "secondary_air_balance_index", "unknown_slots",
    "unpack_state", "wall_area_m2",
]
