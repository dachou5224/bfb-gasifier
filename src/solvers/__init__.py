"""求解器：单 cell fsolve + 多 cell 扫描迭代 / 可选全局 NR。"""

from .cell_solver import solve_cell
from .convergence import InnerConvergence, OuterConvergence
from .global_nr_solver import solve_global_nr
from .outer_loop import OuterLoopResult, run_global_nr_outer_abgleich

__all__ = [
    "solve_cell",
    "solve_global_nr",
    "run_global_nr_outer_abgleich",
    "OuterLoopResult",
    "InnerConvergence",
    "OuterConvergence",
]
