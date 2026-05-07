"""求解器与预计算：对应流程图中 Vorabrechnung、Zellenmodell/NR、Check2、Ausgabe。

**单一实现**位于 ``src.solvers``；本包仅 **聚合重导出**，便于在
``src/gasifier/`` 命名空间内与论文「Module 2/3」对照，而不复制代码。
"""

from src.solvers.cell_solver import solve_cell
from src.solvers.convergence import InnerConvergence, OuterConvergence
from src.solvers.global_nr_solver import solve_global_nr
from src.solvers.outer_loop import OuterLoopResult, run_global_nr_outer_abgleich
from src.solvers.result_builder import finalize_global_nr_result
from src.solvers.vorabrechnung import (
    estimate_axial_T_profile,
    generate_initial_x0,
    refresh_vorabrechnung_for_cells,
)

__all__ = [
    "estimate_axial_T_profile",
    "generate_initial_x0",
    "refresh_vorabrechnung_for_cells",
    "solve_cell",
    "solve_global_nr",
    "run_global_nr_outer_abgleich",
    "OuterLoopResult",
    "InnerConvergence",
    "OuterConvergence",
    "finalize_global_nr_result",
]
