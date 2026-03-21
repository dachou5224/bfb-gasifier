"""求解器：单 cell fsolve + 多 cell 扫描迭代 / 可选全局 NR。"""

from .cell_solver import solve_cell
from .global_nr_solver import solve_global_nr

__all__ = ["solve_cell", "solve_global_nr"]
