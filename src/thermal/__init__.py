"""热解子模型：干燥、热解（DAEM）。"""

from .devolatilization import daem_conversion, devolatilization_rate_for_cell
from .drying import drying_rate_for_cell, solve_drying_CN

__all__ = [
    "daem_conversion",
    "devolatilization_rate_for_cell",
    "drying_rate_for_cell",
    "solve_drying_CN",
]
