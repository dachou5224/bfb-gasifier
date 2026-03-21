"""热力学模块：Gibbs 平衡与可逆反应驱动力。

- equilibrium: K_eq、Q_p 计算（Phase 1）
- gibbs_minimizer: Gibbs 自由焓最小化（Phase 2）
- minor_species: 微量组分分布接口（Phase 2）

Source: BFB_TechSpec_v11 §5, §5.6; Hamel (1999) 附录 A1
"""

from src.thermodynamics.equilibrium import (
    calc_gibbs_driving_force,
    calc_reaction_quotient,
    get_K_eq,
)
from src.thermodynamics.gibbs_minimizer import GibbsMinimizer
from src.thermodynamics.minor_species import (
    solve_minor_species,
    solve_nitrogen_distribution,
    solve_sulfur_distribution,
)

__all__ = [
    "get_K_eq",
    "calc_reaction_quotient",
    "calc_gibbs_driving_force",
    "GibbsMinimizer",
    "solve_sulfur_distribution",
    "solve_nitrogen_distribution",
    "solve_minor_species",
]
