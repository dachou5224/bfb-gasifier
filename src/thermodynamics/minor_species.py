"""微量组分 Gibbs 分布求解接口。

调用 gibbs_minimizer 求解 H2S、SO2、COS、NH3、HCN、NO 等微量组分的平衡分布。
Source: BFB_TechSpec_v11 §5.2, §5.5
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import src.core.species as species_module
from src.thermodynamics.gibbs_minimizer import GibbsMinimizer


class SpeciesDBAdapter:
    """基于 species.py 的 SpeciesDB 适配器。"""

    def get_mu0(self, species: str, T: float) -> float:
        return species_module.gibbs_molar(species, T)

    def get_atom_count(self, species: str, element: str) -> int:
        return species_module.get_atom_count(species, element)


# 硫系候选组分
SULFUR_CANDIDATES: List[str] = ["H2S", "SO2", "COS"]

# 氮系候选组分
NITROGEN_CANDIDATES: List[str] = ["NH3", "HCN", "NO"]

# 默认微量组分求解器实例
_default_minimizer: GibbsMinimizer | None = None


def _get_minimizer() -> GibbsMinimizer:
    global _default_minimizer
    if _default_minimizer is None:
        _default_minimizer = GibbsMinimizer(SpeciesDBAdapter())
    return _default_minimizer


def solve_sulfur_distribution(
    T: float,
    P: float,
    elements: Dict[str, float],
    candidates: List[str] | None = None,
    lambda0: list[float] | None = None,
    ln_N0: float | None = None,
    return_diag: bool = False,
) -> Dict[str, float] | tuple[Dict[str, float], Dict[str, Any]]:
    """求解硫系微量组分（H2S、SO2、COS）的 Gibbs 平衡分布。"""
    minimizer = _get_minimizer()
    cand = candidates or SULFUR_CANDIDATES
    return minimizer.solve(
        T,
        P,
        elements,
        cand,
        lambda0=None if lambda0 is None else np.array(lambda0, dtype=float),
        ln_N0=ln_N0,
        return_diag=return_diag,
    )


def solve_nitrogen_distribution(
    T: float,
    P: float,
    elements: Dict[str, float],
    candidates: List[str] | None = None,
    lambda0: list[float] | None = None,
    ln_N0: float | None = None,
    return_diag: bool = False,
) -> Dict[str, float] | tuple[Dict[str, float], Dict[str, Any]]:
    """求解氮系微量组分（NH3、HCN、NO）的 Gibbs 平衡分布。"""
    minimizer = _get_minimizer()
    cand = candidates or NITROGEN_CANDIDATES
    return minimizer.solve(
        T,
        P,
        elements,
        cand,
        lambda0=None if lambda0 is None else np.array(lambda0, dtype=float),
        ln_N0=ln_N0,
        return_diag=return_diag,
    )


def solve_minor_species(
    T: float,
    P: float,
    elements: Dict[str, float],
    candidates: List[str],
    lambda0: list[float] | None = None,
    ln_N0: float | None = None,
    return_diag: bool = False,
) -> Dict[str, float] | tuple[Dict[str, float], Dict[str, Any]]:
    """通用微量组分 Gibbs 最小化求解。"""
    return _get_minimizer().solve(
        T,
        P,
        elements,
        candidates,
        lambda0=None if lambda0 is None else np.array(lambda0, dtype=float),
        ln_N0=ln_N0,
        return_diag=return_diag,
    )
