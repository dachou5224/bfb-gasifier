"""Arrhenius 速率常数工厂函数（唯一合法入口）。

所有动力学模块必须 import 本模块中的函数，禁止在其他文件中直接写
k = A * exp(-E/RT)。

Source: docs/CLAUDE.md §全局规则 2; specs/04_kinetics.md
"""

from __future__ import annotations

import numpy as np

Rg: float = 8.314  # [J/(mol·K)]


def k_hobbs(k0: float, E: float, T: float) -> float:
    """形式 A：R1–R4（炭气化/燃烧）。

    k = k0 * T * exp(-E / (Rg * T))

    Source: Hobbs et al. (1992); specs/04_kinetics.md R1–R4
    """
    return k0 * T * np.exp(np.clip(-E / (Rg * T), -100.0, 100.0))


def k_standard(k0: float, E: float, T: float) -> float:
    """形式 C：标准 Arrhenius（R6 等）。

    k = k0 * exp(-E / (Rg * T))

    Source: specs/04_kinetics.md
    """
    return k0 * np.exp(np.clip(-E / (Rg * T), -100.0, 100.0))


def k_jensen_r7(A: float, E_T: float, T: float) -> float:
    """形式 B：R7 专用。

    k = (A / T) * exp(-E_T / T)，其中 E_T = E/Rg [K]

    Source: Jensen et al.; specs/04_kinetics.md R7
    """
    return (A / T) * np.exp(np.clip(-E_T / T, -100.0, 100.0))
