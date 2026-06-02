"""Arrhenius 速率常数工厂函数（唯一合法入口）。

所有动力学模块必须 import 本模块中的函数，禁止在其他文件中直接写
``k = A * exp(-E/RT)``。

Source of truth:
- ``docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md``
- Hobbs et al. (1992); Jensen et al.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from src.core.constants import Rg


def _clamp_exp_arg(x: float) -> float:
    val = float(x)
    if val < -100.0:
        return -100.0
    if val > 100.0:
        return 100.0
    return val


def _is_scalar_temperature(T: Any) -> bool:
    return not isinstance(T, np.ndarray)


def k_hobbs(k0: float, E: float, T: float) -> float:
    """形式 A：R1–R4（炭气化/燃烧）。

    k = k0 * T * exp(-E / (Rg * T))

    Source: Hobbs et al. (1992); ``docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md``
    """
    if _is_scalar_temperature(T):
        t = float(T)
        return k0 * t * math.exp(_clamp_exp_arg(-E / (Rg * t)))
    return k0 * T * np.exp(np.clip(-E / (Rg * T), -100.0, 100.0))


def k_standard(k0: float, E: float, T: float) -> float:
    """形式 C：标准 Arrhenius（R6 等）。

    k = k0 * exp(-E / (Rg * T))

    Source: ``docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md``
    """
    if _is_scalar_temperature(T):
        t = float(T)
        return k0 * math.exp(_clamp_exp_arg(-E / (Rg * t)))
    return k0 * np.exp(np.clip(-E / (Rg * T), -100.0, 100.0))


def k_jensen_r7(A: float, E_T: float, T: float) -> float:
    """形式 B：R7 专用。

    k = (A / T) * exp(-E_T / T)，其中 E_T = E/Rg [K]

    Source: Jensen et al.; ``docs/hamel_submodels/06_kinetics_r1_r11_and_equilibrium_driving.md``
    """
    if _is_scalar_temperature(T):
        t = float(T)
        return (A / t) * math.exp(_clamp_exp_arg(-E_T / t))
    return (A / T) * np.exp(np.clip(-E_T / T, -100.0, 100.0))
