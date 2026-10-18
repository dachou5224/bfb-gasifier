"""速率的公共类型。负的有限输入不截断。"""
from __future__ import annotations

import math
from dataclasses import dataclass


class KineticsInputError(ValueError):
    """输入不是有限实数，或使印刷式没有定义。不代替速率下限。"""


@dataclass(frozen=True)
class ReactionRate:
    """原式比速率，以及乘上调用方给出的体积或质量后的进度速率。

    ``specific`` 的基准：R1、R2、R4 为单颗粒 mol/s；R3 为 1/s；
    浓度式为 mol/(m³·s)；R8、R9 为 mol/(s·kg 灰)；
    Corella 为 mol/(s·kg 催化剂)；Serio 为 mol/(m³·s)。
    ``extent_mol_s`` 是有符号反应进度速率，mol/s。正值表示反应式正向。
    """

    specific: float
    extent_mol_s: float


def _finite(name: str, value: object) -> float:
    """有限实数。布尔值不是数。"""
    if isinstance(value, bool):
        raise KineticsInputError(f"{name} 必须是有限实数")
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise KineticsInputError(f"{name} 必须是有限实数") from exc
    if not math.isfinite(number):
        raise KineticsInputError(f"{name} 必须是有限实数")
    return number


def _positive(name: str, value: object) -> float:
    """印刷式要做除法或非整数幂时，自变量必须大于 0。"""
    number = _finite(name, value)
    if number <= 0.0:
        raise KineticsInputError(f"{name} 必须大于 0")
    return number


def _iteration_count(name: str, value: object) -> int:
    number = _finite(name, value)
    if number < 1.0 or number != math.trunc(number):
        raise KineticsInputError(f"{name} 必须是不小于 1 的整数")
    return int(number)


def _real_power(name: str, base: object, exponent: float) -> float:
    """负底的非整数幂没有实值。不把底改成 0。"""
    value = _finite(name, base)
    if value < 0.0 and exponent != math.trunc(exponent):
        raise KineticsInputError(f"{name} 为负，不能取非整数幂")
    return value**exponent


def _sphere_area_m2(diameter_m: float) -> float:
    """与颗粒数所用球体同一外表面，A = π d²。"""
    return math.pi * diameter_m * diameter_m


def _particle_count(mass_kg: float, density_kg_m3: float, diameter_m: float) -> float:
    """N = m / (ρ · π d³ / 6)。质量与密度由调用方传入。"""
    volume_m3 = math.pi * diameter_m**3 / 6.0
    denominator = density_kg_m3 * volume_m3
    if denominator == 0.0:
        raise KineticsInputError("颗粒数的分母为 0")
    return mass_kg / denominator
