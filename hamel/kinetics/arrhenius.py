"""全仓唯一的 Arrhenius 工厂。"""
from __future__ import annotations

import math

from hamel.kinetics.common import KineticsInputError, _finite, _positive


def arrhenius_rate_constant(
    pre_factor: float,
    temperature_k: float,
    *,
    activation_energy_j_mol: float | None = None,
    activation_temperature_k: float | None = None,
    gas_constant_j_mol_k: float | None = None,
) -> float:
    """k = A·exp(−E/(R_g·T))。

    表给 E/R_g（K）时只传 ``activation_temperature_k``。表给 E（J/mol）时传
    ``activation_energy_j_mol`` 与 ``gas_constant_j_mol_k``。两者都不给默认值。
    T 的幂、压力因子不在这里。负活化能保持原号。

    Source: Hamel Eq.5.37, S. 90
    """
    factor = _finite("pre_factor", pre_factor)
    temperature = _positive("temperature_k", temperature_k)
    has_energy = activation_energy_j_mol is not None
    has_ratio = activation_temperature_k is not None
    if has_energy == has_ratio:
        raise KineticsInputError("活化能与 E/R_g 必须且只能给出一个")
    if has_ratio:
        if gas_constant_j_mol_k is not None:
            raise KineticsInputError("给出 E/R_g 时不要再传气体常数")
        ratio = _finite("activation_temperature_k", activation_temperature_k)
        exponent = -ratio / temperature
    else:
        energy = _finite("activation_energy_j_mol", activation_energy_j_mol)
        gas_constant = _positive("gas_constant_j_mol_k", gas_constant_j_mol_k)
        exponent = -energy / (gas_constant * temperature)
    return factor * math.exp(exponent)
