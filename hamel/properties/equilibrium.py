"""附录 A1 的自由焓最小化，以及气相绝热平衡。"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from hamel.properties.data import (
    GAS_SPECIES,
    EquilibriumError,
    PropertyError,
    atoms,
    gas_constant_j_mol_k,
    is_gas,
    standard_pressure_pa,
)
from hamel.properties.thermo import chemical_potential_j_mol, mixture_enthalpy_j

# ADR 0017 第 16 项。SolverConfig 目前没有字段，数值放在这里。
_ETA0 = 0.5
_MAX_HALVINGS = 30
_RESIDUAL_TOL = 1e-10
_MAX_ITERATIONS = 200


@dataclass(frozen=True)
class EquilibriumResult:
    """平衡组成。被移除的组分摩尔分数为 0。"""
    mole_fractions: dict[str, float]
    amounts_mol: dict[str, float]
    total_mol: float
    temperature_k: float
    pressure_pa: float
    iterations: int
    dropped: tuple[str, ...]


def is_converged(max_abs_residual: float, balance_deviation: float) -> bool:
    """式 (A.26) 的每一行和式 (A.34) 必须同时满足。ADR 0017 第 16 项。 Source: Hamel Eq.A.26, S. 166 Source: Hamel Eq.A.34, S. 167"""
    return max_abs_residual <= _RESIDUAL_TOL and balance_deviation <= _RESIDUAL_TOL
def _positive_elements(amounts: dict[str, float]) -> dict[str, float]:
    cleaned = {}
    for element, amount in amounts.items():
        value = float(amount)
        if value < 0:
            raise PropertyError(f"元素 {element} 的量为负")
        if value > 0:
            cleaned[element] = value
    return cleaned
def _prepare(element_amounts: dict[str, float], species: list[str], sources: dict[str, float]):
    """去掉含缺失元素的组分，并删掉该元素的方程。ADR 0017 第 16 项。"""
    present = _positive_elements(element_amounts)
    for name, amount in sources.items():
        if float(amount) < 0:
            raise PropertyError(f"源汇物质 {name} 的量为负")
        for element, count in atoms(name).items():
            present[element] = present.get(element, 0.0) - count * float(amount)
            if present[element] < -1e-9:
                raise EquilibriumError(f"源汇物质带走的 {element} 多于进料")
            if present[element] < 0:
                present[element] = 0.0
    missing = {element for element, amount in present.items() if amount == 0.0}
    # 进料里没有、但源汇也没带来的元素不在 present 中。组分若含这种元素，也视为缺失。
    known = set(present)
    dropped = []
    active = []
    for name in species:
        if name in sources:
            continue
        composition = atoms(name)
        if any(count > 0 and (element not in known or element in missing) for element, count in composition.items()):
            dropped.append(name)
            continue
        active.append(name)
    elements = [element for element, amount in present.items() if amount > 0.0]
    return elements, {element: present[element] for element in elements}, active, tuple(dropped)
def _c_values(names: list[str], temperature_k: float, pressure_pa: float) -> np.ndarray:
    """式 (A.12)：c_j = μ^0/(RT) + ln(P/P_0)。μ^0 已换到 P_0。"""
    gas_constant = gas_constant_j_mol_k()
    pressure_term = math.log(float(pressure_pa) / standard_pressure_pa())
    values = []
    for name in names:
        mu = chemical_potential_j_mol(name, temperature_k)
        values.append(mu / (gas_constant * temperature_k) + pressure_term)
    return np.asarray(values, dtype=float)
def _stoichiometry(names: list[str], elements: list[str]) -> np.ndarray:
    matrix = np.zeros((len(elements), len(names)))
    for column, name in enumerate(names):
        composition = atoms(name)
        for row, element in enumerate(elements):
            matrix[row, column] = composition.get(element, 0)
    return matrix
def _q_matrix(stoichiometry: np.ndarray, b_plus: np.ndarray, source_moles: float) -> tuple[np.ndarray, np.ndarray]:
    """式 (A.26) 下的 q_ij、d_i。指数留在求和之内。ADR 0017 第 15 项。"""
    n_elements, n_species = stoichiometry.shape
    q_matrix = np.zeros_like(stoichiometry)
    d_values = np.zeros(n_elements)
    reference = stoichiometry[0]
    q_matrix[0] = 1.0 + reference * source_moles / b_plus[0]
    d_values[0] = 1.0
    for row in range(1, n_elements):
        q_matrix[row] = stoichiometry[row] - (b_plus[row] / b_plus[0]) * reference
    return q_matrix, d_values
def _total_moles(stoichiometry: np.ndarray, b_plus: np.ndarray, mole_fraction: np.ndarray) -> float:
    """式 (A.33)。"""
    denominator = float(np.sum(stoichiometry @ mole_fraction))
    if denominator <= 0:
        raise EquilibriumError("式 (A.33) 的分母不是正数")
    return float(np.sum(b_plus) / denominator)


def balance_deviation(stoichiometry: np.ndarray, b_plus: np.ndarray, amounts: np.ndarray) -> float:
    """式 (A.34) 的 D。"""
    total = 0.0
    for row, element_amount in enumerate(b_plus):
        total += abs(1.0 - float(stoichiometry[row] @ amounts) / element_amount)
    return total


def initial_lambdas(
    stoichiometry: np.ndarray,
    b_plus: np.ndarray,
    c_values: np.ndarray,
    mole_fraction: np.ndarray,
    total_moles: float,
) -> np.ndarray:
    """式 (A.40)，右端第一项取 b_i^+/N。ADR 0017 第 15 项。 Source: Hamel Eq.A.40, S. 168"""
    weight = mole_fraction
    left = stoichiometry * weight @ stoichiometry.T
    log_term = weight * (1.0 - np.log(weight))
    right = (
        b_plus / total_moles
        - stoichiometry @ log_term
        + stoichiometry @ (weight * c_values)
    )
    return np.linalg.solve(left, right)
def _residuals(q_matrix: np.ndarray, d_values: np.ndarray, mole_fraction: np.ndarray) -> np.ndarray:
    """式 (A.26) 的 F = Σ q exp − d。"""
    return q_matrix @ mole_fraction - d_values
def _exponent_feasible(stoichiometry: np.ndarray, lambdas: np.ndarray, c_values: np.ndarray) -> bool:
    """式 (A.31)。同页写明对应 0 ≤ y_j ≤ 1，故指数 ≤ 0，允许 y_j = 1。"""
    exponent = stoichiometry.T @ lambdas - c_values
    return bool(np.all(exponent <= 1e-8))


def gibbs_minimization(
    temperature_k: float,
    pressure_pa: float,
    element_amounts_mol: dict[str, float],
    species: list[str] | None = None,
    source_amounts_mol: dict[str, float] | None = None,
) -> EquilibriumResult:
    """给定 T、P、元素量与源汇量，返回平衡组成。 Source: Hamel Eq.A.26, S. 166 Source: Hamel Eq.A.15, S. 165"""
    temperature = float(temperature_k)
    pressure = float(pressure_pa)
    if temperature <= 0 or pressure <= 0:
        raise PropertyError("温度和压力必须为正")
    requested = list(GAS_SPECIES if species is None else species)
    for name in requested:
        if not is_gas(name):
            raise PropertyError(f"{name} 不是气相组分，不能进入附录 A1")
    sources = {name: float(amount) for name, amount in (source_amounts_mol or {}).items() if amount}
    elements, b_map, active, dropped = _prepare(element_amounts_mol, requested, sources)
    if not active:
        if elements:
            raise EquilibriumError("有元素量，但含这些元素的组分都被移除了")
        return _assemble(requested, [], None, None, sources, temperature, pressure, 0, dropped)
    b_plus = np.asarray([b_map[element] for element in elements], dtype=float)
    stoichiometry = _stoichiometry(active, elements)
    c_values = _c_values(active, temperature, pressure)
    source_total = float(sum(sources.values()))
    q_matrix, d_values = _q_matrix(stoichiometry, b_plus, source_total)
    mole_fraction = np.full(len(active), 1.0 / len(active))
    total = _total_moles(stoichiometry, b_plus, mole_fraction)
    lambdas = initial_lambdas(stoichiometry, b_plus, c_values, mole_fraction, total)
    iterations = 0
    for iterations in range(1, _MAX_ITERATIONS + 1):
        exponent = stoichiometry.T @ lambdas - c_values
        # 式 (A.26) 使用未截断的指数。超出浮点范围时失败，不 clip。ADR 0025。
        if not np.all(np.isfinite(exponent)) or float(np.max(exponent)) > 700.0:
            raise EquilibriumError("式 (A.26) 的指数超出可计算范围")
        mole_fraction = np.exp(exponent)
        residual = _residuals(q_matrix, d_values, mole_fraction)
        amounts = total_amounts(stoichiometry, b_plus, mole_fraction)
        deviation = balance_deviation(stoichiometry, b_plus, amounts)
        if is_converged(float(np.max(np.abs(residual))), deviation):
            break
        jacobian = (q_matrix * mole_fraction) @ stoichiometry.T
        try:
            step = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError as exc:
            raise EquilibriumError("式 (A.29) 的线性方程奇异") from exc
        eta = _ETA0
        accepted = False
        # 当前点已不满足式 (A.31) 时，减半仍停在不可行点。解释性选择：ADR 0025。
        if not _exponent_feasible(stoichiometry, lambdas, c_values):
            lambdas = lambdas + eta * step
            continue
        for _halving in range(_MAX_HALVINGS + 1):
            trial = lambdas + eta * step
            if _exponent_feasible(stoichiometry, trial, c_values):
                lambdas = trial
                accepted = True
                break
            eta *= 0.5
        if not accepted:
            raise EquilibriumError("式 (A.31) 在减半 30 次后仍不成立")
    else:
        raise EquilibriumError(f"迭代 {_MAX_ITERATIONS} 次后未同时满足式 (A.26) 和式 (A.34)")
    amounts = total_amounts(stoichiometry, b_plus, mole_fraction)
    return _assemble(
        requested, active, mole_fraction, amounts, sources,
        temperature, pressure, iterations, dropped,
    )


def _assemble(requested, active, mole_fraction, amounts, sources, temperature, pressure, iterations, dropped):
    """待平衡量加上给定源汇量。源汇不因零元素移除而被删掉。ADR 0017 第 16 项。"""
    fractions = {name: 0.0 for name in requested}
    amount_map = {name: 0.0 for name in requested}
    total = float(sum(sources.values()))
    if mole_fraction is not None:
        for index, name in enumerate(active):
            fractions[name] = float(mole_fraction[index])
            amount_map[name] = float(amounts[index])
            total += float(amounts[index])
    if total > 0:
        for name, amount in sources.items():
            fractions[name] = amount / total
            amount_map[name] = amount
    return EquilibriumResult(fractions, amount_map, total, temperature, pressure, iterations, dropped)


def total_amounts(stoichiometry: np.ndarray, b_plus: np.ndarray, mole_fraction: np.ndarray) -> np.ndarray:
    """N_j = N·y_j，N 由式 (A.33)。"""
    return _total_moles(stoichiometry, b_plus, mole_fraction) * mole_fraction


def adiabatic_gas_equilibrium(
    feed_enthalpy_j: float,
    pressure_pa: float,
    element_amounts_mol: dict[str, float],
    species: list[str] | None = None,
) -> EquilibriumResult:
    """给定进料焓与压力，求气相绝热平衡。不含凝相。ADR 0017 第 9 项。 只对温度做一维求根，平衡计算与 ``gibbs_minimization`` 相同。"""
    names = list(GAS_SPECIES if species is None else species)
    target = float(feed_enthalpy_j)
    def residual(temperature: float) -> float:
        result = gibbs_minimization(temperature, pressure_pa, element_amounts_mol, names)
        return mixture_enthalpy_j(temperature, result.amounts_mol) - target
    grid = np.linspace(400.0, 2800.0, 25)
    values = []
    for temperature in grid:
        try:
            values.append((float(temperature), residual(float(temperature))))
        except (EquilibriumError, PropertyError):
            values.append((float(temperature), None))
    bracket = None
    for (left_t, left_v), (right_t, right_v) in zip(values, values[1:]):
        if left_v is None or right_v is None or left_v == 0:
            if left_v == 0:
                bracket = (left_t, left_t)
                break
            continue
        if left_v * right_v <= 0:
            bracket = (left_t, right_t)
            break
    if bracket is None:
        raise EquilibriumError("在 400–2800 K 内没有找到绝热温度")
    low, high = bracket
    if low == high:
        return gibbs_minimization(low, pressure_pa, element_amounts_mol, names)
    for _ in range(60):
        mid = 0.5 * (low + high)
        if residual(low) * residual(mid) <= 0:
            high = mid
        else:
            low = mid
    return gibbs_minimization(0.5 * (low + high), pressure_pa, element_amounts_mol, names)
