"""收敛门与外部平衡的数值判据。阈值来自 SolverConfig，不在这里另设。
ADR 0012。外部平衡用 hamel.cases.balance。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.cases.balance import BalanceGap


@dataclass(frozen=True)
class ResidualBlocks:
    """三块的下标和 ADR 0012 的参考流。参考流必须为正。"""

    gas: np.ndarray
    solid: np.ndarray
    energy: np.ndarray
    n_dot_ref: float
    m_dot_ref: float
    q_dot_ref: float


@dataclass(frozen=True)
class OuterMonitor:
    """外层 (a) 要比较的量。固体流含向上、向下、带出、惰性补料和排灰。"""

    temperature_k: np.ndarray
    solid_flow_kg_s: np.ndarray
    inventory_kg: np.ndarray
    gas_flow_mol_s: np.ndarray


def block_ratios(residual, blocks: ResidualBlocks) -> dict[str, float]:
    """各块 ‖r‖∞ 除以该块参考流。"""
    values = np.asarray(residual, dtype=float)
    return {
        "gas": _peak(values, blocks.gas, blocks.n_dot_ref),
        "solid": _peak(values, blocks.solid, blocks.m_dot_ref),
        "energy": _peak(values, blocks.energy, blocks.q_dot_ref),
    }


def inner_met(residual, blocks: ResidualBlocks, tolerance: float) -> bool:
    """ADR 0012 第 1 项。三块都要达到。"""
    return all(ratio <= tolerance for ratio in block_ratios(residual, blocks).values())


def outer_met(current: OuterMonitor, previous: OuterMonitor, *, m_dot_ref, n_dot_ref, m_tot, temperature_tolerance, relative_tolerance) -> bool:
    """ADR 0012 第 2(a) 项。分母是全炉尺度。新旧都为 0 时变化为 0。"""
    temperature = _abs_max(current.temperature_k, previous.temperature_k)
    solid = _abs_max(current.solid_flow_kg_s, previous.solid_flow_kg_s) / m_dot_ref
    inventory = _abs_max(current.inventory_kg, previous.inventory_kg) / m_tot
    gas = _abs_max(current.gas_flow_mol_s, previous.gas_flow_mol_s) / n_dot_ref
    return (
        temperature <= temperature_tolerance
        and solid <= relative_tolerance
        and inventory <= relative_tolerance
        and gas <= relative_tolerance
    )


def arrays_allowed(flows, inventories, temperatures, inert) -> bool:
    """气体出流与存量不小于 0，温度大于 0，惰性存量不小于 0。ADR 0017 第 7、14 项。"""
    groups = [np.asarray(item, dtype=float) for item in (*flows, *inventories, inert)]
    if any(np.any(item < 0.0) or not np.isfinite(item).all() for item in groups):
        return False
    temperature = np.asarray(temperatures, dtype=float)
    return bool(np.isfinite(temperature).all() and np.all(temperature > 0.0))


def balance_met(gap: BalanceGap, inlet_mol: tuple[float, float, float], heat_loss_w: float, q_dot_ref: float, element_tolerance: float, energy_tolerance: float) -> dict:
    """ADR 0012 第 3 项。热损失单独减去，不写进物流焓。进口元素不是正数则该项不通过。"""
    ratios = {}
    closed = True
    for name, delta, inlet in zip(("C", "H", "O"), (gap.C_mol_s, gap.H_mol_s, gap.O_mol_s), inlet_mol):
        if not np.isfinite(inlet) or inlet <= 0.0:
            ratios[name] = None
            closed = False
            continue
        ratios[name] = abs(float(delta)) / float(inlet)
        closed = closed and ratios[name] <= element_tolerance
    energy = abs(float(gap.energy_w) - float(heat_loss_w)) / float(q_dot_ref)
    ratios["energy"] = energy
    return {"ok": bool(closed and energy <= energy_tolerance and np.isfinite(energy)), "ratios": ratios}


def _peak(values, index, reference) -> float:
    chosen = np.asarray(index, dtype=int)
    if chosen.size == 0:
        return 0.0
    return float(np.max(np.abs(values[chosen])) / reference)


def _abs_max(current, previous) -> float:
    left = np.asarray(current, dtype=float).ravel()
    right = np.asarray(previous, dtype=float).ravel()
    if left.shape != right.shape:
        raise ValueError("外层监测量的形状不一致")
    if left.size == 0:
        return 0.0
    return float(np.max(np.abs(left - right)))
