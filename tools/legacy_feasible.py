"""ADR 0017 第 7 项的可行条件：供历史诊断工具复现旧记录。

正式路径已由 ADR 0054 第 2 项的三层试探规则取代（``hamel.solver.newton.TrialRule``），
本模块只服务 T26c–T26l 的运行记录，不参与求解，也不在 ``hamel/`` 中使用。
"""
from __future__ import annotations

import numpy as np

from hamel.solver.vector import unpack_state


def feasible(values, frozen) -> bool:
    """气体出流与存量不小于 0，温度大于 0，惰性存量不小于 0。ADR 0017 第 7、14 项。"""
    state = unpack_state(values, frozen)
    arrays = [state.bed_gas_mol_s, state.freeboard_gas_mol_s, state.bed_char_kg,
              state.freeboard_char_free_kg, state.freeboard_char_hover_kg,
              state.freeboard_ash_free_kg, state.freeboard_ash_hover_kg,
              _inert(state, frozen)]
    if any(np.any(np.asarray(item, dtype=float) < 0.0) for item in arrays):
        return False
    temperature = np.concatenate([state.bed_temperature_k, state.freeboard_temperature_k])
    return bool(np.isfinite(temperature).all() and np.all(temperature > 0.0))


def clause(values, frozen) -> str:
    """与 ``feasible`` 同一组条件，按 ADR 0017 第 7 项的顺序报告第一条。"""
    state = unpack_state(values, frozen)
    groups = (
        ("Ṅ < 0", state.bed_gas_mol_s), ("Ṅ < 0", state.freeboard_gas_mol_s),
        ("m < 0", state.bed_char_kg), ("m < 0", state.freeboard_char_free_kg),
        ("m < 0", state.freeboard_char_hover_kg), ("m < 0", state.freeboard_ash_free_kg),
        ("m < 0", state.freeboard_ash_hover_kg), ("m_inert < 0", _inert(state, frozen)),
    )
    for label, array in groups:
        flat = np.asarray(array, dtype=float).ravel()
        if not np.isfinite(flat).all():
            return f"{label} 的量非有限"
        if np.any(flat < 0.0):
            return f"{label}（下标 {int(np.flatnonzero(flat < 0.0)[0])}）"
    temperature = np.concatenate([state.bed_temperature_k, state.freeboard_temperature_k])
    if not np.isfinite(temperature).all():
        return "T 非有限"
    if np.any(temperature <= 0.0):
        return f"T ≤ 0（下标 {int(np.flatnonzero(temperature <= 0.0)[0])}）"
    return "feasible 不满足"


def _inert(state, frozen):
    return frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
