"""未知向量的编号、名字、打包与解包。"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.properties.data import GAS_SPECIES
from hamel.solver.frozen import AssemblyError, FrozenQuantities
# 工况 1：床区 8×0.15 m，自由空间 8×0.2875 m。二次风在自由空间第 2 格下沿。
CASE1_BED_CELLS = 8
CASE1_FREEBOARD_CELLS = 8
CASE1_SECONDARY_AIR_M = 1.4875
_PHASES_BED = ("suspension", "bubble")

@dataclass(frozen=True)
class UnknownSlot:
    """变量名表的一行。size_class 为 0 起算，原文粒级 k = size_class + 1。"""
    index: int
    name: str
    zone: str
    balance_index: int
    geometric_index: int
    kind: str
    phase: str
    species: str
    size_class: int

@dataclass
class ReactorState:
    """解包后的试探点。不存在的状态保持 0，不进入未知向量。"""
    bed_gas_mol_s: np.ndarray
    freeboard_gas_mol_s: np.ndarray
    bed_temperature_k: np.ndarray
    freeboard_temperature_k: np.ndarray
    bed_char_kg: np.ndarray
    freeboard_char_free_kg: np.ndarray
    freeboard_char_hover_kg: np.ndarray
    freeboard_ash_free_kg: np.ndarray
    freeboard_ash_hover_kg: np.ndarray

def balance_index(geometric_index: int, count: int) -> int:
    """区内自下而上的几何下标，换成衡算用的自上而下下标。下标 0 是最上层。
    Source: Hamel Eq.2.4, S. 16
    """
    _check_index(geometric_index, count)
    return count - 1 - geometric_index

def geometric_index(balance: int, count: int) -> int:
    """衡算下标换成区内自下而上的几何下标。
    Source: Hamel Eq.2.4, S. 16
    """
    _check_index(balance, count)
    return count - 1 - balance

def secondary_air_balance_index(n_freeboard: int = CASE1_FREEBOARD_CELLS) -> int:
    """自由空间第 2 格（自下而上，1.4875 m）的衡算下标。
    Source: Hamel Eq.2.4, S. 16
    """
    if n_freeboard < 2:
        raise AssemblyError("二次风所在的自由空间第 2 格要求至少两个自由空间格子")
    return balance_index(1, n_freeboard)

def unknown_slots(frozen: FrozenQuantities) -> tuple[UnknownSlot, ...]:
    """未知量顺序：自由空间自顶向下，然后床区自顶向下。
    每格先各相各组分出流，再温度，再各粒级的焦炭与（自由空间）灰。
    自由空间的灰在各粒级都是独立未知量，与 m_C 是否为零无关。解释性选择：ADR 0044
    Source: Hamel Eq.2.1, S. 15
    """
    slots: list[UnknownSlot] = []
    _append_zone(slots, frozen, "freeboard", frozen.n_freeboard)
    _append_zone(slots, frozen, "bed", frozen.n_bed)
    return tuple(slots)

def jacobian_layout(slots: tuple[UnknownSlot, ...], frozen: FrozenQuantities) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """节点链：床区各格自下而上，接自由空间各格自下而上；非相邻列按 ADR 0038 的 B 判定。求解手段：ADR 0046。"""
    nodes = np.zeros(len(slots), dtype=int)
    within = np.zeros(len(slots), dtype=int)
    separate = np.zeros(len(slots), dtype=bool)
    seen: dict[int, int] = {}
    for slot in slots:
        node = frozen.n_bed - 1 - slot.balance_index if slot.zone == "bed" else frozen.n_bed + frozen.n_freeboard - 1 - slot.balance_index
        nodes[slot.index] = node
        within[slot.index] = seen.get(node, 0)
        seen[node] = within[slot.index] + 1
        separate[slot.index] = slot.zone == "freeboard" and frozen.hover_closure == "B" and abs(node - frozen.n_bed + 1) > 1
    return nodes, within, separate


def pack_state(state: ReactorState, frozen: FrozenQuantities) -> np.ndarray:
    """按变量名表把状态收成未知向量。
    Source: Hamel Eq.2.4, S. 16
    """
    values = np.empty(len(unknown_slots(frozen)), dtype=float)
    for slot in unknown_slots(frozen):
        values[slot.index] = _read_slot(state, slot)
    if not np.isfinite(values).all():
        raise AssemblyError("未知向量含 NaN 或无穷大")
    return values

def unpack_state(values: np.ndarray, frozen: FrozenQuantities) -> ReactorState:
    """未知向量拆回各数组。不存在的状态为 0。
    Source: Hamel Eq.2.4, S. 16
    """
    slots = unknown_slots(frozen)
    vector = np.asarray(values, dtype=float)
    if vector.shape != (len(slots),) or not np.isfinite(vector).all():
        raise AssemblyError("未知向量的长度必须与变量名表一致，且全是有限数")
    state = _empty_state(frozen)
    for slot, value in zip(slots, vector):
        _write_slot(state, slot, float(value))
    return state

def _check_index(index: int, count: int) -> None:
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)):
        raise AssemblyError("格子下标必须是整数")
    if count < 1 or not 0 <= int(index) < count:
        raise AssemblyError("格子下标超出该区")

def _append_zone(slots: list[UnknownSlot], frozen: FrozenQuantities, zone: str, count: int) -> None:
    for balance in range(count):
        geometric = geometric_index(balance, count)
        phases = ("gas",) if zone == "freeboard" else _PHASES_BED
        for phase in phases:
            for species in GAS_SPECIES:
                _add(slots, zone, balance, geometric, "gas_outflow", phase, species, -1)
        _add(slots, zone, balance, geometric, "temperature", "", "", -1)
        for klass in range(frozen.n_class):
            if zone == "bed":
                _add(slots, zone, balance, geometric, "char", "bed", "", klass)
                continue
            if not frozen.free_absent[balance, klass]:
                _add(slots, zone, balance, geometric, "char", "free", "", klass)
                _add(slots, zone, balance, geometric, "ash", "free", "", klass)
            if frozen.hover[balance, klass]:
                _add(slots, zone, balance, geometric, "char", "hover", "", klass)
                _add(slots, zone, balance, geometric, "ash", "hover", "", klass)

def _add(slots, zone, balance, geometric, kind, phase, species, size_class) -> None:
    label = f"{zone}[b={balance},geo={geometric}].{kind}"
    if phase:
        label += f".{phase}"
    if species:
        label += f".{species}"
    if size_class >= 0:
        label += f".k{size_class + 1}"
    slots.append(UnknownSlot(len(slots), label, zone, balance, geometric, kind, phase, species, size_class))

def _empty_state(frozen: FrozenQuantities) -> ReactorState:
    n_species = len(GAS_SPECIES)
    bed, free, klass = frozen.n_bed, frozen.n_freeboard, frozen.n_class
    zeros = lambda shape: np.zeros(shape, dtype=float)
    return ReactorState(
        zeros((bed, 2, n_species)), zeros((free, n_species)), zeros(bed), zeros(free),
        zeros((bed, klass)), zeros((free, klass)), zeros((free, klass)),
        zeros((free, klass)), zeros((free, klass)),
    )

def _read_slot(state: ReactorState, slot: UnknownSlot) -> float:
    if slot.kind == "gas_outflow":
        if slot.zone == "bed":
            return float(state.bed_gas_mol_s[slot.balance_index, _PHASES_BED.index(slot.phase), GAS_SPECIES.index(slot.species)])
        return float(state.freeboard_gas_mol_s[slot.balance_index, GAS_SPECIES.index(slot.species)])
    if slot.kind == "temperature":
        temperatures = state.bed_temperature_k if slot.zone == "bed" else state.freeboard_temperature_k
        return float(temperatures[slot.balance_index])
    return float(_solid_array(state, slot)[slot.balance_index, slot.size_class])

def _write_slot(state: ReactorState, slot: UnknownSlot, value: float) -> None:
    if slot.kind == "gas_outflow":
        if slot.zone == "bed":
            state.bed_gas_mol_s[slot.balance_index, _PHASES_BED.index(slot.phase), GAS_SPECIES.index(slot.species)] = value
        else:
            state.freeboard_gas_mol_s[slot.balance_index, GAS_SPECIES.index(slot.species)] = value
        return
    if slot.kind == "temperature":
        temperatures = state.bed_temperature_k if slot.zone == "bed" else state.freeboard_temperature_k
        temperatures[slot.balance_index] = value
        return
    _solid_array(state, slot)[slot.balance_index, slot.size_class] = value

def _solid_array(state: ReactorState, slot: UnknownSlot) -> np.ndarray:
    if slot.zone == "bed":
        return state.bed_char_kg
    if slot.kind == "char" and slot.phase == "free":
        return state.freeboard_char_free_kg
    if slot.kind == "char":
        return state.freeboard_char_hover_kg
    if slot.phase == "free":
        return state.freeboard_ash_free_kg
    return state.freeboard_ash_hover_kg
