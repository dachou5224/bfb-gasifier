"""外部元素平衡与总能量缺口。

只使用进出口物流，不使用模型内部的 cell 残差。
元素缺口 = 进口 − 出口，单位 mol/s；能量缺口单位 W。
负的碳缺口表示出口碳多于进口（床顶凭空产生碳时即为负）。
组分到元素的对应是分子式记账，不是反应速率方程。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence, Union

# ADR 0004 相容性检查写明的原子量，kg/kmol。这里换成 kg/mol。
_MOLAR_MASS_KG_PER_MOL = {
    "C": 12.011 / 1000.0,
    "H": 1.008 / 1000.0,
    "O": 15.999 / 1000.0,
}
_ELEMENTS = ("C", "H", "O")
# 本检查器只报告 C/H/O。氮和硫可以出现在物流里，但不计入缺口。
_IGNORED_ELEMENTS = {"N", "S"}
_SPECIES_ATOMS = {
    "CO": {"C": 1, "O": 1},
    "CO2": {"C": 1, "O": 2},
    "H2": {"H": 2},
    "H2O": {"H": 2, "O": 1},
    "CH4": {"C": 1, "H": 4},
    "O2": {"O": 2},
    "N2": {},
}

class BalanceInputError(Exception):
    """物流缺元素、组分未知，或流率为负。不截断、不补零。"""

@dataclass(frozen=True)
class GasStream:
    """组分摩尔流，mol/s。``enthalpy_flow_w`` 须显式给出，单位 W。"""

    mole_flow_mol_s: Mapping[str, float]
    enthalpy_flow_w: float

@dataclass(frozen=True)
class SolidStream:
    """固体质量流。元素质量分数单位 kg/kg，必须含 C、H、O。"""

    mass_flow_kg_s: float
    element_mass_fraction: Mapping[str, float]
    enthalpy_flow_w: float

@dataclass(frozen=True)
class ElementStream:
    """已经是元素摩尔流的物流，mol/s。T1 探针给出的是这一层数字。"""

    mole_flow_mol_s: Mapping[str, float]
    enthalpy_flow_w: float

Stream = Union[GasStream, SolidStream, ElementStream]

@dataclass(frozen=True)
class BalanceGap:
    """进口减出口。``energy_w`` 为总能量缺口。"""

    C_mol_s: float
    H_mol_s: float
    O_mol_s: float
    energy_w: float

def external_balance(inlets: Sequence[Stream], outlets: Sequence[Stream]) -> BalanceGap:
    """计算 C/H/O 缺口（mol/s）与总能量缺口（W）。"""
    inlet = _sum_streams(inlets)
    outlet = _sum_streams(outlets)
    return BalanceGap(
        C_mol_s=inlet[0] - outlet[0],
        H_mol_s=inlet[1] - outlet[1],
        O_mol_s=inlet[2] - outlet[2],
        energy_w=inlet[3] - outlet[3],
    )

def _sum_streams(streams: Sequence[Stream]) -> tuple[float, float, float, float]:
    total = [0.0, 0.0, 0.0, 0.0]
    for stream in streams:
        atoms, energy = _stream_atoms(stream)
        for index, element in enumerate(_ELEMENTS):
            total[index] += atoms[element]
        total[3] += energy
    return total[0], total[1], total[2], total[3]

def _stream_atoms(stream: Stream) -> tuple[dict[str, float], float]:
    if isinstance(stream, GasStream):
        atoms = _from_species(stream.mole_flow_mol_s)
        energy = _energy(stream.enthalpy_flow_w)
    elif isinstance(stream, SolidStream):
        atoms = _from_solid(stream.mass_flow_kg_s, stream.element_mass_fraction)
        energy = _energy(stream.enthalpy_flow_w)
    elif isinstance(stream, ElementStream):
        atoms = _from_elements(stream.mole_flow_mol_s)
        energy = _energy(stream.enthalpy_flow_w)
    else:
        raise BalanceInputError(f"无法识别的物流类型：{type(stream).__name__}")
    return atoms, energy

def _from_species(mole_flow: Mapping[str, float]) -> dict[str, float]:
    atoms = {element: 0.0 for element in _ELEMENTS}
    if not mole_flow:
        raise BalanceInputError("气体物流没有组分")
    for species, flow in mole_flow.items():
        _require_nonnegative(species, flow)
        if species not in _SPECIES_ATOMS:
            raise BalanceInputError(f"未知组分 {species}，不能忽略")
        for element, count in _SPECIES_ATOMS[species].items():
            if element in atoms:
                atoms[element] += float(flow) * count
    return atoms

def _from_solid(mass_flow_kg_s: float, fractions: Mapping[str, float]) -> dict[str, float]:
    _require_nonnegative("mass_flow_kg_s", mass_flow_kg_s)
    missing = [element for element in _ELEMENTS if element not in fractions]
    if missing:
        raise BalanceInputError("固体元素组成缺少 " + "、".join(missing))
    atoms = {}
    for element in _ELEMENTS:
        fraction = fractions[element]
        _require_nonnegative(element, fraction)
        atoms[element] = float(mass_flow_kg_s) * float(fraction) / _MOLAR_MASS_KG_PER_MOL[element]
    return atoms

def _from_elements(mole_flow: Mapping[str, float]) -> dict[str, float]:
    unknown = [name for name in mole_flow if name not in _ELEMENTS and name not in _IGNORED_ELEMENTS]
    if unknown:
        raise BalanceInputError("无法计入的元素 " + "、".join(unknown))
    missing = [element for element in _ELEMENTS if element not in mole_flow]
    if missing:
        raise BalanceInputError("元素物流缺少 " + "、".join(missing))
    atoms = {}
    for element in _ELEMENTS:
        _require_nonnegative(element, mole_flow[element])
        atoms[element] = float(mole_flow[element])
    return atoms

def _energy(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise BalanceInputError("焓流必须显式给出，单位 W")
    number = float(value)
    if number != number or number in (float("inf"), float("-inf")):
        raise BalanceInputError("焓流不是有限值")
    return number

def _require_nonnegative(name: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise BalanceInputError(f"{name} 不是数值")
    number = float(value)
    if number != number or number in (float("inf"), float("-inf")):
        raise BalanceInputError(f"{name} 不是有限值")
    if number < 0:
        raise BalanceInputError(f"{name} 为负，不截断")
