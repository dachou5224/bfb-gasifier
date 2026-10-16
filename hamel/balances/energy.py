"""式 (2.7) 的能量残差，以及式 (2.8) 把一股物流收成能量流。"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from hamel.balances.inputs import (
    ResidualInputError,
    _inflow_from_above,
    _inflow_from_below,
    _require_cell_array,
    _require_float_array,
    _require_scalar,
    _require_shape,
)

@dataclass(frozen=True)
class EnergyResidual:
    """每个 cell 一条能量残差，W。下标 0 为最上层 cell。"""

    residual_w: np.ndarray
    top_upward_out_w: float
    bottom_downward_out_w: float

def energy_residual(
    h_in_ext_w,
    h_rec_w,
    h_up_w,
    h_down_w,
    h_out_sys_w,
    q_wall_w,
    q_heat_exchanger_w,
) -> EnergyResidual:
    """给定每个 cell 的能量流，返回式 (2.7) 的残差。

    各参数形状为 (cell,)，单位 W。``h_up`` 离开本 cell 进入上一层，``h_down`` 离开
    本 cell 进入下一层。缺失的邻居入流为零。最上层的上行出流和最下层的下行出流
    越出床层，作为边界流返回，不判定去向。解释性选择：ADR 0007。这里只用该 ADR 的
    编号方向，不把两端出流指定到自由空间或外排。

    壁面热流和内置换热面热流都按印刷式从残差减去。热流的物理正方向原文未另行规定，
    本函数不把减号改写成别的符号。

    Source: Hamel Eq.2.7, S. 19
    """
    upward = _require_cell_array("h_up_w", h_up_w, 1)
    shape = upward.shape
    downward = _require_shape("h_down_w", h_down_w, shape)
    external = _require_shape("h_in_ext_w", h_in_ext_w, shape)
    recycle = _require_shape("h_rec_w", h_rec_w, shape)
    leaving = _require_shape("h_out_sys_w", h_out_sys_w, shape)
    wall = _require_shape("q_wall_w", q_wall_w, shape)
    exchanger = _require_shape("q_heat_exchanger_w", q_heat_exchanger_w, shape)

    residual = (
        external
        + recycle
        + _inflow_from_below(upward)
        + _inflow_from_above(downward)
        - upward
        - downward
        - leaving
        - wall
        - exchanger
    )
    return EnergyResidual(
        residual_w=residual,
        top_upward_out_w=float(upward[0]),
        bottom_downward_out_w=float(downward[-1]),
    )

def stream_enthalpy_flow(
    gas_mole_flow_mol_s,
    gas_molar_enthalpy_j_per_mol,
    solid_mass_flow_kg_s,
    solid_specific_enthalpy_j_per_kg,
    enthalpy_flow_s_w,
) -> float:
    """把一股物流收成能量流，W。

    气体摩尔流和摩尔焓的形状为 (相, 组分)，单位 mol/s 与 J/mol。每种固体单独给
    一个粒级数组，质量流单位 kg/s，比焓单位 J/kg；不同固体的粒级数不必相同。
    单个一维数组表示只有一种固体。焓由调用方给出，含标准生成焓与否也由调用方决定。

    ``enthalpy_flow_s_w`` 是印刷式最后一项的乘积，单位 W，必填，没有默认值。
    下标 s 的含义尚未裁定（T8 歧义第 4 条），本函数不解释它是哪一种固体。

    Hamel 在式 (2.8) 写明这一加和采用 Wirsum（1998）对“与物流结合的能量流”的组成方式。
    本函数只做加和。

    Source: Hamel Eq.2.8, S. 19
    """
    gas_flow = _require_float_array("gas_mole_flow_mol_s", gas_mole_flow_mol_s, 2)
    gas_enthalpy = _require_shape(
        "gas_molar_enthalpy_j_per_mol", gas_molar_enthalpy_j_per_mol, gas_flow.shape
    )
    mass_flows = _solid_species(solid_mass_flow_kg_s, "solid_mass_flow_kg_s")
    enthalpies = _solid_species(solid_specific_enthalpy_j_per_kg, "solid_specific_enthalpy_j_per_kg")
    if len(mass_flows) != len(enthalpies):
        raise ResidualInputError("固体质量流与比焓的种类数不一致")
    extra = _require_scalar("enthalpy_flow_s_w", enthalpy_flow_s_w)

    total = float(np.sum(gas_flow * gas_enthalpy))
    for index, (mass, enthalpy) in enumerate(zip(mass_flows, enthalpies)):
        if mass.shape != enthalpy.shape:
            raise ResidualInputError(
                f"第 {index} 种固体的质量流形状是 {mass.shape}，比焓形状是 {enthalpy.shape}"
            )
        total += float(np.sum(mass * enthalpy))
    return total + extra

def _solid_species(value: object, name: str) -> list[np.ndarray]:
    """每种固体一个一维粒级数组。矩形二维数组不接受，避免把不同粒级数排齐。"""
    if isinstance(value, np.ndarray) and value.ndim == 1:
        items: Sequence[object] = (value,)
    elif isinstance(value, np.ndarray):
        raise ResidualInputError(f"{name} 应按固体种类分别给出一维粒级数组")
    elif isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ResidualInputError(f"{name} 应按固体种类分别给出一维粒级数组")
    else:
        items = value
    return [_require_float_array(f"{name}[{index}]", item, 1) for index, item in enumerate(items)]
