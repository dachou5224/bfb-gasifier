"""式 (2.6) 的固体质量流残差。只做流量形式，不闭合存量。"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.balances.inputs import (
    _copy_row,
    _inflow_from_above,
    _inflow_from_below,
    _require_cell_array,
    _require_shape,
)

@dataclass(frozen=True)
class SolidMassFlowResidual:
    """一种固体的粒级质量流残差，kg/s。第 0 轴是 cell，第 1 轴是粒级。"""

    residual_kg_s: np.ndarray
    top_upward_out_kg_s: np.ndarray
    bottom_downward_out_kg_s: np.ndarray

def solid_mass_flow_residual(
    m_in_ext_kg_s,
    m_rec_kg_s,
    m_out_sys_kg_s,
    m_reaction_kg_s,
    m_from_class_kg_s,
    m_to_class_kg_s,
    m_up_kg_s,
    m_down_kg_s,
) -> SolidMassFlowResidual:
    """给定一种固体各粒级的质量流，返回式 (2.6) 的残差。

    形状为 (cell, 粒级)，单位 kg/s。多种固体的粒级数可以不同，按固体种类分别调用。
    第 1 轴只放粒级 k，不把不同固体排进同一条轴。``m_from_class`` 是自其它粒级进入
    本粒级的流，对应印刷式的左箭头；``m_to_class`` 是离开本粒级的流，对应右箭头。

    相邻入流由出流数组按编号方向取得。来自下方的上行入流是下方 cell 同一粒级的
    ``m_up``。解释性选择：ADR 0008。来自上方的下行入流是上方 cell 的 ``m_down``。
    最下层没有来自更下方的入流，最上层没有来自更上方的入流。

    最上层的上行出流和最下层的下行出流不进入任何 cell，由边界流字段原样返回。
    解释性选择：ADR 0007。这里只采用该 ADR 的编号方向和“两端出流不自行指定去向”；
    不实现存量闭合，本函数不接收存量，也不接收系数 K。

    Source: Hamel Eq.2.6, S. 17
    """
    upward = _require_cell_array("m_up_kg_s", m_up_kg_s, 2)
    shape = upward.shape
    downward = _require_shape("m_down_kg_s", m_down_kg_s, shape)
    external = _require_shape("m_in_ext_kg_s", m_in_ext_kg_s, shape)
    recycle = _require_shape("m_rec_kg_s", m_rec_kg_s, shape)
    leaving = _require_shape("m_out_sys_kg_s", m_out_sys_kg_s, shape)
    reaction = _require_shape("m_reaction_kg_s", m_reaction_kg_s, shape)
    from_class = _require_shape("m_from_class_kg_s", m_from_class_kg_s, shape)
    to_class = _require_shape("m_to_class_kg_s", m_to_class_kg_s, shape)

    from_below = _inflow_from_below(upward)
    from_above = _inflow_from_above(downward)
    residual = (
        external
        + recycle
        - leaving
        - upward
        - downward
        + from_below
        + from_above
        + reaction
        + from_class
        - to_class
    )
    return SolidMassFlowResidual(
        residual_kg_s=residual,
        top_upward_out_kg_s=_copy_row(upward, 0),
        bottom_downward_out_kg_s=_copy_row(downward, -1),
    )
