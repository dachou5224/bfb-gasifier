"""式 (2.4)、(2.5) 的气体摩尔流残差。各股流由调用方给出。"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.balances.inputs import (
    _copy_row,
    _inflow_from_below,
    _require_cell_array,
    _require_shape,
)

@dataclass(frozen=True)
class GasMoleResidual:
    """悬浮相与气泡相残差，mol/s。第 0 轴是 cell，下标 0 为最上层。"""

    suspension_mol_s: np.ndarray
    bubble_mol_s: np.ndarray
    top_suspension_out_mol_s: np.ndarray
    top_bubble_out_mol_s: np.ndarray

def gas_mole_residual(
    n_in_ext_suspension_mol_s,
    n_rec_suspension_mol_s,
    n_reaction_suspension_mol_s,
    n_out_suspension_mol_s,
    n_in_ext_bubble_mol_s,
    n_rec_bubble_mol_s,
    n_reaction_bubble_mol_s,
    n_out_bubble_mol_s,
    n_exchange_bd_mol_s,
) -> GasMoleResidual:
    """给定各股摩尔流，返回式 (2.4)、(2.5) 的残差。

    数组形状为 (cell, 气体组分)，单位 mol/s。``n_out_*`` 是本 cell 离开该相、
    向上进入上一层 cell 的流。最下层没有来自更下方的入流。最上层这两股出流离开
    床层，由 ``top_*`` 返回，供床层以上耦合；本函数不判定它们的去向。

    ``n_exchange_bd_mol_s`` 只按印刷式进入残差：悬浮相减、气泡相加。它的物理正方向
    尚未裁定（T8 歧义第 3 条）。这里不写浓度差，也不写交换速率。

    下标 0 为最上层，i+1 在下。这个编号方向与 ADR 0007 一致；该 ADR 不覆盖气相衡算的形式。

    Source: Hamel Eq.2.4, S. 16
    Source: Hamel Eq.2.5, S. 16
    """
    suspension_out = _require_cell_array("n_out_suspension_mol_s", n_out_suspension_mol_s, 2)
    shape = suspension_out.shape
    bubble_out = _require_shape("n_out_bubble_mol_s", n_out_bubble_mol_s, shape)
    external_suspension = _require_shape("n_in_ext_suspension_mol_s", n_in_ext_suspension_mol_s, shape)
    recycle_suspension = _require_shape("n_rec_suspension_mol_s", n_rec_suspension_mol_s, shape)
    reaction_suspension = _require_shape(
        "n_reaction_suspension_mol_s", n_reaction_suspension_mol_s, shape
    )
    external_bubble = _require_shape("n_in_ext_bubble_mol_s", n_in_ext_bubble_mol_s, shape)
    recycle_bubble = _require_shape("n_rec_bubble_mol_s", n_rec_bubble_mol_s, shape)
    reaction_bubble = _require_shape("n_reaction_bubble_mol_s", n_reaction_bubble_mol_s, shape)
    exchange = _require_shape("n_exchange_bd_mol_s", n_exchange_bd_mol_s, shape)

    from_below_suspension = _inflow_from_below(suspension_out)
    from_below_bubble = _inflow_from_below(bubble_out)
    suspension = (
        external_suspension
        + recycle_suspension
        + from_below_suspension
        + reaction_suspension
        - suspension_out
        - exchange
    )
    bubble = (
        external_bubble
        + recycle_bubble
        + from_below_bubble
        + reaction_bubble
        - bubble_out
        + exchange
    )
    return GasMoleResidual(
        suspension_mol_s=suspension,
        bubble_mol_s=bubble,
        top_suspension_out_mol_s=_copy_row(suspension_out, 0),
        top_bubble_out_mol_s=_copy_row(bubble_out, 0),
    )
