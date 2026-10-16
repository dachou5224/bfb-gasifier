"""残差输入检查。负值不截断、不报错，按原值进入残差。"""
from __future__ import annotations

import numpy as np

class ResidualInputError(ValueError):
    """形状不一致，或数值里出现 NaN。"""

def _require_float_array(name: str, value: object, ndim: int) -> np.ndarray:
    if isinstance(value, bool):
        raise ResidualInputError(f"{name} 不是数组")
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ResidualInputError(f"{name} 不能转换成浮点数组") from exc
    if array.ndim != ndim:
        raise ResidualInputError(f"{name} 应为 {ndim} 维数组，实际为 {array.ndim} 维")
    if np.isnan(array).any():
        raise ResidualInputError(f"{name} 含 NaN")
    return array

def _require_cell_array(name: str, value: object, ndim: int) -> np.ndarray:
    """第 0 轴是 cell，下标 0 为最上层。至少要有一个 cell。"""
    array = _require_float_array(name, value, ndim)
    if array.shape[0] < 1:
        raise ResidualInputError(f"{name} 至少要有一个 cell")
    return array

def _require_shape(name: str, value: object, shape: tuple[int, ...]) -> np.ndarray:
    array = _require_float_array(name, value, len(shape))
    if array.shape != shape:
        raise ResidualInputError(f"{name} 的形状是 {array.shape}，应为 {shape}")
    return array

def _require_scalar(name: str, value: object) -> float:
    if isinstance(value, bool) or isinstance(value, np.ndarray):
        raise ResidualInputError(f"{name} 必须是标量，单位 W")
    if not isinstance(value, (int, float, np.integer, np.floating)):
        raise ResidualInputError(f"{name} 必须是标量，单位 W")
    number = float(value)
    if number != number:
        raise ResidualInputError(f"{name} 含 NaN")
    return number

def _inflow_from_below(outflow: np.ndarray) -> np.ndarray:
    """cell i 接收下方 cell（下标 i+1）的出流。最下层没有这一项。"""
    inflow = np.zeros_like(outflow)
    if outflow.shape[0] > 1:
        inflow[:-1] = outflow[1:]
    return inflow

def _inflow_from_above(outflow: np.ndarray) -> np.ndarray:
    """cell i 接收上方 cell（下标 i−1）的出流。最上层没有这一项。"""
    inflow = np.zeros_like(outflow)
    if outflow.shape[0] > 1:
        inflow[1:] = outflow[:-1]
    return inflow

def _copy_row(array: np.ndarray, index: int) -> np.ndarray:
    return np.array(array[index], dtype=float, copy=True)
