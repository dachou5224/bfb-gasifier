"""验收判定。每个验收量四种状态之一：通过、不通过、边界待定、原文无。

边界待定：结果与合格区间某一端的距离小于该量的 ``est_error``，区间内外都算。
距离等于读数误差时不算边界待定（章程写的是“小于”）。
任一验收量不通过，工况即为不通过；没有不通过但有边界待定，工况为待定。
§2.1 的前提不满足时，工况直接不通过，不再做边界判定。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

from hamel.cases.targets import AcceptanceTarget

PASS = "通过"
FAIL = "不通过"
BORDERLINE = "边界待定"
ABSENT = "原文无"
PENDING = "待定"

@dataclass(frozen=True)
class Premises:
    """章程 §2.1。三项都要显式传入，不设默认。"""

    residuals_met: bool
    balance_closed: bool
    state_finite_nonnegative: bool

@dataclass(frozen=True)
class QuantityJudgement:
    quantity: str
    status: str
    value: float | None
    distance_to_boundary: float | None

@dataclass(frozen=True)
class CaseJudgement:
    case_id: str
    status: str
    quantities: tuple[QuantityJudgement, ...]
    premises_met: bool
    balance_closed: bool

def judge(
    case_id: str,
    measured: Mapping[str, float],
    targets: Sequence[AcceptanceTarget],
    premises: Premises,
) -> CaseJudgement:
    """判定一个工况。``measured`` 的键与验收量名称一致，缺测不补数。"""
    state_ok = premises.state_finite_nonnegative and _measurements_finite_nonnegative(measured, targets)
    premises_met = bool(premises.residuals_met and premises.balance_closed and state_ok)
    if not premises_met:
        quantities = tuple(_premise_failure(target, measured) for target in targets)
        return CaseJudgement(
            case_id=case_id,
            status=FAIL,
            quantities=quantities,
            premises_met=False,
            balance_closed=bool(premises.balance_closed),
        )
    quantities = tuple(_judge_one(target, measured) for target in targets)
    return CaseJudgement(
        case_id=case_id,
        status=_rollup(quantity.status for quantity in quantities),
        quantities=quantities,
        premises_met=True,
        balance_closed=True,
    )

def distance_to_interval(value: float, lower: float, upper: float) -> float:
    if lower <= value <= upper:
        return min(value - lower, upper - value)
    if value < lower:
        return lower - value
    return value - upper

def _judge_one(target: AcceptanceTarget, measured: Mapping[str, float]) -> QuantityJudgement:
    if not target.present:
        return QuantityJudgement(target.quantity, ABSENT, measured.get(target.quantity), None)
    if target.quantity not in measured:
        return QuantityJudgement(target.quantity, FAIL, None, None)
    value = float(measured[target.quantity])
    distance = distance_to_interval(value, target.lower, target.upper)
    outside = value < target.lower or value > target.upper
    if distance < target.est_error:
        status = BORDERLINE
    elif outside:
        status = FAIL
    else:
        status = PASS
    return QuantityJudgement(target.quantity, status, value, distance)

def _premise_failure(target: AcceptanceTarget, measured: Mapping[str, float]) -> QuantityJudgement:
    """前提失败时不计算到边界的距离，因此不会得到“边界待定”。"""
    value = measured.get(target.quantity)
    if not target.present:
        return QuantityJudgement(target.quantity, ABSENT, value, None)
    return QuantityJudgement(target.quantity, FAIL, value, None)

def _measurements_finite_nonnegative(
    measured: Mapping[str, float],
    targets: Sequence[AcceptanceTarget],
) -> bool:
    for target in targets:
        if not target.present or target.quantity not in measured:
            continue
        value = measured[target.quantity]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return False
        if not math.isfinite(float(value)) or float(value) < 0:
            return False
    return True

def _rollup(statuses: Sequence[str] | object) -> str:
    judged = [status for status in statuses if status != ABSENT]
    if any(status == FAIL for status in judged):
        return FAIL
    if any(status == BORDERLINE for status in judged):
        return PENDING
    if judged and all(status == PASS for status in judged):
        return PASS
    return FAIL
