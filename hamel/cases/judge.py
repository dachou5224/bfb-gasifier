"""验收判定。四种状态：通过、不通过、边界待定、原文无。前提失败则直接不通过。"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

from hamel.cases.targets import AcceptanceTarget, registered_signature

PASS = "通过"
FAIL = "不通过"
BORDERLINE = "边界待定"
ABSENT = "原文无"
PENDING = "待定"

class JudgementError(Exception):
    """目标集合不完整、工况对不上，或测量的单位与基准和目标不一致。"""

@dataclass(frozen=True)
class Measurement:
    """一个验收量的模型结果。单位和气体基准必须与目标一致，不做混比换算。"""

    value: float
    unit: str
    basis: str | None

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
    measured: Mapping[str, Measurement],
    targets: Sequence[AcceptanceTarget],
    premises: Premises,
) -> CaseJudgement:
    """判定一个工况。缺测不补数。目标集合不完整或单位、基准不一致时报错，不给通过。"""
    _require_targets(case_id, targets)
    _require_measurements(measured, targets)
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

def _require_targets(case_id: str, targets: Sequence[AcceptanceTarget]) -> None:
    names = [item.quantity for item in targets]
    if len(names) != len(set(names)):
        raise JudgementError("验收量有重复，不能给出整案结论")
    if any(item.case_id != case_id for item in targets):
        raise JudgementError(f"目标不属于 {case_id}")
    got = set()
    for item in targets:
        unit = item.unit if item.present else None
        basis = item.basis if item.present else None
        got.add((item.quantity, item.present, unit, basis))
    if got != registered_signature(case_id):
        raise JudgementError(f"{case_id} 的验收量与章程登记不一致")

def _require_measurements(measured: Mapping[str, Measurement], targets: Sequence[AcceptanceTarget]) -> None:
    for target in targets:
        if not target.present or target.quantity not in measured:
            continue
        item = measured[target.quantity]
        if not isinstance(item, Measurement):
            raise JudgementError(f"{target.quantity} 必须带单位和基准")
        if item.unit != target.unit or item.basis != target.basis:
            raise JudgementError(f"{target.quantity} 的单位或基准与目标不一致，不能混比")

def _judge_one(target: AcceptanceTarget, measured: Mapping[str, Measurement]) -> QuantityJudgement:
    if not target.present:
        return QuantityJudgement(target.quantity, ABSENT, None, None)
    item = measured.get(target.quantity)
    if item is None:
        return QuantityJudgement(target.quantity, FAIL, None, None)
    value = float(item.value)
    distance = distance_to_interval(value, target.lower, target.upper)
    outside = value < target.lower or value > target.upper
    if distance < target.est_error:
        status = BORDERLINE
    elif outside:
        status = FAIL
    else:
        status = PASS
    return QuantityJudgement(target.quantity, status, value, distance)

def _premise_failure(target: AcceptanceTarget, measured: Mapping[str, Measurement]) -> QuantityJudgement:
    """前提失败时不计算到边界的距离，因此不会得到“边界待定”。"""
    if not target.present:
        return QuantityJudgement(target.quantity, ABSENT, None, None)
    item = measured.get(target.quantity)
    value = None if item is None else item.value
    return QuantityJudgement(target.quantity, FAIL, value, None)

def _measurements_finite_nonnegative(
    measured: Mapping[str, Measurement],
    targets: Sequence[AcceptanceTarget],
) -> bool:
    for target in targets:
        if not target.present or target.quantity not in measured:
            continue
        value = measured[target.quantity].value
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
