"""从 fig_*.json 读取 Hamel 模拟值。验收量以章程 §2.2 的登记为准。"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

_GAS = "major_gas"
_CH4 = "ch4"
_TEMP = "temperature"
_XC = "carbon_conversion"

# CHARTER v0.8 §2.2。VTT 加压 1 的 7 个量都在 fig_7_10.json。
_VTT1 = (
    ("T_FRB_K", _TEMP),
    ("X_C", _XC),
    ("CO", _GAS),
    ("CO2", _GAS),
    ("H2", _GAS),
    ("CH4", _CH4),
    ("H2O", _GAS),
)
# CHARTER §2.2。WSV400 为后备工况：干基气体，原文无 H2O 与碳转化率。
_WSV = (
    ("T_FRB_K", _TEMP),
    ("CO", _GAS),
    ("CO2", _GAS),
    ("H2", _GAS),
    ("CH4", _CH4),
)
_WSV_ABSENT = ("H2O", "X_C")
_WSV_LETTERS = tuple("ABCDEF")
_GAS_NAMES = {"CO", "CO2", "H2", "CH4", "H2O"}

@dataclass(frozen=True)
class AcceptanceTarget:
    """一个验收量。``present`` 为假时表示原文无，没有区间。"""

    case_id: str
    quantity: str
    present: bool
    sim: float | None
    lower: float | None
    upper: float | None
    est_error: float | None
    unit: str | None
    basis: str | None
    source: str

class TargetError(Exception):
    """目标文件缺验收量、缺读数误差，或该工况尚未登记。"""

def half_width(kind: str, sim: float) -> float:
    """章程 §2.2 的半宽。气体取绝对容差与相对容差中较宽者。"""
    if kind == _TEMP:
        return 20.0
    if kind == _GAS:
        return max(0.015, 0.10 * abs(sim))
    if kind == _CH4:
        return max(0.005, 0.20 * abs(sim))
    if kind == _XC:
        return 3.0
    raise TargetError(f"没有容差种类 {kind}")

def _meta(quantity: str, present: bool, gas_basis: str | None) -> tuple:
    """章程登记的（量，是否有原文，单位，气体基准）。"""
    if not present:
        return (quantity, False, None, None)
    if quantity in {"T_FRB_K", "T_exit_K"}:
        return (quantity, True, "K", None)
    if quantity == "X_C":
        return (quantity, True, "%", None)
    if quantity in _GAS_NAMES:
        return (quantity, True, "mol/mol", gas_basis)
    raise TargetError(f"没有登记 {quantity} 的单位")

def registered_signature(case_id: str) -> set[tuple]:
    """已登记工况的验收量集合。未登记则报错，不把图上的点自动算进去。"""
    if case_id == "CASE_VTT_PRESSURIZED_1":
        return {_meta(name, True, "wet") for name, _ in _VTT1}
    letter = _wsv_letter(case_id)
    if letter is None:
        raise TargetError(f"{case_id} 的验收量尚未在章程 §2.2 登记")
    rows = {_meta(name, True, "dry") for name, _ in _WSV}
    rows.update(_meta(name, False, None) for name in _WSV_ABSENT)
    if letter == "A":
        rows.add(_meta("T_exit_K", True, None))
    return rows

def acceptance_targets(case_id: str, digitized_dir: Path | str) -> tuple[AcceptanceTarget, ...]:
    """生成已登记工况的目标。未登记的工况报错，不自行选定量。"""
    folder = Path(digitized_dir)
    if case_id == "CASE_VTT_PRESSURIZED_1":
        return _from_parity(case_id, folder / "fig_7_10.json", _VTT1, (), "wet")
    letter = _wsv_letter(case_id)
    if letter is None:
        raise TargetError(f"{case_id} 的验收量尚未在章程 §2.2 登记")
    targets = list(_from_parity(case_id, folder / "fig_7_17.json", _WSV, _WSV_ABSENT, "dry"))
    if letter == "A":
        targets.append(_exit_temperature(case_id, folder / "fig_7_16.json"))
    return tuple(targets)

def _wsv_letter(case_id: str) -> str | None:
    prefix = "CASE_WSV400_"
    if case_id.startswith(prefix) and case_id[len(prefix):] in _WSV_LETTERS:
        return case_id[len(prefix):]
    return None

def _from_parity(
    case_id: str,
    path: Path,
    quantities: tuple[tuple[str, str], ...],
    absent: tuple[str, ...],
    gas_basis: str,
) -> tuple[AcceptanceTarget, ...]:
    doc = _load(path)
    case = _parity_case(doc, case_id)
    if doc.get("basis_gas") != gas_basis:
        raise TargetError(f"{path.name} 的气体基准不是 {gas_basis}")
    units = doc.get("units", {})
    found = []
    for quantity, kind in quantities:
        pair = case.get(quantity)
        if not isinstance(pair, list) or len(pair) != 2:
            raise TargetError(f"{path.name} 的 {case_id} 没有 {quantity} 的模拟值")
        if quantity not in doc.get("est_error", {}):
            raise TargetError(f"{path.name} 没有 {quantity} 的 est_error")
        _, _, unit, basis = _meta(quantity, True, gas_basis)
        if units.get(quantity) != unit:
            raise TargetError(f"{path.name} 的 {quantity} 单位不是 {unit}")
        sim = float(pair[1])
        found.append(_point(
            case_id, quantity, path.name, present=True, sim=sim, kind=kind,
            error=doc["est_error"][quantity], unit=unit, basis=basis,
        ))
    for quantity in absent:
        found.append(_point(case_id, quantity, path.name, present=False))
    return tuple(found)

def _point(case_id, quantity, source, *, present, sim=None, kind=None, error=None, unit=None, basis=None):
    if not present:
        return AcceptanceTarget(case_id, quantity, False, None, None, None, None, None, None, source)
    width = half_width(kind, sim)
    return AcceptanceTarget(
        case_id, quantity, True, sim, sim - width, sim + width, float(error), unit, basis, source,
    )

def _exit_temperature(case_id: str, path: Path) -> AcceptanceTarget:
    """Bild 7.16 出口端温度，只补充核对 WSV400 A，不替代 T̄_FRB。"""
    doc = _load(path)
    if doc.get("case_id") != case_id:
        raise TargetError(f"{path.name} 不是 {case_id}")
    if "T" not in doc.get("est_error", {}) or doc.get("units", {}).get("T") != "K":
        raise TargetError(f"{path.name} 没有以 K 给出出口温度及其 est_error")
    sim = None
    for series in doc.get("series", []):
        if series.get("kind") == "sim" and series.get("variable") == "T":
            exits = [point for point in series.get("points", []) if point[0] == 1 or point[0] == 1.0]
            if len(exits) != 1:
                raise TargetError(f"{path.name} 出口端温度不是唯一的 ξ=1 点")
            sim = float(exits[0][1])
    if sim is None:
        raise TargetError(f"{path.name} 没有模拟出口温度")
    return _point(
        case_id, "T_exit_K", path.name, present=True, sim=sim, kind=_TEMP,
        error=doc["est_error"]["T"], unit="K", basis=None,
    )

def _parity_case(doc: dict, case_id: str) -> dict:
    for case in doc.get("cases", []):
        if case.get("case_id") == case_id:
            return case
    raise TargetError(f"图中没有 {case_id}")

def _load(path: Path) -> dict:
    if not path.is_file():
        raise TargetError(f"找不到 {path}")
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)
