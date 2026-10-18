"""R1–R11 的计量系数。产物为正，反应物为负，对应一摩尔进度。"""
from __future__ import annotations

from hamel.kinetics.common import KineticsInputError, _finite

_FIXED = {
    "R2": {"C": -1.0, "H2O": -1.0, "CO": 1.0, "H2": 1.0},
    "R3": {"C": -1.0, "CO2": -1.0, "CO": 2.0},
    "R4": {"C": -1.0, "H2": -2.0, "CH4": 1.0},
    "R5": {"CO": -1.0, "O2": -0.5, "CO2": 1.0},
    "R6": {"H2": -1.0, "O2": -0.5, "H2O": 1.0},
    "R7": {"CH4": -1.0, "O2": -1.5, "CO": 1.0, "H2O": 2.0},
    "R8": {"CO": -1.0, "H2O": -1.0, "H2": 1.0, "CO2": 1.0},
    "R9": {"CH4": -1.0, "H2O": -1.0, "CO": 1.0, "H2": 3.0},
}


def reaction_stoichiometry(
    reaction: str,
    phi_c: float | None = None,
    tar: str | None = None,
) -> dict[str, float]:
    """返回 ν_i。R10、R11 只实现苯与萘。十六烷不在本卡范围。

    Source: Hamel Eq.5.8, S. 86
    Source: Hamel Eq.5.52, S. 95
    Source: Hamel Eq.5.55, S. 95
    """
    if reaction in _FIXED:
        return dict(_FIXED[reaction])
    if reaction == "R1":
        phi = _finite("phi_c", phi_c)
        if phi == 0.0:
            raise KineticsInputError("phi_c 不能为 0")
        return {
            "C": -1.0,
            "O2": -1.0 / phi,
            "CO": 2.0 - 2.0 / phi,
            "CO2": 2.0 / phi - 1.0,
        }
    if reaction == "R10" and tar == "C6H6":
        return {"C6H6": -1.0, "O2": -3.0, "CO": 6.0, "H2": 3.0}
    if reaction == "R10" and tar == "C10H8":
        return {"C10H8": -1.0, "O2": -5.0, "CO": 10.0, "H2": 4.0}
    if reaction == "R11" and tar == "C6H6":
        return {"C6H6": -1.0, "H2O": -4.5, "CO": 4.5, "H2": 4.5, "CH4": 1.5}
    if reaction == "R11" and tar == "C10H8":
        return {"C10H8": -1.0, "H2O": -23.0 / 3.0, "CO": 23.0 / 3.0, "H2": 7.0, "CH4": 7.0 / 3.0}
    raise KineticsInputError("反应编号或焦油代表物不在 R1–R11 的实现范围内")
