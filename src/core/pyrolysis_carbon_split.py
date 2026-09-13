"""热解碳分率：Hegermann / 固定 tar。

与 ``calc_drying_pyrolysis_sources`` 同口径。``koks=0`` 且 tar 固定时，
偏碳气相池可能对 {CO,CO2,CH4,H2,H2O} 的 ``An=b`` 不可行；Hamel §4.2
先按 Eq.4.21–4.22 分出 Koks/Teer，禁止事后把剩余 C 划回 Koks。

Source: Hamel (1999) Eq.4.21–4.22, Tabelle 4.3; Hegermann (1990)
"""

from __future__ import annotations


def resolve_pyrolysis_carbon_fractions(
    *,
    fuel_type: str,
    yield_mode: str,
    P_Pa: float,
    T_K: float,
    pyrolysis_tar_carbon_frac: float = 0.2,
) -> tuple[float, float]:
    """返回 ``(f_koks, f_tar)``，均为相对输入 ``nC`` 的碳分率 [-]。

    ``coal`` + ``hegermann``：Eq.4.21 / 4.22。
    其它：``f_koks=0``，``f_tar`` 用调用方固定值。
    """
    mode = str(yield_mode).strip().lower()
    fuel = str(fuel_type).strip().lower()
    if mode == "hegermann" and fuel == "coal":
        from src.thermal.hegermann_yield import xc_koks_hegermann, xc_teer_hegermann

        return (
            float(xc_koks_hegermann(float(P_Pa), float(T_K))),
            float(xc_teer_hegermann(float(P_Pa), float(T_K))),
        )
    if mode not in {"hegermann", "fixed"}:
        raise ValueError(
            f"Unknown pyrolysis_tar_yield_mode={yield_mode!r}; expected 'hegermann' or 'fixed'"
        )
    return 0.0, float(pyrolysis_tar_carbon_frac)
