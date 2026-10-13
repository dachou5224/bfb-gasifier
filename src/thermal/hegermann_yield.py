"""Hegermann (1990) 热解碳分配：Xc,Koks / Xc,Teer（Hamel Eq.4.21–4.22）。

仅适用于 Rheinische Braunkohle（本仓库 ``fuel_type="coal"``）。
压力单位 Pa，温度单位 K；参数见 Hamel Tabelle 4.3。

Source: Hamel (1999) §4.2.1.2, Eq.4.21–4.22, Tabelle 4.3;
        Hegermann (1990) 压力流化床元素衡算。
"""

from __future__ import annotations

import numpy as np

# Tabelle 4.3 — C-Anteil im Koks, Gleichung (4.21)
# x_C,Koks = ξ1 · P^ξ2 · ξ3^(1e-6·P) · T^ξ4
XI_KOKS: tuple[float, float, float, float] = (0.3027, 0.1047, 0.9356, -0.08578)

# Tabelle 4.3 — C-Anteil im Teer, Gleichung (4.22)
# x_C,Teer = ξ1 · P^ξ2 · T^ξ3 + ξ4
XI_TEER: tuple[float, float, float, float] = (-0.2168, 0.04148, 0.1566, 1.333)

# Gültigkeitsbereich Tabelle 4.3
T_WS_MIN_K: float = 1073.0  # [K]
T_WS_MAX_K: float = 1273.0  # [K]
P_WS_MIN_PA: float = 0.1e6  # [Pa]
P_WS_MAX_PA: float = 2.5e6  # [Pa]


def _clip_ws_state(P_Pa: float, T_K: float) -> tuple[float, float]:
    """将 (P,T) 裁到 Tabelle 4.3 有效域（界外用边界值，不外推）。"""
    P = float(np.clip(P_Pa, P_WS_MIN_PA, P_WS_MAX_PA))
    T = float(np.clip(T_K, T_WS_MIN_K, T_WS_MAX_K))
    return P, T


def hegermann_applicability(P_Pa: float, T_K: float) -> dict[str, float | bool | str]:
    """Tabelle 4.3 适用范围。仅 Rheinische Braunkohle（``fuel_type="coal"``）。

    Source: Hamel (1999) §4.2.1.2, Tabelle 4.3
    """
    P_used, T_used = _clip_ws_state(P_Pa, T_K)
    return {
        "fuel": "Rheinische Braunkohle",
        "fuel_type": "coal",
        "source": "Hamel (1999) Tabelle 4.3 / Hegermann (1990)",
        "T_ws_min_K": float(T_WS_MIN_K),
        "T_ws_max_K": float(T_WS_MAX_K),
        "P_ws_min_Pa": float(P_WS_MIN_PA),
        "P_ws_max_Pa": float(P_WS_MAX_PA),
        "T_K": float(T_K),
        "P_Pa": float(P_Pa),
        "T_used_K": float(T_used),
        "P_used_Pa": float(P_used),
        "T_clipped": abs(float(T_K) - float(T_used)) > 1.0e-12,
        "P_clipped": abs(float(P_Pa) - float(P_used)) > 1.0e-9,
    }


def xc_koks_hegermann(P_Pa: float, T_K: float) -> float:
    """焦炭碳份额 x_C,Koks [kg C in Koks / kg C in Brennstoff waf]。

    Hamel (1999) Eq.4.21
    """
    P, T = _clip_ws_state(P_Pa, T_K)
    xi1, xi2, xi3, xi4 = XI_KOKS
    return float(xi1 * (P**xi2) * (xi3 ** (1.0e-6 * P)) * (T**xi4))


def xc_teer_hegermann(P_Pa: float, T_K: float) -> float:
    """焦油碳份额 x_C,Teer [kg C in Teer / kg C in Brennstoff waf]。

    Hamel (1999) Eq.4.22
    """
    P, T = _clip_ws_state(P_Pa, T_K)
    xi1, xi2, xi3, xi4 = XI_TEER
    return float(xi1 * (P**xi2) * (T**xi3) + xi4)


def tar_carbon_frac_of_volatile_c(P_Pa: float, T_K: float) -> float:
    """VM 碳池中进入焦油的碳分率 f_tar = x_C,Teer / (1 - x_C,Koks)。

    在「非焦炭碳最终全部进入挥发分池」的 Hamel 闭合假设下，使
    tar_C / fuel_C = x_C,Teer。返回值 clip 到 [0, 0.9]。
    """
    x_koks = xc_koks_hegermann(P_Pa, T_K)
    x_teer = xc_teer_hegermann(P_Pa, T_K)
    volatile = max(1.0 - x_koks, 1.0e-12)
    f_tar = x_teer / volatile
    return float(np.clip(f_tar, 0.0, 0.9))
