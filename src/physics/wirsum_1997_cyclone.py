"""旋风分级效率对照核。不是 Hamel 原文。

Hamel (1999) 把 Zyklon 建成 1–2 个 cell（Bild 7.7），VTT 算例甚至没有旋风几何。
生产路径用常数 ``freeboard_cyclone_capture_*_frac``。本核给出具名的
Lapple (1951) 切割粒径与分级效率，供阶段 F 后半与常数 η 对照。

``η = 1 / (1 + (d_{50}/d_p)^2)``，
``d_{50} = [9 μ W / (2 π N_e (ρ_p − ρ_g) v_i)]^{1/2}``。

Source: Lapple (1951); Shepherd & Lapple (1939) ``N_e`` 习惯值；
Hamel (1999) §7.1 只提供旋风 cell 拓扑，不含本式。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from src.core.constants import Rg

CYCLONE_MODEL_CONSTANT: str = "constant_eta"
CYCLONE_MODEL_LAPPLE: str = "engineering_lapple_1951"


@dataclass(frozen=True)
class CycloneGeometry:
    """单级旋风入口几何。全程 SI。"""

    inlet_height: float = 0.040  # [m]
    inlet_width: float = 0.011  # [m]
    n_turns: float = 5.0  # [-] Shepherd–Lapple 常用有效圈数


def volumetric_flow(*, N_tot: float, T: float, P: float) -> float:
    """理想气体体积流 [m³/s]。"""
    if T <= 0.0 or P <= 0.0:
        raise ValueError("T, P must be positive SI")
    return float(N_tot) * Rg * float(T) / float(P)


def inlet_velocity(*, N_tot: float, T: float, P: float, geo: CycloneGeometry) -> float:
    """入口气速 ``v_i = Q / (a W)`` [m/s]。"""
    area = float(geo.inlet_height) * float(geo.inlet_width)
    if area <= 0.0:
        raise ValueError("cyclone inlet area must be positive")
    return volumetric_flow(N_tot=N_tot, T=T, P=P) / area


def lapple_cut_size(
    *,
    mu: float,
    inlet_width: float,
    n_turns: float,
    v_i: float,
    rho_p: float,
    rho_g: float,
) -> float:
    """Lapple 切割粒径 ``d_{50}`` [m]。

    Source: Lapple (1951)
    """
    if mu <= 0.0 or inlet_width <= 0.0 or n_turns <= 0.0 or v_i <= 0.0:
        raise ValueError("Lapple cut-size inputs must be positive")
    dens = float(rho_p) - float(rho_g)
    if dens <= 0.0:
        raise ValueError("rho_p must exceed rho_g")
    arg = 9.0 * float(mu) * float(inlet_width) / (
        2.0 * math.pi * float(n_turns) * dens * float(v_i)
    )
    return math.sqrt(arg)


def grade_efficiency_lapple(d_p: float, d50: float) -> float:
    """Lapple 分级效率 ``η ∈ [0, 1]``。``d_p = d_{50}`` 时为 0.5。"""
    if d_p <= 0.0 or d50 <= 0.0:
        raise ValueError("d_p and d50 must be positive")
    ratio = float(d50) / float(d_p)
    return 1.0 / (1.0 + ratio * ratio)


def cyclone_eta(
    *,
    mode: str,
    d_p: float,
    mu: float,
    rho_p: float,
    rho_g: float,
    N_tot: float,
    T: float,
    P: float,
    geo: CycloneGeometry,
    eta_constant: float = 0.95,
) -> tuple[float, float, str]:
    """返回 ``(η, d50_or_nan, model_name)``。常数模式的 ``d50`` 为 NaN。"""
    name = str(mode).strip().lower()
    if name == CYCLONE_MODEL_CONSTANT:
        if not (0.0 <= float(eta_constant) <= 1.0):
            raise ValueError(f"constant eta must be in [0, 1], got {eta_constant}")
        return float(eta_constant), float("nan"), CYCLONE_MODEL_CONSTANT
    if name != CYCLONE_MODEL_LAPPLE:
        raise ValueError(
            f"Unknown cyclone_model={mode!r}; expected "
            f"{CYCLONE_MODEL_CONSTANT!r} or {CYCLONE_MODEL_LAPPLE!r}"
        )
    v_i = inlet_velocity(N_tot=N_tot, T=T, P=P, geo=geo)
    d50 = lapple_cut_size(
        mu=mu,
        inlet_width=geo.inlet_width,
        n_turns=geo.n_turns,
        v_i=v_i,
        rho_p=rho_p,
        rho_g=rho_g,
    )
    return grade_efficiency_lapple(d_p, d50), d50, CYCLONE_MODEL_LAPPLE
