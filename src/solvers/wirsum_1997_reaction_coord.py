"""阶段 F3：离开炭与内层气源汇用同一套 ``rate_R1``。

外层给出 ``m_leave``（轨迹，不重算）。库存 ``M = …`` 与面积 ``6M/(ρ d_p)``
与阶段 E 相同。``K_fb`` 具名为 ``engineering_K_fb_leave``，不是 Hamel 停留。
只做 R1（Hamel 5.8–5.23），不做 R2–R4 / R5–R12。

``r`` 只由 ``char_r1_extent`` 计算；内层残差不得另写 Arrhenius。
生产自由板路径未改。

Source: Hamel (1999) Eq. 5.8–5.23; Wirsum (1997) Eq. 2.21–2.23;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 F
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import MW_C_KG_PER_MOL, gas_diffusivity_correlation
from src.kinetics.char_reactions import rate_R1
from src.physics.wirsum_1997_cyclone import CYCLONE_MODEL_LAPPLE, CycloneGeometry
from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED
from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
from src.solvers.wirsum_1997_bed_freeboard_coord import outer_update_solids
from src.solvers.wirsum_1997_stiffness_full_gas import char_area_m2
from src.solvers.wirsum_1997_stiffness_ladder import (
    ExchangeGeometry,
    classify_newton_stop,
)
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

GAS: tuple[str, ...] = ("O2", "CO", "CO2", "N2")
N_GAS_R: int = 4
I_O2, I_CO, I_CO2, I_N2 = 0, 1, 2, 3
K_FB_MODEL: Literal["engineering_K_fb_leave"] = "engineering_K_fb_leave"
K_FB_DEFAULT: float = 1.0  # [1/s]


@dataclass(frozen=True)
class R1Extent:
    """共享 R1 反应量。全程 SI。"""

    e1: float  # [mol C/s]
    alpha: float  # [-]
    area: float  # [m²]
    C_O2: float  # [mol/m³]
    r_surf: float  # [mol/(m²·s)]


def gas_source_from_extent(ext: R1Extent) -> npt.NDArray[np.float64]:
    """R1 气源 [mol/s]，与阶段 E ``_char_gas_source`` 的 R1 列相同。"""
    e1, a = float(ext.e1), float(ext.alpha)
    src = np.zeros(N_GAS_R, dtype=np.float64)
    src[I_O2] = -a * e1
    src[I_CO] = 2.0 * (1.0 - a) * e1
    src[I_CO2] = (2.0 * a - 1.0) * e1
    return src


def char_r1_extent(
    *,
    T: float,
    P: float,
    N: npt.NDArray[np.float64],
    M_kg: float,
    d_p: float,
    rho_s: float,
) -> R1Extent:
    """唯一允许的 R1 入口。速率只调用生产 ``rate_R1``。"""
    n = np.maximum(np.asarray(N, dtype=np.float64), 0.0)
    ntot = float(np.sum(n))
    y_o2 = float(n[I_O2]) / ntot if ntot > 0.0 else 0.0
    C_O2 = y_o2 * float(P) / (Rg * float(T))  # [mol/m³]
    D_g = gas_diffusivity_correlation(float(T), float(P))  # [m²/s]
    r_surf, alpha = rate_R1(float(T), C_O2, float(d_p), D_g)
    area = char_area_m2(M_kg, rho_s=rho_s, d_p=d_p)
    return R1Extent(
        e1=float(r_surf) * area,
        alpha=float(alpha),
        area=area,
        C_O2=C_O2,
        r_surf=float(r_surf),
    )


def unpack_state(x: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], float]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (N_GAS_R + 1,)
    return z[:N_GAS_R].copy(), float(z[-1])


def r1_fb_residual(
    x: npt.NDArray[np.float64],
    *,
    T: float,
    P: float,
    N_in: npt.NDArray[np.float64],
    m_leave_char: float,
    K_fb: float,
    d_p: float,
    rho_s: float,
) -> npt.NDArray[np.float64]:
    """自由空间单相气 + 离开炭库存。``K_fb`` [1/s]。"""
    N, M = unpack_state(x)
    ext = char_r1_extent(T=T, P=P, N=N, M_kg=M, d_p=d_p, rho_s=rho_s)
    src = gas_source_from_extent(ext)
    F_g = np.asarray(N_in, dtype=np.float64) - N + src
    F_m = float(m_leave_char) - float(K_fb) * M - MW_C_KG_PER_MOL * ext.e1
    return np.concatenate([F_g, np.array([F_m], dtype=np.float64)])


def r1_fb_jacobian(
    x: npt.NDArray[np.float64],
    *,
    eps: float = 1.0e-8,
    **kw: Any,
) -> npt.NDArray[np.float64]:
    """前向差分雅可比。扰动与 ``x`` 同单位。"""
    x0 = np.asarray(x, dtype=np.float64)
    f0 = r1_fb_residual(x0, **kw)
    n = x0.size
    J = np.zeros((n, n), dtype=np.float64)
    for k in range(n):
        xp = x0.copy()
        step = eps * max(abs(float(x0[k])), 1.0)
        xp[k] += step
        J[:, k] = (r1_fb_residual(xp, **kw) - f0) / step
    return J


def default_fb_inlets() -> npt.NDArray[np.float64]:
    """床顶总量示意：空气进自由空间 [mol/s]。"""
    return np.array([0.30, 0.0, 0.0, 10.0], dtype=np.float64)


def leaving_char_from_outer(*, recycle: bool = False) -> float:
    """外层轨迹的离开炭 [kg/s]。细粉、Hamel 扩展、Lapple。"""
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    geo_ex = ExchangeGeometry()
    N_d_vec = np.array([0.2, 0.8], dtype=np.float64)
    N_b_vec = np.array([0.6, 0.2], dtype=np.float64)
    streams, _eta, _d50, _name = outer_update_solids(
        geo=geo,
        N_d_vec=N_d_vec,
        N_b_vec=N_b_vec,
        m_char=0.32,
        m_ash=0.08,
        freeboard_model=FREEBOARD_MODEL_HAMEL_EXTENDED,
        recycle=recycle,
        cyclone_mode=CYCLONE_MODEL_LAPPLE,
        cyclone_geo=CycloneGeometry(),
        T=geo_ex.T,
        P=geo_ex.P,
        N_tot=20.0,
    )
    return float(streams.m_leave[0])


def solve_r1_fb(
    *,
    T: float,
    P: float,
    N_in: npt.NDArray[np.float64],
    m_leave_char: float,
    K_fb: float,
    d_p: float,
    rho_s: float,
) -> dict[str, Any]:
    kw = dict(T=T, P=P, N_in=N_in, m_leave_char=m_leave_char, K_fb=K_fb, d_p=d_p, rho_s=rho_s)
    M0 = max(float(m_leave_char) / max(float(K_fb), 1.0e-12), 1.0e-9)
    x0 = np.concatenate([np.asarray(N_in, dtype=np.float64), np.array([M0])])

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return r1_fb_residual(z, **kw)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return r1_fb_jacobian(z, **kw)

    sol = wirsum_1997_strict_newton(res, jac, x0, max_iter=40)
    stop = classify_newton_stop(sol, model_is_full_rate=True)
    N, M = unpack_state(sol.x)
    ext = char_r1_extent(T=T, P=P, N=N, M_kg=M, d_p=d_p, rho_s=rho_s)
    src = gas_source_from_extent(ext)
    gas_c = float(N[I_CO] + N[I_CO2])
    return {
        **stop,
        "x": sol.x,
        "N": N,
        "M_kg": M,
        "extent": ext,
        "src_mol_s": src,
        "r_outer_mol_s": ext.e1,
        "r_inner_mol_s": ext.e1,
        "r_gap_mol_s": 0.0,
        "carbon_gas_minus_e1": gas_c - ext.e1,
        "K_fb_model": K_FB_MODEL,
    }


def build_stage_f3_payload() -> dict[str, Any]:
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    geo_ex = ExchangeGeometry()
    m_leave = leaving_char_from_outer(recycle=False)
    N_in = default_fb_inlets()
    solved = solve_r1_fb(
        T=geo_ex.T,
        P=geo_ex.P,
        N_in=N_in,
        m_leave_char=m_leave,
        K_fb=K_FB_DEFAULT,
        d_p=geo.d_p,
        rho_s=geo.rho_p,
    )
    ext0 = char_r1_extent(
        T=geo_ex.T, P=geo_ex.P, N=N_in, M_kg=m_leave / K_FB_DEFAULT,
        d_p=geo.d_p, rho_s=geo.rho_p,
    )
    return {
        "source": "Wirsum alignment stage F3 leaving-char R1 coordination",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "units": {"gas": "mol/s", "solid": "kg/s", "holdup": "kg", "r": "mol/s"},
        "m_leave_char_kg_s": m_leave,
        "K_fb_1_per_s": K_FB_DEFAULT,
        "K_fb_model": K_FB_MODEL,
        "newton": {
            "physical_convergence": solved["physical_convergence"],
            "converged_newton": solved["converged_newton"],
            "failure_reason": solved["failure_reason"],
            "residual_l2": solved["residual_l2"],
            "iterations": solved["iterations"],
        },
        "root": {
            "N_mol_s": solved["N"].tolist(),
            "M_kg": solved["M_kg"],
            "e1_mol_s": solved["extent"].e1,
            "alpha": solved["extent"].alpha,
            "C_O2_mol_m3": solved["extent"].C_O2,
            "r_gap_mol_s": solved["r_gap_mol_s"],
            "carbon_gas_minus_e1": solved["carbon_gas_minus_e1"],
        },
        "inlet_extent_e1_mol_s": ext0.e1,
        "inlet_minus_root_e1": ext0.e1 - solved["extent"].e1,
        "reactions": ["R1"],
        "not": ["R2", "R3", "R4", "R5", "R12", "fastox"],
        "note": (
            "r 只来自 char_r1_extent / rate_R1。"
            "K_fb 是工程停留，不是 Hamel。生产路径未改。"
        ),
    }
