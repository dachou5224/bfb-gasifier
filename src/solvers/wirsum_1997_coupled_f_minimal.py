"""Q6：等温、单炭粒径、床内仅 R1 的联立 F。不经 ``Reactor``。

38 未知 / 38 方程。从一个 x 生成一次边流和同一 R1 进度。
旋风 K 冻结；返料无滞留，``m_rez = K_ab M_{c1}`` 直接回床。
自由空间气相 ``N_fb,j − N_d,j − N_b,j = 0``（物种守恒接口，不是两相同速）。

省略项全部 ``omitted, named``。本轮不切生产默认。

Source: Hamel (1999) Eq. 2.4–2.6, 5.8–5.23;
docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q6
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import MW_C_KG_PER_MOL, gas_diffusivity_correlation
from src.kinetics.char_reactions import rate_R1
from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED
from src.solvers.wirsum_1997_bed_freeboard import (
    BedFreeboardGeometry,
    eject_bed_top,
    particle_fate_for_mean_launch,
)
from src.solvers.wirsum_1997_cyclone_cells import CycloneCellK, cyclone_k_implied
from src.solvers.wirsum_1997_stiffness_full_gas import (
    IDX,
    MAJOR,
    N_MAJ,
    char_area_m2,
    fd_jacobian,
)
from src.solvers.wirsum_1997_stiffness_ladder import ExchangeGeometry, classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

N_X: int = 38
N_CYC: int = 2
OMITTED: tuple[str, ...] = (
    "pyrolysis_drying",
    "R2_R12",
    "tar",
    "multi_size_migration",
    "freeboard_reaction",
    "energy_wall_heat",
    "return_leg_holdup",
    "cyclone_K_update_P",
    "full_outer_P",
)


@dataclass(frozen=True)
class CoupledMinimalQ:
    """冻结上下文 q。系数不乘旧库存当流量。"""

    geo_ex: ExchangeGeometry
    geo_fb: BedFreeboardGeometry
    N_in_d: npt.NDArray[np.float64]
    N_in_b: npt.NDArray[np.float64]
    m_feed: float  # [kg/s]
    K_eject: float  # [1/s] 外层线性化，不是原文常数
    w_leave: float  # [-]
    w_return: float  # [-]
    k: CycloneCellK
    D_g: float  # [m²/s]
    H_bed: float  # [m]


def _elem_C(N: npt.NDArray[np.float64]) -> float:
    i = IDX
    return float(N[i["CO"]] + N[i["CO2"]] + N[i["CH4"]])


def _elem_H(N: npt.NDArray[np.float64]) -> float:
    i = IDX
    return float(2.0 * N[i["H2"]] + 2.0 * N[i["H2O"]] + 4.0 * N[i["CH4"]])


def _elem_O(N: npt.NDArray[np.float64]) -> float:
    i = IDX
    return float(N[i["CO"]] + 2.0 * N[i["CO2"]] + N[i["H2O"]] + 2.0 * N[i["O2"]])


def freeze_q(
    *,
    geo_ex: ExchangeGeometry | None = None,
    geo_fb: BedFreeboardGeometry | None = None,
    N_in_d: npt.NDArray[np.float64] | None = None,
    N_in_b: npt.NDArray[np.float64] | None = None,
    m_feed: float = 0.02,
    H_bed: float = 1.0,
    k: CycloneCellK | None = None,
) -> CoupledMinimalQ:
    """从水力/轨迹冻结 ``K_eject`` 与逃逸权重。``m_eject = K_eject M``。"""
    ex = geo_ex or ExchangeGeometry()
    fb = geo_fb or BedFreeboardGeometry(d_p=2.0e-4)
    if N_in_d is None:
        N_in_d = np.array([0.05, 0.02, 0.04, 0.03, 0.01, 0.0, 2.0], dtype=np.float64)
    if N_in_b is None:
        N_in_b = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.40, 1.5], dtype=np.float64)
    ck = k or CycloneCellK()
    m_eject = eject_bed_top(geo=fb, m_char=1.0, m_ash=0.0, model=FREEBOARD_MODEL_HAMEL_EXTENDED)
    M_geom = (
        float(fb.rho_p)
        * max(1.0 - float(fb.eps_d), 0.0)
        * max(1.0 - float(fb.eps_b), 0.0)
        * float(fb.area)
        * float(H_bed)
    )
    K_eject = float(m_eject[0]) / max(M_geom, 1.0e-12)
    fate = particle_fate_for_mean_launch(geo=fb, model=FREEBOARD_MODEL_HAMEL_EXTENDED)
    if fate == "returns":
        w_leave, w_return = 0.0, 1.0
    elif fate == "reaches_next_cell":
        w_leave, w_return = 1.0, 0.0
    else:
        w_leave, w_return = 1.0, 0.0
    D_g = gas_diffusivity_correlation(float(ex.T), float(ex.P))
    return CoupledMinimalQ(
        geo_ex=ex,
        geo_fb=fb,
        N_in_d=np.asarray(N_in_d, dtype=np.float64),
        N_in_b=np.asarray(N_in_b, dtype=np.float64),
        m_feed=float(m_feed),
        K_eject=float(K_eject),
        w_leave=float(w_leave),
        w_return=float(w_return),
        k=ck,
        D_g=float(D_g),
        H_bed=float(H_bed),
    )


def unpack(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (N_X,), f"expected {N_X}, got {z.shape}"
    return {
        "N_d": z[0:7].copy(),
        "N_b": z[7:14].copy(),
        "M": float(z[14]),
        "N_fb": z[15:22].copy(),
        "N_c1": z[22:29].copy(),
        "N_c2": z[29:36].copy(),
        "M_c1": float(z[36]),
        "M_c2": float(z[37]),
    }


def pack(
    *,
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    M: float,
    N_fb: npt.NDArray[np.float64],
    N_c1: npt.NDArray[np.float64],
    N_c2: npt.NDArray[np.float64],
    M_c1: float,
    M_c2: float,
) -> npt.NDArray[np.float64]:
    return np.concatenate(
        [
            np.asarray(N_d, dtype=np.float64),
            np.asarray(N_b, dtype=np.float64),
            np.array([M], dtype=np.float64),
            np.asarray(N_fb, dtype=np.float64),
            np.asarray(N_c1, dtype=np.float64),
            np.asarray(N_c2, dtype=np.float64),
            np.array([M_c1, M_c2], dtype=np.float64),
        ]
    )


def edges_from_x(x: npt.NDArray[np.float64], q: CoupledMinimalQ) -> dict[str, float]:
    """当前 x 与冻结系数唯一决定的边流 [kg/s]。"""
    u = unpack(x)
    m_eject = float(q.K_eject) * float(u["M"])
    m_leave = float(q.w_leave) * m_eject
    m_return = float(q.w_return) * m_eject
    m_auf0 = float(q.k.K_auf) * float(u["M_c1"])
    m_ab0 = float(q.k.K_ab) * float(u["M_c1"])
    m_auf1 = float(q.k.K_auf) * float(u["M_c2"])
    m_ab1 = float(q.k.K_ab) * float(u["M_c2"])
    return {
        "m_eject": m_eject,
        "m_leave": m_leave,
        "m_return": m_return,
        "m_auf0": m_auf0,
        "m_ab0": m_ab0,
        "m_auf1": m_auf1,
        "m_ab1": m_ab1,
        "m_rez": m_ab0,
        "m_fly": m_auf1,
    }


def r1_extent_from_x(x: npt.NDArray[np.float64], q: CoupledMinimalQ) -> dict[str, float]:
    """唯一 R1 进度。只调用生产 ``rate_R1``。仅 dense。"""
    u = unpack(x)
    N_d = u["N_d"]
    ntot = float(np.sum(N_d))
    C_O2 = 0.0
    if abs(ntot) > 1.0e-30:
        y = float(N_d[IDX["O2"]]) / ntot
        C_O2 = y * float(q.geo_ex.P) / (Rg * float(q.geo_ex.T))
    r_surf, alpha = rate_R1(float(q.geo_ex.T), C_O2, float(q.geo_fb.d_p), float(q.D_g))
    area = char_area_m2(u["M"], rho_s=q.geo_fb.rho_p, d_p=q.geo_fb.d_p)
    e1 = float(r_surf) * float(area)
    return {"e1": e1, "alpha": float(alpha), "area": float(area), "C_O2": float(C_O2)}


def r1_gas_source(ext: dict[str, float]) -> npt.NDArray[np.float64]:
    e1, a = float(ext["e1"]), float(ext["alpha"])
    src = np.zeros(N_MAJ, dtype=np.float64)
    src[IDX["O2"]] = -a * e1
    src[IDX["CO"]] = 2.0 * (1.0 - a) * e1
    src[IDX["CO2"]] = (2.0 * a - 1.0) * e1
    return src


def coupled_residual(
    x: npt.NDArray[np.float64],
    q: CoupledMinimalQ,
    *,
    rate_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    """一张 F。边流与 R1 各算一次。``rate_scale`` 只用于延拓。"""
    u = unpack(x)
    ed = edges_from_x(x, q)
    ext = r1_extent_from_x(x, q)
    ext = {**ext, "e1": float(rate_scale) * float(ext["e1"])}
    src = r1_gas_source(ext)
    C_d = (u["N_d"] / max(float(np.sum(u["N_d"])), 1.0e-30)) * (
        float(q.geo_ex.P) / (Rg * float(q.geo_ex.T))
    )
    C_b = (u["N_b"] / max(float(np.sum(u["N_b"])), 1.0e-30)) * (
        float(q.geo_ex.P) / (Rg * float(q.geo_ex.T))
    )
    N_ex = float(q.geo_ex.K_bd) * float(q.geo_ex.V_b) * (C_b - C_d)
    F_d = q.N_in_d - u["N_d"] + N_ex + src
    F_b = q.N_in_b - u["N_b"] - N_ex
    F_M = (
        float(q.m_feed)
        + ed["m_return"]
        + ed["m_rez"]
        - ed["m_eject"]
        - MW_C_KG_PER_MOL * ext["e1"]
    )
    F_fb = u["N_fb"] - u["N_d"] - u["N_b"]
    F_c1g = u["N_fb"] - u["N_c1"]
    F_c2g = u["N_c1"] - u["N_c2"]
    F_Mc1 = ed["m_leave"] + ed["m_ab1"] - ed["m_auf0"] - ed["m_ab0"]
    F_Mc2 = ed["m_auf0"] - ed["m_auf1"] - ed["m_ab1"]
    return np.concatenate(
        [
            F_d,
            F_b,
            np.array([F_M], dtype=np.float64),
            F_fb,
            F_c1g,
            F_c2g,
            np.array([F_Mc1, F_Mc2], dtype=np.float64),
        ]
    )


def atom_gap_from_F(x: npt.NDArray[np.float64], q: CoupledMinimalQ) -> dict[str, float]:
    """任意 x：床气残差减 FB 接口拷贝，再加旋风气残差，应闭合进出口。"""
    F = coupled_residual(x, q)
    u = unpack(x)
    ed = edges_from_x(x, q)
    gas_F = F[0:7] + F[7:14] - F[15:22] + F[22:29] + F[29:36]
    solid_F = float(F[14]) + float(F[36]) + float(F[37])
    N_in = q.N_in_d + q.N_in_b
    return {
        "C": abs(
            _elem_C(gas_F)
            + solid_F / MW_C_KG_PER_MOL
            + _elem_C(u["N_c2"])
            + ed["m_fly"] / MW_C_KG_PER_MOL
            - _elem_C(N_in)
            - q.m_feed / MW_C_KG_PER_MOL
        ),
        "H": abs(_elem_H(gas_F) + _elem_H(u["N_c2"]) - _elem_H(N_in)),
        "O": abs(_elem_O(gas_F) + _elem_O(u["N_c2"]) - _elem_O(N_in)),
        "N": abs(float(gas_F[IDX["N2"]]) + float(u["N_c2"][IDX["N2"]]) - float(N_in[IDX["N2"]])),
        "mass_kg_s": abs(solid_F + ed["m_fly"] + MW_C_KG_PER_MOL * r1_extent_from_x(x, q)["e1"] - q.m_feed),
    }


def default_x0(q: CoupledMinimalQ) -> npt.NDArray[np.float64]:
    N_tot = q.N_in_d + q.N_in_b
    M = max(q.m_feed / max(q.K_eject * q.w_leave + 1.0e-6, 1.0e-6), 0.05)
    m_leave = q.K_eject * q.w_leave * M
    ksum = q.k.K_auf + q.k.K_ab
    Mc = m_leave / max(ksum, 1.0e-12)
    return pack(
        N_d=q.N_in_d.copy(),
        N_b=q.N_in_b.copy(),
        M=M,
        N_fb=N_tot.copy(),
        N_c1=N_tot.copy(),
        N_c2=N_tot.copy(),
        M_c1=Mc,
        M_c2=Mc,
    )


def omitted_named() -> dict[str, str]:
    return {name: "omitted, named" for name in OMITTED}


def q_diagnostics(q: CoupledMinimalQ) -> dict[str, Any]:
    return {
        "K_eject_1_per_s": q.K_eject,
        "K_eject_source": "frozen_linearization_m_eject_over_M_geom",
        "w_leave": q.w_leave,
        "w_return": q.w_return,
        "fate_drag": "hamel_eq386_asymptotic",
        "engineering_K_fb_leave": "not_used",
        "cyclone": cyclone_k_implied(q.k),
        "units": {"gas": "mol/s", "solid": "kg/s", "K": "1/s"},
        "n_unknowns": N_X,
        "omitted": omitted_named(),
        "species": list(MAJOR),
    }


def solve_coupled(
    q: CoupledMinimalQ,
    *,
    x0: npt.NDArray[np.float64] | None = None,
    rate_scale: float = 1.0,
    max_iter: int = 40,
) -> dict[str, Any]:
    """同一张 F 上的 Wirsum strict Newton。不是四套根拼接。"""
    x_start = default_x0(q) if x0 is None else np.asarray(x0, dtype=np.float64)
    scale = np.maximum(np.abs(x_start), 0.05)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return coupled_residual(z, q, rate_scale=rate_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return fd_jacobian(res, z, scale=scale)

    sol = wirsum_1997_strict_newton(res, jac, x_start, max_iter=max_iter)
    stop = classify_newton_stop(sol, model_is_full_rate=abs(rate_scale - 1.0) < 1.0e-15)
    F2 = coupled_residual(sol.x, q, rate_scale=rate_scale)
    return {
        **stop,
        "x": np.asarray(sol.x, dtype=np.float64),
        "F": np.asarray(sol.F, dtype=np.float64),
        "rate_scale": float(rate_scale),
        "q_purity_gap": float(np.max(np.abs(F2 - sol.F))),
        "rank_last": int(sol.rank_last),
        "cond_last": float(sol.cond_last),
        "edges": edges_from_x(sol.x, q),
        "extent": r1_extent_from_x(sol.x, q),
        "atom_gap": atom_gap_from_F(sol.x, q),
    }


def recycle_changes_bed_root(
    *,
    k_lo: CycloneCellK,
    k_hi: CycloneCellK,
) -> dict[str, Any]:
    """同一进料，只改冻结 K_ab，床库存根必须变。"""
    a = solve_coupled(freeze_q(k=k_lo))
    b = solve_coupled(freeze_q(k=k_hi))
    return {
        "M_lo_kg": float(a["x"][14]),
        "M_hi_kg": float(b["x"][14]),
        "gap_kg": abs(float(a["x"][14]) - float(b["x"][14])),
        "both_physical": bool(a["physical_convergence"] and b["physical_convergence"]),
    }


def build_stage_q6_payload() -> dict[str, Any]:
    q = freeze_q()
    out = solve_coupled(q)
    return {
        "source": "Wirsum alignment Q6 isothermal single-size coupled F",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "n_unknowns": N_X,
        "physical_convergence": out["physical_convergence"],
        "residual_l2": out["residual_l2"],
        "q_purity_gap": out["q_purity_gap"],
        "rank_last": out["rank_last"],
        "cond_last": out["cond_last"],
        "atom_gap": out["atom_gap"],
        "edges": out["edges"],
        "extent": out["extent"],
        "q": q_diagnostics(q),
        "omitted": omitted_named(),
        "note": (
            "一张 38×38 F。飞灰与 R1 消耗闭合固体账。"
            "生产 Reactor 未改。省略项见 omitted。"
        ),
    }
