"""阶段 F 后半：外层轨迹更新固体，内层解气/T。不经 ``Reactor``。

内层：床层示踪交换（Hamel Eq. 2.2，闭式解）+ 气体焓定 T。
外层：Hamel 扩展抛射 / 命运分类 + Lapple 旋风。``η`` 随 ``Q(T)`` 变。
反应量本核标 ``named_absent``（0=0），不是 Hamel 炭动力学。

旋风是 ``engineering_lapple_1951``，不是 Hamel Bild 7.7 的 cell 拓扑。
生产 ``simulate_freeboard`` 与常数捕集分数未改。

Source: Hamel (1999) Eq. 2.2; Lapple (1951);
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 F
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import cp_molar
from src.physics.wirsum_1997_cyclone import (
    CYCLONE_MODEL_CONSTANT,
    CYCLONE_MODEL_LAPPLE,
    CycloneGeometry,
    cyclone_eta,
    grade_efficiency_lapple,
)
from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED
from src.solvers.wirsum_1997_bed_freeboard import (
    OUTER_CRITERIA,
    BedFreeboardGeometry,
    Covered,
    InterfaceStreams,
    cyclone_split,
    eject_bed_top,
    gas_interface_gap,
    particle_fate_for_mean_launch,
    route_total_gas,
    solid_edge_gaps,
    split_eject_by_fate,
)
from src.solvers.wirsum_1997_stiffness_ladder import (
    ExchangeGeometry,
    exact_tracer_exchange,
    tracer_exchange_residual,
)

# 气体焓身份用 N2 Cp(T_in)，不是混合物积分。
_CP_SPECIES: str = "N2"


@dataclass(frozen=True)
class InnerGasState:
    """内层气/T。示踪 [mol/s]，温度 [K]。"""

    N_d: float
    N_b: float
    N_fb: float
    T: float
    residual_gas: npt.NDArray[np.float64]
    residual_T: float
    Q_wall: float


def inner_temperature(
    *,
    T_in: float,
    N_tot: float,
    Q_wall: float,
) -> tuple[float, float]:
    """气体焓身份：``N cp (T_in − T) − Q_W = 0``。

    ``cp`` 取 ``cp_molar('N2', T_in)`` [J/(mol·K)]。``Q_W`` 是具名工程热损，
    不是 Hamel Eq. 2.7 全式。返回 ``(T, residual)``。
    """
    if T_in <= 0.0 or N_tot <= 0.0:
        raise ValueError("T_in and N_tot must be positive")
    cp = float(cp_molar(_CP_SPECIES, float(T_in)))  # [J/(mol·K)]
    T = float(T_in) - float(Q_wall) / (float(N_tot) * cp)
    if T <= 0.0:
        raise ValueError(f"Q_wall={Q_wall} drove T={T} K non-physical")
    residual = float(N_tot) * cp * (float(T_in) - T) - float(Q_wall)
    return T, residual


def inner_solve_gas_T(
    *,
    geo_ex: ExchangeGeometry,
    N_d_in: float,
    N_b_in: float,
    N_tot_d: float,
    N_tot_b: float,
    Q_wall: float,
) -> InnerGasState:
    """示踪交换闭式解 + 焓身份 T。不重算抛射。"""
    N = exact_tracer_exchange(
        geo=geo_ex, N_d_in=N_d_in, N_b_in=N_b_in, N_tot_d=N_tot_d, N_tot_b=N_tot_b
    )
    N_d, N_b = float(N[0]), float(N[1])
    N_fb = N_d + N_b
    T, rT = inner_temperature(T_in=geo_ex.T, N_tot=N_tot_d + N_tot_b, Q_wall=Q_wall)
    F = tracer_exchange_residual(
        N, geo=geo_ex, N_d_in=N_d_in, N_b_in=N_b_in, N_tot_d=N_tot_d, N_tot_b=N_tot_b
    )
    return InnerGasState(
        N_d=N_d,
        N_b=N_b,
        N_fb=N_fb,
        T=T,
        residual_gas=F,
        residual_T=rT,
        Q_wall=float(Q_wall),
    )


def outer_update_solids(
    *,
    geo: BedFreeboardGeometry,
    N_d_vec: npt.NDArray[np.float64],
    N_b_vec: npt.NDArray[np.float64],
    m_char: float,
    m_ash: float,
    freeboard_model: str,
    recycle: bool,
    cyclone_mode: str,
    cyclone_geo: CycloneGeometry,
    T: float,
    P: float,
    N_tot: float,
    eta_constant: float = 0.95,
) -> tuple[InterfaceStreams, float, float, str]:
    """轨迹更新固体。``η`` 由内层 ``T, N_tot`` 进入 Lapple。"""
    N_tot_vec = route_total_gas(N_d_vec, N_b_vec)
    m_eject = eject_bed_top(geo=geo, m_char=m_char, m_ash=m_ash, model=freeboard_model)
    fate = particle_fate_for_mean_launch(geo=geo, model=freeboard_model)
    m_ret, m_leave = split_eject_by_fate(m_eject, fate)
    eta, d50, cyc_name = cyclone_eta(
        mode=cyclone_mode,
        d_p=geo.d_p,
        mu=geo.mu_g,
        rho_p=geo.rho_p,
        rho_g=geo.rho_g,
        N_tot=N_tot,
        T=T,
        P=P,
        geo=cyclone_geo,
        eta_constant=eta_constant,
    )
    catch, fly, rec = cyclone_split(m_leave, eta=eta, recycle=recycle)
    streams = InterfaceStreams(
        N_bed_total=N_tot_vec,
        N_fb_in=N_tot_vec.copy(),
        m_eject=m_eject,
        m_return=m_ret,
        m_leave=m_leave,
        m_catch=catch,
        m_flyash=fly,
        m_recycle=rec,
        fate=fate,
        model=freeboard_model,
    )
    return streams, eta, d50, cyc_name


def coordinate_once(
    *,
    geo: BedFreeboardGeometry,
    geo_ex: ExchangeGeometry,
    cyclone_geo: CycloneGeometry,
    N_d_in: float,
    N_b_in: float,
    N_tot_d: float,
    N_tot_b: float,
    N_d_vec: npt.NDArray[np.float64],
    N_b_vec: npt.NDArray[np.float64],
    m_char: float,
    m_ash: float,
    recycle: bool,
    Q_wall: float,
    cyclone_mode: str = CYCLONE_MODEL_LAPPLE,
    freeboard_model: str = FREEBOARD_MODEL_HAMEL_EXTENDED,
) -> dict[str, Any]:
    """一次外层→内层。固体不进残差，气体不重算抛射。"""
    inner = inner_solve_gas_T(
        geo_ex=geo_ex,
        N_d_in=N_d_in,
        N_b_in=N_b_in,
        N_tot_d=N_tot_d,
        N_tot_b=N_tot_b,
        Q_wall=Q_wall,
    )
    streams, eta, d50, cyc_name = outer_update_solids(
        geo=geo,
        N_d_vec=N_d_vec,
        N_b_vec=N_b_vec,
        m_char=m_char,
        m_ash=m_ash,
        freeboard_model=freeboard_model,
        recycle=recycle,
        cyclone_mode=cyclone_mode,
        cyclone_geo=cyclone_geo,
        T=inner.T,
        P=geo_ex.P,
        N_tot=N_tot_d + N_tot_b,
    )
    return {
        "inner": inner,
        "streams": streams,
        "eta": eta,
        "d50_m": d50,
        "cyclone_model": cyc_name,
        "reaction": {
            "model": "named_absent",
            "r_outer_mol_s": 0.0,
            "r_inner_mol_s": 0.0,
            "gap_mol_s": 0.0,
        },
    }


def outer_criterion_coverage_coord() -> dict[str, Covered]:
    """协调核覆盖。全炉能量/热解/库存仍未做。"""
    return {
        "temperature": "covered",
        "total_gas": "covered",
        "phase_split": "not_covered",
        "hydro_exchange": "frozen_context",
        "solid_bidirectional": "covered",
        "solid_holdup": "not_covered",
        "pyrolysis_release": "not_covered",
        "bed_surface_fallback_entrainment": "covered",
        "element_balance": "not_covered",
        "energy": "not_covered",
    }


def _stream_brief(streams: InterfaceStreams) -> dict[str, Any]:
    return {
        "fate": streams.fate,
        "m_eject_kg_s": streams.m_eject.tolist(),
        "m_leave_kg_s": streams.m_leave.tolist(),
        "m_recycle_kg_s": streams.m_recycle.tolist(),
        "m_flyash_kg_s": streams.m_flyash.tolist(),
        "gas_gap_mol_s": gas_interface_gap(streams),
        "solid_gaps": solid_edge_gaps(streams),
    }


def build_stage_f2_payload() -> dict[str, Any]:
    geo = BedFreeboardGeometry(d_p=2.0e-4)  # [m] 离开到旋风
    geo_ex = ExchangeGeometry()
    cyc = CycloneGeometry()
    N_d_vec = np.array([0.2, 0.8], dtype=np.float64)
    N_b_vec = np.array([0.6, 0.2], dtype=np.float64)
    base = dict(
        geo=geo,
        geo_ex=geo_ex,
        cyclone_geo=cyc,
        N_d_in=1.0,
        N_b_in=0.0,
        N_tot_d=10.0,
        N_tot_b=10.0,
        N_d_vec=N_d_vec,
        N_b_vec=N_b_vec,
        m_char=0.32,
        m_ash=0.08,
        freeboard_model=FREEBOARD_MODEL_HAMEL_EXTENDED,
        recycle=True,
    )
    a = coordinate_once(**base, Q_wall=0.0, cyclone_mode=CYCLONE_MODEL_LAPPLE)
    a2 = coordinate_once(**base, Q_wall=0.0, cyclone_mode=CYCLONE_MODEL_LAPPLE)
    b = coordinate_once(**base, Q_wall=6.0e4, cyclone_mode=CYCLONE_MODEL_LAPPLE)
    const = coordinate_once(**base, Q_wall=0.0, cyclone_mode=CYCLONE_MODEL_CONSTANT)
    eta_repeat = abs(float(a["eta"]) - float(a2["eta"]))
    rec_repeat = float(
        np.max(np.abs(a["streams"].m_recycle - a2["streams"].m_recycle))
    )
    # 离开级 200 μm 的 η≈1；T 耦合看切割粒径处：Q_W>0 → d50 变大 → 原 d50 处 η<0.5。
    eta_at_old_cut = grade_efficiency_lapple(float(a["d50_m"]), float(b["d50_m"]))
    return {
        "source": "Wirsum alignment stage F2 bed-freeboard coordination",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "units": {"gas": "mol/s", "solid": "kg/s", "T": "K", "Q_wall": "W", "d50": "m"},
        "inner_Q0": {
            "N_d": a["inner"].N_d,
            "N_b": a["inner"].N_b,
            "N_fb": a["inner"].N_fb,
            "T": a["inner"].T,
            "residual_gas_max": float(np.max(np.abs(a["inner"].residual_gas))),
            "residual_T": a["inner"].residual_T,
        },
        "outer_Q0": {
            **_stream_brief(a["streams"]),
            "eta": a["eta"],
            "d50_m": a["d50_m"],
            "cyclone_model": a["cyclone_model"],
        },
        "repeat_outer_gap": {"eta": eta_repeat, "recycle_kg_s": rec_repeat},
        "inner_Qwall": {
            "T": b["inner"].T,
            "eta_leave": b["eta"],
            "d50_m": b["d50_m"],
            "eta_at_Q0_cut": eta_at_old_cut,
        },
        "constant_eta_minus_lapple_Q0": float(const["eta"] - a["eta"]),
        "reaction": a["reaction"],
        "outer_criteria": outer_criterion_coverage_coord(),
        "outer_criteria_keys": list(OUTER_CRITERIA),
        "note": (
            "内层解示踪与 T；外层用 T 更新 Lapple η。"
            "Hamel 旋风仍是 cell 拓扑，不是本式。生产路径未改。"
        ),
    }
