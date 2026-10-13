"""阶段 F：床–自由空间界面约化核。不经 ``Reactor``，不改生产默认。

先关回流：气体只传 ``N_d+N_b`` 总量；固体由轨迹抛射，回落回床，离开进旋风。
Wirsum exact 只用 Roof 3.59；Hamel 扩展加 Wake 3.54。二者不得混名。
轨迹更新固体；气相界面不把床层相分裂带进自由空间。

外层判据在本核只覆盖已实现的场；温度/热解/能量仍标 ``not_covered``。

Source: Hamel (1999) Eq. 2.2, 3.53–3.54, 3.73; Wirsum (1997) Eq. 3.59;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 F
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from src.physics.wirsum_1997_freeboard import (
    FREEBOARD_MODEL_HAMEL_EXTENDED,
    FREEBOARD_MODEL_WIRSUM_1997_EXACT,
    HAMEL_PHI_S_MEAN_DEFAULT,
    classify_particle_fate_branching,
    drag_z_for_spec,
    eject_size_class_mass_flows,
    resolve_freeboard_model,
)

OUTER_CRITERIA: tuple[str, ...] = (
    "temperature",
    "total_gas",
    "phase_split",
    "hydro_exchange",
    "solid_bidirectional",
    "solid_holdup",
    "pyrolysis_release",
    "bed_surface_fallback_entrainment",
    "element_balance",
    "energy",
)

Covered = Literal["covered", "not_covered", "frozen_context"]


@dataclass(frozen=True)
class BedFreeboardGeometry:
    """单床顶 + 一段自由空间。全程 SI。"""

    area: float = 0.126  # [m²] 约 D=0.4 m
    H_freeboard: float = 2.0  # [m]
    d_p: float = 2.25e-3  # [m]
    rho_p: float = 1400.0  # [kg/m³]
    eps_d: float = 0.45  # [-]
    eps_b: float = 0.13  # [-]
    u_b: float = 0.80  # [m/s]
    d_b: float = 0.21  # [m]
    u_g: float = 1.20  # [m/s]
    rho_g: float = 0.80  # [kg/m³]
    mu_g: float = 4.0e-5  # [Pa·s]
    phi_s: float = HAMEL_PHI_S_MEAN_DEFAULT  # [-] 形状因子均值，未解析分布


@dataclass(frozen=True)
class InterfaceStreams:
    """界面物流。气体 [mol/s]，固体 [kg/s]。"""

    N_bed_total: npt.NDArray[np.float64]
    N_fb_in: npt.NDArray[np.float64]
    m_eject: npt.NDArray[np.float64]
    m_return: npt.NDArray[np.float64]
    m_leave: npt.NDArray[np.float64]
    m_catch: npt.NDArray[np.float64]
    m_flyash: npt.NDArray[np.float64]
    m_recycle: npt.NDArray[np.float64]
    fate: str
    model: str


def route_total_gas(
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """床顶 → 自由空间：只传总量，不传相分裂。"""
    nd = np.asarray(N_d, dtype=np.float64)
    nb = np.asarray(N_b, dtype=np.float64)
    assert nd.shape == nb.shape
    return nd + nb


def eject_bed_top(
    *,
    geo: BedFreeboardGeometry,
    m_char: float,
    m_ash: float,
    model: str,
) -> npt.NDArray[np.float64]:
    """床面抛射 [kg/s]，顺序 ``(char, ash)``。"""
    feed = np.array([max(float(m_char), 0.0), max(float(m_ash), 0.0)], dtype=np.float64)
    d_p = np.array([geo.d_p, geo.d_p], dtype=np.float64)
    return eject_size_class_mass_flows(
        d_p_classes=d_p,
        feed_classes=feed,
        rho_p=geo.rho_p,
        eps_d=geo.eps_d,
        eps_b=geo.eps_b,
        u_b=geo.u_b,
        d_b=geo.d_b,
        area=geo.area,
        mode=model,
    )


def particle_fate_for_mean_launch(*, geo: BedFreeboardGeometry, model: str) -> str:
    """平均抛射 ``u_{p,0}``。``hamel_extended`` 用 Eq. 3.86 极限 z(φ)。

    回床是 ``z≤0``，不是 ``u_p=0``。本段气速冻结，跨 ``u_p=u_g`` 换分支。
    """
    spec = resolve_freeboard_model(model)
    u_p0 = spec.u_p0_over_u_b * float(geo.u_b)
    z_drag = drag_z_for_spec(spec, phi_s=geo.phi_s)
    return classify_particle_fate_branching(
        u_p0,
        u_g=geo.u_g,
        rho_g=geo.rho_g,
        rho_p=geo.rho_p,
        d_p=geo.d_p,
        mu_g=geo.mu_g,
        z=z_drag,
        cell_height=geo.H_freeboard,
    )


def split_eject_by_fate(
    m_eject: npt.NDArray[np.float64],
    fate: str,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """回落回床 vs 离开到旋风。``undecided`` 两边都是 0，不是逃逸证据。"""
    m = np.asarray(m_eject, dtype=np.float64)
    z = np.zeros_like(m)
    if fate == "returns":
        return m.copy(), z
    if fate == "reaches_next_cell":
        return z, m.copy()
    return z.copy(), z.copy()


def cyclone_split(
    m_leave: npt.NDArray[np.float64],
    *,
    eta: float,
    recycle: bool,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """``(catch, flyash, recycle)`` [kg/s]。关回流时 recycle=0。"""
    if not (0.0 <= float(eta) <= 1.0):
        raise ValueError(f"cyclone eta must be in [0, 1], got {eta}")
    leave = np.asarray(m_leave, dtype=np.float64)
    catch = float(eta) * leave
    fly = leave - catch
    rec = catch if recycle else np.zeros_like(leave)
    return catch, fly, rec


def evaluate_interface(
    *,
    geo: BedFreeboardGeometry,
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    m_char: float,
    m_ash: float,
    model: str,
    recycle: bool,
    eta_cyclone: float = 0.95,
) -> InterfaceStreams:
    N_tot = route_total_gas(N_d, N_b)
    m_eject = eject_bed_top(geo=geo, m_char=m_char, m_ash=m_ash, model=model)
    fate = particle_fate_for_mean_launch(geo=geo, model=model)
    m_ret, m_leave = split_eject_by_fate(m_eject, fate)
    catch, fly, rec = cyclone_split(m_leave, eta=eta_cyclone, recycle=recycle)
    return InterfaceStreams(
        N_bed_total=N_tot,
        N_fb_in=N_tot.copy(),
        m_eject=m_eject,
        m_return=m_ret,
        m_leave=m_leave,
        m_catch=catch,
        m_flyash=fly,
        m_recycle=rec,
        fate=fate,
        model=resolve_freeboard_model(model).name,
    )


def gas_interface_gap(streams: InterfaceStreams) -> float:
    """``||N_fb_in − (N_d+N_b)||`` [mol/s]。"""
    return float(np.max(np.abs(streams.N_fb_in - streams.N_bed_total)))


def solid_edge_gaps(streams: InterfaceStreams) -> dict[str, float]:
    """内部边：抛射=接收，回落=床接收；离开=捕集+飞灰。"""
    return {
        "eject_minus_return_leave_kg_s": float(
            np.max(np.abs(streams.m_eject - streams.m_return - streams.m_leave))
        ),
        "leave_minus_catch_fly_kg_s": float(
            np.max(np.abs(streams.m_leave - streams.m_catch - streams.m_flyash))
        ),
        "recycle_kg_s": float(np.max(np.abs(streams.m_recycle))),
    }


def outer_criterion_coverage(*, recycle: bool) -> dict[str, Covered]:
    """本核对外层 ``q=P(x)`` 的覆盖。不是全炉已协调。

    ``element_balance`` 仍是全炉判据，这里只有界面质量恒等式，故 ``not_covered``。
    """
    _ = recycle
    return {
        "temperature": "not_covered",
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


def _stream_record(s: InterfaceStreams) -> dict[str, Any]:
    return {
        "fate": s.fate,
        "model": s.model,
        "m_eject_kg_s": s.m_eject.tolist(),
        "m_return_kg_s": s.m_return.tolist(),
        "m_leave_kg_s": s.m_leave.tolist(),
        "m_catch_kg_s": s.m_catch.tolist(),
        "m_flyash_kg_s": s.m_flyash.tolist(),
        "m_recycle_kg_s": s.m_recycle.tolist(),
        "gas_gap_mol_s": gas_interface_gap(s),
        "solid_gaps": solid_edge_gaps(s),
    }


def build_stage_f_payload() -> dict[str, Any]:
    geo_coarse = BedFreeboardGeometry()  # d_p=2.25 mm，本几何回落
    geo_fine = BedFreeboardGeometry(d_p=2.0e-4)  # [m] 离开到旋风
    N_d = np.array([0.2, 0.8], dtype=np.float64)  # [mol/s] 示踪
    N_b = np.array([0.6, 0.2], dtype=np.float64)
    m_char, m_ash = 0.32, 0.08  # [kg/s] 组成，不是必须等于抛射总量
    closed: dict[str, Any] = {}
    for model in (FREEBOARD_MODEL_WIRSUM_1997_EXACT, FREEBOARD_MODEL_HAMEL_EXTENDED):
        s = evaluate_interface(
            geo=geo_coarse, N_d=N_d, N_b=N_b, m_char=m_char, m_ash=m_ash,
            model=model, recycle=False,
        )
        closed[model] = _stream_record(s)
    fine_off = evaluate_interface(
        geo=geo_fine, N_d=N_d, N_b=N_b, m_char=m_char, m_ash=m_ash,
        model=FREEBOARD_MODEL_HAMEL_EXTENDED, recycle=False,
    )
    fine_on = evaluate_interface(
        geo=geo_fine, N_d=N_d, N_b=N_b, m_char=m_char, m_ash=m_ash,
        model=FREEBOARD_MODEL_HAMEL_EXTENDED, recycle=True,
    )
    return {
        "source": "Wirsum alignment stage F bed-freeboard interface",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "units": {"gas": "mol/s", "solid": "kg/s", "length": "m"},
        "recycle_closed_coarse": closed,
        "fine_leave_d_p_m": geo_fine.d_p,
        "recycle_closed_fine_hamel": _stream_record(fine_off),
        "recycle_open_fine_hamel": _stream_record(fine_on),
        "eject_hamel_minus_wirsum_kg_s": float(
            np.sum(np.array(closed[FREEBOARD_MODEL_HAMEL_EXTENDED]["m_eject_kg_s"]))
            - np.sum(np.array(closed[FREEBOARD_MODEL_WIRSUM_1997_EXACT]["m_eject_kg_s"]))
        ),
        "outer_criteria": outer_criterion_coverage(recycle=False),
        "fate_drag": {
            FREEBOARD_MODEL_WIRSUM_1997_EXACT: "wirsum_1997_cw",
            FREEBOARD_MODEL_HAMEL_EXTENDED: "hamel_eq386_asymptotic",
        },
        "note": (
            "关回流验收界面。粗颗粒回落；细粉离开后旋风才有物流。"
            "Wirsum exact ≠ Hamel 扩展。生产自由板路径未改。"
        ),
    }
