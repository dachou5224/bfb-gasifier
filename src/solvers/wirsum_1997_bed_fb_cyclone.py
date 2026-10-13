"""阶段 F5：床–自由空间界面与 Bild 7.7 旋风共用同一套 ``N_fb``、``m_leave``。

界面只做到离开；旋风 cell 接收这两股流，不再用常数 η / Lapple 切一刀。
轨迹用 Hamel 扩展抛射；``hamel_extended`` 命运用 Eq. 3.86 极限 z(φ)。
回料以 ``m_rez`` 进床残差。生产路径未改。

Source: Hamel (1999) Eq. 2.4, 2.6, §7.1 Bild 7.7;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 F
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED
from src.solvers.wirsum_1997_bed_freeboard import (
    OUTER_CRITERIA,
    BedFreeboardGeometry,
    Covered,
    eject_bed_top,
    particle_fate_for_mean_launch,
    route_total_gas,
    split_eject_by_fate,
)
from src.solvers.wirsum_1997_cyclone_cells import (
    CycloneLoopSpec,
    bed_residual_with_rez,
    solve_loop,
)


@dataclass(frozen=True)
class BedFbInlets:
    """床顶气与抛射组成。气体 [mol/s]，固体组成 [kg/s]。"""

    N_d: npt.NDArray[np.float64]
    N_b: npt.NDArray[np.float64]
    m_char: float
    m_ash: float
    m_feed_char: float


def default_inlets() -> BedFbInlets:
    return BedFbInlets(
        N_d=np.array([0.2, 0.8], dtype=np.float64),
        N_b=np.array([0.6, 0.2], dtype=np.float64),
        m_char=0.32,
        m_ash=0.08,
        m_feed_char=0.32,
    )


def interface_to_cyclone(
    *,
    geo: BedFreeboardGeometry,
    inn: BedFbInlets,
    freeboard_model: str = FREEBOARD_MODEL_HAMEL_EXTENDED,
) -> dict[str, Any]:
    """床顶 → 自由空间。到 ``m_leave`` 为止，不切旋风。"""
    N_fb_vec = route_total_gas(inn.N_d, inn.N_b)
    m_eject = eject_bed_top(
        geo=geo, m_char=inn.m_char, m_ash=inn.m_ash, model=freeboard_model
    )
    fate = particle_fate_for_mean_launch(geo=geo, model=freeboard_model)
    m_return, m_leave = split_eject_by_fate(m_eject, fate)
    return {
        "N_fb_vec": N_fb_vec,
        "N_fb": float(np.sum(N_fb_vec)),
        "m_eject": m_eject,
        "m_return": m_return,
        "m_leave": m_leave,
        "fate": fate,
        "model": freeboard_model,
        "eject_minus_return_leave": float(np.max(np.abs(m_eject - m_return - m_leave))),
    }


def couple(
    *,
    geo: BedFreeboardGeometry,
    inn: BedFbInlets,
    recycle: bool,
    n_cyc: int = 2,
    freeboard_model: str = FREEBOARD_MODEL_HAMEL_EXTENDED,
) -> dict[str, Any]:
    """同一 ``N_fb``、``m_leave[char]`` 交给旋风 cell。"""
    face = interface_to_cyclone(geo=geo, inn=inn, freeboard_model=freeboard_model)
    spec = CycloneLoopSpec(
        n_cyc=n_cyc,
        recycle=recycle,
        N_fb=float(face["N_fb"]),
        m_leave=float(face["m_leave"][0]),
        m_feed=float(inn.m_feed_char),
    )
    cyc = solve_loop(spec)
    return {
        "face": face,
        "spec": spec,
        "cyc": cyc,
        "N_fb_handoff_gap": abs(float(cyc["N_mol_s"][0]) - float(face["N_fb"])),
        "m_leave_handoff_gap": abs(float(spec.m_leave) - float(face["m_leave"][0])),
        "F_bed": bed_residual_with_rez(spec, float(cyc["m_rez_kg_s"])),
    }


def outer_criterion_coverage_coupled() -> dict[str, Covered]:
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


def _brief(out: dict[str, Any]) -> dict[str, Any]:
    cyc = out["cyc"]
    return {
        "fate": out["face"]["fate"],
        "N_fb": out["face"]["N_fb"],
        "m_leave_char_kg_s": float(out["face"]["m_leave"][0]),
        "N_fb_handoff_gap": out["N_fb_handoff_gap"],
        "m_leave_handoff_gap": out["m_leave_handoff_gap"],
        "eject_minus_return_leave": out["face"]["eject_minus_return_leave"],
        "leave_minus_fly_under": cyc["leave_minus_fly_under"],
        "m_rez_kg_s": cyc["m_rez_kg_s"],
        "F_bed": out["F_bed"],
        "dF_bed_dMret": cyc["dF_bed_dMret"],
        "physical_convergence": cyc["physical_convergence"],
        "residual_l2": cyc["residual_l2"],
    }


def build_stage_f5_payload() -> dict[str, Any]:
    inn = default_inlets()
    fine = BedFreeboardGeometry(d_p=2.0e-4)
    coarse = BedFreeboardGeometry()
    fine_off = couple(geo=fine, inn=inn, recycle=False)
    fine_on = couple(geo=fine, inn=inn, recycle=True)
    coarse_off = couple(geo=coarse, inn=inn, recycle=False)
    return {
        "source": "Wirsum alignment stage F5 bed-FB-cyclone handoff",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "units": {"gas": "mol/s", "solid": "kg/s"},
        "topology": "bed → freeboard → cyclone(2) → return → bed",
        "fine_recycle_closed": _brief(fine_off),
        "fine_recycle_open": _brief(fine_on),
        "coarse_recycle_closed": _brief(coarse_off),
        "outer_criteria": outer_criterion_coverage_coupled(),
        "outer_criteria_keys": list(OUTER_CRITERIA),
        "not": ["engineering_lapple_1951", "constant_eta", "evaluate_interface"],
        "note": (
            "N_fb 与 m_leave 只算一次，旋风 cell 原样接收。"
            "粗粒回落则旋风固相入口为 0。生产路径未改。"
        ),
    }
