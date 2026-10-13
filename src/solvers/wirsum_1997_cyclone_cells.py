"""Hamel Bild 7.7：旋风 1–2 cell + 返料。不经 ``Reactor``，不用 Lapple。

旋风与返料是和床、自由空间同一套衡算，不是分级效率曲线。
气相：Eq. 2.4，退化 ``N_b=0``、``ex_bd=0``、``r=0``。
固相：Eq. 2.6，``m = K M``（物种自己的库存）。``K_auf/K_ab`` 是冻结上下文，
论文未给出旋风 K 的计算式，不得用 Lapple / 常数 η 冒充原文。

回料：``m_rez`` 进入床残差；关回流时 ``m_rez=0``，下行是 ``m_aus``。
这是 Jacobian Nebenelemente 的来源（Hamel §2.3）。

Source: Hamel (1999) Eq. 2.4, 2.6, §2.3, §7.1 Bild 7.7;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 F
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED
from src.solvers.wirsum_1997_bed_freeboard import (
    BedFreeboardGeometry,
    eject_bed_top,
    particle_fate_for_mean_launch,
    split_eject_by_fate,
)
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

K_MODEL: Literal["frozen_context"] = "frozen_context"


def cyclone_k_implied(k: CycloneCellK) -> dict[str, Any]:
    """冻结 K 的隐含分流与停留。不是分离规律，也不是状态更新 P。"""
    ksum = float(k.K_auf) + float(k.K_ab)
    if ksum <= 0.0:
        raise ValueError("K_auf + K_ab must be positive")
    return {
        "K_source": k.model,
        "implied_split_down": float(k.K_ab) / ksum,
        "implied_split_up": float(k.K_auf) / ksum,
        "tau_cyc_s": 1.0 / ksum,
        "not_a_separation_law": True,
        "not_a_state_update_P": True,
    }


@dataclass(frozen=True)
class CycloneCellK:
    """冻结输运系数 [1/s]。不是 Hamel 给出的计算公式。"""

    K_auf: float = 0.25  # [1/s] 随气上行
    K_ab: float = 0.75  # [1/s] 下行
    K_return: float = 1.0  # [1/s] 返料腿只下行
    model: str = K_MODEL


@dataclass(frozen=True)
class CycloneLoopSpec:
    """Bild 7.7 拓扑。``n_cyc`` 论文允许自选。"""

    n_cyc: int = 2
    recycle: bool = False
    N_fb: float = 1.8  # [mol/s] 床顶总量，单相进入旋风
    m_leave: float = 0.83  # [kg/s] 自由空间离开炭，轨迹侧已算
    m_feed: float = 0.32  # [kg/s] 床给料，只用于床残差
    k: CycloneCellK = CycloneCellK()


def _assert_n_cyc(n: int) -> None:
    if n not in (1, 2):
        raise ValueError(f"n_cyc must be 1 or 2 (Bild 7.7 / frei wählbar), got {n}")


def n_unknowns(n_cyc: int) -> int:
    """``N`` 每格 1 个 + ``M`` 每旋风格 1 个 + 返料 ``M``。"""
    _assert_n_cyc(n_cyc)
    return 2 * n_cyc + 1


def unpack(x: npt.NDArray[np.float64], n_cyc: int) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (n_unknowns(n_cyc),)
    N = z[:n_cyc].copy()
    M = z[n_cyc:2 * n_cyc].copy()
    return N, M, float(z[-1])


def cyclone_gas_residual(
    N: npt.NDArray[np.float64],
    *,
    N_fb: float,
) -> npt.NDArray[np.float64]:
    """Eq. 2.4，单相、无交换、无反应、无 ``N_rez``。

    格 0：``N_zu = N_fb``。其后：上游出流作轴向入流。
    """
    n = N.size
    F = np.zeros(n, dtype=np.float64)
    for i in range(n):
        N_zu = float(N_fb) if i == 0 else 0.0
        N_in = 0.0 if i == 0 else float(N[i - 1])
        F[i] = N_zu + N_in - float(N[i])
    return F


def cyclone_solid_residual(
    M: npt.NDArray[np.float64],
    M_ret: float,
    *,
    spec: CycloneLoopSpec,
) -> npt.NDArray[np.float64]:
    """Eq. 2.6。格间 ``m_auf = K_auf M``、``m_ab = K_ab M``。返料只收床底格下行。"""
    n = M.size
    K_auf = float(spec.k.K_auf)
    K_ab = float(spec.k.K_ab)
    m_auf = K_auf * M
    m_ab = K_ab * M
    F = np.zeros(n, dtype=np.float64)
    for i in range(n):
        m_zu = float(spec.m_leave) if i == 0 else 0.0
        m_auf_in = 0.0 if i == 0 else float(m_auf[i - 1])
        m_ab_in = 0.0 if i + 1 >= n else float(m_ab[i + 1])
        F[i] = m_zu + m_auf_in + m_ab_in - float(m_auf[i]) - float(m_ab[i])
    m_under = float(m_ab[0])
    F_ret = m_under - float(spec.k.K_return) * float(M_ret)
    return np.concatenate([F, np.array([F_ret], dtype=np.float64)])


def recycle_streams(M_ret: float, spec: CycloneLoopSpec) -> tuple[float, float]:
    """返料出流 [kg/s]：开回流进 ``m_rez``，关回流进 ``m_aus``。"""
    m = float(spec.k.K_return) * float(M_ret)
    if spec.recycle:
        return m, 0.0
    return 0.0, m


def bed_residual_with_rez(spec: CycloneLoopSpec, m_rez: float) -> float:
    """床固相残差里显式含 ``m_rez``。不是求解后再回写。"""
    return float(spec.m_feed) + float(m_rez) - float(spec.m_leave)


def loop_residual(x: npt.NDArray[np.float64], spec: CycloneLoopSpec) -> npt.NDArray[np.float64]:
    """旋风气/固 + 返料。床 ``F`` 单独报告，不混进本核未知量。"""
    N, M, M_ret = unpack(x, spec.n_cyc)
    F_g = cyclone_gas_residual(N, N_fb=spec.N_fb)
    F_s = cyclone_solid_residual(M, M_ret, spec=spec)
    return np.concatenate([F_g, F_s])


def loop_jacobian(x: npt.NDArray[np.float64], spec: CycloneLoopSpec, *, eps: float = 1.0e-8) -> npt.NDArray[np.float64]:
    x0 = np.asarray(x, dtype=np.float64)
    f0 = loop_residual(x0, spec)
    n = x0.size
    J = np.zeros((n, n), dtype=np.float64)
    for k in range(n):
        xp = x0.copy()
        step = eps * max(abs(float(x0[k])), 1.0)
        xp[k] += step
        J[:, k] = (loop_residual(xp, spec) - f0) / step
    return J


def d_bed_d_Mret(spec: CycloneLoopSpec) -> float:
    """``∂F_bed / ∂M_ret``。开回流为 ``K_return``，关回流为 0。"""
    return float(spec.k.K_return) if spec.recycle else 0.0


def leaving_char_from_trajectory() -> float:
    """自由空间离开炭 [kg/s]。Hamel 扩展抛射 + 命运，不经旋风经验式。"""
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    m_eject = eject_bed_top(geo=geo, m_char=0.32, m_ash=0.08, model=FREEBOARD_MODEL_HAMEL_EXTENDED)
    fate = particle_fate_for_mean_launch(geo=geo, model=FREEBOARD_MODEL_HAMEL_EXTENDED)
    _ret, leave = split_eject_by_fate(m_eject, fate)
    return float(leave[0])


def solve_loop(spec: CycloneLoopSpec) -> dict[str, Any]:
    n = spec.n_cyc
    N0 = np.full(n, spec.N_fb, dtype=np.float64)
    M0 = np.full(n, spec.m_leave / max(spec.k.K_auf + spec.k.K_ab, 1.0e-12), dtype=np.float64)
    Mr0 = spec.m_leave / max(spec.k.K_return, 1.0e-12)
    x0 = np.concatenate([N0, M0, np.array([Mr0])])

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return loop_residual(z, spec)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return loop_jacobian(z, spec)

    sol = wirsum_1997_strict_newton(res, jac, x0, max_iter=20)
    stop = classify_newton_stop(sol, model_is_full_rate=True)
    N, M, M_ret = unpack(sol.x, n)
    m_rez, m_aus = recycle_streams(M_ret, spec)
    fly = float(spec.k.K_auf) * float(M[-1])
    under = float(spec.k.K_ab) * float(M[0])
    return {
        **stop,
        "N_mol_s": N,
        "M_cyc_kg": M,
        "M_ret_kg": M_ret,
        "m_rez_kg_s": m_rez,
        "m_aus_kg_s": m_aus,
        "m_flyash_kg_s": fly,
        "m_underflow_kg_s": under,
        "leave_minus_fly_under": float(spec.m_leave) - fly - under,
        "F_bed": bed_residual_with_rez(spec, m_rez),
        "dF_bed_dMret": d_bed_d_Mret(spec),
        "K_model": spec.k.model,
        **cyclone_k_implied(spec.k),
    }


def build_stage_f_cyclone_payload() -> dict[str, Any]:
    m_leave = leaving_char_from_trajectory()
    off = solve_loop(CycloneLoopSpec(n_cyc=2, recycle=False, m_leave=m_leave))
    on = solve_loop(CycloneLoopSpec(n_cyc=2, recycle=True, m_leave=m_leave))
    one = solve_loop(CycloneLoopSpec(n_cyc=1, recycle=False, m_leave=m_leave))
    return {
        "source": "Wirsum alignment Hamel Bild 7.7 cyclone cells",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "units": {"gas": "mol/s", "solid": "kg/s", "holdup": "kg", "K": "1/s"},
        "topology": "bed → freeboard → cyclone(n) → return → bed",
        "n_cyc_frei_waehlbar": [1, 2],
        "K_model": K_MODEL,
        "K_implied": cyclone_k_implied(CycloneCellK()),
        "m_leave_char_kg_s": m_leave,
        "recycle_closed": {
            "physical_convergence": off["physical_convergence"],
            "residual_l2": off["residual_l2"],
            "m_rez_kg_s": off["m_rez_kg_s"],
            "m_aus_kg_s": off["m_aus_kg_s"],
            "leave_minus_fly_under": off["leave_minus_fly_under"],
            "dF_bed_dMret": off["dF_bed_dMret"],
            "N_mol_s": off["N_mol_s"].tolist(),
        },
        "recycle_open": {
            "physical_convergence": on["physical_convergence"],
            "m_rez_kg_s": on["m_rez_kg_s"],
            "F_bed": on["F_bed"],
            "dF_bed_dMret": on["dF_bed_dMret"],
            "N_mol_s": on["N_mol_s"].tolist(),
        },
        "one_cell": {
            "physical_convergence": one["physical_convergence"],
            "leave_minus_fly_under": one["leave_minus_fly_under"],
            "N_mol_s": one["N_mol_s"].tolist(),
        },
        "not": ["engineering_lapple_1951", "constant_eta", "R1"],
        "note": (
            "旋风是 cell 衡算。K 是冻结上下文，不是论文公式。"
            "Lapple 对照核不参与本残差。"
        ),
    }
