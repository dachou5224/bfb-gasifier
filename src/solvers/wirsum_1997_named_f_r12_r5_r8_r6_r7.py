"""F_E0+R12+R5+R8+R6 数学根上接入 R7（CH4 蒸汽重整）。不改 ``rate_R7``，不经 ``Reactor``。

R12/R5/R8/R6 保持完整速率。R7 用 ``r7_scale`` 延拓。
仓库 R7 = thesis R9。Hamel Eq. 5.51：``R_m = k p_CH4`` [mol/(s·kg_ash)]，
体积源 ``r_V = ρ_ash R_m``，只进密相。气泡与无灰自由空间为零。
C_CH4=0 时密相保留 CH4 右导数。不做 Startwert。

Source: Hamel (1999) Eq. 5.50–5.51; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.kinetics.arrhenius import k_standard
from src.kinetics.gas_reactions import R7_E_J, R7_k0_Pa, rate_R7
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    rho_ash_used,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_r12 import _conc_and_dN, named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import freeze_r12_r5_r8_q
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import (
    freeze_r12_r5_r8_r6_q,
    jacobian_r12_r5_r8_plus_scaled_r6,
    named_r12_r5_r8_r6_newton_step,
    solve_named_r12_r5_r8_then_r6,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

# CH4 + H2O → CO + 3 H2。Ref: Hamel (1999) Eq. 5.50。
NU_R7_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R7_MAJOR[IDX["CH4"]] = -1.0
NU_R7_MAJOR[IDX["H2O"]] = -1.0
NU_R7_MAJOR[IDX["CO"]] = 1.0
NU_R7_MAJOR[IDX["H2"]] = 3.0
_T_FLOOR: float = 300.0  # [K] 与 rate_R7 相同

_R7_CONTINUATION: tuple[float, ...] = (
    0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0
)


def freeze_r12_r5_r8_r6_r7_q() -> NamedStageQ:
    return freeze_stage_q(2, enabled_rxns=frozenset({"R12", "R5", "R8", "R6", "R7"}))


def residual_r12_r5_r8_r6_plus_scaled_r7(
    x: npt.NDArray[np.float64],
    nq_r6: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r7_scale: float,
) -> npt.NDArray[np.float64]:
    """R12+R5+R8+R6 完整速率，R7 乘 ``r7_scale``。"""
    z = np.asarray(x, dtype=np.float64)
    f6 = coupled_residual_named(z, nq_r6, rate_scale=1.0)
    s = float(r7_scale)
    if abs(s) < 1.0e-30:
        return f6
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    return f6 + s * (both - f6)


def _r7_r_and_dC(
    T: float, C: npt.NDArray[np.float64], *, rho_ash: float
) -> tuple[float, npt.NDArray[np.float64]]:
    """密相体积速率 ``r_V = ρ_ash k p_CH4``。ρ_ash<=0 时速率与导数均为 0。"""
    dr = np.zeros(N_MAJ, dtype=np.float64)
    rho = max(float(rho_ash), 0.0)
    if rho <= 0.0:
        return 0.0, dr
    t = max(float(T), _T_FLOOR)
    k_mass = float(k_standard(R7_k0_Pa, float(R7_E_J), t))
    dr[IDX["CH4"]] = rho * k_mass * Rg * t
    ch4 = float(C[IDX["CH4"]])
    if ch4 <= 0.0:
        return 0.0, dr
    r = float(rate_R7(
        t, ch4, float(C[IDX["H2O"]]), float(C[IDX["CO"]]), float(C[IDX["H2"]]),
        rho_ash_kg_m3=rho,
    ))
    return r, dr


def _r7_dr_dT(T: float, C: npt.NDArray[np.float64], *, rho_ash: float) -> float:
    """冻结摩尔分率。p_CH4 = y_CH4 P，与 T 无关；r ∝ k(T) ρ_ash。"""
    rho = max(float(rho_ash), 0.0)
    ch4 = float(C[IDX["CH4"]])
    if rho <= 0.0 or ch4 <= 0.0:
        return 0.0
    t = max(float(T), _T_FLOOR)
    r = float(rate_R7(
        t, ch4, float(C[IDX["H2O"]]), float(C[IDX["CO"]]), float(C[IDX["H2"]]),
        rho_ash_kg_m3=rho,
    ))
    dlnk = float(R7_E_J) / (Rg * t * t)
    return float(r * dlnk)


def r7_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """密相 R7 对 N_d/T/M_ash 的解析块。气泡与无灰为 0。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    T = float(z[-1])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    s = float(rate_scale)
    if abs(s) < 1.0e-30:
        return J
    u = unpack_named(z, nq.stage) if n >= 41 else {"M_ash": float(nq.rho_ash) * V_d}
    rho = rho_ash_used(nq, u) if n >= 41 else float(nq.rho_ash)
    C, dCdN = _conc_and_dN(z[0:N_MAJ], T, P)
    r, drdC = _r7_r_and_dC(T, C, rho_ash=rho)
    dxi_dN = (s * float(V_d)) * (drdC @ dCdN)
    J[0:N_MAJ, 0:N_MAJ] += NU_R7_MAJOR[:, None] * dxi_dN
    dxi_dT = (s * float(V_d)) * _r7_dr_dT(T, C, rho_ash=rho)
    J[0:N_MAJ, -1] += NU_R7_MAJOR * dxi_dT
    if n >= 41:
        r_unit = float(rate_R7(
            T, float(C[IDX["CH4"]]), float(C[IDX["H2O"]]),
            float(C[IDX["CO"]]), float(C[IDX["H2"]]),
            rho_ash_kg_m3=1.0,
        ))
        J[0:N_MAJ, 39] += NU_R7_MAJOR * (s * r_unit)
    _ = r
    return J


def jacobian_r12_r5_r8_r6_plus_scaled_r7(
    x: npt.NDArray[np.float64],
    nq_r6: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r7_scale: float,
) -> npt.NDArray[np.float64]:
    J6 = jacobian_r12_r5_r8_plus_scaled_r6(
        x, freeze_r12_r5_r8_q(), nq_r6, r6_scale=1.0
    )
    s = float(r7_scale)
    if abs(s) < 1.0e-30:
        return J6
    return J6 + s * r7_analytic_jacobian(x, nq_all, rate_scale=1.0)


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _newton_at_r7_scale(
    nq_r6: NamedStageQ,
    nq_all: NamedStageQ,
    x: npt.NDArray[np.float64],
    r7_scale: float,
    max_iter: int,
) -> Any:
    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_r12_r5_r8_r6_plus_scaled_r7(
            z, nq_r6, nq_all, r7_scale=r7_scale
        )

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_r12_r5_r8_r6_plus_scaled_r7(
            z, nq_r6, nq_all, r7_scale=r7_scale
        )

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_r8_r6_newton_step,
    )


def r7_source_snapshot(
    x: npt.NDArray[np.float64],
    nq_r6: NamedStageQ,
    nq_all: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f6 = coupled_residual_named(z, nq_r6, rate_scale=1.0)
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    dF = both - f6
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_CH4_d": float(dF[IDX["CH4"]]),
        "dF_H2O_d": float(dF[IDX["H2O"]]),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_H2_d": float(dF[IDX["H2"]]),
        "dF_CH4_b": float(dF[7 + IDX["CH4"]]),
        "dF_H2_b": float(dF[7 + IDX["H2"]]),
        "dF_energy": float(dF[-1]),
        "N_d_CH4": float(z[IDX["CH4"]]),
        "N_d_H2": float(z[IDX["H2"]]),
        "N_d_O2": float(z[IDX["O2"]]),
        "N_d_H2O": float(z[IDX["H2O"]]),
        "T_K": float(z[-1]),
    }


def solve_named_r12_r5_r8_r6_then_r7(
    *,
    max_iter: int = 20,
    continuation: bool = True,
    r7_scale: float | None = None,
    r6: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E0+R12+R5+R8+R6 数学根做 R7 延拓。"""
    nq_r6 = freeze_r12_r5_r8_r6_q()
    nq_all = freeze_r12_r5_r8_r6_r7_q()
    if r6 is None:
        r6 = solve_named_r12_r5_r8_then_r6(max_iter=max_iter, apply_startwert=True)
    x_start = np.asarray(r6["x"], dtype=np.float64).copy()
    snap0 = r7_source_snapshot(x_start, nq_r6, nq_all)
    scales = (float(r7_scale),) if r7_scale is not None else (
        _R7_CONTINUATION if continuation else (1.0,)
    )
    n = n_unknowns(nq_all.stage)
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    for s in scales:
        last = _newton_at_r7_scale(nq_r6, nq_all, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_r12_r5_r8_r6_plus_scaled_r7(last.x, nq_r6, nq_all, r7_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        history.append(
            {
                "r7_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "N_d_CH4": float(last.x[IDX["CH4"]]),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
                "N_d_CO": float(last.x[IDX["CO"]]),
                "N_d_H2O": float(last.x[IDX["H2O"]]),
            }
        )
        if int(last.rank_last) == n and bool(stop_s["converged_newton"]):
            x_cur = np.asarray(last.x, dtype=np.float64)

    assert last is not None
    stop = classify_newton_stop(
        last, model_is_full_rate=abs(float(scales[-1]) - 1.0) < 1.0e-15
    )
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "T_K": float(last.x[-1]),
        "T_r6_K": float(r6["T_K"]),
        "energy_W": float(energy_row_named(last.x, nq_all)),
        "f_name": "F_E0+R12+R5+R8+R6+R7",
        "enabled_rxns": ("R12", "R5", "R8", "R6", "R7"),
        "r6_physical_convergence": bool(r6["physical_convergence"]),
        "r7_at_r6_root": snap0,
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
    }
