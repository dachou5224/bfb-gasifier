"""F_E0+R12+R5 数学根上接入 R8（WGSR）。不改 ``rate_R8``，不经 ``Reactor``。

R12/R5 保持完整速率。R8 用 ``r8_scale`` 延拓，只进密相（气泡无灰）。
压力项内部 atm：``P_atm = P / P0``。平衡由 ``y_CO - y_eq,CO`` 驱动。
在 R12+R5 根上 R8 源 ~0.04 mol/s；H2 仍被剩余 O2 经 R12 烧掉，T 几乎不动。

Source: Hamel (1999) Eq. 5.45–5.48; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import P0, Rg
from src.kinetics.arrhenius import k_standard
from src.kinetics.gas_reactions import (
    R8_E_Rg,
    R8_a_R8,
    R8_k0,
    R8_p_corr_T,
    _clamp_exp_arg,
    rate_R8,
    wgsr_equilibrium_constant,
    wgsr_equilibrium_y_co,
)
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
)
from src.solvers.wirsum_1997_named_f_r12 import freeze_r12_only_q, named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5 import (
    freeze_r12_r5_q,
    jacobian_r12_plus_scaled_r5,
    named_r12_r5_newton_step,
    solve_named_r12_then_r5,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

# CO + H2O ⇌ CO2 + H2。Ref: Hamel (1999) Eq. 5.45–5.46。
NU_R8_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R8_MAJOR[IDX["CO"]] = -1.0
NU_R8_MAJOR[IDX["H2O"]] = -1.0
NU_R8_MAJOR[IDX["CO2"]] = 1.0
NU_R8_MAJOR[IDX["H2"]] = 1.0
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]
_Y_H2O_FLOOR: float = 1.0e-30  # [-] 与 wgsr_equilibrium_y_co 一致
_T_FLOOR: float = 300.0  # [K]

_R8_CONTINUATION: tuple[float, ...] = (
    0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0
)


def freeze_r12_r5_r8_q() -> NamedStageQ:
    return freeze_stage_q(2, enabled_rxns=frozenset({"R12", "R5", "R8"}))


def residual_r12_r5_plus_scaled_r8(
    x: npt.NDArray[np.float64],
    nq_r5: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r8_scale: float,
) -> npt.NDArray[np.float64]:
    """R12+R5 完整速率，R8 乘 ``r8_scale``。"""
    z = np.asarray(x, dtype=np.float64)
    f5 = coupled_residual_named(z, nq_r5, rate_scale=1.0)
    s = float(r8_scale)
    if abs(s) < 1.0e-30:
        return f5
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    return f5 + s * (both - f5)


def _y_and_dydN(
    N: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    n = np.asarray(N, dtype=np.float64)
    ntot = float(np.sum(n))
    if abs(ntot) < _NTOT_FLOOR:
        return np.zeros(N_MAJ), np.zeros((N_MAJ, N_MAJ))
    y = n / ntot
    dydN = np.zeros((N_MAJ, N_MAJ), dtype=np.float64)
    for i in range(N_MAJ):
        for k in range(N_MAJ):
            dydN[i, k] = (((1.0 if i == k else 0.0) * ntot) - n[i]) / (ntot * ntot)
    return y, dydN


def _r8_k_mass(T: float, P: float) -> tuple[float, float]:
    """``k_mass`` [1/s per kg_ash] 与 ``∂ln k_mass / ∂T`` [1/K]。"""
    t = max(float(T), _T_FLOOR)
    k_arr = float(k_standard(R8_k0, float(R8_E_Rg) * Rg, t))
    p_atm = max(float(P), 0.0) / P0
    p_exp = 0.5 - p_atm / 250.0
    exp_corr = math.exp(_clamp_exp_arg(-8.91 + float(R8_p_corr_T) / t))
    k_mass = float(R8_a_R8) * k_arr * (p_atm ** p_exp) * exp_corr
    dlnk_dT = (float(R8_E_Rg) - float(R8_p_corr_T)) / (t * t)
    return k_mass, dlnk_dT


def r8_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """密相 R8 对 N_d/T 的解析块。气泡 ρ_ash=0，不贡献。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    T = float(z[-1])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    rho = float(nq.rho_ash)
    if rho <= 0.0 or abs(float(rate_scale)) < 1.0e-30:
        return J
    y, dydN = _y_and_dydN(z[0:N_MAJ])
    y_co = max(float(y[IDX["CO"]]), 0.0)
    y_h2o = max(float(y[IDX["H2O"]]), 0.0)
    y_co2 = max(float(y[IDX["CO2"]]), 0.0)
    y_h2 = max(float(y[IDX["H2"]]), 0.0)
    k_eq = float(wgsr_equilibrium_constant(T))
    y_eq = float(wgsr_equilibrium_y_co(y_h2o, y_co2, y_h2, k_eq))
    k_mass, dlnk_dT = _r8_k_mass(T, P)
    r = float(rate_R8(T, P, y_co, y_h2o, y_co2, y_h2, rho))
    drive = y_co - y_eq
    ddrive = np.zeros(N_MAJ, dtype=np.float64)
    ddrive[IDX["CO"]] = 1.0
    denom = max(k_eq * max(y_h2o, _Y_H2O_FLOOR), _Y_H2O_FLOOR)
    ddrive[IDX["CO2"]] -= y_h2 / denom
    ddrive[IDX["H2"]] -= y_co2 / denom
    ddrive[IDX["H2O"]] += y_co2 * y_h2 / (denom * max(y_h2o, _Y_H2O_FLOOR))
    drdN = (rho * k_mass) * (ddrive @ dydN)
    sV = float(rate_scale) * float(V_d)
    J[0:N_MAJ, 0:N_MAJ] += NU_R8_MAJOR[:, None] * (sV * drdN)
    dyeq_dT = y_eq * (4019.0 / (max(T, _T_FLOOR) ** 2))
    dr_dT = rho * (dlnk_dT * k_mass * drive + k_mass * (-dyeq_dT))
    J[0:N_MAJ, -1] += NU_R8_MAJOR * (sV * dr_dT)
    _ = r
    return J


def jacobian_r12_r5_plus_scaled_r8(
    x: npt.NDArray[np.float64],
    nq_r5: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r8_scale: float,
) -> npt.NDArray[np.float64]:
    J5 = jacobian_r12_plus_scaled_r5(x, freeze_r12_only_q(), nq_r5, r5_scale=1.0)
    s = float(r8_scale)
    if abs(s) < 1.0e-30:
        return J5
    return J5 + s * r8_analytic_jacobian(x, nq_all, rate_scale=1.0)


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _newton_at_r8_scale(
    nq_r5: NamedStageQ,
    nq_all: NamedStageQ,
    x: npt.NDArray[np.float64],
    r8_scale: float,
    max_iter: int,
) -> Any:
    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_r12_r5_plus_scaled_r8(z, nq_r5, nq_all, r8_scale=r8_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_r12_r5_plus_scaled_r8(z, nq_r5, nq_all, r8_scale=r8_scale)

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_newton_step,
    )


def r8_source_snapshot(
    x: npt.NDArray[np.float64],
    nq_r5: NamedStageQ,
    nq_all: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f5 = coupled_residual_named(z, nq_r5, rate_scale=1.0)
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    dF = both - f5
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_H2_d": float(dF[IDX["H2"]]),
        "dF_CO2_d": float(dF[IDX["CO2"]]),
        "dF_H2O_d": float(dF[IDX["H2O"]]),
        "dF_CO_b": float(dF[7 + IDX["CO"]]),
        "dF_H2_b": float(dF[7 + IDX["H2"]]),
        "dF_energy": float(dF[-1]),
        "N_d_CO": float(z[IDX["CO"]]),
        "N_d_H2": float(z[IDX["H2"]]),
        "N_d_H2O": float(z[IDX["H2O"]]),
        "N_d_O2": float(z[IDX["O2"]]),
        "T_K": float(z[-1]),
    }


def solve_named_r12_r5_then_r8(
    *,
    max_iter: int = 20,
    continuation: bool = True,
    r8_scale: float | None = None,
    r5: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E0+R12+R5 数学根做 R8 延拓。"""
    nq_r5 = freeze_r12_r5_q()
    nq_all = freeze_r12_r5_r8_q()
    if r5 is None:
        r5 = solve_named_r12_then_r5(max_iter=max_iter, apply_startwert=True)
    x_start = np.asarray(r5["x"], dtype=np.float64).copy()
    snap0 = r8_source_snapshot(x_start, nq_r5, nq_all)
    scales = (float(r8_scale),) if r8_scale is not None else (
        _R8_CONTINUATION if continuation else (1.0,)
    )
    n = n_unknowns(nq_all.stage)
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    for s in scales:
        last = _newton_at_r8_scale(nq_r5, nq_all, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_r12_r5_plus_scaled_r8(last.x, nq_r5, nq_all, r8_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        history.append(
            {
                "r8_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_CO": float(last.x[IDX["CO"]]),
                "N_d_CO2": float(last.x[IDX["CO2"]]),
                "N_d_H2O": float(last.x[IDX["H2O"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
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
        "T_r5_K": float(r5["T_K"]),
        "energy_W": float(energy_row_named(last.x, nq_all)),
        "f_name": "F_E0+R12+R5+R8",
        "enabled_rxns": ("R12", "R5", "R8"),
        "r5_physical_convergence": bool(r5["physical_convergence"]),
        "r8_at_r5_root": snap0,
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
    }
