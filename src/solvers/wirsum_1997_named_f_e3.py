"""F_E3：在均相+R2–R4 数学根上打开灰库存。不经 ``Reactor``，不切生产。

40 维根嵌入 41 维。``M_ash`` 用 Eq. 2.6 闭合 ``m_feed = f_fly K_e M`` 做 Startwert。
灰不参与 R1–R4。耦合只经密相 R7/R8（``ρ_ash = M_ash / V_d``）和飞灰焓。
冻结 ``rho_ash=75`` 与闭合库存差两个数量级；不得把 75 kg/m³ 当成根。

打开飞灰焓后 ``E ≈ -H_ash``。能量行用 ``dT = -E / (∂E/∂T)`` 做 Startwertwahl。
R7/R8 的灰密度从冻结 75 kg/m³ 延拓到 ``M_ash/V_d``，不改速率式。

旋风 K 与 ``hamel_outer_refresh`` 仍冻结。

Source: Hamel (1999) Eq. 2.6, §2.1 Startwertwahl; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import SOLID_H_MODE_INTEGRAL, h_solid_mass
from src.kinetics.gas_reactions import rate_R7, rate_R8
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import IDX_COARSE
from src.solvers.wirsum_1997_named_f import (
    N_SPECIES_BLOCK,
    RHO_ASH_Q_KG_M3,
    NamedStageQ,
    ash_feed_kg_s,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
    homog_sources,
    n_unknowns,
    pack_named,
    q_at_T,
    rho_ash_used,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_continuation import (
    bound_feasible_direction,
    shrink_scale,
)
from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import NU_R8_MAJOR, r8_analytic_jacobian
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import NU_R7_MAJOR, r7_analytic_jacobian
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import named_r12_r5_r8_r6_newton_step
from src.solvers.wirsum_1997_named_f_r234 import (
    HOMOG_ON,
    R234_ON,
    freeze_homog_r234_q,
    jacobian_homog_plus_scaled_r234,
    solve_named_homog_then_r234,
)
from src.solvers.wirsum_1997_named_f_registry import f_name
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ASH_HF: float = -1.2e6  # [J/kg] 与 energy_row_named 相同
_T_FD: float = 1.0  # [K] 灰焓温度差商
_IDX_ASH: int = 39
_ENABLED: frozenset[str] = HOMOG_ON | R234_ON
_T_FLOOR: float = 1.0  # [K]
# ρ：s=0 冻结 q 灰密度，s=1 为 M_ash/V_d。中间步校正 R8 组成，不是改速率。
_E3_RHO_CONTINUATION: tuple[float, ...] = (
    0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0
)
_S_MIN: float = 1.0e-3  # [-] 灰密度延拓最小步


def freeze_e3_q() -> NamedStageQ:
    return freeze_stage_q(3, enabled_rxns=_ENABLED)


def mash_eq26_kg(nq: NamedStageQ) -> float:
    """Eq. 2.6 灰库存闭合 [kg]：``M = m_feed / (f_fly K_e)``。"""
    k = nq.energy.species.k
    k_sum = float(k.K_auf) + float(k.K_ab)
    fly = float(k.K_auf) / max(k_sum, 1.0e-30)
    k_e = float(nq.energy.species.K_eject[IDX_COARSE])
    return float(ash_feed_kg_s(nq)) / max(fly * k_e, 1.0e-30)


def _ash_fly_coeff(nq: NamedStageQ) -> float:
    """``f_fly K_e`` [1/s]。"""
    k = nq.energy.species.k
    k_sum = float(k.K_auf) + float(k.K_ab)
    fly = float(k.K_auf) / max(k_sum, 1.0e-30)
    return fly * float(nq.energy.species.K_eject[IDX_COARSE])


def embed_e2_root_in_e3(
    x40: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    *,
    M_ash: float | None = None,
) -> npt.NDArray[np.float64]:
    """40 维根嵌入 F_E3。默认 Eq. 2.6 闭合库存。"""
    z = np.asarray(x40, dtype=np.float64)
    m = float(mash_eq26_kg(nq3) if M_ash is None else M_ash)
    return pack_named(z[:39], float(z[-1]), 3, M_ash=m)


def _x40_from_e3(x41: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], float]:
    z = np.asarray(x41, dtype=np.float64)
    return np.concatenate([z[:39], z[-1:]]), float(z[_IDX_ASH])


def _lift_j40(J40: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    j = np.asarray(J40, dtype=np.float64)
    J = np.zeros((41, 41), dtype=np.float64)
    J[:39, :39] = j[:39, :39]
    J[:39, -1] = j[:39, -1]
    J[-1, :39] = j[-1, :39]
    J[-1, -1] = j[-1, -1]
    return J


def _h_ash(T: float) -> float:
    return float(
        h_solid_mass("ash", float(T), h_f=_ASH_HF, mode=SOLID_H_MODE_INTEGRAL)
    )


def _rho_homog(z: npt.NDArray[np.float64], nq3: NamedStageQ, rho_scale: float) -> float:
    """密相 R7/R8 灰密度 [kg/m³]。``s=0`` 为冻结 q，``s=1`` 为 ``M_ash/V_d``。"""
    V_d = float(nq3.energy.species.geo_ex.V_d)
    rho_true = float(z[_IDX_ASH]) / max(V_d, 1.0e-12)
    s = float(rho_scale)
    return (1.0 - s) * float(RHO_ASH_Q_KG_M3) + s * rho_true


def residual_named_e3(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    *,
    rho_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    """物理残差上把均相源从 ``M/V`` 换到 ``ρ(s)``。``s=1`` 即 ``coupled_residual_named``。"""
    z = np.asarray(x, dtype=np.float64)
    F = coupled_residual_named(z, nq3, rate_scale=1.0)
    s = float(rho_scale)
    if abs(s - 1.0) < 1.0e-15:
        return F
    V_d = float(nq3.energy.species.geo_ex.V_d)
    rho_true = float(z[_IDX_ASH]) / max(V_d, 1.0e-12)
    rho_h = _rho_homog(z, nq3, s)
    T = float(z[-1])
    qT = q_at_T(nq3.energy.species, T)
    sp = z[:N_SPECIES_BLOCK]
    on = nq3.enabled_rxns
    sd_t, sb_t = homog_sources(sp, qT, rate_scale=1.0, rho_ash=rho_true, enabled=on)
    sd_h, sb_h = homog_sources(sp, qT, rate_scale=1.0, rho_ash=rho_h, enabled=on)
    out = F.copy()
    out[0:7] = out[0:7] + sd_h - sd_t
    out[7:14] = out[7:14] + sb_h - sb_t
    return out


def jacobian_named_e3(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    *,
    rho_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    """40 维均相+R234 雅可比上嵌入灰行/列。R7/R8 用 ``ρ(s)``。

    ``s=1`` 时 ``ρ = M_ash/V_d``，与物理残差一致。
    """
    z = np.asarray(x, dtype=np.float64)
    x40, m_ash = _x40_from_e3(z)
    nq2 = freeze_homog_r234_q()
    s = float(rho_scale)
    rho_h = _rho_homog(z, nq3, s)
    J40 = jacobian_homog_plus_scaled_r234(x40, nq2, nq2, r234_scale=1.0)
    nq_rho = replace(nq2, rho_ash=rho_h)
    J40 = (
        J40
        - r8_analytic_jacobian(x40, nq2, rate_scale=1.0)
        + r8_analytic_jacobian(x40, nq_rho, rate_scale=1.0)
        - r7_analytic_jacobian(x40, nq2, rate_scale=1.0)
        + r7_analytic_jacobian(x40, nq_rho, rate_scale=1.0)
    )
    J = _lift_j40(J40)
    coeff = _ash_fly_coeff(nq3)
    J[_IDX_ASH, _IDX_ASH] = -coeff
    T = float(z[-1])
    h_a = _h_ash(T)
    J[-1, _IDX_ASH] -= coeff * h_a
    h_a2 = _h_ash(T + _T_FD)
    J[-1, -1] -= coeff * float(m_ash) * ((h_a2 - h_a) / _T_FD)
    y = z[0:N_MAJ]
    ntot = float(np.sum(y))
    if abs(ntot) >= 1.0e-30:
        y = y / ntot
    else:
        y = np.zeros(N_MAJ)
    P = float(nq3.energy.species.geo_ex.P)
    r_unit = float(
        rate_R8(
            T, P,
            float(y[IDX["CO"]]), float(y[IDX["H2O"]]),
            float(y[IDX["CO2"]]), float(y[IDX["H2"]]),
            1.0,
        )
    )
    J[0:N_MAJ, _IDX_ASH] += s * NU_R8_MAJOR * r_unit
    C_d = (float(P) / (Rg * max(T, 1.0))) * y
    r7_unit = float(
        rate_R7(
            T,
            float(C_d[IDX["CH4"]]), float(C_d[IDX["H2O"]]),
            float(C_d[IDX["CO"]]), float(C_d[IDX["H2"]]),
            rho_ash_kg_m3=1.0,
        )
    )
    J[0:N_MAJ, _IDX_ASH] += s * NU_R7_MAJOR * r7_unit
    return J


def apply_energy_t_startwert(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
) -> dict[str, Any]:
    """打开飞灰焓后的温度 Startwert：``dT = -E / (∂E/∂T)``。

    在 ``ρ=75`` 残差上取步，不改物种。灰行已由 Eq. 2.6 闭合。

    Source: Hamel (1999) §2.1 Startwertwahl; Wirsum (1997) §2.3.3
    """
    z = np.asarray(x, dtype=np.float64).copy()
    F = residual_named_e3(z, nq3, rho_scale=0.0)
    J = jacobian_named_e3(z, nq3, rho_scale=0.0)
    e0 = float(F[-1])
    je = float(J[-1, -1])
    dT = 0.0  # [K]
    applied = False
    if abs(je) >= 1.0e-30:
        dT = -e0 / je
        t1 = float(z[-1]) + dT
        if t1 > _T_FLOOR:
            z[-1] = t1
            applied = True
    F1 = residual_named_e3(z, nq3, rho_scale=0.0)
    return {
        "x": z,
        "applied": bool(applied),
        "dT_K": float(dT),
        "E_before_W": e0,
        "E_after_W": float(F1[-1]),
        "dE_dT_W_per_K": je,
    }


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _e3_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    dx, rank, cond = named_r12_r5_r8_r6_newton_step(J, F, x)
    p, _changed = bound_feasible_direction(x, dx)
    return p, rank, cond


def _newton_at_rho_scale(
    nq3: NamedStageQ,
    x: npt.NDArray[np.float64],
    rho_scale: float,
    max_iter: int,
) -> Any:
    s = float(rho_scale)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_named_e3(z, nq3, rho_scale=s)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_named_e3(z, nq3, rho_scale=s)

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=_e3_newton_step,
    )


def e3_source_snapshot(
    x40: npt.NDArray[np.float64],
    nq2: NamedStageQ,
    nq3: NamedStageQ,
    x41: npt.NDArray[np.float64],
) -> dict[str, Any]:
    f2 = coupled_residual_named(x40, nq2, rate_scale=1.0)
    f3 = coupled_residual_named(x41, nq3, rate_scale=1.0)
    u = unpack_named(x41, 3)
    return {
        "F40_l2": float(np.linalg.norm(f2)),
        "F41_l2": float(np.linalg.norm(f3)),
        "F_ash": float(f3[_IDX_ASH]),
        "F_energy": float(f3[-1]),
        "F_gas_l2": float(np.linalg.norm(f3[:39])),
        "M_ash": float(u["M_ash"]),
        "rho_ash": float(rho_ash_used(nq3, u)),
        "rho_q": float(nq2.rho_ash),
        "T_K": float(x41[-1]),
    }


def solve_named_e3_from_r234(
    *,
    max_iter: int = 20,
    r234: dict[str, Any] | None = None,
    continuation: bool = True,
    rho_scale: float | None = None,
) -> dict[str, Any]:
    """从均相+R2–R4 根做 F_E3。``M_ash`` 用 Eq. 2.6 闭合。

    先做能量 T Startwert，再把密相灰密度从冻结 75 kg/m³ 延拓到 ``M/V``。
    """
    nq2 = freeze_homog_r234_q()
    nq3 = freeze_e3_q()
    if r234 is None:
        r234 = solve_named_homog_then_r234(max_iter=max_iter)
    x40 = np.asarray(r234["x"], dtype=np.float64)
    x0 = embed_e2_root_in_e3(x40, nq3)
    snap0 = e3_source_snapshot(x40, nq2, nq3, x0)
    sw_e = apply_energy_t_startwert(x0, nq3)
    x_start = np.asarray(sw_e["x"], dtype=np.float64)
    if rho_scale is not None:
        planned = [float(rho_scale)]
    elif continuation:
        planned = list(_E3_RHO_CONTINUATION)
    else:
        planned = [1.0]
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    s_ok = 0.0
    i = 0
    while i < len(planned):
        s = float(planned[i])
        last = _newton_at_rho_scale(nq3, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_named_e3(last.x, nq3, rho_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        passed = bool(stop_s["physical_convergence"]) if abs(s - 1.0) < 1.0e-15 else True
        history.append(
            {
                "rho_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "physical_convergence": bool(stop_s["physical_convergence"]),
                "N_d_CO": float(last.x[IDX["CO"]]),
                "N_d_H2": float(last.x[IDX["H2"]]),
            }
        )
        if abs(s - 1.0) >= 1.0e-15:
            x_cur = np.asarray(last.x, dtype=np.float64)
            s_ok = s
            i += 1
            continue
        if passed:
            x_cur = np.asarray(last.x, dtype=np.float64)
            s_ok = s
            i += 1
            continue
        s_mid = shrink_scale(s_ok, s)
        if abs(s - s_ok) < float(_S_MIN) or abs(s_mid - s) < 1.0e-18:
            break
        planned.insert(i, s_mid)

    assert last is not None
    F_phys = coupled_residual_named(last.x, nq3, rate_scale=1.0)
    last.F = np.asarray(F_phys, dtype=np.float64)
    stop = classify_newton_stop(last, model_is_full_rate=True)
    info = named_r12_scaled_rank(jacobian_named_e3(last.x, nq3), last.F, last.x)
    u = unpack_named(last.x, 3)
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": float(last.x[-1]),
        "T_r234_K": float(r234["T_K"]),
        "M_ash": float(u["M_ash"]),
        "M_ash_eq26": float(mash_eq26_kg(nq3)),
        "rho_ash": float(rho_ash_used(nq3, u)),
        "energy_W": float(energy_row_named(last.x, nq3)),
        "f_name": f_name(3),
        "n_unknowns": n_unknowns(3),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "r234_physical_convergence": bool(r234["physical_convergence"]),
        "at_embed": snap0,
        "energy_t_startwert": {k: v for k, v in sw_e.items() if k != "x"},
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "bound_feasible_direction": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
        "production_hamel_outer_refresh": False,
    }
