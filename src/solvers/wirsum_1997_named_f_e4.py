"""F_E4：在 F_E3 数学根上打开 DAEM ``X(T)``。不经 ``Reactor``，不切生产。

未知量仍是 41。残差走具名外层 ``outer_P_named``：``n_py(T)``、焦油由
DAEM ``X(T)`` 与热解分配更新。这不是生产 ``hamel_outer_refresh``。

F_E3 冻结热解在 geo 1100 K、``X=1``。E3 根上 DAEM 已饱和（X≈1），
主要变化是 Hegermann 产物随 T 重分配，嵌入点能量行约数百瓦。
``n_py`` 增量用 ``F_E3 + s (H_E4 - F_E3)`` 延拓。密相出口残氧不够
烧掉 Hegermann 重分配的 H2；把 Δs 的 H2 增量按反向 WGSR
（H2+CO2→CO+H2O）接到 ``N_d``，不写入 H2/CH4。不改速率式。

雅可比：N/M_ash 用 F_E3 块；T 列对复合残差 ``H(x)=F(x;P_DAEM(x))`` 差商。
默认耦合是 ``composite_P``，不是生产内外层冻结 q 再刷新。``frozen_q`` 接口保留作对照。

Source: Hamel (1999) Eq. 4.10–4.12, §2.1 Startwertwahl; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    daem_X_at_T,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_continuation import (
    CouplingMode,
    bound_feasible_direction,
    shrink_scale,
)
from src.solvers.wirsum_1997_named_f_e3 import (
    freeze_e3_q,
    jacobian_named_e3,
    solve_named_e3_from_r234,
)
from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import named_r12_r5_r8_r6_newton_step
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_named_f_registry import f_name
from src.solvers.wirsum_1997_stiffness_full_gas import IDX
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON
_T_FD: float = 1.0  # [K]
_T_FLOOR: float = 1.0  # [K]
_S_MIN: float = 2.5e-4  # [-]
_TOL_MERIT: float = 1.0e-8  # 与 physical L2 同量级，少做无效二次迭代
# 0.20 仍能从 E3 一跳；0.40 已进反向 WGSR 支；其后 Δs≤0.08。
_E4_DAEM_CONTINUATION: tuple[float, ...] = (
    0.20, 0.40, 0.48, 0.56, 0.64, 0.72, 0.80, 0.88, 0.96, 1.0
)


def freeze_e4_q() -> NamedStageQ:
    return freeze_stage_q(4, enabled_rxns=_ENABLED)


def residual_named_e4(
    x: npt.NDArray[np.float64],
    nq4: NamedStageQ,
    *,
    coupling: CouplingMode = "composite_P",
) -> npt.NDArray[np.float64]:
    """默认 ``H(x)=F(x;P_DAEM(x))``。``frozen_q`` 不调用 ``outer_P_named``。"""
    z = np.asarray(x, dtype=np.float64)
    nq = outer_P_named(z, nq4) if coupling == "composite_P" else nq4
    return coupled_residual_named(z, nq, rate_scale=1.0)


def jacobian_named_e4(
    x: npt.NDArray[np.float64],
    nq4: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """N/M 用 E3 块（n_py 对 N 无显式依赖）；T 列含 DAEM。"""
    z = np.asarray(x, dtype=np.float64)
    nq = outer_P_named(z, nq4)
    J = jacobian_named_e3(z, nq)
    xp = z.copy()
    xp[-1] = float(z[-1]) + _T_FD
    dF = (residual_named_e4(xp, nq4) - residual_named_e4(z, nq4)) / _T_FD
    J[:, -1] = dF
    return J


def residual_e3_plus_scaled_e4(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    nq4: NamedStageQ,
    *,
    daem_scale: float,
    coupling: CouplingMode = "composite_P",
) -> npt.NDArray[np.float64]:
    """``F_E3 + s (H_E4 - F_E3)``。``s=0`` 退回冻结热解，``s=1`` 为 ``H(x;P_DAEM)``。"""
    z = np.asarray(x, dtype=np.float64)
    f3 = coupled_residual_named(z, nq3, rate_scale=1.0)
    s = float(daem_scale)
    if abs(s) < 1.0e-30:
        return f3
    f4 = residual_named_e4(z, nq4, coupling=coupling)
    if abs(s - 1.0) < 1.0e-15:
        return f4
    return f3 + s * (f4 - f3)


def _fd_jacobian(
    residual_fn: Any,
    x: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """延拓中间步用差商 J。n_py 对 N 无显式依赖，避免 1 s 级解析块。"""
    z = np.asarray(x, dtype=np.float64)
    f0 = residual_fn(z)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        xp = z.copy()
        if i == n - 1:
            d = _T_FD
        else:
            d = max(1.0e-12, 1.0e-6 * abs(float(z[i])))
        xp[i] = float(xp[i]) + d
        J[:, i] = (residual_fn(xp) - f0) / d
    return J


def _sync_cyclone_gas(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    out = np.asarray(z, dtype=np.float64)
    ntot = out[0:7] + out[7:14]
    out[16:23] = ntot
    out[23:30] = ntot
    out[30:37] = ntot
    return out


def apply_pyro_reverse_wgsr_startwert_e4(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    nq4: NamedStageQ,
    *,
    daem_scale_step: float,
) -> dict[str, Any]:
    """把 ``Δs`` 的热解 H2 增量按反向 WGSR 接到密相 CO/H2O。

    ξ8 = Δs · max(Δn_py,H2, 0)。CO += Δs Δn_CO + ξ8，CO2 += Δs Δn_CO2 − ξ8，
    H2O += Δs Δn_H2O + ξ8。不写 H2（R12/R8 会炸），不减 CH4（``N_d`` 远小于
    |Δn_CH4|）。CO2 不够时削 ξ8。同步 FB/旋风气。不改 F(x)。

    Source: Hamel (1999) §2.1 Startwertwahl; Eq. 5.45–5.46
    """
    z = np.asarray(x, dtype=np.float64).copy()
    ds = float(daem_scale_step)
    n_py0 = np.asarray(nq3.energy.species.n_py_d, dtype=np.float64)
    n_py1 = np.asarray(outer_P_named(z, nq4).energy.species.n_py_d, dtype=np.float64)
    dn = n_py1 - n_py0
    n_d = z[0:7].copy()
    xi8 = ds * max(float(dn[IDX["H2"]]), 0.0)  # [mol/s]
    d = np.zeros(7, dtype=np.float64)
    d[IDX["CO"]] = ds * float(dn[IDX["CO"]]) + xi8
    d[IDX["CO2"]] = ds * float(dn[IDX["CO2"]]) - xi8
    d[IDX["H2O"]] = ds * float(dn[IDX["H2O"]]) + xi8
    n_co2 = float(n_d[IDX["CO2"]]) + float(d[IDX["CO2"]])
    if n_co2 < 0.0:
        extra = -n_co2  # [mol/s]
        d[IDX["CO2"]] += extra
        d[IDX["CO"]] -= extra
        d[IDX["H2O"]] -= extra
        xi8 = max(xi8 - extra, 0.0)
    z[0:7] = np.maximum(n_d + d, 0.0)
    z = _sync_cyclone_gas(z)
    return {
        "x": z,
        "xi8_mol_s": float(xi8),
        "daem_scale_step": ds,
        "dn_py": dn,
    }


def apply_energy_t_startwert_e4(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    nq4: NamedStageQ,
    *,
    daem_scale: float,
) -> dict[str, Any]:
    """打开 DAEM 焓后的温度 Startwert：``dT = -E / (∂E/∂T)``。

    Source: Hamel (1999) §2.1 Startwertwahl; Wirsum (1997) §2.3.3
    """
    z = np.asarray(x, dtype=np.float64).copy()
    s = float(daem_scale)
    F = residual_e3_plus_scaled_e4(z, nq3, nq4, daem_scale=s)
    xp = z.copy()
    xp[-1] = float(z[-1]) + _T_FD
    Fp = residual_e3_plus_scaled_e4(xp, nq3, nq4, daem_scale=s)
    je = float(Fp[-1] - F[-1]) / _T_FD
    e0 = float(F[-1])
    dT = 0.0  # [K]
    applied = False
    if abs(je) >= 1.0e-30:
        dT = -e0 / je
        t1 = float(z[-1]) + dT
        if t1 > _T_FLOOR:
            z[-1] = t1
            applied = True
    F1 = residual_e3_plus_scaled_e4(z, nq3, nq4, daem_scale=s)
    return {
        "x": z,
        "applied": bool(applied),
        "dT_K": float(dT),
        "E_before_W": e0,
        "E_after_W": float(F1[-1]),
        "dE_dT_W_per_K": je,
        "daem_scale": s,
    }


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _e4_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    dx, rank, cond = named_r12_r5_r8_r6_newton_step(J, F, x)
    p, _changed = bound_feasible_direction(x, dx)
    return p, rank, cond


def _newton_at_daem_scale(
    nq3: NamedStageQ,
    nq4: NamedStageQ,
    x: npt.NDArray[np.float64],
    daem_scale: float,
    max_iter: int,
    *,
    live_p: bool = False,
    analytic_j: bool = False,
) -> Any:
    s = float(daem_scale)
    nq_eval = nq4 if live_p else outer_P_named(x, nq4)
    coupling: CouplingMode = "composite_P" if live_p else "frozen_q"

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_e3_plus_scaled_e4(
            z, nq3, nq_eval, daem_scale=s, coupling=coupling,
        )

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        if live_p and abs(s - 1.0) < 1.0e-15:
            return jacobian_named_e4(z, nq4)
        if analytic_j:
            J = jacobian_named_e3(z, nq3)
            xp = z.copy()
            xp[-1] = float(z[-1]) + _T_FD
            J[:, -1] = (res(xp) - res(z)) / _T_FD
            return J
        return _fd_jacobian(res, z)

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=_e4_newton_step, tol_merit=_TOL_MERIT,
    )


def e4_source_snapshot(
    x: npt.NDArray[np.float64],
    nq3: NamedStageQ,
    nq4: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f3 = coupled_residual_named(z, nq3, rate_scale=1.0)
    nqP = outer_P_named(z, nq4)
    f4 = coupled_residual_named(z, nqP, rate_scale=1.0)
    dF = f4 - f3
    t = float(z[-1])
    n_py0 = np.asarray(nq3.energy.species.n_py_d, dtype=np.float64)
    n_py1 = np.asarray(nqP.energy.species.n_py_d, dtype=np.float64)
    return {
        "F3_l2": float(np.linalg.norm(f3)),
        "F4_l2": float(np.linalg.norm(f4)),
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_energy": float(dF[-1]),
        "dF_H2_d": float(dF[IDX["H2"]]),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_CH4_d": float(dF[IDX["CH4"]]),
        "n_py_sum_frozen": float(np.sum(n_py0)),
        "n_py_sum_daem": float(np.sum(n_py1)),
        "X_vm": float(daem_X_at_T(t, nq4.tau_daem_s)),
        "T_K": t,
    }


def solve_named_e4_from_e3(
    *,
    max_iter: int = 40,
    e3: dict[str, Any] | None = None,
    continuation: bool = True,
    daem_scale: float | None = None,
) -> dict[str, Any]:
    """从 F_E3 根做 F_E4。反向 WGSR 接热解 H2 增量，再延拓 DAEM。"""
    nq3 = freeze_e3_q()
    nq4 = freeze_e4_q()
    if e3 is None:
        e3 = solve_named_e3_from_r234(max_iter=max_iter)
    x0 = np.asarray(e3["x"], dtype=np.float64).copy()
    snap0 = e4_source_snapshot(x0, nq3, nq4)
    sw0 = apply_energy_t_startwert_e4(x0, nq3, nq4, daem_scale=1.0)
    history: list[dict[str, Any]] = []
    x_cur = x0
    last = None
    s_ok = 0.0
    if daem_scale is not None:
        planned = [float(daem_scale)]
    elif continuation:
        planned = list(_E4_DAEM_CONTINUATION)
    else:
        planned = [1.0]

    def _record(s: float, xi8: float, last: Any, stop_s: dict[str, Any]) -> None:
        history.append(
            {
                "daem_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "xi8_mol_s": float(xi8),
                "rank_last": int(last.rank_last),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "physical_convergence": bool(stop_s["physical_convergence"]),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
                "N_d_CO": float(last.x[IDX["CO"]]),
            }
        )

    def _advance(s: float) -> tuple[Any, bool]:
        ds = float(s) - float(s_ok)
        pyro = apply_pyro_reverse_wgsr_startwert_e4(
            x_cur, nq3, nq4, daem_scale_step=ds,
        )
        analytic = s >= 0.08 - 1.0e-15
        nxt = _newton_at_daem_scale(
            nq3, nq4, pyro["x"], s, max_iter,
            live_p=analytic, analytic_j=analytic,
        )
        stop_s = classify_newton_stop(nxt, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        _record(s, float(pyro["xi8_mol_s"]), nxt, stop_s)
        return nxt, float(stop_s["residual_l2"]) < 1.0e-6

    i = 0
    while i < len(planned):
        s = float(planned[i])
        last, passed = _advance(s)
        if passed:
            x_cur = np.asarray(last.x, dtype=np.float64)
            s_ok = s
            i += 1
            continue
        if abs(s - s_ok) < float(_S_MIN):
            break
        s_mid = shrink_scale(s_ok, s)
        if abs(s_mid - s) < 1.0e-18:
            break
        planned.insert(i, s_mid)

    assert last is not None
    F_phys = residual_named_e4(last.x, nq4)
    last.F = np.asarray(F_phys, dtype=np.float64)
    stop = classify_newton_stop(last, model_is_full_rate=True)
    info = named_r12_scaled_rank(jacobian_named_e4(last.x, nq4), last.F, last.x)
    nqP = outer_P_named(last.x, nq4)
    u = unpack_named(last.x, 4)
    t = float(last.x[-1])
    n_py = np.asarray(nqP.energy.species.n_py_d, dtype=np.float64)
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": t,
        "T_e3_K": float(e3["T_K"]),
        "M_ash": float(u["M_ash"]),
        "X_vm": float(daem_X_at_T(t, nq4.tau_daem_s)),
        "n_py_d": n_py,
        "n_py_sum": float(np.sum(n_py)),
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": f_name(4),
        "n_unknowns": n_unknowns(4),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "e3_physical_convergence": bool(e3["physical_convergence"]),
        "at_e3_root": snap0,
        "energy_t_startwert": {k: v for k, v in sw0.items() if k != "x"},
        "continuation": history,
        "pyro_reverse_wgsr_startwert": True,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "bound_feasible_direction": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
        "production_hamel_outer_refresh": False,
        "tau_daem_s": float(nq4.tau_daem_s),
        "coupling_mode": "composite_P",
    }
