"""F_E6 数学根上接入 R10（焦油氧化）。不改 ``rate_R10``，不经 ``Reactor``。

R1–R8/R12 与 DAEM、灰、``Q_W`` 保持完整。R10 用 ``r10_scale`` 延拓：
F(x;s) = F_{E6}(x) + s (F_{E6+R10}(x) - F_{E6}(x)).

Hamel Eq. 5.59：r = k10 exp(-E/RT) T P^{0.3} C_tar^{0.5} C_O2，走 ``k_hobbs``。
两相同一式。具名核没有焦油相间交换，C_tar 取全格摩尔分率。

E6 根上 ξ ~ 10^5 mol/s。从 E6 根做密网格延拓，不把 holdup 上的 O₂=0
当作 s=1 初值：气泡入口仍供氧，O₂=0 且 H₂/CO>0 会撞上 R12/R5 硬源。
定义域 x[:-1]≥0；越界试探缩短 α，残差本身不 clip。
细炭碰到 M=0 时只修正搜索方向（bound_feasible_direction），不改 F。
C_O2=0 但 C_tar>0 时保留 O₂ 右导数（氧一阶）。C_tar^{0.5} 在 0 处发散，不写有限右导数。

Source: Hamel (1999) Eq. 5.59–5.61, §2.1; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    R10_FUEL_TYPE,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
    r10_c_tar,
    r10_nu_major,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_e6 import (
    freeze_e6_q,
    jacobian_named_e6,
    residual_named_e6,
    solve_named_e6_from_e5,
)
from src.solvers.wirsum_1997_named_f_continuation import (
    AdvanceMode,
    bound_feasible_direction,
    certify_named_full_rate_endpoint,
    named_in_domain,
    named_linear_step,
    stage_passed,
)
from src.solvers.wirsum_1997_named_f_r12 import _conc_and_dN, apply_named_r12_startwert, named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_startwert_r10 import apply_named_r10_startwert
from src.solvers.wirsum_1997_startwert_r5 import apply_named_r5_startwert
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON | frozenset({"R10"})
_IDX_TAR1: int = 40
_IDX_TAR2: int = 41
_T_FD: float = 1.0  # [K]
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]

_R10_CONTINUATION: tuple[float, ...] = (
    1.0e-8, 1.0e-6, 3.0e-6, 1.0e-2, 0.1, 1.0
)
_R10_STARTWERT_CONTINUATION: tuple[float, ...] = (1.0e-4, 1.0e-2, 1.0e-1, 1.0)


def freeze_e6_r10_q(*, Q_W: float = 0.0) -> NamedStageQ:
    return freeze_stage_q(6, Q_W=float(Q_W), enabled_rxns=_ENABLED)


def residual_e6_plus_scaled_r10(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r10_scale: float,
) -> npt.NDArray[np.float64]:
    """F_E6 完整速率，R10 乘 ``r10_scale``。不 clip。"""
    z = np.asarray(x, dtype=np.float64)
    f6 = residual_named_e6(z, nq6)
    s = float(r10_scale)
    if abs(s) < 1.0e-30:
        return f6
    both = residual_named_e6(z, nq_all)
    return f6 + s * (both - f6)


def r10_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """R10 对 N_d/N_b/N_tar 的解析块。C_O2=0 时保留 O₂ 右导数（氧一阶）。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    if "R10" not in nq.enabled_rxns or int(nq.stage) < 6:
        return J
    s = float(rate_scale)
    if abs(s) < 1.0e-30:
        return J
    T = float(z[-1])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    n1 = max(float(z[_IDX_TAR1]), 0.0)
    n2 = max(float(z[_IDX_TAR2]), 0.0)
    n_tar = n1 + n2
    u = unpack_named(z, nq.stage)
    C_tar = r10_c_tar(u["N_d"], u["N_b"], n_tar, T, P)
    if C_tar <= 0.0:
        return J
    k_eff = float(rate_R10(T, C_tar, 1.0, P, R10_FUEL_TYPE)) / (C_tar ** 0.5)
    nu = r10_nu_major()
    sto = get_lumped_tar_stoichiometry("R10", R10_FUEL_TYPE)
    nu_t1 = float(sto.get("TAR1", 0.0))
    nu_t2 = float(sto.get("TAR2", 0.0))
    n_gas = float(np.sum(u["N_d"]) + np.sum(u["N_b"]))
    n_tot = n_gas + n_tar
    c_tot = float(P) / (Rg * float(T))
    if n_tot < _NTOT_FLOOR:
        return J
    dCtar_dnt = (n_gas / (n_tot * n_tot)) * c_tot
    dCtar_dNj = -(n_tar / (n_tot * n_tot)) * c_tot
    sqrt_ct = C_tar ** 0.5
    blocks = ((0, V_d), (7, V_b))
    for off, V in blocks:
        C, dCdN = _conc_and_dN(z[off:off + N_MAJ], T, P)
        c_o2 = float(C[IDX["O2"]])
        dr_dCo2 = k_eff * sqrt_ct
        dr_dCtar = (0.5 * k_eff / sqrt_ct) * max(c_o2, 0.0)
        dxi_this = (s * float(V)) * (
            dr_dCo2 * dCdN[IDX["O2"], :] + dr_dCtar * dCtar_dNj
        )
        dxi_other = (s * float(V)) * (dr_dCtar * dCtar_dNj)
        dxi_tar = (s * float(V)) * (dr_dCtar * dCtar_dnt)
        other = 7 if off == 0 else 0
        J[off:off + N_MAJ, off:off + N_MAJ] += nu[:, None] * dxi_this
        J[off:off + N_MAJ, other:other + N_MAJ] += nu[:, None] * dxi_other
        J[off:off + N_MAJ, _IDX_TAR1] += nu * dxi_tar
        J[off:off + N_MAJ, _IDX_TAR2] += nu * dxi_tar
        J[_IDX_TAR1, off:off + N_MAJ] -= nu_t1 * dxi_this
        J[_IDX_TAR2, off:off + N_MAJ] -= nu_t2 * dxi_this
        J[_IDX_TAR1, other:other + N_MAJ] -= nu_t1 * dxi_other
        J[_IDX_TAR2, other:other + N_MAJ] -= nu_t2 * dxi_other
        J[_IDX_TAR1, _IDX_TAR1] -= nu_t1 * dxi_tar
        J[_IDX_TAR1, _IDX_TAR2] -= nu_t1 * dxi_tar
        J[_IDX_TAR2, _IDX_TAR1] -= nu_t2 * dxi_tar
        J[_IDX_TAR2, _IDX_TAR2] -= nu_t2 * dxi_tar
    return J


def jacobian_e6_plus_scaled_r10(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r10_scale: float,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    J = jacobian_named_e6(z, nq6)
    s = float(r10_scale)
    if abs(s) >= 1.0e-30:
        J = J + r10_analytic_jacobian(z, nq_all, rate_scale=s)
        xp = z.copy()
        xp[-1] = float(z[-1]) + _T_FD
        dF = (
            residual_e6_plus_scaled_r10(xp, nq6, nq_all, r10_scale=s)
            - residual_e6_plus_scaled_r10(z, nq6, nq_all, r10_scale=s)
        ) / _T_FD
        J[:, -1] = dF
    return J


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return named_in_domain(z)


def _r10_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    dx, rank, cond, _diag = named_linear_step(
        J, F, x, extra_diag_species=("H2", "O2", "CO", "CH4"),
    )
    p, _changed = bound_feasible_direction(x, dx)
    return p, rank, cond


def r10_source_snapshot(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_all: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f6 = residual_named_e6(z, nq6)
    both = residual_named_e6(z, nq_all)
    dF = both - f6
    u = unpack_named(z, 6)
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_O2_d": float(dF[IDX["O2"]]),
        "dF_O2_b": float(dF[7 + IDX["O2"]]),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_H2_d": float(dF[IDX["H2"]]),
        "dF_tar1": float(dF[_IDX_TAR1]),
        "dF_tar2": float(dF[_IDX_TAR2]),
        "dF_energy": float(dF[-1]),
        "N_d_O2": float(u["N_d"][IDX["O2"]]),
        "N_b_O2": float(u["N_b"][IDX["O2"]]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "T_K": float(u["T"]),
    }


def _newton_at_r10_scale(
    nq6: NamedStageQ,
    nq_all: NamedStageQ,
    x: npt.NDArray[np.float64],
    r10_scale: float,
    max_iter: int,
) -> Any:
    x0 = np.asarray(x, dtype=np.float64)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_e6_plus_scaled_r10(z, nq6, nq_all, r10_scale=r10_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_e6_plus_scaled_r10(z, nq6, nq_all, r10_scale=r10_scale)

    return wirsum_1997_strict_newton(
        res, jac, x0, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=_r10_newton_step,
    )


def solve_named_e6_then_r10(
    *,
    Q_W: float = 0.0,
    max_iter: int = 20,
    continuation: bool = True,
    r10_scale: float | None = None,
    apply_startwert: bool = False,
    e6: dict[str, Any] | None = None,
    advance_mode: AdvanceMode = "heuristic_search",
    s_min: float = 1.0e-16,
) -> dict[str, Any]:
    """从 F_E6 根做 R10 延拓。默认启发式初值搜索，不是逐根跟踪。"""
    nq6 = freeze_e6_q(Q_W=Q_W)
    nq_all = freeze_e6_r10_q(Q_W=Q_W)
    if e6 is None:
        e6 = solve_named_e6_from_e5(Q_W=Q_W, max_iter=max_iter)
    x_start = np.asarray(e6["x"], dtype=np.float64).copy()
    snap_root = r10_source_snapshot(x_start, nq6, nq_all)
    sw10 = None
    if apply_startwert:
        packed = apply_named_r10_startwert(x_start)
        x_start = packed["x"]
        sw10 = {"dense": packed["dense"], "bubble": packed["bubble"]}
        x_start = apply_named_r12_startwert(x_start)["x"]
        x_start = apply_named_r5_startwert(x_start)["x"]
    snap0 = r10_source_snapshot(x_start, nq6, nq_all)
    if r10_scale is not None:
        planned = [float(r10_scale)]
    elif continuation:
        planned = list(
            _R10_STARTWERT_CONTINUATION if apply_startwert else _R10_CONTINUATION
        )
    else:
        planned = [1.0]
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    s_ok = 0.0
    branch_tracked = True
    i = 0
    while i < len(planned):
        s = float(planned[i])
        last = _newton_at_r10_scale(nq6, nq_all, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_e6_plus_scaled_r10(last.x, nq6, nq_all, r10_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        passed = stage_passed(
            stop_s, last.x, advance_mode=advance_mode,
            accepted_alphas=list(last.accepted_alphas),
        )
        history.append(
            {
                "r10_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "solver_stop_reason": stop_s["failure_reason"],
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "stage_passed": bool(passed),
                "N_d_O2": float(last.x[IDX["O2"]]),
                "N_b_O2": float(last.x[7 + IDX["O2"]]),
                "N_tar1": float(last.x[_IDX_TAR1]),
                "N_tar2": float(last.x[_IDX_TAR2]),
            }
        )
        if passed:
            x_cur = np.asarray(last.x, dtype=np.float64)
            s_ok = s
            i += 1
            continue
        branch_tracked = False
        if advance_mode == "root_correcting":
            s_mid = 0.5 * (s_ok + s)
            if abs(s - s_ok) < float(s_min) or abs(s_mid - s) < 1.0e-18:
                break
            planned.insert(i, s_mid)
            continue
        i += 1

    assert last is not None
    s_last = float(history[-1]["r10_scale"]) if history else 1.0
    stop = classify_newton_stop(
        last, model_is_full_rate=abs(s_last - 1.0) < 1.0e-15
    )
    info = named_r12_scaled_rank(
        jacobian_e6_plus_scaled_r10(last.x, nq6, nq_all, r10_scale=1.0),
        last.F, last.x,
    )
    nqP = outer_P_named(last.x, nq_all)
    u = unpack_named(last.x, 6)
    cert = certify_named_full_rate_endpoint(last.x, nq_all, coupling="composite_P")
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": float(last.x[-1]),
        "T_e6_K": float(e6["T_K"]),
        "M_ash": float(u["M_ash"]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "N_d_O2": float(u["N_d"][IDX["O2"]]),
        "N_b_O2": float(u["N_b"][IDX["O2"]]),
        "Q_W": float(nq_all.Q_W),
        "Q_W_is_fitted_to_1100K": False,
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": "F_E6+R10",
        "n_unknowns": n_unknowns(6),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "e6_physical_convergence": bool(e6["physical_convergence"]),
        "r10_at_e6_root": snap_root,
        "r10_after_startwert": snap0,
        "startwert_r10": sw10,
        "apply_startwert": bool(apply_startwert),
        "advance_mode": str(advance_mode),
        "branch_tracked": bool(
            branch_tracked
            and bool(history)
            and abs(float(history[-1]["r10_scale"]) - 1.0) < 1.0e-15
        ),
        "solver_stop_reason": stop["failure_reason"],
        "endpoint_cert": {k: v for k, v in cert.items() if k != "F"},
        "endpoint_certified": bool(cert["endpoint_certified"]),
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
        "production_hamel_outer_refresh": False,
        "r10_tar_oxidation": True,
        "r11_tar_conversion": False,
        "x_nonneg_projection": False,
        "coupling_mode": "composite_P",
        "bound_feasible_direction": True,
    }
