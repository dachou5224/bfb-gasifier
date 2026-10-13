"""F_E0+R12+R5+R8 数学根上接入 R6（CH4 氧化）。不改 ``rate_R6``，不经 ``Reactor``。

R12/R5/R8 保持完整速率。R6 用 ``r6_scale`` 延拓：
F(x;s) = F_{E0+R12+R5+R8}(x) + s (F_{E0+R12+R5+R8+R6}(x) - F_{E0+R12+R5+R8}(x)).

仓库 R6 = thesis R7。Hamel Eq. 5.40：r = k_ch,R7 C_CH4 C_O2；
k = (A/T) exp(-15700/T)，走 ``k_jensen_r7``。两相同一式。

R8 根上 CH4∩O2 给出 10²–10³ mol/s 源。默认同相 R6 Startwert；
R6 产 CO，已开的 R5 会再硬，所以随后再跑 R5 Startwert。
C_CH4=0 但 O2>0 时保留 CH4 右导数（级数 1）。

Source: Hamel (1999) Eq. 5.40–5.41, §2.1; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.kinetics.arrhenius import k_jensen_r7
from src.kinetics.gas_reactions import R6_E_Rg, R6_k0, rate_R6
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
)
from src.solvers.wirsum_1997_named_f_r12 import (
    _conc_and_dN,
    named_r12_scaled_rank,
)
from src.solvers.wirsum_1997_named_f_r12_r5 import freeze_r12_r5_q
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import (
    freeze_r12_r5_r8_q,
    jacobian_r12_r5_plus_scaled_r8,
    solve_named_r12_r5_then_r8,
)
from src.solvers.wirsum_1997_startwert_r5 import apply_named_r5_startwert
from src.solvers.wirsum_1997_startwert_r6 import apply_named_r6_startwert
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ, NU_R6
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

# CH4 + 1.5 O2 → CO + 2 H2O。Ref: Hamel (1999) Eq. 5.40。
NU_R6_MAJOR: npt.NDArray[np.float64] = np.asarray(NU_R6, dtype=np.float64)
_T_FLOOR: float = 300.0  # [K] 与 rate_R6 相同

_R6_CONTINUATION: tuple[float, ...] = (
    0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0
)
# Startwert 后 s=0 会按 F_E0+R12+R5+R8 把 CH4/O2 再混回去。
_R6_STARTWERT_CONTINUATION: tuple[float, ...] = (1.0e-4, 1.0e-2, 1.0e-1, 1.0)


def freeze_r12_r5_r8_r6_q() -> NamedStageQ:
    return freeze_stage_q(2, enabled_rxns=frozenset({"R12", "R5", "R8", "R6"}))


def residual_r12_r5_r8_plus_scaled_r6(
    x: npt.NDArray[np.float64],
    nq_r8: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r6_scale: float,
) -> npt.NDArray[np.float64]:
    """R12+R5+R8 完整速率，R6 乘 ``r6_scale``。"""
    z = np.asarray(x, dtype=np.float64)
    f8 = coupled_residual_named(z, nq_r8, rate_scale=1.0)
    s = float(r6_scale)
    if abs(s) < 1.0e-30:
        return f8
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    return f8 + s * (both - f8)


def _r6_r_and_dC(
    T: float, C: npt.NDArray[np.float64]
) -> tuple[float, npt.NDArray[np.float64]]:
    """与 ``rate_R6`` 同式。

    C_CH4=0 但 O2>0 时，CH4 右导数为 k C_O2（级数 1），不得取 0。
    C_O2=0 但 CH4>0 时，O2 右导数为 k C_CH4。
    """
    ch4 = float(C[IDX["CH4"]])
    o2 = float(C[IDX["O2"]])
    dr = np.zeros(N_MAJ, dtype=np.float64)
    t = max(float(T), _T_FLOOR)
    k = float(k_jensen_r7(R6_k0, float(R6_E_Rg), t))
    if o2 > 0.0:
        dr[IDX["CH4"]] = k * o2
    if ch4 > 0.0:
        dr[IDX["O2"]] = k * ch4
    if ch4 <= 0.0 or o2 <= 0.0:
        return 0.0, dr
    r = float(rate_R6(t, ch4, o2))
    return r, dr


def _r6_dr_dT(T: float, C: npt.NDArray[np.float64]) -> float:
    """冻结摩尔分率。k ∝ T^{-1} exp(-E_T/T)，C ∝ T^{-1}。"""
    ch4 = float(C[IDX["CH4"]])
    o2 = float(C[IDX["O2"]])
    if ch4 <= 0.0 or o2 <= 0.0:
        return 0.0
    t = max(float(T), _T_FLOOR)
    r = float(rate_R6(t, ch4, o2))
    dlnk = -1.0 / t + float(R6_E_Rg) / (t * t)
    return float(r * (dlnk - 2.0 / t))


def r6_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """R6 源项对 N_d/N_b/T 的解析块。两相同一 de Souza-Santos 式。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    T = float(z[-1])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    s = float(rate_scale)
    if abs(s) < 1.0e-30:
        return J
    for off, V in ((0, V_d), (7, V_b)):
        C, dCdN = _conc_and_dN(z[off:off + N_MAJ], T, P)
        r, drdC = _r6_r_and_dC(T, C)
        dxi_dN = (s * float(V)) * (drdC @ dCdN)
        J[off:off + N_MAJ, off:off + N_MAJ] += NU_R6_MAJOR[:, None] * dxi_dN
        dxi_dT = (s * float(V)) * _r6_dr_dT(T, C)
        J[off:off + N_MAJ, -1] += NU_R6_MAJOR * dxi_dT
        _ = r
    return J


def jacobian_r12_r5_r8_plus_scaled_r6(
    x: npt.NDArray[np.float64],
    nq_r8: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r6_scale: float,
) -> npt.NDArray[np.float64]:
    J8 = jacobian_r12_r5_plus_scaled_r8(x, freeze_r12_r5_q(), nq_r8, r8_scale=1.0)
    s = float(r6_scale)
    if abs(s) < 1.0e-30:
        return J8
    return J8 + s * r6_analytic_jacobian(x, nq_all, rate_scale=1.0)


def named_r12_r5_r8_r6_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    """与 R12/R5 同一套尺度。亏秩时 H2/O2/CO/CH4 补对角（修正 Newton）。"""
    from src.solvers.wirsum_1997_named_f_continuation import named_linear_step

    dx, rank, cond, _diag = named_linear_step(
        J, F, x, extra_diag_species=("H2", "O2", "CO", "CH4"),
    )
    return dx, rank, cond


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _newton_at_r6_scale(
    nq_r8: NamedStageQ,
    nq_all: NamedStageQ,
    x: npt.NDArray[np.float64],
    r6_scale: float,
    max_iter: int,
) -> Any:
    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_r12_r5_r8_plus_scaled_r6(z, nq_r8, nq_all, r6_scale=r6_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_r12_r5_r8_plus_scaled_r6(z, nq_r8, nq_all, r6_scale=r6_scale)

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_r8_r6_newton_step,
    )


def r6_source_snapshot(
    x: npt.NDArray[np.float64],
    nq_r8: NamedStageQ,
    nq_all: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f8 = coupled_residual_named(z, nq_r8, rate_scale=1.0)
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    dF = both - f8
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_CH4_d": float(dF[IDX["CH4"]]),
        "dF_O2_d": float(dF[IDX["O2"]]),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_H2O_d": float(dF[IDX["H2O"]]),
        "dF_CH4_b": float(dF[7 + IDX["CH4"]]),
        "dF_O2_b": float(dF[7 + IDX["O2"]]),
        "dF_CO_b": float(dF[7 + IDX["CO"]]),
        "dF_energy": float(dF[-1]),
        "N_d_CH4": float(z[IDX["CH4"]]),
        "N_d_O2": float(z[IDX["O2"]]),
        "N_d_CO": float(z[IDX["CO"]]),
        "N_d_H2": float(z[IDX["H2"]]),
        "N_b_CH4": float(z[7 + IDX["CH4"]]),
        "N_b_O2": float(z[7 + IDX["O2"]]),
        "N_b_CO": float(z[7 + IDX["CO"]]),
        "T_K": float(z[-1]),
    }


def solve_named_r12_r5_r8_then_r6(
    *,
    max_iter: int = 20,
    continuation: bool = True,
    r6_scale: float | None = None,
    apply_startwert: bool = True,
    r8: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E0+R12+R5+R8 数学根做 R6 延拓。默认先做同相 CH4∩O2 Startwert。"""
    nq_r8 = freeze_r12_r5_r8_q()
    nq_all = freeze_r12_r5_r8_r6_q()
    if r8 is None:
        r8 = solve_named_r12_r5_then_r8(max_iter=max_iter)
    x_start = np.asarray(r8["x"], dtype=np.float64).copy()
    snap_root = r6_source_snapshot(x_start, nq_r8, nq_all)
    sw6 = None
    sw5 = None
    if apply_startwert:
        packed6 = apply_named_r6_startwert(x_start)
        x_start = packed6["x"]
        sw6 = {"dense": packed6["dense"], "bubble": packed6["bubble"]}
        packed5 = apply_named_r5_startwert(x_start)
        x_start = packed5["x"]
        sw5 = {"dense": packed5["dense"], "bubble": packed5["bubble"]}
    snap0 = r6_source_snapshot(x_start, nq_r8, nq_all)
    if r6_scale is not None:
        scales = (float(r6_scale),)
    elif continuation:
        scales = _R6_STARTWERT_CONTINUATION if apply_startwert else _R6_CONTINUATION
    else:
        scales = (1.0,)
    n = n_unknowns(nq_all.stage)
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    for s in scales:
        last = _newton_at_r6_scale(nq_r8, nq_all, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_r12_r5_r8_plus_scaled_r6(last.x, nq_r8, nq_all, r6_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        history.append(
            {
                "r6_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "N_d_CH4": float(last.x[IDX["CH4"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
                "N_d_CO": float(last.x[IDX["CO"]]),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_H2O": float(last.x[IDX["H2O"]]),
                "N_b_CH4": float(last.x[7 + IDX["CH4"]]),
                "N_b_O2": float(last.x[7 + IDX["O2"]]),
            }
        )
        if int(last.rank_last) == n and bool(stop_s["converged_newton"]):
            x_cur = np.asarray(last.x, dtype=np.float64)

    assert last is not None
    stop = classify_newton_stop(
        last, model_is_full_rate=abs(float(scales[-1]) - 1.0) < 1.0e-15
    )
    startwert_abandoned = False
    startwert_history = None
    # Hegermann 焦炭给料提高床温后，同相 Startwert 会把 x 拉离刚硬的 R8 根，
    # 从 s=1e-4 无法收回。不改 rate_R6；改从 R8 根做无 Startwert 延拓。
    if (
        bool(apply_startwert)
        and r6_scale is None
        and continuation
        and not bool(stop["physical_convergence"])
    ):
        startwert_abandoned = True
        startwert_history = history
        x_cur = np.asarray(r8["x"], dtype=np.float64).copy()
        scales = _R6_CONTINUATION
        history = []
        last = None
        for s in scales:
            last = _newton_at_r6_scale(nq_r8, nq_all, x_cur, s, max_iter)
            stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
            J = jacobian_r12_r5_r8_plus_scaled_r6(last.x, nq_r8, nq_all, r6_scale=s)
            info = named_r12_scaled_rank(J, last.F, last.x)
            history.append(
                {
                    "r6_scale": float(s),
                    "residual_l2": float(stop_s["residual_l2"]),
                    "T_K": float(last.x[-1]),
                    "rank_last": int(last.rank_last),
                    "raw_rank": int(info["raw_rank"]),
                    "equilibrated_rank": int(info["equilibrated_rank"]),
                    "failure_reason": stop_s["failure_reason"],
                    "converged_newton": bool(stop_s["converged_newton"]),
                    "N_d_CH4": float(last.x[IDX["CH4"]]),
                    "N_d_O2": float(last.x[IDX["O2"]]),
                    "N_d_CO": float(last.x[IDX["CO"]]),
                    "N_d_H2": float(last.x[IDX["H2"]]),
                    "N_d_H2O": float(last.x[IDX["H2O"]]),
                    "N_b_CH4": float(last.x[7 + IDX["CH4"]]),
                    "N_b_O2": float(last.x[7 + IDX["O2"]]),
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
        "T_r8_K": float(r8["T_K"]),
        "energy_W": float(energy_row_named(last.x, nq_all)),
        "f_name": "F_E0+R12+R5+R8+R6",
        "enabled_rxns": ("R12", "R5", "R8", "R6"),
        "r8_physical_convergence": bool(r8["physical_convergence"]),
        "r6_at_r8_root": snap_root,
        "r6_after_startwert": snap0,
        "startwert_r6": sw6,
        "startwert_r5_after_r6": sw5,
        "apply_startwert": bool(apply_startwert),
        "startwert_abandoned": bool(startwert_abandoned),
        "startwert_continuation": startwert_history,
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
    }
