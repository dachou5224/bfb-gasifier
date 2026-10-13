"""F_E0+R12 数学根上接入 R5。不改 ``rate_R5_*`` / ``rate_R12``，不经 ``Reactor``。

R12 保持完整速率。R5 用 ``rate_scale`` 延拓：
F(x;s) = F_{E0+R12}(x) + s (F_{E0+R12+R5}(x) - F_{E0+R12}(x)).

气泡相 Eq. 5.36：C_CO C_O2^{0.3} C_H2O^{0.5}；悬浮相 Table 5.3 分压式。
C_i=0 时：CO 一阶保留右导数；O2/H2O 分数幂在 0 处发散，不写有限右导数。

无 Startwert 时气泡仍是 CO∩O2，Hottel 源 ~10^3 mol/s，延拓只能接到
``r5_scale≤1e-6``。默认在 R12 根上做同相 R5 Startwert；C_CO=0 时保留
CO 右导数。完整速率可落到数学根（T 升到 ~1100 K 量级），仍不是炉温验收。

Source: Hamel (1999) Eq. 5.36–5.37, Table 5.3, §2.1; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import P0, Rg
from src.kinetics.arrhenius import k_standard
from src.kinetics.gas_reactions import (
    R5_bubble_E_Rg,
    R5_bubble_beta,
    R5_bubble_gamma,
    R5_bubble_k0,
    R5_suspension_k_T05,
    rate_R5_bubble,
    rate_R5_suspension,
)
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
)
from src.solvers.wirsum_1997_named_f_r12 import (
    _conc_and_dN,
    freeze_r12_only_q,
    jacobian_named_r12,
    named_r12_scaled_rank,
    solve_named_r12_only,
)
from src.solvers.wirsum_1997_named_f_continuation import (
    AdvanceMode,
    certify_named_full_rate_endpoint,
    shrink_scale,
    stage_passed,
)
from src.solvers.wirsum_1997_startwert_r5 import apply_named_r5_startwert
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

# CO + 1/2 O2 → CO2。Ref: Hamel (1999) Table 5.3 计量。
NU_R5_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R5_MAJOR[IDX["CO"]] = -1.0
NU_R5_MAJOR[IDX["O2"]] = -0.5
NU_R5_MAJOR[IDX["CO2"]] = 1.0
_T_SUSP_FLOOR: float = 300.0  # [K] 与 rate_R5_suspension 相同
_DUMMY_Y: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)

_R5_CONTINUATION: tuple[float, ...] = (
    0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4, 1.0
)
# Startwert 后 s=0 会按 F_E0+R12 把 CO/O2 再混回去；从有限 R5 速率起步。
_R5_STARTWERT_CONTINUATION: tuple[float, ...] = (1.0e-4, 1.0e-2, 1.0e-1, 1.0)


def freeze_r12_r5_q() -> NamedStageQ:
    return freeze_stage_q(2, enabled_rxns=frozenset({"R12", "R5"}))


def residual_r12_plus_scaled_r5(
    x: npt.NDArray[np.float64],
    nq12: NamedStageQ,
    nq_both: NamedStageQ,
    *,
    r5_scale: float,
) -> npt.NDArray[np.float64]:
    """R12 完整速率，R5 乘 ``r5_scale``。s=0 退回 F_E0+R12。"""
    z = np.asarray(x, dtype=np.float64)
    f12 = coupled_residual_named(z, nq12, rate_scale=1.0)
    s = float(r5_scale)
    if abs(s) < 1.0e-30:
        return f12
    both = coupled_residual_named(z, nq_both, rate_scale=1.0)
    return f12 + s * (both - f12)


def _r5_bubble_r_and_dC(
    T: float, C: npt.NDArray[np.float64]
) -> tuple[float, npt.NDArray[np.float64]]:
    """与 ``rate_R5_bubble`` 同式。Hamel 气泡：C_CO C_O2^{0.3} C_H2O^{0.5}。

    C_CO=0 但氧化剂仍在：CO 一阶，右导数有限，必须保留。
    C_O2=0 或 C_H2O=0：分数幂导数在 0 处发散，**不能**套用 R12 那种有限右导数；
    此处 r=0 且分数幂方向导数记为无有限值（分量保持 0）。
    """
    co = float(C[IDX["CO"]])
    o2 = float(C[IDX["O2"]])
    h2o = float(C[IDX["H2O"]])
    dr = np.zeros(N_MAJ, dtype=np.float64)
    if o2 <= 0.0 or h2o <= 0.0:
        return 0.0, dr
    e_j = float(R5_bubble_E_Rg) * Rg  # [J/mol]
    k = float(k_standard(R5_bubble_k0, e_j, float(T)))
    beta = float(R5_bubble_beta)
    gamma = float(R5_bubble_gamma)
    dr[IDX["CO"]] = k * (o2**beta) * (h2o**gamma)
    if co <= 0.0:
        return 0.0, dr
    r = float(rate_R5_bubble(T, co, o2, h2o, P0, _DUMMY_Y, None))
    dr[IDX["O2"]] = k * co * beta * (o2 ** (beta - 1.0)) * (h2o**gamma)
    dr[IDX["H2O"]] = k * co * (o2**beta) * gamma * (h2o ** (gamma - 1.0))
    return r, dr


def _r5_suspension_r_and_dC(
    T: float, C: npt.NDArray[np.float64]
) -> tuple[float, npt.NDArray[np.float64]]:
    """与 ``rate_R5_suspension`` 同式。p_i = C_i Rg T / P0。

    CO 一阶：C_CO=0 时右导数保留。O2/H2O 为 0.5 次：零处导数发散，不写有限右导数。
    """
    co = float(C[IDX["CO"]])
    o2 = float(C[IDX["O2"]])
    h2o = float(C[IDX["H2O"]])
    dr = np.zeros(N_MAJ, dtype=np.float64)
    if o2 <= 0.0 or h2o <= 0.0:
        return 0.0, dr
    t = max(float(T), _T_SUSP_FLOOR)
    k = float(R5_suspension_k_T05) * (t ** -0.5)
    alpha = (Rg * t) / P0  # [atm·m³/mol]
    pref = k * (alpha**2)
    dr[IDX["CO"]] = pref * (o2**0.5) * (h2o**0.5)
    if co <= 0.0:
        return 0.0, dr
    r = float(rate_R5_suspension(t, co, o2, h2o, P0, _DUMMY_Y, None))
    dr[IDX["O2"]] = pref * co * 0.5 * (o2 ** -0.5) * (h2o**0.5)
    dr[IDX["H2O"]] = pref * co * (o2**0.5) * 0.5 * (h2o ** -0.5)
    return r, dr


def _r5_bubble_dr_dT(T: float, C: npt.NDArray[np.float64]) -> float:
    """冻结摩尔分率。r ∝ k(T) T^{-1.8}。"""
    co = float(C[IDX["CO"]])
    o2 = float(C[IDX["O2"]])
    h2o = float(C[IDX["H2O"]])
    if co <= 0.0 or o2 <= 0.0 or h2o <= 0.0:
        return 0.0
    t = float(T)
    r = float(rate_R5_bubble(t, co, o2, h2o, P0, _DUMMY_Y, None))
    order = 1.0 + float(R5_bubble_beta) + float(R5_bubble_gamma)
    dlnk = float(R5_bubble_E_Rg) / (t * t)
    return float(r * (dlnk - order / t))


def _r5_suspension_dr_dT(T: float, C: npt.NDArray[np.float64]) -> float:
    """冻结 y 时 p_i 与 T 无关，r ∝ T^{-0.5}。"""
    co = float(C[IDX["CO"]])
    o2 = float(C[IDX["O2"]])
    h2o = float(C[IDX["H2O"]])
    if co <= 0.0 or o2 <= 0.0 or h2o <= 0.0:
        return 0.0
    t = max(float(T), _T_SUSP_FLOOR)
    r = float(rate_R5_suspension(t, co, o2, h2o, P0, _DUMMY_Y, None))
    return float(r * (-0.5 / t))


def r5_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """R5 源项对 N_d/N_b/T 的解析块。密相悬浮式，气泡 Hottel。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    T = float(z[-1])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    s = float(rate_scale)
    specs = (
        (0, V_d, _r5_suspension_r_and_dC, _r5_suspension_dr_dT),
        (7, V_b, _r5_bubble_r_and_dC, _r5_bubble_dr_dT),
    )
    for off, V, r_fn, dt_fn in specs:
        C, dCdN = _conc_and_dN(z[off:off + N_MAJ], T, P)
        r, drdC = r_fn(T, C)
        dxi_dN = (s * float(V)) * (drdC @ dCdN)
        J[off:off + N_MAJ, off:off + N_MAJ] += NU_R5_MAJOR[:, None] * dxi_dN
        dxi_dT = (s * float(V)) * dt_fn(T, C)
        J[off:off + N_MAJ, -1] += NU_R5_MAJOR * dxi_dT
        _ = r
    return J


def jacobian_r12_plus_scaled_r5(
    x: npt.NDArray[np.float64],
    nq12: NamedStageQ,
    nq_both: NamedStageQ,
    *,
    r5_scale: float,
) -> npt.NDArray[np.float64]:
    J12 = jacobian_named_r12(x, nq12, rate_scale=1.0)
    s = float(r5_scale)
    if abs(s) < 1.0e-30:
        return J12
    return J12 + s * r5_analytic_jacobian(x, nq_both, rate_scale=1.0)


def named_r12_r5_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    """与 R12 同一套尺度。亏秩时 H2/O2/CO 补对角（修正 Newton，非解析导数）。"""
    from src.solvers.wirsum_1997_named_f_continuation import named_linear_step

    dx, rank, cond, _diag = named_linear_step(
        J, F, x, extra_diag_species=("H2", "O2", "CO"),
    )
    return dx, rank, cond


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _newton_at_r5_scale(
    nq12: NamedStageQ,
    nq_both: NamedStageQ,
    x: npt.NDArray[np.float64],
    r5_scale: float,
    max_iter: int,
) -> Any:
    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_r12_plus_scaled_r5(z, nq12, nq_both, r5_scale=r5_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_r12_plus_scaled_r5(z, nq12, nq_both, r5_scale=r5_scale)

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_newton_step,
    )


def r5_source_snapshot(
    x: npt.NDArray[np.float64],
    nq12: NamedStageQ,
    nq_both: NamedStageQ,
) -> dict[str, Any]:
    """R12 根上 R5 源项有多大。不改 x。"""
    z = np.asarray(x, dtype=np.float64)
    f12 = coupled_residual_named(z, nq12, rate_scale=1.0)
    both = coupled_residual_named(z, nq_both, rate_scale=1.0)
    dF = both - f12
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_O2_d": float(dF[IDX["O2"]]),
        "dF_CO2_d": float(dF[IDX["CO2"]]),
        "dF_CO_b": float(dF[7 + IDX["CO"]]),
        "dF_O2_b": float(dF[7 + IDX["O2"]]),
        "dF_energy": float(dF[-1]),
        "N_d_CO": float(z[IDX["CO"]]),
        "N_d_O2": float(z[IDX["O2"]]),
        "N_d_H2O": float(z[IDX["H2O"]]),
        "N_d_H2": float(z[IDX["H2"]]),
        "N_b_CO": float(z[7 + IDX["CO"]]),
        "N_b_O2": float(z[7 + IDX["O2"]]),
        "T_K": float(z[-1]),
    }


def solve_named_r12_then_r5(
    *,
    max_iter: int = 20,
    continuation: bool = True,
    r5_scale: float | None = None,
    apply_startwert: bool = True,
    advance_mode: AdvanceMode = "root_correcting",
    s_min: float = 1.0e-16,
    r12: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E0+R12 数学根做 R5 延拓。默认逐根校正。"""
    nq12 = freeze_r12_only_q()
    nq_both = freeze_r12_r5_q()
    if r12 is None:
        r12 = solve_named_r12_only(max_iter=max_iter)
    x_start = np.asarray(r12["x"], dtype=np.float64).copy()
    snap_root = r5_source_snapshot(x_start, nq12, nq_both)
    sw = None
    if apply_startwert:
        packed = apply_named_r5_startwert(x_start)
        x_start = packed["x"]
        sw = {"dense": packed["dense"], "bubble": packed["bubble"]}
    snap0 = r5_source_snapshot(x_start, nq12, nq_both)
    if r5_scale is not None:
        planned = [float(r5_scale)]
    elif continuation:
        planned = list(
            _R5_STARTWERT_CONTINUATION if apply_startwert else _R5_CONTINUATION
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
        last = _newton_at_r5_scale(nq12, nq_both, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_r12_plus_scaled_r5(last.x, nq12, nq_both, r5_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        passed = stage_passed(
            stop_s, last.x, advance_mode=advance_mode,
            accepted_alphas=list(last.accepted_alphas),
        )
        history.append(
            {
                "r5_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "equilibrated_condition": float(info["equilibrated_condition"]),
                "solver_stop_reason": stop_s["failure_reason"],
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "stage_passed": bool(passed),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
                "N_d_CO": float(last.x[IDX["CO"]]),
                "N_d_CO2": float(last.x[IDX["CO2"]]),
                "N_d_H2O": float(last.x[IDX["H2O"]]),
                "N_b_CO": float(last.x[7 + IDX["CO"]]),
                "N_b_O2": float(last.x[7 + IDX["O2"]]),
            }
        )
        if passed:
            x_cur = np.asarray(last.x, dtype=np.float64)
            s_ok = s
            i += 1
            continue
        branch_tracked = False
        if advance_mode == "root_correcting":
            s_mid = shrink_scale(s_ok, s)
            if abs(s - s_ok) < float(s_min) or abs(s_mid - s) < 1.0e-18:
                break
            planned.insert(i, s_mid)
            continue
        i += 1

    assert last is not None
    s_last = float(history[-1]["r5_scale"]) if history else 1.0
    stop = classify_newton_stop(
        last, model_is_full_rate=abs(s_last - 1.0) < 1.0e-15
    )
    cert = certify_named_full_rate_endpoint(last.x, nq_both, coupling="frozen_q")
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "T_K": float(last.x[-1]),
        "T_r12_K": float(r12["T_K"]),
        "energy_W": float(energy_row_named(last.x, nq_both)),
        "f_name": "F_E0+R12+R5",
        "enabled_rxns": ("R12", "R5"),
        "r12_physical_convergence": bool(r12["physical_convergence"]),
        "r5_at_r12_root": snap_root,
        "r5_after_startwert": snap0,
        "startwert": sw,
        "apply_startwert": bool(apply_startwert),
        "advance_mode": str(advance_mode),
        "branch_tracked": bool(branch_tracked),
        "solver_stop_reason": stop["failure_reason"],
        "endpoint_certified": bool(cert["endpoint_certified"]),
        "endpoint_cert": {k: v for k, v in cert.items() if k != "F"},
        "coupling_mode": "frozen_q",
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
    }
