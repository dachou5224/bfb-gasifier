"""仅 R12（Hamel thesis R6）接入具名 40 维核。不改 ``rate_R12``，不经 ``Reactor``。

Hamel PDF 55 / 印刷 91 Eq. 5.38–5.39 已复核：式子与生产函数一致。
本模块处理同相 H2∩O2 的 Startwert、C_H2^{3/2} 解析雅可比，
以及 H2 耗尽流形上的微量组分尺度 / 行列均衡（亏秩时 H2 补对角，修正 Newton）。
``rate_scale<1`` 只用于延拓。

Source: Hamel (1999) Eq. 5.38–5.39, §2.1 Startwertwahl; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.kinetics.arrhenius import k_jensen_h2
from src.kinetics.gas_reactions import R12_A, R12_E_Rg, R12_h2_order, R12_o2_order
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    default_x0_named,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    solve_inner_named,
    unpack_named,
)
from src.solvers.wirsum_1997_startwert_r12 import apply_same_phase_r12_startwert
from src.solvers.wirsum_1997_stiffness_full_gas import (
    IDX,
    N_MAJ,
    fd_jacobian,
    row_col_equilibrate,
    scale_linear_system,
)
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton
from src.core.solid_balance_reduced import jacobian_diagnostics

# H2 + 1/2 O2 → H2O。Ref: Hamel (1999) §5.2.2 计量，印刷 91。
NU_R12_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R12_MAJOR[IDX["H2"]] = -1.0
NU_R12_MAJOR[IDX["O2"]] = -0.5
NU_R12_MAJOR[IDX["H2O"]] = 1.0
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]
# 微量组分尺度：H2 耗尽后 |N_H2|~1e-9，不能再用 0.05 mol/s 当地板。
_X_FLOOR_GAS: float = 1.0e-12  # [mol/s]
_X_FLOOR_SOLID: float = 0.05  # [kg]
_X_FLOOR_T: float = 50.0  # [K]
_F_SCALE_GAS: float = 1.0  # [mol/s]
_F_SCALE_SOLID: float = 1.0e-3  # [kg/s]
_F_SCALE_ENERGY: float = 1.0e5  # [W]
# 解析对角：仅当均衡后仍亏秩时加在 H2 行，单位 [1/s] 量级的无量纲 1。
_TRACE_DIAG: float = 1.0


def freeze_r12_only_q() -> NamedStageQ:
    return freeze_stage_q(2, enabled_rxns=frozenset({"R12"}))


def freeze_r12_inert_q(nq: NamedStageQ) -> NamedStageQ:
    """同一份 q，关掉均相。用于 FD 无 R12 部分。"""
    return NamedStageQ(
        energy=nq.energy,
        stage=0,
        Q_W=nq.Q_W,
        rho_ash=nq.rho_ash,
        tau_daem_s=nq.tau_daem_s,
        enabled_rxns=frozenset(),
    )


def apply_named_r12_startwert(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    """同相 R12 Startwert，并同步 FB/旋风气接口。不改 F(x) 定义。"""
    z = np.asarray(x, dtype=np.float64).copy()
    dense = apply_same_phase_r12_startwert(z[0:7])
    bubble = apply_same_phase_r12_startwert(z[7:14])
    n_tot = z[0:7] + z[7:14]
    z[16:23] = n_tot
    z[23:30] = n_tot
    z[30:37] = n_tot
    return {"x": z, "dense": dense, "bubble": bubble}


def _conc_and_dN(
    N: npt.NDArray[np.float64], T: float, P: float
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    n = np.asarray(N, dtype=np.float64)
    ntot = float(np.sum(n))
    c_tot = float(P) / (Rg * float(T))  # [mol/m³]
    if abs(ntot) < _NTOT_FLOOR:
        return np.zeros(N_MAJ, dtype=np.float64), np.zeros((N_MAJ, N_MAJ), dtype=np.float64)
    C = (n / ntot) * c_tot
    dCdN = np.zeros((N_MAJ, N_MAJ), dtype=np.float64)
    for i in range(N_MAJ):
        for k in range(N_MAJ):
            dCdN[i, k] = c_tot * (((1.0 if i == k else 0.0) * ntot) - n[i]) / (ntot * ntot)
    return C, dCdN


def _r12_r_and_dC(
    T: float, C: npt.NDArray[np.float64]
) -> tuple[float, npt.NDArray[np.float64]]:
    """与 ``rate_R12`` 同式。Hamel Eq. 5.38：r = k C_H2^{3/2} C_O2。

    C_H2>0 且 C_O2=0：r=0，但 ∂r/∂C_O2 |_{0+} = k C_H2^{3/2} ≠ 0（氧一阶）。
    C_H2=0：H2 级数 3/2，∂r/∂C_H2 |_{0+} = 0。不得把两种边界写成同一套零导数。
    """
    k = float(k_jensen_h2(R12_A, R12_E_Rg, float(T)))
    h2 = float(C[IDX["H2"]])
    o2 = float(C[IDX["O2"]])
    dr = np.zeros(N_MAJ, dtype=np.float64)
    if h2 > 0.0 and o2 > 0.0:
        r = k * (h2**R12_h2_order) * (o2**R12_o2_order)
        dr[IDX["H2"]] = k * R12_h2_order * (h2 ** (R12_h2_order - 1.0)) * (o2**R12_o2_order)
        dr[IDX["O2"]] = k * (h2**R12_h2_order) * R12_o2_order * (o2 ** (R12_o2_order - 1.0))
        return float(r), dr
    if h2 > 0.0:
        dr[IDX["O2"]] = k * (h2**R12_h2_order) * R12_o2_order
    return 0.0, dr


def _r12_dr_dT(T: float, C: npt.NDArray[np.float64]) -> float:
    """冻结摩尔分率时 ``∂r/∂T`` [mol/(m³·s·K)]。"""
    h2 = float(C[IDX["H2"]])
    o2 = float(C[IDX["O2"]])
    if h2 <= 0.0 or o2 <= 0.0:
        return 0.0
    t = float(T)
    k = float(k_jensen_h2(R12_A, R12_E_Rg, t))
    r = k * (h2**R12_h2_order) * (o2**R12_o2_order)
    dlnk = 1.5 / t + float(R12_E_Rg) / (t * t)
    # C ∝ 1/T
    return float(r * (dlnk - (R12_h2_order + R12_o2_order) / t))


def r12_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """R12 源项对 N_d/N_b/T 的解析块。能量行无显式 R12。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    T = float(z[-1])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    s = float(rate_scale)
    blocks = ((0, V_d), (7, V_b))
    for off, V in blocks:
        C, dCdN = _conc_and_dN(z[off:off + N_MAJ], T, P)
        r, drdC = _r12_r_and_dC(T, C)
        dxi_dN = (s * float(V)) * (drdC @ dCdN)
        J[off:off + N_MAJ, off:off + N_MAJ] += NU_R12_MAJOR[:, None] * dxi_dN
        dxi_dT = (s * float(V)) * _r12_dr_dT(T, C)
        J[off:off + N_MAJ, -1] += NU_R12_MAJOR * dxi_dT
        _ = r
    return J


def named_r12_linear_scales(
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """具名 ``D_x``/``D_F``。气体地板 1e-12 mol/s，能量 1e5 W。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    xs = np.maximum(np.abs(z), _X_FLOOR_GAS)
    fs = np.full(n, _F_SCALE_GAS, dtype=np.float64)
    for i in (14, 15, 37, 38):
        if i < n - 1:
            xs[i] = max(abs(float(z[i])), _X_FLOOR_SOLID)
            fs[i] = _F_SCALE_SOLID
    if n >= 41:
        xs[39] = max(abs(float(z[39])), _X_FLOOR_SOLID)
        fs[39] = _F_SCALE_SOLID
    xs[-1] = max(abs(float(z[-1])), _X_FLOOR_T)
    fs[-1] = _F_SCALE_ENERGY
    return xs, fs


def named_r12_scaled_rank(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> dict[str, Any]:
    """原始秩 vs 微量组分尺度 + 行列均衡后的秩。"""
    xs, fs = named_r12_linear_scales(x)
    raw = jacobian_diagnostics(J)
    J_hat, F_hat = scale_linear_system(J, F, x_scale=xs, F_scale=fs)
    hat = jacobian_diagnostics(J_hat)
    J_eq, dr, dc = row_col_equilibrate(J_hat)
    eq = jacobian_diagnostics(J_eq)
    return {
        "raw_rank": int(raw.rank),
        "raw_condition": float(raw.cond),
        "explicit_scale_rank": int(hat.rank),
        "equilibrated_rank": int(eq.rank),
        "equilibrated_condition": float(eq.cond),
        "n": int(J.shape[0]),
        "x_scale": xs,
        "F_scale": fs,
        "J_hat": J_hat,
        "F_hat": F_hat,
        "J_eq": J_eq,
        "dr": dr,
        "dc": dc,
    }


def named_r12_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    """``J Δx = -F``：微量组分尺度 + 行列均衡。亏秩时给 H2 补对角（修正 Newton）。"""
    from src.solvers.wirsum_1997_named_f_continuation import named_linear_step

    dx, rank, cond, _diag = named_linear_step(
        J, F, x, extra_diag_species=("H2",),
    )
    return dx, rank, cond


def jacobian_named_r12(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """无 R12 的 FD 雅可比 + R12 解析块。避免流形上差分 C_H2^{3/2}。"""
    z = np.asarray(x, dtype=np.float64)
    nq0 = freeze_r12_inert_q(nq)
    scale = np.maximum(np.abs(z), 0.05)
    scale[-1] = max(float(scale[-1]), 50.0)
    J0 = fd_jacobian(lambda w: coupled_residual_named(w, nq0), z, scale=scale)
    return J0 + r12_analytic_jacobian(z, nq, rate_scale=rate_scale)


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


_R12_CONTINUATION: tuple[float, ...] = (
    0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4, 1.0
)


def _newton_at_scale(
    nq: NamedStageQ,
    x: npt.NDArray[np.float64],
    rate_scale: float,
    max_iter: int,
) -> Any:
    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return coupled_residual_named(z, nq, rate_scale=rate_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_named_r12(z, nq, rate_scale=rate_scale)

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_newton_step,
    )


def solve_named_r12_only(
    *,
    x0: npt.NDArray[np.float64] | None = None,
    from_e0_root: bool = True,
    apply_startwert: bool = False,
    rate_scale: float | None = None,
    continuation: bool = True,
    max_iter: int = 20,
) -> dict[str, Any]:
    """F_E0 + 仅 R12。默认从 F_E0 根做 rate_scale 延拓，混合雅可比。

    ``apply_startwert`` 只改 x0，不改残差。完整速率仍是 ``rate_scale=1``。
    微量组分尺度下可在 H2 耗尽流形上求得数学根；T 远低于 1100 K，
    不得写成炉温收敛。
    """
    nq = freeze_r12_only_q()
    if x0 is None and from_e0_root:
        e0 = solve_inner_named(freeze_stage_q(0), max_iter=40)
        x_start = np.asarray(e0["x"], dtype=np.float64).copy()
        e0_ok = bool(e0["physical_convergence"])
        t_e0 = float(e0["T_K"])
    elif x0 is None:
        x_start = default_x0_named(nq)
        e0_ok = False
        t_e0 = float(x_start[-1])
    else:
        x_start = np.asarray(x0, dtype=np.float64).copy()
        e0_ok = False
        t_e0 = float(x_start[-1])
    sw = None
    if apply_startwert:
        packed = apply_named_r12_startwert(x_start)
        x_start = packed["x"]
        sw = {"dense": packed["dense"], "bubble": packed["bubble"]}

    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    scales = (float(rate_scale),) if rate_scale is not None else (
        _R12_CONTINUATION if continuation else (1.0,)
    )
    n = n_unknowns(nq.stage)
    for s in scales:
        last = _newton_at_scale(nq, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        history.append(
            {
                "rate_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
            }
        )
        if int(last.rank_last) == n:
            x_cur = np.asarray(last.x, dtype=np.float64)

    assert last is not None
    stop = classify_newton_stop(last, model_is_full_rate=abs(float(scales[-1]) - 1.0) < 1.0e-15)
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "T_K": float(last.x[-1]),
        "T_e0_K": t_e0,
        "energy_W": float(energy_row_named(last.x, nq)),
        "f_name": "F_E0+R12",
        "enabled_rxns": ("R12",),
        "from_e0_root": bool(from_e0_root),
        "e0_root_physical": e0_ok,
        "startwert": sw,
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
    }


def r12_x0_explosion_snapshot() -> dict[str, Any]:
    """默认 x0（含密相 O2 种子）上的 R12 爆炸对照。"""
    nq = freeze_r12_only_q()
    x = default_x0_named(nq)
    F0 = coupled_residual_named(x, nq)
    packed = apply_named_r12_startwert(x)
    F1 = coupled_residual_named(packed["x"], nq)
    u = unpack_named(x, nq.stage)
    return {
        "n_unknowns": n_unknowns(nq.stage),
        "x0_l2": float(np.linalg.norm(F0)),
        "after_startwert_l2": float(np.linalg.norm(F1)),
        "N_d_H2": float(u["N_d"][IDX["H2"]]),
        "N_d_O2": float(u["N_d"][IDX["O2"]]),
        "N_b_O2": float(u["N_b"][IDX["O2"]]),
        "dense_xi12_mol_s": float(packed["dense"]["xi12_mol_s"]),
        "not_a_furnace_convergence_claim": True,
    }
