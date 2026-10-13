"""阶段 E：最小联立刚性阶梯。不经 ``Reactor``，不改生产默认。

阶梯（后一级必须回到完整生产速率，延拓不是最终模型）：

1. ``isothermal_solid_transport`` — 已有 holdup 开链，线性
2. ``isothermal_two_phase_exchange`` — 示踪 + Hamel Eq. 2.2 交换，线性
3. ``isothermal_r12_exchange`` — 交换 + 生产 ``rate_R12``（Hamel 5.38）
4. ``full_gas`` / ``char`` / ``energy`` — ``wirsum_1997_stiffness_full_gas.py``

Newton 用 ``wirsum_1997_strict_newton``。``raw_l2`` 是 Wirsum Eq. 2.23；
``scaled_rms`` 是工程对照。线搜索失败、奇异、触 ``λ_min`` 不得标成物理收敛。
禁止在残差里 clip，禁止恢复 fastox。

Source: Hamel (1999) Eq. 2.2, 5.38–5.39; Wirsum (1997) §2.3.3;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 E
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.solid_balance_reduced import (
    canonical_open_chain_K,
    holdup_jacobian,
    holdup_residual,
    jacobian_diagnostics,
)
from src.core.species import N_GAS
from src.kinetics.arrhenius import k_jensen_h2
from src.kinetics.gas_reactions import R12_A, R12_E_Rg, R12_h2_order, R12_o2_order, rate_R12
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

RungName = Literal[
    "isothermal_solid_transport",
    "isothermal_two_phase_exchange",
    "isothermal_r12_exchange",
    "full_gas",
    "char",
    "energy",
]

SPECIES: tuple[str, ...] = ("H2", "O2", "H2O", "N2")
N_SPEC: int = 4
IDX_H2: int = 0
IDX_O2: int = 1
IDX_H2O: int = 2
IDX_N2: int = 3
NU_R12: npt.NDArray[np.float64] = np.array([-1.0, -0.5, 1.0, 0.0], dtype=np.float64)
_Y_DUMMY = np.zeros(N_GAS, dtype=np.float64)
_UNIMPLEMENTED = frozenset()
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]


@dataclass(frozen=True)
class ExchangeGeometry:
    """单格等温双相几何。全程 SI。"""

    T: float = 1100.0  # [K]
    P: float = 2.5e6  # [Pa]
    K_bd: float = 6.57  # [1/s]
    V_b: float = 0.10  # [m³]
    V_d: float = 0.32  # [m³]


def implemented_rungs() -> tuple[str, ...]:
    return (
        "isothermal_solid_transport",
        "isothermal_two_phase_exchange",
        "isothermal_r12_exchange",
        "full_gas",
        "char",
        "energy",
    )


def assert_rung_implemented(rung: str) -> None:
    if rung not in implemented_rungs():
        raise ValueError(f"未知 rung={rung!r}")


_PHYSICAL_L2: float = 1.0e-6  # 与 F 同单位的欧氏阈值；小步长停机不得越过它


def classify_newton_stop(result: Any, *, model_is_full_rate: bool) -> dict[str, Any]:
    """把 Newton 停止原因与「物理收敛」分开。

    Wirsum 核可因 ``||αΔx||`` 过小停机；那只是 ``small_step``，不是物理根。
    """
    l2 = float(np.linalg.norm(result.F))
    if result.converged and l2 < _PHYSICAL_L2:
        reason = "newton_merit"
    elif result.converged:
        reason = "small_step"
    else:
        reason = str(result.failure_reason or "")
    return {
        "converged_newton": bool(result.converged),
        "failure_reason": reason,
        "model_is_full_rate": bool(model_is_full_rate),
        "physical_convergence": bool(model_is_full_rate) and l2 < _PHYSICAL_L2,
        "iterations": int(result.iterations),
        "residual_l2": l2,
        "accepted_alphas": [float(a) for a in result.accepted_alphas],
        "merit_kind": str(result.merit_kind),
    }


def linear_backward_error(
    J: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
) -> float:
    """相对后向误差 ``||F|| / (||J||_F ||x||)``，F 与 Jx 同单位。"""
    jn = float(np.linalg.norm(J, ord="fro"))
    xn = float(np.linalg.norm(x))
    fn = float(np.linalg.norm(F))
    return fn / max(jn * xn, 1.0e-30)


def newton_step_backward_error(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    dx: npt.NDArray[np.float64],
) -> float:
    """``||J Δx + F|| / (||J||_F ||Δx|| + ||F||)``。"""
    resid = J @ dx + F
    denom = float(np.linalg.norm(J, ord="fro")) * float(np.linalg.norm(dx)) + float(
        np.linalg.norm(F)
    )
    return float(np.linalg.norm(resid)) / max(denom, 1.0e-30)


def directional_derivative_error(
    residual_fn: Any,
    jacobian_fn: Any,
    x: npt.NDArray[np.float64],
    direction: npt.NDArray[np.float64],
    eps: float,
) -> float:
    """``||(F(x+εv)-F(x))/ε - J v|| / max(||J v||, 1)``。"""
    v = np.asarray(direction, dtype=np.float64)
    vn = float(np.linalg.norm(v))
    if vn == 0.0:
        raise ValueError("direction 不能为零")
    v = v / vn
    x = np.asarray(x, dtype=np.float64)
    Jv = jacobian_fn(x) @ v
    fd = (residual_fn(x + float(eps) * v) - residual_fn(x)) / float(eps)
    return float(np.linalg.norm(fd - Jv)) / max(float(np.linalg.norm(Jv)), 1.0)


# --- rung 1: solid transport -------------------------------------------------

def _transport_problem() -> dict[str, npt.NDArray[np.float64]]:
    K_up, K_down = canonical_open_chain_K()
    m_feed = np.array([1.0, 0.0, 0.0], dtype=np.float64)  # [kg/s]
    r_dot = np.zeros(3, dtype=np.float64)  # [kg/s]
    return {"K_up": K_up, "K_down": K_down, "m_feed": m_feed, "r_dot": r_dot}


def transport_residual(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    p = _transport_problem()
    return holdup_residual(M, k0_anchor="none", **p)


def transport_jacobian(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    p = _transport_problem()
    return holdup_jacobian(M, k0_anchor="none", **p)


# --- rung 2: tracer exchange -------------------------------------------------

def _c_tot(geo: ExchangeGeometry) -> float:
    return float(geo.P) / (Rg * float(geo.T))  # [mol/m³]


def _exchange_coeff(geo: ExchangeGeometry) -> float:
    """``K_bd V_b c_tot`` [mol/s]。"""
    return float(geo.K_bd) * float(geo.V_b) * _c_tot(geo)


def tracer_exchange_residual(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_d_in: float,
    N_b_in: float,
    N_tot_d: float,
    N_tot_b: float,
) -> npt.NDArray[np.float64]:
    """示踪摩尔流残差 [mol/s]。``N_tot`` 由惰性补齐，不进未知量。"""
    N_d, N_b = float(x[0]), float(x[1])
    beta = _exchange_coeff(geo)
    N_ex = beta * (N_b / float(N_tot_b) - N_d / float(N_tot_d))
    return np.array([N_d_in - N_d + N_ex, N_b_in - N_b - N_ex], dtype=np.float64)


def tracer_exchange_jacobian(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_tot_d: float,
    N_tot_b: float,
) -> npt.NDArray[np.float64]:
    _ = x
    beta = _exchange_coeff(geo)
    d_ex_d = -beta / float(N_tot_d)
    d_ex_b = beta / float(N_tot_b)
    return np.array(
        [[-1.0 + d_ex_d, d_ex_b], [-d_ex_d, -1.0 - d_ex_b]],
        dtype=np.float64,
    )


def exact_tracer_exchange(
    *,
    geo: ExchangeGeometry,
    N_d_in: float,
    N_b_in: float,
    N_tot_d: float,
    N_tot_b: float,
) -> npt.NDArray[np.float64]:
    """线性示踪交换的闭式解 [mol/s]。"""
    beta = _exchange_coeff(geo)
    a = 1.0 + beta / float(N_tot_d)
    b = -beta / float(N_tot_b)
    c = -beta / float(N_tot_d)
    d = 1.0 + beta / float(N_tot_b)
    det = a * d - b * c
    N_d = (d * float(N_d_in) - b * float(N_b_in)) / det
    N_b = (-c * float(N_d_in) + a * float(N_b_in)) / det
    return np.array([N_d, N_b], dtype=np.float64)


# --- rung 3: R12 + exchange --------------------------------------------------

def _phase_conc(N: npt.NDArray[np.float64], geo: ExchangeGeometry) -> npt.NDArray[np.float64]:
    ntot = float(np.sum(N))
    if abs(ntot) < _NTOT_FLOOR:
        return np.zeros(N_SPEC, dtype=np.float64)
    return (np.asarray(N, dtype=np.float64) / ntot) * _c_tot(geo)


def _dC_dN(N: npt.NDArray[np.float64], geo: ExchangeGeometry) -> npt.NDArray[np.float64]:
    n = np.asarray(N, dtype=np.float64)
    ntot = float(np.sum(n))
    J = np.zeros((N_SPEC, N_SPEC), dtype=np.float64)
    if abs(ntot) < _NTOT_FLOOR:
        return J
    c = _c_tot(geo)
    for i in range(N_SPEC):
        for k in range(N_SPEC):
            J[i, k] = c * (((1.0 if i == k else 0.0) * ntot) - n[i]) / (ntot * ntot)
    return J


def _r12_rate(T: float, C: npt.NDArray[np.float64], P: float) -> float:
    return float(rate_R12(float(T), float(C[IDX_H2]), float(C[IDX_O2]), float(P), _Y_DUMMY, None))


def _r12_rate_and_dC(T: float, C: npt.NDArray[np.float64]) -> tuple[float, npt.NDArray[np.float64]]:
    """与 ``rate_R12`` 同式；导数在 C≤0 的分量取 0。"""
    k = float(k_jensen_h2(R12_A, R12_E_Rg, float(T)))
    h2 = float(C[IDX_H2])
    o2 = float(C[IDX_O2])
    h2c = max(h2, 0.0)
    o2c = max(o2, 0.0)
    r = k * (h2c**R12_h2_order) * (o2c**R12_o2_order)
    dr = np.zeros(N_SPEC, dtype=np.float64)
    if h2 > 0.0 and o2 > 0.0:
        dr[IDX_H2] = k * R12_h2_order * (h2 ** (R12_h2_order - 1.0)) * (o2**R12_o2_order)
        dr[IDX_O2] = k * (h2**R12_h2_order) * R12_o2_order * (o2 ** (R12_o2_order - 1.0))
    return float(r), dr


def unpack_r12_state(x: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (2 * N_SPEC,), "R12 状态必须是 8 元摩尔流"
    return z[:N_SPEC].copy(), z[N_SPEC:].copy()


def r12_exchange_residual(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_in_d: npt.NDArray[np.float64],
    N_in_b: npt.NDArray[np.float64],
    rate_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    """交换 + R12。``rate_scale<1`` 只用于延拓，不是物理模型。"""
    N_d, N_b = unpack_r12_state(x)
    C_d = _phase_conc(N_d, geo)
    C_b = _phase_conc(N_b, geo)
    N_ex = float(geo.K_bd) * float(geo.V_b) * (C_b - C_d)  # [mol/s]
    xi_d = float(rate_scale) * _r12_rate(geo.T, C_d, geo.P) * float(geo.V_d)
    xi_b = float(rate_scale) * _r12_rate(geo.T, C_b, geo.P) * float(geo.V_b)
    F_d = N_in_d - N_d + N_ex + NU_R12 * xi_d
    F_b = N_in_b - N_b - N_ex + NU_R12 * xi_b
    return np.concatenate([F_d, F_b])


def r12_exchange_jacobian(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_in_d: npt.NDArray[np.float64],
    N_in_b: npt.NDArray[np.float64],
    rate_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    _ = (N_in_d, N_in_b)
    N_d, N_b = unpack_r12_state(x)
    C_d = _phase_conc(N_d, geo)
    C_b = _phase_conc(N_b, geo)
    dCd = _dC_dN(N_d, geo)
    dCb = _dC_dN(N_b, geo)
    beta_v = float(geo.K_bd) * float(geo.V_b)
    dNex_dNd = -beta_v * dCd
    dNex_dNb = beta_v * dCb
    _r_d, drdC_d = _r12_rate_and_dC(geo.T, C_d)
    _r_b, drdC_b = _r12_rate_and_dC(geo.T, C_b)
    dxi_d = float(rate_scale) * float(geo.V_d) * (drdC_d @ dCd)
    dxi_b = float(rate_scale) * float(geo.V_b) * (drdC_b @ dCb)
    J = np.zeros((2 * N_SPEC, 2 * N_SPEC), dtype=np.float64)
    J[:N_SPEC, :N_SPEC] = -np.eye(N_SPEC) + dNex_dNd + NU_R12[:, None] * dxi_d
    J[:N_SPEC, N_SPEC:] = dNex_dNb
    J[N_SPEC:, :N_SPEC] = -dNex_dNd
    J[N_SPEC:, N_SPEC:] = -np.eye(N_SPEC) - dNex_dNb + NU_R12[:, None] * dxi_b
    return J


def default_r12_inlets() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """H₂ 只进 dense，O₂ 只进 bubble：必须经交换才能反应。"""
    N_in_d = np.array([0.50, 0.0, 0.0, 10.0], dtype=np.float64)  # [mol/s]
    N_in_b = np.array([0.0, 0.30, 0.0, 10.0], dtype=np.float64)  # [mol/s]
    return N_in_d, N_in_b


def r12_initial_guesses(
    N_in_d: npt.NDArray[np.float64],
    N_in_b: npt.NDArray[np.float64],
) -> dict[str, npt.NDArray[np.float64]]:
    separated = np.concatenate([N_in_d, N_in_b])
    mixed = 0.5 * np.concatenate([N_in_d + N_in_b, N_in_d + N_in_b])
    converted = np.array(
        [0.20, 1.0e-12, 0.30, 10.0, 1.0e-12, 1.0e-12, 0.0, 10.0],
        dtype=np.float64,
    )
    return {"separated_inlet": separated, "mixed_overlap": mixed, "stoich_converted": converted}


# --- payload -----------------------------------------------------------------

def build_stage_e_payload() -> dict[str, Any]:
    """阶梯证据。``solved`` 只表示本核跑完，不是全炉收敛。"""
    from src.solvers.wirsum_1997_stiffness_full_gas import build_full_gas_char_energy_section
    geo = ExchangeGeometry()
    N_d_in_tr, N_b_in_tr = 1.0, 0.0  # [mol/s]
    N_tot_d, N_tot_b = 10.0, 10.0  # [mol/s]

    M0 = np.array([0.05, 0.05, 0.05], dtype=np.float64)  # [kg]
    tr = wirsum_1997_strict_newton(transport_residual, transport_jacobian, M0, max_iter=20)
    J_tr = transport_jacobian(tr.x)
    dx_tr = np.linalg.solve(J_tr, -transport_residual(M0))
    transport = {
        **classify_newton_stop(tr, model_is_full_rate=True),
        "rung": "isothermal_solid_transport",
        "units": {"state": "kg", "residual": "kg/s"},
        "backward_error_at_root": linear_backward_error(J_tr, tr.x, tr.F),
        "newton_step_backward_error_at_x0": newton_step_backward_error(
            transport_jacobian(M0), transport_residual(M0), dx_tr
        ),
        "jacobian_rank": int(jacobian_diagnostics(J_tr).rank),
    }

    def _ex_res(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return tracer_exchange_residual(
            x, geo=geo, N_d_in=N_d_in_tr, N_b_in=N_b_in_tr, N_tot_d=N_tot_d, N_tot_b=N_tot_b
        )

    def _ex_jac(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return tracer_exchange_jacobian(x, geo=geo, N_tot_d=N_tot_d, N_tot_b=N_tot_b)

    x_ex0 = np.array([0.2, 0.8], dtype=np.float64)
    exact = exact_tracer_exchange(
        geo=geo, N_d_in=N_d_in_tr, N_b_in=N_b_in_tr, N_tot_d=N_tot_d, N_tot_b=N_tot_b
    )
    raw = wirsum_1997_strict_newton(_ex_res, _ex_jac, x_ex0, merit_kind="raw_l2", max_iter=20)
    scaled = wirsum_1997_strict_newton(
        _ex_res,
        _ex_jac,
        x_ex0,
        merit_kind="scaled_rms",
        scale=np.array([1.0, 1.0], dtype=np.float64),
        max_iter=20,
    )
    exchange = {
        "rung": "isothermal_two_phase_exchange",
        "units": {"state": "mol/s", "residual": "mol/s"},
        "exact_N_mol_s": exact.tolist(),
        "raw_l2": classify_newton_stop(raw, model_is_full_rate=True),
        "scaled_rms": classify_newton_stop(scaled, model_is_full_rate=True),
        "max_abs_raw_vs_exact": float(np.max(np.abs(raw.x - exact))),
        "F_xx_max_abs": float(np.max(np.abs(_ex_res(raw.x) - _ex_res(raw.x)))),
        "dir_deriv_eps": {
            "1e-3": directional_derivative_error(_ex_res, _ex_jac, raw.x, np.array([1.0, -0.3]), 1e-3),
            "1e-5": directional_derivative_error(_ex_res, _ex_jac, raw.x, np.array([1.0, -0.3]), 1e-5),
            "1e-7": directional_derivative_error(_ex_res, _ex_jac, raw.x, np.array([1.0, -0.3]), 1e-7),
        },
    }

    N_in_d, N_in_b = default_r12_inlets()
    guesses = r12_initial_guesses(N_in_d, N_in_b)

    def _make_r12(scale: float) -> tuple[Any, Any]:
        def _res(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return r12_exchange_residual(
                x, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, rate_scale=scale
            )

        def _jac(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return r12_exchange_jacobian(
                x, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, rate_scale=scale
            )

        return _res, _jac

    continuation: dict[str, Any] = {}
    x_cont = guesses["separated_inlet"]
    for scale in (0.0, 1.0e-12, 1.0e-8, 1.0e-4, 1.0):
        res_fn, jac_fn = _make_r12(scale)
        sol = wirsum_1997_strict_newton(res_fn, jac_fn, x_cont, max_iter=40, lambda_min=1.0 / 4096.0)
        key = "scale_0" if scale == 0.0 else f"scale_{scale:.0e}"
        continuation[key] = {
            **classify_newton_stop(sol, model_is_full_rate=abs(scale - 1.0) < 1.0e-15),
            "rate_scale": float(scale),
        }
        if sol.converged:
            x_cont = sol.x

    uniqueness: dict[str, Any] = {}
    probe_scale = 1.0e-8
    res_u, jac_u = _make_r12(probe_scale)
    xs: list[npt.NDArray[np.float64]] = []
    for name, x0 in guesses.items():
        sol = wirsum_1997_strict_newton(res_u, jac_u, x0, max_iter=40)
        uniqueness[name] = classify_newton_stop(sol, model_is_full_rate=False)
        uniqueness[name]["rate_scale"] = probe_scale
        if sol.converged:
            xs.append(sol.x)
    if len(xs) >= 2:
        uniqueness["max_abs_x_spread"] = float(np.max(np.abs(xs[0] - xs[1])))
        if len(xs) == 3:
            uniqueness["max_abs_x_spread"] = float(
                max(np.max(np.abs(a - b)) for a, b in ((xs[0], xs[1]), (xs[0], xs[2]), (xs[1], xs[2])))
            )
    else:
        uniqueness["max_abs_x_spread"] = None

    res1, jac1 = _make_r12(1.0)
    x_full = x_cont if continuation["scale_1e+00"]["converged_newton"] else guesses["separated_inlet"]
    sol_pert = wirsum_1997_strict_newton(
        res1, jac1, x_full * 1.01 + 1.0e-8, max_iter=40
    )
    uniqueness_full = {
        "from_continuation": {
            "converged_newton": continuation["scale_1e+00"]["converged_newton"],
            "physical_convergence": continuation["scale_1e+00"]["physical_convergence"],
            "failure_reason": continuation["scale_1e+00"]["failure_reason"],
            "residual_l2": continuation["scale_1e+00"]["residual_l2"],
        },
        "from_perturbation": classify_newton_stop(sol_pert, model_is_full_rate=True),
        "max_abs_x_vs_continuation": (
            float(np.max(np.abs(sol_pert.x - x_full))) if sol_pert.converged else None
        ),
        "direct_separated_at_full_rate": classify_newton_stop(
            wirsum_1997_strict_newton(res1, jac1, guesses["separated_inlet"], max_iter=20),
            model_is_full_rate=True,
        ),
    }
    x_fd = np.array([0.25, 1.0e-4, 0.05, 10.0, 0.05, 0.10, 0.02, 10.0], dtype=np.float64)
    dir_r12 = {
        "1e-3": directional_derivative_error(res1, jac1, x_fd, np.ones(8), 1e-3),
        "1e-5": directional_derivative_error(res1, jac1, x_fd, np.ones(8), 1e-5),
        "1e-7": directional_derivative_error(res1, jac1, x_fd, np.ones(8), 1e-7),
    }

    payload = {
        "source": "Wirsum alignment stage E stiffness ladder",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "implemented_rungs": list(implemented_rungs()),
        "not_implemented_rungs": sorted(_UNIMPLEMENTED),
        "units": {"gas": "mol/s", "solid": "kg", "T": "K", "P": "Pa", "volume": "m3"},
        "geometry": {
            "T_K": geo.T,
            "P_Pa": geo.P,
            "K_bd_1_per_s": geo.K_bd,
            "V_b_m3": geo.V_b,
            "V_d_m3": geo.V_d,
        },
        "transport": transport,
        "exchange": exchange,
        "r12": {
            "rung": "isothermal_r12_exchange",
            "rate": "production rate_R12 Hamel 5.38",
            "not": ["fastox", "extent_limiter", "physical_convergence_on_continuation"],
            "inlet_dense_mol_s": N_in_d.tolist(),
            "inlet_bubble_mol_s": N_in_b.tolist(),
            "continuation": continuation,
            "uniqueness_at_scale_1e-8": uniqueness,
            "uniqueness_at_full_rate": uniqueness_full,
            "dir_deriv_at_overlap": dir_r12,
        },
        "note": "rate_scale<1 是延拓，不是 Hamel 模型。能量不含固体瞬时 cp 积分。",
    }
    later = build_full_gas_char_energy_section()
    payload.update(later)
    return payload
