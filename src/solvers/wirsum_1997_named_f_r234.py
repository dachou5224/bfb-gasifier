"""均相数学根上接入 R2/R3/R4（F_E1 炭捆包）。不改 ``rate_R2/R3/R4``，不经 ``Reactor``。

均相 R12/R5/R8/R6/R7 保持完整速率。R2–R4 用 ``r234_scale`` 延拓，只进密相。
R1 已在 F_E0，这里不重复。

R2：``R_C = A_C k_ch,R2 C_H2O``，Eq. 5.24–5.25，``k_hobbs``。
R3：加氢气化 Eq. 5.34–5.35，``k_hobbs``；均相根上 H2 已尽，源 ~10^{-13}。
R4：Weeda L-H Boudouard Eq. 5.26–5.33，``k_standard``；禁止改 k_b/k_c 符号。

均相根上 R2 ~0.77 mol/s、R4 ~0.018 mol/s。不做气体同相 Startwert。

Source: Hamel (1999) Eq. 5.24–5.35; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import MW_C_KG_PER_MOL
from src.kinetics.arrhenius import k_hobbs, k_standard
from src.kinetics.char_reactions import (
    R2_E_Rg,
    R2_k0,
    R3_E_Rg,
    R3_k0,
    R4_M_C,
    R4_kb_E,
    R4_kb_k0,
    R4_kc_E,
    R4_kc_k0,
    R4_kf_E,
    R4_kf_k0,
    rate_R2,
    rate_R3,
    rate_R4,
)
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import IDX_COARSE, IDX_FINE
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
)
from src.solvers.wirsum_1997_named_f_r12 import _conc_and_dN, named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import (
    freeze_r12_r5_r8_r6_q,
    named_r12_r5_r8_r6_newton_step,
)
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import (
    freeze_r12_r5_r8_r6_r7_q,
    jacobian_r12_r5_r8_r6_plus_scaled_r7,
    solve_named_r12_r5_r8_r6_then_r7,
)
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ, char_area_m2
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

HOMOG_ON: frozenset[str] = frozenset({"R12", "R5", "R8", "R6", "R7"})
R234_ON: frozenset[str] = frozenset({"R2", "R3", "R4"})
_M_IDX: tuple[int, int] = (14, 15)  # M_fine, M_coarse
_T_FLOOR: float = 300.0  # [K]
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]
_MW_C: float = float(MW_C_KG_PER_MOL)  # [kg/mol]

NU_R2_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R2_MAJOR[IDX["CO"]] = 1.0
NU_R2_MAJOR[IDX["H2"]] = 1.0
NU_R2_MAJOR[IDX["H2O"]] = -1.0
NU_R3_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R3_MAJOR[IDX["CH4"]] = 1.0
NU_R3_MAJOR[IDX["H2"]] = -2.0
NU_R4_MAJOR: npt.NDArray[np.float64] = np.zeros(N_MAJ, dtype=np.float64)
NU_R4_MAJOR[IDX["CO"]] = 2.0
NU_R4_MAJOR[IDX["CO2"]] = -1.0

_R234_CONTINUATION: tuple[float, ...] = (
    0.0, 1.0e-16, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0
)


def freeze_homog_r234_q() -> NamedStageQ:
    return freeze_stage_q(2, enabled_rxns=HOMOG_ON | R234_ON)


def residual_homog_plus_scaled_r234(
    x: npt.NDArray[np.float64],
    nq_h: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r234_scale: float,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    fh = coupled_residual_named(z, nq_h, rate_scale=1.0)
    s = float(r234_scale)
    if abs(s) < 1.0e-30:
        return fh
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    return fh + s * (both - fh)


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


def _add_extent(
    J: npt.NDArray[np.float64],
    nu: npt.NDArray[np.float64],
    de_dN: npt.NDArray[np.float64],
    de_dM: float,
    de_dT: float,
    m_idx: int,
) -> None:
    J[0:N_MAJ, 0:N_MAJ] += nu[:, None] * de_dN
    J[0:N_MAJ, m_idx] += nu * float(de_dM)
    J[0:N_MAJ, -1] += nu * float(de_dT)
    J[m_idx, 0:N_MAJ] -= _MW_C * de_dN
    J[m_idx, m_idx] -= _MW_C * float(de_dM)
    J[m_idx, -1] -= _MW_C * float(de_dT)


def _r2_k(T: float) -> float:
    t = max(float(T), _T_FLOOR)
    return float(k_hobbs(R2_k0, float(R2_E_Rg) * Rg, t))


def _r3_k(T: float) -> float:
    t = max(float(T), _T_FLOOR)
    return float(k_hobbs(R3_k0, float(R3_E_Rg) * Rg, t))


def r234_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """R2/R3/R4 对密相 N、两档 M、T 的解析块。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    s = float(rate_scale)
    if abs(s) < 1.0e-30:
        return J
    T = float(z[-1])
    t = max(T, _T_FLOOR)
    P = float(nq.energy.species.geo_ex.P)
    rho = float(nq.energy.species.rho_p)
    d_p = nq.energy.species.d_p
    C, dCdN = _conc_and_dN(z[0:N_MAJ], T, P)
    y, dydN = _y_and_dydN(z[0:N_MAJ])
    k2 = _r2_k(t)
    k3 = _r3_k(t)
    dlnk2 = 1.0 / t + float(R2_E_Rg) / (t * t)
    dlnk3 = 1.0 / t + float(R3_E_Rg) / (t * t)
    h2o = float(C[IDX["H2O"]])
    h2 = float(C[IDX["H2"]])
    p_co2 = float(y[IDX["CO2"]]) * P
    p_co = float(y[IDX["CO"]]) * P
    kf = float(k_standard(R4_kf_k0, R4_kf_E, t))
    kb = float(k_standard(R4_kb_k0, R4_kb_E, t))
    kc = float(k_standard(R4_kc_k0, R4_kc_E, t))
    denom = 1.0 + kb * max(p_co2, 0.0) + kc * max(p_co, 0.0)
    r_star = float(rate_R4(t, p_co2, p_co))
    d_rstar_dp_co2 = kf * (1.0 + kc * max(p_co, 0.0)) / (denom * denom)
    d_rstar_dp_co = -kf * max(p_co2, 0.0) * kc / (denom * denom)
    dp_co2_dN = P * dydN[IDX["CO2"]]
    dp_co_dN = P * dydN[IDX["CO"]]
    d_rstar_dN = d_rstar_dp_co2 * dp_co2_dN + d_rstar_dp_co * dp_co_dN
    dkf = kf * (float(R4_kf_E) / (Rg * t * t))
    dkb = kb * (float(R4_kb_E) / (Rg * t * t))
    dkc = kc * (float(R4_kc_E) / (Rg * t * t))
    d_rstar_dT = (
        dkf * max(p_co2, 0.0) * denom
        - kf * max(p_co2, 0.0) * (dkb * max(p_co2, 0.0) + dkc * max(p_co, 0.0))
    ) / (denom * denom)
    for bin_i, m_idx in ((IDX_FINE, _M_IDX[0]), (IDX_COARSE, _M_IDX[1])):
        M = float(z[m_idx])  # [kg]
        dp = float(d_p[bin_i])
        area = char_area_m2(M, rho_s=rho, d_p=dp)
        darea_dM = 6.0 / max(float(rho) * dp, 1.0e-30)
        # R2：C_H2O=0 时保留右导数 k*area。
        de2_dC = np.zeros(N_MAJ, dtype=np.float64)
        de2_dC[IDX["H2O"]] = s * k2 * area
        de2_dN = de2_dC @ dCdN
        e2 = s * float(rate_R2(t, max(h2o, 0.0), dp, 1.0)) * area
        de2_dM = s * k2 * max(h2o, 0.0) * darea_dM
        de2_dT = e2 * (dlnk2 - 1.0 / t) if h2o > 0.0 else 0.0
        _add_extent(J, NU_R2_MAJOR, de2_dN, de2_dM, de2_dT, m_idx)
        # R3：C_H2=0 时保留右导数。
        de3_dC = np.zeros(N_MAJ, dtype=np.float64)
        de3_dC[IDX["H2"]] = s * k3 * area
        de3_dN = de3_dC @ dCdN
        e3 = s * float(rate_R3(t, max(h2, 0.0), dp, 1.0)) * area
        de3_dM = s * k3 * max(h2, 0.0) * darea_dM
        de3_dT = e3 * (dlnk3 - 1.0 / t) if h2 > 0.0 else 0.0
        _add_extent(J, NU_R3_MAJOR, de3_dN, de3_dM, de3_dT, m_idx)
        # R4：extent = r* M / M_C。
        de4_dN = (s * M / R4_M_C) * d_rstar_dN
        de4_dM = s * r_star / R4_M_C
        de4_dT = (s * M / R4_M_C) * d_rstar_dT
        _add_extent(J, NU_R4_MAJOR, de4_dN, de4_dM, de4_dT, m_idx)
    return J


def jacobian_homog_plus_scaled_r234(
    x: npt.NDArray[np.float64],
    nq_h: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r234_scale: float,
) -> npt.NDArray[np.float64]:
    Jh = jacobian_r12_r5_r8_r6_plus_scaled_r7(
        x, freeze_r12_r5_r8_r6_q(), nq_h, r7_scale=1.0
    )
    s = float(r234_scale)
    if abs(s) < 1.0e-30:
        return Jh
    return Jh + s * r234_analytic_jacobian(x, nq_all, rate_scale=1.0)


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def _newton_at_r234_scale(
    nq_h: NamedStageQ,
    nq_all: NamedStageQ,
    x: npt.NDArray[np.float64],
    r234_scale: float,
    max_iter: int,
) -> Any:
    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_homog_plus_scaled_r234(
            z, nq_h, nq_all, r234_scale=r234_scale
        )

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_homog_plus_scaled_r234(
            z, nq_h, nq_all, r234_scale=r234_scale
        )

    return wirsum_1997_strict_newton(
        res, jac, x, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_r8_r6_newton_step,
    )


def r234_source_snapshot(
    x: npt.NDArray[np.float64],
    nq_h: NamedStageQ,
    nq_all: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    fh = coupled_residual_named(z, nq_h, rate_scale=1.0)
    both = coupled_residual_named(z, nq_all, rate_scale=1.0)
    dF = both - fh
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_H2_d": float(dF[IDX["H2"]]),
        "dF_H2O_d": float(dF[IDX["H2O"]]),
        "dF_CO2_d": float(dF[IDX["CO2"]]),
        "dF_CH4_d": float(dF[IDX["CH4"]]),
        "dF_M_fine": float(dF[14]),
        "dF_M_coarse": float(dF[15]),
        "dF_energy": float(dF[-1]),
        "N_d_H2": float(z[IDX["H2"]]),
        "N_d_O2": float(z[IDX["O2"]]),
        "N_d_H2O": float(z[IDX["H2O"]]),
        "N_d_CO": float(z[IDX["CO"]]),
        "M_fine": float(z[14]),
        "M_coarse": float(z[15]),
        "T_K": float(z[-1]),
    }


def solve_named_homog_then_r234(
    *,
    max_iter: int = 20,
    continuation: bool = True,
    r234_scale: float | None = None,
    homog: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E0+均相数学根做 R2/R3/R4 延拓。"""
    nq_h = freeze_r12_r5_r8_r6_r7_q()
    nq_all = freeze_homog_r234_q()
    if homog is None:
        homog = solve_named_r12_r5_r8_r6_then_r7(max_iter=max_iter)
    x_start = np.asarray(homog["x"], dtype=np.float64).copy()
    snap0 = r234_source_snapshot(x_start, nq_h, nq_all)
    scales = (float(r234_scale),) if r234_scale is not None else (
        _R234_CONTINUATION if continuation else (1.0,)
    )
    n = n_unknowns(nq_all.stage)
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    for s in scales:
        last = _newton_at_r234_scale(nq_h, nq_all, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_homog_plus_scaled_r234(last.x, nq_h, nq_all, r234_scale=s)
        info = named_r12_scaled_rank(J, last.F, last.x)
        history.append(
            {
                "r234_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "N_d_H2": float(last.x[IDX["H2"]]),
                "N_d_O2": float(last.x[IDX["O2"]]),
                "N_d_CO": float(last.x[IDX["CO"]]),
                "N_d_H2O": float(last.x[IDX["H2O"]]),
                "N_d_CO2": float(last.x[IDX["CO2"]]),
                "M_fine": float(last.x[14]),
                "M_coarse": float(last.x[15]),
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
        "T_homog_K": float(homog["T_K"]),
        "energy_W": float(energy_row_named(last.x, nq_all)),
        "f_name": "F_E0+R12+R5+R8+R6+R7+R2+R3+R4",
        "enabled_rxns": ("R12", "R5", "R8", "R6", "R7", "R2", "R3", "R4"),
        "homog_physical_convergence": bool(homog["physical_convergence"]),
        "r234_at_homog_root": snap0,
        "continuation": history,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
    }
