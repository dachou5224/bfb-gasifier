"""F_E6+R10 数学根上接入 R11（焦油转化）。不改 ``rate_R11_*``，不经 ``Reactor``。

R1–R8/R10/R12 与 DAEM、灰、``Q_W`` 保持完整。R11 用 ``r11_scale`` 延拓：
F(x;s) = F_{E6+R10}(x) + s (F_{E6+R10+R11}(x) - F_{E6+R10}(x)).

气泡 Table 5.4：Serio 一级 ``k(T) C_tar``。悬浮：Corella ``k(T) rho_cat C_tar``，
``rho_cat = M_ash/V_d``。速率不含 C_H2O；H2O 只在计量里。
C_tar=0 时一级右导数 ``k`` 保留。不升到 F_E7 自由板。

定义域 x[:-1]≥0；越界试探缩短 α，残差不 clip。
每档用气泡相 ``ξ=Δs·N_tar`` 预测子，多出的 H2 在气泡里走反向 WGSR，不写入 ``N_d``。

Source: Hamel (1999) Table 5.4, §5.2.6, §2.1; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.kinetics.arrhenius import k_standard
from src.kinetics.tar_reactions import (
    R11_CORELLA_E,
    R11_CORELLA_K0,
    R11_SERIO_E,
    R11_SERIO_K0,
    get_lumped_tar_stoichiometry,
)
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    N_SPECIES_BLOCK,
    R10_FUEL_TYPE,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
    r10_c_tar,
    r11_nu_major,
    rho_ash_used,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_continuation import (
    AdvanceMode,
    bound_feasible_direction,
    certify_named_full_rate_endpoint,
    named_in_domain,
    named_linear_step,
    shrink_scale,
    stage_passed,
)
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q, residual_named_e6
from src.solvers.wirsum_1997_named_f_r10 import (
    freeze_e6_r10_q,
    jacobian_e6_plus_scaled_r10,
    solve_named_e6_then_r10,
)
from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_startwert_r11 import apply_scaled_bubble_r11_predictor
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON | frozenset({"R10", "R11"})
_IDX_TAR1: int = 40
_IDX_TAR2: int = 41
_IDX_ASH: int = int(N_SPECIES_BLOCK)
_T_FD: float = 1.0  # [K]
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]

_R11_CONTINUATION: tuple[float, ...] = (
    1.0e-8, 1.0e-6, 1.0e-4, 0.01, 0.02, 0.03, 0.08, 0.16, 0.32, 0.64, 1.0
)
_R11_STARTWERT_CONTINUATION: tuple[float, ...] = (
    0.08, 0.16, 0.32, 0.64, 1.0
)


def freeze_e6_r10_r11_q(*, Q_W: float = 0.0) -> NamedStageQ:
    return freeze_stage_q(6, Q_W=float(Q_W), enabled_rxns=_ENABLED)


def residual_r10_plus_scaled_r11(
    x: npt.NDArray[np.float64],
    nq_r10: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r11_scale: float,
) -> npt.NDArray[np.float64]:
    """F_E6+R10 完整速率，R11 乘 ``r11_scale``。不 clip。"""
    z = np.asarray(x, dtype=np.float64)
    f10 = residual_named_e6(z, nq_r10)
    s = float(r11_scale)
    if abs(s) < 1.0e-30:
        return f10
    both = residual_named_e6(z, nq_all)
    return f10 + s * (both - f10)


def r11_analytic_jacobian(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """R11 对 N/N_tar/M_ash 的解析块。C_tar=0 时保留一级右导数。"""
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    if "R11" not in nq.enabled_rxns or int(nq.stage) < 6:
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
    rho_cat = rho_ash_used(nq, u)
    k_b = float(k_standard(R11_SERIO_K0, R11_SERIO_E, T))
    k_d_per_rho = float(k_standard(R11_CORELLA_K0, R11_CORELLA_E, T))
    k_d = k_d_per_rho * float(rho_cat)
    nu = r11_nu_major()
    sto = get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)
    nu_t1 = float(sto.get("TAR1", 0.0))
    nu_t2 = float(sto.get("TAR2", 0.0))
    n_gas = float(np.sum(u["N_d"]) + np.sum(u["N_b"]))
    n_tot = n_gas + n_tar
    c_tot = float(P) / (Rg * float(T))
    if n_tot < _NTOT_FLOOR:
        return J
    dCtar_dnt = (n_gas / (n_tot * n_tot)) * c_tot
    dCtar_dNj = -(n_tar / (n_tot * n_tot)) * c_tot
    ks = ((0, V_d, k_d), (7, V_b, k_b))
    for off, V, k_r in ks:
        dr_dCtar = float(k_r)
        dxi_gas = (s * float(V)) * dr_dCtar * dCtar_dNj
        dxi_tar = (s * float(V)) * dr_dCtar * dCtar_dnt
        other = 7 if off == 0 else 0
        dxi_this = np.full(N_MAJ, dxi_gas, dtype=np.float64)
        dxi_other = np.full(N_MAJ, dxi_gas, dtype=np.float64)
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
    dxi_d_dM = s * k_d_per_rho * float(C_tar)
    J[0:N_MAJ, _IDX_ASH] += nu * dxi_d_dM
    J[_IDX_TAR1, _IDX_ASH] -= nu_t1 * dxi_d_dM
    J[_IDX_TAR2, _IDX_ASH] -= nu_t2 * dxi_d_dM
    return J


def jacobian_r10_plus_scaled_r11(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_r10: NamedStageQ,
    nq_all: NamedStageQ,
    *,
    r11_scale: float,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    J = jacobian_e6_plus_scaled_r10(z, nq6, nq_r10, r10_scale=1.0)
    s = float(r11_scale)
    if abs(s) >= 1.0e-30:
        J = J + r11_analytic_jacobian(z, nq_all, rate_scale=s)
        xp = z.copy()
        xp[-1] = float(z[-1]) + _T_FD
        dF = (
            residual_r10_plus_scaled_r11(xp, nq_r10, nq_all, r11_scale=s)
            - residual_r10_plus_scaled_r11(z, nq_r10, nq_all, r11_scale=s)
        ) / _T_FD
        J[:, -1] = dF
    return J


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return named_in_domain(z)


def _r11_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    dx, rank, cond, _diag = named_linear_step(
        J, F, x, extra_diag_species=("H2", "O2", "CO", "CH4"),
    )
    p, _changed = bound_feasible_direction(x, dx)
    return p, rank, cond


def r11_source_snapshot(
    x: npt.NDArray[np.float64],
    nq_r10: NamedStageQ,
    nq_all: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f10 = residual_named_e6(z, nq_r10)
    both = residual_named_e6(z, nq_all)
    dF = both - f10
    u = unpack_named(z, 6)
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_H2O_d": float(dF[IDX["H2O"]]),
        "dF_H2O_b": float(dF[7 + IDX["H2O"]]),
        "dF_CO_d": float(dF[IDX["CO"]]),
        "dF_H2_d": float(dF[IDX["H2"]]),
        "dF_CH4_d": float(dF[IDX["CH4"]]),
        "dF_tar1": float(dF[_IDX_TAR1]),
        "dF_tar2": float(dF[_IDX_TAR2]),
        "dF_energy": float(dF[-1]),
        "N_d_H2O": float(u["N_d"][IDX["H2O"]]),
        "N_b_H2O": float(u["N_b"][IDX["H2O"]]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "T_K": float(u["T"]),
    }


def _newton_at_r11_scale(
    nq6: NamedStageQ,
    nq_r10: NamedStageQ,
    nq_all: NamedStageQ,
    x: npt.NDArray[np.float64],
    r11_scale: float,
    max_iter: int,
) -> Any:
    x0 = np.asarray(x, dtype=np.float64)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_r10_plus_scaled_r11(z, nq_r10, nq_all, r11_scale=r11_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_r10_plus_scaled_r11(
            z, nq6, nq_r10, nq_all, r11_scale=r11_scale
        )

    return wirsum_1997_strict_newton(
        res, jac, x0, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=_r11_newton_step,
    )


def solve_named_r10_then_r11(
    *,
    Q_W: float = 0.0,
    max_iter: int = 20,
    continuation: bool = True,
    r11_scale: float | None = None,
    apply_startwert: bool = True,
    r10: dict[str, Any] | None = None,
    e6: dict[str, Any] | None = None,
    advance_mode: AdvanceMode | None = None,
    s_min: float = 1.0e-16,
) -> dict[str, Any]:
    """从 F_E6+R10 根做 R11 延拓。默认 Startwert + 逐根校正。"""
    nq6 = freeze_e6_q(Q_W=Q_W)
    nq_r10 = freeze_e6_r10_q(Q_W=Q_W)
    nq_all = freeze_e6_r10_r11_q(Q_W=Q_W)
    if r10 is None:
        r10 = solve_named_e6_then_r10(Q_W=Q_W, max_iter=max_iter, e6=e6)
    if advance_mode is None:
        advance_mode = "root_correcting" if apply_startwert else "heuristic_search"
    x_start = np.asarray(r10["x"], dtype=np.float64).copy()
    snap_root = r11_source_snapshot(x_start, nq_r10, nq_all)
    sw11 = None
    if r11_scale is not None:
        planned = [float(r11_scale)]
    elif continuation:
        planned = list(
            _R11_STARTWERT_CONTINUATION if apply_startwert else _R11_CONTINUATION
        )
    else:
        planned = [1.0]
    if apply_startwert:
        packed = apply_scaled_bubble_r11_predictor(x_start, ds=float(planned[0]))
        sw11 = {
            "predictor": {k: v for k, v in packed.items() if k != "x"},
            "dense": {"xi11_mol_s": 0.0},
            "bubble": {
                "xi11_mol_s": float(packed["xi11_mol_s"]),
                "N_tar1_after": float(packed["x"][_IDX_TAR1]),
                "N_tar2_after": float(packed["x"][_IDX_TAR2]),
            },
        }
        snap0 = r11_source_snapshot(packed["x"], nq_r10, nq_all)
    else:
        snap0 = r11_source_snapshot(x_start, nq_r10, nq_all)
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    s_ok = 0.0
    branch_tracked = True
    i = 0
    while i < len(planned):
        s = float(planned[i])
        x_newton = x_cur
        if apply_startwert and s > s_ok:
            pred = apply_scaled_bubble_r11_predictor(x_cur, ds=s - s_ok)
            x_newton = pred["x"]
        last = _newton_at_r11_scale(nq6, nq_r10, nq_all, x_newton, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_r10_plus_scaled_r11(
            last.x, nq6, nq_r10, nq_all, r11_scale=s
        )
        info = named_r12_scaled_rank(J, last.F, last.x)
        passed = stage_passed(
            stop_s, last.x, advance_mode=advance_mode,
            accepted_alphas=list(last.accepted_alphas),
        )
        history.append(
            {
                "r11_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "solver_stop_reason": stop_s["failure_reason"],
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "stage_passed": bool(passed),
                "N_d_H2O": float(last.x[IDX["H2O"]]),
                "N_b_H2O": float(last.x[7 + IDX["H2O"]]),
                "N_tar1": float(last.x[_IDX_TAR1]),
                "N_tar2": float(last.x[_IDX_TAR2]),
            }
        )
        if passed:
            x_cur = np.asarray(last.x, dtype=np.float64)
            s_ok = s
            i += 1
            continue
        if advance_mode == "root_correcting":
            s_mid = shrink_scale(s_ok, s)
            if abs(s - s_ok) < float(s_min) or abs(s_mid - s) < 1.0e-18:
                branch_tracked = False
                break
            planned.insert(i, s_mid)
            continue
        branch_tracked = False
        i += 1

    assert last is not None
    s_last = float(history[-1]["r11_scale"]) if history else 1.0
    stop = classify_newton_stop(
        last, model_is_full_rate=abs(s_last - 1.0) < 1.0e-15
    )
    info = named_r12_scaled_rank(
        jacobian_r10_plus_scaled_r11(
            last.x, nq6, nq_r10, nq_all, r11_scale=1.0
        ),
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
        "T_r10_K": float(r10["T_K"]),
        "M_ash": float(u["M_ash"]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "N_d_H2O": float(u["N_d"][IDX["H2O"]]),
        "N_b_H2O": float(u["N_b"][IDX["H2O"]]),
        "N_d_CH4": float(u["N_d"][IDX["CH4"]]),
        "N_d_H2": float(u["N_d"][IDX["H2"]]),
        "Q_W": float(nq_all.Q_W),
        "Q_W_is_fitted_to_1100K": False,
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": "F_E6+R10+R11",
        "n_unknowns": n_unknowns(6),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "r10_physical_convergence": bool(r10.get("physical_convergence", False)),
        "r10_endpoint_certified": bool(r10.get("endpoint_certified", False)),
        "r11_at_r10_root": snap_root,
        "r11_after_startwert": snap0,
        "startwert_r11": sw11,
        "apply_startwert": bool(apply_startwert),
        "advance_mode": str(advance_mode),
        "branch_tracked": bool(
            branch_tracked
            and bool(history)
            and abs(float(history[-1]["r11_scale"]) - 1.0) < 1.0e-15
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
        "r11_tar_conversion": True,
        "x_nonneg_projection": False,
        "coupling_mode": "composite_P",
        "bound_feasible_direction": True,
        "freeboard_reaction": False,
    }
