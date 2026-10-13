"""F_E6+R10+R11 数学根上打开 F_E7 自由板均相。不经 ``Reactor``。

未知量仍是 43。床内 R1–R8/R10/R11/R12 保持完整。自由板是单槽 CSTR：
``N_fb = N_d + N_b + ξ_fb``，``ξ_fb`` 只用均相 R5/R6/R7/R8/R12，气泡动力学，
``ρ_ash=0``（R8 关闭）。不含 R10/R11，也不搬生产 ``_freeboard_reaction_step``。

延拓 ``F(x;s)=F_{E6+R10+R11}(x)+s(F_{E7}-F_{E6+R10+R11})``。
混合 ``N_fb`` 上的燃料 ∩ O₂ 用 Startwert 只改 ``x_0``。残差不 clip。

浓度与速率用状态 ``T=x[-1]``（``q_at_T``），不是冻结的 geo 1100 K。
不升 F_E8（``K_eject_in_this_P``），不开放 ``hamel_outer_refresh``。

Source: Hamel (1999) Eq. 2.2, 5.36–5.41; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_coupled_f_energy import q_at_T
from src.solvers.wirsum_1997_named_f import (
    N_SPECIES_BLOCK,
    NamedStageQ,
    energy_row_named,
    fb_gas_source,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
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
from src.solvers.wirsum_1997_named_f_e4 import residual_named_e4
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q
from src.solvers.wirsum_1997_named_f_r11 import (
    freeze_e6_r10_r11_q,
    jacobian_r10_plus_scaled_r11,
    solve_named_r10_then_r11,
)
from src.solvers.wirsum_1997_named_f_r12 import (
    NU_R12_MAJOR,
    _conc_and_dN,
    _r12_r_and_dC,
    named_r12_scaled_rank,
)
from src.solvers.wirsum_1997_named_f_r12_r5 import NU_R5_MAJOR, _r5_bubble_r_and_dC
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import NU_R6_MAJOR, _r6_r_and_dC
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import NU_R7_MAJOR, _r7_r_and_dC
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_named_f_registry import f_name, omitted_named
from src.solvers.wirsum_1997_startwert_e7 import apply_named_e7_fb_startwert
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, N_MAJ
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON | frozenset({"R10", "R11"})
_IDX_TAR1: int = 40
_IDX_TAR2: int = 41
_T_FD: float = 1.0  # [K]

_E7_CONTINUATION: tuple[float, ...] = (
    1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 0.1, 1.0
)
# Startwert 后小 s 会把 N_d+N_b 的微量氧加回 N_fb，R12 右导数再次硬源。
_E7_STARTWERT_CONTINUATION: tuple[float, ...] = (1.0,)
_MAX_SHRINKS: int = 8


def freeze_e7_q(*, Q_W: float = 0.0) -> NamedStageQ:
    return freeze_stage_q(7, Q_W=float(Q_W), enabled_rxns=_ENABLED)


def residual_named_e7(
    x: npt.NDArray[np.float64],
    nq7: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """``F_E7(x; P_DAEM(x))``，含自由板均相源。"""
    return residual_named_e4(x, nq7)


def residual_r11_plus_scaled_fb(
    x: npt.NDArray[np.float64],
    nq_r11: NamedStageQ,
    nq7: NamedStageQ,
    *,
    fb_scale: float,
) -> npt.NDArray[np.float64]:
    """F_E6+R10+R11 完整速率，自由板均相乘 ``fb_scale``。不 clip。"""
    z = np.asarray(x, dtype=np.float64)
    f11 = residual_named_e4(z, nq_r11)
    s = float(fb_scale)
    if abs(s) < 1.0e-30:
        return f11
    both = residual_named_e7(z, nq7)
    return f11 + s * (both - f11)


def fb_source_jacobian(
    x: npt.NDArray[np.float64],
    nq7: NamedStageQ,
    *,
    rate_scale: float,
) -> npt.NDArray[np.float64]:
    """``-ξ_fb`` 对 ``N_fb`` 的解析块。气泡均相；C_O2=0 时保留 R12 氧右导数。

    T 列由完整残差差商覆盖。R8 在自由板 ``ρ_ash=0``，不写入。
    """
    z = np.asarray(x, dtype=np.float64)
    n = int(z.size)
    J = np.zeros((n, n), dtype=np.float64)
    s = float(rate_scale)
    if abs(s) < 1.0e-30 or int(nq7.stage) < 7:
        return J
    T = float(z[-1])
    P = float(nq7.energy.species.geo_ex.P)
    V_fb = float(nq7.energy.species.geo_fb.area) * float(
        nq7.energy.species.geo_fb.H_freeboard
    )
    on = nq7.enabled_rxns
    off = 16
    C, dCdN = _conc_and_dN(z[off:off + N_MAJ], T, P)
    blocks: list[tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] = []
    if "R12" in on:
        _r, dr = _r12_r_and_dC(T, C)
        blocks.append((NU_R12_MAJOR, dr))
    if "R5" in on:
        _r, dr = _r5_bubble_r_and_dC(T, C)
        blocks.append((NU_R5_MAJOR, dr))
    if "R6" in on:
        _r, dr = _r6_r_and_dC(T, C)
        blocks.append((NU_R6_MAJOR, dr))
    if "R7" in on:
        _r, dr = _r7_r_and_dC(T, C, rho_ash=0.0)
        if float(np.linalg.norm(dr)) > 0.0:
            blocks.append((NU_R7_MAJOR, dr))
    for nu, drdC in blocks:
        dxi_dN = (s * float(V_fb)) * (drdC @ dCdN)
        J[off:off + N_MAJ, off:off + N_MAJ] -= nu[:, None] * dxi_dN
    return J


def jacobian_r11_plus_scaled_fb(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_r10: NamedStageQ,
    nq_r11: NamedStageQ,
    nq7: NamedStageQ,
    *,
    fb_scale: float,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    J = jacobian_r10_plus_scaled_r11(z, nq6, nq_r10, nq_r11, r11_scale=1.0)
    s = float(fb_scale)
    if abs(s) >= 1.0e-30:
        J = J + fb_source_jacobian(z, nq7, rate_scale=s)
        xp = z.copy()
        xp[-1] = float(z[-1]) + _T_FD
        dF = (
            residual_r11_plus_scaled_fb(xp, nq_r11, nq7, fb_scale=s)
            - residual_r11_plus_scaled_fb(z, nq_r11, nq7, fb_scale=s)
        ) / _T_FD
        J[:, -1] = dF
    return J


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return named_in_domain(z)


def _e7_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    dx, rank, cond, _diag = named_linear_step(
        J, F, x,
        extra_diag_species=("H2", "O2", "CO", "CH4"),
        phase_offsets=(0, 7, 16),
    )
    p, _changed = bound_feasible_direction(x, dx)
    return p, rank, cond


def e7_source_snapshot(
    x: npt.NDArray[np.float64],
    nq_r11: NamedStageQ,
    nq7: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f11 = residual_named_e4(z, nq_r11)
    f7 = residual_named_e7(z, nq7)
    dF = f7 - f11
    u = unpack_named(z, nq7.stage)
    T = float(u["T"])
    nqP = outer_P_named(z, nq7)
    qT = q_at_T(nqP.energy.species, T)
    src = fb_gas_source(
        z[:N_SPECIES_BLOCK], qT, rate_scale=1.0, enabled=nq7.enabled_rxns
    )
    return {
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_fb_l2": float(np.linalg.norm(dF[16:23])),
        "dF_energy": float(dF[-1]),
        "src_O2": float(src[IDX["O2"]]),
        "src_H2": float(src[IDX["H2"]]),
        "src_CH4": float(src[IDX["CH4"]]),
        "src_CO": float(src[IDX["CO"]]),
        "N_fb_O2": float(u["N_fb"][IDX["O2"]]),
        "N_fb_H2": float(u["N_fb"][IDX["H2"]]),
        "N_fb_CH4": float(u["N_fb"][IDX["CH4"]]),
        "N_fb_CO": float(u["N_fb"][IDX["CO"]]),
        "T_K": T,
        "V_fb_m3": float(qT.geo_fb.area) * float(qT.geo_fb.H_freeboard),
    }


def _newton_at_fb_scale(
    nq6: NamedStageQ,
    nq_r10: NamedStageQ,
    nq_r11: NamedStageQ,
    nq7: NamedStageQ,
    x: npt.NDArray[np.float64],
    fb_scale: float,
    max_iter: int,
) -> Any:
    x0 = np.asarray(x, dtype=np.float64)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_r11_plus_scaled_fb(z, nq_r11, nq7, fb_scale=fb_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_r11_plus_scaled_fb(
            z, nq6, nq_r10, nq_r11, nq7, fb_scale=fb_scale
        )

    return wirsum_1997_strict_newton(
        res, jac, x0, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=_e7_newton_step,
    )


def solve_named_e7_from_r11(
    *,
    Q_W: float = 0.0,
    max_iter: int = 20,
    continuation: bool = True,
    fb_scale: float | None = None,
    apply_startwert: bool = True,
    r11: dict[str, Any] | None = None,
    r10: dict[str, Any] | None = None,
    e6: dict[str, Any] | None = None,
    advance_mode: AdvanceMode | None = None,
    s_min: float = 1.0e-16,
) -> dict[str, Any]:
    """从 F_E6+R10+R11 根做 F_E7 延拓。默认 FB Startwert + 逐根校正。"""
    nq6 = freeze_e6_q(Q_W=Q_W)
    nq_r10 = freeze_e6_r10_q(Q_W=Q_W)
    nq_r11 = freeze_e6_r10_r11_q(Q_W=Q_W)
    nq7 = freeze_e7_q(Q_W=Q_W)
    if r11 is None:
        r11 = solve_named_r10_then_r11(Q_W=Q_W, max_iter=max_iter, r10=r10, e6=e6)
    if advance_mode is None:
        advance_mode = "root_correcting" if apply_startwert else "heuristic_search"
    x_start = np.asarray(r11["x"], dtype=np.float64).copy()
    snap_root = e7_source_snapshot(x_start, nq_r11, nq7)
    sw7 = None
    if apply_startwert:
        packed = apply_named_e7_fb_startwert(x_start)
        x_start = packed["x"]
        sw7 = {
            "xi12_mol_s": packed["xi12_mol_s"],
            "xi5_mol_s": packed["xi5_mol_s"],
            "xi6_mol_s": packed["xi6_mol_s"],
        }
    snap0 = e7_source_snapshot(x_start, nq_r11, nq7)
    cert0 = certify_named_full_rate_endpoint(
        x_start, nq7, coupling="composite_P"
    )
    if fb_scale is not None:
        planned = [float(fb_scale)]
    elif continuation:
        planned = list(
            _E7_STARTWERT_CONTINUATION if apply_startwert else _E7_CONTINUATION
        )
    else:
        planned = [1.0]
    history: list[dict[str, Any]] = []
    x_cur = x_start
    last = None
    s_ok = 0.0
    branch_tracked = True
    shrinks = 0
    i = 0
    while i < len(planned):
        s = float(planned[i])
        last = _newton_at_fb_scale(nq6, nq_r10, nq_r11, nq7, x_cur, s, max_iter)
        stop_s = classify_newton_stop(last, model_is_full_rate=abs(s - 1.0) < 1.0e-15)
        J = jacobian_r11_plus_scaled_fb(
            last.x, nq6, nq_r10, nq_r11, nq7, fb_scale=s
        )
        info = named_r12_scaled_rank(J, last.F, last.x)
        passed = stage_passed(
            stop_s, last.x, advance_mode=advance_mode,
            accepted_alphas=list(last.accepted_alphas),
        )
        if abs(s - 1.0) < 1.0e-15:
            cert_s = certify_named_full_rate_endpoint(
                last.x, nq7, coupling="composite_P"
            )
            if cert_s["endpoint_certified"]:
                passed = True
            elif apply_startwert and cert0["endpoint_certified"]:
                last.x = np.asarray(x_start, dtype=np.float64)
                last.F = residual_named_e7(last.x, nq7)
                stop_s = classify_newton_stop(last, model_is_full_rate=True)
                J = jacobian_r11_plus_scaled_fb(
                    last.x, nq6, nq_r10, nq_r11, nq7, fb_scale=s
                )
                info = named_r12_scaled_rank(J, last.F, last.x)
                passed = True
        history.append(
            {
                "fb_scale": float(s),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "solver_stop_reason": stop_s["failure_reason"],
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "stage_passed": bool(passed),
                "N_fb_O2": float(last.x[16 + IDX["O2"]]),
                "N_fb_H2": float(last.x[16 + IDX["H2"]]),
                "N_fb_CH4": float(last.x[16 + IDX["CH4"]]),
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
            # Startwert 后失败若缩回 s≪1，混合约束会把微量氧加回 N_fb。
            if apply_startwert and abs(s_ok) < 1.0e-18:
                branch_tracked = False
                break
            if shrinks >= _MAX_SHRINKS:
                branch_tracked = False
                break
            s_mid = shrink_scale(s_ok, s)
            if abs(s - s_ok) < float(s_min) or abs(s_mid - s) < 1.0e-18:
                branch_tracked = False
                break
            planned.insert(i, s_mid)
            shrinks += 1
            continue
        branch_tracked = False
        i += 1

    assert last is not None
    s_last = float(history[-1]["fb_scale"]) if history else 1.0
    stop = classify_newton_stop(
        last, model_is_full_rate=abs(s_last - 1.0) < 1.0e-15
    )
    info = named_r12_scaled_rank(
        jacobian_r11_plus_scaled_fb(
            last.x, nq6, nq_r10, nq_r11, nq7, fb_scale=1.0
        ),
        last.F, last.x,
    )
    nqP = outer_P_named(last.x, nq7)
    u = unpack_named(last.x, 7)
    cert = certify_named_full_rate_endpoint(last.x, nq7, coupling="composite_P")
    omitted = omitted_named(7)
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": float(last.x[-1]),
        "T_r11_K": float(r11["T_K"]),
        "M_ash": float(u["M_ash"]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "N_fb_O2": float(u["N_fb"][IDX["O2"]]),
        "N_fb_H2": float(u["N_fb"][IDX["H2"]]),
        "N_fb_CH4": float(u["N_fb"][IDX["CH4"]]),
        "N_d_O2": float(u["N_d"][IDX["O2"]]),
        "N_b_O2": float(u["N_b"][IDX["O2"]]),
        "Q_W": float(nq7.Q_W),
        "Q_W_is_fitted_to_1100K": False,
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": f_name(7),
        "n_unknowns": n_unknowns(7),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "r11_physical_convergence": bool(r11.get("physical_convergence", False)),
        "r11_endpoint_certified": bool(r11.get("endpoint_certified", False)),
        "e7_at_r11_root": snap_root,
        "e7_after_startwert": snap0,
        "startwert_e7": sw7,
        "apply_startwert": bool(apply_startwert),
        "advance_mode": str(advance_mode),
        "branch_tracked": bool(
            branch_tracked
            and bool(history)
            and abs(float(history[-1]["fb_scale"]) - 1.0) < 1.0e-15
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
        "freeboard_reaction": True,
        "freeboard_r10_r11": False,
        "production_freeboard_step": False,
        "x_nonneg_projection": False,
        "coupling_mode": "composite_P",
        "bound_feasible_direction": True,
        "omitted_named": omitted,
        "K_eject_in_this_P": omitted.get("K_eject_in_this_P"),
    }
