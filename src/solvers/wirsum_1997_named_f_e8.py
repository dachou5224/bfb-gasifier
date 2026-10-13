"""F_E7 猜测上打开 F_E8：有限水力 ``P(x;b)``、两相焦油与 ``K_eject``。不经 ``Reactor``。

未知量 45。E7 的 43 维根只是起点，不是尺度 0 的根：气泡焦油初值为 0 时
新交换行和各相焦油浓度已经把残差抬到约 0.05 mol/s。先求尺度 0，再插值水力。
延拓 ``F(x;s)=F(x;(1-s)P_{E7}+s P_{E8})``。``(hydro,fate)=(1,0)`` 不是完整 E8。
旋风 ``k`` 冻结。未判定不是逃逸。不开放 ``hamel_outer_refresh``。

Source: Hamel (1999) Eq. 2.2, 2.4–2.5, 3.13, 3.50, 3.72–3.73, 5.59;
Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import time

import numpy as np
import numpy.typing as npt

from src.physics.wirsum_1997_freeboard import (
    drag_equilibrium_speed_m_s,
    drag_z_for_spec,
    resolve_freeboard_model,
    FREEBOARD_MODEL_HAMEL_EXTENDED,
)
from src.solvers.wirsum_1997_named_f_p2 import bed_face_exchange, gas_velocity_and_density
from src.solvers.wirsum_1997_named_f import (
    GEO_FB_RHO_G_T_REF_K,
    NamedStageQ,
    apply_tar_reaction_transport_rows,
    char_face_flows_kg_s,
    coupled_residual_named,
    inventories_from_sum_fraction,
    linear_tar_row_map,
    original_residual_from_tar_rows,
    tar_consumption_fractions,
    tar_sum_and_fraction,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
    overlay_eject_from_p,
    inner_hydro_frozen_residual,
    refresh_pyrolysis_named,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_continuation import (
    AdvanceMode,
    bound_feasible_direction,
    certify_named_full_rate_endpoint,
    diagnose_linear_step,
    ENDPOINT_L2,
    named_endpoint_blocks,
    named_in_domain,
    named_linear_step,
    shrink_scale,
    stage_passed,
)
from src.solvers.wirsum_1997_named_f_e4 import residual_named_e4
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q
from src.solvers.wirsum_1997_named_f_e7 import (
    freeze_e7_q,
    residual_named_e7,
    solve_named_e7_from_r11,
)
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q
from src.solvers.wirsum_1997_named_f_r11 import freeze_e6_r10_r11_q
from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_named_f_registry import f_name, n_unknowns, omitted_named
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, fd_jacobian
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import (
    WirsumStrictSolveResult,
    wirsum_1997_strict_newton,
)

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON | frozenset({"R10", "R11"})
_T_FD: float = 1.0  # [K]
_E8_HYDRO_CONTINUATION: tuple[float, ...] = (
    0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
)
_MAX_SHRINKS: int = 48
_N_PICARD: int = 6
_TRACE_ABS: float = 1.0e-8  # [与 x 同单位] 低于此用微量步长
_TRACE_EPS_FLOOR: float = 1.0e-20  # [与 x 同单位] 微量差商下限


def freeze_e8_q(*, Q_W: float = 0.0) -> NamedStageQ:
    return freeze_stage_q(8, Q_W=float(Q_W), enabled_rxns=_ENABLED)


def as_e8_state(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """E7 的 43 维状态补上两个气泡焦油未知量，初值 0。已是 45 维则原样返回。

    悬浮相焦油沿用 E6 的两个槽。热解给料进悬浮相，气泡只靠交换和反应。
    """
    z = np.asarray(x, dtype=np.float64).reshape(-1)
    n8 = n_unknowns(8)
    if int(z.size) == n8:
        return z.copy()
    if int(z.size) == n_unknowns(7):
        return np.concatenate(
            [z[:-1], np.zeros(2, dtype=np.float64), z[-1:]]
        )
    raise ValueError(f"F_E8 expected {n8} or {n_unknowns(7)}, got {z.size}")


def residual_named_e8(
    x: npt.NDArray[np.float64],
    nq8: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """``F_E8(x; P_DAEM,K_eject(x))``。焦油两相衡算，交换与主气体同号。"""
    return residual_named_e4(as_e8_state(x), nq8)


def _lerp(a: float, b: float, s: float) -> float:
    return (1.0 - s) * float(a) + s * float(b)


def blend_hydro_named_q(
    q7: NamedStageQ,
    q8: NamedStageQ,
    *,
    hydro_scale: float,
    fate_scale: float | None = None,
) -> NamedStageQ:
    """``q`` 分字段插值。``hydro_scale`` 管气速/密度/``K_bd``；``fate_scale`` 管抛射。

    旋风 ``k`` 不混合。``fate_scale is None`` 时与 ``hydro_scale`` 同步。
    """
    s_g = float(hydro_scale)
    s_f = s_g if fate_scale is None else float(fate_scale)
    if s_g <= 0.0 and s_f <= 0.0:
        return overlay_eject_from_p(q8, q7)
    if s_g >= 1.0 and s_f >= 1.0:
        return q8
    sp7 = q7.energy.species
    sp8 = q8.energy.species
    g7 = sp7.geo_fb
    g8 = sp8.geo_fb
    geo_fb = replace(
        g8,
        u_g=_lerp(g7.u_g, g8.u_g, s_g),
        rho_g=_lerp(g7.rho_g, g8.rho_g, s_g),
        u_b=_lerp(g7.u_b, g8.u_b, s_g),
    )
    geo_ex = replace(
        sp8.geo_ex,
        K_bd=_lerp(sp7.geo_ex.K_bd, sp8.geo_ex.K_bd, s_g),
    )
    w_leave = (1.0 - s_f) * np.asarray(sp7.w_leave, dtype=np.float64) + s_f * np.asarray(
        sp8.w_leave, dtype=np.float64
    )
    w_return = (1.0 - s_f) * np.asarray(sp7.w_return, dtype=np.float64) + s_f * np.asarray(
        sp8.w_return, dtype=np.float64
    )
    K_eject = (1.0 - s_f) * np.asarray(sp7.K_eject, dtype=np.float64) + s_f * np.asarray(
        sp8.K_eject, dtype=np.float64
    )
    D_g = _lerp(sp7.D_g, sp8.D_g, s_g)
    sp = replace(
        sp8,
        K_eject=K_eject,
        w_leave=w_leave,
        w_return=w_return,
        geo_fb=geo_fb,
        geo_ex=geo_ex,
        D_g=D_g,
        k=sp8.k,
    )
    hydro = q8.hydro
    if hydro is not None:
        w_u7 = np.maximum(
            0.0,
            1.0
            - np.asarray(sp7.w_leave, dtype=np.float64)
            - np.asarray(sp7.w_return, dtype=np.float64),
        )
        hydro = replace(
            hydro,
            u_g=float(geo_fb.u_g),
            rho_g=float(geo_fb.rho_g),
            u_b=float(geo_fb.u_b),
            K_bd=float(geo_ex.K_bd),
            D_g=float(D_g),
            w_leave=w_leave,
            w_return=w_return,
            w_undecided=(1.0 - s_f) * w_u7 + s_f * np.asarray(
                hydro.w_undecided, dtype=np.float64
            ),
        )
    return replace(q8, energy=replace(q8.energy, species=sp), hydro=hydro)


def residual_e7_plus_scaled_e8(
    x: npt.NDArray[np.float64],
    nq7: NamedStageQ,
    nq8: NamedStageQ,
    *,
    hydro_scale: float,
    fate_scale: float | None = None,
) -> npt.NDArray[np.float64]:
    """F(x; 分字段插值 q)。不 clip。"""
    z = np.asarray(x, dtype=np.float64)
    q7 = outer_P_named(z, nq7)
    q8 = outer_P_named(z, nq8)
    q = blend_hydro_named_q(
        q7, q8, hydro_scale=float(hydro_scale), fate_scale=fate_scale
    )
    return coupled_residual_named(z, q, rate_scale=1.0)


def jacobian_named_e8(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_r10: NamedStageQ,
    nq_r11: NamedStageQ,
    nq7: NamedStageQ,
    nq8: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """非 T 列在冻结 ``q=P(x)`` 下对 ``F(x;q)`` 差商；T 列对 ``H(x)`` 差商。

    链底 ``freeze_r12_only_q`` 的 ``w_leave`` 不得进入 N/M 块。
    ``nq6``/``nq_r10``/``nq_r11``/``nq7`` 仅保持调用签名。
    """
    _ = (nq6, nq_r10, nq_r11, nq7)
    z = np.asarray(x, dtype=np.float64)
    nqP = outer_P_named(z, nq8)
    scale = np.maximum(np.abs(z), 0.05)
    scale[-1] = max(float(scale[-1]), 50.0)

    def f_frozen(w: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return coupled_residual_named(w, nqP, rate_scale=1.0)

    J = fd_jacobian(f_frozen, z, scale=scale)
    f0 = residual_named_e8(z, nq8)
    xp = z.copy()
    xp[-1] = float(z[-1]) + _T_FD
    J[:, -1] = (residual_named_e8(xp, nq8) - f0) / _T_FD
    return J


def jacobian_e7_plus_scaled_e8(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
    nq_r10: NamedStageQ,
    nq_r11: NamedStageQ,
    nq7: NamedStageQ,
    nq8: NamedStageQ,
    *,
    hydro_scale: float,
    fate_scale: float | None = None,
) -> npt.NDArray[np.float64]:
    """非 T 列在冻结 ``q(s)`` 下差商；T 列对 ``F(x;q(s(x)))`` 差商。"""
    _ = (nq6, nq_r10, nq_r11)
    z = np.asarray(x, dtype=np.float64)
    s = float(hydro_scale)
    q7 = outer_P_named(z, nq7)
    q8 = outer_P_named(z, nq8)
    q = blend_hydro_named_q(q7, q8, hydro_scale=s, fate_scale=fate_scale)
    scale = np.maximum(np.abs(z), 0.05)
    scale[-1] = max(float(scale[-1]), 50.0)

    def f_frozen(w: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return coupled_residual_named(w, q, rate_scale=1.0)

    J = fd_jacobian(f_frozen, z, scale=scale)
    xp = z.copy()
    xp[-1] = float(z[-1]) + _T_FD
    J[:, -1] = (
        residual_e7_plus_scaled_e8(
            xp, nq7, nq8, hydro_scale=s, fate_scale=fate_scale
        )
        - residual_e7_plus_scaled_e8(
            z, nq7, nq8, hydro_scale=s, fate_scale=fate_scale
        )
    ) / _T_FD
    return J


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return named_in_domain(z)


def _energy_T_step(x: npt.NDArray[np.float64], q: NamedStageQ) -> npt.NDArray[np.float64]:
    """能量行对 T 的一维牛顿。只改 ``x[-1]``，不改 F。"""
    z = np.asarray(x, dtype=np.float64).copy()
    f0 = coupled_residual_named(z, q, rate_scale=1.0)
    zp = z.copy()
    zp[-1] = float(z[-1]) + _T_FD
    dEdT = coupled_residual_named(zp, q, rate_scale=1.0)[-1] - f0[-1]
    if abs(float(dEdT)) > 1.0e-12:
        z[-1] = float(z[-1]) - float(f0[-1]) / float(dEdT)
    return z


def apply_named_e8_exchange_startwert(
    x: npt.NDArray[np.float64],
    q: NamedStageQ,
    *,
    n_inner: int = 8,
) -> npt.NDArray[np.float64]:
    """相间交换 Startwert 试探。延拓不调用。

    把 ``N_ex=K_bd V_b(C_b-C_d)`` 配回旧通量的有限转移会改浓度，
    均相源的增量大于被消掉的交换残差。这里留下的一维转移同样只压 ``F_d``。

    Source: Hamel (1999) Eq. 2.2, 3.50
    """
    z = np.asarray(x, dtype=np.float64).copy()
    for _ in range(int(n_inner)):
        f = coupled_residual_named(z, q, rate_scale=1.0)
        fd = np.asarray(f[0:7], dtype=np.float64)
        if float(np.linalg.norm(fd)) < 1.0e-10:
            break
        jt = np.zeros((7, 7), dtype=np.float64)
        for j in range(7):
            h = max(1.0e-8, 1.0e-6 * abs(float(z[j])))
            if float(z[j]) - h < 0.0:
                h = 0.5 * float(z[j]) if float(z[j]) > 0.0 else 0.0
            if h <= 0.0:
                continue
            zp = z.copy()
            zp[j] = float(zp[j]) - h
            zp[7 + j] = float(zp[7 + j]) + h
            fp = coupled_residual_named(zp, q, rate_scale=1.0)
            jt[:, j] = (np.asarray(fp[0:7], dtype=np.float64) - fd) / h
        try:
            t = np.linalg.solve(jt, -fd)
        except np.linalg.LinAlgError:
            break
        for j in range(7):
            if float(z[j]) - float(t[j]) < 0.0:
                t[j] = 0.9 * float(z[j])
            if float(z[7 + j]) + float(t[j]) < 0.0:
                t[j] = -0.9 * float(z[7 + j])
        z[0:7] = z[0:7] - t
        z[7:14] = z[7:14] + t
        z[0:14] = np.maximum(z[0:14], 0.0)
    return z


_E8_STEP_EVIDENCE: dict[str, Any] = {}
_DIAG_SKIP: frozenset[str] = frozenset(
    {"F", "J_modified", "J_eq", "J_hat", "x_scale", "F_scale", "dr", "dc"}
)


def apply_e8_certified_rollback(
    last: Any,
    x0: npt.NDArray[np.float64],
    F0: npt.NDArray[np.float64],
    cert: dict[str, Any],
    cert0: dict[str, Any],
) -> tuple[Any, dict[str, Any], str]:
    """未认证且旧根已认证时退回旧根。不改 F 定义。"""
    if (not bool(cert["endpoint_certified"])) and bool(cert0["endpoint_certified"]):
        last.x = np.asarray(x0, dtype=np.float64)
        last.F = np.asarray(F0, dtype=np.float64)
        return last, cert0, "e7_certified_root"
    return last, cert, "newton"


def _e8_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    dx, rank, cond, diag = named_linear_step(
        J, F, x,
        extra_diag_species=("H2", "O2", "CO", "CH4"),
        phase_offsets=(0, 7, 16),
    )
    j_mod = diag.get("J_modified")
    p, changed = bound_feasible_direction(x, dx)
    feas = diagnose_linear_step(
        J, F, p,
        used_trace_diag=bool(diag["used_trace_diag"]),
        J_modified=j_mod,
    )
    raw_only = {k: diag[k] for k in diag if k not in _DIAG_SKIP}
    _E8_STEP_EVIDENCE.clear()
    _E8_STEP_EVIDENCE.update(
        {
            "raw_J": {
                "raw_rank": int(diag["raw_rank"]),
                "equilibrated_rank_original": int(diag["equilibrated_rank_original"]),
                "linear_residual_raw_l2_unclipped": float(diag["linear_residual_raw_l2"]),
                "directional_phi_dot_unclipped": float(diag["directional_phi_dot"]),
            },
            "modified_matrix": {
                "used_trace_diag": bool(diag["used_trace_diag"]),
                "trace_diag": float(diag["trace_diag"]),
                "equilibrated_rank_used": int(diag["equilibrated_rank_used"]),
                "linear_residual_modified_l2_unclipped": (
                    None if j_mod is None
                    else float(diag.get("linear_residual_modified_l2", np.nan))
                ),
                "linear_residual_modified_l2_feasible": feas.get(
                    "linear_residual_modified_l2"
                ),
            },
            "feasible_direction": {
                "direction_clipped": bool(changed),
                "linear_residual_raw_l2": float(feas["linear_residual_raw_l2"]),
                "directional_phi_dot": float(feas["directional_phi_dot"]),
                "p": np.asarray(p, dtype=np.float64).tolist(),
            },
            "rollback": {"returned_state_source": "pending"},
            "dx_l2": float(np.linalg.norm(dx)),
            "p_l2": float(np.linalg.norm(p)),
            "diag_scalar": raw_only,
        }
    )
    return p, rank, cond


def _fd_jacobian_trace(residual_fn, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """差商。``|x_i|<1e-8`` 时步长与该微量同阶。

    全局 ``fd_jacobian`` 的绝对步长约 ``5e-10``。气泡 ``N_O2~1e-14`` 上
    这条割线比 ``R12`` 对 ``C_O2`` 的真导数小约两个数量级，``Jp=-F`` 仍成立，
    但沿该方向的 ``F`` 上升。不改 ``jacobian_named_e8`` 的 P0 装配。
    """
    z = np.asarray(x, dtype=np.float64)
    f0 = np.asarray(residual_fn(z), dtype=np.float64)
    jac = np.zeros((z.size, z.size), dtype=np.float64)
    for i in range(int(z.size)):
        xi = float(z[i])
        if i == z.size - 1:
            eps = max(1.0e-6, 1.0e-8 * abs(xi))  # [K]
        elif abs(xi) < _TRACE_ABS:
            eps = max(_TRACE_EPS_FLOOR, 1.0e-3 * max(abs(xi), 1.0e-30))
        else:
            eps = max(1.0e-10, 1.0e-8 * max(abs(xi), 0.05))
        xp = z.copy()
        xp[i] = xi + eps
        jac[:, i] = (np.asarray(residual_fn(xp), dtype=np.float64) - f0) / eps
    return jac


def _e8_active_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], int, float]:
    """非正未知量从线性步拿掉。正的微量氧保留：差商已按它缩放。不改 F。"""
    jm = np.asarray(J, dtype=np.float64).copy()
    fm = np.asarray(F, dtype=np.float64).copy()
    z = np.asarray(x, dtype=np.float64)
    frozen = np.zeros(z.size, dtype=bool)
    frozen[:-1] = z[:-1] <= 0.0
    # 残差还在的零分量要打开。自由板 O2 贴零时 R12 仍要 ~1e-21 mol/s 才能配平入口。
    for i in range(int(z.size) - 1):
        if frozen[i] and abs(float(F[i])) > 1.0e-12:
            frozen[i] = False
    for i in np.where(frozen)[0]:
        jm[i, :] = 0.0
        jm[:, i] = 0.0
        jm[i, i] = 1.0
        fm[i] = 0.0
    dx = np.linalg.solve(jm, -fm)
    dx[frozen] = 0.0
    p = dx
    for i in range(int(p.size) - 1):
        if z[i] + p[i] < 0.0:
            p[i] = -0.5 * float(z[i]) if float(z[i]) > 0.0 else 0.0
    return p, int(z.size), 0.0


_GAS_ENERGY_FREE: tuple[int, ...] = tuple(
    i for i in list(range(14)) + list(range(16, 37)) if i not in (21, 28, 35)
)


def _newton_gas_energy_reserved(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
    *,
    n_iter: int = 8,
) -> npt.NDArray[np.float64]:
    """床层气体与下游气体一起走，温度列不进这些方程。

    全步牛顿在水力尺度约 0.98 处用 ``ΔT`` 去消交换残差，真实 ``F`` 上升。
    这里温度只配能量行。下游氧（自由板、两级旋风）保持不动：那里的 R12
    对微量氧是硬的。差商只建一次。

    这一步不改焦油。焦油行因此按水力尺度累积，大约每 ``10^{-5}`` 留下
    ``6×10^{-7} mol/s``。把它并进同一线性系统会经 R10 把气泡氧列的条件数
    抬到约 ``10^{16}``，方向不再下降。
    """
    z = np.asarray(x, dtype=np.float64).copy()

    def res(w: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return coupled_residual_named(w, q, rate_scale=1.0)

    f0 = res(z)
    jac = np.zeros((z.size, z.size), dtype=np.float64)
    for i in _GAS_ENERGY_FREE + (z.size - 1,):
        xi = float(z[i])
        if i == z.size - 1:
            eps = max(1.0e-6, 1.0e-8 * abs(xi))
        elif abs(xi) < _TRACE_ABS:
            eps = max(_TRACE_EPS_FLOOR, 1.0e-3 * max(abs(xi), 1.0e-30))
        else:
            eps = max(1.0e-10, 1.0e-8 * max(abs(xi), 0.05))
        xp = z.copy()
        xp[i] = xi + eps
        jac[:, i] = (res(xp) - f0) / eps
    frozen = [i for i in range(z.size - 1) if i not in _GAS_ENERGY_FREE]
    for _k in range(int(n_iter)):
        f = res(z)
        phi = float(np.linalg.norm(f))
        if phi < 5.0e-7:
            return z
        jm = jac.copy()
        jm[:-1, -1] = 0.0
        fm = f.copy()
        for i in frozen:
            jm[i, :] = 0.0
            jm[:, i] = 0.0
            jm[i, i] = 1.0
            fm[i] = 0.0
        p = np.linalg.solve(jm, -fm)
        for i in frozen:
            p[i] = 0.0
        for i in range(int(p.size) - 1):
            if z[i] + p[i] < 0.0:
                p[i] = -0.5 * float(z[i]) if float(z[i]) > 0.0 else 0.0
        moved = False
        alpha = 1.0
        while alpha >= 1.0 / 64.0:
            trial = z + alpha * p
            if float(np.min(trial[:-1])) < 0.0 or float(trial[-1]) <= 0.0:
                alpha *= 0.5
                continue
            if float(np.linalg.norm(res(trial))) < phi:
                z = trial
                moved = True
                break
            alpha *= 0.5
        if not moved:
            return z
    return z


def residual_block_scale(n: int) -> npt.NDArray[np.float64]:
    """与内层最小二乘相同的分块尺度。气体/焦油 0.01 mol/s，固体 0.001 kg/s，能量 1 W。"""
    scale = np.full(int(n), 0.01, dtype=np.float64)
    for i in (14, 15, 37, 38, 39, 46, 47):
        if i < scale.size:
            scale[i] = 0.001
    scale[-1] = 1.0
    if scale.size >= 46:
        scale[-2] = 1.0
    return scale


def residual_block_merit(f: npt.NDArray[np.float64]) -> float:
    """无量纲残差。只用于比较候选步，不代替原单位分块门。"""
    z = np.asarray(f, dtype=np.float64)
    return float(np.linalg.norm(z / residual_block_scale(int(z.size))))


class _InnerBudget(Exception):
    """内层 J/F 预算用尽。不是物理根，也不是差分探针。"""


class _StepBudget(Exception):
    """线性子问题碰到校正器的截止时刻。不是物理停滞。"""


def _stop_at_deadline(deadline_at: float | None):
    """SLSQP 与 trust-constr 共用。回调抛出后不再继续纯代数迭代。"""

    def callback(xk: npt.NDArray[np.float64], *_args: Any) -> None:
        if deadline_at is not None and time.monotonic() >= float(deadline_at):
            raise _StepBudget()

    return callback


def _scaled_tar_row(
    z: npt.NDArray[np.float64],
    qq: NamedStageQ,
    scale: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    raw = apply_tar_reaction_transport_rows(
        z, qq, coupled_residual_named(z, qq, rate_scale=1.0),
    )
    return raw / scale


# 首个受控试验值，不是物理常数。约束的是组合方向上的沉降余量变化。
SETTLING_MARGIN_ETA = 0.2
_EJECT_SUBSTANTIAL_KG_S = 1.0e-8


def freeboard_settling_margin_m(
    z: npt.NDArray[np.float64],
    b: Any,
) -> npt.NDArray[np.float64]:
    """均匀场沉降余量 ``m = v_t - u_g`` [m/s]。``m>0`` 表示颗粒终速向下。

    ``v_t`` 随该状态的自由空间密度和组成变化，不只随表观气速变化。
    """
    state = np.asarray(z, dtype=np.float64)
    u_g, rho_g, _n, _mw = gas_velocity_and_density(
        state[16:23], np.zeros(7, dtype=np.float64),
        T=float(state[-1]), P=float(b.P), area=float(b.area),
    )
    spec = resolve_freeboard_model(FREEBOARD_MODEL_HAMEL_EXTENDED)
    z_drag = drag_z_for_spec(spec, phi_s=float(b.phi_s))
    margin = np.zeros(int(b.d_p.size), dtype=np.float64)
    for i, d_p in enumerate(b.d_p):
        u_eq = drag_equilibrium_speed_m_s(
            u_g=u_g, u_p=0.0, rho_g=rho_g, rho_p=float(b.rho_p),
            d_p=float(d_p), mu_g=float(b.mu_g), z=z_drag,
        )
        margin[i] = -float(u_eq)
    return margin


def _settling_margin_limits(
    z: npt.NDArray[np.float64],
    q: NamedStageQ,
    unpack,
    n_var: int,
    eta: float | None,
) -> list[tuple[npt.NDArray[np.float64], float]]:
    """``|g·v| ≤ eta |m|``。入口 m 取当前目标参数下的状态，不写进残差。"""
    if eta is None or q.reference is None or int(z.size) < 23:
        return []
    b = q.reference
    m0 = freeboard_settling_margin_m(z, b)
    mass = np.asarray(z[14:16], dtype=np.float64)
    eject = np.asarray(q.energy.species.K_eject, dtype=np.float64) * mass
    limits: list[tuple[npt.NDArray[np.float64], float]] = []
    h = 1.0e-4
    for i, m_i in enumerate(m0):
        if i >= int(eject.size) or float(eject[i]) < _EJECT_SUBSTANTIAL_KG_S:
            continue
        if abs(float(m_i)) < 1.0e-6:
            continue
        grad = np.zeros(int(n_var), dtype=np.float64)
        for j in range(int(n_var)):
            vp = np.zeros(int(n_var), dtype=np.float64)
            vm = vp.copy()
            vp[j] = h
            vm[j] = -h
            mp = freeboard_settling_margin_m(unpack(vp), b)
            mm = freeboard_settling_margin_m(unpack(vm), b)
            grad[j] = (float(mp[i]) - float(mm[i])) / (2.0 * h)
        limits.append((grad, float(eta) * abs(float(m_i))))
    return limits


def _trajectory_log_state_index(index: int, n_state: int) -> bool:
    """这一列会改变抛射速度或自由空间气速，从而改变飞行时间。

    床层气量与床温进入床面气泡速度。自由空间气量和自由空间温度进入表观气速。
    旋风气量和焦油比例不在这条轨迹里。
    """
    i = int(index)
    n = int(n_state)
    if i < 14 or 16 <= i < 23:
        return True
    n_t = 2 if n >= 46 else 1
    return i >= n - n_t


def _add_exchange_chain(
    jac: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
    q: NamedStageQ,
    ids: npt.NDArray[np.int64],
    n_log: int,
    log_ids: npt.NDArray[np.int64],
    scale: npt.NDArray[np.float64],
    pyro_cache: dict[float, NamedStageQ],
    skip_columns: set[int] | None = None,
) -> tuple[float, float]:
    """把 ``K_eject`` 和 ``K_bd`` 对对数未知量的链式项加进内层雅可比。

    内层差商冻结这两个系数。这里只对床层气量和床温差分床面公式，不重积分轨迹。
    抛射系数进入所有使用它的质量和能量行。返回残差次数和床面公式次数。
    ``skip_columns`` 里的列已经用完整 ``P`` 差商替换，不再叠加本项。
    """
    if q.reference is None or int(q.stage) < 8 or n_log <= 0:
        return 0.0, 0.0
    skipped = set() if skip_columns is None else set(skip_columns)
    b = q.reference
    qq = _pyrolysis_q_cached(q, z, pyro_cache)
    K0, kbd0, *_rest = bed_face_exchange(z, b, T_bed_K=_bed_temperature(z))
    slopes_k = np.zeros((int(K0.size), n_log), dtype=np.float64)
    slopes_b = np.zeros(n_log, dtype=np.float64)
    bed_ids = set(range(14))
    t_bed = int(z.size - 2) if z.size >= 46 else int(z.size - 1)
    n_face = 1.0
    for j in range(n_log):
        index = int(log_ids[j])
        if index not in bed_ids and index != t_bed:
            continue
        n_face += 2.0
        h = 1.0e-7 if index == t_bed else 1.0e-4
        zp = z.copy()
        zm = z.copy()
        zp[index] *= float(np.exp(h))
        zm[index] *= float(np.exp(-h))
        kp, bp, *_ = bed_face_exchange(zp, b, T_bed_K=_bed_temperature(zp))
        km, bm, *_ = bed_face_exchange(zm, b, T_bed_K=_bed_temperature(zm))
        slopes_k[:, j] = (kp - km) / (2.0 * h)
        slopes_b[j] = (bp - bm) / (2.0 * h)
    species = qq.energy.species
    columns = []
    for i in range(int(K0.size)):
        step = max(1.0e-8, 1.0e-4 * abs(float(K0[i])))
        plus = np.asarray(K0, dtype=np.float64).copy()
        minus = plus.copy()
        plus[i] += step
        minus[i] -= step
        columns.append((
            _row_with_exchange(z, qq, species, plus, float(kbd0), scale)
            - _row_with_exchange(z, qq, species, minus, float(kbd0), scale)
        ) / (2.0 * step))
    step_b = max(1.0e-8, 1.0e-4 * abs(float(kbd0)))
    row_b = (
        _row_with_exchange(z, qq, species, K0, float(kbd0) + step_b, scale)
        - _row_with_exchange(z, qq, species, K0, float(kbd0) - step_b, scale)
    ) / (2.0 * step_b)
    for j, index in enumerate(ids):
        if j >= n_log or int(j) in skipped:
            continue
        column = row_b * float(slopes_b[j])
        for i, drow in enumerate(columns):
            column = column + drow * float(slopes_k[i, j])
        jac[:, j] = jac[:, j] + column
    return 2.0 * float(K0.size) + 2.0, n_face


def _bed_temperature(z: npt.NDArray[np.float64]) -> float:
    if int(z.size) >= 46:
        return float(z[-2])
    return float(z[-1])


def _row_with_exchange(
    z: npt.NDArray[np.float64],
    qq: NamedStageQ,
    species: Any,
    K_eject: npt.NDArray[np.float64],
    K_bd: float,
    scale: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    geo = replace(species.geo_ex, K_bd=float(K_bd))
    patched = replace(species, K_eject=np.asarray(K_eject, dtype=np.float64), geo_ex=geo)
    qk = replace(qq, energy=replace(qq.energy, species=patched))
    return _scaled_tar_row(z, qk, scale)


def _gate_vector(n: int) -> npt.NDArray[np.float64]:
    """原始残差的门尺度。物料行 1e-6，能量行 1 W。与 ``progress_phi`` 同一度量。"""
    gate = np.full(int(n), 1.0e-6, dtype=np.float64)
    gate[-1] = 1.0
    if int(n) >= 46:
        gate[-2] = 1.0
    return gate


def _material_blocks(n: int) -> tuple[npt.NDArray[np.int64], ...]:
    gas = list(range(min(14, n)))
    if n >= 37:
        gas.extend(range(16, 37))
    elif n >= 23:
        gas.extend(range(16, 23))
    solid = [i for i in (14, 15, 37, 38, 39, 46, 47) if i < n - 1]
    tar = list(range(40, 46)) if n >= 46 else []
    return (
        np.asarray(gas, dtype=np.int64),
        np.asarray(solid, dtype=np.int64),
        np.asarray(tar, dtype=np.int64),
    )


_ND_FEASIBLE = 1.0e-8  # [-] 无量纲约束容差。物料和能量各自用门尺度，不共用 1e-8 的有量纲平方。


def material_protection_margin(block_norm: float, tol: float = 1.0e-6) -> float:
    """``1 - ||F_block||^2 / tol^2``。块范数等于门时为 0，十倍门时为 -99。"""
    ratio = float(block_norm) / float(tol)
    return 1.0 - ratio * ratio


def _gate_progress_step(
    matrix: npt.NDArray[np.float64],
    entrance_row: npt.NDArray[np.float64],
    scale: npt.NDArray[np.float64],
    lo: npt.NDArray[np.float64],
    hi: npt.NDArray[np.float64],
    quality: dict[str, Any] | None = None,
    *,
    deadline_at: float | None = None,
    energy_guard_W: float = 0.0,
    linear_abs_limits: list[tuple[npt.NDArray[np.float64], float]] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], int]:
    """在原始门尺度上解有界线性步。变量和约束都无量纲化，约束带解析导数。

    返回状态：1 表示子问题已解好；2 表示可行下降但求解器未成功，只是非精确步；
    0 表示退回零步。状态 2 不能称为已解好，也不能称为没有下降方向。
    """
    from scipy.optimize import lsq_linear, minimize

    n = int(entrance_row.size)
    gate = _gate_vector(n)
    jacobian = np.column_stack([
        original_residual_from_tar_rows(matrix[:, j] * scale) for j in range(matrix.shape[1])
    ])
    base = original_residual_from_tar_rows(entrance_row * scale)
    design = jacobian / gate[:, None]
    target = base / gate
    phi0 = float(np.linalg.norm(target))
    nv = int(design.shape[1])
    report = quality if quality is not None else {}
    if nv == 0:
        report.update({"solved": True, "status": 0, "message": "no variables", "max_nd_violation": 0.0})
        return np.zeros(0, dtype=np.float64), base, 0

    def cost(v: npt.NDArray[np.float64]) -> float:
        gap = design @ v + target
        return 0.5 * float(gap @ gap)

    def grad(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return design.T @ (design @ v + target)

    v_margins: list[Any] = []
    n_t = 2 if n >= 46 else 1
    def _energy_limit(level: float) -> float:
        if abs(float(level)) <= 1.0:
            return max(1.0e-3, 1.0 - max(float(energy_guard_W), 0.0))
        return abs(float(level))

    for k in range(n - n_t, n):
        level = float(base[k])
        limit = _energy_limit(level)
        row = np.asarray(jacobian[k], dtype=np.float64).copy()
        v_margins.append(lambda v, row=row, level=level, limit=limit: (
            limit - (level + float(row @ v))
        ) / limit)
        v_margins.append(lambda v, row=row, level=level, limit=limit: (
            limit + (level + float(row @ v))
        ) / limit)
    for block in _material_blocks(n):
        if int(block.size) == 0 or float(np.linalg.norm(base[block])) >= 1.0e-6:
            continue
        block_base = np.asarray(base[block], dtype=np.float64).copy()
        block_j = np.asarray(jacobian[block], dtype=np.float64).copy()
        v_margins.append(lambda v, block_base=block_base, block_j=block_j: (
            material_protection_margin(float(np.linalg.norm(block_base + block_j @ v)))
        ))
    for margin_grad, limit in list(linear_abs_limits or []):
        row = np.asarray(margin_grad, dtype=np.float64).copy()
        cap = max(float(limit), 0.0)
        v_margins.append(lambda v, row=row, cap=cap: cap - float(row @ v))
        v_margins.append(lambda v, row=row, cap=cap: cap + float(row @ v))
    bounds = (
        np.asarray(lo, dtype=np.float64).copy(),
        np.asarray(hi, dtype=np.float64).copy(),
    )
    order = bounds[0] > bounds[1]
    bounds[0][order] = bounds[1][order]
    span = np.maximum(np.maximum(np.abs(bounds[0]), np.abs(bounds[1])), 1.0e-12)
    free = lsq_linear(design, -target, bounds=bounds, method="bvls", tol=1.0e-12, max_iter=200)
    start = np.asarray(free.x, dtype=np.float64)

    def feasible(v: npt.NDArray[np.float64]) -> bool:
        return all(float(margin(v)) >= -_ND_FEASIBLE for margin in v_margins)

    alpha = 1.0
    while alpha > 1.0e-8 and not feasible(alpha * start):
        alpha *= 0.5
    start = alpha * start
    if not feasible(start):
        start = np.zeros(nv, dtype=np.float64)
    obj_norm = max(cost(start), 1.0)
    y_lo = bounds[0] / span
    y_hi = bounds[1] / span
    constraints = []
    constraint_hessians: list[Any] = []
    for k in range(n - n_t, n):
        level = float(base[k])
        limit = _energy_limit(level)
        row = np.asarray(jacobian[k], dtype=np.float64).copy()
        constraints.append({
            "type": "ineq",
            "fun": lambda y, row=row, level=level, limit=limit: (
                limit - (level + float(row @ (span * y)))
            ) / limit,
            "jac": lambda y, row=row, limit=limit: -(row * span) / limit,
        })
        constraint_hessians.append(None)
        constraints.append({
            "type": "ineq",
            "fun": lambda y, row=row, level=level, limit=limit: (
                limit + (level + float(row @ (span * y)))
            ) / limit,
            "jac": lambda y, row=row, limit=limit: (row * span) / limit,
        })
        constraint_hessians.append(None)
    for block in _material_blocks(n):
        if int(block.size) == 0 or float(np.linalg.norm(base[block])) >= 1.0e-6:
            continue
        block_base = np.asarray(base[block], dtype=np.float64).copy()
        block_j = np.asarray(jacobian[block], dtype=np.float64).copy()

        def material_fun(y, block_base=block_base, block_j=block_j):
            return material_protection_margin(float(np.linalg.norm(block_base + block_j @ (span * y))))

        def material_jac(y, block_base=block_base, block_j=block_j):
            phys = block_base + block_j @ (span * y)
            # d(1-||phys||^2/tol^2)/dy = -2/tol^2 * phys^T * block_j * diag(span)
            return (-2.0 / 1.0e-12) * (block_j.T @ phys) * span

        constraints.append({"type": "ineq", "fun": material_fun, "jac": material_jac})

        def material_hess(y, block_j=block_j):
            gram = block_j.T @ block_j
            return (-2.0 / 1.0e-12) * (span[:, None] * gram * span[None, :])

        constraint_hessians.append(material_hess)
    for margin_grad, limit in list(linear_abs_limits or []):
        row = np.asarray(margin_grad, dtype=np.float64).copy()
        cap = max(float(limit), 0.0)
        constraints.append({
            "type": "ineq",
            "fun": lambda y, row=row, cap=cap: cap - float(row @ (span * y)),
            "jac": lambda y, row=row: -(row * span),
        })
        constraint_hessians.append(None)
        constraints.append({
            "type": "ineq",
            "fun": lambda y, row=row, cap=cap: cap + float(row @ (span * y)),
            "jac": lambda y, row=row: row * span,
        })
        constraint_hessians.append(None)

    def cost_y(y: npt.NDArray[np.float64]) -> float:
        return cost(span * y) / obj_norm

    def grad_y(y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return grad(span * y) * span / obj_norm

    timed_out = deadline_at is not None and time.monotonic() >= float(deadline_at)
    fit = None
    if not timed_out:
        try:
            fit = minimize(
                cost_y,
                start / span,
                jac=grad_y,
                bounds=list(zip(y_lo, y_hi)),
                constraints=constraints,
                method="SLSQP",
                tol=1.0e-10,
                callback=_stop_at_deadline(deadline_at),
                options={"maxiter": 200, "ftol": 1.0e-12},
            )
        except _StepBudget:
            timed_out = True
    if fit is not None and not bool(fit.success) and constraints and not timed_out:
        from scipy.optimize import NonlinearConstraint

        def stacked(y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return np.array([float(row["fun"](y)) for row in constraints], dtype=np.float64)

        def stacked_jac(y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return np.vstack([np.atleast_1d(row["jac"](y)) for row in constraints])

        def objective_hess(y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            gram = design.T @ design
            return (span[:, None] * gram * span[None, :]) / obj_norm

        def constraint_hess(y: npt.NDArray[np.float64], weights: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            total = np.zeros((int(y.size), int(y.size)), dtype=np.float64)
            for weight, hess in zip(np.atleast_1d(weights), constraint_hessians):
                if hess is not None:
                    total = total + float(weight) * hess(y)
            return total

        cone = NonlinearConstraint(stacked, 0.0, np.inf, jac=stacked_jac, hess=constraint_hess)
        try:
            second = minimize(
                cost_y,
                np.asarray(fit.x, dtype=np.float64),
                jac=grad_y,
                hess=objective_hess,
                bounds=list(zip(y_lo, y_hi)),
                constraints=cone,
                method="trust-constr",
                callback=_stop_at_deadline(deadline_at),
                options={"maxiter": 200, "verbose": 0, "gtol": 1.0e-8},
            )
        except _StepBudget:
            timed_out = True
            second = None
        if second is not None:
            cand = span * np.asarray(second.x, dtype=np.float64)
            current = span * np.asarray(fit.x, dtype=np.float64)
            report["fallback_success"] = bool(second.success)
            report["fallback_status"] = int(second.status)
            report["fallback_message"] = str(second.message)
            report["fallback_feasible"] = bool(feasible(cand))
            if feasible(cand) and (bool(second.success) or cost(cand) < cost(current)):
                fit = second
    if fit is None:
        fit = type("_Missed", (), {
            "success": False,
            "status": -1,
            "message": "deadline",
            "nit": 0,
            "x": start / span,
        })()
    attempted = np.minimum(np.maximum(span * np.asarray(fit.x, dtype=np.float64), bounds[0]), bounds[1])
    use = attempted if feasible(attempted) else start
    phi = float(np.linalg.norm(design @ use + target))
    descent = phi < phi0 - 1.0e-9
    if phi > phi0 + 1.0e-9 or not descent:
        use = np.zeros(nv, dtype=np.float64)
        phi = phi0
        descent = False
    solved = (
        bool(fit.success)
        and feasible(attempted)
        and float(np.linalg.norm(design @ attempted + target)) <= phi0 + 1.0e-9
        and max((max(0.0, -float(margin(attempted))) for margin in v_margins), default=0.0) <= _ND_FEASIBLE
    )
    if solved:
        use = attempted
        phi = float(np.linalg.norm(design @ use + target))
        status = 1
    elif descent:
        status = 2
    else:
        status = 0
        use = np.zeros(nv, dtype=np.float64)
        phi = phi0
    violation = max((max(0.0, -float(margin(use))) for margin in v_margins), default=0.0)
    gy = grad(use) * span
    y = use / span
    projected = y - np.minimum(np.maximum(y - gy, y_lo), y_hi)
    active = int(np.sum(
        (use <= bounds[0] + 1.0e-8 * span) | (use >= bounds[1] - 1.0e-8 * span)
    ))
    report.update({
        "success": bool(fit.success),
        "status": int(fit.status),
        "message": str(fit.message),
        "nit": int(fit.nit),
        "objective": float(phi),
        "objective_at_zero": float(phi0),
        "max_nd_violation": float(violation),
        "active_bound_count": active,
        "projected_gradient_norm": float(np.linalg.norm(projected)),
        "solved": bool(solved),
        "inexact_descent": bool(status == 2),
        "deadline_hit": bool(timed_out),
    })
    return use, base + jacobian @ use, status


_CYCLONE_MASS_INDEX: tuple[int, int] = (37, 38)


def cyclone_equilibrium_mass_kg(
    inflow_kg_s: float,
    k_auf: float,
    k_ab: float,
) -> npt.NDArray[np.float64]:
    """两级旋风对入流的平衡库存 [kg]。排出系数冻结。

    ``0 = L + K_ab M_2 - (K_auf + K_ab) M_1``
    ``0 = K_auf M_1 - (K_auf + K_ab) M_2``

    Source: Hamel (1999) Eq. 2.6
    """
    inflow = float(inflow_kg_s)
    k_up = float(k_auf)
    k_down = float(k_ab)
    denom = k_up * k_up + k_up * k_down + k_down * k_down
    mass = np.zeros(2, dtype=np.float64)
    if inflow > 0.0 and denom > 0.0:
        mass[0] = inflow * (k_up + k_down) / denom
        mass[1] = inflow * k_up / denom
    return mass


def cyclone_reactivation_layout(
    mass_kg: npt.NDArray[np.float64],
    equilibrium_kg: npt.NDArray[np.float64],
    *,
    k_auf: float,
    k_ab: float,
    inflow_kg_s: float,
    log_span: float,
) -> dict[str, npt.NDArray[np.float64]]:
    """入流尺度超出当前对数步长时，把旋风库存改成有尺度的线性坐标。

    ``M = max(base + scale v, 0)``，``v = 0`` 是当前库存。尺度取平衡库存和当前
    库存中较大的一个。入流和现有排出都低于固体门时保持对数坐标。返回的下标、
    基准、尺度和边界一一对应。
    """
    mass = np.asarray(mass_kg, dtype=np.float64)
    equilibrium = np.asarray(equilibrium_kg, dtype=np.float64)
    k_sum = float(k_auf) + float(k_ab)
    span = max(float(log_span), 0.0)
    picked: list[int] = []
    base: list[float] = []
    scale: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    for slot, index in enumerate(_CYCLONE_MASS_INDEX):
        if slot >= int(mass.size) or slot >= int(equilibrium.size):
            break
        current = max(float(mass[slot]), 0.0)
        target = max(float(equilibrium[slot]), 0.0)
        live = (
            float(inflow_kg_s) >= float(ENDPOINT_L2)
            or k_sum * current >= float(ENDPOINT_L2)
        )
        if not live or max(current, target) <= 0.0:
            continue
        separated = current <= 0.0 or target <= 0.0 or abs(float(np.log(target / current))) > span
        if not separated:
            continue
        unit = max(current, target)
        step = (target - current) / unit
        picked.append(int(index))
        base.append(current)
        scale.append(unit)
        lower.append(-current / unit)
        upper.append(max(step, 0.0) + span)
    return {
        "index": np.asarray(picked, dtype=np.int64),
        "base": np.asarray(base, dtype=np.float64),
        "scale": np.asarray(scale, dtype=np.float64),
        "lo": np.asarray(lower, dtype=np.float64),
        "hi": np.asarray(upper, dtype=np.float64),
    }


def _cyclone_reactivation_for_state(
    z: npt.NDArray[np.float64],
    q: NamedStageQ,
    log_span: float,
) -> dict[str, npt.NDArray[np.float64]]:
    """由当前反应后入流和冻结排出系数构造旋风重激活坐标。"""
    empty = {
        "index": np.zeros(0, dtype=np.int64),
        "base": np.zeros(0, dtype=np.float64),
        "scale": np.zeros(0, dtype=np.float64),
        "lo": np.zeros(0, dtype=np.float64),
        "hi": np.zeros(0, dtype=np.float64),
    }
    if int(z.size) <= _CYCLONE_MASS_INDEX[-1]:
        return empty
    species = q.energy.species
    k = getattr(species, "k", None)
    if k is None:
        return empty
    _eject, _returned, leave = char_face_flows_kg_s(z, q)
    inflow = float(np.sum(leave))
    equilibrium = cyclone_equilibrium_mass_kg(inflow, float(k.K_auf), float(k.K_ab))
    return cyclone_reactivation_layout(
        np.asarray(z[_CYCLONE_MASS_INDEX[0]:_CYCLONE_MASS_INDEX[-1] + 1], dtype=np.float64),
        equilibrium,
        k_auf=float(k.K_auf),
        k_ab=float(k.K_ab),
        inflow_kg_s=inflow,
        log_span=float(log_span),
    )


def _e8_log_scaled_correct(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
    *,
    deadline_s: float | None = None,
    max_nfev: int = 800,
    ftol: float = 1.0e-12,
    temperature_factor: float | None = None,
    log_span: float = 40.0,
    ratio_span: float = 1.0,
    max_j: int | None = None,
    stats: dict[str, Any] | None = None,
    residual_offset: npt.NDArray[np.float64] | None = None,
    temperature_radius_scale: float = 1.0,
    couple_trajectory_jacobian: bool = False,
    energy_guard_W: float = 0.0,
    settling_margin_eta: float | None = None,
) -> npt.NDArray[np.float64]:
    """正坐标上的缩放最小二乘。不是 Hamel 规定的牛顿法。

    已为正的未知量走 ``x = x_ref exp(v)``，含低于 ``1e-12`` mol/s 的自由空间氧。
    两级旋风库存在入流尺度超出本步对数跨度时改走 ``M = max(base + scale v, 0)``，
    不再同时保留该变量的对数列。尺度来自冻结排出系数和当前反应后入流。
    标量配平只在自由空间氧仍非正时构造初值，联合步之后不再覆盖。
    差商在对数坐标上取固定步长。
    残差按气体/焦油 0.01 mol/s、固体 0.001 kg/s、能量 1 W 缩放，只用于求解。
    ``couple_trajectory_jacobian`` 为真时，影响抛射速度或自由空间气速的对数列
    改用完整 ``P`` 差商，并且不再叠加床面交换链式项。其余列仍冻结飞行时间。
    """
    z0 = _seed_downstream_o2(q, np.asarray(x, dtype=np.float64).copy())
    f_ids = coupled_residual_named(z0, q, rate_scale=1.0)
    a1, a2 = tar_consumption_fractions() if int(z0.size) >= 46 else (0.0, 0.0)
    pairs = []
    if int(z0.size) >= 46:
        for i1, i2 in ((40, 41), (42, 43), (44, 45)):
            if float(z0[i1]) < 0.0 or float(z0[i2]) < 0.0:
                continue
            total, frac = tar_sum_and_fraction(float(z0[i1]), float(z0[i2]), a1, a2)
            frac = min(max(frac, -a1), a2)
            pairs.append((i1, i2, total, frac))
    paired = {i for i1, i2, _s, _f in pairs for i in (i1, i2)}
    cyclone = _cyclone_reactivation_for_state(z0, q, log_span)
    react_ids = np.asarray(cyclone["index"], dtype=np.int64)
    react_base = np.asarray(cyclone["base"], dtype=np.float64)
    react_scale = np.asarray(cyclone["scale"], dtype=np.float64)
    react_set = {int(i) for i in react_ids}
    log_ids = np.asarray([
        i for i in _e8_positive_log_ids(z0)
        if int(i) not in paired and int(i) not in react_set
    ], dtype=np.int64)
    zero_ids = np.asarray(
        [i for i in _linear_bubble_zero_ids(z0, f_ids) if int(i) not in paired and int(i) not in react_set],
        dtype=np.int64,
    )
    n_log = int(log_ids.size)
    n_pair = len(pairs)
    n_react = int(react_ids.size)
    ids = np.arange(n_log + 2 * n_pair + int(zero_ids.size) + n_react, dtype=np.int64)
    base = z0[log_ids].copy() if n_log else np.zeros(0, dtype=np.float64)
    react_at = n_log + 2 * n_pair + int(zero_ids.size)
    scale = residual_block_scale(int(z0.size))
    hi_lin = max(1.0, float(np.sum(np.maximum(z0, 0.0))))
    t0 = time.monotonic()
    calls = {"F": 0, "J": 0, "P": 0, "P_wall_s": 0.0, "chain": 0.0, "bed_face": 0.0}
    cached = None if stats is None else stats.get("jacobian_cache")
    entrance_row = apply_tar_reaction_transport_rows(
        z0, q, inner_hydro_frozen_residual(z0, q),
    ) / scale
    state = {
        "v": np.zeros(int(ids.size), dtype=np.float64),
        "best_v": np.zeros(int(ids.size), dtype=np.float64),
        "best_merit": float(np.linalg.norm(entrance_row)),
        "improved": False,
        "in_jac": False,
    }
    status = "linear_step"
    pyro_cache: dict[float, NamedStageQ] = {}

    def unpack(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        vv = np.asarray(v, dtype=np.float64)
        z = z0.copy()
        if n_log:
            z[log_ids] = base * np.exp(vv[:n_log])
        for k, (i1, i2, total, frac) in enumerate(pairs):
            v_s = float(vv[n_log + 2 * k])
            v_f = float(vv[n_log + 2 * k + 1])
            amount = total * float(np.exp(v_s)) if total > 0.0 else max(v_s, 0.0)
            n1, n2 = inventories_from_sum_fraction(amount, frac + v_f, a1, a2)
            z[i1] = n1
            z[i2] = n2
        if zero_ids.size:
            z[zero_ids] = np.maximum(vv[n_log + 2 * n_pair:react_at], 0.0)
        if n_react:
            z[react_ids] = np.maximum(react_base + react_scale * vv[react_at:react_at + n_react], 0.0)
        return z

    def fun(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        calls["F"] += 1
        if deadline_s is not None and time.monotonic() - t0 > float(deadline_s):
            raise TimeoutError("log corrector budget")
        if not state["in_jac"]:
            vv = np.asarray(v, dtype=np.float64).copy()
            state["v"] = vv
        z = unpack(v)
        qq = _pyrolysis_q_cached(q, z, pyro_cache)
        raw = apply_tar_reaction_transport_rows(
            z, qq, coupled_residual_named(z, qq, rate_scale=1.0),
        )
        row = raw / scale
        if not state["in_jac"]:
            merit = float(np.linalg.norm(row))
            if merit < float(state["best_merit"]):
                state["best_merit"] = merit
                state["best_v"] = np.asarray(v, dtype=np.float64).copy()
                state["improved"] = True
        return row

    def full_coupled_column(v_plus: npt.NDArray[np.float64], v_minus: npt.NDArray[np.float64], h_step: float) -> npt.NDArray[np.float64]:
        """完整 ``P`` 的对数列。飞行时间和床面系数都随状态更新，不再另加链式项。"""
        def row_at(vv: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            if deadline_s is not None and time.monotonic() - t0 > float(deadline_s):
                raise TimeoutError("log corrector budget")
            calls["F"] += 1
            calls["P"] += 1
            t_p = time.monotonic()
            zz = unpack(vv)
            qq = outer_P_named(zz, q)
            calls["P_wall_s"] += time.monotonic() - t_p
            raw = apply_tar_reaction_transport_rows(
                zz, qq, coupled_residual_named(zz, qq, rate_scale=1.0),
            )
            return raw / scale

        return (row_at(v_plus) - row_at(v_minus)) / (2.0 * float(h_step))

    def jac(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        if max_j is not None and calls["J"] >= int(max_j):
            raise _InnerBudget()
        calls["J"] += 1
        out = np.zeros((z0.size, int(ids.size)), dtype=np.float64)
        vv = np.asarray(v, dtype=np.float64)
        state["in_jac"] = True
        coupled_columns: set[int] = set()
        try:
            for j in range(n_log):
                i = int(log_ids[j])
                temperature = int(i) == z0.size - 1 or (z0.size >= 46 and int(i) == z0.size - 2)
                if couple_trajectory_jacobian and _trajectory_log_state_index(i, int(z0.size)):
                    # 1e-5 与 1e-6 在 s=0.678 的自由空间温度列上同号且相对差小于 1%。
                    h = 1.0e-5 if temperature else 1.0e-4
                    vp = vv.copy()
                    vm = vv.copy()
                    vp[j] += h
                    vm[j] -= h
                    out[:, j] = full_coupled_column(vp, vm, h)
                    coupled_columns.add(int(j))
                    continue
                h = 1.0e-7 if temperature else 1.0e-4
                vp = vv.copy()
                vm = vv.copy()
                vp[j] += h
                vm[j] -= h
                out[:, j] = (fun(vp) - fun(vm)) / (2.0 * h)
            for k, (_i1, _i2, total, frac) in enumerate(pairs):
                j_s = n_log + 2 * k
                j_f = j_s + 1
                if total > 0.0:
                    h_s = 1.0e-4
                    vp = vv.copy()
                    vm = vv.copy()
                    vp[j_s] += h_s
                    vm[j_s] -= h_s
                    out[:, j_s] = (fun(vp) - fun(vm)) / (2.0 * h_s)
                else:
                    h_s = 1.0e-8
                    vp = vv.copy()
                    vp[j_s] += h_s
                    out[:, j_s] = (fun(vp) - fun(vv)) / h_s
                h_f = 1.0e-4
                if frac - h_f < -a1:
                    vp = vv.copy()
                    vp[j_f] += h_f
                    out[:, j_f] = (fun(vp) - fun(vv)) / h_f
                elif frac + h_f > a2:
                    vm = vv.copy()
                    vm[j_f] -= h_f
                    out[:, j_f] = (fun(vv) - fun(vm)) / h_f
                else:
                    vp = vv.copy()
                    vm = vv.copy()
                    vp[j_f] += h_f
                    vm[j_f] -= h_f
                    out[:, j_f] = (fun(vp) - fun(vm)) / (2.0 * h_f)
            for j, _i in enumerate(zero_ids):
                col = n_log + 2 * n_pair + j
                vp = vv.copy()
                vp[col] += 1.0e-4
                out[:, col] = (fun(vp) - fun(vv)) / 1.0e-4
            for j in range(n_react):
                col = react_at + j
                h_r = 1.0e-4
                vp = vv.copy()
                vp[col] += h_r
                down = float(react_base[j] + react_scale[j] * (float(vv[col]) - h_r))
                if down > 0.0:
                    vm = vv.copy()
                    vm[col] -= h_r
                    out[:, col] = (fun(vp) - fun(vm)) / (2.0 * h_r)
                else:
                    out[:, col] = (fun(vp) - fun(vv)) / h_r
        finally:
            state["in_jac"] = False
        n_res, n_face = _add_exchange_chain(
            out, unpack(vv), q, ids, n_log, log_ids, scale, pyro_cache,
            skip_columns=coupled_columns,
        )
        calls["chain"] = float(calls.get("chain", 0.0)) + n_res
        calls["bed_face"] = float(calls.get("bed_face", 0.0)) + n_face
        return out

    lo = np.full(int(ids.size), -float(log_span), dtype=np.float64)
    hi = np.full(int(ids.size), float(log_span), dtype=np.float64)
    temp_span = float(log_span) * (1.0e-7 / 1.0e-4) * float(temperature_radius_scale)
    n_t = 2 if z0.size >= 46 else 1
    for j, i in enumerate(log_ids):
        if int(i) >= z0.size - n_t:
            lo[j] = -temp_span
            hi[j] = temp_span
    for k, (_i1, _i2, total, frac) in enumerate(pairs):
        j_s = n_log + 2 * k
        j_f = j_s + 1
        if total > 0.0:
            lo[j_s] = -float(log_span)
            hi[j_s] = float(log_span)
        else:
            lo[j_s] = 0.0
            hi[j_s] = hi_lin
        lo[j_f] = max(-a1 - frac, -float(ratio_span))
        hi[j_f] = min(a2 - frac, float(ratio_span))
    if zero_ids.size:
        lo[n_log + 2 * n_pair:react_at] = 0.0
        hi[n_log + 2 * n_pair:react_at] = hi_lin
    if n_react:
        lo[react_at:react_at + n_react] = np.asarray(cyclone["lo"], dtype=np.float64)
        hi[react_at:react_at + n_react] = np.asarray(cyclone["hi"], dtype=np.float64)
    v_out = state["best_v"]
    checkpoint = "entrance"
    linear_status = -1
    try:
        if deadline_s is not None and float(deadline_s) <= 0.0:
            raise TimeoutError("log corrector budget")
        if (
            isinstance(cached, dict)
            and int(cached.get("n", -1)) == int(z0.size)
            and bool(cached.get("trajectory_coupled")) == bool(couple_trajectory_jacobian)
            and np.asarray(cached["ids"]).shape == ids.shape
            and np.array_equal(cached["ids"], ids)
            and np.array_equal(np.asarray(cached.get("log_state_ids")), log_ids)
            and np.array_equal(np.asarray(cached.get("react_state_ids")), react_ids)
        ):
            matrix = np.asarray(cached["J"], dtype=np.float64)
            entrance_row = np.asarray(cached["entrance_row"], dtype=np.float64)
        else:
            matrix = jac(np.zeros(int(ids.size), dtype=np.float64))
            if stats is not None:
                stats["jacobian_cache"] = {
                    "n": int(z0.size),
                    "trajectory_coupled": bool(couple_trajectory_jacobian),
                    "ids": np.asarray(ids, dtype=np.int64).copy(),
                    "log_state_ids": np.asarray(log_ids, dtype=np.int64).copy(),
                    "react_state_ids": np.asarray(react_ids, dtype=np.int64).copy(),
                    "J": matrix,
                    "entrance_row": np.asarray(entrance_row, dtype=np.float64).copy(),
                }
        if residual_offset is not None:
            entrance_row = entrance_row + linear_tar_row_map(
                np.asarray(residual_offset, dtype=np.float64),
            ) / scale
        quality: dict[str, Any] = {}
        deadline_at = None if deadline_s is None else t0 + float(deadline_s)
        margin_limits = _settling_margin_limits(
            unpack(np.zeros(int(ids.size), dtype=np.float64)),
            q, unpack, int(ids.size), settling_margin_eta,
        )
        v_out, pred_raw, linear_status = _gate_progress_step(
            matrix, entrance_row, scale, lo, hi, quality, deadline_at=deadline_at,
            energy_guard_W=float(energy_guard_W), linear_abs_limits=margin_limits,
        )
        checkpoint = "linear_step"
        if quality.get("deadline_hit"):
            status = "deadline"
        if stats is not None:
            stats["predicted_raw"] = pred_raw
            stats["objective"] = "gate_ratio_phi"
            stats["linear_status"] = linear_status
            stats["active_bound_count"] = int(quality.get("active_bound_count", 0))
            stats["v_abs_max"] = float(np.max(np.abs(v_out))) if v_out.size else 0.0
            stats["linear_v"] = np.asarray(v_out, dtype=np.float64).copy()
            stats["subproblem_solved"] = bool(quality.get("solved"))
            stats["solver_success"] = quality.get("success")
            stats["solver_status"] = quality.get("status")
            stats["solver_message"] = quality.get("message")
            stats["solver_nit"] = quality.get("nit")
            stats["max_nd_violation"] = quality.get("max_nd_violation")
            stats["projected_gradient_norm"] = quality.get("projected_gradient_norm")
            stats["inexact_descent"] = bool(quality.get("inexact_descent"))
            stats["fallback_success"] = quality.get("fallback_success")
            stats["fallback_status"] = quality.get("fallback_status")
            stats["fallback_message"] = quality.get("fallback_message")
            stats["fallback_feasible"] = quality.get("fallback_feasible")
            stats["deadline_hit"] = bool(quality.get("deadline_hit"))
            temp_js = [j for j, i in enumerate(log_ids) if int(i) >= z0.size - n_t]
            width = np.maximum(hi - lo, 1.0e-12)
            stats["temperature_bound_active"] = bool(any(
                float(v_out[j]) <= float(lo[j]) + 1.0e-8 * float(width[j])
                or float(v_out[j]) >= float(hi[j]) - 1.0e-8 * float(width[j])
                for j in temp_js
            )) if v_out.size else False
    except TimeoutError:
        v_out = np.asarray(state["best_v"], dtype=np.float64)
        status = "deadline"
        checkpoint = "best_inner_candidate" if state["improved"] else "entrance"
    except _InnerBudget:
        v_out = np.asarray(state["best_v"], dtype=np.float64)
        status = "jacobian_budget"
        checkpoint = "best_inner_candidate" if state["improved"] else "entrance"
    finally:
        if stats is not None:
            stats["fd_residual_calls"] = float(calls["F"])
            stats["chain_residual_calls"] = float(calls["chain"])
            stats["bed_face_calls"] = float(calls["bed_face"])
            stats["P_calls"] = float(stats.get("P_calls", 0.0)) + float(calls["P"])
            stats["trajectory_P_wall_s"] = float(stats.get("trajectory_P_wall_s", 0.0)) + float(calls["P_wall_s"])
            stats["F_calls"] = float(stats.get("F_calls", 0.0)) + float(calls["F"]) + float(calls["chain"])
            stats["J_calls"] = float(stats.get("J_calls", 0.0)) + float(calls["J"])
            stats["inner_status"] = status
            stats["checkpoint"] = checkpoint
            stats["linear_status"] = linear_status
            stats["inner_merit_entrance"] = float(np.linalg.norm(entrance_row))
            stats["inner_merit_best"] = float(state["best_merit"])
            stats["inner_merit_best_is_a_candidate"] = bool(state["improved"])
            stats["linear_zero_ids"] = [int(i) for i in zero_ids]
            stats["cyclone_reactivation_ids"] = [int(i) for i in react_ids]
            stats["cyclone_reactivation_columns"] = [int(react_at + j) for j in range(n_react)]
            stats["cyclone_reactivation_scale_kg"] = [float(v) for v in react_scale]
    return unpack(v_out)


def _pyrolysis_q_cached(
    q: NamedStageQ,
    z: npt.NDArray[np.float64],
    cache: dict[float, NamedStageQ],
) -> NamedStageQ:
    """热解映射只随床温变化。同一温度的差商列复用，不重算 DAEM。"""
    if int(q.stage) < 4:
        return q
    index = z.size - 2 if int(q.stage) >= 9 else z.size - 1
    T = float(z[index])
    hit = cache.get(T)
    if hit is None:
        hit = refresh_pyrolysis_named(z, q)
        cache[T] = hit
    return hit


def _pyrolysis_temperature_column(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
    h_K: float = 1.0,
) -> npt.NDArray[np.float64]:
    """固定状态，只让热解映射随床温变化的残差列 [与 F 同单位 / K]。"""
    z = np.asarray(x, dtype=np.float64).copy()
    index = z.size - 2 if int(q.stage) >= 9 else z.size - 1
    T = float(z[index])

    def at(temp: float) -> npt.NDArray[np.float64]:
        zz = z.copy()
        zz[index] = float(temp)
        return coupled_residual_named(z, refresh_pyrolysis_named(zz, q), rate_scale=1.0)

    return (at(T + float(h_K)) - at(T - float(h_K))) / (2.0 * float(h_K))


def _linear_undeveloped_bubble_tar_ids(
    x: npt.NDArray[np.float64],
    f: npt.NDArray[np.float64],
) -> npt.NDArray[np.intp]:
    """气泡焦油仍在 1e-12 mol/s 以下、且该行残差明显不为零时走线性非负坐标。

    对数跨度 4 无法把 1e-18 mol/s 抬到残差所要求的 1e-5 mol/s。
    不把自由空间氧从对数坐标里拿掉。
    """
    z = np.asarray(x, dtype=np.float64)
    row = np.asarray(f, dtype=np.float64)
    ids = [
        i for i in (42, 43)
        if i < int(z.size) and 0.0 <= float(z[i]) <= 1.0e-12 and abs(float(row[i])) > 1.0e-9
    ]
    return np.asarray(ids, dtype=np.intp)


def _linear_bubble_zero_ids(
    x: npt.NDArray[np.float64],
    f: npt.NDArray[np.float64],
) -> npt.NDArray[np.intp]:
    """气泡气体中仍为零、且该行残差不为零的下标。初值保持 0，由联合步移动。"""
    z = np.asarray(x, dtype=np.float64)
    row = np.asarray(f, dtype=np.float64)
    ids = [
        i for i in range(7, min(14, int(z.size)))
        if float(z[i]) <= 0.0 and abs(float(row[i])) > 0.0
    ]
    return np.asarray(ids, dtype=np.intp)


def _e8_positive_log_ids(x: npt.NDArray[np.float64]) -> npt.NDArray[np.intp]:
    """严格为正的未知量都进入对数坐标，不按绝对阈值剔除自由空间氧。"""
    return np.flatnonzero(np.asarray(x, dtype=np.float64) > 0.0)


def _seed_downstream_o2(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """氧仍非正、且对应行需要修正时，给一个正初值。已为正的值保持不动。"""
    z = np.asarray(x, dtype=np.float64).copy()
    if int(z.size) < 36:
        return z
    if float(z[21]) <= 0.0:
        z = _balance_freeboard_o2(q, z)
        if float(z[21]) <= 0.0:
            f21 = float(coupled_residual_named(z, q, rate_scale=1.0)[21])
            if abs(f21) > 1.0e-12:
                z[21] = 1.0e-40
    f = coupled_residual_named(z, q, rate_scale=1.0)
    if float(z[28]) <= 0.0 and abs(float(f[28])) > 1.0e-9:
        z[28] = float(z[21]) if float(z[21]) > 0.0 else 1.0e-40
    if float(z[35]) <= 0.0 and abs(float(f[35])) > 1.0e-9:
        z[35] = float(z[28]) if float(z[28]) > 0.0 else 1.0e-40
    return z


def _balance_freeboard_o2(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """只动自由板氧，使该行过零。床层氧和旋风氧保持原未知量。

    精确 0 时自由板均相源关掉，行残差等于床层氧之和的相反数。
    正的流量一出现，R5/R12 就以极陡的斜率配平这一行，根在 1e-35 mol/s 附近。
    """
    z = np.asarray(x, dtype=np.float64).copy()
    if int(z.size) < 22:
        return z

    def row(v: float) -> float:
        w = z.copy()
        w[21] = float(v)
        return float(coupled_residual_named(w, q, rate_scale=1.0)[21])

    if float(z[21]) > 0.0:
        return z
    lo, hi = 1.0e-60, 1.0e-28
    f_lo, f_hi = row(lo), row(hi)
    expands = 0
    while f_lo * f_hi > 0.0 and expands < 8:
        if abs(f_lo) <= abs(f_hi):
            lo *= 1.0e-4
        else:
            hi *= 1.0e4
        f_lo, f_hi = row(lo), row(hi)
        expands += 1
    if f_lo == 0.0:
        z[21] = lo
        return z
    if f_lo * f_hi > 0.0:
        return z
    flo = f_lo
    for _ in range(48):
        mid = (lo * hi) ** 0.5
        fm = row(mid)
        if flo * fm <= 0.0:
            hi = mid
        else:
            lo = mid
            flo = fm
    z[21] = (lo * hi) ** 0.5
    return z


def _newton_at_frozen_q(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
    max_iter: int,
) -> Any:
    z0 = np.asarray(x, dtype=np.float64)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return coupled_residual_named(z, q, rate_scale=1.0)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return _fd_jacobian_trace(res, z)

    return wirsum_1997_strict_newton(
        res, jac, z0, max_iter=max_iter, tol_merit=1.0e-9, tol_step=0.0,
        trial_in_domain=_in_domain,
        linear_solve_fn=_e8_active_step,
    )


def _picard_at_hydro_scale(
    nq7: NamedStageQ,
    nq8: NamedStageQ,
    x: npt.NDArray[np.float64],
    hydro_scale: float,
    max_iter: int,
    *,
    fate_scale: float | None = None,
) -> Any:
    """冻结 ``q(s)`` 求根，再按 ``P(x)`` 刷新。复合 ``H`` 写入 ``last.F``。

    ``fate_scale is None`` 时与 ``hydro_scale`` 同步。先走水力时传 ``0``。
    """
    z = np.asarray(x, dtype=np.float64).copy()
    s = float(hydro_scale)
    last = None
    h = residual_e7_plus_scaled_e8(
        z, nq7, nq8, hydro_scale=s, fate_scale=fate_scale
    )
    for p_i in range(_N_PICARD):
        q = blend_hydro_named_q(
            outer_P_named(z, nq7),
            outer_P_named(z, nq8),
            hydro_scale=s,
            fate_scale=fate_scale,
        )
        # 对数坐标先走。原坐标牛顿会在事后裁剪里把步做成上升方向。
        z_c = _e8_log_scaled_correct(q, z)
        h_c = residual_e7_plus_scaled_e8(
            z_c, nq7, nq8, hydro_scale=s, fate_scale=fate_scale
        )
        if bool(named_endpoint_blocks(h_c)["passed"]) or (
            residual_block_merit(h_c) < residual_block_merit(h)
        ):
            z = z_c
            h = h_c
        if bool(named_endpoint_blocks(h)["passed"]):
            last = WirsumStrictSolveResult(
                x=z, F=h, converged=True, iterations=0,
                residual_l2_history=[float(np.linalg.norm(h))],
                failure_reason="block_gates",
            )
        else:
            last = _newton_at_frozen_q(q, z, max_iter)
            z_n = np.asarray(last.x, dtype=np.float64)
            h_n = residual_e7_plus_scaled_e8(
                z_n, nq7, nq8, hydro_scale=s, fate_scale=fate_scale
            )
            if residual_block_merit(h_n) < residual_block_merit(h):
                z = z_n
                h = h_n
        if not bool(named_endpoint_blocks(h)["passed"]):
            z_b = _newton_gas_energy_reserved(q, z)
            h_b = residual_e7_plus_scaled_e8(
                z_b, nq7, nq8, hydro_scale=s, fate_scale=fate_scale
            )
            if residual_block_merit(h_b) < residual_block_merit(h):
                z = z_b
                h = h_b
        last.x = z
        last.F = h
        if bool(named_endpoint_blocks(h)["passed"]):
            break
        if int(last.iterations) <= 1 and p_i > 0 and (
            not bool(named_endpoint_blocks(h)["passed"])
        ) and float(np.linalg.norm(h)) > 1.0e-3:
            break
    assert last is not None
    return last


def e8_source_snapshot(
    x: npt.NDArray[np.float64],
    nq7: NamedStageQ,
    nq8: NamedStageQ,
) -> dict[str, Any]:
    z8 = as_e8_state(np.asarray(x, dtype=np.float64))
    z7 = np.concatenate([z8[:-3], z8[-1:]])
    f7 = residual_named_e7(z7, nq7)
    f8 = residual_named_e8(z8, nq8)
    d_gas = f8[:14] - f7[:14]
    dF_energy = float(f8[-1] - f7[-1])
    q7 = outer_P_named(z7, nq7)
    q8 = outer_P_named(z8, nq8)
    return {
        "dF_l2": float(np.linalg.norm(np.concatenate([d_gas, [dF_energy]]))),
        "dF_energy": dF_energy,
        "F7_l2": float(np.linalg.norm(f7)),
        "F8_l2": float(np.linalg.norm(f8)),
        "T_K": float(z8[-1]),
        "T_ref_K": float(GEO_FB_RHO_G_T_REF_K),
        "rho_g7": float(q7.energy.species.geo_fb.rho_g),
        "rho_g8": float(q8.energy.species.geo_fb.rho_g),
        "K_eject7": np.asarray(q7.energy.species.K_eject, dtype=np.float64).tolist(),
        "K_eject8": np.asarray(q8.energy.species.K_eject, dtype=np.float64).tolist(),
        "w_leave7": np.asarray(q7.energy.species.w_leave, dtype=np.float64).tolist(),
        "w_leave8": np.asarray(q8.energy.species.w_leave, dtype=np.float64).tolist(),
        "cyclone_k_unchanged": q8.energy.species.k == nq8.energy.species.k,
    }


def solve_named_e8_from_e7(
    *,
    Q_W: float = 0.0,
    max_iter: int = 35,
    continuation: bool = True,
    hydro_scale: float | None = None,
    e7: dict[str, Any] | None = None,
    r11: dict[str, Any] | None = None,
    r10: dict[str, Any] | None = None,
    e6: dict[str, Any] | None = None,
    advance_mode: AdvanceMode | None = None,
    s_min: float = 1.0e-16,
    deadline_s: float | None = None,
) -> dict[str, Any]:
    """从 F_E7 根做 F_E8 水力延拓。默认 ``Q_W=0``。"""
    nq6 = freeze_e6_q(Q_W=Q_W)
    nq_r10 = freeze_e6_r10_q(Q_W=Q_W)
    nq_r11 = freeze_e6_r10_r11_q(Q_W=Q_W)
    nq7 = freeze_e7_q(Q_W=Q_W)
    nq8 = freeze_e8_q(Q_W=Q_W)
    if e7 is None:
        e7 = solve_named_e7_from_r11(
            Q_W=Q_W, max_iter=max_iter, r11=r11, r10=r10, e6=e6
        )
    if advance_mode is None:
        advance_mode = "root_correcting"
    x_e7 = np.asarray(e7["x"], dtype=np.float64).copy()
    snap0 = e8_source_snapshot(x_e7, nq7, nq8)
    x0 = as_e8_state(x_e7)
    # 正的气泡焦油只是避开 sqrt(C_tar) 奇点的初值，不是根，也不是给料。
    if int(x0.size) >= 45:
        if float(x0[42]) <= 0.0:
            x0[42] = 1.2e-5
        if float(x0[43]) <= 0.0:
            x0[43] = 1.8e-5
    cert0 = certify_named_full_rate_endpoint(x0, nq8, coupling="composite_P")
    # (hydro_scale, fate_scale)。None 表示与水力同步。
    # 默认先把气速/密度/K_bd 开到 1（命运保持 E7 的 w=[1,0]），再开 Eq. 3.72。
    if hydro_scale is not None:
        targets: list[tuple[float, float | None]] = [(float(hydro_scale), None)]
    elif continuation:
        targets = [(0.0, 0.0)] + [(float(s), 0.0) for s in _E8_HYDRO_CONTINUATION]
        targets += [(1.0, float(s)) for s in _E8_HYDRO_CONTINUATION]
    else:
        targets = [(1.0, 1.0)]
    history: list[dict[str, Any]] = []
    x_cur = x0.copy()
    last = None
    h_ok = 0.0
    f_ok = 0.0
    branch_tracked = True
    shrinks = 0
    i = 0
    timed_out = False
    t_deadline = None if deadline_s is None else time.monotonic() + float(deadline_s)
    while i < len(targets):
        if t_deadline is not None and time.monotonic() > t_deadline:
            timed_out = True
            break
        s_h, s_f = targets[i]
        at_end = abs(s_h - 1.0) < 1.0e-15 and (
            s_f is None or abs(float(s_f) - 1.0) < 1.0e-15
        )
        last = _picard_at_hydro_scale(
            nq7, nq8, x_cur, s_h, max_iter, fate_scale=s_f
        )
        stop_s = classify_newton_stop(last, model_is_full_rate=at_end)
        J = jacobian_e7_plus_scaled_e8(
            last.x, nq6, nq_r10, nq_r11, nq7, nq8,
            hydro_scale=s_h, fate_scale=s_f,
        )
        info = named_r12_scaled_rank(J, last.F, last.x)
        passed = stage_passed(
            stop_s, last.x, advance_mode=advance_mode,
            accepted_alphas=list(last.accepted_alphas),
        ) or bool(named_endpoint_blocks(last.F)["passed"])
        if at_end:
            cert_s = certify_named_full_rate_endpoint(
                last.x, nq8, coupling="composite_P"
            )
            if cert_s["endpoint_certified"]:
                passed = True
        history.append(
            {
                "hydro_scale": float(s_h),
                "fate_scale": None if s_f is None else float(s_f),
                "residual_l2": float(stop_s["residual_l2"]),
                "T_K": float(last.x[-1]),
                "rank_last": int(last.rank_last),
                "raw_rank": int(info["raw_rank"]),
                "equilibrated_rank": int(info["equilibrated_rank"]),
                "solver_stop_reason": stop_s["failure_reason"],
                "failure_reason": stop_s["failure_reason"],
                "converged_newton": bool(stop_s["converged_newton"]),
                "stage_passed": bool(passed),
            }
        )
        if passed:
            x_cur = np.asarray(last.x, dtype=np.float64)
            h_ok = float(s_h)
            f_ok = float(s_h if s_f is None else s_f)
            i += 1
            continue
        if advance_mode == "root_correcting":
            if shrinks >= _MAX_SHRINKS:
                branch_tracked = False
                break
            shrink_hydro = s_f is None or abs(float(s_h) - h_ok) > 1.0e-15
            if shrink_hydro:
                s_mid = shrink_scale(h_ok, s_h)
                stalled = (
                    abs(s_h - h_ok) < float(s_min) or abs(s_mid - s_h) < 1.0e-18
                )
                nxt = (s_mid, s_f)
            else:
                s_mid = shrink_scale(f_ok, float(s_f))
                stalled = (
                    abs(float(s_f) - f_ok) < float(s_min)
                    or abs(s_mid - float(s_f)) < 1.0e-18
                )
                nxt = (s_h, s_mid)
            if stalled:
                branch_tracked = False
                break
            targets.insert(i, nxt)
            shrinks += 1
            continue
        branch_tracked = False
        i += 1

    assert last is not None
    s_last = float(history[-1]["hydro_scale"]) if history else 1.0
    fate_last = history[-1]["fate_scale"] if history else 1.0
    full_rate = abs(s_last - 1.0) < 1.0e-15 and (
        fate_last is None or abs(float(fate_last) - 1.0) < 1.0e-15
    )
    last.F = residual_named_e8(last.x, nq8)
    step_ev = dict(_E8_STEP_EVIDENCE)
    cert = certify_named_full_rate_endpoint(last.x, nq8, coupling="composite_P")
    f0 = residual_named_e8(x0, nq8)
    last, cert, returned_state_source = apply_e8_certified_rollback(
        last, x0, f0, cert, cert0,
    )
    if "rollback" not in step_ev:
        step_ev["rollback"] = {}
    step_ev["rollback"]["returned_state_source"] = returned_state_source
    stop = classify_newton_stop(last, model_is_full_rate=full_rate)
    info = named_r12_scaled_rank(
        jacobian_named_e8(last.x, nq6, nq_r10, nq_r11, nq7, nq8),
        last.F, last.x,
    )
    nqP = outer_P_named(last.x, nq8)
    u = unpack_named(last.x, 8)
    omitted = omitted_named(8)
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": float(last.x[-1]),
        "T_e7_K": float(e7["T_K"]),
        "M_ash": float(u["M_ash"]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "N_fb_O2": float(u["N_fb"][IDX["O2"]]),
        "K_eject": np.asarray(nqP.energy.species.K_eject, dtype=np.float64).tolist(),
        "w_leave": np.asarray(nqP.energy.species.w_leave, dtype=np.float64).tolist(),
        "rho_g": float(nqP.energy.species.geo_fb.rho_g),
        "Q_W": float(nq8.Q_W),
        "Q_W_is_fitted_to_1100K": False,
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": f_name(8),
        "n_unknowns": n_unknowns(8),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "e7_physical_convergence": bool(e7.get("physical_convergence", False)),
        "e7_endpoint_certified": bool(e7.get("endpoint_certified", False)),
        "e8_at_e7_root": snap0,
        "solver_stop_reason": stop["failure_reason"],
        "endpoint_cert": {k: v for k, v in cert.items() if k != "F"},
        "endpoint_certified": bool(cert["endpoint_certified"]),
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
        "production_hamel_outer_refresh": False,
        "r10_tar_oxidation": True,
        "r11_tar_conversion": True,
        "freeboard_reaction": True,
        "K_eject_in_this_P": True,
        "finite_hydro_P": True,
        "cyclone_K_update_P": False,
        "x_nonneg_projection": False,
        "coupling_mode": "composite_P",
        "advance_mode": str(advance_mode),
        "branch_tracked": bool(
            branch_tracked
            and bool(history)
            and full_rate
            and bool(history[-1]["stage_passed"])
        ),
        "continuation": history,
        "timed_out": bool(timed_out),
        "bound_feasible_direction": True,
        "returned_state_source": returned_state_source,
        "linear_step_last": step_ev,
        "omitted_named": omitted,
        "T_ref_K": float(GEO_FB_RHO_G_T_REF_K),
    }
