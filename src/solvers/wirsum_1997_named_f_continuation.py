"""具名核速率延拓的共享推进 / 诊断 / 端点认证。不是 Hamel 方程。

Source: docs/wirsum_1997_named_f_rate_continuation.md;
docs/audit_reports/rate_continuation_review_2026-09-15.md §7
"""

from __future__ import annotations

from typing import Any, Callable, Literal

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    outer_P_named,
)
from src.solvers.wirsum_1997_named_f_r12 import (
    _TRACE_DIAG,
    IDX,
    named_r12_scaled_rank,
)

AdvanceMode = Literal["root_correcting", "heuristic_search"]
CouplingMode = Literal["composite_P", "frozen_q"]

ENDPOINT_L2: float = 1.0e-6
_LAMBDA_MIN: float = 1.0 / 4096.0
_SOLID_IDX: tuple[int, ...] = (14, 15, 37, 38)


def named_endpoint_blocks(F: npt.NDArray[np.float64]) -> dict[str, float | bool]:
    """原单位分块门。气体/固体/焦油 L2 < ``ENDPOINT_L2``，能量绝对值 < 1 W。

    不把 mol/s、kg/s 和 W 合成一个无权 L2 再和 ``1e-6`` 比较。
    """
    f = np.asarray(F, dtype=np.float64)
    n = int(f.size)
    gas = f[0:14]
    if n >= 37:
        gas = np.concatenate([gas, f[16:37]])
    elif n >= 23:
        gas = np.concatenate([gas, f[16:23]])
    solid_idx = [i for i in _SOLID_IDX if i < n - 1]
    if n >= 41:
        solid_idx.append(39)
    if n >= 50:
        solid_idx.extend([46, 47])
    solid = f[solid_idx] if solid_idx else np.zeros(0, dtype=np.float64)
    if n >= 48:
        tar = f[40:46]
    elif n >= 45:
        tar = f[40:44]
    elif n >= 43:
        tar = f[40:42]
    else:
        tar = np.zeros(0, dtype=np.float64)
    if n >= 46:
        energy_rows = f[-2:]
        energy = float(np.sum(energy_rows))
        energy_ok = bool(np.max(np.abs(energy_rows)) < 1.0)
    else:
        energy = float(f[-1])
        energy_ok = bool(abs(energy) < 1.0)
    gas_l2 = float(np.linalg.norm(gas))
    solid_l2 = float(np.linalg.norm(solid))
    tar_l2 = float(np.linalg.norm(tar))
    passed = bool(
        gas_l2 < ENDPOINT_L2
        and solid_l2 < ENDPOINT_L2
        and tar_l2 < ENDPOINT_L2
        and energy_ok
    )
    return {
        "passed": passed,
        "gas_l2": gas_l2,
        "solid_l2": solid_l2,
        "tar_l2": tar_l2,
        "energy_W": energy,
    }


def named_in_domain(x: npt.NDArray[np.float64]) -> bool:
    z = np.asarray(x, dtype=np.float64)
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def bound_feasible_direction(
    x: npt.NDArray[np.float64],
    dx: npt.NDArray[np.float64],
    *,
    lambda_min: float = _LAMBDA_MIN,
) -> tuple[npt.NDArray[np.float64], bool]:
    """把会在最小 α 下越界的分量从方向里拿掉。改的是步，不是 F。"""
    z = np.asarray(x, dtype=np.float64)
    p = np.asarray(dx, dtype=np.float64).copy()
    changed = False
    a_min = float(lambda_min)
    for i in range(int(p.size) - 1):
        if p[i] >= 0.0:
            continue
        if float(z[i]) <= 0.0 or float(z[i]) + a_min * float(p[i]) < 0.0:
            p[i] = 0.0
            changed = True
    return p, changed


def diagnose_linear_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    dx: npt.NDArray[np.float64],
    *,
    used_trace_diag: bool,
    J_modified: npt.NDArray[np.float64] | None = None,
) -> dict[str, Any]:
    """原 J 与补对角后的线性残差。补对角不是等价缩放。"""
    j = np.asarray(J, dtype=np.float64)
    f = np.asarray(F, dtype=np.float64)
    p = np.asarray(dx, dtype=np.float64)
    r_raw = j @ p + f
    out: dict[str, Any] = {
        "used_trace_diag": bool(used_trace_diag),
        "trace_diag": float(_TRACE_DIAG) if used_trace_diag else 0.0,
        "linear_residual_raw_l2": float(np.linalg.norm(r_raw)),
        "directional_phi_dot": float(np.dot(f, j @ p)),
        "minus_f_l2_sq": float(-np.dot(f, f)),
    }
    if J_modified is not None:
        r_mod = np.asarray(J_modified, dtype=np.float64) @ p + f
        out["linear_residual_modified_l2"] = float(np.linalg.norm(r_mod))
    return out


def named_linear_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    x: npt.NDArray[np.float64],
    *,
    extra_diag_species: tuple[str, ...] = (),
    phase_offsets: tuple[int, ...] = (0, 7),
) -> tuple[npt.NDArray[np.float64], int, float, dict[str, Any]]:
    """微量组分尺度 + 行列均衡。亏秩时才补对角，并记入诊断。"""
    j0 = np.asarray(J, dtype=np.float64)
    info0 = named_r12_scaled_rank(j0, F, x)
    n = int(info0["n"])
    used = False
    j_mod: npt.NDArray[np.float64] | None = None
    info = info0
    rank = int(info["equilibrated_rank"])
    if rank < n:
        used = True
        j_mod = j0.copy()
        species = extra_diag_species if extra_diag_species else ("H2",)
        for off in phase_offsets:
            for sp in species:
                i = int(off) + IDX[sp]
                if 0 <= i < n - 1:
                    j_mod[i, i] += _TRACE_DIAG
        info = named_r12_scaled_rank(j_mod, F, x)
        rank = int(info["equilibrated_rank"])
        if rank < n:
            raise np.linalg.LinAlgError("named scaled Jacobian singular after trace diag")
    w = np.linalg.solve(info["J_eq"], -(info["dr"] * info["F_hat"]))
    dx = info["x_scale"] * (info["dc"] * w)
    diag = diagnose_linear_step(
        j0, F, dx, used_trace_diag=used, J_modified=j_mod,
    )
    diag["raw_rank"] = int(info0["raw_rank"])
    diag["equilibrated_rank_original"] = int(info0["equilibrated_rank"])
    diag["equilibrated_rank_used"] = rank
    diag["J_modified"] = j_mod
    return dx, rank, float(info["equilibrated_condition"]), diag


def stage_passed(
    stop: dict[str, Any],
    x: npt.NDArray[np.float64],
    *,
    advance_mode: AdvanceMode,
    accepted_alphas: list[float],
) -> bool:
    if not named_in_domain(x):
        return False
    if advance_mode == "root_correcting":
        return float(stop["residual_l2"]) < ENDPOINT_L2
    return bool(accepted_alphas) or bool(stop["converged_newton"])


def shrink_scale(s_ok: float, s_fail: float) -> float:
    return 0.5 * (float(s_ok) + float(s_fail))


def certify_named_full_rate_endpoint(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    coupling: CouplingMode = "composite_P",
) -> dict[str, Any]:
    """不经延拓差值、不经 clip，直接装配目标 F。"""
    z = np.asarray(x, dtype=np.float64)
    in_dom = named_in_domain(z)
    nq_eval = outer_P_named(z, nq) if coupling == "composite_P" else nq
    f = coupled_residual_named(z, nq_eval, rate_scale=1.0)
    blocks = named_endpoint_blocks(f)
    n = int(f.size)
    gas = f[0:14]
    if n >= 37:
        gas = np.concatenate([gas, f[16:37]])
    elif n >= 23:
        gas = np.concatenate([gas, f[16:23]])
    solid_idx = [i for i in _SOLID_IDX if i < n - 1]
    if n >= 41:
        solid_idx.append(39)
    solid = f[solid_idx] if solid_idx else np.zeros(0, dtype=np.float64)
    tar = f[40:44] if n >= 45 else (f[40:42] if n >= 43 else np.zeros(0, dtype=np.float64))
    energy = float(f[-1])
    f_frozen = coupled_residual_named(z, nq, rate_scale=1.0)
    q_delta = float(np.linalg.norm(f - f_frozen))
    gas_l2 = float(np.linalg.norm(gas))
    solid_l2 = float(np.linalg.norm(solid))
    tar_l2 = float(np.linalg.norm(tar))
    certified = bool(in_dom and blocks["passed"])
    return {
        "endpoint_certified": certified,
        "in_domain": in_dom,
        "coupling": str(coupling),
        "residual_l2": float(np.linalg.norm(f)),
        "gas_mol_s_l2": gas_l2,
        "gas_mol_s_linf": float(np.max(np.abs(gas))) if gas.size else 0.0,
        "solid_kg_s_l2": solid_l2,
        "solid_kg_s_linf": float(np.max(np.abs(solid))) if solid.size else 0.0,
        "tar_mol_s_l2": tar_l2,
        "energy_W": energy,
        "energy_W_abs": abs(energy),
        "q_consistency_l2": q_delta,
        "n_unknowns": n,
        "enabled_rxns": tuple(sorted(nq.enabled_rxns)),
        "rate_scale": 1.0,
        "clip_applied": False,
        "F": f,
    }
