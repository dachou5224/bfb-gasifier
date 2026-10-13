"""P3/E8 共用的外层控制。基准只提供初值，不各自复制循环。

内层是一次有界线性步，雅可比含冻结水力项以及 K_eject、K_bd 的链式项。
接受看气体、固体、焦油的门比值。能量行低于 1 W 时允许变大，不能越过 1 W 还变差。
线性步本身不是根。

Source: docs/audit_reports/plan_ab_diagnosis_2026-09-28.md 工作包 1
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    inner_hydro_frozen_residual,
    outer_P_named,
)
from src.solvers.wirsum_1997_named_f_continuation import ENDPOINT_L2, named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e8 import (
    _EJECT_SUBSTANTIAL_KG_S,
    _e8_log_scaled_correct,
    freeboard_settling_margin_m,
    residual_block_merit,
)


def _count_residual(
    z: npt.NDArray[np.float64],
    q: NamedStageQ,
    stats: dict[str, Any] | None,
) -> npt.NDArray[np.float64]:
    if stats is not None:
        stats["F_calls"] = float(stats.get("F_calls", 0.0)) + 1.0
    return coupled_residual_named(z, q, rate_scale=1.0)


def bubble_source_zeros(
    x: npt.NDArray[np.float64],
    f: npt.NDArray[np.float64],
) -> list[int]:
    """气泡气体里仍为零、且该行残差不为零的下标。含 H2 与 CH4。"""
    z = np.asarray(x, dtype=np.float64)
    row = np.asarray(f, dtype=np.float64)
    if z.size < 14:
        return []
    return [
        i for i in range(7, 14)
        if float(z[i]) <= 0.0 and abs(float(row[i])) > 0.0
    ]


def joint_open_zeros(
    x: npt.NDArray[np.float64],
    f: npt.NDArray[np.float64],
    indices: list[int],
    q: NamedStageQ,
    stats: dict[str, Any] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], list[int]]:
    """同时给一组零变量赋值。只在完整复合残差下降时接受。

    单行残差不下降不能把变量留在边界上。计量耦合时要看联合步。
    """
    z = np.asarray(x, dtype=np.float64).copy()
    row = np.asarray(f, dtype=np.float64).copy()
    if not indices:
        return z, row, []
    phi0 = progress_phi(row)
    positive = float(np.sum(np.maximum(z, 0.0)))
    scale0 = 0.01 * positive if positive > 0.0 else 1.0e-8
    weights = np.array([max(abs(float(row[i])), 1.0e-12) for i in indices], dtype=np.float64)
    weights = weights / float(np.sum(weights))
    best_m = phi0
    best_z = z
    best_f = row
    opened: list[int] = []
    for factor in (0.01, 0.1, 1.0):
        trial = z.copy()
        for i, w in zip(indices, weights):
            trial[i] = scale0 * factor * float(w)
        f_try = _count_residual(trial, q, stats)
        phi_try = progress_phi(f_try)
        if phi_try < best_m:
            best_m = phi_try
            best_z = trial
            best_f = f_try
            opened = list(indices)
    return best_z, best_f, opened


_TRACE_O2_INDEX = (5, 12, 21)  # 床层悬浮相、气泡、自由空间 [mol/s]


def oxygen_rows_needing_restart(
    f: npt.NDArray[np.float64],
    predicted: npt.NDArray[np.float64] | None,
) -> list[int]:
    """氧行主导未过门的气体残差，且当前线性步修不动这些行。

    不按氧流量的绝对阈值触发。气体块已经低于门时返回空。
    """
    row = np.asarray(f, dtype=np.float64)
    n = int(row.size)
    gas = list(range(min(14, n)))
    if n > 16:
        gas.extend(range(16, min(37, n)))
    if not gas:
        return []
    gas_l2 = float(np.linalg.norm(row[gas]))
    if gas_l2 <= float(ENDPOINT_L2):
        return []
    present = [i for i in _TRACE_O2_INDEX if i < n]
    if not present:
        return []
    oxygen_l2 = float(np.linalg.norm(row[present]))
    if oxygen_l2 * oxygen_l2 < 0.5 * gas_l2 * gas_l2:
        return []
    pred = None if predicted is None else np.asarray(predicted, dtype=np.float64)
    chosen: list[int] = []
    for index in present:
        if abs(float(row[index])) < 0.2 * gas_l2:
            continue
        if pred is not None and pred.shape == row.shape:
            correction = abs(float(pred[index]) - float(row[index]))
            if correction >= 0.1 * abs(float(row[index])):
                continue
        chosen.append(int(index))
    return chosen


def bracket_positive_scalar(
    level: float,
    row_at,
    *,
    max_samples: int = 8,
) -> dict[str, Any]:
    """在正未知量的对数尺度上找唯一变号点。找不到就不抬到固定地板。

    ``row_at(level)`` 返回该未知量自己的原始残差。采样次数有上限。
    """
    amount = float(level)
    if not np.isfinite(amount) or amount <= 0.0:
        return {"status": "nonpositive", "samples": []}
    samples: list[tuple[float, float]] = []
    sample_limit = int(max_samples) + 20

    def evaluate(value: float) -> float | None:
        if len(samples) >= sample_limit:
            return None
        residual = float(row_at(float(value)))
        if not np.isfinite(residual):
            samples.append((float(value), residual))
            return None
        samples.append((float(value), residual))
        return residual

    origin = evaluate(amount)
    if origin is None:
        return {"status": "discontinuous", "samples": samples}
    for decade in (4.0, -4.0, 8.0, -8.0, 16.0, -16.0, 32.0, -32.0):
        if len(samples) >= int(max_samples):
            break
        evaluate(amount * float(10.0 ** decade))
    finite = [(value, residual) for value, residual in samples if np.isfinite(residual)]
    if len(finite) < len(samples):
        return {"status": "discontinuous", "samples": samples}
    finite.sort(key=lambda item: item[0])
    changes: list[tuple[float, float]] = []
    for (left, left_row), (right, right_row) in zip(finite, finite[1:]):
        if left_row == 0.0:
            return {"status": "root", "level": float(left), "samples": samples}
        if left_row * right_row < 0.0:
            changes.append((float(left), float(right)))
    if len(changes) == 0:
        return {"status": "no_bracket", "samples": samples}
    if len(changes) > 1:
        return {"status": "multiple_roots", "samples": samples}
    low, high = changes[0]
    from scipy.optimize import brentq

    def on_log(exponent: float) -> float:
        value = evaluate(float(10.0 ** exponent))
        if value is None:
            raise RuntimeError("scalar sample budget")
        return float(value)

    try:
        exponent = brentq(on_log, np.log10(low), np.log10(high), xtol=1.0e-8, maxiter=16)
    except (RuntimeError, ValueError):
        finite_now = [item for item in samples if np.isfinite(item[1])]
        if not finite_now:
            return {"status": "no_bracket", "samples": samples}
        best_level, best_row = min(finite_now, key=lambda item: abs(item[1]))
        if abs(best_row) <= max(1.0e-12, 1.0e-6 * abs(origin)):
            return {"status": "root", "level": float(best_level), "samples": samples}
        return {"status": "no_bracket", "samples": samples}
    return {"status": "root", "level": float(10.0 ** exponent), "samples": samples}


def trace_oxygen_candidate(
    x: npt.NDArray[np.float64],
    q: NamedStageQ,
    rows: list[int],
    stats: dict[str, Any] | None = None,
) -> tuple[npt.NDArray[np.float64] | None, list[dict[str, Any]]]:
    """冻结其余未知量，只为列出的正氧行找标量根。任一失败则不改状态。"""
    z = np.asarray(x, dtype=np.float64)
    trial = z.copy()
    report: list[dict[str, Any]] = []
    for index in rows:
        def row_at(level: float, index: int = int(index)) -> float:
            w = z.copy()
            w[index] = float(level)
            if stats is not None:
                stats["F_calls"] = float(stats.get("F_calls", 0.0)) + 1.0
            return float(coupled_residual_named(w, q, rate_scale=1.0)[index])

        found = bracket_positive_scalar(float(z[index]), row_at)
        found["index"] = int(index)
        report.append(found)
        if found.get("status") != "root":
            return None, report
        trial[index] = float(found["level"])
    return trial, report


def reactivate_zero_unknowns(
    q: NamedStageQ,
    x: npt.NDArray[np.float64],
    stats: dict[str, Any] | None = None,
) -> tuple[npt.NDArray[np.float64], list[int]]:
    """打开仍为零、但方程尚未平衡的未知量。整组接受后立刻更新残差。"""
    z = np.asarray(x, dtype=np.float64).copy()
    f = _count_residual(z, q, stats)
    opened: list[int] = []
    groups: list[tuple[slice, np.ndarray]] = []
    n_in_b = np.asarray(q.energy.species.N_in_b, dtype=np.float64)
    if z.size >= 14 and float(np.max(z[7:14])) <= 0.0 and float(np.max(n_in_b)) > 0.0:
        base = np.maximum(n_in_b, 0.0)
        for factor in (0.05, 0.2, 1.0):
            groups.append((slice(7, 14), base * factor))
    n_tar = q.energy.n_tar
    tar_fb = np.array(
        [float(n_tar.get("TAR1", 0.0)), float(n_tar.get("TAR2", 0.0))],
        dtype=np.float64,
    )
    if z.size >= 46 and float(np.max(z[44:46])) <= 0.0 and float(np.max(tar_fb)) > 0.0:
        groups.append((slice(44, 46), np.maximum(tar_fb, 0.0)))
    for sl, feed in groups:
        trial = z.copy()
        trial[sl] = feed
        f_try = _count_residual(trial, q, stats)
        phi_try = progress_phi(f_try)
        if phi_try < progress_phi(f):
            z = trial
            f = f_try
            opened.extend(range(sl.start, sl.stop))
    z, f, joint = joint_open_zeros(z, f, bubble_source_zeros(z, f), q, stats)
    opened.extend(i for i in joint if i not in opened)
    for sl in (slice(7, 14), slice(44, 46)):
        if z.size < sl.stop:
            continue
        block = z[sl]
        if float(np.max(block)) <= 0.0:
            continue
        floor = 0.01 * float(np.sum(np.maximum(block, 0.0)))
        for i in range(sl.start, sl.stop):
            if float(z[i]) > 0.0:
                continue
            best_v = 0.0
            best_abs = abs(float(f[i]))
            best_f = f
            for factor in (0.1, 1.0, 10.0):
                trial = z.copy()
                trial[i] = floor * factor
                f_try = _count_residual(trial, q, stats)
                row = abs(float(f_try[i]))
                if row < best_abs:
                    best_abs = row
                    best_v = float(trial[i])
                    best_f = f_try
            if best_v > 0.0:
                z[i] = best_v
                f = best_f
                opened.append(int(i))
    return z, opened


TEMPERATURE_DOMAIN = "temporary_solver_domain_250_2500_K"
_T_LO = 250.0  # [K] 临时求解域，不是物性适用域
_T_HI = 2500.0  # [K]


_ENERGY_GATE_W = 1.0  # [W] 与分块门相同，不是新的松弛
_RADIUS_MAX = 1.0


def material_gate_ratios(f: npt.NDArray[np.float64]) -> tuple[float, float, float]:
    """气体、固体、焦油的 L2 除以 ``1e-6``。能量不放进这个三元组。"""
    blocks = _block_report(f)
    tol = float(ENDPOINT_L2)
    return (
        float(blocks["gas_l2"]) / tol,
        float(blocks["solid_l2"]) / tol,
        float(blocks["tar_l2"]) / tol,
    )


def progress_ratios(f: npt.NDArray[np.float64]) -> tuple[float, ...]:
    """气体、固体、焦油和各能量行的门比值。能量用 1 W，和认证门同一尺度。"""
    gas, solid, tar = material_gate_ratios(f)
    n_t = 2 if int(np.asarray(f).size) >= 46 else 1
    energy = tuple(
        abs(float(v)) / _ENERGY_GATE_W for v in np.asarray(f, dtype=np.float64)[-n_t:]
    )
    return (gas, solid, tar, *energy)


def progress_phi(f: npt.NDArray[np.float64]) -> float:
    """门比值的总体违反度。毫瓦级能量在这里远小于尚未过门的物料块。"""
    ratios = progress_ratios(f)
    return float(np.sqrt(sum(float(v) * float(v) for v in ratios)))


def step_improves_progress(
    f_old: npt.NDArray[np.float64],
    f_new: npt.NDArray[np.float64],
) -> bool:
    """总体违反下降。已经低于门的块不能被抬过门；门上的块允许交换。"""
    old = progress_ratios(f_old)
    new = progress_ratios(f_new)
    if not (np.all(np.isfinite(old)) and np.all(np.isfinite(new))):
        return False
    for before, after in zip(old, new):
        if float(before) < 1.0 and float(after) > 1.0 + 1.0e-9:
            return False
    return progress_phi(f_new) < progress_phi(f_old) - 1.0e-9


def step_reduces_material_violation(
    f_old: npt.NDArray[np.float64],
    f_new: npt.NDArray[np.float64],
) -> bool:
    """三块物料门比值都不升，且最大违反下降。已低于 1 W 的能量行允许变大。"""
    old = material_gate_ratios(f_old)
    new = material_gate_ratios(f_new)
    if any(float(n) > float(o) + 1.0e-9 for n, o in zip(new, old)):
        return False
    if max(new) >= max(old) - 1.0e-9:
        return False
    n_t = 2 if int(np.asarray(f_old).size) >= 46 else 1
    old_e = np.asarray(f_old, dtype=np.float64)[-n_t:]
    new_e = np.asarray(f_new, dtype=np.float64)[-n_t:]
    for a, b in zip(old_e, new_e):
        if abs(float(b)) >= _ENERGY_GATE_W and abs(float(b)) > abs(float(a)) + 1.0e-9:
            return False
    return True


def _block_report(f: npt.NDArray[np.float64]) -> dict[str, Any]:
    blocks = named_endpoint_blocks(f)
    return {
        k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v))
        for k, v in blocks.items()
    }


def _state_admissible(
    x: npt.NDArray[np.float64],
    f: npt.NDArray[np.float64],
    hydro: Any,
) -> bool:
    """全步、回溯和零量打开共用：有限、在临时域内、输运闭合。"""
    cert = certify_controlled_state(x, f, hydro)
    return bool(cert["finite"] and cert["in_domain"] and cert["transport_closed"])


def transport_is_closed(hydro: Any) -> bool:
    """未判定份额还在时，输运没有闭合。不能把这种 q 当成认证稳态。"""
    if hydro is None:
        return False
    und = np.asarray(hydro.w_undecided, dtype=np.float64)
    return bool(np.all(und <= 1.0e-12))


def state_in_temporary_domain(x: npt.NDArray[np.float64]) -> bool:
    """250–2500 K 是临时求解域。通过它不等于物理认证。"""
    z = np.asarray(x, dtype=np.float64)
    if not np.all(np.isfinite(z)):
        return False
    n_t = 2 if z.size >= 46 else 1
    temps = z[-n_t:]
    flows = z[:-n_t] if n_t < z.size else z[:0]
    return bool(
        np.all(temps >= _T_LO) and np.all(temps <= _T_HI) and np.all(flows >= 0.0)
    )


def certify_controlled_state(
    x: npt.NDArray[np.float64],
    f: npt.NDArray[np.float64],
    hydro: Any = None,
) -> dict[str, Any]:
    """成功出口共用。残差小但出域、非有限或输运未闭合时不通过。"""
    z = np.asarray(x, dtype=np.float64)
    row = np.asarray(f, dtype=np.float64)
    finite = bool(np.all(np.isfinite(z)) and np.all(np.isfinite(row)))
    in_domain = finite and state_in_temporary_domain(z)
    blocks = _block_report(row) if finite else _block_report(np.full(z.shape, np.inf))
    closed = transport_is_closed(hydro)
    n_t = 2 if row.size >= 46 else 1
    energy_rows = [float(v) for v in row[-n_t:]] if finite else []
    passed = bool(finite and in_domain and closed and blocks["passed"])
    return {
        "passed": passed,
        "finite": finite,
        "in_domain": in_domain,
        "transport_closed": closed,
        "temperature_domain": TEMPERATURE_DOMAIN,
        "blocks": blocks,
        "energy_rows_W": energy_rows,
        "gate_ratio": {
            "gas": float(blocks["gas_l2"]) / 1.0e-6,
            "solid": float(blocks["solid_l2"]) / 1.0e-6,
            "tar": float(blocks["tar_l2"]) / 1.0e-6,
            "energy_rows": [abs(v) / 1.0 for v in energy_rows],
        },
    }


_BACKTRACK_ALPHAS = (0.5, 0.25, 0.1, 0.01, 0.001)


def log_blend_state(
    accepted: npt.NDArray[np.float64],
    proposed: npt.NDArray[np.float64],
    alpha: float,
) -> npt.NDArray[np.float64]:
    """正分量沿对数步缩到 ``alpha``，其余分量线性缩。"""
    base = np.asarray(accepted, dtype=np.float64)
    prop = np.asarray(proposed, dtype=np.float64)
    trial = base.copy()
    positive = (base > 0.0) & (prop > 0.0)
    trial[positive] = base[positive] * (prop[positive] / base[positive]) ** float(alpha)
    trial[~positive] = base[~positive] + float(alpha) * (prop[~positive] - base[~positive])
    return trial


def progress_backtrack(
    nq: NamedStageQ,
    accepted: npt.NDArray[np.float64],
    proposed: npt.NDArray[np.float64],
    f_accepted: npt.NDArray[np.float64],
    stats: dict[str, Any],
    extra_ok=None,
) -> dict[str, Any] | None:
    """全步把已过门的块顶出门时，沿同一对数方向找仍满足该规则的最大分数步。"""
    for alpha in (0.5, 0.25, 0.1):
        trial = log_blend_state(accepted, proposed, alpha)
        if not state_in_temporary_domain(trial):
            continue
        q_try = outer_P_named(trial, nq)
        stats["P_calls"] = float(stats.get("P_calls", 0.0)) + 1.0
        f_try = _count_residual(trial, q_try, stats)
        cert = certify_controlled_state(trial, f_try, q_try.hydro)
        if extra_ok is not None and not extra_ok(trial):
            continue
        if _state_admissible(trial, f_try, q_try.hydro) and step_improves_progress(f_accepted, f_try):
            return {
                "alpha": float(alpha),
                "x": trial,
                "q": q_try,
                "F": f_try,
                "cert": cert,
                "phi": progress_phi(f_try),
            }
    return None


def composite_backtrack(
    nq: NamedStageQ,
    accepted: npt.NDArray[np.float64],
    proposed: npt.NDArray[np.float64],
    f_accepted: npt.NDArray[np.float64],
    stats: dict[str, Any],
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """全步被拒后，沿同一方向用完整 P 试分数步。不把 0.25 写成固定松弛。"""
    base = np.asarray(accepted, dtype=np.float64)
    prop = np.asarray(proposed, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    chosen: dict[str, Any] | None = None
    for alpha in _BACKTRACK_ALPHAS:
        trial = base.copy()
        positive = (base > 0.0) & (prop > 0.0)
        trial[positive] = base[positive] * (prop[positive] / base[positive]) ** float(alpha)
        trial[~positive] = base[~positive] + float(alpha) * (prop[~positive] - base[~positive])
        if not state_in_temporary_domain(trial):
            rows.append({"alpha": float(alpha), "status": "out_of_domain"})
            continue
        q_try = outer_P_named(trial, nq)
        stats["P_calls"] = float(stats.get("P_calls", 0.0)) + 1.0
        f_try = _count_residual(trial, q_try, stats)
        merit = residual_block_merit(f_try)
        cert = certify_controlled_state(trial, f_try, q_try.hydro)
        improved = step_reduces_material_violation(f_accepted, f_try)
        row = {
            "alpha": float(alpha),
            "merit": merit,
            "material_ratios": material_gate_ratios(f_try),
            "passed": bool(cert["passed"]),
            "transport_closed": bool(cert["transport_closed"]),
            "material_improved": improved,
        }
        rows.append(row)
        if not cert["transport_closed"]:
            continue
        if bool(cert["passed"]) or improved:
            rank = max(material_gate_ratios(f_try))
            if chosen is None or rank < float(chosen["rank"]):
                chosen = {
                    "alpha": float(alpha),
                    "x": trial,
                    "q": q_try,
                    "F": f_try,
                    "merit": merit,
                    "cert": cert,
                    "rank": rank,
                }
    stats["backtrack"] = rows
    return chosen, rows


def shrink_trust(log_span: float, temperature_factor: float) -> tuple[float, float]:
    """拒绝后缩小对数步和温度倍数。下一次候选因此不同。"""
    span = max(0.25, float(log_span) * 0.5)
    factor = max(1.05, 1.0 + (float(temperature_factor) - 1.0) * 0.5)
    return span, factor


def _finish_oxygen_restart(
    nq: NamedStageQ,
    accepted: npt.NDArray[np.float64],
    q: NamedStageQ,
    f: npt.NDArray[np.float64],
    stats: dict[str, Any],
    *,
    deadline_s: float,
    temperature_radius_scale: float,
    couple_trajectory_jacobian: bool,
    settling_margin_eta: float | None,
    energy_guard: bool,
) -> dict[str, Any] | None:
    """氧标量根只作候选初值。联合校正通过后才返回新状态。"""
    rows = oxygen_rows_needing_restart(f, stats.get("predicted_raw"))
    stats["oxygen_restart_used"] = True
    stats.pop("jacobian_cache", None)
    if not rows:
        stats["oxygen_restart"] = {"status": "not_oxygen"}
        return None
    candidate, report = trace_oxygen_candidate(accepted, q, rows, stats)
    stats["oxygen_restart"] = {
        "status": "candidate" if candidate is not None else "declined",
        "rows": rows,
        "report": [
            {key: value for key, value in item.items() if key != "samples"}
            for item in report
        ],
    }
    if candidate is None:
        return None
    if deadline_s <= 0.0:
        return None
    inner = solve_named_controlled(
        nq,
        candidate,
        deadline_s=float(deadline_s),
        temperature_radius_scale=temperature_radius_scale,
        couple_trajectory_jacobian=couple_trajectory_jacobian,
        settling_margin_eta=settling_margin_eta,
        energy_guard=energy_guard,
        oxygen_restart=False,
    )
    for key in ("F_calls", "J_calls", "P_calls"):
        stats[key] = float(stats.get(key, 0.0)) + float(inner["stats"].get(key, 0.0))
    if not bool(inner["passed"]):
        stats["oxygen_restart"]["status"] = "joint_rejected"
        return None
    restarted = np.asarray(inner["x"], dtype=np.float64)
    q_new = outer_P_named(restarted, nq)
    stats["P_calls"] = float(stats.get("P_calls", 0.0)) + 1.0
    f_new = _count_residual(restarted, q_new, stats)
    cert = certify_controlled_state(restarted, f_new, q_new.hydro)
    if not bool(cert["passed"]):
        stats["oxygen_restart"]["status"] = "cold_rejected"
        return None
    stats["oxygen_restart"]["status"] = "accepted"
    return {"x": restarted, "q": q_new, "F": f_new, "cert": cert}


def solve_named_controlled(
    nq: NamedStageQ,
    x0: npt.NDArray[np.float64],
    *,
    deadline_s: float,
    temperature_radius_scale: float = 1.0,
    couple_trajectory_jacobian: bool = True,
    settling_margin_eta: float | None = None,
    energy_guard: bool = False,
    oxygen_restart: bool = True,
) -> dict[str, Any]:
    """有限内层预算。接受时 accepted、q、F、merit 一起换；拒绝时这四者都不换。"""
    t0 = time.monotonic()
    accepted = np.asarray(x0, dtype=np.float64).copy()
    stats: dict[str, Any] = {"F_calls": 0.0, "J_calls": 0.0, "P_calls": 0.0, "outers": 0.0, "rejects": 0.0}
    opened: list[int] = []
    reason = "deadline"
    log_span = 0.5
    ratio_span = 1.0
    energy_guard_W = 0.0
    temp_factor = 1.5
    ftol = 1.0e-4
    max_j = 4
    history: list[dict[str, Any]] = []
    q = outer_P_named(accepted, nq)
    stats["P_calls"] += 1.0
    f = _count_residual(accepted, q, stats)
    best = residual_block_merit(f)
    cert = certify_controlled_state(accepted, f, q.hydro)
    while time.monotonic() - t0 < float(deadline_s):
        if bool(cert["passed"]):
            reason = "block_gates"
            break
        x_try, just = reactivate_zero_unknowns(q, accepted, stats)
        opened.extend(i for i in just if i not in opened)
        if just:
            q_try = outer_P_named(x_try, nq)
            stats["P_calls"] += 1.0
            f_open = _count_residual(x_try, q_try, stats)
            merit_open = residual_block_merit(f_open)
            opened_ok = (
                _state_admissible(x_try, f_open, q_try.hydro)
                and step_improves_progress(f, f_open)
            )
            if opened_ok:
                accepted = x_try
                q = q_try
                f = f_open
                best = merit_open
                cert = certify_controlled_state(accepted, f, q.hydro)
                stats["jacobian_cache"] = None
                if bool(cert["passed"]):
                    reason = "block_gates"
                    break
        remain = float(deadline_s) - (time.monotonic() - t0)
        if remain <= 0.0:
            break
        inner_limit = min(20.0, remain)
        proposed = _e8_log_scaled_correct(
            q, accepted, deadline_s=inner_limit, max_nfev=max(8, max_j), max_j=max_j, ftol=ftol,
            temperature_factor=temp_factor, log_span=log_span, ratio_span=ratio_span, stats=stats,
            temperature_radius_scale=temperature_radius_scale,
            couple_trajectory_jacobian=bool(couple_trajectory_jacobian),
            energy_guard_W=float(energy_guard_W) if energy_guard else 0.0,
            settling_margin_eta=settling_margin_eta,
        )
        stats["outers"] += 1.0
        stats["log_span"] = float(log_span)
        stats["j_budget"] = float(max_j)
        pred = stats.get("predicted_raw")
        q_new = outer_P_named(proposed, nq)
        stats["P_calls"] += 1.0
        f_new = _count_residual(proposed, q_new, stats)
        merit_new = residual_block_merit(f_new)
        cert_new = certify_controlled_state(proposed, f_new, q_new.hydro)
        phi0 = progress_phi(f)
        phi_act = progress_phi(f_new)
        phi_pred = progress_phi(np.asarray(pred, dtype=np.float64)) if pred is not None else phi0
        pred_drop = phi0 - phi_pred
        act_drop = phi0 - phi_act
        rho = float(act_drop / pred_drop) if pred_drop > 1.0e-12 else -1.0
        def margin_ok(state: npt.NDArray[np.float64]) -> bool:
            if settling_margin_eta is None or nq.reference is None:
                return True
            m0 = freeboard_settling_margin_m(accepted, nq.reference)
            m1 = freeboard_settling_margin_m(state, nq.reference)
            eject = np.asarray(q.energy.species.K_eject, dtype=np.float64) * np.asarray(accepted[14:16], dtype=np.float64)
            for i, m_i in enumerate(m0):
                if i >= int(eject.size) or float(eject[i]) < _EJECT_SUBSTANTIAL_KG_S:
                    continue
                if abs(float(m1[i]) - float(m_i)) > float(settling_margin_eta) * abs(float(m_i)) + 1.0e-6:
                    return False
            return True

        admissible = _state_admissible(proposed, f_new, q_new.hydro) and margin_ok(proposed)
        improved = step_improves_progress(f, f_new)
        stats["acceptance_target"] = "gate_ratio_phi"
        stats["material_ratios"] = material_gate_ratios(f)
        stats["progress_phi"] = phi0
        stats["predicted_phi"] = phi_pred
        stats["actual_phi"] = phi_act
        stats["rho"] = rho
        stats["log_span"] = float(log_span)
        on_bound = int(stats.get("active_bound_count", 0)) > 0
        n_t = 2 if int(np.asarray(f_new).size) >= 46 else 1
        take = admissible and (bool(cert_new["passed"]) or (improved and rho >= 0.25 and pred_drop > 1.0e-12))
        f_full = np.asarray(f_new, dtype=np.float64).copy()
        backtrack_alpha = None
        crossed = False
        if not take:
            old_ratios = progress_ratios(f)
            new_ratios = progress_ratios(f_new)
            crossed = any(
                float(before) < 1.0 and float(after) > 1.0
                for before, after in zip(old_ratios, new_ratios)
            )
        if not take and crossed and bool(stats.get("subproblem_solved")):
            partial = progress_backtrack(nq, accepted, proposed, f, stats, extra_ok=margin_ok)
            if partial is not None:
                proposed = partial["x"]
                q_new = partial["q"]
                f_new = partial["F"]
                cert_new = partial["cert"]
                phi_act = float(partial["phi"])
                merit_new = residual_block_merit(f_new)
                backtrack_alpha = float(partial["alpha"])
                take = True
        pred_vec = np.asarray(pred, dtype=np.float64) if pred is not None else None
        phi_pred_acc = phi_pred
        rho_acc = rho
        if backtrack_alpha is not None and pred_vec is not None and pred_vec.shape == np.asarray(f).shape:
            f_lin = np.asarray(f, dtype=np.float64) + float(backtrack_alpha) * (pred_vec - np.asarray(f, dtype=np.float64))
            phi_pred_acc = progress_phi(f_lin)
            drop_acc = phi0 - phi_pred_acc
            rho_acc = float((phi0 - phi_act) / drop_acc) if drop_acc > 1.0e-12 else -1.0
        if energy_guard and pred_vec is not None and pred_vec.shape == f_full.shape:
            n_obs = 2 if int(f_full.size) >= 46 else 1
            err = max(abs(float(f_full[-n_obs + i]) - float(pred_vec[-n_obs + i])) for i in range(n_obs))
            energy_guard_W = min(0.2, max(energy_guard_W, float(err)))
        history.append({
            "log_span": float(log_span),
            "ratio_span": float(ratio_span),
            "phi": phi0,
            "full_step_predicted_phi": phi_pred,
            "full_step_actual_phi": float(stats["actual_phi"]),
            "full_step_rho": rho,
            "predicted_phi": phi_pred_acc,
            "actual_phi": phi_act,
            "rho": rho_acc,
            "accepted": bool(take),
            "passed": bool(cert_new["passed"]),
            "material_ratios": material_gate_ratios(f_new),
            "linear_status": stats.get("linear_status"),
            "subproblem_solved": stats.get("subproblem_solved"),
            "solver_status": stats.get("solver_status"),
            "max_nd_violation": stats.get("max_nd_violation"),
            "projected_gradient_norm": stats.get("projected_gradient_norm"),
            "active_bound_count": stats.get("active_bound_count"),
            "temperature_bound_active": bool(stats.get("temperature_bound_active")),
            "predicted_material_ratios": (
                list(material_gate_ratios(np.asarray(pred, dtype=np.float64)))
                if pred is not None else None
            ),
            "energy_rows_W": [float(v) for v in np.asarray(f_new, dtype=np.float64)[-n_t:]],
            "solid_rows_kg_s": [
                float(np.asarray(f_new, dtype=np.float64)[i])
                for i in (14, 15, 37, 38, 39, 46, 47)
                if i < int(np.asarray(f_new).size)
            ],
            "cyclone_mass_kg": [
                float(np.asarray(proposed, dtype=np.float64)[i])
                for i in (37, 38)
                if i < int(np.asarray(proposed).size)
            ],
            "backtrack_alpha": backtrack_alpha,
        })
        if take:
            accepted = proposed
            q = q_new
            f = f_new
            best = merit_new
            cert = cert_new
            stats["jacobian_cache"] = None
            stats["material_ratios"] = material_gate_ratios(f)
            if backtrack_alpha is None and rho > 0.75 and on_bound:
                log_span = min(_RADIUS_MAX, float(log_span) * 2.0)
                ratio_span = min(1.0, float(ratio_span) * 2.0)
            if bool(cert["passed"]):
                reason = "block_gates"
                break
        else:
            stats["composite_merit"] = merit_new
            if pred_drop <= 1.0e-12:
                stats["direction"] = "no_predicted_decrease"
            elif not admissible:
                stats["direction"] = "inadmissible"
            elif rho < 0.25:
                stats["direction"] = "model_mismatch"
            else:
                stats["direction"] = "progress_rejected"
            log_span *= 0.5
            ratio_span *= 0.5
            noise = max(1.0e-8, 1.0e-6 * max(phi0, 1.0))
            if float(log_span) < 1.0e-3 and float(ratio_span) < 1.0e-3 and (
                pred_drop <= noise or act_drop <= 0.0
            ):
                if not bool(stats.get("subproblem_solved")):
                    reason = "subproblem_not_solved"
                elif stats["direction"] == "no_predicted_decrease" or pred_drop <= noise:
                    reason = "no_descent_direction"
                elif on_bound and act_drop <= 0.0:
                    reason = "constraint_stall"
                else:
                    reason = "model_mismatch"
                if (
                    reason == "no_descent_direction"
                    and oxygen_restart
                    and not stats.get("oxygen_restart_used")
                ):
                    remain_restart = float(deadline_s) - (time.monotonic() - t0)
                    restarted = _finish_oxygen_restart(
                        nq, accepted, q, f, stats,
                        deadline_s=min(20.0, remain_restart),
                        temperature_radius_scale=temperature_radius_scale,
                        couple_trajectory_jacobian=bool(couple_trajectory_jacobian),
                        settling_margin_eta=settling_margin_eta,
                        energy_guard=energy_guard,
                    )
                    if restarted is not None:
                        accepted = restarted["x"]
                        q = restarted["q"]
                        f = restarted["F"]
                        cert = restarted["cert"]
                        best = residual_block_merit(f)
                        reason = "block_gates"
                break
        if stats.get("inner_status") == "deadline" and remain <= inner_limit + 1.0e-9 and time.monotonic() - t0 >= float(deadline_s) - 1.0e-6:
            reason = "deadline"
            break
    stats.pop("jacobian_cache", None)
    stats.pop("predicted_raw", None)
    stats["history"] = history
    elapsed = time.monotonic() - t0
    cert = certify_controlled_state(accepted, f, q.hydro)
    if bool(cert["passed"]):
        reason = "block_gates"
    return {
        "x": accepted,
        "F": f,
        "q_matches_x": True,
        "blocks": cert["blocks"],
        "merit": best,
        "opened_zeros": opened,
        "reason": reason,
        "elapsed_s": elapsed,
        "stats": stats,
        "passed": bool(cert["passed"]),
        "transport_closed": bool(cert["transport_closed"]),
        "temperature_domain": TEMPERATURE_DOMAIN,
        "gate_ratio": cert["gate_ratio"],
    }


def advance_parameter_factor(
    current: float,
    target: float,
    step: float,
    *,
    passed: bool,
) -> tuple[float, float, bool]:
    """参数延拓走一步。未通过时因子不动，步长减半；步长落到 0.01 以下就停止。

    停止时当前因子不是目标因子。不能把未认证的目标写成已到达。
    """
    span = abs(float(target) - float(current))
    if span <= 1.0e-12:
        return float(current), float(step), False
    if passed:
        direction = 1.0 if float(target) > float(current) else -1.0
        moved = float(current) + direction * min(float(step), span)
        grown = min(abs(float(target) - 1.0), float(step) * 2.0)
        return moved, max(grown, 1.0e-6), False
    shrunk = float(step) * 0.5
    return float(current), shrunk, shrunk < 0.01
