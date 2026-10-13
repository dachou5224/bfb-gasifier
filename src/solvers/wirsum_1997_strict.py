"""Wirsum (1997) §2.3.3 最小阻尼 Newton。

只做：

1. 组装 ``F(x)`` 与 ``J(x)``
2. 解 ``J Δx = -F``（Wirsum Eq. 2.21）
3. 试探 ``α = 1, 1/2, 1/4, …``
4. 仅当 merit 严格下降时接受（原文 Eq. 2.23 为 ``||F||₂``）

不含温度 fence、fastox、syngas floor、split/energy merit、
refresh-aware accept 或 residual 内 clip。

生产 ``solve_global_nr`` 仍是工程路径；本模块是 reduced / A/B 基准。

Source: Wirsum (1997) §2.3.3 Eq. 2.21–2.23; docs/wirsum_1997_nr_solver_notes.md
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

import numpy as np
import numpy.typing as npt

from src.core.solid_balance_reduced import jacobian_diagnostics

NR_SOLVER_MODE_ENGINEERING: str = "engineering"
NR_SOLVER_MODE_WIRSUM_1997_STRICT: str = "wirsum_1997_strict"
_NR_SOLVER_MODES = frozenset(
    {NR_SOLVER_MODE_ENGINEERING, NR_SOLVER_MODE_WIRSUM_1997_STRICT}
)

MeritKind = Literal["raw_l2", "scaled_rms"]
JacobianRefresh = Literal["fresh", "frozen"]

ResidualFn = Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]]
JacobianFn = Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]]
# 返回 Δx、用于停机的秩、条件数。缺省走原始 J 的 SVD 秩。
LinearSolveFn = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    tuple[npt.NDArray[np.float64], int, float],
]


def resolve_nr_solver_mode(configured: str) -> str:
    """校验 ``nr_solver_mode``。"""
    name = str(configured).strip().lower()
    if name not in _NR_SOLVER_MODES:
        raise ValueError(
            f"Unknown nr_solver_mode={configured!r}; "
            f"expected {NR_SOLVER_MODE_ENGINEERING!r} or {NR_SOLVER_MODE_WIRSUM_1997_STRICT!r}"
        )
    return name


def apply_wirsum_1997_strict_config(cfg: Any) -> None:
    """关掉可由 config 开关的工程桥，并禁止 K=0 库存锚。

    不改 ``solve_global_nr`` 内部的缩放 RMS merit；Eq. 2.23 的 ``||F||₂``
    只由 ``wirsum_1997_strict_newton`` 实现。
    不改 ``thesis_vorab_sources_single_shot`` 与 ``vorab_outer_source_policy``。
    """
    cfg.nr_solver_mode = NR_SOLVER_MODE_WIRSUM_1997_STRICT
    cfg.holdup_k0_anchor = "none"
    cfg.solid_transport_inventory_mode = "species_own_inventory"
    cfg.solid_transport_bottom_no_penetration = True
    cfg.nr_step_model_thesis = "newton"
    cfg.nr_prefer_full_step_thesis = True
    cfg.nr_enforce_gas_phase_split_convergence_thesis = False
    cfg.nr_line_search_gas_phase_split_merit_thesis = False
    cfg.nr_line_search_energy_merit_thesis = False
    cfg.nr_refresh_hydrodynamics_on_accepted_step_thesis = False
    cfg.nr_line_search_fast_oxidation_projection_thesis = False
    cfg.nr_line_search_per_phase_fastox_thesis = False
    cfg.nr_clip_same_phase_oxidizer_fuel_step_thesis = False
    cfg.nr_clip_allow_oxidizer_into_fuel_thesis = False
    cfg.nr_line_search_syngas_collapse_guard_thesis = False
    cfg.nr_line_search_syngas_floor_ratio_thesis = 0.0
    cfg.nr_line_search_refresh_aware_accept_thesis = False
    cfg.nr_line_search_split_constrained_step_thesis = False
    cfg.nr_clip_solid_ftb_relative_floor_thesis = False
    cfg.nr_nonmonotone_enabled_thesis = False


def eval_merit(
    F: npt.NDArray[np.float64],
    *,
    kind: MeritKind,
    scale: npt.NDArray[np.float64] | None,
) -> float:
    """线搜索标量。``raw_l2`` 是 Wirsum Eq. 2.23；``scaled_rms`` 是工程对照。"""
    F = np.asarray(F, dtype=np.float64)
    if kind == "raw_l2":
        return float(np.linalg.norm(F))
    if kind == "scaled_rms":
        if scale is None:
            raise ValueError("scaled_rms 需要 scale [与 F 同形]")
        scale_arr = np.asarray(scale, dtype=np.float64)
        assert scale_arr.shape == F.shape, "scale 与 F 维度不一致"
        F_hat = F / np.maximum(scale_arr, 1.0e-30)
        return float(np.sqrt(np.mean(F_hat**2))) if F_hat.size else 0.0
    raise ValueError(f"未知 merit kind={kind!r}")


@dataclass
class WirsumStrictSolveResult:
    """一次 strict Newton 的轨迹与终止原因。"""

    x: npt.NDArray[np.float64]
    F: npt.NDArray[np.float64]
    converged: bool
    iterations: int
    accepted_alphas: list[float] = field(default_factory=list)
    residual_l2_history: list[float] = field(default_factory=list)
    merit_history: list[float] = field(default_factory=list)
    n_jac_builds: int = 0
    n_residual_evals: int = 0
    cond_last: float = float("inf")
    n_zero_rows_last: int = 0
    n_zero_cols_last: int = 0
    rank_last: int = 0
    failure_reason: str = ""
    merit_kind: str = "raw_l2"
    jacobian_refresh: str = "fresh"

    def as_dict(self) -> dict[str, Any]:
        return {
            "converged": bool(self.converged),
            "iterations": int(self.iterations),
            "accepted_alphas": [float(a) for a in self.accepted_alphas],
            "residual_l2_history": [float(v) for v in self.residual_l2_history],
            "merit_history": [float(v) for v in self.merit_history],
            "n_jac_builds": int(self.n_jac_builds),
            "n_residual_evals": int(self.n_residual_evals),
            "cond_last": float(self.cond_last),
            "n_zero_rows_last": int(self.n_zero_rows_last),
            "n_zero_cols_last": int(self.n_zero_cols_last),
            "rank_last": int(self.rank_last),
            "failure_reason": str(self.failure_reason),
            "merit_kind": str(self.merit_kind),
            "jacobian_refresh": str(self.jacobian_refresh),
            "x": np.asarray(self.x, dtype=np.float64).tolist(),
            "F": np.asarray(self.F, dtype=np.float64).tolist(),
            "residual_l2_final": float(np.linalg.norm(self.F)),
        }


def wirsum_1997_strict_newton(
    residual_fn: ResidualFn,
    jacobian_fn: JacobianFn,
    x0: npt.ArrayLike,
    *,
    merit_kind: MeritKind = "raw_l2",
    jacobian_refresh: JacobianRefresh = "fresh",
    scale: npt.ArrayLike | None = None,
    max_iter: int = 25,
    lambda_min: float = 1.0 / 4096.0,
    tol_merit: float = 1.0e-10,
    tol_step: float = 1.0e-12,
    trial_in_domain: Callable[[npt.NDArray[np.float64]], bool] | None = None,
    linear_solve_fn: LinearSolveFn | None = None,
) -> WirsumStrictSolveResult:
    """对任意 ``F(x)`` 跑 Wirsum 二分阻尼 Newton。

    ``scale`` 仅 ``merit_kind='scaled_rms'`` 使用，单位与 ``F`` 相同。
    ``trial_in_domain`` 为 None 时保持原行为；给定时越界试探缩短 α，
    且不调用 ``residual_fn``。线搜索 merit 仍是 ``raw_l2`` 或具名 scaled_rms。
    ``linear_solve_fn`` 只替换 ``J Δx = -F`` 与秩判定；不改 merit。
    """
    if merit_kind not in {"raw_l2", "scaled_rms"}:
        raise ValueError(f"未知 merit_kind={merit_kind!r}")
    if jacobian_refresh not in {"fresh", "frozen"}:
        raise ValueError(f"未知 jacobian_refresh={jacobian_refresh!r}")

    x = np.asarray(x0, dtype=np.float64).copy()
    scale_arr = None if scale is None else np.asarray(scale, dtype=np.float64)
    if scale_arr is not None:
        assert scale_arr.shape == x.shape, "scale 与 x 维度不一致"

    F = np.asarray(residual_fn(x), dtype=np.float64)
    n_res = 1
    assert F.shape == x.shape, "F(x) 必须与 x 同形"
    merit = eval_merit(F, kind=merit_kind, scale=scale_arr)
    l2_hist = [float(np.linalg.norm(F))]
    merit_hist = [float(merit)]
    alphas: list[float] = []
    n_jac = 0
    J_frozen: npt.NDArray[np.float64] | None = None
    cond_last = float("inf")
    n_zero_rows = 0
    n_zero_cols = 0
    rank_last = 0
    failure = ""
    converged = bool(merit < float(tol_merit))

    if jacobian_refresh == "frozen":
        J_frozen = np.asarray(jacobian_fn(x), dtype=np.float64)
        n_jac += 1

    it = 0
    if not converged:
        for it in range(1, int(max_iter) + 1):
            if jacobian_refresh == "fresh":
                J = np.asarray(jacobian_fn(x), dtype=np.float64)
                n_jac += 1
            else:
                assert J_frozen is not None
                J = J_frozen
            assert J.shape == (x.size, x.size), "J 必须是 n×n"
            diag = jacobian_diagnostics(J)
            n_zero_rows = int(diag.n_zero_rows)
            n_zero_cols = int(diag.n_zero_cols)
            if linear_solve_fn is None:
                cond_last = float(diag.cond)
                rank_last = int(diag.rank)
                if diag.rank < diag.n:
                    failure = "singular_jacobian"
                    break
                try:
                    dx = np.linalg.solve(J, -F)
                except np.linalg.LinAlgError:
                    failure = "singular_jacobian"
                    break
            else:
                try:
                    dx, rank_last, cond_last = linear_solve_fn(J, F, x)
                    rank_last = int(rank_last)
                    cond_last = float(cond_last)
                except np.linalg.LinAlgError:
                    failure = "singular_jacobian"
                    rank_last = int(diag.rank)
                    cond_last = float(diag.cond)
                    break
                if rank_last < int(J.shape[0]):
                    failure = "singular_jacobian"
                    break

            accepted = False
            alpha = 1.0
            F_trial = F
            x_trial = x
            merit_trial = merit
            while alpha >= float(lambda_min) * 0.5:
                if alpha < float(lambda_min):
                    break
                x_trial = x + alpha * dx
                if trial_in_domain is not None and not trial_in_domain(x_trial):
                    alpha *= 0.5
                    continue
                F_trial = np.asarray(residual_fn(x_trial), dtype=np.float64)
                n_res += 1
                merit_trial = eval_merit(F_trial, kind=merit_kind, scale=scale_arr)
                if np.isfinite(merit_trial) and merit_trial < merit:
                    accepted = True
                    break
                alpha *= 0.5

            if not accepted:
                failure = "line_search"
                break

            x = x_trial
            F = F_trial
            merit = float(merit_trial)
            alphas.append(float(alpha))
            l2_hist.append(float(np.linalg.norm(F)))
            merit_hist.append(float(merit))
            if merit < float(tol_merit) or float(np.linalg.norm(alpha * dx)) < float(tol_step):
                converged = True
                failure = ""
                break
        else:
            failure = "max_iter"

    return WirsumStrictSolveResult(
        x=x,
        F=F,
        converged=converged,
        iterations=int(it if not (converged and it == 0) else 0),
        accepted_alphas=alphas,
        residual_l2_history=l2_hist,
        merit_history=merit_hist,
        n_jac_builds=n_jac,
        n_residual_evals=n_res,
        cond_last=cond_last,
        n_zero_rows_last=n_zero_rows,
        n_zero_cols_last=n_zero_cols,
        rank_last=rank_last,
        failure_reason="" if converged else failure,
        merit_kind=str(merit_kind),
        jacobian_refresh=str(jacobian_refresh),
    )


def run_wirsum_strict_ab(
    residual_fn: ResidualFn,
    jacobian_fn: JacobianFn,
    x0: npt.ArrayLike,
    *,
    scale: npt.ArrayLike | None = None,
    max_iter: int = 40,
    lambda_min: float = 1.0 / 4096.0,
    tol_merit: float = 1.0e-10,
) -> dict[str, dict[str, Any]]:
    """四组 A/B：raw/scaled × fresh/frozen。"""
    variants: tuple[tuple[str, MeritKind, JacobianRefresh], ...] = (
        ("A_raw_fresh", "raw_l2", "fresh"),
        ("B_scaled_fresh", "scaled_rms", "fresh"),
        ("C_raw_frozen", "raw_l2", "frozen"),
        ("D_scaled_frozen", "scaled_rms", "frozen"),
    )
    out: dict[str, dict[str, Any]] = {}
    for name, merit_kind, refresh in variants:
        if merit_kind == "scaled_rms" and scale is None:
            raise ValueError("A/B scaled 分支需要 scale")
        result = wirsum_1997_strict_newton(
            residual_fn,
            jacobian_fn,
            x0,
            merit_kind=merit_kind,
            jacobian_refresh=refresh,
            scale=scale,
            max_iter=max_iter,
            lambda_min=lambda_min,
            tol_merit=tol_merit,
        )
        out[name] = result.as_dict()
    return out


def _slim_ab_result(result: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "converged",
        "iterations",
        "accepted_alphas",
        "residual_l2_history",
        "merit_history",
        "n_jac_builds",
        "cond_last",
        "n_zero_rows_last",
        "rank_last",
        "failure_reason",
        "merit_kind",
        "jacobian_refresh",
        "residual_l2_final",
    )
    slim = {key: result[key] for key in keys}
    cond = slim["cond_last"]
    if isinstance(cond, float) and not np.isfinite(cond):
        slim["cond_last"] = None
    return slim


def build_reduced_holdup_ab_payload() -> dict[str, Any]:
    """P2.1 三格 holdup A/B：线性 / 非线性 MM / 行加权 / K=0。"""
    from src.core.solid_balance_reduced import (
        N_CELLS,
        canonical_open_chain_K,
        holdup_jacobian,
        holdup_residual,
    )

    K_up, K_down = canonical_open_chain_K()
    m_feed = np.zeros(N_CELLS, dtype=np.float64)
    m_feed[0] = 1.0  # [kg/s]
    r_dot = np.zeros(N_CELLS, dtype=np.float64)
    k_sat = np.array([0.50, 0.20, 0.10], dtype=np.float64)  # [kg/s]
    M_ref = np.ones(N_CELLS, dtype=np.float64)  # [kg]
    x0 = np.array([0.05, 0.05, 0.05], dtype=np.float64)  # [kg]
    scale = np.ones(N_CELLS, dtype=np.float64)  # [kg/s]
    weights = np.array([1.0, 1.0e3, 1.0e-3], dtype=np.float64)

    def linear_residual(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return holdup_residual(
            M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
        )

    def linear_jacobian(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return holdup_jacobian(
            M, m_feed=m_feed, K_up=K_up, K_down=K_down, r_dot=r_dot, k0_anchor="none"
        )

    def mm_residual(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return linear_residual(M) - k_sat * M / (1.0 + M / M_ref)

    def mm_jacobian(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        J = linear_jacobian(M).copy()
        J[np.diag_indices(N_CELLS)] -= k_sat * M_ref / (M_ref + M) ** 2
        return J

    def weighted_residual(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return weights * mm_residual(M)

    def weighted_jacobian(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return weights[:, None] * mm_jacobian(M)

    zero_K = np.zeros(N_CELLS, dtype=np.float64)

    def k0_residual(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return holdup_residual(
            M, m_feed=m_feed, K_up=zero_K, K_down=zero_K, r_dot=r_dot, k0_anchor="none"
        )

    def k0_jacobian(M: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return holdup_jacobian(
            M, m_feed=m_feed, K_up=zero_K, K_down=zero_K, r_dot=r_dot, k0_anchor="none"
        )

    k0 = wirsum_1997_strict_newton(k0_residual, k0_jacobian, np.ones(N_CELLS, dtype=np.float64))
    return {
        "problem": {
            "name": "3-cell holdup reduced",
            "x0_kg": x0.tolist(),
            "feed_kg_s": m_feed.tolist(),
            "k_sat_kg_s": k_sat.tolist(),
            "row_weights": weights.tolist(),
        },
        "linear_K_positive": {
            key: _slim_ab_result(val)
            for key, val in run_wirsum_strict_ab(
                linear_residual, linear_jacobian, x0, scale=scale, max_iter=20
            ).items()
        },
        "nonlinear_mm": {
            key: _slim_ab_result(val)
            for key, val in run_wirsum_strict_ab(
                mm_residual, mm_jacobian, x0, scale=scale, max_iter=40
            ).items()
        },
        "nonlinear_mm_row_weighted": {
            key: _slim_ab_result(val)
            for key, val in run_wirsum_strict_ab(
                weighted_residual, weighted_jacobian, x0, scale=weights, max_iter=40
            ).items()
        },
        "K_zero_no_constitutive_outflow": _slim_ab_result(k0.as_dict()),
    }
