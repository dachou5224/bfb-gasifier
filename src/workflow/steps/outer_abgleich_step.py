"""流程图：Check2 外层 Abgleich（Vorabrechnung ↔ 单元模型对齐）。

薄封装 ``run_global_nr_outer_abgleich``，便于流程图审计与导入路径固定。
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from src.solvers.outer_loop import OuterLoopResult, run_global_nr_outer_abgleich

if TYPE_CHECKING:
    from src.core.cell import Cell
    from src.solvers.convergence import OuterConvergence, OuterMeritEarlyStop


def run_outer_abgleich_for_global_nr(
    *,
    cells: list["Cell"],
    outer_max: int,
    total_inner_budget: int,
    tol: float,
    jacobian_mode: str,
    jacobian_lag: int,
    tol_rms: float,
    refresh_fn: Callable[[bool], None],
    snapshot_signature_fn: Callable[[], str],
    solve_inner_fn: Callable[[int, float, str, int], dict],
    pre_inner_fn: Callable[[], dict[str, bool | float] | None] | None = None,
    post_inner_holdup_fn: Callable[[], bool] | None = None,
    pre_inner_on_continuation: bool = False,
    outer_convergence: "OuterConvergence | None" = None,
    gas_metrics_fn: Callable[[], dict[str, float]] | None = None,
    check2_bridge_fn: Callable[[], dict[str, bool | float]] | None = None,
    strict_check1_before_refresh: bool = False,
    check1_energy_regress_refresh_ratio: float = 0.0,
    outer_context_fn: Callable[..., None] | None = None,
    merit_early_stop: "OuterMeritEarlyStop | None" = None,
    hot_restart_after_y_co: bool = False,
    hot_restart_y_co_threshold: float = 0.14,
    hot_restart_bonus_inner_budget: int = 0,
    hot_restart_bonus_outer_slots: int = 0,
    y_co_top_fn: Callable[[], float] | None = None,
    hot_restart_clear_runtime_fn: Callable[[], None] | None = None,
) -> OuterLoopResult:
    """外层对齐循环：刷新 Vorabrechnung → 相分配预投影 → 内层 NR → Check2 对齐。"""
    return run_global_nr_outer_abgleich(
        cells=cells,
        outer_max=outer_max,
        total_inner_budget=total_inner_budget,
        tol=tol,
        jacobian_mode=jacobian_mode,
        jacobian_lag=jacobian_lag,
        tol_rms=tol_rms,
        refresh_fn=refresh_fn,
        snapshot_signature_fn=snapshot_signature_fn,
        solve_inner_fn=solve_inner_fn,
        pre_inner_fn=pre_inner_fn,
        post_inner_holdup_fn=post_inner_holdup_fn,
        pre_inner_on_continuation=pre_inner_on_continuation,
        outer_convergence=outer_convergence,
        gas_metrics_fn=gas_metrics_fn,
        check2_bridge_fn=check2_bridge_fn,
        strict_check1_before_refresh=strict_check1_before_refresh,
        check1_energy_regress_refresh_ratio=check1_energy_regress_refresh_ratio,
        outer_context_fn=outer_context_fn,
        merit_early_stop=merit_early_stop,
        hot_restart_after_y_co=hot_restart_after_y_co,
        hot_restart_y_co_threshold=hot_restart_y_co_threshold,
        hot_restart_bonus_inner_budget=hot_restart_bonus_inner_budget,
        hot_restart_bonus_outer_slots=hot_restart_bonus_outer_slots,
        y_co_top_fn=y_co_top_fn,
        hot_restart_clear_runtime_fn=hot_restart_clear_runtime_fn,
    )
