"""流程图：Check2 外层 Abgleich（Vorabrechnung ↔ 单元模型对齐）。

薄封装 ``run_global_nr_outer_abgleich``，便于流程图审计与导入路径固定。
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from src.solvers.outer_loop import OuterLoopResult, run_global_nr_outer_abgleich

if TYPE_CHECKING:
    from src.core.cell import Cell
    from src.solvers.convergence import OuterConvergence


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
    outer_convergence: "OuterConvergence | None" = None,
) -> OuterLoopResult:
    """外层对齐循环：刷新 Vorabrechnung → 内层 NR → 检查轴向温度对齐。"""
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
        outer_convergence=outer_convergence,
    )
