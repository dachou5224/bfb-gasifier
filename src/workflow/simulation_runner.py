"""薄编排：Config → Vorabrechnung → Zellenmodell / global NR → outer Abgleich → Result。

与流程图对应关系见 ``describe_hamel_flowchart_mapping`` 及
``docs/gasifier_model_flowchart_trilingual.mmd``。
"""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from src.core.reactor import Reactor, ReactorConfig


class HamelFlowStep(str, Enum):
    """流程图节点 → 代码入口（命名固定，便于审计）。"""

    CONFIG = "config"  # ReactorConfig / Module 1
    INIT = "init"  # Reactor.__init__ → _build_cells, 赋初值在 _solve_global_nr 内
    VORABRECHNUNG = "vorabrechnung"  # vorabrechnung.estimate_axial_T_profile, generate_initial_x0, refresh
    CELL_MODEL_NR = "cell_model_nr"  # global_nr_solver.solve_global_nr（全域联立）
    OUTER_ABGLEICH = "outer_abgleich"  # solvers.outer_loop.run_global_nr_outer_abgleich
    OUTPUT = "output"  # result_builder.finalize_global_nr_result


def describe_hamel_flowchart_mapping() -> dict[str, str]:
    """流程图节点到 Python 模块的静态映射（无源码论文反演时的主线索引）。"""
    return {
        "Konfiguration / Read Input": "src.core.reactor.ReactorConfig",
        "Initialisierung / Initialize": "src.core.reactor.Reactor._build_cells",
        "Belegung der Startwerte": "src.solvers.vorabrechnung.generate_initial_x0",
        "Vorabrechnung (Strömung, Trocknung/Pyrolyse)": "src.solvers.vorabrechnung.refresh_vorabrechnung_for_cells",
        "Zellenmodell / Bilanzgleichungen / Newton-Raphson": "src.solvers.global_nr_solver.solve_global_nr",
        "Abgleich Vorabrechnung ↔ Zellenmodell (Check2)": "src.solvers.outer_loop.run_global_nr_outer_abgleich",
        "Ausgabe": "src.solvers.result_builder.finalize_global_nr_result",
        # 流程图编排层（与 mmd 顺序一致的可导入步骤；算法仍在 solvers/core）
        "Workflow · INIT+PRECALC step": "src.workflow.steps.init_precalc_step.run_init_and_precalc_for_global_nr",
        "Workflow · NR inner (Check1) builder": "src.workflow.steps.nr_inner_step.build_global_nr_inner_solve_fn",
        "Workflow · OUTER_ABGLEICH (Check2) step": "src.workflow.steps.outer_abgleich_step.run_outer_abgleich_for_global_nr",
        "Freeboard explicit graph (closure → cells)": "src.core.freeboard_bridge.refresh_explicit_freeboard_transport_from_closure",
        "Side-block explicit graph (top chain → cyclone/return)": "src.core.side_block_bridge.initialize_explicit_side_block_states",
        "NR connectivity graph neighborhood": "src.core.connectivity_graph.affected_nr_residual_cells",
        "Cell residual pipeline": "src.core.cell_residual_pipeline.execute_cell_residual_pipeline",
        "Cell · mole fractions / concentrations / enthalpy helpers": "src.core.cell_thermo_helpers.gas_mole_fractions",
    }


def hamel_global_nr_pipeline_order() -> list[str]:
    """主求解链路的静态顺序（用于契约测试与审计脚本）。"""
    return [
        "src.workflow.steps.init_precalc_step.run_init_and_precalc_for_global_nr",
        "src.workflow.steps.nr_inner_step.build_global_nr_inner_solve_fn",
        "src.solvers.outer_loop.run_global_nr_outer_abgleich",
        "src.solvers.result_builder.finalize_global_nr_result",
    ]


def run_global_nr_workflow(
    reactor: "Reactor",
    *,
    max_global_iter: int = 20,
    tol_global: float = 1e-4,
    verbose: bool = False,
    nr_init_strategy: str | None = None,
    nr_jacobian_strategy: str | None = None,
    nr_linear_solver_backend: str | None = None,
    nr_jacobian_lag_steps: int | None = None,
) -> dict[str, Any]:
    """执行与 Hamel 主线一致的 global NR 求解（薄封装，内部即 ``Reactor.solve`` NR 路径）。

    步骤顺序：:py:data:`HamelFlowStep` 与 ``describe_hamel_flowchart_mapping()`` 一致。
    """
    # 单一入口，避免业务代码散落直接调 _solve_global_nr
    return reactor.solve(
        max_global_iter=max_global_iter,
        tol_global=tol_global,
        solver="global_nr",
        verbose=verbose,
        nr_init_strategy=nr_init_strategy,
        nr_jacobian_strategy=nr_jacobian_strategy,
        nr_linear_solver_backend=nr_linear_solver_backend,
        nr_jacobian_lag_steps=nr_jacobian_lag_steps,
    )


def build_reactor_and_run(
    config: "ReactorConfig",
    **solve_kwargs: Any,
) -> dict[str, Any]:
    """由配置构建 Reactor 并运行 global NR workflow。"""
    from src.core.reactor import Reactor

    return run_global_nr_workflow(Reactor(config), **solve_kwargs)
