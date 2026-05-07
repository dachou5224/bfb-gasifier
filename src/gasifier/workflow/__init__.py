"""流程编排层：与 ``docs/gasifier_model_flowchart_trilingual.mmd`` 主路径一致。

实现与枚举定义在 ``src.workflow``；本包仅作 **命名空间别名**，避免在
``src/gasifier/`` 下缺失「workflow」目录时与文档/流程图分层不一致。
"""

from src.workflow.simulation_runner import (
    HamelFlowStep,
    build_reactor_and_run,
    describe_hamel_flowchart_mapping,
    hamel_global_nr_pipeline_order,
    run_global_nr_workflow,
)

__all__ = [
    "HamelFlowStep",
    "build_reactor_and_run",
    "describe_hamel_flowchart_mapping",
    "hamel_global_nr_pipeline_order",
    "run_global_nr_workflow",
]
