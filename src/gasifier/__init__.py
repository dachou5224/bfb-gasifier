"""与 Hamel 流程图分层的命名空间（domain / models / correlations / workflow / solvers）。

- **主执行与验证路径（canonical）**：``src.core`` + ``src.solvers`` + ``src.workflow``。
- **本包内** ``workflow``、``solvers`` 子目录为上述模块的 **重导出别名**，与
  ``docs/gasifier_model_flowchart_trilingual.mmd`` 及
  ``describe_hamel_flowchart_mapping()`` 对照时目录结构一致，**不重复实现**。
"""

__all__ = []
