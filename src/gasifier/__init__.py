"""与 Hamel 流程图分层的命名空间（domain / models / correlations / workflow / solvers）。

- **主执行与验证路径（canonical）**：``src.core`` + ``src.solvers`` + ``src.workflow``。
- **本包内** ``workflow``、``solvers`` 子目录为上述模块的 **重导出别名**。
- **``models/*Model`` / ``domain.GlobalState``** 为 Phase0 演示（``DeprecationWarning``），
  不得与 ``src.core.connectivity`` 等 production 实现并行演进；活跃子包见 ``models.correlations``。
"""

__all__ = []
