"""全局物理常数（集中管理，禁止在各函数内硬编码）。

Source: docs/CLAUDE.md §全局规则 3
"""

Rg: float = 8.314          # [J/(mol·K)] 通用气体常数
g: float = 9.81            # [m/s²] 重力加速度
P0: float = 101_325.0      # [Pa] 标准大气压（加压修正基准）
n_b: float = 2.7           # [-] 气体交换因子（Hilligardt 常数）
T_REF: float = 298.15      # [K] 热力学参考温度
