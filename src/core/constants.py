"""全局物理常数（集中管理，禁止在各函数内硬编码）。

Source: docs/CLAUDE.md §全局规则 3
"""

Rg: float = 8.314          # [J/(mol·K)] 通用气体常数
g: float = 9.81            # [m/s²] 重力加速度
P0: float = 101_325.0      # [Pa] 标准大气压（加压修正基准）
BAR_PA: float = 100_000.0  # [Pa] 1 bar
MMHG_TO_PA: float = 133.32236842105263  # [Pa/mmHg] Antoine 方程常用换算
P0_HAMEL: float = 101_300.0  # [Pa] Hamel (1999) / DIN-era reference pressure for pressure corrections
n_b: float = 2.7           # [-] 气体交换因子（Hilligardt 常数）
n_B_coalescence: float = 5.0  # [-] Briens et al. (1988); simultaneous coalescing bubbles in Eq.3.41
T_REF: float = 298.15      # [K] 热力学参考温度
