# 求解器审计与稳定性报告 (Solver Audit & Stability Report)

## 1. 科学审计重构 (Scientific Audit Refactoring)
日期：2026-03-21

为了确保物理模型的严谨性，已对全代码库进行了以下重构：
- **单位制统一**：强制执行 SI 单位制 (Pa, K, m, mol, kg, s)。修复了 `gibbs_minimizer.py` 等文件中硬编码的 P0 和 Rg。
- **动力学工厂化**：所有反应速率 (R1–R11) 现在必须调用 `src/kinetics/arrhenius.py` 中的工厂函数 (`k_hobbs`, `k_standard`, `k_jensen_r7`)。严禁在反应逻辑中直接使用 `np.exp`，以确保数值溢出保护（`np.clip`）的全局一致性。
- **常数管理**：所有物理常数统一从 `src/core/constants.py` 导入。

## 2. 求解器稳定性检查 (Solver Stability Check)

### Gauss-Seidel (GS) 扫描求解器
- **状态**：**生产级稳定 (Production Ready)**。
- **表现**：在 Table 2 LU 工况下收敛速度快（约 4s），对 1200K 的粗略初值不敏感。
- **收敛特征**：利用了流体流动的方向性（从底至顶），通过 `least_squares` 逐格解决非线性方程，鲁棒性极高。

### 全局 Newton-Raphson (NR) 求解器
- **状态**：**调试中 / 辅助验证 (Debugging / Warm-start only)**。
- **已知问题**：
    - **初值敏感性**：直接从 1200K 开始全局迭代会导致残差在早期发散或陷入高值平台（如 1607K 的非物理状态）。
    - **Jacobian 性能**：有限差分构造稀疏矩阵在大规模 cell 时较慢。
- **已实施修复**：
    - 修正了有限差分步长 $h = \epsilon (|x| + 10^{-12})$，解决了微量组分导致的导数失效问题。
    - 重构了阻尼线搜索 (Damped Newton)，确保残差范数单调递减。

## 3. 收敛策略建议 (Convergence Strategy)
为了平衡鲁棒性与高精度，推荐执行 **"GS-then-NR"** 策略：
1. 使用 GS 扫描求解器快速获得一个物理合理的粗略剖面。
2. 将 GS 的结果作为初值 (Warm Start) 传递给全局 NR，利用其二阶收敛特性进行最后的高精度精修。
