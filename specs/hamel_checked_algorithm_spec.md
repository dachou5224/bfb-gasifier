# Hamel 原版算法对照规格（已核验版）

说明：本文件仅保留“已通过原版论文 PDF 对照核验”的算法条目。  
核验源：`docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf`。

## 1. 程序结构与求解

- `Bild 2.2`：程序结构包含配置、初始化、Vorabrechnung、Zellenmodell、Newton-Raphson、外层 Abgleich。  
- `Eq.2.9`：非线性方程组采用 Newton-Raphson，原文同时说明可用“gedämpfte”变体。  
- 原文还明确：函数矩阵为块三对角，非邻接流引入 Nebenelemente（如旋风回流）。

实现对照：
- `src/solvers/global_nr_solver.py`
- `src/solvers/vorabrechnung.py`
- `src/core/reactor.py`

## 2. 守恒方程

- `Eq.2.1` / `Eq.2.2`：两相气体守恒。  
- `Eq.2.7`：能量守恒。  
- 文中明确焓流包含标准生成焓（与反应热自动并入守恒一致）。

实现对照：
- `src/core/cell_balances.py`
- `src/core/connectivity.py`

## 3. 水力学主链

- `Eq.3.11`、`Eq.3.12`：`u_d` 相关闭合。  
- `Eq.3.41`、`Eq.3.42`、`Eq.3.43`、`Eq.3.44`：加压气泡链。  
- `Eq.3.50`：气泡-悬浮相传质系数主式。  

实现对照：
- `src/physics/bubble_dynamics.py`
- `src/physics/mass_transfer.py`
- `src/core/cell_hydrodynamics.py`

## 4. 干燥与热解（Chapter 4）

- `Eq.4.2`、`Eq.4.4`、`Eq.4.6`、`Eq.4.9`：干燥边界与换热。  
- `Eq.4.10`~`Eq.4.12`：DAEM 与径向积分。  

实现对照：
- `src/thermal/drying.py`
- `src/thermal/devolatilization.py`

## 5. 动力学与平衡

- `Eq.5.36`：CO 氧化速率表达。  
- `Eq.5.45`、`Eq.5.46`、`Eq.5.47`：WGSR 驱动力与平衡修正。  
- `Eq.5.51`：甲烷重整。  
- `Eq.5.59`：焦油氧化主式。  
- `Table 5.4`：焦油裂解/催化参数来源。  

实现对照：
- `src/kinetics/gas_reactions.py`
- `src/kinetics/tar_reactions.py`
- `src/thermodynamics/equilibrium.py`

## 6. A1（气相 Gibbs 平衡）

- 附录 A1 可定位 `A.27`~`A.34` 等 Newton 求解链。  
- 作为“同构级”验收，需逐式映射实现，不得只做功能等价声明。

实现对照：
- `src/thermodynamics/gibbs_minimizer.py`
- `src/thermodynamics/gibbs_hamel_reduced.py`

## 7. 当前明确不一致项（保留为整改项）

1. `R10 (Eq.5.59)`：原文参数单位记号为 Pa 指数项；当前实现存在 `Pa -> bar` 转换与 `k10` 工程缩放。  
2. 默认水力学配置并非论文加压链主配置（需 `thesis_mode` 才切换）。  
3. Jacobian 目前以 FD 近似为主，尚未完成“逐项同构证明”。

这些差异在未整改前，不得宣称“与原版论文完全一致”。
