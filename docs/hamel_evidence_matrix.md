# Hamel 原版算法证据矩阵（清理后）

口径：只保留“原版论文 PDF 可定位 + 源码可定位”的证据对照。  
主源：`docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf`。

## 1) 求解结构

| 项 | 论文锚点 | 实现锚点 | 判定 |
|---|---|---|---|
| 全局 NR 主线 | Kapitel 2.3, Eq.2.9, Bild 2.2 | `src/solvers/global_nr_solver.py` | 部分一致（结构一致，离散未证同构） |
| 预算-求解层级 | Bild 2.2 | `src/solvers/vorabrechnung.py`, `src/core/reactor.py` | 部分一致（存在工程分支） |
| 连接矩阵与侧块 | Kapitel 2（Verbindungsmatrix/Nebenelemente） | `src/core/connectivity.py`, `src/solvers/structured_jacobian.py` | 部分一致（模板化实现，未逐式同构） |

## 2) 守恒方程

| 项 | 论文锚点 | 实现锚点 | 判定 |
|---|---|---|---|
| 两相气相守恒 | Eq.2.1, Eq.2.2 | `src/core/cell_balances.py` | 一致（方程结构） |
| 能量守恒 | Eq.2.7 | `src/core/cell_balances.py` | 部分一致（工程闭合项待严格同构） |

## 3) 水力学与相间交换

| 项 | 论文锚点 | 实现锚点 | 判定 |
|---|---|---|---|
| `u_d` 与可见气泡率链 | Eq.3.11, Eq.3.12 | `src/core/cell_hydrodynamics.py` | 部分一致（默认路径非论文主链） |
| 加压气泡链 | Eq.3.41, Eq.3.42, Eq.3.43, Eq.3.44 | `src/physics/bubble_dynamics.py` | 部分一致（公式存在，默认配置偏离） |
| `K_bd` | Eq.3.50 | `src/physics/mass_transfer.py` | 一致（主式存在） |

## 4) 干燥与热解

| 项 | 论文锚点 | 实现锚点 | 判定 |
|---|---|---|---|
| 干燥 | Eq.4.2, Eq.4.4, Eq.4.6, Eq.4.9 | `src/thermal/drying.py` | 部分一致（数值法工程替代） |
| DAEM | Eq.4.10~4.12 | `src/thermal/devolatilization.py` | 一致（方程与参数映射） |

## 5) 动力学与平衡

| 项 | 论文锚点 | 实现锚点 | 判定 |
|---|---|---|---|
| 气相主反应 | Eq.5.36, Eq.5.45~5.47, Eq.5.51 | `src/kinetics/gas_reactions.py` | 部分一致（编号与工程标签并存） |
| 焦油氧化/重整 | Eq.5.59, Table 5.4 | `src/kinetics/tar_reactions.py` | 不一致（R10 单位口径与参数缩放偏离论文） |
| A1 Gibbs | Appendix A1 | `src/thermodynamics/gibbs_minimizer.py` | 部分一致（求解链存在，未证同构） |

## 6) 使用边界

1. 本矩阵不再引用已清理历史文档，仅保留论文与代码锚点。  
2. “部分一致/不一致”条目，不能写成“已完全复现”。  
3. 同构验收以 `docs/hamel_original_algorithm_reconstruction.md` 与 `specs/hamel_checked_algorithm_spec.md` 联合为准。
