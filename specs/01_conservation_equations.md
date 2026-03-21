# 01 守恒方程

**来源**：Hamel & Krumm (2001) Section 2.1；Hamel (1999) 技术报告

> **说明**：下文为**单 cell 内**守恒式的文献形式；**全炉联立求解**（全局 NR、分块 Jacobian、Vorabrechnung、连接矩阵）见 **`docs/validation_gap_analysis.md` §0** 与 **`specs/00_implementation_scope.md`**。

## 1. 相间气体交换

$$\dot{N}_{ex,bd,j,i} = K_{bd,i} \cdot V_{b,i} \cdot (C_{j,b,i} - C_{j,d,i})$$

## 2.1 气相摩尔守恒

**悬浮相（d）：**
$$0 = \dot{N}_{zu,d,j,i} + \dot{N}_{rez,d,j,i} + \dot{N}_{d,j,i+1} + \dot{N}_{r,d,j,i} - \dot{N}_{d,j,i} - \dot{N}_{ex,bd,j,i}$$

**气泡相（b）：**
$$0 = \dot{N}_{zu,b,j,i} + \dot{N}_{rez,b,j,i} + \dot{N}_{b,j,i+1} + \dot{N}_{r,b,j,i} - \dot{N}_{b,j,i} + \dot{N}_{ex,bd,j,i}$$

> $\dot{N}_{ex,bd}$ 在悬浮相为负，在气泡相为正，两相使用同一数值。

**反应源项 $\dot{N}_{r}$**：含 R1–R8、R10、R11 动力学源项。当 `Cell.use_gibbs_minor=True` 时，H2S、NH3 等微量组分由 Gibbs 最小化求解，其松弛源项亦计入 $\dot{N}_{r,d}$（见 `specs/04_kinetics.md`）。

## 2.2 固相质量守恒（含粒径类）

$$0 = \sum_k \left[\dot{m}_{zu} + \dot{m}_{rez} - \dot{m}_{aus} + \dot{m}_{auf,i+1} - \dot{m}_{auf,i} + \dot{m}_{ab,i-1} - \dot{m}_{ab,i} + \dot{m}_{r} + \dot{m}_{left} - \dot{m}_{right}\right]_{j,k,i}$$

## 2.3 全局能量守恒

$$0 = \dot{H}_{zu,i} + \dot{H}_{rez,i} + \dot{H}_{auf,i+1} + \dot{H}_{ab,i-1} - \dot{H}_{auf,i} - \dot{H}_{ab,i} - \dot{H}_{aus,i} - \dot{Q}_{W,i} - \dot{Q}_{WU,i}$$

焓值必须包含标准生成焓，使反应热自动体现在守恒中。
