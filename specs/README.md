# `specs/` 方程规格说明

本目录下的 **`01`–`04` 与 `species.md`** 整理 **文献（Hamel & Krumm / Hamel 1999）中的控制方程、符号与物性约定**，供实现与测试对照。

## 与 `src/` 的关系

- **`specs/`**：回答「**方程长什么样**」— 守恒式、源项接口、水力学/动力学/干燥热解公式。
- **数值求解策略**（全局 NR vs 扫描、`fsolve`、Vorabrechnung、连接矩阵、循环在 Jacobian 内/外）**不在**此重复，见：
  - **`docs/validation_gap_analysis.md` §0**（论文 vs 代码分层对照）
  - **`specs/00_implementation_scope.md`**（一页摘要）
  - **`docs/BFB_TechSpec_v11.md` §5.5–§5.6**（目录与 Bild 2.2 差异）
- **论文程序结构（Bild 2.2）与单格摩尔守恒（Bild 2.3）**：Mermaid 流程图已写在根目录 **`README.md`**、**`docs/techspec.md`**、**`docs/BFB_TechSpec_v11.md` §0** 与 **`docs/hamel_dissertation_vs_python_architecture.md` §0**（与 `01_conservation_equations.md` 中 Eq. 3.50 / 两相摩尔守恒一致）。

## 文件索引

| 文件 | 内容 |
|------|------|
| `01_conservation_equations.md` | 气/固/能量守恒，相间交换 |
| `02_hydrodynamics.md` | 最小流化、气泡、传质等 |
| `03_drying_devolatilization.md` | 干燥与 DAEM 热解 |
| `04_kinetics.md` | R1–R11 与 Gibbs 修正角色 |
| `species.md` | 物种列表与约定 |
