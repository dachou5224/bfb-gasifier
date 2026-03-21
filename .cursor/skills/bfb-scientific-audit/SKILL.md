---
name: bfb-scientific-audit
description: BFB 一维模型代码审查专用技能。检查 SI 单位、arrhenius.py 工厂函数、变量单位注释与 Hamel (1999) 方程溯源。在 Validate 阶段或执行 /review 类审查前使用。
---

# bfb-scientific-audit（科学建模审计）

## 何时启用

- **Execution 的 Validate 阶段**：合并或宣称「完成」前，对本次改动做科学一致性门禁。
- **使用 gstack `/review`（或任意深度代码审查）之前**：先跑本清单，避免单位与动力学形式类问题漏到通用审查中。

## 职责清单（必须逐项核对）

### 1. 单位制（SI）

- 全链路是否为 **SI**：**Pa, K, m, mol, kg, s**（见 `docs/CLAUDE.md` §全局编程规则）。
- **已知例外**：R8 WGSR 速率式中压力项为 **atm** 时，是否在函数**内部**完成 `P_atm = P_Pa / 101325.0`，对外仍传 **Pa**。
- 新增/修改的 **函数参数与关键局部变量** 是否在定义处或 docstring 中标注单位（可与 `.cursor/rules/scientific-modeling-core.mdc` 对齐）。

### 2. Arrhenius 与反应速率

- 所有需 Arrhenius 形式的速率常数是否**仅**通过 `src/kinetics/arrhenius.py` 中的工厂函数构造：
  - `k_hobbs`（R1–R4 等形式 A）
  - `k_standard`（形式 C）
  - `k_jensen_r7`（R7 等形式 B）
- **禁止**：在其它文件中手写 `k = A*exp(-E/(R*T))` 实现 R1–R4 等与规范冲突的形式（见 `docs/CLAUDE.md` §禁止事项）。

### 3. 物理常数

- `Rg`、`g`、`P0`、`n_b` 等是否来自**单一约定**（见 `docs/CLAUDE.md`），避免在多处硬编码不一致。

### 4. 文献与方程溯源

- 对 **气泡动力学、相间传质、守恒方程、关键速率式** 等，是否在 docstring 或紧邻注释中标明：
  - **Hamel (1999)** 方程编号（如 Eq. 3.50、4.4/4.5），或
  - **Hamel & Krumm (2001)** / `specs/*.md` 章节。
- 若实现与论文/规格有**已知差异**（如全局 NR vs GS），是否指向 `docs/validation_gap_analysis.md` 或 `docs/hamel_dissertation_vs_python_architecture.md` 相应章节。

### 5. 交叉引用

- 数量级与单位约定是否与 **`docs/source_units_audit.md`** 一致（若本次改动触及输入输出或导出量）。

## 输出格式建议

审查结束后用简短表格汇总：

| 检查项 | 结果（通过/待修/不适用） | 备注（文件:行） |
|--------|-------------------------|-----------------|
| SI + 例外 R8 | | |
| arrhenius 工厂 | | |
| 单位注释 | | |
| Hamel/规格引用 | | |

## 相关文件

- `docs/CLAUDE.md`（全局规则与 Phase）
- `src/kinetics/arrhenius.py`
- `.cursor/rules/scientific-modeling-core.mdc`
- `docs/source_units_audit.md`
