# T4：第 5 章反应动力学的独立转录与比对

| 字段 | 内容 |
|---|---|
| 执行者 | codex |
| 前置任务 | T2（需要原文索引） |
| 分支 | `task/T4-transcribe-kinetics` |
| 时间上限 | 2 个工作日 |
| 状态 | 草稿 |

## 目标

用“独立重做再比对”的方式，确认反应动力学（风险最高的上游输入）的每个常数、级数、单位和适用范围，为迁入 `hamel/kinetics/` 准备审计依据。

## 范围内

1. **独立转录阶段**（完成之前不得打开 `docs/hamel_submodels/` 和归档代码）：只依据 `refs/hamel/` 语料，转录 Hamel 第 5 章全部反应（异相 R1–R4、均相 R5–R11，编号以 Hamel 为准，见 `docs/adr/0001-reaction-numbering.md`）的：
   - 计量式；
   - 速率表达式；
   - 前因子、活化能、各组分反应级数及其单位；
   - 温度、压力适用范围；
   - Hamel 引用的原始文献（按章程 §4.2 标注）；
   - 每项附 PDF 页、印刷页、式号、德文关键原句。

   同时核实 `docs/adr/CANDIDATES.md` 中交给 T4 的条目（C5–C7、C9），结论写进转录文件。

   输出到 `refs/hamel/transcription/ch5_kinetics.md`。
2. **比对阶段**：将转录结果与以下两者逐项比对：
   - `docs/hamel_submodels/` 中对应的转录；
   - 归档代码中的实现（tag `archive-2026-10-12` 下的 `src/kinetics/` 等），包括单位换算和编号映射（对照 `GLOSSARY.md` 中的旧编号表）。
3. 输出差异表：项目、Hamel 原文值（页码）、旧转录值、旧代码值、判定（一致 / 旧转录错 / 旧代码错 / 原文有歧义）。
4. 对每个“原文有歧义”的项，写出可选解释及其影响，作为 ADR 候选，不自行裁定。

## 范围外

- 不写 `hamel/` 下的任何代码。
- 不修改归档代码。
- 不裁定歧义。

## 原文答卷

本卡的转录文件即为答卷，须覆盖 `hamel-source-check` 的六个问题。

## 交付物

- `refs/hamel/transcription/ch5_kinetics.md`
- `refs/hamel/transcription/ch5_kinetics_diff.md`

## 验收证据

- 每个反应的每个参数都有页码和原句。
- 目标所有者抽查约 10%，重点是活化能、前因子、单位和反应级数。

## 停止条件

时间上限到期，或原文某处无法辨认时，停下并列出待人工辨认的位置。
