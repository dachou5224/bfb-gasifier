# T6：护栏工具——hamel-source-check skill、引用校验与偏离检查

| 字段 | 内容 |
|---|---|
| 执行者 | skill 文本：总司令起草；两个脚本：cursor |
| 审查方 | codex（脚本与总司令起草的 skill 文本） |
| 前置任务 | T0 合并；引用校验脚本依赖 T2 的 `index.json` 结构 |
| 分支 | `task/T6-guardrails` |
| 时间上限 | 1 个工作日 |
| 状态 | 待开始；`drift_check.py` 须在 T9 之前完成 |

## 目标

把章程 §4 与 §9 从“提醒”变成可执行、可机械检查的约束。

## 范围内

1. **`hamel-source-check` skill**（总司令起草），同时放在 `.claude/skills/hamel-source-check/SKILL.md` 和 `.cursor/skills/hamel-source-check/SKILL.md`。规定原文答卷的固定问题：
   1. 涉及哪些式号、页码、表号；
   2. 每个变量的含义、单位、适用范围；
   3. 德文关键原句；
   4. 依赖哪些上游量；
   5. 原文没有说明、需要裁定的地方；
   6. 与归档旧代码的差异。

   答卷沿用 grilling 的编号格式，由总司令逐题审查。
2. **引用校验脚本** `tools/check_citations.py`（cursor）：
   - 扫描 `hamel/` 下所有 `Source: Hamel Eq.x.y, S. N` 引用，确认式号存在于 `refs/hamel/index.json`、页码一致；
   - 核心模块（`balances/`、`kinetics/`、`hydrodynamics/`、`pyrolysis/`、`freeboard/`）的公开函数缺少引用时报错。
3. **偏离检查脚本** `tools/drift_check.py`（cursor），输出章程 §9 的各项信号：
   - 工作树是否干净；
   - 本卡相对共同祖先（`git merge-base <基准> HEAD`）的 `hamel/` 净增行数，以及 `hamel/` 总行数是否超过上限；新文件名与新增配置字段同样以共同祖先为比较基点；
   - 配置字段数量的变化，以及新增字段是否有出处：求解器字段须在正式 ADR（`docs/adr/NNNN-*.md`）中按完整标识符出现；物理字段在正式 ADR 中出现，**或者**满足 CHARTER v0.6 §9 的例外：是 `data/cases/*.json` 中某个字段的键名，且该字段同时满足——
     - `status` 为 `source` 或 `derived`；
     - `source` 中含可识别的印刷页码 `S. <数字>`（格式见 T7 任务卡）；
     - `status` 为 `derived` 时，`formula` 非空。
     不满足任一条件（如 `status=missing`、`source` 为“待补”等占位文字或无页码、`derived` 缺 `formula`）即不享受例外，仍须正式 ADR；
   - 新文件名是否匹配阶段后缀模式（`_v\d+`、`_r\d+_r\d+`、`_stage\d+` 等）；
   - `docs/STATUS.md` 中工况指标表最近一次变化的日期；
   - skill 副本一致性（2026-10-14 增补）：`.claude/skills/<名称>/SKILL.md` 与 `.cursor/skills/<名称>/SKILL.md` 内容不一致时触发；只存在于其中一边时也触发；两边都没有则不检查。
4. 两个脚本都加入默认测试集，或作为 PR 检查的一部分。

## 范围外

- 不引入 CI 平台配置以外的新依赖。
- 不写业务代码。

## 交付物

- 两份 SKILL.md；`tools/check_citations.py`、`tools/drift_check.py` 及其测试。

## 验收证据

- 用构造的违规样例（缺引用、假式号、阶段后缀文件名、未登记的配置字段）各触发一次，脚本都能报出。
- 配置出处的例外判定，正反样例各至少一个：`source` 无页码、`source` 为占位文字、`derived` 缺 `formula`、`status=missing`，均须报出；`source` 带 `S. <数字>` 的 `source` 字段，以及带页码且有 `formula` 的 `derived` 字段，可以豁免；与输入文件同名的**求解器**字段仍须报出。

## 停止条件

时间上限到期，提交已完成的部分。
