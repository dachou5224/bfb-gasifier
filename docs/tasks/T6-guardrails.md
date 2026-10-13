# T6：护栏工具——hamel-source-check skill、引用校验与偏离检查

| 字段 | 内容 |
|---|---|
| 执行者 | 总司令起草 skill 文本；cursor 写两个脚本 |
| 前置任务 | T0 合并；引用校验脚本依赖 T2 的 `index.json` 结构 |
| 分支 | `task/T6-guardrails` |
| 时间上限 | 1 个工作日 |
| 状态 | 草稿 |

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
   - 相对 main 的净增行数（按任务卡分支）与 `hamel/` 总行数是否超过上限；
   - 配置字段数量的变化，以及新增字段是否有对应 ADR（按字段名在 `docs/adr/` 中检索）；
   - 新文件名是否匹配阶段后缀模式（`_v\d+`、`_r\d+_r\d+`、`_stage\d+` 等）；
   - `docs/STATUS.md` 中工况指标表最近一次变化的日期。
4. 两个脚本都加入默认测试集，或作为 PR 检查的一部分。

## 范围外

- 不引入 CI 平台配置以外的新依赖。
- 不写业务代码。

## 交付物

- 两份 SKILL.md；`tools/check_citations.py`、`tools/drift_check.py` 及其测试。

## 验收证据

- 用构造的违规样例（缺引用、假式号、阶段后缀文件名、未登记的配置字段）各触发一次，脚本都能报出。

## 停止条件

时间上限到期，提交已完成的部分。
