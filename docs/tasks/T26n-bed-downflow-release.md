# T26n：实现 ADR 0053（床区整体向下流补入干燥与热解失重和随排灰离开的燃料），定位余项

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，在主仓库目录内切任务分支，不建 worktree） |
| 审查方 | codex（代码类与运行记录，CHARTER §8） |
| 前置任务 | ADR 0053 已批准；T26m（PR #141）已合并 |
| 分支 | `task/T26n-bed-downflow-release` |
| 时间上限 | 人工工作 1 个工作日；运行墙钟按 ADR 0037。`hamel/` 净增超过 40 行即停下上报；全部代码净增不超过 120 行（不含测试） |

## 来由

ADR 0009 递推的最底格闭合差约 7,8 g/s（约为燃料量的 60 %），按 ADR 0009 有效域第 5 条与 ADR 0012 第 4 项，即使求解收敛也通不过 CHARTER §2.1。研究底稿 `docs/research/solver-path-alternation-o2.md` 第 5 节表明，闭合差几乎全部来自两项实现漏项：

- 床区干燥与热解的失重没有计入反应汇；
- 随排灰离开的未转化燃料没有计入外排。

两项之外，还余约 2,5·10⁻⁷ kg/s 未定位。ADR 0053 据此修正。

本卡发生在目标所有者 2026-11-02 决定的冻结期内：第一次内层收敛之前冻结全部 F 项与闭合。本卡只补齐 ADR 0009 已写明的项，不改闭合。

## 执行要求（给 cursor）

- **新会话，在主仓库目录 `bfb-gasifier/` 内工作，不建 worktree**：确认工作树干净后，`git fetch origin && git switch -c task/T26n-bed-downflow-release origin/main`。
- **不在仓库目录以外创建或写入任何文件**：临时文件放 `.scratch/`，报告放 `reports/t26n-validation/`。
- **开工前必读**：本卡、ADR 0053、ADR 0009 第 2 项与有效域第 5 条、ADR 0014（对 ADR 0009 的 VTT 适用修订）、ADR 0017 第 13 项、ADR 0052 第 1 项第 2 步、研究底稿第 5 节；`hamel/solver/precalc.py`（`bed_downflow`、`downward_mass_flow`、`_pyro_flows`、`_precalculate` 中递推与 `_release` 的先后）、`hamel/solver/initial.py`（`recompute_flows`）、`hamel/solver/driver.py`（`_precalc`、`solid_feedback`）、`hamel/pyrolysis/release.py`（`ReleaseDistribution`、`_closure`、`_state_closure`）。
- 会话结束时工作树干净，开 PR，同一 PR 内更新 STATUS 的 T26n 一行。

## 范围内

### A. 实现（ADR 0053 第 1–3 项）

1. `bed_downflow` 的反应汇：在焦炭消耗之外，减去各床格干燥释放的水分、热解释放的气体与焦油的质量流。取 `ReleaseDistribution` 中与床区气相源项相同的那组量，不按产率重算。
2. 最底一格的外排：在排灰的灰、焦炭之外，加上随排灰离开的未转化燃料（湿燃料与干燥后燃料的出口流）。
3. 求值所用的状态：
   - 预计算中取上一次预计算的释放分布，与现行反应汇、排灰的回馈方式相同，第一次预计算取 0；
   - `initial.recompute_flows` 取当前预计算的释放分布；
   - 两处调用同一函数，不复制公式。

### B. 定位余项（ADR 0053 第 4 项）

4. 在工况 1 同一状态上，按床格列出整体固体的全部进出项，与 `release` 模块自身的质量闭合对照，找出约 2,5·10⁻⁷ kg/s 的来源。
5. 若是另一处实现漏项，在本卡内修正并在 PR 中说明；若需要新的解释性选择，停下上报。

### C. 报告（只报告）

6. 用 `tools/alternation_diagnostic.py`（研究工具，A5 余额为负时取负值，其余与正式写法相同）在修正后的代码上运行初值构造。
   - 报告写入 `reports/t26n-validation/`：每次交替的最底格闭合差（kg/s，以及除以 ṁ_ref 的比值）、两项补入量、各界面 Ṁ_ab 的修正前后对比、余项定位结果。
   - 工具若须适配修正后的接口，只改工具本身。
7. **判据（预先登记）**：每次交替收敛固相上的最底格闭合差满足 ADR 0017 第 13 项（|残差|/ṁ_ref ≤ 1·10⁻⁶）；否则逐项列出余项并停止。
   - 不运行主求解的内层，不进入整批重跑。

## 范围外（明确禁止）

- 不改 ADR 0009 的闭合、ADR 0014 A5、ADR 0052 的求解路径、冻结值的定义、阻尼、可行条件与收敛门。
- 不在残差或递推中加 clip、下限、开关；不把闭合差覆写为 0。
- 不改 A5 在预计算中的符号检查（待求解路径诊断后裁定）。

## 交付物

- A 的代码与测试；B 的定位结果；C 的报告；`docs/STATUS.md` 中 T26n 一行。

## 验收证据

- `pytest` 默认集通过且 5 分钟内；`check_citations.py`、`drift_check.py` 通过。
- 测试按 ADR 0053“实现约束”：
  - 补入两项后小例子的闭合差为零到舍入；
  - 去掉任一项，闭合差恰为该项；
  - 预计算与 `recompute_flows` 在同一状态下的 Ṁ_ab 逐位相同。
- `hamel/` 净增不超过 40 行；全部代码不超过 120 行（不含测试）。

## 停止条件

判据未达、到达时间或行数上限、发现需要新解释性选择的余项，或需要越出本卡范围时，停下并在 PR / STATUS 中写明。
