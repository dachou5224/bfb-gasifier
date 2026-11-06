# T26q：`hamel/` 行数清理（行为逐位不变）

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，在主仓库目录内切任务分支，不建 worktree） |
| 审查方 | codex（代码类，CHARTER §8） |
| 前置任务 | T26p 已合并（PR #154） |
| 分支 | `task/T26q-hamel-line-budget` |
| 时间上限 | 人工工作 1 个工作日。目标：`hamel/` 净减 ≥ 150 行；净减不足 100 行时停下上报（列出已查而不能删的项及原因）。`tools/` 净增不超过 +150 行 |

## 来由

T26p 合并后 `hamel/` 为 10 886 行，距 11 000 行上限只剩 114 行（CHARTER；目标所有者 2026-11-02 决定第 3 项保留上限）。下一项 ADR 0056（初值构造）及后续工作都要用行数。目标所有者 2026-11-05 决定另立本卡腾出行数。

本卡**只重构，不改行为**：所有数值结果逐位不变。

## 执行要求（给 cursor）

- 新会话，在主仓库目录 `bfb-gasifier/` 内工作，不建 worktree：确认工作树干净后，`git fetch origin && git switch -c task/T26q-hamel-line-budget origin/main`。
- 不在仓库目录以外创建或写入任何文件；临时文件放 `.scratch/`。
- 会话结束时工作树干净，开 PR，同一 PR 内更新 STATUS 的 T26q 一行。

## 范围内

### A. 动手前：基准快照（在 origin/main 上，先于任何改动）

1. 写 `tools/t26q_golden.py`，在工况 1 上生成并保存到 `.scratch/t26q/golden.npz`（不入库）：
   - 第一次预计算的全部冻结量数组（`FrozenQuantities` 各字段）；
   - 初值状态上的完整残差向量；
   - 以 φ = 0,3 的工况副本跑初值构造（ADR 0054 第 3 项）的每次交替：固相牛顿迭代数、收敛固相向量、排灰与 w_f 更新记录（E、E₀、E₁、评价次数）；
   - **全程 φ = 0,3**：上述初值构造结束后，经正式路径（`driver.solve_case`，工况副本 φ = 0,3，`tools` 中外推 Corella，`outer_round_limit = 1`）进入外层第 1 轮，记录内层前 3 次迭代的归一残差范数、y 上的步向量与阻尼试探记录（达到 3 次即用期限或包装函数中止，不改正式路径）。不做 φ = 0,3 → 1 的切换（那是 ADR 0056 的内容，PR #155 审查）；
   - **φ = 1 的初值构造拒绝终态**（工况 1 原文件）：`InitialStateError.diagnostics` 的全部内容，含 #154 补齐的被拒解五组存量数组、逐格惰性物、排灰、归一残差、`min_stock_kg` 与第 3 层 detail。
   - 快照所用的基准提交（origin/main 的 SHA）写进 PR 与 `golden.npz` 的元数据。
2. 改动后用同一工具重算，与快照**逐位比较**（`np.array_equal`，不是容差）；任何不等即停下上报。

### B. 候选（逐项查，能删则删，不能删写明原因）

1. **合并两处“子集牛顿”**：`initial.solve_solids` 与 `inert_balance.free_ash_state` 都是“固定其余未知量、只取一部分行与未知量、调用 `damped_newton`”。抽成一个共用函数（选行、嵌入、缩放、门、负存量检查），两处调用。
2. **只供报告或工具的函数移出 `hamel/`**：
   - `residual.particle_temperatures`（只有 `tools/` 使用）移到 `tools/`；
   - 初值构造**成功交替记录**（`initial._record`）中的 ADR 0052 遗留报告项 `step_history`、`volume_ratio`、`min_stock_kg` 可删，`_volume_ratio` 与 `SolidSolve.step_history` 随之移到 `tools/` 或删去，并改 `tools/initial_solids_diagnostic.py`、`tools/t26l_validate.py` 的导入；
   - **保留**被拒解的证据（PR #154 审查第 4 项补齐，PR #155 审查）：`_rejected()` 写入 `InitialStateError.diagnostics` 的五组存量数组、逐格惰性物、排灰、归一残差、`min_stock_kg` 与第 3 层 detail 一项不少。`_minimum_stock` 仍被 `_rejected()` 使用，留在 `hamel/`（可在 `hamel/` 内等价合并），`hamel/` 不得反向导入 `tools/`。
3. **只被测试引用的函数**：`hydrodynamics.bed.sphericity_diameter_product`、`ergun_archimedes`：若对应 Hamel 式号的计算已在正式路径中由别的函数完成，则删去并把测试改为测正式函数；若它们是某式的唯一实现，保留并在 PR 写明。
4. `newton._dump`（失败点 npz）：保留功能，但查能否与 `_note` 共用数组整理，压缩重复。
5. `inert_balance.py`：`_text`、`_Exhausted` 与报告字典的组装能否简化；`InertBalanceError` 能否复用现有错误类型。
6. 其余自查：`driver.py`、`precalc.py`、`initial.py` 中只为已删功能留下的参数、分支与说明（T26o、T26p 后的残留）。

### C. 规则

- **不改**：方程、残差、闭合、冻结值、收敛门与容差、求解手段（ADR 0054、0055 的行为）、来源引用（`Source:`、`解释性选择：ADR …`）。
- **不得用格式手段凑数**：不把多行合并成长行、不删 docstring 中的来源与 ADR 引用、不删有效注释。行数以 `tools/drift_check.py` 的计数为准；PR 中按文件列出净变化与每项删改理由。
- 移到 `tools/` 的函数不得被 `hamel/` 的正式路径调用。

### D. 验证

- 第 A 项逐位比较全部相等（报告写入 PR）。
- 默认 pytest、`tools/check_citations.py`、`tools/drift_check.py`；已知失败 `test_cli_skip_render_does_not_rewrite_png` 照旧报告。
- 被移出或删去的函数，其测试改写或删除逐条写进 PR。

## 范围外（明确禁止）

- 任何行为变化；ADR 0056 的内容（初值构造）；O4。
- 改 `tools/` 中历史运行记录或重生成已有报告。

## 停止条件

- 逐位比较出现差异；
- 净减不足 100 行；
- 某项删改需要改变行为或范围外条目。

## 审查记录

- PR #155 审查（codex）：两项 P2——基准中 φ 切换无唯一定义（改为全程 φ = 0,3，另记 φ = 1 初值构造拒绝终态）；被拒解证据不得删除（只删成功交替记录中的旧报告项）。复审通过。
- **2026-11-05：总司令批准**，可派发。
- **2026-11-06：`tools/` 超限，目标所有者批准。** 本卡预算 `tools/` 净增不超过 +150，实际 +278（`tools/t26q_golden.py` +231，`tools/run_prep.py` +43，`tools/initial_solids_diagnostic.py` +4）。超限来自第 A 项要求写的基准工具，卡上预算未给它留足行数，属任务卡的疏漏。保留该工具，ADR 0056 实现中 φ ≤ 0,3 逐位不变的验证可复用；`hamel/` 不导入它。PR #157 审查（codex）的 P1 由此关闭。
- PR #157 审查（codex）：两轮。第 1 轮一项 P1（`tools/` 超预算，见上条）。第 2 轮（`98d51c9`）通过；审查方核对了执行者的 407 项相等日志，未重跑长快照。
- 总司令独立复核（2026-11-06）：在 `.scratch` 的导出副本中，用 `a08bced` 的代码重新生成基准，再用 PR 代码（`c4822fd`，其后只有文档提交）比较：407 个 `data.*` 数组 `np.array_equal` 全部相等，`missing`、`extra`、`different` 均为空。重新生成的基准与执行者的 `golden.npz` 也逐位相等。
- **2026-11-06：总司令批准**。`hamel/` 10 886 → 10 733，净减 153；第 B1 项（合并子集牛顿）未做，理由见 PR，接受。
