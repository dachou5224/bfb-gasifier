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
   - 上述初值构造结束状态上的 `damped_newton` 内层前 3 次迭代（φ = 1 的工况，`tools` 中外推 Corella），每次的残差范数与步向量。
2. 改动后用同一工具重算，与快照**逐位比较**（`np.array_equal`，不是容差）；任何不等即停下上报。

### B. 候选（逐项查，能删则删，不能删写明原因）

1. **合并两处“子集牛顿”**：`initial.solve_solids` 与 `inert_balance.free_ash_state` 都是“固定其余未知量、只取一部分行与未知量、调用 `damped_newton`”。抽成一个共用函数（选行、嵌入、缩放、门、负存量检查），两处调用。
2. **只供报告或工具的函数移出 `hamel/`**：`residual.particle_temperatures`（只有 `tools/` 使用）；`initial._volume_ratio`、`initial._minimum_stock`、`SolidSolve.step_history` 及其报告项（ADR 0052 遗留，只被 `tools/initial_solids_diagnostic.py`、`tools/t26l_validate.py` 使用）——移到 `tools/` 并改这些工具的导入，正式路径的记录里不再带。
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

（待审查）
