# T26u：`hamel/` 行数清理（第二次，行为逐位不变）

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，在主仓库目录内切任务分支，不建 worktree） |
| 审查方 | codex（代码类，CHARTER §8） |
| 前置任务 | T26t 已合并（PR #167） |
| 分支 | `task/T26u-hamel-line-budget` |
| 时间上限 | 人工工作 1 个工作日。目标：`hamel/` 净减 ≥ 250 行。净减不足 200 行时停下上报，并列出已查过但不能删的项及原因。`tools/` 与 `reports/` 脚本合计净增不超过 +300 行，其中搬入 `tools/` 的模块按搬入后的行数计 |

## 来由

T26t 合并后 `hamel/` 为 10 966 行，距 11 000 行上限（CHARTER v0.14 §6.2）只剩 34 行。工况 1、φ = 1 已经收敛，后续还有 O4（Corella 外推范围）与第 1 关验收要做，都需要行数。目标所有者 2026-11-18 同意另立本卡腾出行数。

本卡**只重构，不改行为**：所有数值结果逐位不变。

## 执行要求（给 cursor）

- 新会话，在主仓库目录 `bfb-gasifier/` 内工作，不建 worktree：确认工作树干净后，`git fetch origin && git switch -c task/T26u-hamel-line-budget origin/main`。
- 不在仓库目录以外创建或写入任何文件；临时文件放 `.scratch/`。
- 会话结束时工作树干净，开 PR，同一 PR 内更新 STATUS 的 T26u 一行。

## 范围内

### A. 动手前：基准快照（在 origin/main 上，先于任何改动）

1. 用现有的 `tools/t26q_golden.py` 在 origin/main 上生成 `.scratch/t26u/golden.npz`（不入库）。基准 SHA 写进 PR 和 `golden.npz` 的元数据。这一步覆盖第一次预计算、初值残差、φ = 0,3 初值构造、外层第 1 轮前 3 次内层迭代，以及 φ = 1 初值构造的拒绝终态。若该工具在现行主线上不能原样运行，只做最小适配，并在 PR 写明。
2. **边界求积路径的补充基准**：用 `.scratch/t26t/` 下已有的两份入参，在 origin/main 上各做一次预计算，保存全部冻结量数组：
   - `precalc_args_r17.pkl`：第 17 轮入参，含第三种签名与端点节点；
   - `final_args.pkl`：T26t 第 D.3 项收敛解的 2(b) 入参。

   同时保存 `account_freeboard` 的 `boundary` 记录。补充基准可以加进 `tools/t26q_golden.py` 的一个新模式，计入 `tools/` 预算。
3. 改动后用同一工具重算，与快照**逐位比较**：用 `np.array_equal`，不用容差；`missing`、`extra`、`different` 都必须为空。任何不等即停下上报。

### B. 候选（逐项查，能删则删，不能删写明原因）

1. **`hamel/solver/campaign.py` 移到 `tools/`**（主要来源，约 −222 行）。这个模块是 ADR 0030 敏感性运行的调度器：登记 S1–S19b 运行表、生成工况副本、起子进程、写报告。它本身已经导入 `tools.driver_report`，`hamel/` 的正式求解路径（`driver.solve_case` 及其下游）不调用它。
   - 移到 `tools/campaign.py`（或并入 `tools/driver_report.py`，择一，PR 写明理由），运行表内容一字不改。
   - 改写全部导入方，包括 `tools/driver_report.py`、`tools/t26c_run.py`、`tools/t26d_validate.py`、`tools/t26e_validate.py`、`tools/t26f_precheck.py`、`tools/t26g_validate.py`，以及 `tests/test_driver.py`、`tests/test_run_prep.py`、`tests/test_bed_ash_split.py`、`tests/test_t26c_diagnosis.py`、`tests/test_sensitivity_params.py` 中的 import 与 monkeypatch 路径。用 `rg "solver.campaign|solver import .*campaign"` 确认没有遗漏。
   - 文档中的 `campaign.py`、`hamel/solver/campaign.py` 路径只改**现行**文档：ADR 0030 冻结清单，以及 STATUS 中“下一步”一类的指引。已接受的 ADR 和历史任务卡正文不改（它们记录当时的路径）；需要时在 ADR 0030 加一行路径变更注记。
   - 搬走后 `hamel/` 不再有对 `tools/` 的**新增**导入。导入检查用 `rg -n '^\s*(from tools\b|import tools\b)' hamel`（含缩进导入），PR 中给出完整输出：应当**只剩**下面一条具名例外。
   - **具名例外（保留，不改）**：`hamel/solver/driver.py` 在 `solve_case` 收敛并通过外部平衡、非负与壁面散热之后，缩进导入 `from tools.run_prep import attach_converged` 并调用它，把记账、闭式比对与加密写进 `outcome.diagnostics`（T26a 起如此，T26q 把相关函数放进 `tools/run_prep.py`）。它在求解结束后执行、不影响求解与收敛判定，但属于 `solve_case` 返回内容的一部分。本卡**不迁移、不删除**这条调用，收敛后的报告内容保持不变；是否迁移另立卡讨论。
2. **`hamel/solver/sensitivity.py` 中只供报告的部分**：逐个函数查调用方。只被 `tools/` 或测试使用、正式路径（residual、precalc、driver）不用的函数，随第 1 项移到 `tools/`；正式路径用到的（如 `diameter_at`、`section_volume_m3`、`pack_reported`）保留。
3. **T26t 留下的压行**：`hamel/freeboard/accounting.py` 中 `_sum_ledgers` 与 `_check_closure` 被写成单行的 `for …: …` 和 `if …: raise …`（commit `96282c8`），按第 C 项规则恢复为正常的多行写法。这会增加几行，计入本卡的净变化。
4. **T26t 实现中的冗余**：`accounting.account_freeboard` 两阶段累加中，只为第一轮写法保留的变量或分支（例如 `plans` 中不再使用的字段、`lookup` 与 `sequence` 的重复构造）能否等价简化。必须保证第 A.2 项逐位相同。
5. **T26q 留下未做的项**：
   - 第 B1 项“合并子集牛顿”（`initial.solve_solids` 与 `inert_balance.free_ash_state`）再查一次，能等价合并则做，不能则写明原因。
   - `hydrodynamics.bed.sphericity_diameter_product`、`ergun_archimedes` 只被测试引用，T26q 曾以“是该式唯一实现”为由保留；若仍是唯一实现，保留并在 PR 写明。
6. **其余自查**：`driver.py`、`precalc.py`、`newton.py`、`residual.py` 中只为已撤销或已替代的功能留下的参数、分支、配置字段和说明。例如 ADR 0057、0058 之前的两态诊断残留；u_min = 0、n_b 某一取值等仅供对照的路径**不算残留**，S19a/S19b 与 ADR 0058 的对照仍要用它们。

### C. 规则

- **不改**：方程、残差、闭合、冻结值、收敛门与容差、求解手段（ADR 0054–0058 的行为）、来源引用（`Source:`、`解释性选择：ADR …`）、`SolverConfig` 字段的名称与默认值。
- **不得用格式手段凑数**：不把多行合并成长行或单行复合语句，不删 docstring 中的来源与 ADR 引用，不删有效注释。行数以 `tools/drift_check.py` 的计数为准；PR 中按文件列出净变化，并逐项写明删改理由。
- 移到 `tools/` 的代码不得被 `hamel/` 的正式路径调用；除第 B.1 项列出的 `driver.py` → `tools.run_prep.attach_converged` 具名例外外，`hamel/` 不得导入 `tools/`。不得以删除或改写该例外调用来满足导入检查。

### D. 验证

- 第 A 项逐位比较全部相等，比较日志摘要写进 PR。
- 默认 pytest、`tools/check_citations.py`、`tools/drift_check.py`。已知失败 `test_cli_skip_render_does_not_rewrite_png` 照旧报告。
- 搬移或删去的函数，其测试的改写或删除逐条写进 PR。
- 不重跑外层完整求解（第 A 项已覆盖预计算与内层路径；外层 8000 s 级运行不属本卡）。

## 范围外（明确禁止）

- 任何行为变化；O4；ADR 0058 的“A…C…A”两次边界局限。
- 改动 `tools/` 中历史运行记录，或重新生成已有报告。
- 修改已接受 ADR、历史任务卡的正文（第 B.1 项允许的 ADR 0030 路径注记除外）。
- 迁移或删除 `driver.py` 对 `tools.run_prep.attach_converged` 的收敛后调用。

## 原文答卷

不适用（不涉及方程、常数或数据的内容变化）。

## 交付物

- 代码改动与测试改写。
- 第 A 项比较结果（PR 正文）、按文件的行数净变化（`hamel/`、`tools/`、`reports/` 分列）、每项候选的处置（删 / 移 / 保留及理由）。
- STATUS 的 T26u 一行。

## 停止条件

到达时间上限，或需要越出范围时，停下并在 PR／`docs/STATUS.md` 中写明做到哪里、卡在哪里、需要什么决策。以下情况也停下上报：

- 逐位比较出现差异；
- `hamel/` 净减不足 200 行；
- `tools/` 与 `reports/` 脚本合计净增超过 +300 行；
- 某项删改需要改变行为、改配置字段，或触及范围外条目。

## 审查记录

- PR #169 审查（codex）第 1 轮（`8e21fa4`）一项 P2：`driver.py` 收敛后缩进导入 `tools.run_prep.attach_converged`，卡上导入检查漏掉缩进导入且未说明这条既有依赖的处置。已修订：导入检查改为含缩进的写法，这条调用列为具名例外，本卡保留不改、不得以删除它满足检查。
- PR #169 复审（codex）第 2 轮（`24081ac`）通过。总司令批准本卡（2026-11-18），可派发。
- PR #170 审查（codex）第 1 轮（`8b5d8fd`）一项 P2：`t26q_golden` 新增的 `boundary` 模式未接硬超时（CHARTER §7）。已改（`c632f95`），超时报 TIMEOUT、不写不完整基准。
- PR #170 复审（codex）第 2 轮（`c632f95`）一项 P2：超时修复后 `tools/` 净增 +316，超过本卡 +300 停止线且 PR 与 STATUS 记录未更新。执行者选择等价化简回到预算内（`d53ef91`），未请求放宽。
- PR #170 复审（codex）第 3 轮（`d53ef91`）通过。总司令独立复核：用 `85576a7` 生成边界求积补充基准（第 17 轮与收敛解两份预计算，160 个数组），再用 PR 代码（`8b5d8fd`，其后只改 `tools/` 与文档）比较，`np.array_equal` 全部相等，missing/extra/different 为空；`campaign.py` 搬移只改子进程模块名与 `_root()` 层级；导入检查只剩 `driver.py` → `tools.run_prep.attach_converged` 具名例外。`hamel/` 10 966 → 10 746（净减 220，未达 250 目标、高于 200 停止线；B2–B6 未删理由接受），`tools/` 净增 +299（上限 +300）。总司令批准（2026-11-19）。
